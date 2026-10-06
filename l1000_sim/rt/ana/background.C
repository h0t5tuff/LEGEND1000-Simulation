//  RT backgrounds: Tl-208 and Bi-214 decaying in the re-entrant tube wall, as the Ge array sees them. one run per
//  chain, reweighted to the tube as built and to every steel (top) / Cu / EFCu (bottom) design with seams L1, L2
//    root -l -b -q 'ana/background.C("output/tl208*.root=Tl208,output/bi214*.root=Bi214")'
//  no detector model: the Ge energy is what Geant4 deposited, and the only cut is M1 (one detector fired)

#include "../sim/rt.h"
#include "TRegexp.h"
#include <sstream>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  // MAJORANA's background window [Haufe et al., arXiv:2209.10592]: 1950-2350 keV minus 10 keV around Qbb, 2103.5
  // (Tl208 single escape), 2118.5 and 2204.1 keV (Bi214) [Arnquist et al., arXiv:2207.07638]. the BI is its counts / 360 keV
  const double winLo = 1950, winHi = 2350, winGap[4] = {2039.0, 2103.5, 2118.5, 2204.1}, winGapHalf = 5.0;
  bool inWindow(double e) { if (e < winLo || e > winHi) return false; for (double g : winGap) if (fabs(e - g) <= winGapHalf) return false; return true; }
  const double winWidth = winHi - winLo - 4 * 2 * winGapHalf; // 360 keV
  const double m1_keV = 5.0;                  // a detector fired above this
  const double geMass_kg = 1000.0;            // the array the BI is per kg of
  const double bgGoal = 1e-5;                 // LEGEND-1000 target [cts/(keV kg yr)] at Qbb
  const double rtBudget = 1e-5;               // the tube's share of it: a design passes when its MC 90% bound after M1 is below this
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000;        // spectrum [keV]
  const int    eBins = 200;
  const int    NDB = 10;                      // depth slabs of the hit table
  const int    NZ = 625;                      // ~1 cm depth bins: where each section's decays lie

  const char  *mat[3] = {"steel", "Cu", "EFCu"};      // top to bottom
  const double density[3] = {7900.0, 8960.0, 8930.0}; // kg/m^3
  struct Act { double v; bool ul; };                  // specific activity [uBq/kg], only an upper limit?
  const Act act[3][2] = {{{2500, false}, {1000, false}},    // steel: Bi214, Tl208. Ralph's materialMix
                         {{1, false}, {1, false}},          // Cu: the same source
                         {{0.19, false}, {0.077, true}}};   // EFCu: Edgar's survival_BI.py radioassay
  int iso(const std::string &s) { return s.find("Bi") != std::string::npos ? 0 : 1; }
  double activity(const std::string &s, int m) { return act[m][iso(s)].v * 1e-6; } // -> Bq/kg
  bool upperLimit(const std::string &s, int m) { return act[m][iso(s)].ul; }
}

//-------------------------------------------------------------------------------
//  2. One run, reduced to what a reweighting needs:
struct Hit { Long64_t ev; float e; bool m1; }; // one detector's energy in one decay; ev: the decay's slot among those with a hit

struct Run
{
  std::string file, iso;
  Long64_t nsim = 0;
  long nFiles = 1, nBad = 0;            // files merged in; decays outside the wall (a confinement bug)
  long nPv[3] = {0, 0, 0};              // decays per physical volume: 0 mother (EFCu), 1 OFHC, 2 SS
  long nZPv[3][cfg::NZ] = {{0}};        // the same, per ~1 cm of depth
  double mcDensity[3] = {0, 0, 0};      // simulated decays per m^3 in each physical volume, measured
  long dN[cfg::NDB] = {0}, dH[cfg::NDB] = {0}, dW[cfg::NDB] = {0}; // per depth slab: decays, hits, window hits
  std::vector<float> depth;             // per decay WITH a Ge hit: mm below the top of the tube
  std::vector<signed char> pv;          // the physical volume it was drawn in
  std::vector<short> win0, win1;        // its hits in the window: no cut, M1
  std::vector<Hit> hits;
  long nHit[2] = {0, 0};                // hits: no cut, M1
};

static int slabOf(double dep, double L1, double L2) { return dep < L1 ? 0 : (dep < L2 ? 1 : 2); } // a design's three slabs, from the top
static int bin(double dep, double tot, int n) { return std::min(n - 1, std::max(0, (int)(dep / tot * n))); }

static Run readRun(const std::string &file, const std::string &iso, const RT &rt)
{
  Run r;
  r.file = file;
  r.iso = iso;
  auto f = TFile::Open(file.c_str());
  auto d = (f && !f->IsZombie()) ? (TDirectory *)f->Get("stp") : nullptr;
  auto vtx = d ? (TTree *)d->Get("vtx") : nullptr;
  if (!vtx) { printf("ERROR: no stp/vtx in %s\n", file.c_str()); return r; }
  const Long64_t n = r.nsim = vtx->GetEntries();
  std::vector<float> depth(n, -1);
  std::vector<short> nFired(n, 0), win0(n, 0), win1(n, 0);
  std::vector<signed char> pv(n, -1);
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it()) // each Ge detector's steps summed per decay
  {
    auto t = rtIsGermanium(k->GetName()) ? (TTree *)d->Get(k->GetName()) : nullptr;
    if (!t || !t->GetBranch("edep_in_keV")) continue;
    std::map<int, float> sum;
    rtScan(t, {"evtid", "edep_in_keV"}, [&](const double *v) { sum[(int)v[0]] += v[1]; });
    for (auto &p : sum)
      if (p.first >= 0 && p.first < n && p.second > cfg::m1_keV) { nFired[p.first]++; r.hits.push_back({p.first, p.second, false}); }
  }
  rtScan(vtx, {"evtid", "xloc_in_m", "yloc_in_m", "zloc_in_m"}, [&](const double *v) {
    int ev = (int)v[0]; // rows are not in event order with -t: index by evtid
    if (ev < 0 || ev >= n) return;
    double x = v[1] * 1000, y = v[2] * 1000, z = v[3] * 1000;
    depth[ev] = rt.zTop - z;
    pv[ev] = rt.physAt(sqrt(x * x + y * y), z); // the volume remage drew it from
    if (pv[ev] < 0) r.nBad++;
  });
  const double tot = rt.zTop - rt.zBottom;
  for (auto &h : r.hits)
  {
    h.m1 = nFired[h.ev] == 1;
    r.nHit[0]++;
    r.nHit[1] += h.m1;
    if (cfg::inWindow(h.e)) { win0[h.ev]++; win1[h.ev] += h.m1; }
  }
  for (Long64_t i = 0; i < n; i++) // every decay, as counts
  {
    r.dN[bin(depth[i], tot, cfg::NDB)]++;
    r.dW[bin(depth[i], tot, cfg::NDB)] += win0[i];
    if (pv[i] >= 0) { r.nPv[pv[i]]++; r.nZPv[pv[i]][bin(depth[i], tot, cfg::NZ)]++; }
  }
  std::vector<Long64_t> slot(n, -1); // a decay with a hit gets one slot, and its hits point there
  for (auto &h : r.hits)
  {
    r.dH[bin(depth[h.ev], tot, cfg::NDB)]++;
    if (slot[h.ev] < 0)
    {
      slot[h.ev] = r.depth.size();
      r.depth.push_back(depth[h.ev]);
      r.pv.push_back(pv[h.ev]);
      r.win0.push_back(win0[h.ev]);
      r.win1.push_back(win1[h.ev]);
    }
    h.ev = slot[h.ev];
  }
  delete f;
  return r;
}

//-------------------------------------------------------------------------------
//  3. Combine the runs: files of one isotope are one run split into jobs
static std::vector<std::string> expand(const std::string &pat) // "output/tl208*.root" -> every match, sorted
{
  if (pat.find_first_of("*?") == std::string::npos) return {pat};
  size_t sl = pat.find_last_of('/');
  std::string dir = sl == std::string::npos ? "." : pat.substr(0, sl), base = pat.substr(sl + 1);
  TRegexp re(base.c_str(), kTRUE); // kTRUE: a shell wildcard
  std::vector<std::string> out;
  void *d = gSystem->OpenDirectory(dir.c_str());
  for (const char *e; d && (e = gSystem->GetDirEntry(d));)
  {
    TString f(e);
    Ssiz_t len = 0;
    if (f.Index(re, &len) == 0 && len == f.Length()) out.push_back(sl == std::string::npos ? f.Data() : dir + "/" + f.Data());
  }
  if (d) gSystem->FreeDirectory(d);
  std::sort(out.begin(), out.end());
  if (out.empty()) printf("no file matches %s\n", pat.c_str());
  return out;
}

void background(const char *runs = "output/tl208*.root=Tl208,output/bi214*.root=Bi214", const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  if (!rt.ok) return;
  std::vector<Run> R; // one per isotope
  std::stringstream list(runs);
  for (std::string tok; std::getline(list, tok, ',');) // "a.root=Tl208,b.root=Bi214"
  {
    size_t eq = tok.find('=');
    std::string iso = eq == std::string::npos ? "Tl208" : tok.substr(eq + 1);
    if (iso.find("Tl") == std::string::npos && iso.find("Bi") == std::string::npos) { printf("ERROR: unknown isotope '%s'\n", iso.c_str()); continue; }
    for (auto &f : expand(tok.substr(0, eq)))
    {
      Run r = readRun(f, iso, rt);
      if (!r.nsim) continue;
      Run *into = nullptr;
      for (auto &m : R) if (m.iso == iso) into = &m;
      if (!into) { R.push_back(std::move(r)); continue; }
      Long64_t off = into->depth.size(); // slots renumbered, so nothing counts twice
      for (auto h : r.hits) { h.ev += off; into->hits.push_back(h); }
      into->depth.insert(into->depth.end(), r.depth.begin(), r.depth.end());
      into->pv.insert(into->pv.end(), r.pv.begin(), r.pv.end());
      into->win0.insert(into->win0.end(), r.win0.begin(), r.win0.end());
      into->win1.insert(into->win1.end(), r.win1.begin(), r.win1.end());
      for (int k = 0; k < 2; k++) into->nHit[k] += r.nHit[k];
      for (int p = 0; p < 3; p++) { into->nPv[p] += r.nPv[p]; for (int b = 0; b < cfg::NZ; b++) into->nZPv[p][b] += r.nZPv[p][b]; }
      for (int b = 0; b < cfg::NDB; b++) { into->dN[b] += r.dN[b]; into->dH[b] += r.dH[b]; into->dW[b] += r.dW[b]; }
      into->nsim += r.nsim; into->nBad += r.nBad; into->nFiles++; into->file += "+" + r.file;
    }
  }
  if (R.empty()) { printf("no usable runs in '%s'\n", runs); return; }

  const double depthTot = rt.zTop - rt.zBottom, L1ks = rt.zTop - rt.seamSS, L2ks = rt.zTop - rt.seamOFHC; // KS seams as depths
  const double UL90 = 2.30, binW = (cfg::eHi - cfg::eLo) / cfg::eBins;
  double physV[3]; // m^3: mother, OFHC shell, SS shell
  for (int p = 0; p < 3; p++) physV[p] = rt.physVolume(p) * 1e-9;
  for (auto &r : R) for (int p = 0; p < 3; p++) r.mcDensity[p] = r.nPv[p] / physV[p]; // remage fills the mother sparser than its daughters: measure, never assume
  // one simulated decay in volume p, where material m sits = rho x A x 1 yr real decays per m^3, over the MC density, per kg of Ge
  auto weight = [&](const Run &r, int m, int p) { return p < 0 || r.mcDensity[p] <= 0 ? 0.0 : cfg::density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[p] * cfg::geMass_kg); };
  auto massOf = [&](double d0, double d1, int m) { return d1 > d0 ? cfg::density[m] * rt.wallVolume(rt.zTop - d1, rt.zTop - d0) * 1e-9 : 0.0; };
  auto ulN = [](long n) { const double t[11] = {2.30, 3.89, 5.32, 6.68, 7.99, 9.27, 10.53, 11.77, 12.99, 14.21, 15.41}; return n <= 10 ? t[n] : n + 1.28 * sqrt((double)n) + 1; }; // Poisson 90% upper limit
  printf("geometry : %s\nRT wall  : %.5f m^3 over %.3f m\n", rt.file.c_str(), rt.wallVolume(rt.zBottom, rt.zTop) * 1e-9, depthTot / 1000);
  printf("window   : %.0f-%.0f keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's BEW): %.0f keV\n", cfg::winLo, cfg::winHi, cfg::winWidth);
  printf("cut      : M1 (exactly one detector above %.0f keV). no detector model: deposited energy\n\n", cfg::m1_keV);

  for (auto &r : R)
  {
//-------------------------------------------------------------------------------
//  4. Hits, per run:
    long w0 = 0, w1 = 0;
    for (size_t i = 0; i < r.win0.size(); i++) { w0 += r.win0[i]; w1 += r.win1[i]; }
    printf("=== %s   (%s, %lld decays)\n", r.iso.c_str(), r.file.c_str(), r.nsim);
    printf("[4] hits: %ld, after M1 %ld; in the window: %ld, after M1 %ld%s\n", r.nHit[0], r.nHit[1], w0, w1, r.nFiles > 1 ? Form("   (merged from %ld files)", r.nFiles) : "");
    if (r.nBad) printf("    WARNING: %ld decays outside the wall - a confinement bug\n", r.nBad);

//-------------------------------------------------------------------------------
//  5. Normalisation, per physical volume as built:
    const int asBuilt[3] = {2, 1, 0}; // mother EFCu, OFHC shell Cu, SS shell steel
    const char *pvName[3] = {"mother", "OFHC", "SS"};
    printf("[5] normalisation:  %-7s %-6s %9s %8s %10s %9s %11s %14s\n", "volume", "mat", "V [m^3]", "M [kg]", "A [Bq]", "MC dec.", "MC per m^3", "decays/yr per");
    for (int p = 0; p < 3; p++)
    {
      int m = asBuilt[p];
      double mass = physV[p] * cfg::density[m], A = mass * cfg::activity(r.iso, m);
      printf("                    %-7s %-6s %9.5f %8.1f %10.3e %9ld %11.0f %14.3g\n", pvName[p], cfg::mat[m], physV[p], mass, A, r.nPv[p], r.mcDensity[p],
             r.nPv[p] ? A * cfg::secPerYear / r.nPv[p] : 0.0);
    }
    double dAvg = 0.5 * (r.mcDensity[1] + r.mcDensity[2]);
    printf("    sampling density mother / shells = %.3f   (1.000 would be uniform; the weights use the measured density)\n", dAvg > 0 ? r.mcDensity[0] / dAvg : 0.0);

//-------------------------------------------------------------------------------
//  6. Where the hits come from, and the statistics the seams need:
    printf("[6] depth [m]          decays     hits    hits/decay   window\n");
    for (int b = 0; b < cfg::NDB; b++)
      printf("    %5.2f..%-8.2f %10ld %8ld %13.2e %8ld%s\n", depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000, r.dN[b], r.dH[b],
             r.dN[b] ? (double)r.dH[b] / r.dN[b] : 0.0, r.dW[b], b == cfg::NDB - 1 ? "   <- nearest the detectors" : "");
    printf("    a slab with no window hit is bounded at 2.30 x one decay's weight. to bound steel there at 10%% of the budget:\n");
    for (int p = 2; p >= 0; p--)
    {
      double b = UL90 * weight(r, 0, p) / cfg::winWidth;
      printf("      steel over %-7s %9ld decays now: bound %.1e, needs %.1e decays in it\n", pvName[p], r.nPv[p], b, r.nPv[p] * b / (0.1 * cfg::rtBudget));
    }
    printf("\n");
  }

//-------------------------------------------------------------------------------
//  7. Background index as built, per chain and section:
  printf("[7] background index as built [cts/(keV kg yr)]: BI +- MC statistics; a section with no window hit gets its 90%% limit\n");
  printf("    %-6s %-6s %9s   %-26s %-26s\n", "chain", "sect", "win hits", "no cut", "M1");
  double gB[2] = {0, 0}, gV[2] = {0, 0};
  bool gUL = false;
  for (auto &r : R)
  {
    double tB[2] = {0, 0}, tV[2] = {0, 0};
    bool tUL = false; // the total rests on an upper-limit activity
    for (int s = 0; s < 3; s++) // the slabs as built: steel, Cu, EFCu
    {
      double B[2] = {0, 0}, V[2] = {0, 0}, wbar = 0;
      long n[2] = {0, 0}, nsec = 0;
      const double edge[4] = {0, L1ks, L2ks, depthTot};
      for (int p = 0; p < 3; p++) // mean weight of a decay in this section, for the limit
        for (int b = bin(edge[s], depthTot, cfg::NZ); b < std::min(cfg::NZ, (int)ceil(edge[s + 1] / depthTot * cfg::NZ)); b++)
        { nsec += r.nZPv[p][b]; wbar += r.nZPv[p][b] * weight(r, s, p) / cfg::winWidth; }
      wbar = nsec ? wbar / nsec : 0;
      for (size_t i = 0; i < r.depth.size(); i++)
      {
        if (slabOf(r.depth[i], L1ks, L2ks) != s) continue;
        double w = weight(r, s, r.pv[i]) / cfg::winWidth;
        short k[2] = {r.win0[i], r.win1[i]};
        for (int c = 0; c < 2; c++) { B[c] += k[c] * w; V[c] += k[c] * w * w; n[c] += k[c]; }
      }
      const bool ul = cfg::upperLimit(r.iso, s);
      std::string cell[2];
      for (int c = 0; c < 2; c++)
      {
        cell[c] = n[c] ? Form("%s%.2e +- %.1e", ul ? "<" : "", B[c], sqrt(V[c])) : Form("< %.2e (90%%)", UL90 * wbar);
        if (n[c]) { tB[c] += B[c]; tV[c] += V[c]; }
      }
      tUL |= ul && n[1];
      printf("    %-6s %-6s %4ld/%-4ld   %-26s %-26s%s\n", s ? "" : r.iso.c_str(), cfg::mat[s], n[1], n[0], cell[0].c_str(), cell[1].c_str(), ul ? "  activity is an upper limit" : "");
    }
    printf("    %-6s %-6s %9s   %-26s %-26s\n", "", "tube", "", Form("%s%.2e +- %.1e", tUL ? "<" : "", tB[0], sqrt(tV[0])), Form("%s%.2e +- %.1e", tUL ? "<" : "", tB[1], sqrt(tV[1])));
    gUL |= tUL;
    for (int c = 0; c < 2; c++) { gB[c] += tB[c]; gV[c] += tV[c]; }
  }
  printf("    ALL CHAINS          %s%.2e +- %.1e (%s%.2f x goal)   M1: %s%.2e +- %.1e (%s%.2f x goal)\n", gUL ? "<" : "", gB[0], sqrt(gV[0]), gUL ? "<" : "", gB[0] / cfg::bgGoal,
         gUL ? "<" : "", gB[1], sqrt(gV[1]), gUL ? "<" : "", gB[1] / cfg::bgGoal);
  printf("    sections with no window hit are left out of the totals; win hits reads M1/no cut\n");

//-------------------------------------------------------------------------------
//  8. Three-material designs: the MC alone, after M1, judged by its 90% upper bound
  // per slab: n window hits count UL(n) x their mean weight; none count 2.30 x the heaviest decay in the slab
  auto judge = [&](double L1, double L2, double &bi) {
    const double edge[4] = {0, L1, L2, depthTot};
    double up = 0;
    bi = 0;
    for (int s = 0; s < 3; s++)
    {
      double B = 0, wmax = 0;
      long n = 0;
      int b0 = bin(edge[s], depthTot, cfg::NZ), b1 = std::min(cfg::NZ, (int)ceil(edge[s + 1] / depthTot * cfg::NZ));
      for (auto &r : R)
      {
        for (int p = 0; p < 3; p++)
          for (int b = b0; b < b1; b++)
            if (r.nZPv[p][b]) { wmax = std::max(wmax, weight(r, s, p) / cfg::winWidth); break; }
        for (size_t i = 0; i < r.depth.size(); i++)
          if (r.win1[i] && slabOf(r.depth[i], L1, L2) == s) { B += r.win1[i] * weight(r, s, r.pv[i]) / cfg::winWidth; n += r.win1[i]; }
      }
      bi += B;
      up += n ? ulN(n) * B / n : UL90 * wmax;
    }
    return up;
  };
  auto passes = [&](double L1, double L2) { double bi; return judge(L1, L2, bi) <= cfg::rtBudget; };
  auto deepestSteel = [&](double L2) { double l1 = -1; for (double L1 = 10; L1 <= L2 - 10 && passes(L1, L2); L1 += 10) l1 = L1; return l1; }; // 1 cm grid, always some Cu
  printf("\n[8] designs: steel to L1, Cu to L2, EFCu below. MC alone, after M1; passes if its 90%% bound <= the budget %.0e\n", cfg::rtBudget);
  printf("    %-7s %9s %9s %10s %9s %11s %11s %9s\n", "L2 [m]", "EFCu [kg]", "steel to", "steel [kg]", "Cu [kg]", "BI", "90% bound", "x budget");
  std::vector<double> rows = {L2ks};
  for (double L2 = 3000; L2 < depthTot; L2 += 250) rows.push_back(L2);
  std::sort(rows.begin(), rows.end());
  double L2best = -1, L1best = -1;
  for (double L2 : rows)
  {
    double l1 = L2 == L2ks ? L1ks : deepestSteel(L2), bi = 0;
    if (l1 < 0) { printf("    %-7.2f %9.0f %9s   no steel can be shown to pass yet\n", L2 / 1000, massOf(L2, depthTot, 2), "-"); continue; }
    double up = judge(l1, L2, bi);
    printf("    %-7.2f %9.0f %7.2f m %10.0f %9.0f %11.2e %11.2e %9.2f%s\n", L2 / 1000, massOf(L2, depthTot, 2), l1 / 1000, massOf(0, l1, 0),
           massOf(l1, L2, 1), bi, up, up / cfg::rtBudget, L2 == L2ks ? "   <- KS as built" : "");
    if (L2 != L2ks && up <= cfg::rtBudget && L2 > L2best) { L2best = L2; L1best = l1; }
  }
  if (L2best > 0) printf("    least EFCu on this grid, then most steel: L2 %.2f m, L1 %.2f m\n", L2best / 1000, L1best / 1000);
  else printf("    no design with steel passes on the MC alone: [6] says how many decays its slabs need\n");
  printf("    masses in kg; steel includes the lid, which sits in the top slab\n");

//-------------------------------------------------------------------------------
//  9. Draw: the spectrum as built, and the 90% bound of every design
  TH1D *sp[2];
  const char *cutName[2] = {"no cut", "M1"};
  for (int c = 0; c < 2; c++) sp[c] = new TH1D(Form("sp%d", c), "RT background in the germanium, as built;energy [keV];cts / (keV kg yr)", cfg::eBins, cfg::eLo, cfg::eHi);
  for (auto &r : R)
    for (auto &h : r.hits)
    {
      double w = weight(r, slabOf(r.depth[h.ev], L1ks, L2ks), r.pv[h.ev]) / binW;
      sp[0]->Fill(h.e, w);
      if (h.m1) sp[1]->Fill(h.e, w);
    }
  const int NS = 40;
  auto hMap = new TH2D("hMap", "MC 90% bound after M1, x budget;L1 steel/Cu seam [m];L2 Cu/EFCu seam [m]", NS, 0, depthTot / 1000, NS, 0, depthTot / 1000);
  for (int i = 1; i <= NS; i++)
    for (int j = i; j <= NS; j++) // only L1 <= L2 is a design
    {
      double bi;
      hMap->SetBinContent(i, j, judge(hMap->GetXaxis()->GetBinCenter(i) * 1000, hMap->GetYaxis()->GetBinCenter(j) * 1000, bi) / cfg::rtBudget);
    }
  gStyle->SetOptStat(0);
  auto c = new TCanvas("c_bkg", "", 1400, 560);
  c->Divide(2, 1);
  c->cd(1)->SetLogy();
  gPad->SetGrid();
  gPad->SetLeftMargin(0.14);
  auto leg = new TLegend(0.6, 0.78, 0.88, 0.88);
  const int col[2] = {kAzure + 2, kOrange + 7};
  for (int k = 0; k < 2; k++)
  {
    sp[k]->SetLineColor(col[k]);
    sp[k]->SetLineWidth(2);
    sp[k]->Draw(k ? "HIST SAME" : "HIST");
    leg->AddEntry(sp[k], cutName[k], "l");
  }
  leg->Draw();
  c->cd(2)->SetLogz();
  gPad->SetRightMargin(0.15);
  hMap->SetMinimum(0.5 * hMap->GetMinimum(0));
  hMap->Draw("COLZ");
  auto star = new TMarker(L1ks / 1000, L2ks / 1000, 29);
  star->SetMarkerColor(kRed + 1);
  star->SetMarkerSize(2.2);
  star->Draw();
  std::string tag;
  for (auto &r : R) tag += (tag.empty() ? "" : "_") + r.iso;
  std::string png = rtOut(tag + "_background.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s\n", png.c_str());
}
