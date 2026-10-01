//  how much background the re-entrant tube puts into the germanium, and whether the choice and placement of its materials matter
//    root -l -b -q 'ana/sim_background.C("output/tl208.root=Tl208,output/bi214.root=Bi214")'
//    root -l -b -q 'ana/sim_background.C("output/tl208_*.root=Tl208,output/bi214_*.root=Bi214")'   job arrays, merged per isotope
//  each run is weighted by its nuclide's activity and its own measured sampling density, then the chains are summed
//  output/<run>_psd.csv from sim_psd.py supplies active energy and PSD. without it the study stops at M1 + argon and says so

#include "../geom/rt.h"
#include "TRegexp.h"
#include <sstream>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  const double qbb = 2039.0, roiHalf = 55.0;  // keV, the neutrinoless double beta window
  const double m1_keV = 5.0;                  // a detector counts as fired above this
  const double lar_keV = 20.0;                // argon veto: fires above this deposit. 20 keV is where the 4 PE cut sits
  const double psdCut = -1.80;                // AoE_class cut in Edgar's analysis. his aoe_class_paras.yaml says -0.83
  const double geMass_kg = 1000.0;            // the array being protected
  const double exposureYr = 10.0;             // live time to quote real-life event counts over
  const double decaysPerSec = 1600.0;         // this laptop, -t 8, remage v0.26.0: 1M decays in ~10 min. only turns a required N into a wall time
  const double bgGoal = 1e-5;                 // LEGEND-1000 target [cts/(keV kg yr)] in the ROI
  const double targetRel = 0.10;              // [7] sizes a run for this relative error on the BI AFTER cuts
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000;        // spectrum range [keV]
  const int    eBins = 200;
  const int    NDB = 10;                      // depth slabs for the reach profile

  const char  *mat[3] = {"steel", "Cu", "EFCu"};      // ordered top to bottom
  const double density[3] = {7900.0, 8960.0, 8930.0}; // kg/m^3
  struct Act { double v, dv; bool ul; };              // specific activity [uBq/kg], its uncertainty, only an upper limit?
  const Act act[3][2] = {{{2500, 0, false}, {1000, 0, false}},     // steel: Bi214, Tl208. Ralph's materialMix, no uncertainty quoted
                         {{1, 0, false}, {1, 0, false}},           // Cu: the same source
                         {{0.19, 0.10, false}, {0.077, 0, true}}}; // EFCu: Edgar's survival_BI.py radioassay
  int iso(const std::string &s) { return s.find("Bi") != std::string::npos ? 0 : 1; }
  double activity(const std::string &s, int m) { return act[m][iso(s)].v * 1e-6; } // -> Bq/kg
  double relUnc(const std::string &s, int m) { const Act &a = act[m][iso(s)]; return a.v > 0 ? a.dv / a.v : 0; }
  bool upperLimit(const std::string &s, int m) { return act[m][iso(s)].ul; }

  const double cdr[2][3] = {{21.0, 21.0, 19.0}, {1.2, 31.0, 16.0}};         // RE vessel (EFCu) survival in the ROI [%]: AC, PSD, PSD|AC. CDR
  const double edgar[2][3] = {{16.94, 18.12, 12.27}, {1.10, 30.32, 22.72}}; // the same from Edgar's remage (survival_BI.py)
}

//-------------------------------------------------------------------------------
//  2. One run, reduced to what a reweighting needs:
struct Hit // one detector's response to one decay
{
  Long64_t ev;     // its decay: the evtid while a file is read, then its slot among the run's decays with a hit
  float e;         // active energy if sim_psd.py ran, else the summed deposit
  bool m1, lar;    // passes M1, passes the argon veto
  signed char psd; // 1 passes PSD, 0 fails, -1 unknown (sim_psd.py not run)
};

struct Run
{
  std::string file, iso;
  Long64_t nsim = 0;
  int ndet = 0;
  bool hasPsd = false;                // detector response supplied by sim_psd.py
  long nPv[3] = {0, 0, 0};            // every decay counts, but only as numbers. per physical volume: 0 mother, 1 OFHC, 2 SS
  long nSlabPv[3][3] = {{0}};         // the same, per as-built section (steel, Cu, EFCu)
  double mcDensity[3] = {0, 0, 0};    // simulated decays per m^3 in each physical volume, measured
  long nBad = 0, nFiles = 1;          // decays outside the wall (a confinement bug); files merged in
  std::vector<float> depth;           // one slot per decay WITH a Ge hit, so memory follows the hits: mm below the top of the tube
  std::vector<signed char> pv;        // physical volume it was drawn in: 0 mother, 1 OFHC, 2 SS, -1 none
  std::vector<short> roiRaw, roiCut;  // ROI hits: before cuts, after the full chain
  std::vector<Long64_t> roiIdx;       // the slots with an ROI hit: the only decays a design can change
  std::vector<Hit> hits;
  long nHit[4] = {0, 0, 0, 0};        // no cuts, M1, M1 + argon, M1 + argon + PSD
  long dN[cfg::NDB] = {0}, dH[cfg::NDB] = {0}, dR[cfg::NDB] = {0}; // per depth slab: decays, hits, ROI hits
};

static int slabOf(double dep, double L1, double L2) { return dep < L1 ? 0 : (dep < L2 ? 1 : 2); } // steel, Cu, EFCu from the top
static int depthBin(double dep, double tot) { return std::min(cfg::NDB - 1, std::max(0, (int)(dep / tot * cfg::NDB))); }

static Run readRun(const std::string &file, const std::string &iso, const RT &rt)
{
  Run r;
  r.file = file;
  r.iso = iso;
  auto f = TFile::Open(file.c_str());
  auto d = (f && !f->IsZombie()) ? (TDirectory *)f->Get("stp") : nullptr;
  auto vtx = d ? (TTree *)d->Get("vtx") : nullptr;
  if (!vtx)
  {
    printf("ERROR: no stp/vtx in %s\n", file.c_str());
    return r;
  }
  const Long64_t n = r.nsim = vtx->GetEntries();
  std::vector<float> eLAr(n, 0), depth(n, -1); // per decay, but for this file only (10^7 on NERSC): the run keeps just the decays with a hit
  std::vector<short> nFired(n, 0), roiRaw(n, 0), roiCut(n, 0);
  std::vector<signed char> pv(n, -1);

  std::map<std::pair<int, std::string>, std::pair<float, float>> psd; // sim_psd.py: (event, detector) -> (active energy, AoE_class)
  std::string base = file.substr(file.find_last_of('/') + 1);
  std::ifstream in("output/" + base.substr(0, base.rfind(".root")) + "_psd.csv");
  for (std::string line; std::getline(in, line);)
  {
    std::string c[4];
    std::stringstream ss(line);
    for (auto &x : c) std::getline(ss, x, ',');
    if (isdigit(c[0][0])) psd[{std::stoi(c[0]), c[1]}] = {std::stof(c[2]), std::stof(c[3])}; // the header line starts with a letter
  }
  r.hasPsd = !psd.empty();

  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it()) // germanium: each detector's deposits summed per decay
  {
    auto t = rtIsGermanium(k->GetName()) ? (TTree *)d->Get(k->GetName()) : nullptr;
    if (!t || !t->GetBranch("edep_in_keV"))
      continue;
    r.ndet++;
    std::map<int, float> sum;
    rtScan(t, {"evtid", "edep_in_keV"}, [&](const double *v) { sum[(int)v[0]] += v[1]; });
    for (auto &p : sum)
    {
      if (p.first < 0 || p.first >= n || p.second <= 0)
        continue;
      float e = p.second;
      signed char pass = -1;
      if (r.hasPsd)
      {
        auto h = psd.find({p.first, k->GetName()});
        if (h == psd.end())
          continue; // sim_psd.py dropped it: every deposit sat in the dead layer
        e = h->second.first;
        pass = h->second.second > cfg::psdCut ? 1 : 0;
      }
      if (e <= cfg::m1_keV)
        continue;
      nFired[p.first]++; // a detector fired: M1 counts these
      r.hits.push_back({p.first, e, false, false, pass});
    }
  }
  if (auto tl = (TTree *)d->Get("undergroundlar"))
    rtScan(tl, {"evtid", "edep_in_keV"}, [&](const double *v) { if (v[0] >= 0 && v[0] < n) eLAr[(int)v[0]] += v[1]; });
  rtScan(vtx, {"evtid", "xloc_in_m", "yloc_in_m", "zloc_in_m"}, [&](const double *v) {
    int ev = (int)v[0]; // rows are NOT in event order with -t 8, so index by evtid or every decay gets someone else's vertex
    if (ev < 0 || ev >= n)
      return;
    double x = v[1] * 1000, y = v[2] * 1000, z = v[3] * 1000;
    depth[ev] = rt.zTop - z;
    pv[ev] = rt.physAt(sqrt(x * x + y * y), z); // the volume remage sampled it from
    if (pv[ev] < 0) r.nBad++;
  });

  const double tot = rt.zTop - rt.zBottom, L1 = rt.zTop - rt.seamSS, L2 = rt.zTop - rt.seamOFHC;
  for (auto &h : r.hits)
  {
    h.m1 = (nFired[h.ev] == 1);
    h.lar = (eLAr[h.ev] <= cfg::lar_keV);
    bool full = h.m1 && h.lar && (!r.hasPsd || h.psd == 1); // the whole chain this run can apply
    r.nHit[0]++;
    r.nHit[1] += h.m1;
    r.nHit[2] += h.m1 && h.lar;
    r.nHit[3] += r.hasPsd && full;
    if (fabs(h.e - cfg::qbb) <= cfg::roiHalf) { roiRaw[h.ev]++; roiCut[h.ev] += full; }
  }
  for (Long64_t i = 0; i < n; i++) // every decay, reduced to counts
  {
    r.dN[depthBin(depth[i], tot)]++;
    r.dR[depthBin(depth[i], tot)] += roiRaw[i];
    if (pv[i] >= 0) { r.nPv[pv[i]]++; r.nSlabPv[slabOf(depth[i], L1, L2)][pv[i]]++; }
  }
  std::vector<int> slot(n, -1); // a decay with a hit gets one slot, and each of its hits points there
  for (auto &h : r.hits)
  {
    r.dH[depthBin(depth[h.ev], tot)]++;
    int &sl = slot[h.ev];
    if (sl < 0)
    {
      sl = r.depth.size();
      r.depth.push_back(depth[h.ev]);
      r.pv.push_back(pv[h.ev]);
      r.roiRaw.push_back(roiRaw[h.ev]);
      r.roiCut.push_back(roiCut[h.ev]);
    }
    h.ev = sl;
  }
  delete f; // hundreds of job files: close each one
  return r;
}

//-------------------------------------------------------------------------------
//  3. Combine the runs:
static std::vector<std::string> expand(const std::string &pat) // "output/tl208_*.root" -> every match, sorted
{
  if (pat.find_first_of("*?") == std::string::npos) return {pat};
  size_t sl = pat.find_last_of('/');
  std::string dir = sl == std::string::npos ? "." : pat.substr(0, sl), base = pat.substr(sl + 1); // npos+1 == 0
  TRegexp re(base.c_str(), kTRUE); // kTRUE: shell wildcard, not a regex
  std::vector<std::string> out;
  void *d = gSystem->OpenDirectory(dir.c_str());
  for (const char *e; d && (e = gSystem->GetDirEntry(d));)
  {
    TString f(e);
    Ssiz_t len = 0;
    if (f.Index(re, &len) == 0 && len == f.Length()) // the whole name must match
      out.push_back(sl == std::string::npos ? f.Data() : dir + "/" + f.Data());
  }
  if (d) gSystem->FreeDirectory(d);
  std::sort(out.begin(), out.end());
  if (out.empty()) printf("no file matches %s\n", pat.c_str());
  return out;
}

void sim_background(const char *runs = "output/tl208.root=Tl208,output/bi214.root=Bi214", const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  if (!rt.ok)
    return;
  std::vector<Run> R; // one per isotope
  auto append = [](auto &a, const auto &b) { a.insert(a.end(), b.begin(), b.end()); };
  std::stringstream list(runs);
  for (std::string tok; std::getline(list, tok, ',');) // "a.root=Tl208,b.root=Bi214": split on commas, then each on '='
  {
    size_t eq = tok.find('=');
    std::string iso = eq == std::string::npos ? "Tl208" : tok.substr(eq + 1);
    for (auto &f : tok.empty() ? std::vector<std::string>() : expand(tok.substr(0, eq)))
    {
      Run r = readRun(f, iso, rt);
      if (!r.nsim)
        continue;
      Run *into = nullptr;
      for (auto &m : R)
        if (m.iso == iso) into = &m;
      if (!into) { R.push_back(std::move(r)); continue; }
      // files of one isotope are one run split into jobs: stitch them, slots renumbered, or it counts once per file
      Long64_t off = into->depth.size();
      for (auto h : r.hits) { h.ev += off; into->hits.push_back(h); }
      append(into->depth, r.depth);
      append(into->pv, r.pv);
      append(into->roiRaw, r.roiRaw);
      append(into->roiCut, r.roiCut);
      for (int k = 0; k < 4; k++) into->nHit[k] += r.nHit[k];
      for (int p = 0; p < 3; p++) { into->nPv[p] += r.nPv[p]; for (int k = 0; k < 3; k++) into->nSlabPv[k][p] += r.nSlabPv[k][p]; }
      for (int b = 0; b < cfg::NDB; b++) { into->dN[b] += r.dN[b]; into->dH[b] += r.dH[b]; into->dR[b] += r.dR[b]; }
      if (into->hasPsd != r.hasPsd)
        printf("WARNING: %s files disagree on sim_psd.py output - run it on every file\n", iso.c_str());
      into->hasPsd = into->hasPsd && r.hasPsd;
      into->nsim += r.nsim;
      into->nBad += r.nBad;
      into->nFiles++;
      into->file += "+" + r.file;
    }
  }
  if (R.empty())
  {
    printf("no usable runs in '%s'\n", runs);
    return;
  }

  const double depthTot = rt.zTop - rt.zBottom, L1ks = rt.zTop - rt.seamSS, L2ks = rt.zTop - rt.seamOFHC; // KS seams as depths
  const double binW = (cfg::eHi - cfg::eLo) / cfg::eBins, UL90 = 2.30, roiW = 2 * cfg::roiHalf;
  const int ORD[3] = {0, 1, 2};
  double physV[3]; // m^3 of the mother, the OFHC shell, the SS shell
  for (int p = 0; p < 3; p++) physV[p] = rt.physVolume(p) * 1e-9;
  for (auto &r : R)
  {
    for (int p = 0; p < 3; p++) r.mcDensity[p] = r.nPv[p] / physV[p]; // remage fills the mother ~20% sparser than its daughters: measure, never assume
    for (Long64_t i = 0; i < (Long64_t)r.roiRaw.size(); i++)
      if (r.roiRaw[i]) r.roiIdx.push_back(i); // after cuts is a subset of before
  }
  printf("geometry  : %s\n", rt.file.c_str());
  printf("RT wall   : %.5f m^3 over %.3f m\n", rt.wallVolume(rt.zBottom, rt.zTop) * 1e-9, depthTot / 1000);
  printf("cuts      : M1 (> %.0f keV) + argon veto (<= %.0f keV) + PSD (AoE_class > %.2f)   ROI %.0f +/- %.0f keV\n",
         cfg::m1_keV, cfg::lar_keV, cfg::psdCut, cfg::qbb, cfg::roiHalf);
  for (auto &r : R)
    printf("            %-6s %s\n", r.iso.c_str(), r.hasPsd ? "detector response from sim_psd.py: active energy + PSD"
                                                          : "NO sim_psd.py output: summed deposit, and the chain stops at M1 + argon");
  printf("\n");

  // one simulated decay = rho x A x 1 yr real decays per m^3 of material m, over the measured MC density of the
  // volume pv it was drawn in, per kg of Ge. m is what the design puts at that depth, so a design is arithmetic
  auto weight = [&](const Run &r, int m, int pv) {
    if (pv < 0 || r.mcDensity[pv] <= 0) return 0.0;
    return cfg::density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[pv] * cfg::geMass_kg);
  };
  auto rate = [&](double L1, double L2, const int *ord, bool cut, long *nmc) { // ROI rate over every run, for seams L1, L2 and a material order
    double sum = 0;
    long n = 0;
    for (auto &r : R)
      for (Long64_t i : r.roiIdx) // the seam scan calls this ~1800 times: never loop over every decay
      {
        int k = cut ? r.roiCut[i] : r.roiRaw[i];
        double w = weight(r, ord[slabOf(r.depth[i], L1, L2)], r.pv[i]) / roiW;
        sum += k * w;
        n += k;
      }
    if (nmc) *nmc = n;
    return sum;
  };
  auto xg = [&](double v) { return Form(v / cfg::bgGoal < 10 ? "%.2f" : "%.0f", v / cfg::bgGoal); }; // in units of the goal

  for (auto &r : R)
  {
//-------------------------------------------------------------------------------
//  4. Cut ladder, per run:
    printf("=== %s   (%s, %lld decays, %d Ge tables) ===\n", r.iso.c_str(), r.file.c_str(), r.nsim, r.ndet);
    printf("[4] hits surviving each cut: no cuts %ld -> M1 %ld -> +argon %ld", r.nHit[0], r.nHit[1], r.nHit[2]);
    if (r.hasPsd)
      printf(" -> +PSD %ld", r.nHit[3]);
    printf("   (argon veto x%.1f)\n", r.nHit[2] ? (double)r.nHit[1] / r.nHit[2] : 0.0);

//-------------------------------------------------------------------------------
//  5. Normalisation, per physical volume as built:
    const int asBuilt[3] = {2, 1, 0}; // mother is EFCu, the OFHC shell Cu, the SS shell steel
    const char *pvName[3] = {"mother", "OFHC", "SS"};
    double T = cfg::exposureYr * cfg::secPerYear;
    printf("[5] normalisation:  %-7s %-6s %9s %8s %10s %11s %9s %11s %10s\n", "volume", "mat", "V [m^3]", "M [kg]", "A [Bq]",
           Form("N in %.0f yr", cfg::exposureYr), "MC dec.", "MC per m^3", "per MC dec");
    for (int p = 0; p < 3; p++)
    {
      int m = asBuilt[p];
      double mass = physV[p] * cfg::density[m], act = mass * cfg::activity(r.iso, m);
      printf("                    %-7s %-6s %9.5f %8.1f %10.3e %11.3e %9ld %11.0f %10.2f\n", pvName[p], cfg::mat[m], physV[p],
             mass, act, T * act, r.nPv[p], r.mcDensity[p], r.nPv[p] ? T * act / r.nPv[p] : 0.0);
    }
    double dAvg = 0.5 * (r.mcDensity[1] + r.mcDensity[2]);
    printf("      sampling density mother / shells = %.3f   (1.000 would be uniform; weights use the measured value)\n",
           dAvg > 0 ? r.mcDensity[0] / dAvg : 0.0);
    if (r.nBad)
      printf("      WARNING: %ld decays outside the wall - a confinement bug\n", r.nBad);
    if (r.nFiles > 1)
      printf("      merged from %ld files\n", r.nFiles);

//-------------------------------------------------------------------------------
//  6. Reach vs depth:
    printf("[6] reach vs depth:  %-16s %9s %8s %13s %9s\n", "depth [m]", "decays", "hits", "hits/decay", "ROI");
    for (int b = 0; b < cfg::NDB; b++)
      printf("      %6.2f..%-9.2f %9ld %8ld %13.2e %9ld%s\n", depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000,
             r.dN[b], r.dH[b], r.dN[b] ? (double)r.dH[b] / r.dN[b] : 0.0, r.dR[b], b == cfg::NDB - 1 ? "   <- nearest the detectors" : "");

//-------------------------------------------------------------------------------
//  7. Sizing for the target after cuts:
    long roiB = 0, roiA = 0; // the BI after cuts needs 1/p^2 SURVIVING ROI hits
    for (Long64_t i : r.roiIdx) { roiB += r.roiRaw[i]; roiA += r.roiCut[i]; }
    int is = cfg::iso(r.iso);
    double surv = r.hasPsd ? cfg::edgar[is][0] * cfg::edgar[is][2] / 1e4 : cfg::edgar[is][0] / 100; // Edgar's Combined = AC x PSD|AC. no PSD: AC alone
    double eff = roiA >= 10 ? (double)roiA / r.nsim : (roiB ? (double)roiB : UL90) / r.nsim * surv;  // measured once 10 survive, else before cuts x Edgar
    bool lower = roiA < 10 && roiB == 0; // no ROI hit at all yet: the rate is a limit, N a lower bound
    double need = 1.0 / (cfg::targetRel * cfg::targetRel) / eff;
    printf("[7] sizing for %.0f%% after cuts: %s%.1e decays (%.0f h here)   ROI hits %ld before / %ld after cuts\n",
           100 * cfg::targetRel, lower ? ">" : "", need, need / cfg::decaysPerSec / 3600, roiB, roiA);
    if (roiA < 10)
      printf("      after-cut hits per decay %s%.1e = %s/%lld before cuts x Edgar's Combined %.2f%%, until 10 survive\n",
             lower ? "< " : "", eff, roiB ? Form("%ld", roiB) : "2.3 (90% CL)", r.nsim, 100 * surv);
    int dead = 0; // the slabs from the top down that no hit came from
    long nDead = 0;
    while (dead < cfg::NDB && r.dH[dead] == 0 && r.dN[dead] > 500) nDead += r.dN[dead++];
    double fr = (double)nDead / r.nsim;
    if (fr > 0.05)
      printf("      %.0f%% of decays landed above %.2f m depth and reached nothing\n", 100 * fr, depthTot * dead / cfg::NDB / 1000);
    printf("\n");
  }

//-------------------------------------------------------------------------------
//  8. Designs, summed over every run:
  struct { const char *name; double L1, L2; int ord[3]; } design[] = {
      {"KS as built: steel / Cu / EFCu", L1ks, L2ks, {0, 1, 2}},
      {"all EFCu", 0, 0, {2, 2, 2}},
      {"all Cu", 0, 0, {1, 1, 1}},
      {"all steel", 0, 0, {0, 0, 0}},
      {"KS, steel reaching twice as far down", 2 * L1ks, L2ks, {0, 1, 2}},
      {"KS, no steel (Cu down to L2)", 0, L2ks, {0, 1, 2}},
      {"SAME slabs, order reversed: EFCu / Cu / steel", L1ks, L2ks, {2, 1, 0}}};
  printf("[8] ROI rate [cts/(keV kg yr)], summed over %zu run(s), target %.0e:\n", R.size(), cfg::bgGoal);
  printf("      %-44s %8s %13s %8s %13s %8s\n", "design", "MC hits", "before cuts", "x goal", "after cuts", "x goal");
  for (auto &q : design)
  {
    long n0, n1;
    double r0 = rate(q.L1, q.L2, q.ord, false, &n0), r1 = rate(q.L1, q.L2, q.ord, true, &n1);
    printf("      %-44s %3ld/%-4ld %12.2e %8s", q.name, n1, n0, r0, xg(r0));
    if (n1 > 0)
      printf(" %12.2e %8s\n", r1, xg(r1));
    else
    {
      double ul = n0 > 0 ? r0 * UL90 / n0 : 0.0; // nothing survived: a 90% CL limit from the before-cut weight
      printf(" %12s %8s\n", Form("< %.1e", ul), Form("< %s", xg(ul)));
    }
  }

//-------------------------------------------------------------------------------
//  9. Survival in the ROI, per chain and section:
  auto surv = [&](long pass, long of) -> std::string { // percent, binomial error; a limit when nothing survives
    if (of <= 0) return "-";
    if (pass == 0) return Form("< %.1f", 100.0 * UL90 / of);
    double f = (double)pass / of;
    return Form("%.1f +- %.1f", 100 * f, 100 * sqrt(f * (1 - f) / of));
  };
  const char *secName[4] = {"steel", "Cu", "EFCu", "tube"};
  std::vector<std::array<long, 4>> tube; // per chain, whole tube: N0, AC, PSD, AC + PSD
  printf("\n[9] survival in the ROI, counted in hits as Edgar does [%%]:\n");
  printf("      %-6s %-6s %5s   %-16s %-16s %-16s %-16s\n", "chain", "sect", "N0", "AC (M1+argon)", "PSD", "PSD | AC", "Combined");
  for (auto &r : R)
  {
    long n[4][4] = {{0}}; // [section, 3 = whole tube][N0, AC, PSD, AC + PSD]
    for (auto &h : r.hits)
      if (fabs(h.e - cfg::qbb) <= cfg::roiHalf)
        for (int k : {slabOf(r.depth[h.ev], L1ks, L2ks), 3})
        {
          n[k][0]++;
          n[k][1] += h.m1 && h.lar;
          n[k][2] += h.psd == 1;
          n[k][3] += h.m1 && h.lar && h.psd == 1;
        }
    for (int k = 0; k < 4; k++)
    {
      printf("      %-6s %-6s %5ld   %-16s", k ? "" : r.iso.c_str(), secName[k], n[k][0], surv(n[k][1], n[k][0]).c_str());
      if (r.hasPsd) // PSD|AC is a step: of the hits AC kept, how many PSD keeps too
        printf(" %-16s %-16s %-16s\n", surv(n[k][2], n[k][0]).c_str(), surv(n[k][3], n[k][1]).c_str(), surv(n[k][3], n[k][0]).c_str());
      else
        printf(" (PSD columns need sim_psd.py)\n");
    }
    tube.push_back({n[3][0], n[3][1], n[3][2], n[3][3]});
  }

//-------------------------------------------------------------------------------
//  10. Background index, per chain and section:
  struct Sum { double w = 0, w2 = 0, wbar = 0; long n = 0; }; // one section: summed weight, its square for Poisson, mean weight of every decay there, ROI hits
  struct Tot { double v = 0, va = 0, vs = 0; long n = 0; bool ul = false; }; // a total: value, radioassay and MC variances, hits, rests on a limit?
  auto cell = [&](const Sum &q, double rel, bool ul) -> std::string {
    if (q.n == 0) return Form("< %.2e", UL90 * q.wbar); // nothing seen: 90% CL on the MC alone
    return Form("%s%.2e +-%.1e +-%.1e", ul ? "<" : "", q.w, ul ? 0.0 : q.w * rel, sqrt(q.w2));
  };
  auto add = [&](Tot &t, const Sum &q, double rel, bool ul) {
    if (q.n == 0) return; // saw nothing: only bounded (its own row shows the limit), so it stays out of the total
    t.v += q.w; t.va += pow(q.w * rel, 2); t.vs += q.w2; t.n += q.n; t.ul |= ul; // radioassay and MC are independent: variances add
  };
  auto show = [&](const Tot &t, bool goal) -> std::string {
    if (t.n == 0) return "no ROI hit";
    std::string out = Form("%s%.2e +-%.1e +-%.1e", t.ul ? "<" : "", t.v, sqrt(t.va), sqrt(t.vs));
    return goal ? out + Form("  (%s%s x goal)", t.ul ? "<" : "", xg(t.v)) : out;
  };
  printf("\n[10] background index, KS as built [cts/(keV kg yr)] = BI +- radioassay +- MC statistics:\n");
  printf("      %-6s %-6s %9s   %-34s %-34s\n", "chain", "sect", "ROI hits", "before cuts", "after cuts");
  Tot gB, gA;
  for (auto &r : R)
  {
    Tot tB, tA;
    for (int sec = 0; sec < 3; sec++)
    {
      double rel = cfg::relUnc(r.iso, sec);
      bool ul = cfg::upperLimit(r.iso, sec);
      Sum qb, qa; // before and after cuts
      long nsec = 0;
      for (int p = 0; p < 3; p++) // mean weight over every decay in the section, from the counts
      {
        nsec += r.nSlabPv[sec][p];
        qb.wbar += r.nSlabPv[sec][p] * weight(r, sec, p) / roiW;
      }
      qa.wbar = qb.wbar = nsec ? qb.wbar / nsec : 0;
      for (Long64_t i : r.roiIdx)
      {
        if (slabOf(r.depth[i], L1ks, L2ks) != sec)
          continue;
        double w = weight(r, sec, r.pv[i]) / roiW; // per decay: the slab can mix volumes
        if (r.roiRaw[i]) { qb.w += r.roiRaw[i] * w; qb.w2 += r.roiRaw[i] * w * w; qb.n += r.roiRaw[i]; }
        if (r.roiCut[i]) { qa.w += r.roiCut[i] * w; qa.w2 += r.roiCut[i] * w * w; qa.n += r.roiCut[i]; }
      }
      printf("      %-6s %-6s %4ld/%-4ld   %-34s %-34s%s\n", sec ? "" : r.iso.c_str(), secName[sec], qa.n, qb.n,
             cell(qb, rel, ul).c_str(), cell(qa, rel, ul).c_str(), ul ? "  activity is an upper limit" : "");
      for (Tot *t : {&tB, &gB}) add(*t, qb, rel, ul);
      for (Tot *t : {&tA, &gA}) add(*t, qa, rel, ul);
    }
    printf("      %-6s %-6s %9s   %-34s %-34s\n", "", "tube", "", show(tB, false).c_str(), show(tA, false).c_str());
  }
  printf("      %-13s %9s   %-34s %-34s\n", "ALL CHAINS", "", show(gB, true).c_str(), show(gA, true).c_str());

//-------------------------------------------------------------------------------
//  11. Draw: spectrum, seam scan, survival
  const int NS = 60;
  auto hScan = new TH2D("hScan", "ROI rate before cuts, all chains;L1 steel/Cu boundary [m];L2 Cu/EFCu boundary [m]",
                        NS, 0, depthTot / 1000, NS, 0, depthTot / 1000);
  for (int i = 1; i <= NS; i++)
    for (int j = i; j <= NS; j++) // only L1 <= L2 is a design
      hScan->SetBinContent(i, j, rate(hScan->GetXaxis()->GetBinCenter(i) * 1000, hScan->GetYaxis()->GetBinCenter(j) * 1000, ORD, false, nullptr));
  hScan->SetMinimum(0.5 * hScan->GetMinimum(0)); // log z would otherwise cut at max/1000 and blank every design with EFCu at the bottom

  const char *cutLabel[4] = {"no cuts", "M1", "M1 + argon veto", "M1 + argon veto + PSD"};
  const int cutColor[4] = {kAzure + 2, kOrange + 7, kRed + 1, kViolet + 1};
  TH1D *sp[4]; // spectra weighted to cts/(keV kg yr) with the as-built design, summed over runs
  for (int k = 0; k < 4; k++)
    sp[k] = new TH1D(Form("sp%d", k), "RT background in the germanium;energy [keV];cts / (keV kg yr)", cfg::eBins, cfg::eLo, cfg::eHi);
  bool anyPsd = false;
  for (auto &r : R)
  {
    anyPsd |= r.hasPsd;
    for (auto &h : r.hits)
    {
      double w = weight(r, ORD[slabOf(r.depth[h.ev], L1ks, L2ks)], r.pv[h.ev]) / binW;
      bool pass[4] = {true, h.m1, h.m1 && h.lar, h.m1 && h.lar && h.psd == 1};
      for (int k = 0; k < 4; k++)
        if (pass[k]) sp[k]->Fill(h.e, w);
    }
  }
  auto dot = [](auto *o, int col, int style, double size) { o->SetMarkerColor(col); o->SetMarkerStyle(style); o->SetMarkerSize(size); };
  gStyle->SetOptStat(0);
  auto c = new TCanvas("c_bkg", "", 1280, 560);
  c->Divide(2, 1);
  c->cd(1)->SetLogy();
  gPad->SetGrid();
  gPad->SetLeftMargin(0.15);
  sp[0]->GetYaxis()->SetTitleOffset(1.6);
  auto leg = new TLegend(0.50, 0.72, 0.88, 0.88);
  leg->SetTextSize(0.032);
  for (int k = 0; k < (anyPsd ? 4 : 3); k++)
  {
    sp[k]->SetLineColor(cutColor[k]);
    sp[k]->SetLineWidth(2);
    sp[k]->Draw(k ? "HIST SAME" : "HIST");
    leg->AddEntry(sp[k], cutLabel[k], "l");
  }
  leg->Draw();
  c->cd(2)->SetLogz();
  gPad->SetRightMargin(0.17);
  hScan->Draw("COLZ");
  auto star = new TMarker(L1ks / 1000, L2ks / 1000, 29); // the KS design
  dot(star, kRed + 1, 29, 2.2);
  star->Draw();
  std::string tag;
  for (auto &r : R) tag += (tag.empty() ? "" : "_") + r.iso;
  std::string png = rtOut(tag + "_background.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s\n", png.c_str());

  auto cs = new TCanvas("c_surv", "", 560 * R.size(), 480); // survival in the ROI against the RE-vessel references
  cs->Divide(R.size(), 1);
  const char *cutName[4] = {"AC (M1+argon)", "PSD", "PSD | AC", "Combined"};
  for (size_t i = 0; i < R.size(); i++)
  {
    cs->cd(i + 1)->SetLogy();
    gPad->SetGrid(0, 1);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    int is = cfg::iso(R[i].iso);
    auto bars = new TH1D(Form("cdr%zu", i), Form("%s from the RT wall;;survival in ROI [%%]", R[i].iso.c_str()), 4, 0, 4); // the CDR
    auto ge = new TGraph();       // Edgar's remage for the same component
    auto gm = new TGraphErrors(); // this simulation: whole tube, binomial errors, each cut against its own denominator
    std::vector<TArrow *> lim;    // nothing survived: a 90% CL limit, drawn as an arrow down from it
    const long *n = tube[i].data();
    long pass[4] = {n[1], n[2], n[3], n[3]}, of[4] = {n[0], n[0], n[1], n[0]};
    for (int k = 0; k < 4; k++)
    {
      bars->SetBinContent(k + 1, k < 3 ? cfg::cdr[is][k] : cfg::cdr[is][0] * cfg::cdr[is][2] / 100); // combined = AC x PSD|AC
      bars->GetXaxis()->SetBinLabel(k + 1, cutName[k]);
      ge->SetPoint(k, k + 0.5, k < 3 ? cfg::edgar[is][k] : cfg::edgar[is][0] * cfg::edgar[is][2] / 100);
      if ((k > 0 && !R[i].hasPsd) || of[k] <= 0) // PSD points need sim_psd.py
        continue;
      double fr = (double)pass[k] / of[k];
      if (pass[k] == 0) // a log axis cannot show zero. above 100% a limit says nothing
        lim.push_back(new TArrow(k + 0.5, std::min(100.0, 100.0 * UL90 / of[k]), k + 0.5, std::min(100.0, 100.0 * UL90 / of[k]) / 4, 0.02, "|>"));
      else
      {
        gm->SetPoint(gm->GetN(), k + 0.5, 100 * fr);
        gm->SetPointError(gm->GetN() - 1, 0, 100 * sqrt(fr * (1 - fr) / of[k]));
      }
    }
    bars->SetFillColor(kAzure - 9);
    bars->SetLineColor(kAzure + 2);
    bars->SetBarWidth(0.5);
    bars->SetBarOffset(0.25);
    bars->SetMinimum(0.05);
    bars->SetMaximum(3000); // headroom above 100% for the legend
    bars->GetXaxis()->SetLabelSize(0.05);
    bars->Draw("BAR");
    dot(ge, kGray + 2, 24, 1.6);
    ge->Draw("P SAME");
    dot(gm, kRed + 1, 20, 1.4);
    gm->SetLineColor(kRed + 1);
    if (gm->GetN()) gm->Draw("P SAME");
    for (auto a : lim)
    {
      a->SetLineColor(kRed + 1);
      a->SetFillColor(kRed + 1);
      a->SetLineWidth(2);
      a->Draw();
    }
    auto lg = new TLegend(0.36, 0.72, 0.89, 0.89);
    lg->SetTextSize(0.035);
    lg->AddEntry(bars, "CDR - RE vessel (EFCu)", "f");
    lg->AddEntry(ge, "Edgar remage - RE Cu", "p");
    lg->AddEntry(gm, Form("this sim - KS tube (%ld ROI hits)", n[0]), "pe");
    if (!lim.empty()) lg->AddEntry(lim[0], "this sim - none survived (90% CL)", "l");
    lg->Draw();
  }
  std::string spng = rtOut(tag + "_survival.png");
  cs->SaveAs(spng.c_str());
  printf("wrote %s\n", spng.c_str());
}
