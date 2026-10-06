//  the cheapest re-entrant tube: always steel (top) / Cu / EFCu (bottom), with the steel/Cu seam L1 as deep and the
//  Cu/EFCu seam L2 as deep as the tube's background budget allows. one run per chain, reweighted to every (L1, L2)
//    root -l -b -q 'ana/background.C("output/tl208.root=Tl208,output/bi214.root=Bi214")'
//    root -l -b -q 'ana/background.C("output/tl208_*.root=Tl208,output/bi214_*.root=Bi214")'   job arrays, merged per isotope
//  each run is weighted by its nuclide's activity and its own measured sampling density, then the chains are summed
//  output/<run>_response.csv from response.py supplies active energy and PSD. without it the study stops at M1 + argon and says so

#include "../geom/rt.h"
#include "TRegexp.h"
#include <sstream>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  // the background window of MAJORANA [Haufe et al., arXiv:2209.10592]: 1950-2350 keV, minus 10 keV around Qbb and the
  // 2103.5 (Tl208 single escape), 2118.5 and 2204.1 keV (Bi214) lines [Arnquist et al., arXiv:2207.07638]. 360 keV of flat
  // continuum, 3.3x the 110 keV ROI around Qbb it replaced: the BI is its counts / 360 keV
  const double winLo = 1950, winHi = 2350, winGap[4] = {2039.0, 2103.5, 2118.5, 2204.1}, winGapHalf = 5.0;
  bool inWindow(double e) { if (e < winLo || e > winHi) return false; for (double g : winGap) if (fabs(e - g) <= winGapHalf) return false; return true; }
  const double winWidth = winHi - winLo - 4 * 2 * winGapHalf; // 360 keV
  const double reachEnd = 4400;               // mm depth. above it a decay's chance of a Ge hit is fitted as exp(-distance / lambda) and
                                              // extrapolated, so steel can be judged where the MC saw nothing; below it the MC is used as it is
  const double line[2] = {1764.5, 2614.5};    // keV: the strongest line of each chain (Bi214, Tl208), the tracer MAJORANA used per detector
  const double m1_keV = 5.0;                  // a detector counts as fired above this
  const double lar_keV = 20.0;                // argon veto: fires above this deposit. 20 keV is where the 4 PE cut sits
  const double psdCut = -1.80;                // AoE_class cut in Edgar's analysis. his aoe_class_paras.yaml says -0.83
  const double geMass_kg = 1000.0;            // the array being protected
  const double exposureYr = 10.0;             // live time to quote real-life event counts over
  const double decaysPerSec = 1600.0;         // this laptop, -t 8, remage v0.26.0: 1M decays in ~10 min. only turns a required N into a wall time
  const double bgGoal = 1e-5;                 // LEGEND-1000 target [cts/(keV kg yr)] at Qbb
  const double rtBudget = 1e-5;               // the tube's share of that goal: a design passes when its projected 90% BI is below this
  const double targetRel = 0.10;              // [7] sizes a run for this relative error on the BI AFTER cuts
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000;        // spectrum range [keV]
  const int    eBins = 200;
  const int    NDB = 10;                      // depth slabs for the reach profile
  const int    NZ = 625;                      // ~1 cm depth bins: decays per section, the reach fit, the projection

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
  short det;       // index into Run::dets
  float e;         // active energy if response.py ran, else the summed deposit
  bool m1, lar;    // passes M1, passes the argon veto
  signed char psd; // 1 passes PSD, 0 fails, -1 unknown (response.py not run)
};

struct Run
{
  std::string file, iso;
  Long64_t nsim = 0;
  int ndet = 0;
  std::vector<std::string> dets;      // Ge detector names, V<string><position>
  bool hasPsd = false;                // detector response supplied by response.py
  long nPv[3] = {0, 0, 0};            // every decay counts, but only as numbers. per physical volume: 0 mother, 1 OFHC, 2 SS
  long nRowPv[4][3] = {{0}};          // the same, per as-built row (steel, Cu, EFCu wall, EFCu head)
  long nZPv[3][cfg::NZ] = {{0}};      // the same, per ~1 cm of depth
  double mcDensity[3] = {0, 0, 0};    // simulated decays per m^3 in each physical volume, measured
  double lam = 0, P0 = 0, lamHi = 0, P0Hi = 0; // the fitted reach: P(hit) = P0 exp((depth - reachEnd) / lam), and its 90% far-reaching variant
  double fWin = 0, fWinHi = 0, fCut = 0, fCutHi = 0; // window hits per decay with a hit, before / after cuts, and 90% upper values
  long nBad = 0, nFiles = 1;          // decays outside the wall (a confinement bug); files merged in
  std::vector<float> depth;           // one slot per decay WITH a Ge hit, so memory follows the hits: mm below the top of the tube
  std::vector<signed char> pv;        // physical volume it was drawn in: 0 mother, 1 OFHC, 2 SS, -1 none
  std::vector<short> winRaw, winCut;  // hits in the background window: before cuts, after the full chain
  std::vector<Long64_t> winIdx;       // the slots with a window hit: the only decays a design can change
  std::vector<Hit> hits;
  long nHit[4] = {0, 0, 0, 0};        // no cuts, M1, M1 + argon, M1 + argon + PSD
  long dN[cfg::NDB] = {0}, dH[cfg::NDB] = {0}, dR[cfg::NDB] = {0}; // per depth slab: decays, hits, window hits
};

static int slabOf(double dep, double L1, double L2) { return dep < L1 ? 0 : (dep < L2 ? 1 : 2); } // a design's three slabs, from the top
static int rowOf(double dep, double L1, double L2, double Lh) { return dep < Lh ? slabOf(dep, L1, L2) : 3; } // as built, with the bottom head split off the EFCu
const int   rowMat[4] = {0, 1, 2, 2};                                  // the material of each row
const char *rowName[5] = {"steel", "Cu", "EFCu-w", "EFCu-h", "tube"}; // EFCu-w: the EFCu wall above the head. EFCu-h: the bottom head
static int depthBin(double dep, double tot) { return std::min(cfg::NDB - 1, std::max(0, (int)(dep / tot * cfg::NDB))); }
static int zBin(double dep, double tot) { return std::min(cfg::NZ - 1, std::max(0, (int)(dep / tot * cfg::NZ))); }

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
  std::vector<short> nFired(n, 0), winRaw(n, 0), winCut(n, 0);
  std::vector<signed char> pv(n, -1);

  std::map<std::pair<int, std::string>, std::pair<float, float>> psd; // response.py: (event, detector) -> (active energy, AoE_class)
  std::string base = file.substr(file.find_last_of('/') + 1);
  std::ifstream in("output/" + base.substr(0, base.rfind(".root")) + "_response.csv");
  r.hasPsd = (bool)in; // the file is there: response.py ran, even if a far section gave it no Ge hit to write
  for (std::string line; std::getline(in, line);)
  {
    std::string c[4];
    std::stringstream ss(line);
    for (auto &x : c) std::getline(ss, x, ',');
    if (isdigit(c[0][0])) psd[{std::stoi(c[0]), c[1]}] = {std::stof(c[2]), std::stof(c[3])}; // the header line starts with a letter
  }

  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it()) // germanium: each detector's deposits summed per decay
  {
    auto t = rtIsGermanium(k->GetName()) ? (TTree *)d->Get(k->GetName()) : nullptr;
    if (!t || !t->GetBranch("edep_in_keV"))
      continue;
    r.ndet++;
    r.dets.push_back(k->GetName());
    short id = r.dets.size() - 1;
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
          continue; // response.py dropped it: every deposit sat in the dead layer
        e = h->second.first;
        pass = h->second.second > cfg::psdCut ? 1 : 0;
      }
      if (e <= cfg::m1_keV)
        continue;
      nFired[p.first]++; // a detector fired: M1 counts these
      r.hits.push_back({p.first, id, e, false, false, pass});
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

  const double tot = rt.zTop - rt.zBottom, L1 = rt.zTop - rt.seamSS, L2 = rt.zTop - rt.seamOFHC, Lh = rt.zTop - rt.zHead;
  for (auto &h : r.hits)
  {
    h.m1 = (nFired[h.ev] == 1);
    h.lar = (eLAr[h.ev] <= cfg::lar_keV);
    bool full = h.m1 && h.lar && (!r.hasPsd || h.psd == 1); // the whole chain this run can apply
    r.nHit[0]++;
    r.nHit[1] += h.m1;
    r.nHit[2] += h.m1 && h.lar;
    r.nHit[3] += r.hasPsd && full;
    if (cfg::inWindow(h.e)) { winRaw[h.ev]++; winCut[h.ev] += full; }
  }
  for (Long64_t i = 0; i < n; i++) // every decay, reduced to counts
  {
    r.dN[depthBin(depth[i], tot)]++;
    r.dR[depthBin(depth[i], tot)] += winRaw[i];
    if (pv[i] >= 0) { r.nPv[pv[i]]++; r.nRowPv[rowOf(depth[i], L1, L2, Lh)][pv[i]]++; r.nZPv[pv[i]][zBin(depth[i], tot)]++; }
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
      r.winRaw.push_back(winRaw[h.ev]);
      r.winCut.push_back(winCut[h.ev]);
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

void background(const char *runs = "output/tl208.root=Tl208,output/bi214.root=Bi214", const char *gdml = "")
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
    if (iso.find("Tl") == std::string::npos && iso.find("Bi") == std::string::npos)
    {
      printf("ERROR: unknown isotope '%s' - cfg::act only has Tl208 and Bi214\n", iso.c_str());
      continue;
    }
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
      std::vector<short> map(r.dets.size()); // detector indices by name, in case a file lists them in another order
      for (size_t k = 0; k < r.dets.size(); k++)
      {
        auto it = std::find(into->dets.begin(), into->dets.end(), r.dets[k]);
        map[k] = it != into->dets.end() ? it - into->dets.begin() : (into->dets.push_back(r.dets[k]), into->dets.size() - 1);
      }
      for (auto h : r.hits) { h.ev += off; h.det = map[h.det]; into->hits.push_back(h); }
      append(into->depth, r.depth);
      append(into->pv, r.pv);
      append(into->winRaw, r.winRaw);
      append(into->winCut, r.winCut);
      for (int k = 0; k < 4; k++) into->nHit[k] += r.nHit[k];
      for (int p = 0; p < 3; p++) { into->nPv[p] += r.nPv[p]; for (int k = 0; k < 4; k++) into->nRowPv[k][p] += r.nRowPv[k][p]; for (int b = 0; b < cfg::NZ; b++) into->nZPv[p][b] += r.nZPv[p][b]; }
      for (int b = 0; b < cfg::NDB; b++) { into->dN[b] += r.dN[b]; into->dH[b] += r.dH[b]; into->dR[b] += r.dR[b]; }
      if (into->hasPsd != r.hasPsd)
        printf("WARNING: %s files disagree on response.py output - run it on every file\n", iso.c_str());
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

  const double depthTot = rt.zTop - rt.zBottom, L1ks = rt.zTop - rt.seamSS, L2ks = rt.zTop - rt.seamOFHC, Lhks = rt.zTop - rt.zHead; // KS seams and head as depths
  const double binW = (cfg::eHi - cfg::eLo) / cfg::eBins, UL90 = 2.30, winW = cfg::winWidth; // the BI is window counts per keV of window
  const int ORD[3] = {0, 1, 2};
  double physV[3]; // m^3 of the mother, the OFHC shell, the SS shell
  for (int p = 0; p < 3; p++) physV[p] = rt.physVolume(p) * 1e-9;
  for (auto &r : R)
  {
    for (int p = 0; p < 3; p++) r.mcDensity[p] = r.nPv[p] / physV[p]; // remage fills the mother ~20% sparser than its daughters: measure, never assume
    for (Long64_t i = 0; i < (Long64_t)r.winRaw.size(); i++)
      if (r.winRaw[i]) r.winIdx.push_back(i); // after cuts is a subset of before
  }
  printf("geometry  : %s\n", rt.file.c_str());
  printf("RT wall   : %.5f m^3 over %.3f m\n", rt.wallVolume(rt.zBottom, rt.zTop) * 1e-9, depthTot / 1000);
  printf("cuts      : M1 (> %.0f keV) + argon veto (<= %.0f keV) + PSD (AoE_class > %.2f)\n", cfg::m1_keV, cfg::lar_keV, cfg::psdCut);
  printf("window    : %.0f-%.0f keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's BEW): %.0f keV\n",
         cfg::winLo, cfg::winHi, cfg::winWidth);
  for (auto &r : R)
    printf("            %-6s %s\n", r.iso.c_str(), r.hasPsd ? "detector response from response.py: active energy + PSD"
                                                          : "NO response.py output: summed deposit, and the chain stops at M1 + argon");
  printf("\n");

  // one simulated decay = rho x A x 1 yr real decays per m^3 of material m, over the measured MC density of the
  // volume pv it was drawn in, per kg of Ge. m is what the design puts at that depth, so a design is arithmetic
  auto weight = [&](const Run &r, int m, int pv) {
    if (pv < 0 || r.mcDensity[pv] <= 0) return 0.0;
    return cfg::density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[pv] * cfg::geMass_kg);
  };
  // a design, honest about what the MC did not see: the measured window rate, and a 90% CL bound summed slab by slab.
  // a slab with n window hits is bounded at UL(n) x its mean hit weight; one with none at 2.30 x the heaviest decay in it
  auto ulN = [](long n) { const double t[11] = {2.30, 3.89, 5.32, 6.68, 7.99, 9.27, 10.53, 11.77, 12.99, 14.21, 15.41}; return n <= 10 ? t[n] : n + 1.28 * sqrt((double)n) + 1; };
  struct Eval { double bi = 0, up = 0; long n = 0; };
  auto judge = [&](double L1, double L2, const int *ord, bool cut) {
    Eval e;
    const double edge[4] = {0, L1, L2, depthTot};
    for (int s = 0; s < 3; s++)
    {
      double B = 0, wmax = 0;
      long n = 0;
      int b0 = (int)(edge[s] / depthTot * cfg::NZ), b1 = std::min(cfg::NZ, (int)ceil(edge[s + 1] / depthTot * cfg::NZ));
      for (auto &r : R)
      {
        for (int p = 0; p < 3; p++)
          for (int b = b0; b < b1; b++)
            if (r.nZPv[p][b]) { wmax = std::max(wmax, weight(r, ord[s], p) / winW); break; } // decays from section p lie in this slab
        for (Long64_t i : r.winIdx) // never loop over every decay: the scan calls this ~1800 times
          if (slabOf(r.depth[i], L1, L2) == s)
          {
            int k = cut ? r.winCut[i] : r.winRaw[i];
            B += k * weight(r, ord[s], r.pv[i]) / winW;
            n += k;
          }
      }
      e.bi += B;
      e.n += n;
      e.up += n ? ulN(n) * B / n : UL90 * wmax; // an empty slab has wmax 0
    }
    return e;
  };
  // the reach: how likely a decay at depth d is to put a hit (any energy) in the Ge. above reachEnd it falls through the
  // argon as exp(-distance / lambda); a Poisson fit of that tail, profiled over its normalisation, judges designs where the
  // MC saw few or no hits, i.e. where steel would go. lamHi is the one-sided 90% value (dlnL = 0.82): the farther reach
  const double dz = depthTot / cfg::NZ;
  const int bEnd = (int)(cfg::reachEnd / dz);
  std::vector<double> Vz(cfg::NZ); // wall volume per ~1 cm of depth [m^3]
  for (int b = 0; b < cfg::NZ; b++) Vz[b] = rt.wallVolume(rt.zTop - (b + 1) * dz, rt.zTop - b * dz) * 1e-9;
  auto survOf = [&](const Run &r) { int is = cfg::iso(r.iso); return r.hasPsd ? cfg::edgar[is][0] * cfg::edgar[is][2] / 1e4 : cfg::edgar[is][0] / 100; }; // Edgar's Combined
  for (auto &r : R)
  {
    const int G = 10, K = bEnd / G; // 10 cm fit bins, whole ones only
    std::vector<double> N(K, 0), n(K, 0), d(K);
    for (int k = 0; k < K; k++)
    {
      d[k] = (k + 0.5) * G * dz;
      for (int b = k * G; b < (k + 1) * G; b++) for (int p = 0; p < 3; p++) N[k] += r.nZPv[p][b];
    }
    long nWin = 0, nWinCut = 0;
    for (size_t i = 0; i < r.depth.size(); i++)
    {
      int k = zBin(r.depth[i], depthTot) / G;
      if (k < K) n[k]++;
      nWin += r.winRaw[i];
      nWinCut += r.winCut[i];
    }
    auto lnL = [&](double lam, double &P0) {
      double sn = 0, se = 0, l = 0;
      for (int k = 0; k < K; k++) { sn += n[k]; se += N[k] * exp((d[k] - cfg::reachEnd) / lam); }
      P0 = se > 0 ? sn / se : 0; // the normalisation that maximises lnL for this lambda
      for (int k = 0; k < K; k++)
      {
        double mu = N[k] * P0 * exp((d[k] - cfg::reachEnd) / lam);
        if (mu > 0) l += n[k] * log(mu) - mu;
      }
      return l;
    };
    double best = -1e300, P;
    for (double lam = 30; lam <= 1500; lam += 1) { double l = lnL(lam, P); if (l > best) { best = l; r.lam = lam; r.P0 = P; } }
    for (double lam = r.lam; lam <= 1500 && lnL(lam, P) >= best - 0.82; lam += 1) { r.lamHi = lam; r.P0Hi = P; }
    long nh = r.depth.size();
    r.fWin = nh ? (double)nWin / nh : 0;
    r.fWinHi = nh ? ulN(nWin) / nh : 0;
    r.fCut = nWinCut >= 10 ? (double)nWinCut / nh : r.fWin * survOf(r); // measured once 10 survive, until then Edgar's Combined survival
    r.fCutHi = nWinCut >= 10 ? ulN(nWinCut) / nh : r.fWinHi * survOf(r);
  }
  // a design's projected BI: the MC where it has statistics (below reachEnd), the fitted reach above it. hi = its 90% value
  auto project = [&](double L1, double L2, const int *ord, bool cut, bool hi) {
    double bi = 0, var = 0;
    for (auto &r : R)
    {
      for (Long64_t i : r.winIdx)
        if (r.depth[i] >= cfg::reachEnd)
        {
          int k = cut ? r.winCut[i] : r.winRaw[i];
          double w = weight(r, ord[slabOf(r.depth[i], L1, L2)], r.pv[i]) / winW;
          bi += k * w;
          var += k * w * w;
        }
      double lam = hi ? r.lamHi : r.lam, P0 = hi ? r.P0Hi : r.P0, f = cut ? (hi ? r.fCutHi : r.fCut) : (hi ? r.fWinHi : r.fWin);
      for (int b = 0; b < bEnd; b++) // decays per year x the reach x the window fraction
      {
        double dd = (b + 0.5) * dz;
        int m = ord[slabOf(dd, L1, L2)];
        bi += cfg::density[m] * cfg::activity(r.iso, m) * cfg::secPerYear * Vz[b] * f * P0 * exp((dd - cfg::reachEnd) / lam) / (cfg::geMass_kg * winW);
      }
    }
    return hi ? bi + 1.28 * sqrt(var) : bi;
  };
  auto massOf = [&](double d0, double d1, int m) { return d1 > d0 ? cfg::density[m] * rt.wallVolume(rt.zTop - d1, rt.zTop - d0) * 1e-9 : 0.0; };
  auto passes = [&](double L1, double L2) { return project(L1, L2, ORD, false, true) <= cfg::rtBudget; };
  auto deepestSteel = [&](double L2) { // the deepest steel/Cu seam that still passes for this L2, on a 1 cm grid. -1: not even 1 cm of steel
    double l1 = -1;
    for (double L1 = 10; L1 <= L2 - 10 && passes(L1, L2); L1 += 10) l1 = L1; // always some Cu between the seams: three materials
    return l1;
  };
  auto xg = [&](double v) { double x = v / cfg::bgGoal; return Form(x < 0.01 ? "%.0e" : (x < 10 ? "%.2f" : "%.0f"), x); };     // in units of the goal
  auto xb = [&](double v) { double x = v / cfg::rtBudget; return Form(x < 0.01 ? "%.0e" : (x < 10 ? "%.2f" : "%.0f"), x); }; // in units of the tube's budget

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
    printf("[6] reach vs depth:  %-16s %9s %8s %13s %9s\n", "depth [m]", "decays", "hits", "hits/decay", "window");
    for (int b = 0; b < cfg::NDB; b++)
      printf("      %6.2f..%-9.2f %9ld %8ld %13.2e %9ld%s\n", depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000,
             r.dN[b], r.dH[b], r.dN[b] ? (double)r.dH[b] / r.dN[b] : 0.0, r.dR[b], b == cfg::NDB - 1 ? "   <- nearest the detectors" : "");
    printf("      reach above %.1f m: P(hit) = %.2e x exp((depth - %.1f m) / lambda), lambda = %.0f mm (90%%: %.0f mm)\n",
           cfg::reachEnd / 1000, r.P0, cfg::reachEnd / 1000, r.lam, r.lamHi);
    printf("      window hits per decay with a hit: %.2e before cuts (90%%: %.2e), %.2e after%s\n", r.fWin, r.fWinHi, r.fCut,
           r.fCut == r.fWin * survOf(r) ? " (before x Edgar's Combined)" : "");

//-------------------------------------------------------------------------------
//  7. Statistics for the seams:
    // a slab of material m over section p that saw no window hit is only bounded, at 2.30 x one decay's weight. that bound
    // falls as 1/decays: this is how many decays each section needs to bound such a slab at 10% of the budget
    printf("[7] statistics: decays each section needs before the MC alone bounds a slab over it, with no window hit, at 10%% of the budget\n");
    printf("      %-8s %9s   %-26s %-26s\n", "section", "decays", "steel there: bound, needs", "Cu there: bound, needs");
    for (int p = 2; p >= 0; p--) // top to bottom
    {
      std::string c[2];
      for (int m = 0; m < 2; m++)
      {
        double b = UL90 * weight(r, m, p) / winW;
        c[m] = r.nPv[p] ? Form("%.1e, %.1e", b, r.nPv[p] * b / (0.1 * cfg::rtBudget)) : "-";
      }
      printf("      %-8s %9ld   %-26s %-26s\n", pvName[p], r.nPv[p], c[0].c_str(), c[1].c_str());
    }
    long winB = 0, winA = 0; // the BI after cuts needs 1/p^2 SURVIVING window hits
    for (Long64_t i : r.winIdx) { winB += r.winRaw[i]; winA += r.winCut[i]; }
    double surv = survOf(r); // Edgar's Combined = AC x PSD|AC. no PSD: AC alone
    double eff = winA >= 10 ? (double)winA / r.nsim : (winB ? (double)winB : UL90) / r.nsim * surv;  // measured once 10 survive, else before cuts x Edgar
    double need = 1.0 / (cfg::targetRel * cfg::targetRel) / eff;
    printf("      %.0f%% on the as-built BI after cuts: %s%.1e decays (%.0f h here); window hits %ld before / %ld after cuts%s\n",
           100 * cfg::targetRel, winA < 10 && winB == 0 ? ">" : "", need, need / cfg::decaysPerSec / 3600, winB, winA,
           winA < 10 ? Form(", rate from before cuts x Edgar's Combined %.2f%%", 100 * surv) : "");
    int dead = 0; // the slabs from the top down that no hit came from
    long nDead = 0;
    while (dead < cfg::NDB && r.dH[dead] == 0 && r.dN[dead] > 500) nDead += r.dN[dead++];
    if (nDead)
      printf("      no hit at all came from above %.2f m depth (%.0f%% of the decays): there the fitted reach stands in\n",
             depthTot * dead / cfg::NDB / 1000, 100.0 * nDead / r.nsim);
    printf("\n");
  }

//-------------------------------------------------------------------------------
//  8. Three-material designs, summed over every run:
  // the cheapest tube in this order: least EFCu first (the deepest L2 for which some steel still passes), then most steel
  double L2best = 0, L2ksSteel = 0;
  for (double L2 = 20; L2 <= depthTot && passes(10, L2); L2 += 10) L2best = L2;       // even 1 cm of steel must pass
  for (double L2 = L1ks; L2 <= depthTot && passes(L1ks, L2); L2 += 10) L2ksSteel = L2;  // keeping KS's steel
  double L1best = deepestSteel(L2best), L1ksEf = deepestSteel(L2ks);
  struct { const char *name; double L1, L2; } design[] = {
      {"KS as built", L1ks, L2ks},
      {"KS EFCu, the most steel", L1ksEf, L2ks},
      {"KS steel, the least EFCu", L1ks, L2ksSteel},
      {"the least EFCu, then the most steel", L1best, L2best}};
  printf("[8] three-material designs [cts/(keV kg yr)], all chains, budget %.0e. projected = MC below %.1f m + the fitted reach above;\n"
         "    MC 90%% = the MC alone, unseen slabs at 2.30 x their heaviest decay: the proof, once statistics allow\n", cfg::rtBudget, cfg::reachEnd / 1000);
  printf("      %-36s %6s %6s %6s %6s %6s %10s %10s %8s %11s %9s\n", "design", "L1 [m]", "L2 [m]", "steel", "Cu", "EFCu", "projected",
         "proj 90%", "x budget", "after: 90%", "MC 90% x");
  for (auto &q : design)
  {
    if (q.L1 < 0) { printf("      %-36s no steel passes at this L2\n", q.name); continue; }
    double pc = project(q.L1, q.L2, ORD, false, false), ph = project(q.L1, q.L2, ORD, false, true), pa = project(q.L1, q.L2, ORD, true, true);
    printf("      %-36s %6.2f %6.2f %6.0f %6.0f %6.0f %10.2e %10.2e %8s %11.2e %9s\n", q.name, q.L1 / 1000, q.L2 / 1000, massOf(0, q.L1, 0),
           massOf(q.L1, q.L2, 1), massOf(q.L2, depthTot, 2), pc, ph, xb(ph), pa, xb(judge(q.L1, q.L2, ORD, false).up));
  }
  printf("      masses in kg; steel includes the lid, which sits in the top slab\n");

//-------------------------------------------------------------------------------
//  9. The trade-off: for each Cu/EFCu seam L2, the deepest steel that still passes
  std::vector<double> rows = {L2ks, L2best, L2ksSteel};
  for (double L2 = 3000; L2 < depthTot; L2 += 250) rows.push_back(L2);
  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  printf("\n[9] the trade-off: for each L2, the deepest steel/Cu seam L1 whose projected 90%% BI is under the budget (before cuts: cuts only lower it)\n");
  printf("      %-7s %9s %9s %10s %9s %10s %8s %11s %9s\n", "L2 [m]", "EFCu [kg]", "steel to", "steel [kg]", "Cu [kg]", "proj 90%",
         "x budget", "after: 90%", "MC 90% x");
  for (double L2 : rows)
  {
    double l1 = deepestSteel(L2);
    if (l1 < 0)
    {
      printf("      %-7.2f %9.0f %9s   no steel passes: EFCu must reach higher\n", L2 / 1000, massOf(L2, depthTot, 2), "-");
      continue;
    }
    double ph = project(l1, L2, ORD, false, true);
    printf("      %-7.2f %9.0f %7.2f m %10.0f %9.0f %10.2e %8s %11.2e %9s%s\n", L2 / 1000, massOf(L2, depthTot, 2), l1 / 1000, massOf(0, l1, 0),
           massOf(l1, L2, 1), ph, xb(ph), project(l1, L2, ORD, true, true), xb(judge(l1, L2, ORD, false).up),
           L2 == L2ks ? "   <- KS EFCu" : (L2 == L2best ? "   <- least EFCu" : (L2 == L2ksSteel ? "   <- least EFCu with KS steel" : "")));
  }
  printf("      KS: steel to %.2f m (%.0f kg), EFCu from %.2f m (%.0f kg)\n", L1ks / 1000, massOf(0, L1ks, 0), L2ks / 1000, massOf(L2ks, depthTot, 2));

//-------------------------------------------------------------------------------
//  10. Survival in the background window, per chain and section:
  auto surv = [&](long pass, long of) -> std::string { // percent, binomial error; a limit when nothing survives
    if (of <= 0) return "-";
    if (pass == 0) return Form("< %.1f", 100.0 * UL90 / of);
    double f = (double)pass / of;
    return Form("%.1f +- %.1f", 100 * f, 100 * sqrt(f * (1 - f) / of));
  };
  std::vector<std::array<long, 4>> tube; // per chain, whole tube: N0, AC, PSD, AC + PSD
  printf("\n[10] survival in the background window, counted in hits as Edgar does [%%]:\n");
  printf("      %-6s %-6s %5s   %-16s %-16s %-16s %-16s\n", "chain", "sect", "N0", "AC (M1+argon)", "PSD", "PSD | AC", "Combined");
  for (auto &r : R)
  {
    long n[5][4] = {{0}}; // [row, 4 = whole tube][N0, AC, PSD, AC + PSD]
    for (auto &h : r.hits)
      if (cfg::inWindow(h.e))
        for (int k : {rowOf(r.depth[h.ev], L1ks, L2ks, Lhks), 4})
        {
          n[k][0]++;
          n[k][1] += h.m1 && h.lar;
          n[k][2] += h.psd == 1;
          n[k][3] += h.m1 && h.lar && h.psd == 1;
        }
    for (int k = 0; k < 5; k++)
    {
      printf("      %-6s %-6s %5ld   %-16s", k ? "" : r.iso.c_str(), rowName[k], n[k][0], surv(n[k][1], n[k][0]).c_str());
      if (r.hasPsd) // PSD|AC is a step: of the hits AC kept, how many PSD keeps too
        printf(" %-16s %-16s %-16s\n", surv(n[k][2], n[k][0]).c_str(), surv(n[k][3], n[k][1]).c_str(), surv(n[k][3], n[k][0]).c_str());
      else
        printf(" (PSD columns need response.py)\n");
    }
    tube.push_back({n[4][0], n[4][1], n[4][2], n[4][3]});
  }

//-------------------------------------------------------------------------------
//  11. Background index, per chain and section:
  struct Sum { double w = 0, w2 = 0, wbar = 0; long n = 0; }; // one section: summed weight, its square for Poisson, mean weight of every decay there, window hits
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
    if (t.n == 0) return "no window hit";
    std::string out = Form("%s%.2e +-%.1e +-%.1e", t.ul ? "<" : "", t.v, sqrt(t.va), sqrt(t.vs));
    return goal ? out + Form("  (%s%s x goal)", t.ul ? "<" : "", xg(t.v)) : out;
  };
  printf("\n[11] background index, KS as built [cts/(keV kg yr)] = BI +- radioassay +- MC statistics:\n");
  printf("      %-6s %-6s %9s   %-34s %-34s\n", "chain", "sect", "win hits", "before cuts", "after cuts");
  Tot gB, gA;
  struct Pure { std::string iso; int sec; Sum b, a; };
  std::vector<Pure> pure; // kept for [12]
  for (auto &r : R)
  {
    Tot tB, tA;
    for (int sec = 0; sec < 4; sec++)
    {
      int m = rowMat[sec];
      double rel = cfg::relUnc(r.iso, m);
      bool ul = cfg::upperLimit(r.iso, m);
      Sum qb, qa; // before and after cuts
      long nsec = 0;
      for (int p = 0; p < 3; p++) // mean weight over every decay in the section, from the counts
      {
        nsec += r.nRowPv[sec][p];
        qb.wbar += r.nRowPv[sec][p] * weight(r, m, p) / winW;
      }
      qa.wbar = qb.wbar = nsec ? qb.wbar / nsec : 0;
      for (Long64_t i : r.winIdx)
      {
        if (rowOf(r.depth[i], L1ks, L2ks, Lhks) != sec)
          continue;
        double w = weight(r, m, r.pv[i]) / winW; // per decay: the row can mix volumes
        if (r.winRaw[i]) { qb.w += r.winRaw[i] * w; qb.w2 += r.winRaw[i] * w * w; qb.n += r.winRaw[i]; }
        if (r.winCut[i]) { qa.w += r.winCut[i] * w; qa.w2 += r.winCut[i] * w * w; qa.n += r.winCut[i]; }
      }
      printf("      %-6s %-6s %4ld/%-4ld   %-34s %-34s%s\n", sec ? "" : r.iso.c_str(), rowName[sec], qa.n, qb.n,
             cell(qb, rel, ul).c_str(), cell(qa, rel, ul).c_str(), ul ? "  activity is an upper limit" : "");
      for (Tot *t : {&tB, &gB}) add(*t, qb, rel, ul);
      for (Tot *t : {&tA, &gA}) add(*t, qa, rel, ul);
      pure.push_back({r.iso, sec, qb, qa});
    }
    printf("      %-6s %-6s %9s   %-34s %-34s\n", "", "tube", "", show(tB, false).c_str(), show(tA, false).c_str());
  }
  printf("      %-13s %9s   %-34s %-34s\n", "ALL CHAINS", "", show(gB, true).c_str(), show(gA, true).c_str());

//-------------------------------------------------------------------------------
//  12. Purity requirement, per chain and section:
  auto allowed = [&](const Sum &q, double A) -> std::string { // BI is proportional to activity, so the activity that gives the goal is A x goal / BI
    if (q.n == 0) return Form("> %.2g", A * cfg::bgGoal / (UL90 * q.wbar)); // no window hit: the MC can only vouch for activities up to this
    return Form("%.2g +- %.1g", A * cfg::bgGoal / q.w, A * cfg::bgGoal / q.w * sqrt(q.w2) / q.w);
  };
  printf("\n[12] purity requirement: specific activity [uBq/kg] at which a section ALONE gives the goal:\n");
  printf("      %-6s %-6s %9s   %-16s %-16s\n", "chain", "sect", "assumed", "before cuts", "after cuts");
  for (auto &q : pure)
  {
    double A = cfg::activity(q.iso, rowMat[q.sec]) * 1e6;
    bool blind = q.b.n == 0 && cfg::bgGoal < UL90 * q.b.wbar; // no hit, and one hit would already be over the goal
    printf("      %-6s %-6s %9s   %-16s %-16s%s\n", q.sec ? "" : q.iso.c_str(), rowName[q.sec],
           Form("%s%g", cfg::upperLimit(q.iso, rowMat[q.sec]) ? "< " : "", A), allowed(q.b, A).c_str(), allowed(q.a, A).c_str(),
           blind ? "  assumed activity is above what this run can exclude: more decays needed" : "");
  }

//-------------------------------------------------------------------------------
//  13. Per detector, as built: the line each chain gives, the tracer MAJORANA used to locate a source [arXiv:2209.10592]
  std::vector<TH2D *> hDet; // string x position map of the line rate, one per chain
  printf("\n[13] per detector, as built, before cuts: counts per year in each chain's strongest line (+-5 keV), the tube's fingerprint in data\n");
  for (auto &r : R)
  {
    int is = cfg::iso(r.iso);
    std::vector<double> rate(r.dets.size(), 0), raw(r.dets.size(), 0), byPos(9, 0);
    for (auto &h : r.hits)
      if (fabs(h.e - cfg::line[is]) <= 5)
      {
        double cy = weight(r, ORD[slabOf(r.depth[h.ev], L1ks, L2ks)], r.pv[h.ev]) * cfg::geMass_kg; // real counts per year in the array
        rate[h.det] += cy;
        raw[h.det]++;
      }
    auto hm = new TH2D(Form("hDet%d", is), Form("%s, %.1f keV line from the tube as built;string;position (1 = top)", r.iso.c_str(), cfg::line[is]),
                       42, 0.5, 42.5, 8, 0.5, 8.5);
    double tot = 0, totRaw = 0;
    for (size_t k = 0; k < r.dets.size(); k++)
    {
      int st = atoi(r.dets[k].substr(1, 2).c_str()), ps = atoi(r.dets[k].substr(3, 2).c_str()); // V<string><position>
      hm->Fill(st, 9 - ps, rate[k]); // position 1 is the top of the string: draw it at the top
      if (ps >= 1 && ps <= 8) byPos[ps] += rate[k];
      tot += rate[k];
      totRaw += raw[k];
    }
    for (int j = 1; j <= 8; j++) hm->GetYaxis()->SetBinLabel(j, Form("%d", 9 - j)); // row j holds position 9 - j
    hDet.push_back(hm);
    printf("      %s %.1f keV: %.2e counts/yr in the whole array, from %.0f MC events\n", r.iso.c_str(), cfg::line[is], tot, totRaw);
    printf("        by position, top to bottom [counts/yr]:");
    for (int ps = 1; ps <= 8; ps++) printf(" %.1e", byPos[ps]);
    std::vector<size_t> order(r.dets.size());
    for (size_t k = 0; k < order.size(); k++) order[k] = k;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return rate[a] > rate[b]; });
    printf("\n        brightest detectors:");
    for (int k = 0; k < 5 && k < (int)order.size(); k++) printf("  %s %.1e (%.0f MC)", r.dets[order[k]].c_str(), rate[order[k]], raw[order[k]]);
    printf("\n");
  }

//-------------------------------------------------------------------------------
//  14. Draw: spectrum, seam scan, trade-off, survival, raw vs weighted per chain x material, per-detector map
  const int NS = 60;
  auto hScan = new TH2D("hScan", "projected 90% BI before cuts, x budget;L1 steel/Cu seam [m];L2 Cu/EFCu seam [m]",
                        NS, 0, depthTot / 1000, NS, 0, depthTot / 1000);
  for (int i = 1; i <= NS; i++)
    for (int j = i; j <= NS; j++) // only L1 <= L2 is a design
      hScan->SetBinContent(i, j, project(hScan->GetXaxis()->GetBinCenter(i) * 1000, hScan->GetYaxis()->GetBinCenter(j) * 1000, ORD, false, true) / cfg::rtBudget);
  auto gEdge = new TGraph(); // the deepest steel that passes, per L2: the edge of the allowed region in the scan
  auto gTrade = new TGraph(); // the same as masses: steel kg against EFCu kg
  for (double L2 = 100; L2 <= depthTot; L2 += 50)
  {
    double l1 = deepestSteel(L2);
    if (l1 < 0) continue;
    gEdge->SetPoint(gEdge->GetN(), l1 / 1000, L2 / 1000);
    gTrade->SetPoint(gTrade->GetN(), massOf(L2, depthTot, 2), massOf(0, l1, 0));
  }
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
  auto c = new TCanvas("c_bkg", "", 1920, 560);
  c->Divide(3, 1);
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
  gEdge->SetLineColor(kRed + 1);
  gEdge->SetLineWidth(3);
  if (gEdge->GetN()) gEdge->Draw("L SAME");
  auto star = new TMarker(L1ks / 1000, L2ks / 1000, 29); // the KS design
  dot(star, kRed + 1, 29, 2.2);
  star->Draw();
  c->cd(3);
  gPad->SetGrid();
  gPad->SetLeftMargin(0.15);
  double kgMax = 2.5 * massOf(L2ks, depthTot, 2), stMax = 1.3 * massOf(0, depthTot, 0); // the decision is around KS
  auto fr = gPad->DrawFrame(0, 0, kgMax, stMax);
  fr->SetTitle("the cheapest tubes: most steel for each EFCu mass;EFCu [kg];steel [kg]");
  fr->GetYaxis()->SetTitleOffset(1.6);
  gTrade->SetLineColor(kRed + 1);
  gTrade->SetLineWidth(3);
  if (gTrade->GetN()) gTrade->Draw("L SAME");
  auto leg3 = new TLegend(0.35, 0.74, 0.88, 0.88);
  leg3->SetTextSize(0.032);
  leg3->AddEntry(gTrade, "projected 90% BI = budget", "l");
  struct { double L1, L2; int col, style; const char *lab; } pt[3] = {{L1ks, L2ks, kRed + 1, 29, "KS"}, {L1best, L2best, kGreen + 2, 20, "least EFCu, then most steel"},
                                                                    {L1ks, L2ksSteel, kAzure + 2, 21, "KS steel, least EFCu"}};
  for (auto &q : pt)
  {
    if (q.L1 < 0) continue;
    auto m = new TMarker(massOf(q.L2, depthTot, 2), massOf(0, q.L1, 0), q.style);
    dot(m, q.col, q.style, 2.0);
    m->Draw();
    leg3->AddEntry(m, q.lab, "p");
  }
  leg3->Draw();
  std::string tag;
  for (auto &r : R) tag += (tag.empty() ? "" : "_") + r.iso;
  std::string png = rtOut(tag + "_background.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s\n", png.c_str());

  auto cs = new TCanvas("c_surv", "", 560 * R.size(), 480); // survival in the window against the RE-vessel references
  cs->Divide(R.size(), 1);
  const char *cutName[4] = {"AC (M1+argon)", "PSD", "PSD | AC", "Combined"};
  for (size_t i = 0; i < R.size(); i++)
  {
    cs->cd(i + 1)->SetLogy();
    gPad->SetGrid(0, 1);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    int is = cfg::iso(R[i].iso);
    auto bars = new TH1D(Form("cdr%zu", i), Form("%s from the RT wall;;survival [%%]", R[i].iso.c_str()), 4, 0, 4); // the CDR
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
      if ((k > 0 && !R[i].hasPsd) || of[k] <= 0) // PSD points need response.py
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
    lg->AddEntry(gm, Form("this sim - KS tube (%ld window hits)", n[0]), "pe");
    if (!lim.empty()) lg->AddEntry(lim[0], "this sim - none survived (90% CL)", "l");
    lg->Draw();
  }
  std::string spng = rtOut(tag + "_survival.png");
  cs->SaveAs(spng.c_str());
  printf("wrote %s\n", spng.c_str());

  // per chain x material, as built: raw MC counts (solid, left axis) beside the weighted rate (dashed, right axis).
  // a single weighted bin can rest on one raw event: the raw curve shows how many sit under each
  const char *colName[4] = {"steel", "Cu", "EFCu", "the whole tube"};
  auto cm = new TCanvas("c_mat", "", 520 * 4, 430 * R.size());
  cm->Divide(4, R.size());
  for (size_t i = 0; i < R.size(); i++)
    for (int col = 0; col < 4; col++) // the three slabs as built, then all of them
    {
      const Run &r = R[i];
      const char *stage[2] = {"M1", r.hasPsd ? "M1+AC+PSD" : "M1+AC"};
      TH1D *raw[2], *wtd[2];
      for (int k = 0; k < 2; k++)
      {
        raw[k] = new TH1D(Form("raw%zu%d%d", i, col, k), "", cfg::eBins, cfg::eLo, cfg::eHi);
        wtd[k] = new TH1D(Form("wtd%zu%d%d", i, col, k), "", cfg::eBins, cfg::eLo, cfg::eHi);
      }
      for (auto &h : r.hits)
      {
        int sl = slabOf(r.depth[h.ev], L1ks, L2ks);
        if (col < 3 && sl != col) continue;
        double w = weight(r, ORD[sl], r.pv[h.ev]) / binW;
        bool pass[2] = {h.m1, h.m1 && h.lar && (!r.hasPsd || h.psd == 1)};
        for (int k = 0; k < 2; k++)
          if (pass[k]) { raw[k]->Fill(h.e); wtd[k]->Fill(h.e, w); }
      }
      cm->cd(i * 4 + col + 1)->SetLogy();
      gPad->SetGrid();
      gPad->SetLeftMargin(0.13);
      gPad->SetRightMargin(0.15);
      std::string title = Form("%s in %s", r.iso.c_str(), colName[col]);
      if (raw[0]->GetEntries() == 0) // nothing came from here: say so instead of drawing an empty frame
      {
        gPad->SetLogy(0);
        gPad->DrawFrame(cfg::eLo, 0, cfg::eHi, 1)->SetTitle(Form("%s;energy [keV];MC counts", title.c_str()));
        auto t = new TLatex(0.5 * (cfg::eLo + cfg::eHi), 0.5, Form("no hit in %lld decays", r.nsim));
        t->SetTextAlign(22);
        t->SetTextSize(0.06);
        t->Draw();
        continue;
      }
      double k = raw[0]->Integral() / wtd[0]->Integral(); // maps the weighted rate onto the counts axis; the right axis undoes it
      double top = 1;
      for (int s = 0; s < 2; s++)
      {
        wtd[s]->Scale(k);
        top = std::max({top, raw[s]->GetMaximum(), wtd[s]->GetMaximum()});
      }
      double y0 = 0.5, y1 = 40 * top; // headroom above the highest line for the legend
      auto fr = gPad->DrawFrame(cfg::eLo, y0, cfg::eHi, y1);
      fr->SetTitle(Form("%s;energy [keV];MC counts", title.c_str()));
      auto band = new TBox(cfg::winLo, y0, cfg::winHi, y1); // the background window (its four 10 keV gaps are too thin to draw)
      band->SetFillColorAlpha(kRed, 0.12);
      band->Draw();
      auto lg = new TLegend(0.14, 0.74, 0.84, 0.89);
      lg->SetNColumns(2);
      lg->SetTextSize(0.032);
      const int col2[2] = {kAzure + 2, kRed + 1};
      for (int s = 0; s < 2; s++)
      {
        raw[s]->SetLineColor(col2[s]);
        raw[s]->SetLineWidth(2);
        raw[s]->Draw("HIST SAME");
        wtd[s]->SetLineColor(col2[s]);
        wtd[s]->SetLineStyle(2);
        wtd[s]->Draw("HIST SAME");
        lg->AddEntry(raw[s], Form("%s raw", stage[s]), "l");
        lg->AddEntry(wtd[s], Form("%s weighted", stage[s]), "l");
      }
      lg->AddEntry(band, "BI window", "f");
      lg->Draw();
      auto ax = new TGaxis(cfg::eHi, y0, cfg::eHi, y1, y0 / k, y1 / k, 510, "+LG"); // plain values: ROOT maps them onto the log pad
      ax->SetTitle("weighted: cts / (keV kg yr)");
      ax->SetLabelSize(0.035);
      ax->SetTitleSize(0.035);
      ax->SetTitleOffset(1.4);
      ax->Draw();
    }
  std::string mpng = rtOut(tag + "_spectra.png");
  cm->SaveAs(mpng.c_str());
  printf("wrote %s\n", mpng.c_str());

  auto cd = new TCanvas("c_det", "", 760 * hDet.size(), 420);
  cd->Divide(hDet.size(), 1);
  for (size_t i = 0; i < hDet.size(); i++)
  {
    cd->cd(i + 1)->SetLogz();
    gPad->SetRightMargin(0.16);
    hDet[i]->SetMinimum(0.5 * hDet[i]->GetMinimum(0));
    hDet[i]->GetZaxis()->SetTitle("counts / yr");
    hDet[i]->Draw("COLZ");
  }
  std::string dpng = rtOut(tag + "_detectors.png");
  cd->SaveAs(dpng.c_str());
  printf("wrote %s\n", dpng.c_str());
}
