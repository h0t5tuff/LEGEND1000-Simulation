//  how much background the re-entrant tube puts into the germanium, and whether the
//  choice and placement of its materials matters
//    root -l -b -q 'sim_background.C("output/tl208.root=Tl208")'
//    root -l -b -q 'sim_background.C("output/tl208.root=Tl208,output/bi214.root=Bi214")'
//    root -l -b -q 'sim_background.C("output/tl208_*.root=Tl208,output/bi214_*.root=Bi214")'   job arrays, merged per isotope
//  takes a comma-separated list of run=isotope. each run is weighted by its own nuclide's
//  activity and its own event count, then the rates are summed: the background is the sum
//  over decay chains, and the 1e-5 goal is measured against that total
//  if sim_psd.py has been run, its output/<run>_psd.csv supplies the detector response:
//  active energy after the dead layer, and AoE_class for the PSD cut. without it the study
//  stops at M1 + argon veto and says so

#include "rt_geom.h"
#include "TRegexp.h"
#include <sstream>
#include <vector>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  const double qbb = 2039.0, roiHalf = 55.0;  // keV, the neutrinoless double beta window
  const double m1_keV = 5.0;                  // a detector counts as fired above this
  const double lar_keV = 20.0;                // argon veto: fires above this deposit. 20 keV is where the 4 PE cut sits
  const double geMass_kg = 1000.0;            // the array being protected
  const double exposureYr = 10.0;             // live time to quote real-life event counts over
  const double decaysPerSec = 50.0;           // this laptop, -t 8. only used to turn a required N into a wall time
  const double bgGoal = 1e-5;                 // LEGEND-1000 target [cts/(keV kg yr)] in the ROI
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000;        // spectrum range [keV]
  const int    eBins = 200;
  const int    NDB = 10;                      // depth slabs for the reach profile

  const char  *mat[3] = {"steel", "Cu", "EFCu"};             // ordered top to bottom
  const double density[3] = {7900.0, 8960.0, 8930.0};        // kg/m^3
  const double psdCut = -1.80;                // AoE_class cut in Edgar's analysis. his aoe_class_paras.yaml says -0.83

  // specific activity [uBq/kg], its uncertainty, and whether it is only an upper limit.
  // rows steel / Cu / EFCu, columns Bi214 (U238 chain) / Tl208 (Th232 chain)
  //   steel, Cu : Ralph's materialMix, no uncertainty quoted
  //   EFCu      : Edgar's survival_BI.py radioassay; Tl208 there is an upper limit
  struct Act { double v, dv; bool ul; };
  const Act act[3][2] = {{{2500, 0, false}, {1000, 0, false}},
                         {{1, 0, false}, {1, 0, false}},
                         {{0.19, 0.10, false}, {0.077, 0, true}}};
  int iso(const std::string &s) { return s.find("Bi") != std::string::npos ? 0 : 1; }
  double activity(const std::string &s, int m) { return act[m][iso(s)].v * 1e-6; } // -> Bq/kg
  double relUnc(const std::string &s, int m) { const Act &a = act[m][iso(s)]; return a.v > 0 ? a.dv / a.v : 0; }
  bool upperLimit(const std::string &s, int m) { return act[m][iso(s)].ul; }

  // RE vessel (EFCu) survival in the ROI [%] as AC, PSD, PSD|AC, from survival_BI.py
  const double cdr[2][3] = {{21.0, 21.0, 19.0}, {1.2, 31.0, 16.0}};          // CDR
  const double edgar[2][3] = {{16.94, 18.12, 12.27}, {1.10, 30.32, 22.72}}; // Edgar's remage
}

static bool rt_is_germanium_tree(const TString &n) // NtupleUseVolumeName true -> "V0101"
{
  if (n.Length() != 5 || n[0] != 'V')
    return false;
  for (int i = 1; i < 5; i++)
    if (!isdigit(n[i]))
      return false;
  return true;
}

//-------------------------------------------------------------------------------
//  2. One run, reduced to what a reweighting needs:
struct Hit // one detector's response to one decay
{
  int ev;
  short det;
  float e;             // active energy if sim_psd.py ran, else the summed deposit
  bool m1, lar;        // passes M1, passes the argon veto
  signed char psd;     // 1 passes PSD, 0 fails, -1 unknown (sim_psd.py not run)
};

struct Run
{
  std::string file, iso;
  Long64_t nsim = 0;
  int ndet = 0;
  bool hasPsd = false;                // detector response supplied by sim_psd.py
  std::vector<std::string> dets;
  std::vector<float> depth;           // mm below the top of the tube, indexed by evtid
  std::vector<signed char> pv;        // physical volume it was sampled in: 0 mother, 1 OFHC, 2 SS, -1 none
  double mcDensity[3] = {0, 0, 0};    // simulated decays per m^3 in each physical volume, measured
  long nBad = 0, nFiles = 1;          // decays outside the wall (a confinement bug); files merged in
  std::vector<short> roiRaw, roiCut;  // ROI hits per decay: before cuts, after the full chain
  std::vector<Hit> hits;
  long nHit[4] = {0, 0, 0, 0};        // no cuts, M1, M1 + argon, M1 + argon + PSD
  long dN[cfg::NDB] = {0}, dH[cfg::NDB] = {0}, dR[cfg::NDB] = {0};
  bool ok = false;
};

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
  r.nsim = vtx->GetEntries();
  std::vector<float> eLAr(r.nsim, 0);
  std::vector<short> nFired(r.nsim, 0);
  r.depth.assign(r.nsim, -1);
  r.roiRaw.assign(r.nsim, 0);
  r.roiCut.assign(r.nsim, 0);

  // germanium: each detector's deposits summed per decay
  struct Resp { int ev; short det; float e; };
  std::vector<Resp> resp;
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it())
  {
    TString nm = k->GetName();
    if (!rt_is_germanium_tree(nm))
      continue;
    auto t = (TTree *)d->Get(nm);
    if (!t || !t->GetBranch("edep_in_keV"))
      continue;
    short id = r.dets.size();
    r.dets.push_back(nm.Data());
    r.ndet++;
    TTreeFormula fe("fe", "evtid", t), fd("fd", "edep_in_keV", t);
    std::map<int, float> perEvent;
    for (Long64_t j = 0; j < t->GetEntries(); j++)
    {
      t->GetEntry(j);
      perEvent[(int)fe.EvalInstance()] += fd.EvalInstance();
    }
    for (auto &p : perEvent)
      if (p.first >= 0 && p.first < r.nsim && p.second > 0)
        resp.push_back({p.first, id, p.second});
  }

  // the detector response from sim_psd.py, joined on (event, detector)
  std::string base = file.substr(file.find_last_of('/') + 1);
  base = base.substr(0, base.rfind(".root"));
  std::map<std::pair<int, std::string>, std::pair<float, float>> psd; // -> (active energy, AoE_class)
  std::ifstream in("output/" + base + "_psd.csv");
  if (in)
  {
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line))
    {
      std::stringstream ss(line);
      std::string ev, det, e, a;
      std::getline(ss, ev, ',');
      std::getline(ss, det, ',');
      std::getline(ss, e, ',');
      std::getline(ss, a, ',');
      psd[{std::stoi(ev), det}] = {std::stof(e), std::stof(a)};
    }
    r.hasPsd = !psd.empty();
  }
  for (auto &p : resp)
  {
    float e = p.e;
    signed char pass = -1;
    if (r.hasPsd)
    {
      auto h = psd.find({p.ev, r.dets[p.det]});
      if (h == psd.end())
        continue; // sim_psd.py dropped it: every deposit sat in the dead layer
      e = h->second.first;
      pass = h->second.second > cfg::psdCut ? 1 : 0;
    }
    if (e <= cfg::m1_keV)
      continue;
    nFired[p.ev]++; // a detector fired: M1 counts these
    r.hits.push_back({p.ev, p.det, e, false, false, pass});
  }

  if (auto tl = (TTree *)d->Get("undergroundlar"))
  {
    TTreeFormula fe("fe", "evtid", tl), fd("fd", "edep_in_keV", tl);
    for (Long64_t j = 0; j < tl->GetEntries(); j++)
    {
      tl->GetEntry(j);
      int ev = (int)fe.EvalInstance();
      if (ev >= 0 && ev < r.nsim)
        eLAr[ev] += fd.EvalInstance();
    }
  }
  // vtx rows are NOT in event order: with -t 8 remage writes them as threads finish, so
  // row 7 can hold evtid 2242. Index by evtid or every decay gets someone else's vertex.
  r.pv.assign(r.nsim, -1);
  TTreeFormula vx("vx", "xloc_in_m", vtx), vy("vy", "yloc_in_m", vtx), vz("vz", "zloc_in_m", vtx),
      ve("ve", "evtid", vtx);
  for (Long64_t i = 0; i < r.nsim; i++)
  {
    vtx->GetEntry(i);
    int ev = (int)ve.EvalInstance();
    if (ev < 0 || ev >= r.nsim)
      continue;
    double x = vx.EvalInstance() * 1000, y = vy.EvalInstance() * 1000, z = vz.EvalInstance() * 1000;
    r.depth[ev] = rt.zTop - z;
    r.pv[ev] = rt.physAt(sqrt(x * x + y * y), z); // the volume remage sampled it from
    if (r.pv[ev] < 0) r.nBad++;
  }

  const double depthTot = rt.zTop - rt.zBottom;
  for (auto &h : r.hits)
  {
    h.m1 = (nFired[h.ev] == 1);
    h.lar = (eLAr[h.ev] <= cfg::lar_keV);
    bool full = h.m1 && h.lar && (!r.hasPsd || h.psd == 1); // the whole chain this run can apply
    bool inRoi = fabs(h.e - cfg::qbb) <= cfg::roiHalf;
    r.nHit[0]++;
    if (h.m1) r.nHit[1]++;
    if (h.m1 && h.lar) r.nHit[2]++;
    if (r.hasPsd && full) r.nHit[3]++;
    if (inRoi) { r.roiRaw[h.ev]++; if (full) r.roiCut[h.ev]++; }
  }
  for (Long64_t i = 0; i < r.nsim; i++)
  {
    int b = std::min(cfg::NDB - 1, std::max(0, (int)(r.depth[i] / depthTot * cfg::NDB)));
    r.dN[b]++;
    r.dR[b] += r.roiRaw[i];
  }
  for (auto &h : r.hits)
    r.dH[std::min(cfg::NDB - 1, std::max(0, (int)(r.depth[h.ev] / depthTot * cfg::NDB)))]++;
  r.ok = true;
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
void sim_background(const char *runs = "output/tl208.root=Tl208,output/bi214.root=Bi214",
                    const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  if (!rt.ok)
    return;

  std::vector<Run> R; // split "a.root=Tl208,b.root=Bi214" on commas, then each on '='
  std::string all(runs);
  for (size_t p = 0; p <= all.size();)
  {
    size_t c = all.find(',', p);
    std::string tok = all.substr(p, c == std::string::npos ? std::string::npos : c - p);
    p = (c == std::string::npos) ? all.size() + 1 : c + 1;
    if (tok.empty())
      continue;
    size_t eq = tok.find('=');
    std::string file = eq == std::string::npos ? tok : tok.substr(0, eq);
    std::string iso = eq == std::string::npos ? "Tl208" : tok.substr(eq + 1);
    for (auto &f : expand(file))
    {
      Run r = readRun(f, iso, rt);
      if (r.ok)
        R.push_back(std::move(r));
    }
  }
  if (R.empty())
  {
    printf("no usable runs in '%s'\n", runs);
    return;
  }

  // files of one isotope are one run split into jobs: stitch them, events renumbered, or it counts once per file
  std::vector<Run> M;
  for (auto &r : R)
  {
    Run *into = nullptr;
    for (auto &m : M)
      if (m.iso == r.iso) into = &m;
    if (!into) { M.push_back(std::move(r)); continue; }
    int off = into->nsim;
    for (auto h : r.hits) { h.ev += off; into->hits.push_back(h); }
    into->depth.insert(into->depth.end(), r.depth.begin(), r.depth.end());
    into->pv.insert(into->pv.end(), r.pv.begin(), r.pv.end());
    into->roiRaw.insert(into->roiRaw.end(), r.roiRaw.begin(), r.roiRaw.end());
    into->roiCut.insert(into->roiCut.end(), r.roiCut.begin(), r.roiCut.end());
    for (int k = 0; k < 4; k++) into->nHit[k] += r.nHit[k];
    for (int b = 0; b < cfg::NDB; b++) { into->dN[b] += r.dN[b]; into->dH[b] += r.dH[b]; into->dR[b] += r.dR[b]; }
    if (into->hasPsd != r.hasPsd)
      printf("WARNING: %s files disagree on sim_psd.py output - run it on every file\n", r.iso.c_str());
    into->hasPsd = into->hasPsd && r.hasPsd;
    into->nsim += r.nsim;
    into->nBad += r.nBad;
    into->nFiles++;
    into->file += "+" + r.file;
  }
  R = std::move(M);

  const double mm3 = 1e-9;
  const double vTot = rt.wallVolume(rt.zBottom, rt.zTop) * mm3;
  double physV[3]; // m^3 of the mother, the OFHC shell, the SS shell
  for (int p = 0; p < 3; p++) physV[p] = rt.physVolume(p) * mm3;
  // remage fills the mother ~19% sparser than its two daughters, so measure each volume's density, never assume it
  for (auto &r : R)
  {
    long n[3] = {0, 0, 0};
    for (Long64_t i = 0; i < r.nsim; i++)
      if (r.pv[i] >= 0) n[(int)r.pv[i]]++;
    for (int p = 0; p < 3; p++) r.mcDensity[p] = n[p] / physV[p];
  }
  const double depthTot = rt.zTop - rt.zBottom;
  const double L1ks = rt.zTop - rt.seamSS, L2ks = rt.zTop - rt.seamOFHC;
  const double binW = (cfg::eHi - cfg::eLo) / cfg::eBins;
  const double UL90 = 2.30;

  printf("geometry  : %s\n", rt.file.c_str());
  printf("RT wall   : %.5f m^3 over %.3f m\n", vTot, depthTot / 1000);
  printf("cuts      : M1 (> %.0f keV) + argon veto (<= %.0f keV) + PSD (AoE_class > %.2f)   ROI %.0f +/- %.0f keV\n",
         cfg::m1_keV, cfg::lar_keV, cfg::psdCut, cfg::qbb, cfg::roiHalf);
  for (auto &r : R)
    printf("            %-6s %s\n", r.iso.c_str(), r.hasPsd
               ? "detector response from sim_psd.py: active energy + PSD"
               : "NO sim_psd.py output: summed deposit, and the chain stops at M1 + argon");
  printf("\n");

  // one simulated decay = rho x A x T real decays per m^3 of material m, over the measured MC density
  // of the volume pv it was drawn in. m is what the design puts at that depth, so a design is arithmetic
  auto weight = [&](const Run &r, int m, int pv) {
    if (pv < 0 || r.mcDensity[pv] <= 0) return 0.0;
    return cfg::density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[pv] * cfg::geMass_kg);
  };
  auto slabOf = [&](double dep, double L1, double L2) { return dep < L1 ? 0 : (dep < L2 ? 1 : 2); };
  // ROI rate summed over every run, for a design given by boundaries and material order
  auto rate = [&](double L1, double L2, const int *ord, bool cut, double *err, long *nmc) {
    double sum = 0, var = 0;
    long n = 0;
    for (auto &r : R)
    {
      const std::vector<short> &roi = cut ? r.roiCut : r.roiRaw;
      for (Long64_t i = 0; i < r.nsim; i++)
      {
        if (!roi[i])
          continue;
        double w = weight(r, ord[slabOf(r.depth[i], L1, L2)], r.pv[i]) / (2 * cfg::roiHalf);
        sum += roi[i] * w;
        var += roi[i] * w * w; // Poisson: var(N w) = N w^2
        n += roi[i];
      }
    }
    if (err) *err = sqrt(var);
    if (nmc) *nmc = n;
    return sum;
  };
  const int ORD[3] = {0, 1, 2};

//-------------------------------------------------------------------------------
//  4. Per run: cuts, normalisation, reach, sizing
  for (auto &r : R)
  {
    printf("=== %s   (%s, %lld decays, %d Ge tables) ===\n", r.iso.c_str(), r.file.c_str(), r.nsim, r.ndet);
    printf("[a] hits surviving each cut: no cuts %ld -> M1 %ld -> +argon %ld", r.nHit[0], r.nHit[1], r.nHit[2]);
    if (r.hasPsd)
      printf(" -> +PSD %ld", r.nHit[3]);
    printf("   (argon veto x%.1f)\n", r.nHit[2] ? (double)r.nHit[1] / r.nHit[2] : 0.0);

    // per physical volume, as built. MC decays and their density are measured, not assumed
    const int asBuilt[3] = {2, 1, 0}; // mother is EFCu, the OFHC shell Cu, the SS shell steel
    const char *pvName[3] = {"mother", "OFHC", "SS"};
    double T = cfg::exposureYr * cfg::secPerYear;
    printf("[b] normalisation:  %-7s %-6s %9s %8s %10s %11s %9s %11s %10s\n", "volume", "mat",
           "V [m^3]", "M [kg]", "A [Bq]", Form("N in %.0f yr", cfg::exposureYr), "MC dec.",
           "MC per m^3", "per MC dec");
    for (int p = 0; p < 3; p++)
    {
      int m = asBuilt[p];
      double mass = physV[p] * cfg::density[m], act = mass * cfg::activity(r.iso, m);
      long nMC = (long)(r.mcDensity[p] * physV[p] + 0.5);
      printf("                    %-7s %-6s %9.5f %8.1f %10.3e %11.3e %9ld %11.0f %10.2f\n", pvName[p],
             cfg::mat[m], physV[p], mass, act, T * act, nMC, r.mcDensity[p], nMC ? T * act / nMC : 0.0);
    }
    double dAvg = 0.5 * (r.mcDensity[1] + r.mcDensity[2]);
    printf("      sampling density mother / shells = %.3f   (1.000 would be uniform; weights use the measured value)\n",
           dAvg > 0 ? r.mcDensity[0] / dAvg : 0.0);
    if (r.nBad)
      printf("      WARNING: %ld decays outside the wall - a confinement bug\n", r.nBad);
    if (r.nFiles > 1)
      printf("      merged from %ld files\n", r.nFiles);

    printf("[c] reach vs depth:  %-16s %9s %8s %13s %9s\n", "depth [m]", "decays", "hits", "hits/decay", "ROI");
    for (int b = 0; b < cfg::NDB; b++)
      printf("      %6.2f..%-9.2f %9ld %8ld %13.2e %9ld%s\n",
             depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000,
             r.dN[b], r.dH[b], r.dN[b] ? (double)r.dH[b] / r.dN[b] : 0.0, r.dR[b],
             b == cfg::NDB - 1 ? "   <- nearest the detectors" : "");

    long roiT = 0;
    for (Long64_t i = 0; i < r.nsim; i++) roiT += r.roiRaw[i];
    bool measured = roiT > 0;
    double eff = measured ? (double)roiT / r.nsim : UL90 / r.nsim;
    printf("[d] sizing: ROI hits per decay %s%.3e (%ld in %lld)", measured ? "" : "< ", eff, roiT, r.nsim);
    for (double p : {0.10, 0.05})
      printf("   %2.0f%%: %s%.1e decays (%.0f h)", 100 * p, measured ? "" : ">",
             1.0 / (p * p) / eff, 1.0 / (p * p) / eff / cfg::decaysPerSec / 3600);
    printf("\n");
    int lastEmpty = -1;
    for (int b = 0; b < cfg::NDB; b++) { if (r.dH[b] == 0 && r.dN[b] > 500) lastEmpty = b; else break; }
    if (lastEmpty >= 0)
    {
      long nDead = 0;
      for (int b = 0; b <= lastEmpty; b++) nDead += r.dN[b];
      double fr = (double)nDead / r.nsim;
      if (fr > 0.05)
        printf("      %.0f%% of decays landed above %.2f m depth and reached nothing: confining the\n"
               "      source below that buys %.1fx the statistics for the same wall time\n",
               100 * fr, depthTot * (lastEmpty + 1) / cfg::NDB / 1000, 1.0 / (1.0 - fr));
    }
    printf("\n");
  }

//-------------------------------------------------------------------------------
//  5. Designs, summed over every run:
  printf("[e] ROI rate [cts/(keV kg yr)], summed over %zu run(s), target %.0e:\n", R.size(), cfg::bgGoal);
  auto xg = [&](double v) { return Form(v / cfg::bgGoal < 10 ? "%.2f" : "%.0f", v / cfg::bgGoal); };
  printf("      %-44s %8s %13s %8s %13s %8s\n", "design", "MC hits", "before cuts", "x goal",
         "after cuts", "x goal");
  struct Scen { const char *name; double L1, L2; int ord[3]; };
  Scen sc[] = {
      {"KS as built: steel / Cu / EFCu", L1ks, L2ks, {0, 1, 2}},
      {"all EFCu", 0, 0, {2, 2, 2}},
      {"all Cu", 0, 0, {1, 1, 1}},
      {"all steel", 0, 0, {0, 0, 0}},
      {"KS, steel reaching twice as far down", 2 * L1ks, L2ks, {0, 1, 2}},
      {"KS, no steel (Cu down to L2)", 0, L2ks, {0, 1, 2}},
      {"SAME slabs, order reversed: EFCu / Cu / steel", L1ks, L2ks, {2, 1, 0}}};
  for (auto &q : sc)
  {
    double e0, e1;
    long n0, n1;
    double r0 = rate(q.L1, q.L2, q.ord, false, &e0, &n0);
    double r1 = rate(q.L1, q.L2, q.ord, true, &e1, &n1);
    printf("      %-44s %3ld/%-4ld %12.2e %8s", q.name, n1, n0, r0, xg(r0));
    if (n1 > 0)
      printf(" %12.2e %8s\n", r1, xg(r1));
    else
    {
      double ul = n0 > 0 ? r0 * UL90 / n0 : 0.0;
      printf(" %12s %8s\n", Form("< %.1e", ul), Form("< %s", xg(ul)));
    }
  }
  if (R.size() > 1)
  {
    printf("\n      per chain, as built:");
    for (auto &r : R)
    {
      double s = 0;
      long n = 0;
      for (Long64_t i = 0; i < r.nsim; i++)
        if (r.roiRaw[i])
        {
          s += r.roiRaw[i] * weight(r, ORD[slabOf(r.depth[i], L1ks, L2ks)], r.pv[i]) / (2 * cfg::roiHalf);
          n += r.roiRaw[i];
        }
      if (n > 0)
        printf("   %s %.2e", r.iso.c_str(), s);
      else
      {
        // no ROI hit from this chain. quote a 90% CL limit using the mean weight of the
        // hits it DID make, so a chain that is merely unmeasured is not reported as zero
        double wbar = 0;
        for (size_t i = 0; i < r.hits.size(); i++)
          wbar += weight(r, ORD[slabOf(r.depth[r.hits[i].ev], L1ks, L2ks)], r.pv[r.hits[i].ev]);
        wbar = r.hits.empty() ? 0 : wbar / r.hits.size();
        printf("   %s < %.2e", r.iso.c_str(), UL90 * wbar / (2 * cfg::roiHalf));
      }
    }
    printf("   (a '<' chain made no ROI hit yet, it is unmeasured, not absent)\n");
  }

//-------------------------------------------------------------------------------
//  6. The background study, per section and per chain:
  auto inRoi = [&](float e) { return fabs(e - cfg::qbb) <= cfg::roiHalf; };
  auto surv = [&](long pass, long of) -> std::string { // percent, binomial error; a limit when nothing survives
    if (of <= 0) return "-";
    if (pass == 0) return Form("< %.1f", 100.0 * UL90 / of);
    double f = (double)pass / of;
    return Form("%.1f +- %.1f", 100 * f, 100 * sqrt(f * (1 - f) / of));
  };
  const char *secName[4] = {"steel", "Cu", "EFCu", "tube"};
  struct Surv { long n[4]; };          // whole tube: N0, AC, PSD, AC+PSD
  std::vector<Surv> tube;

  printf("\n[f] survival in the ROI, counted in hits as Edgar does [%%]:\n");
  printf("      %-6s %-6s %5s   %-16s %-16s %-16s %-16s\n", "chain", "sect", "N0",
         "AC (M1+argon)", "PSD", "PSD | AC", "Combined");
  for (auto &r : R)
  {
    long n[4][4] = {{0}}; // [section, 3 = whole tube][N0, AC, PSD, AC + PSD]
    for (auto &h : r.hits)
    {
      if (!inRoi(h.e))
        continue;
      int sec = slabOf(r.depth[h.ev], L1ks, L2ks);
      for (int k : {sec, 3})
      {
        n[k][0]++;
        if (h.m1 && h.lar) n[k][1]++;
        if (h.psd == 1) n[k][2]++;
        if (h.m1 && h.lar && h.psd == 1) n[k][3]++;
      }
    }
    for (int k = 0; k < 4; k++)
    {
      printf("      %-6s %-6s %5ld   %-16s", k ? "" : r.iso.c_str(), secName[k], n[k][0],
             surv(n[k][1], n[k][0]).c_str());
      if (r.hasPsd) // PSD|AC is a step: of the hits AC kept, how many PSD keeps too
        printf(" %-16s %-16s %-16s\n", surv(n[k][2], n[k][0]).c_str(),
               surv(n[k][3], n[k][1]).c_str(), surv(n[k][3], n[k][0]).c_str());
      else
        printf(" (PSD columns need sim_psd.py)\n");
    }
    tube.push_back({{n[3][0], n[3][1], n[3][2], n[3][3]}});
  }

  // background index: every ROI hit weighted by its section's activity.
  // radioassay and MC statistics are independent, so they add in quadrature
  printf("\n[g] background index, KS as built [cts/(keV kg yr)] = BI +- radioassay +- MC statistics:\n");
  printf("      %-6s %-6s %9s   %-34s %-34s\n", "chain", "sect", "ROI hits", "before cuts", "after cuts");
  // one section's ROI contribution: summed weight, summed weight^2 for Poisson, hits behind it,
  // and the mean weight of every decay there, which sets the limit when it saw nothing
  struct Sum { double w = 0, w2 = 0, wbar = 0; long n = 0; };
  auto cell = [&](const Sum &q, double rel, bool ul) -> std::string {
    if (q.n == 0) return Form("< %.2e", UL90 * q.wbar); // nothing seen: 90% CL on the MC alone
    return Form("%s%.2e +-%.1e +-%.1e", ul ? "<" : "", q.w, ul ? 0.0 : q.w * rel, sqrt(q.w2));
  };
  // a total: summed value, its two variances, how many hits back it, whether any part of it
  // rests on an upper-limit activity, and which sections saw nothing and are only bounded
  struct Tot { double v = 0, va = 0, vs = 0; long n = 0; bool ul = false; std::string blind; };
  auto add = [&](Tot &t, const Sum &q, double rel, bool ul, const char *name) {
    if (q.n == 0) { t.blind += (t.blind.empty() ? "" : ", ") + std::string(name); return; }
    t.v += q.w; t.va += pow(q.w * rel, 2); t.vs += q.w2; t.n += q.n; t.ul |= ul;
  };
  auto show = [&](const Tot &t, bool goal) -> std::string {
    if (t.n == 0) return "no ROI hit";
    std::string out = Form("%s%.2e +-%.1e +-%.1e", t.ul ? "<" : "", t.v, sqrt(t.va), sqrt(t.vs));
    if (goal) out += Form("  (%s%s x goal)", t.ul ? "<" : "", xg(t.v));
    return out;
  };
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
      for (Long64_t i = 0; i < r.nsim; i++)
      {
        if (slabOf(r.depth[i], L1ks, L2ks) != sec)
          continue;
        double w = weight(r, sec, r.pv[i]) / (2 * cfg::roiHalf); // per decay: the slab can mix volumes
        qb.wbar += w; nsec++;
        if (r.roiRaw[i]) { qb.w += r.roiRaw[i] * w; qb.w2 += r.roiRaw[i] * w * w; qb.n += r.roiRaw[i]; }
        if (r.roiCut[i]) { qa.w += r.roiCut[i] * w; qa.w2 += r.roiCut[i] * w * w; qa.n += r.roiCut[i]; }
      }
      qa.wbar = qb.wbar = nsec ? qb.wbar / nsec : 0;
      printf("      %-6s %-6s %4ld/%-4ld   %-34s %-34s%s\n", sec ? "" : r.iso.c_str(), secName[sec],
             qa.n, qb.n, cell(qb, rel, ul).c_str(), cell(qa, rel, ul).c_str(),
             ul ? "  activity is an upper limit" : "");
      std::string tag2 = r.iso + " " + secName[sec];
      add(tB, qb, rel, ul, secName[sec]); add(tA, qa, rel, ul, secName[sec]);
      add(gB, qb, rel, ul, tag2.c_str()); add(gA, qa, rel, ul, tag2.c_str());
    }
    printf("      %-6s %-6s %9s   %-34s %-34s\n", "", "tube", "", show(tB, false).c_str(), show(tA, false).c_str());
  }
  printf("      %-13s %9s   %-34s %-34s\n", "ALL CHAINS", "", show(gB, true).c_str(), show(gA, true).c_str());
  if (!gB.blind.empty())
    printf("      not in the totals, no ROI hit so only bounded (limits above): %s\n", gB.blind.c_str());
  if (gA.n == 0 && gB.n > 0)
    printf("      after cuts: none of the %ld ROI hits survived, so there is no after-cut value yet\n", gB.n);
  printf("      radioassay dominates wherever its term is larger than the MC one: more simulation\n"
         "      cannot shrink it, only a better activity measurement can\n");

//-------------------------------------------------------------------------------
//  7. Scan both boundaries, and draw:
  const int NS = 60;
  auto hScan = new TH2D("hScan", "ROI rate before cuts, all chains;L1 steel/Cu boundary [m];L2 Cu/EFCu boundary [m]",
                        NS, 0, depthTot / 1000, NS, 0, depthTot / 1000);
  for (int i = 1; i <= NS; i++)
    for (int j = 1; j <= NS; j++)
    {
      double L1 = hScan->GetXaxis()->GetBinCenter(i) * 1000, L2 = hScan->GetYaxis()->GetBinCenter(j) * 1000;
      if (L1 <= L2)
        hScan->SetBinContent(i, j, rate(L1, L2, ORD, false, nullptr, nullptr));
    }

  // spectra weighted to cts/(keV kg yr) with the as-built design, summed over runs
  auto hAll = new TH1D("hAll", "", cfg::eBins, cfg::eLo, cfg::eHi);
  auto hM1 = new TH1D("hM1", "", cfg::eBins, cfg::eLo, cfg::eHi);
  auto hCut = new TH1D("hCut", "", cfg::eBins, cfg::eLo, cfg::eHi);
  auto hPsd = new TH1D("hPsd", "", cfg::eBins, cfg::eLo, cfg::eHi);
  bool anyPsd = false;
  for (auto &r : R)
  {
    anyPsd |= r.hasPsd;
    for (auto &h : r.hits)
    {
      double w = weight(r, ORD[slabOf(r.depth[h.ev], L1ks, L2ks)], r.pv[h.ev]) / binW;
      hAll->Fill(h.e, w);
      if (h.m1) hM1->Fill(h.e, w);
      if (h.m1 && h.lar) hCut->Fill(h.e, w);
      if (h.m1 && h.lar && h.psd == 1) hPsd->Fill(h.e, w);
    }
  }

  gStyle->SetOptStat(0);
  auto c = new TCanvas("c_bkg", "", 1280, 560);
  c->Divide(2, 1);
  c->cd(1)->SetLogy();
  gPad->SetGrid();
  gPad->SetLeftMargin(0.15);
  hAll->SetTitle("RT background in the germanium;energy [keV];cts / (keV kg yr)");
  hAll->GetYaxis()->SetTitleOffset(1.6);
  struct { TH1D *h; int col; const char *lab; } sp[4] = {
      {hAll, kAzure + 2, "no cuts"}, {hM1, kOrange + 7, "M1"}, {hCut, kRed + 1, "M1 + argon veto"},
      {hPsd, kViolet + 1, "M1 + argon veto + PSD"}};
  auto leg = new TLegend(0.50, 0.72, 0.88, 0.88);
  leg->SetTextSize(0.032);
  for (int i = 0; i < (anyPsd ? 4 : 3); i++)
  {
    sp[i].h->SetLineColor(sp[i].col);
    sp[i].h->SetLineWidth(2);
    sp[i].h->Draw(i ? "HIST SAME" : "HIST");
    leg->AddEntry(sp[i].h, sp[i].lab, "l");
  }
  leg->Draw();
  c->cd(2)->SetLogz();
  gPad->SetRightMargin(0.17);
  hScan->Draw("COLZ");
  auto mk = new TMarker(L1ks / 1000, L2ks / 1000, 29);
  mk->SetMarkerColor(kRed + 1);
  mk->SetMarkerSize(2.2);
  mk->Draw();

  std::string tag;
  for (auto &r : R) tag += (tag.empty() ? "" : "_") + r.iso;
  std::string png = rtOut(tag + "_background.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s\n", png.c_str());

  // survival in the ROI against the RE-vessel references: the slide's chart, for the RT
  auto cs = new TCanvas("c_surv", "", 560 * R.size(), 480);
  cs->Divide(R.size(), 1);
  const char *cutName[4] = {"AC (M1+argon)", "PSD", "PSD | AC", "Combined"};
  for (size_t i = 0; i < R.size(); i++)
  {
    cs->cd(i + 1)->SetLogy();
    gPad->SetGrid(0, 1);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    int is = cfg::iso(R[i].iso);
    auto bars = new TH1D(Form("cdr%zu", i), Form("%s from the RT wall;;survival in ROI [%%]", R[i].iso.c_str()), 4, 0, 4);
    double ref[2][4];
    for (int k = 0; k < 3; k++) { ref[0][k] = cfg::cdr[is][k]; ref[1][k] = cfg::edgar[is][k]; }
    for (int j = 0; j < 2; j++) ref[j][3] = ref[j][0] * ref[j][2] / 100; // combined = AC x PSD|AC
    for (int k = 0; k < 4; k++)
    {
      bars->SetBinContent(k + 1, ref[0][k]);
      bars->GetXaxis()->SetBinLabel(k + 1, cutName[k]);
    }
    bars->SetFillColor(kAzure - 9);
    bars->SetLineColor(kAzure + 2);
    bars->SetBarWidth(0.5);
    bars->SetBarOffset(0.25);
    bars->SetMinimum(0.05);
    bars->SetMaximum(300);
    bars->GetXaxis()->SetLabelSize(0.05);
    bars->Draw("BAR");

    auto ge = new TGraph(); // Edgar's remage for the same component
    for (int k = 0; k < 4; k++) ge->SetPoint(k, k + 0.5, ref[1][k]);
    ge->SetMarkerStyle(24);
    ge->SetMarkerSize(1.6);
    ge->SetMarkerColor(kGray + 2);
    ge->Draw("P SAME");

    // this simulation: whole tube, binomial errors, each cut against its own denominator
    const long *n = tube[i].n;
    long pass[4] = {n[1], n[2], n[3], n[3]}, of[4] = {n[0], n[0], n[1], n[0]};
    auto gm = new TGraphErrors();
    for (int k = 0; k < 4; k++)
    {
      if (k > 0 && !R[i].hasPsd) break;       // PSD bars need sim_psd.py
      if (of[k] <= 0 || pass[k] <= 0) continue; // a log axis cannot show zero
      double f = (double)pass[k] / of[k];
      int m = gm->GetN();
      gm->SetPoint(m, k + 0.5, 100 * f);
      gm->SetPointError(m, 0, 100 * sqrt(f * (1 - f) / of[k]));
    }
    gm->SetMarkerStyle(20);
    gm->SetMarkerSize(1.4);
    gm->SetMarkerColor(kRed + 1);
    gm->SetLineColor(kRed + 1);
    gm->Draw("P SAME");

    auto lg = new TLegend(0.40, 0.76, 0.89, 0.89);
    lg->SetTextSize(0.035);
    lg->AddEntry(bars, "CDR - RE vessel (EFCu)", "f");
    lg->AddEntry(ge, "Edgar remage - RE Cu", "p");
    lg->AddEntry(gm, Form("this sim - KS tube (%ld ROI hits)", n[0]), "pe");
    lg->Draw();
  }
  std::string spng = rtOut(tag + "_survival.png");
  cs->SaveAs(spng.c_str());
  printf("wrote %s\n", spng.c_str());
}
