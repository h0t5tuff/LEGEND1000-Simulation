#include "TRegexp.h"
#include <fstream>
#include <sstream>
#include <set>
#include <unordered_map>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  // the ROI: one window around Qbb, Qbb +- 50 keV (1989-2089 keV), short of the Tl208 single-escape peak at 2103.5 keV.
  // the BI is its counts / 100 keV
  const double Qbb = 2039.0, roiHalf = 50.0;
  const double winLo = Qbb - roiHalf, winHi = Qbb + roiHalf;
  bool inWindow(double e) { return e >= winLo && e <= winHi; }
  const double winWidth = winHi - winLo; // 100 keV
  // the cuts of the LEGEND-1000 simulation production (legend1000-metadata simprod/config, l1000dsg01: tier/evt, pars/geds/psdcuts)
  const double m1_keV = 25.0;  // a detector fired above this; M1: exactly one did
  const double lar_pe = 4.0;   // the LAr veto rejects a decay whose SiPMs see this many photoelectrons, summed
  const double psd_low = -1.5; // the production's A/E cut: keep a hit whose A/E classifier is above this (low side)
  // LEGEND tunes it in data to keep 90% of the Tl208 double-escape peak (DEP); here that is only checked (checkAoE)
  const double depKeep = 0.90, peakHalf = 2.5, sideLo = 5.0, sideHi = 10.0; // keV: peak +-2.5, sidebands 5-10 off it
  const int NC = 4;                                                         // cut levels: none, M1, M1 + LAr, M1 + LAr + A/E
  const char *cutName[NC] = {"no cut", "M1", "M1 + LAr", "M1 + LAr + A/E"};
  const double bgGoal = 1e-5; // LEGEND-1000 goal [counts/(keV kg yr)] at Qbb, for scale
  const double simTol = 0.10; // a design passes when its BI after all cuts is at most this much above the tube's as built
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000; // spectrum [keV]
  const int eBins = 200;
  const int NDB = 10; // depth slabs of the DEPTH table

  const char *mat[3] = {"SS", "Cu", "EFCu"}; // top to bottom; their densities are the GDML's (output/tube.txt)
  // chain activities [uBq/kg], {232Th, 238U}: steel from "RT MaterialMix, ralph.pdf" (p. 5), OFHC Cu from MAJORANA's assay
  // paper (Phys. Rev. C 110, 055804, Table II), EFCu from the LEGEND-1000 CDR's Table XVII
  const double chain[3][2] = {{1000, 2500}, {1.10, 1.37}, {0.37, 0.19}};
  const double branch[2] = {0.3594, 1.0};                                                           // Tl208 is 35.94% of the 232Th chain (Nuclear Data Sheets for A = 208); Bi214 100% of the 238U chain
  int iso(const std::string &s) { return s.find("Bi") != std::string::npos ? 1 : 0; }               // 0 Tl208, 1 Bi214
  double activity(const std::string &s, int m) { return chain[m][iso(s)] * branch[iso(s)] * 1e-6; } // Bq/kg of the simulated nuclide
}

//-------------------------------------------------------------------------------
//  2. The tube, from output/tube.txt (geom/l1000geom/l1000-autopeel.py --table reads it off the GDML): z in mm
struct Tube
{
  double zTop = 0, zBottom = 0, seamCu = 0, seamSS = 0, density[3] = {0, 0, 0}, geMass = 0; // density: SS, Cu, EFCu [kg/m^3]
  int nGe = 0;
  std::vector<double> area;                                                                  // the wall's cross-section [mm^2] per 0.5 mm, from the top down
  int sectionAt(double z) const { return z >= seamSS ? 0 : (z >= seamCu ? 1 : 2); }          // 0 SS, 1 Cu, 2 EFCu
  double volume(double d0, double d1) const                                                  // [m^3] the wall between two depths below the top [mm]
  {
    double v = 0;
    for (long i = std::max(0L, lround(2 * d0)); i < std::min((long)area.size(), lround(2 * d1)); i++)
      v += 0.5 * area[i];
    return v * 1e-9;
  }
};

static Tube readTube(const char *file)
{
  Tube t;
  std::ifstream in(file);
  in >> t.zTop >> t.zBottom >> t.seamCu >> t.seamSS >> t.density[0] >> t.density[1] >> t.density[2] >> t.geMass >> t.nGe;
  for (double a; in >> a;)
    t.area.push_back(a);
  return t;
}

template <class F>
static void scanTree(TTree *t, const std::vector<const char *> &cols, F f) // the named columns as doubles, row by row
{
  std::vector<TTreeFormula *> col;
  for (auto c : cols)
    col.push_back(new TTreeFormula(c, c, t));
  std::vector<double> v(cols.size());
  for (Long64_t i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    for (size_t k = 0; k < col.size(); k++)
      v[k] = col[k]->EvalInstance();
    f(v.data());
  }
  for (auto c : col)
    delete c;
}

//-------------------------------------------------------------------------------
//  3. The runs of one isotope, reduced to what a reweighting needs:
struct Hit // one detector's energy in one decay; ev: the decay's slot; lar: the decay passes the LAr veto, larMap: the same with the map alone
{
  Long64_t ev;
  float e;
  bool m1, lar, larMap, psd;
  bool pass(int c) const { return c == 0 || (m1 && (c == 1 || (lar && (c == 2 || psd)))); }
  bool passMap(int k) const { return m1 && larMap && (k == 0 || psd); } // M1 + LAr (k 0), + A/E (k 1), the LAr veto off the map alone
};

struct Run
{
  std::string iso;
  Long64_t nsim = 0;
  long nPv[3] = {0, 0, 0};                                            // decays per section: SS, Cu, EFCu
  double mcDensity[3] = {0, 0, 0};                                    // simulated decays per m^3 in each section, measured
  long dN[cfg::NDB] = {0}, dW1[cfg::NDB] = {0}, dWA[cfg::NDB] = {0}; // per depth slab: decays, ROI counts after M1 and after all cuts
  double S = 0;                                                       // the cut survival near the detectors (DESIGNS)
  std::vector<float> depth;                                           // per decay WITH a Ge hit: mm below the top of the tube
  std::vector<signed char> pv;                                        // its section
  std::vector<short> win[cfg::NC];                                    // its hits in the ROI, per cut level
  std::vector<short> winMap[2];                                       // ...and with the LAr veto read off the map alone: M1 + LAr, + A/E
  std::vector<size_t> roi;                                            // the decays with an ROI count after M1: the only ones a design changes
  std::vector<Hit> hits;
};

static int slabOf(double dep, double L1, double L2) { return dep < L1 ? 0 : (dep < L2 ? 1 : 2); } // a design's three slabs, from the top
static int bin(double dep, double tot, int n) { return std::min(n - 1, std::max(0, (int)(dep / tot * n))); }

static void readFile(const std::string &file, Run &r, const Tube &t) // one hit file, added to its isotope's run
{
  auto f = TFile::Open(file.c_str());
  auto vtx = (f && !f->IsZombie()) ? (TTree *)f->Get("vtx") : nullptr;
  auto tg = vtx ? (TTree *)f->Get("geds") : nullptr, tl = vtx ? (TTree *)f->Get("lar") : nullptr;
  if (!tg || !tl)
  {
    printf("ERROR: no vtx, geds or lar tree in %s\n", file.c_str());
    delete f;
    return;
  }
  r.nsim += vtx->GetEntries();
  const Long64_t n = (Long64_t)vtx->GetMaximum("evtid") + 1; // the event numbers, up to it: some may have no vertex, no decay
  std::vector<float> depth(n, -1), pe(n, 0), peMap(n, 0), aoe;
  std::vector<signed char> pv(n, -1);
  std::vector<short> nFired(n, 0);
  std::vector<Hit> hits;
  scanTree(tg, {"evtid", "energy", "aoe_class"}, [&](const double *v) { // one row per detector hit, the response applied
    if (v[0] < 0 || v[0] >= n || v[1] <= cfg::m1_keV)
      return;
    nFired[(int)v[0]]++;
    hits.push_back({(Long64_t)v[0], (float)v[1], false, false, false, false});
    aoe.push_back(v[2]);
  });
  scanTree(tl, {"evtid", "pe", "pe_map"}, [&](const double *v)
           { if (v[0] >= 0 && v[0] < n) { pe[(int)v[0]] += v[1]; peMap[(int)v[0]] += v[2]; } });
  scanTree(vtx, {"evtid", "zloc"}, [&](const double *v) { // [m]; rows are not in event order with -t: index by evtid
    if (v[0] >= 0 && v[0] < n)
      depth[(int)v[0]] = t.zTop - 1000 * v[1], pv[(int)v[0]] = t.sectionAt(1000 * v[1]);
  });
  std::vector<std::vector<short>> win(cfg::NC, std::vector<short>(n, 0)), winMap(2, std::vector<short>(n, 0));
  for (size_t i = 0; i < hits.size(); i++)
  {
    Hit &h = hits[i];
    h.m1 = nFired[h.ev] == 1;
    h.lar = pe[h.ev] < cfg::lar_pe;
    h.larMap = peMap[h.ev] < cfg::lar_pe;
    h.psd = aoe[i] > cfg::psd_low; // an undefined A/E (NaN) fails
    for (int c = 0; c < cfg::NC; c++)
      win[c][h.ev] += h.pass(c) && cfg::inWindow(h.e);
    for (int k = 0; k < 2; k++)
      winMap[k][h.ev] += h.passMap(k) && cfg::inWindow(h.e);
  }
  const double tot = t.zTop - t.zBottom;
  for (Long64_t i = 0; i < n; i++) // every decay
    if (depth[i] >= 0)
    {
      const int b = bin(depth[i], tot, cfg::NDB);
      r.dN[b]++, r.dW1[b] += win[1][i], r.dWA[b] += win[cfg::NC - 1][i], r.nPv[pv[i]]++;
    }
  std::vector<Long64_t> slot(n, -1); // a decay with a hit gets one slot in the run, and its hits point there
  for (auto h : hits)
  {
    if (slot[h.ev] < 0)
    {
      slot[h.ev] = r.depth.size();
      r.depth.push_back(depth[h.ev]);
      r.pv.push_back(pv[h.ev]);
      for (int c = 0; c < cfg::NC; c++)
        r.win[c].push_back(win[c][h.ev]);
      for (int k = 0; k < 2; k++)
        r.winMap[k].push_back(winMap[k][h.ev]);
    }
    h.ev = slot[h.ev];
    r.hits.push_back(h);
  }
  delete f;
}

static std::vector<std::string> expand(const std::string &pat) // "output/tl208*.root" -> every match, sorted
{
  if (pat.find_first_of("*?") == std::string::npos)
    return {pat};
  size_t sl = pat.find_last_of('/');
  std::string dir = sl == std::string::npos ? "." : pat.substr(0, sl), base = pat.substr(sl + 1);
  TRegexp re(base.c_str(), kTRUE); // a shell wildcard
  std::vector<std::string> out;
  void *d = gSystem->OpenDirectory(dir.c_str());
  for (const char *e; d && (e = gSystem->GetDirEntry(d));)
  {
    TString f(e);
    Ssiz_t len = 0;
    if (f.Index(re, &len) == 0 && len == f.Length())
      out.push_back(sl == std::string::npos ? f.Data() : dir + "/" + f.Data());
  }
  if (d)
    gSystem->FreeDirectory(d);
  std::sort(out.begin(), out.end());
  if (out.empty())
    printf("no file matches %s\n", pat.c_str());
  return out;
}

//-------------------------------------------------------------------------------
//  4. The A/E cut, checked on the Tl208 runs' M1 hits against the LEGEND-1000 pCDR (p. 71): what it keeps of the single-site
//     DEP (net of the continuum, from the sidebands; LEGEND tunes the cut to 90%) and of the ROI (the pCDR: reduced 2-3 x)
static void checkAoE(const std::vector<std::string> &files)
{
  const double dep = 1592.5;
  std::vector<float> peak, side, roi;
  for (auto &file : files)
  {
    auto f = TFile::Open(file.c_str());
    auto t = (f && !f->IsZombie()) ? (TTree *)f->Get("geds") : nullptr;
    std::unordered_map<Long64_t, short> nFired;
    if (t)
      scanTree(t, {"evtid", "energy"}, [&](const double *v)
               { if (v[1] > cfg::m1_keV) nFired[(Long64_t)v[0]]++; });
    if (t)
      scanTree(t, {"evtid", "energy", "aoe_class"}, [&](const double *v)
               {
        if (v[1] <= cfg::m1_keV || nFired[(Long64_t)v[0]] != 1) return;
        float a = std::isnan(v[2]) ? -1e30 : v[2]; // an undefined A/E fails
        if (cfg::inWindow(v[1])) roi.push_back(a);
        double d = fabs(v[1] - dep);
        if (d < cfg::peakHalf) peak.push_back(a);
        else if (d > cfg::sideLo && d < cfg::sideHi) side.push_back(a); });
    delete f;
  }
  const double w = cfg::peakHalf / (cfg::sideHi - cfg::sideLo); // sideband counts per peak width
  auto above = [](const std::vector<float> &x)
  { return (double)std::count_if(x.begin(), x.end(), [](float a) { return a > cfg::psd_low; }); };
  const double net = peak.size() - w * side.size(), roiKept = roi.empty() ? 0.0 : above(roi) / roi.size();
  printf("A/E CUT  classifier > %.2f (the production's): it keeps %.2f of the DEP (+- %.2f, %.0f net; LEGEND tunes it to %.2f)\n",
         cfg::psd_low, (above(peak) - w * above(side)) / net, sqrt(peak.size() + w * w * side.size()) / net, net, cfg::depKeep);
  printf("    and %.2f of the ROI: a factor %.1f (the pCDR: 2-3)\n\n", roiKept, roiKept > 0 ? 1 / roiKept : 0.0);
}

void ana(const char *runs = "output/output_root/tl208*_hit.root=Tl208,output/output_root/bi214*_hit.root=Bi214", const char *tubeFile = "output/tube.txt")
{
  gErrorIgnoreLevel = kWarning; // no "png file has been created" notices in the printout
  const Tube rt = readTube(tubeFile);
  if (rt.area.empty())
  {
    printf("ERROR: no %s: run geom/l1000geom/l1000-autopeel.py --table %s first\n", tubeFile, tubeFile);
    return;
  }
  std::vector<std::pair<std::string, std::string>> pats; // "a_hit.root=Tl208,b_hit.root=Bi214": (pattern, isotope)
  std::stringstream list(runs);
  for (std::string tok; std::getline(list, tok, ',');)
  {
    size_t eq = tok.find('=');
    pats.push_back({tok.substr(0, eq), eq == std::string::npos ? "Tl208" : tok.substr(eq + 1)});
  }
  std::vector<std::string> tlFiles; // every Tl208 file, for the A/E check
  for (auto &[pat, iso] : pats)
    if (cfg::iso(iso) == 0)
      for (auto &f : expand(pat))
        tlFiles.push_back(f);
  checkAoE(tlFiles);
  std::map<int, Run> R; // one per isotope: 0 Tl208, 1 Bi214
  for (auto &[pat, iso] : pats)
    for (auto &f : expand(pat))
      R[cfg::iso(iso)].iso = iso, readFile(f, R[cfg::iso(iso)], rt);

  const double depthTot = rt.zTop - rt.zBottom, L1built = rt.zTop - rt.seamSS, L2built = rt.zTop - rt.seamCu; // the seams as built, as depths
  const double binW = (cfg::eHi - cfg::eLo) / cfg::eBins;
  const int NC = cfg::NC, ALL = NC - 1;
  const double edge[4] = {0, L1built, L2built, depthTot};
  double physV[3]; // m^3: SS, Cu, EFCu, as built
  for (int p = 0; p < 3; p++)
    physV[p] = rt.volume(edge[p], edge[p + 1]);
  for (auto &[key, r] : R)
    for (int p = 0; p < 3; p++)
      r.mcDensity[p] = r.nPv[p] / physV[p]; // remage fills the mother sparser than its daughters: measure, never assume
  // one simulated decay in section p, where material m sits = rho x A x 1 yr real decays per m^3, over the MC density, per kg of Ge
  auto weight = [&](const Run &r, int m, int p)
  { return p < 0 || r.mcDensity[p] <= 0 ? 0.0 : rt.density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[p] * rt.geMass); };
  auto massOf = [&](double d0, double d1, int m)
  { return d1 > d0 ? rt.density[m] * rt.volume(d0, d1) : 0.0; };
  printf("INPUTS  ROI %.0f-%.0f keV (%.0f keV); cuts: M1 (one detector above %.0f keV), LAr veto (< %.0f PE), A/E > %.2f; germanium %.1f kg\n\n",
         cfg::winLo, cfg::winHi, cfg::winWidth, cfg::m1_keV, cfg::lar_pe, cfg::psd_low, rt.geMass);

  //-------------------------------------------------------------------------------
  //  5. The background index as built, per chain and section: the README's as-built table
  printf("BI AS BUILT  [counts/(keV kg yr)], +- simulation statistics\n");
  printf("    %-6s %-5s %8s %10s %10s %9s  %-30s %10s  %s\n", "chain", "sect", "mass", "real/yr", "simulated", "w", "ROI counts: none/M1/+LAr/+A/E",
         "counts/yr", "BI, all cuts");
  double gB[NC] = {0}, gV[NC] = {0}, mB[2] = {0}, mV[2] = {0}; // all chains, per cut level; and with the LAr veto off the map alone
  for (auto &[key, r] : R)
    for (int s = 0; s < 3; s++)
    {
      const double wt = weight(r, s, s) / cfg::winWidth, real = physV[s] * rt.density[s] * cfg::activity(r.iso, s) * cfg::secPerYear;
      long n[NC] = {0};
      for (size_t i = 0; i < r.depth.size(); i++)
        if (r.pv[i] == s)
        {
          for (int c = 0; c < NC; c++)
            n[c] += r.win[c][i];
          for (int k = 0; k < 2; k++)
            mB[k] += r.winMap[k][i] * wt, mV[k] += r.winMap[k][i] * wt * wt;
        }
      for (int c = 0; c < NC; c++)
        gB[c] += n[c] * wt, gV[c] += n[c] * wt * wt;
      printf("    %-6s %-5s %6.1f kg %10.4g %10ld %9.2g  %-30s %10.2g  %.2e +- %.1e\n", r.iso.c_str(), cfg::mat[s], physV[s] * rt.density[s], real,
             r.nPv[s], r.nPv[s] ? real / r.nPv[s] : 0.0, Form("%ld/%ld/%ld/%ld", n[0], n[1], n[2], n[3]), n[ALL] * wt * cfg::winWidth * rt.geMass,
             n[ALL] * wt, sqrt((double)n[ALL]) * wt);
    }
  printf("    all chains:");
  for (int c = 0; c < NC; c++)
    printf("  %s %.2e +- %.1e;", cfg::cutName[c], gB[c], sqrt(gV[c]));
  printf("\n    with the LAr veto off the optical map alone (no light from the argon beyond its edge): M1 + LAr %.2e, all cuts %.2e\n\n", mB[0], mB[1]);

  //-------------------------------------------------------------------------------
  //  6. Where the ROI counts come from, by depth: the README's depth table
  for (auto &[key, r] : R)
  {
    printf("DEPTH  %s  per slab below the tube's top: decays, ROI counts after M1, after all cuts\n", r.iso.c_str());
    for (int b = 0; b < cfg::NDB; b++)
      printf("    %5.2f-%-5.2f m %12ld %8ld %6ld\n", depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000, r.dN[b], r.dW1[b], r.dWA[b]);
  }

  //-------------------------------------------------------------------------------
  //  7. Designs against the tube as built: steel to L1, Cu to L2, EFCu below, for less EFCu (L2 deeper) and more steel (L1 deeper)
  // a design differs from the tube as built only where a slab changed material, so its BI after all cuts is the as-built's plus
  // that change, dB. Too few counts survive all cuts away from the detectors to count dB on them, so it is counted on the
  // changed slabs' ROI counts after M1 times the cut survival (the share of them all cuts keep, see DEPTH), measured where the
  // change is: near the detectors, below farEdge, per chain; above it, both chains together
  const double biBuilt = gB[ALL], maxRise = cfg::simTol * biBuilt;
  const int nFar = 7; // the slabs above 4.37 m: above the array
  const double farEdge = depthTot * nFar / cfg::NDB;
  long fA = 0, fM = 0;
  for (auto &[key, r] : R)
  {
    long a = 0, m = 0;
    for (int b = 0; b < cfg::NDB; b++)
      (b < nFar ? fA : a) += r.dWA[b], (b < nFar ? fM : m) += r.dW1[b];
    r.S = m ? (double)a / m : 0;
  }
  const double sFar = fM ? (double)fA / fM : 0;
  for (auto &[key, r] : R)
    for (size_t i = 0; i < r.depth.size(); i++)
      if (r.win[1][i])
        r.roi.push_back(i);
  auto change = [&](double L1, double L2)
  {
    double dB = 0;
    for (auto &[key, r] : R)
      for (size_t i : r.roi)
      {
        const int mb = slabOf(r.depth[i], L1built, L2built), md = slabOf(r.depth[i], L1, L2);
        if (mb != md)
          dB += (r.depth[i] < farEdge ? sFar : r.S) * r.win[1][i] * (weight(r, md, r.pv[i]) - weight(r, mb, r.pv[i])) / cfg::winWidth;
      }
    return dB;
  };
  auto deepestSteel = [&](double L2) { // 1 cm at a time below the steel as built, always some Cu
    double l1 = -1;
    for (double L1 = L1built; L1 <= L2 - 10 && change(L1, L2) <= maxRise; L1 += 10)
      l1 = L1;
    return l1;
  };
  auto split = [&](double l1, double L2)
  { return std::string(Form("%4.0f : %3.0f : %3.0f kg", massOf(0, l1, 0), massOf(l1, L2, 1), massOf(L2, depthTot, 2))); };
  printf("\nDESIGNS  SS down to L1, Cu down to L2, EFCu below; a design passes when its BI stays within %.0f%% of the as-built's (%.2e)\n",
         100 * cfg::simTol, biBuilt);
  printf("    cut survival: below %.2f m", farEdge / 1000);
  for (auto &[key, r] : R)
    printf(" %s %.2e", r.iso.c_str(), r.S);
  printf("; above it, both chains, %ld of %ld = %.2e\n", fA, fM, sFar);
  printf("    %-6s %6s | %-10s %-24s %9s\n", "L2 [m]", "EFCu", "deepest L1", "SS : Cu : EFCu", "BI change");
  std::vector<double> rows = {L2built};
  for (double L2 = 4100; L2 <= 5000; L2 += 100)
    if (L2 > L2built)
      rows.push_back(L2);
  for (double L2 = 5250; L2 < depthTot; L2 += 250)
    rows.push_back(L2);
  double best[2] = {-1, -1}; // L1, L2 of the least EFCu, then the most steel
  for (double L2 : rows)
  {
    const double l1 = deepestSteel(L2);
    if (l1 < 0)
    {
      printf("    %-6.2f %3.0f kg | the BI rises more than 10%% at any L1\n", L2 / 1000, massOf(L2, depthTot, 2));
      continue;
    }
    printf("    %-6.2f %3.0f kg | %8.2f m %-24s %+8.0f%%\n", L2 / 1000, massOf(L2, depthTot, 2), l1 / 1000, split(l1, L2).c_str(), 100 * change(l1, L2) / biBuilt);
    if (L2 > best[1] || (L2 == best[1] && l1 > best[0]))
      best[0] = l1, best[1] = L2;
  }
  printf("    best: L2 %.2f m, L1 %.2f m, %s, BI %+.0f%%; as built: %s\n", best[1] / 1000, best[0] / 1000, split(best[0], best[1]).c_str(),
         100 * change(best[0], best[1]) / biBuilt, split(L1built, L2built).c_str());

  //-------------------------------------------------------------------------------
  //  8. The tube's spectrum as built, per cut level, for the slide below
  TH1D *sp[NC];
  for (int c = 0; c < NC; c++)
    sp[c] = new TH1D(Form("sp%d", c), "Re-entrant tube background, as built;energy [keV];counts / (keV kg yr)", cfg::eBins, cfg::eLo, cfg::eHi);
  for (auto &[key, r] : R)
    for (auto &h : r.hits)
    {
      double wt = weight(r, slabOf(r.depth[h.ev], L1built, L2built), r.pv[h.ev]) / binW;
      for (int c = 0; c < NC; c++)
        if (h.pass(c))
          sp[c]->Fill(h.e, wt);
    }
  gStyle->SetOptStat(0);
  const int col[NC] = {kAzure + 2, kOrange + 7, kRed + 1, kViolet + 1};
  for (int k = 0; k < NC; k++)
    sp[k]->SetLineColor(col[k]);
  sp[0]->SetMinimum(0.1 * std::max(1e-12, sp[ALL]->GetMinimum(0) > 0 ? sp[ALL]->GetMinimum(0) : sp[2]->GetMinimum(0)));

  //-------------------------------------------------------------------------------
  //  9. Two slide figures, one message each: the spectrum with the ROI and the BI marked; the design map near the tube as built
  auto sci = [](double x) { // "9.8 #times 10^{-9}"
    const int e = (int)floor(log10(x));
    return std::string(Form("%.1f #times 10^{%d}", x / pow(10, e), e));
  };
  auto line = [](double x0, double y0, double x1, double y1, int col, int width, int style = 1) {
    auto l = new TLine(x0, y0, x1, y1);
    l->SetLineColor(col);
    l->SetLineWidth(width);
    l->SetLineStyle(style);
    l->Draw();
  };
  auto label = [](double x, double y, const std::string &s, int col, double size, int align = 12) {
    auto t = new TLatex(x, y, s.c_str());
    t->SetTextColor(col);
    t->SetTextSize(size);
    t->SetTextFont(42);
    t->SetTextAlign(align);
    t->Draw();
  };
  { // the spectrum
    auto cs = new TCanvas("c_slide_sp", "", 1300, 760);
    cs->SetLogy();
    cs->SetGridy();
    cs->SetMargin(0.09, 0.03, 0.11, 0.08);
    const char *what[NC] = {"no cut", "M1: one detector only", "+ LAr veto", "+ A/E: all cuts"};
    sp[0]->SetTitle("What the as-built tube leaves in the detectors;energy [keV];counts / (keV kg yr)");
    sp[0]->SetMaximum(5 * sp[0]->GetMaximum());
    sp[0]->Draw("HIST");
    gPad->Update();
    const double y0 = pow(10, gPad->GetUymin()), y1 = pow(10, gPad->GetUymax());
    const int roiCol = TColor::GetColor("#2e7d32");
    auto b = new TBox(cfg::winLo, y0, cfg::winHi, y1); // the ROI, one band
    b->SetFillColorAlpha(roiCol, 0.2);
    b->SetLineWidth(0);
    b->Draw();
    line(cfg::Qbb, y0, cfg::Qbb, y1, roiCol, 2);
    auto leg = new TLegend(0.10, 0.80, 0.50, 0.91);
    leg->SetNColumns(2);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.029);
    for (int k = 0; k < NC; k++)
    {
      sp[k]->SetLineWidth(3);
      sp[k]->Draw("HIST SAME");
      leg->AddEntry(sp[k], what[k], "l");
    }
    leg->Draw();
    line(cfg::eLo, cfg::bgGoal, cfg::eHi, cfg::bgGoal, kGray + 2, 2, 2); // the goal, dashed
    label(cfg::eHi - 20, 1.6 * cfg::bgGoal, "LEGEND-1000 goal: 10^{-5}", kGray + 2, 0.034, 31);
    line(cfg::winLo, biBuilt, cfg::winHi, biBuilt, kBlack, 4); // the BI over the ROI
    label(cfg::winHi + 15, biBuilt / 3, "#splitline{BI = " + sci(biBuilt) + "}{purple, averaged over the ROI}", kBlack, 0.032, 13); // below the spectra
    label(cfg::winHi + 6, pow(10, gPad->GetUymax() - 0.06 * (gPad->GetUymax() - gPad->GetUymin())), Form("ROI: Q_{#beta#beta} #pm %.0f keV", cfg::roiHalf), roiCol, 0.036);
    label(2614.5 - 12, 0.6 * sp[0]->GetBinContent(sp[0]->FindBin(2614.5)), "Tl-208, 2614 keV", kGray + 3, 0.032, 32);
    cs->SaveAs("output/slide_spectrum.png");
  }
  { // the design map, near the tube as built
    auto hz = new TH2D("hZoom", ";steel ends at depth L1 [m];EFCu starts at depth L2 [m];BI / BI as built", 48, 1.6, 3.2, 36, 3.9, 4.8);
    for (int i = 1; i <= 48; i++)
      for (int j = 1; j <= 36; j++)
        hz->SetBinContent(i, j, 1 + change(hz->GetXaxis()->GetBinCenter(i) * 1000, hz->GetYaxis()->GetBinCenter(j) * 1000) / biBuilt);
    auto cd = new TCanvas("c_slide_designs", "", 1100, 800);
    cd->SetMargin(0.10, 0.15, 0.11, 0.05);
    gStyle->SetPalette(kBird);
    hz->SetMinimum(0.9);
    hz->SetMaximum(1.3);
    hz->SetContour(40);
    hz->GetZaxis()->SetTitleOffset(1.3);
    hz->Draw("COLZ");
    auto leg = new TLegend(0.12, 0.80, 0.35, 0.92);
    for (auto [x, y, col, name] : {std::tuple{L1built, L2built, kRed + 1, "as built"}, std::tuple{best[0], best[1], (int)kBlack, "best design"}})
    {
      auto m = new TMarker(x / 1000, y / 1000, 29); // a star
      m->SetMarkerColor(col);
      m->SetMarkerSize(3.5);
      m->Draw();
      leg->AddEntry(m, name, "p");
    }
    leg->Draw();
    cd->SaveAs("output/slide_designs.png");
  }
  { // 10. Every design for LEGEND.app's rt (app/): a 2 cm grid through the seams as built, the BI over the as-built one,
    //     with the masses -> output/designs.json
    std::vector<double> g1, g2; // L1, L2 [mm]
    for (double v = fmod(L1built, 20); v <= depthTot; v += 20)
      g1.push_back(v);
    for (double v = fmod(L2built, 20); v <= depthTot; v += 20)
      g2.push_back(v);
    auto arr = [](size_t n, auto f) { // "[f(0),f(1),...]"
      std::string s = "[";
      for (size_t i = 0; i < n; i++)
        s += (i ? "," : "") + f(i);
      return s + "]";
    };
    auto num = [](const char *fmt, double x) { return std::string(Form(fmt, x)); };
    std::ofstream("output/designs.json")
        << "{\"depthTot\":" << depthTot / 1000 << ",\"L1built\":" << L1built / 1000 << ",\"L2built\":" << L2built / 1000
        << ",\"biBuilt\":" << num("%.4g", biBuilt) << ",\"tol\":" << cfg::simTol
        << ",\"l1\":" << arr(g1.size(), [&](size_t i) { return num("%.3f", g1[i] / 1000); })
        << ",\"l2\":" << arr(g2.size(), [&](size_t j) { return num("%.3f", g2[j] / 1000); })
        << ",\"mSS\":" << arr(g1.size(), [&](size_t i) { return num("%.1f", massOf(0, g1[i], 0)); })
        << ",\"mCu1\":" << arr(g1.size(), [&](size_t i) { return num("%.1f", massOf(0, g1[i], 1)); })
        << ",\"mCu2\":" << arr(g2.size(), [&](size_t j) { return num("%.1f", massOf(0, g2[j], 1)); })
        << ",\"mEF\":" << arr(g2.size(), [&](size_t j) { return num("%.1f", massOf(g2[j], depthTot, 2)); })
        << ",\"ratio\":" << arr(g1.size(), [&](size_t i) { return arr(g2.size(), [&](size_t j) { // null: no Cu left
             return g1[i] >= g2[j] ? std::string("null") : num("%.4g", (biBuilt + change(g1[i], g2[j])) / biBuilt); }); })
        << "}";
  }
  printf("\nwrote output/slide_spectrum.png, output/slide_designs.png and output/designs.json\n");
}

//-------------------------------------------------------------------------------
//  11. Which decays can reach Qbb (README step 2): sim/sim.mac's whole Th-232 and U-238 chains in the EFCu detector holders, no detector
//      response. Each detector's steps within 10 us of its first are one hit; the spectrum of the hits, its lines named.
//      root -l -b -q -e '.L ana/ana.C' -e 'holders()'
static Long64_t holderHits(const char *file, TH1D *h) // fills h, returns the chains
{
  TFile f(file);
  auto stp = f.IsZombie() ? nullptr : f.GetDirectory("stp");
  auto vtx = stp ? (TTree *)stp->Get("vtx") : nullptr;
  if (!vtx)
  {
    printf("ERROR: no stp/vtx in %s\n", file);
    return 0;
  }
  long nHit = 0;
  std::set<std::string> done; // a tree written in parts has several keys
  for (auto k : *stp->GetListOfKeys())
  {
    const std::string det = k->GetName();
    if (det[0] != 'V' || !done.insert(det).second) // the germanium detectors: V00101Z ...
      continue;
    std::vector<std::tuple<double, double, double>> s; // evtid, time [ns], energy [keV]
    scanTree((TTree *)stp->Get(det.c_str()), {"evtid", "time_in_ns", "edep_in_keV"}, [&](const double *v)
             { s.push_back({v[0], v[1], v[2]}); });
    std::sort(s.begin(), s.end()); // the threads write their events interleaved
    double ev = -1, t0 = 0, sum = 0;
    for (auto [e, t, E] : s)
    {
      if (e != ev || t - t0 > 1e4)
      {
        if (ev >= 0)
          h->Fill(sum), nHit++;
        ev = e, t0 = t, sum = 0;
      }
      sum += E;
    }
    if (ev >= 0)
      h->Fill(sum), nHit++;
  }
  printf("%s: %lld chains, %ld hits\n", file, vtx->GetEntries(), nHit);
  return vtx->GetEntries();
}

void holders()
{
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  const char *file[2] = {"output/holders_th232.root", "output/holders_u238.root"};
  const char *chain[2] = {"Th-232", "U-238"}, *fill[2] = {"#2f74c8", "#e8862a"}, *ink[2] = {"#1d4f91", "#b35900"};
  struct Line // a line, its label centre [keV] and foot [log fraction of the frame]
  {
    double e;
    const char *keV, *name;
    int ch;
    double xLab, yLab;
  };
  const Line lines[] = {{238.6, "238", "Pb-212", 0, 200, 0.90}, {351.9, "352", "Pb-214", 1, 390, 0.79}, {583.2, "583", "Tl-208", 0, 590, 0.90},
                        {609.3, "609", "Bi-214", 1, 790, 0.79}, {911.2, "911", "Ac-228", 0, 1000, 0.90}, {1120.3, "1120", "Bi-214", 1, 1220, 0.72},
                        {1764.5, "1764", "Bi-214", 1, 1764, 0.70}, {2447.9, "2448", "Bi-214", 1, 2448, 0.50}, {2614.5, "2614", "Tl-208", 0, 2614, 0.74}};
  const double eMax = 3000, band = 2 * sqrt(0.5 + 0.001 * cfg::Qbb) / 2.35482; // Qbb +- 2 sigma, the tube's resolution (step 3)
  TH1D *h[2];
  Long64_t nCh[2] = {0, 0};
  for (int k = 0; k < 2; k++)
  {
    h[k] = new TH1D(Form("hh%d", k), ";energy per detector per event [keV];counts / 5 keV", 600, 0, eMax);
    nCh[k] = holderHits(file[k], h[k]);
  }
  if (!nCh[0] && !nCh[1])
    return;

  auto c = new TCanvas("c_holders", "", 1400, 760);
  c->SetMargin(0.09, 0.02, 0.10, 0.08);
  c->SetLogy();
  const double yMin = 0.5, yMax = 50 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum());
  auto ndc = [&](double x, double frac, const char *txt) { // text in the pad's own coordinates: subscripts misplace on a log axis
    auto t = new TLatex(c->GetLeftMargin() + x / eMax * (1 - c->GetLeftMargin() - c->GetRightMargin()),
                        c->GetBottomMargin() + frac * (1 - c->GetBottomMargin() - c->GetTopMargin()), txt);
    t->SetNDC();
    t->SetTextFont(62);
    t->SetTextSize(0.026);
    return t;
  };
  h[0]->SetMinimum(yMin);
  h[0]->SetMaximum(yMax);
  h[0]->GetYaxis()->SetTitleOffset(1.1);
  h[0]->Draw("axis");
  auto leg = new TLegend(0.79, 0.82, 0.975, 0.91);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  for (int k = 0; k < 2; k++)
    if (nCh[k])
    {
      const int col = TColor::GetColor(fill[k]);
      h[k]->SetFillColorAlpha(col, 0.55);
      h[k]->SetLineColor(col);
      h[k]->Draw("hist same");
      leg->AddEntry(h[k], Form("%s chain (%.0f k)", chain[k], nCh[k] / 1e3), "f");
    }
  const int roiCol = TColor::GetColor("#2e7d32"); // the Qbb +- 2 sigma band: 2.7 keV, a line at this scale
  auto b = new TBox(cfg::Qbb - band, yMin, cfg::Qbb + band, yMax);
  b->SetFillColorAlpha(roiCol, 0.8);
  b->SetLineWidth(0);
  b->Draw();
  auto roiText = ndc(cfg::Qbb + 15, 0.62, Form("#splitline{Q_{#beta#beta} #pm 2#sigma}{%.1f-%.1f keV}", cfg::Qbb - band, cfg::Qbb + band));
  roiText->SetTextAlign(11);
  roiText->SetTextColor(roiCol);
  roiText->Draw();
  gPad->RedrawAxis();
  leg->Draw();
  auto title = new TLatex(0.5, 0.95, "Radioactivity of the EFCu detector holders: Th-232 and U-238 chains");
  title->SetNDC();
  title->SetTextAlign(22);
  title->SetTextSize(0.04);
  title->Draw();
  for (const auto &L : lines) // each line, in its chain's colour
  {
    const int col = TColor::GetColor(ink[L.ch]);
    const double top = h[L.ch]->GetBinContent(h[L.ch]->FindBin(L.e)), yLab = yMin * pow(yMax / yMin, L.yLab);
    auto a = new TArrow(L.xLab, yLab / 1.35, L.e, top * 1.25, 0.008, "|>");
    a->SetLineColor(col);
    a->SetFillColor(col);
    a->SetLineWidth(2);
    a->Draw();
    auto t = ndc(L.xLab, L.yLab, Form("#splitline{%s keV}{%s (%s)}", L.keV, L.name, chain[L.ch]));
    t->SetTextAlign(21);
    t->SetTextColor(col);
    t->Draw();
  }
  c->SaveAs("output/holders_spectrum.png");
}
