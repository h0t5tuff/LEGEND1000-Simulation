//  RT backgrounds, everything after the simulation: the tube read from the GDML and checked, the runs read and checked,
//  the background as built, and every steel (top) / Cu / EFCu (bottom) design with seams L1, L2, by reweighting
//    root -l -b -q 'ana/background.C("output/tl208*.root=Tl208,output/bi214*.root=Bi214")'
//  energies are what Geant4 deposited, no detector model. cuts: M1, and an argon veto on the energy deposited in the UGLAr

#include "TRegexp.h"
#include <fstream>
#include <sstream>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  // MAJORANA's background window [Haufe et al., arXiv:2209.10592]: 1950-2350 keV minus 10 keV around Qbb, 2103.5 (Tl208
  // single escape), 2118.5 and 2204.1 keV (Bi214) [Arnquist et al., arXiv:2207.07638]. the BI is its counts / 360 keV
  const double winLo = 1950, winHi = 2350, winGap[4] = {2039.0, 2103.5, 2118.5, 2204.1}, winGapHalf = 5.0;
  bool inWindow(double e) { if (e < winLo || e > winHi) return false; for (double g : winGap) if (fabs(e - g) <= winGapHalf) return false; return true; }
  const double winWidth = winHi - winLo - 4 * 2 * winGapHalf; // 360 keV
  const double m1_keV = 5.0;                  // a detector fired above this
  const double lar_keV = 50.0;                // the argon veto rejects a decay leaving more than this in the UGLAr: ~4 photoelectrons (Ralph's slides, p. 4)
  const double geMass_kg = 1000.0;            // the array the BI is per kg of
  const double bgGoal = 1e-5;                 // LEGEND-1000 target [cts/(keV kg yr)] at Qbb
  const double rtBudget = 1e-5;               // the tube's share of it: a design passes when its 90% bound after all cuts is below this
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000;        // spectrum [keV]
  const int    eBins = 200;
  const int    NDB = 10;                      // depth slabs of the hit table
  const int    NZ = 625;                      // ~1 cm depth bins: where each section's decays lie

  const char  *mat[3] = {"steel", "Cu", "EFCu"};      // top to bottom
  const double density[3] = {7900.0, 8960.0, 8930.0}; // kg/m^3
  // chain activities [uBq/kg], {232Th, 238U}, from Ralph's MaterialMix slides (23 Jun 2026): steel from Bernhard (p. 5),
  // OFHC Cu from the MAJORANA assay paper (p. 15), EFCu from M. Green, CD-1 (p. 5)
  const double chain[3][2] = {{1000, 2500}, {1.1, 1.3}, {0.37, 0.19}};
  const double branch[2] = {0.3594, 1.0};     // decays of the simulated nuclide per chain decay: Tl208 is 35.94% of 232Th (Bi212 branch), Bi214 every 238U
  int iso(const std::string &s) { return s.find("Bi") != std::string::npos ? 1 : 0; } // 0 Tl208, 1 Bi214
  double activity(const std::string &s, int m) { return chain[m][iso(s)] * branch[iso(s)] * 1e-6; } // Bq/kg of the simulated nuclide
}

//-------------------------------------------------------------------------------
//  2. The tube, read from the GDML:
struct Outline // the (r, z) corners of one GDML polycone
{
  std::string name;
  std::vector<double> r, z;
  double rmin = 0, rmax = 0, zmin = 0, zmax = 0;
  int size() const { return (int)r.size(); }
  double radiusAt(double zq) const // the outermost radius at height zq, -1 outside
  {
    double best = -1;
    for (int i = 0; i + 1 < size(); i++)
    {
      if ((zq < z[i] && zq < z[i + 1]) || (zq > z[i] && zq > z[i + 1])) continue;
      double rr = fabs(z[i + 1] - z[i]) < 1e-12 ? std::max(r[i], r[i + 1]) : r[i] + (zq - z[i]) / (z[i + 1] - z[i]) * (r[i + 1] - r[i]);
      best = std::max(best, rr);
    }
    return best;
  }
  bool contains(double rq, double zq) const // even-odd ray cast
  {
    bool in = false;
    for (int i = 0, j = size() - 1; i < size(); j = i++)
      if ((z[i] > zq) != (z[j] > zq) && rq < r[i] + (zq - z[i]) / (z[j] - z[i]) * (r[j] - r[i])) in = !in;
    return in;
  }
  // inside or within tol of it (near), at least tol inside (deep). nudged in r AND z: on the lid's shallow cone a float32 z error of 0.2 um is 40 um in r
  bool near(double rq, double zq, double t) const { return contains(rq, zq) || contains(rq + t, zq) || contains(fabs(rq - t), zq) || contains(rq, zq + t) || contains(rq, zq - t); }
  bool deep(double rq, double zq, double t) const { return contains(rq, zq) && contains(rq + t, zq) && contains(fabs(rq - t), zq) && contains(rq, zq + t) && contains(rq, zq - t); }
  double areaAt(double zq) const // cross-section [mm^2]: the crossings pair into rings, so the argon around the lid counts right
  {
    std::vector<double> x;
    for (int i = 0, j = size() - 1; i < size(); j = i++)
      if ((z[i] > zq) != (z[j] > zq)) x.push_back(r[i] + (zq - z[i]) / (z[j] - z[i]) * (r[j] - r[i]));
    std::sort(x.begin(), x.end());
    double a = 0;
    for (size_t k = 0; k + 1 < x.size(); k += 2) a += x[k + 1] * x[k + 1] - x[k] * x[k];
    return M_PI * a;
  }
  double volume(double dz = 0.5) const { double v = 0; for (double zz = zmin; zz < zmax; zz += dz) v += areaAt(zz + 0.5 * std::min(dz, zmax - zz)) * std::min(dz, zmax - zz); return v; }
};

static Outline readOutline(const std::string &file, const char *solid) // one <genericPolycone> out of the GDML text
{
  Outline o;
  o.name = solid;
  std::ifstream in(file.c_str());
  const std::string wanted = std::string("name=\"") + solid + "\"";
  bool on = false;
  for (std::string line; std::getline(in, line);)
  {
    if (!on) { on = line.find("<genericPolycone") != std::string::npos && line.find(wanted) != std::string::npos; continue; }
    if (line.find("</genericPolycone>") != std::string::npos) break;
    double rr, zz;
    if (sscanf(line.c_str(), " <rzpoint r=\"%lf\" z=\"%lf\"", &rr, &zz) == 2) { o.r.push_back(rr); o.z.push_back(zz); }
  }
  if (o.size())
  {
    o.rmin = *std::min_element(o.r.begin(), o.r.end()); o.rmax = *std::max_element(o.r.begin(), o.r.end());
    o.zmin = *std::min_element(o.z.begin(), o.z.end()); o.zmax = *std::max_element(o.z.begin(), o.z.end());
  }
  return o;
}

struct RT
{
  std::string file;
  Outline wall, argon, ofhcOuter, ofhcInner, ssOuter, ssInner; // the tube (EFCu mother), the UGLAr inside it, the two shells
  double zBottom = 0, zTop = 0, seamOFHC = 0, seamSS = 0, zHead = 0; // zHead: the widest point, where the head meets the barrel
  bool ok = false;
  const char *sectionAt(double zq) const { return zq < seamOFHC ? "EFCu" : (zq < seamSS ? "OFHC" : "SS"); }
  const char *materialAt(double rq, double zq, double t = 1e-3) const // 1 um of slack: a float32 vertex lands off the surface it was drawn on
  {
    if (!wall.near(rq, zq, t)) return "outside";
    if (ofhcOuter.near(rq, zq, t) && !ofhcInner.deep(rq, zq, t)) return "OFHC";
    if (ssOuter.near(rq, zq, t) && !ssInner.deep(rq, zq, t)) return "SS";
    return argon.deep(rq, zq, t) ? "LAr" : "EFCu";
  }
  int physAt(double rq, double zq) const { const char *m = materialAt(rq, zq); return !strcmp(m, "OFHC") ? 1 : (!strcmp(m, "SS") ? 2 : (!strcmp(m, "EFCu") ? 0 : -1)); } // 0 mother, 1 OFHC, 2 SS
  double wallVolume(double z0, double z1, double dz = 0.5) const // [mm^3]: the tube's solid minus the argon in it
  {
    z0 = std::max(std::min(z0, z1), zBottom); z1 = std::min(std::max(z0, z1), zTop);
    double v = 0;
    for (double z = z0; z < z1; z += dz) { double s = std::min(dz, z1 - z); v += (wall.areaAt(z + 0.5 * s) - argon.areaAt(z + 0.5 * s)) * s; }
    return v;
  }
  double physVolume(int p) const // [mm^3]: the shells are outer minus inner bound, the mother is the rest of the wall
  {
    double ofhc = ofhcOuter.volume() - ofhcInner.volume(), ss = ssOuter.volume() - ssInner.volume();
    return p == 1 ? ofhc : (p == 2 ? ss : wallVolume(zBottom, zTop) - ofhc - ss);
  }
};

static std::string findFile(const std::string &name) // in the run directory or up to 5 levels above it
{
  std::string p = name;
  for (int up = 0; up < 6; up++, p = "../" + p) if (!gSystem->AccessPathName(p.c_str())) return p;
  return "";
}

static RT loadRT(const std::string &gdml)
{
  RT rt;
  rt.file = findFile(gdml);
  if (rt.file.empty()) { printf("ERROR: %s not found\n", gdml.c_str()); return rt; }
  rt.wall = readOutline(rt.file, "reentrancetube");
  rt.argon = readOutline(rt.file, "undergroundlar");
  rt.ofhcOuter = readOutline(rt.file, "ofhc_cu_outer_bound");
  rt.ofhcInner = readOutline(rt.file, "ofhc_cu_inner_bound");
  rt.ssOuter = readOutline(rt.file, "ss_316l_outer_bound");
  rt.ssInner = readOutline(rt.file, "ss_316l_inner_bound");
  if (!rt.wall.size() || !rt.ofhcOuter.size() || !rt.ssOuter.size()) { printf("ERROR: %s has no re-entrant tube\n", rt.file.c_str()); return rt; }
  rt.zBottom = rt.wall.zmin; rt.zTop = rt.wall.zmax;
  rt.seamOFHC = rt.ofhcOuter.zmin; rt.seamSS = rt.ssOuter.zmin;
  rt.zHead = rt.wall.z[std::max_element(rt.wall.r.begin(), rt.wall.r.end()) - rt.wall.r.begin()];
  rt.ok = true;
  return rt;
}

template <class F> static void scanTree(TTree *t, const std::vector<const char *> &cols, F f) // the named columns as doubles, row by row
{
  std::vector<TTreeFormula *> col;
  for (auto c : cols) col.push_back(new TTreeFormula(c, c, t));
  std::vector<double> v(cols.size());
  for (Long64_t i = 0; i < t->GetEntries(); i++) { t->GetEntry(i); for (size_t k = 0; k < col.size(); k++) v[k] = col[k]->EvalInstance(); f(v.data()); }
  for (auto c : col) delete c;
}

static std::string outPath(const std::string &name) { gSystem->mkdir("output", kTRUE); return "output/" + name; }

//-------------------------------------------------------------------------------
//  3. Check the tube: one clean surface, the endcap flush with the barrel, and how it differs from the unmodified geometry
static void checkTube(const RT &s)
{
  const Outline &w = s.wall;
  const std::vector<double> &r = w.r, &z = w.z;
  const int n = w.size();
  auto turn = [](double ar, double az, double br, double bz, double cr, double cz) { double v = (bz - az) * (cr - br) - (br - ar) * (cz - bz); return fabs(v) < 1e-9 ? 0 : (v > 0 ? 1 : 2); };
  int nAxis = 0, nDup = 0, nSpike = 0, nCross = 0;
  for (int i = 0; i < n; i++)
  {
    int i2 = (i + 1) % n;
    if (r[i] == 0) nAxis++;                                                   // only the bottom tip and the top sit on the axis
    for (int j = 0; j < i; j++) if (r[i] == r[j] && z[i] == z[j]) nDup++;     // the same point twice
    if (i > 0 && i < n - 1)
    {
      double ar = r[i] - r[i - 1], az = z[i] - z[i - 1], br = r[i + 1] - r[i], bz = z[i + 1] - z[i], m = sqrt((ar * ar + az * az) * (br * br + bz * bz));
      if (m > 0 && (ar * br + az * bz) / m < -0.99) nSpike++;                 // the outline folds straight back
    }
    for (int j = i + 2; j < n; j++)                                           // segments that are not neighbours must not cross
    {
      int j2 = (j + 1) % n;
      if (i == 0 && j == n - 1) continue;
      if (turn(r[i], z[i], r[i2], z[i2], r[j], z[j]) != turn(r[i], z[i], r[i2], z[i2], r[j2], z[j2]) &&
          turn(r[j], z[j], r[j2], z[j2], r[i], z[i]) != turn(r[j], z[j], r[j2], z[j2], r[i2], z[i2])) nCross++;
    }
  }
  const double tol = 0.05, zBarrel = 0.5 * (s.seamOFHC + s.seamSS), rBarrel = w.radiusAt(zBarrel), step = w.rmax - rBarrel;
  printf("[3] the tube in %s: z %.1f .. %.1f mm (%.3f m), %d points, seams EFCu | %.1f | OFHC | %.1f | SS\n", s.file.c_str(), s.zBottom, s.zTop,
         (s.zTop - s.zBottom) / 1000, n, s.seamOFHC, s.seamSS);
  printf("    clean surface: on-axis points %d (expect 2), duplicates %d, spikes %d, self-intersections %d (expect 0)\n", nAxis, nDup, nSpike, nCross);
  printf("    endcap: radius %.4f at z %.1f, barrel %.4f: step %+.4f mm%s\n", w.rmax, s.zHead, rBarrel, step, fabs(step) <= tol ? "" : "  <- NOT flush");
  printf("    shells against the barrel: OFHC %+.4f, SS %+.4f mm\n", s.ofhcOuter.rmax - rBarrel, s.ssOuter.rmax - rBarrel);
  printf("    wall thickness [mm]:");
  for (double zq : {s.zHead, 0.5 * (s.zHead + s.seamOFHC), s.seamOFHC - 10, zBarrel, 0.5 * (s.seamSS + s.zTop)})
    printf("  %s %.2f (z %.0f)", s.sectionAt(zq), w.radiusAt(zq) - s.argon.radiusAt(zq), zq);
  printf("\n");
  RT m0;
  std::string mint = findFile("l1000.gdml");
  if (!mint.empty()) m0 = loadRT(mint);
  if (m0.ok)
    printf("    against the unmodified %s: bottom %+.1f, top %+.1f, length %+.1f, seams %+.1f / %+.1f, barrel radius %+.1f, endcap radius %+.1f mm\n",
           mint.c_str(), s.zBottom - m0.zBottom, s.zTop - m0.zTop, (s.zTop - s.zBottom) - (m0.zTop - m0.zBottom), s.seamOFHC - m0.seamOFHC,
           s.seamSS - m0.seamSS, rBarrel - m0.wall.radiusAt(0.5 * (m0.seamOFHC + m0.seamSS)), w.rmax - m0.wall.rmax);
  bool ok = nAxis == 2 && !nDup && !nSpike && !nCross && fabs(step) <= tol && fabs(s.ofhcOuter.rmax - rBarrel) <= tol && fabs(s.ssOuter.rmax - rBarrel) <= tol;
  printf("    tube: %s\n\n", ok ? "PASS" : "FAIL");

  auto c = new TCanvas("c_rt", "", 1250, 780); // whole tube | endcap | endcap/barrel junction
  c->Divide(3, 1);
  struct { const char *title; double z0, z1, r0, r1; } pad[] = {{"whole tube", s.zBottom - 150, s.zTop + 150, 0, w.rmax * 1.12},
      {"endcap bottom profile", s.zBottom - 30, s.zHead + 250, 0, w.rmax * 1.06}, {"endcap / barrel junction", s.zHead - 120, s.zHead + 160, rBarrel - 2.5, w.rmax + 1.0}};
  auto line = [](const Outline &o, int col, int style) { auto g = new TGraph(o.size(), o.r.data(), o.z.data()); g->SetLineColor(col); g->SetLineStyle(style); g->SetMarkerColor(col); g->SetMarkerStyle(20); g->SetMarkerSize(0.3); return g; };
  for (int i = 0; i < 3; i++)
  {
    c->cd(i + 1)->SetGrid();
    gPad->SetLeftMargin(0.17);
    auto fr = gPad->DrawFrame(pad[i].r0, pad[i].z0, pad[i].r1, pad[i].z1);
    fr->SetTitle(Form("%s;r [mm];z [mm]", pad[i].title));
    fr->GetYaxis()->SetTitleOffset(1.7);
    fr->GetXaxis()->SetNdivisions(506);
    if (m0.ok) line(m0.wall, kGray + 1, 2)->Draw("L SAME");  // unmodified, dashed grey
    line(w, kAzure + 2, 1)->Draw(i ? "LP SAME" : "L SAME");   // zoomed pads show the corner points too
    if (i == 0) for (double v : {s.seamOFHC, s.seamSS}) { auto l = new TLine(pad[0].r0, v, pad[0].r1, v); l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw(); }
    if (i == 2)
    {
      for (double v : {rBarrel, w.rmax}) { auto l = new TLine(v, pad[2].z0, v, pad[2].z1); l->SetLineStyle(2); l->SetLineColor(v == w.rmax ? kRed + 1 : kGreen + 2); l->Draw(); }
      auto t = new TLatex(rBarrel - 2.2, s.zHead + 45, Form("step = %+.4f mm", step));
      t->SetTextColor(kRed + 1); t->SetTextSize(0.04); t->Draw();
    }
  }
  c->SaveAs(outPath("tube.png").c_str());
}

//-------------------------------------------------------------------------------
//  4. One run, reduced to what a reweighting needs:
struct Hit { Long64_t ev; float e; bool m1, lar; }; // one detector's energy in one decay; ev: the decay's slot; lar: the decay passes the argon veto

struct Run
{
  std::string file, iso;
  Long64_t nsim = 0;
  long nFiles = 1, nMat[5] = {0, 0, 0, 0, 0}; // decays in EFCu, OFHC, SS, the argon, outside the wall
  long nPv[3] = {0, 0, 0};              // decays per section: 0 mother (EFCu), 1 OFHC, 2 SS
  long nZPv[3][cfg::NZ] = {{0}};        // the same, per ~1 cm of depth
  double mcDensity[3] = {0, 0, 0};      // simulated decays per m^3 in each section, measured
  long dN[cfg::NDB] = {0}, dH[cfg::NDB] = {0}, dW[cfg::NDB] = {0}; // per depth slab: decays, hits, window hits
  std::vector<float> depth;             // per decay WITH a Ge hit: mm below the top of the tube
  std::vector<signed char> pv;          // its section
  std::vector<short> win[3];            // its hits in the window: no cut, M1, M1 + argon veto
  std::vector<Hit> hits;
  long nHit[3] = {0, 0, 0};             // hits: no cut, M1, M1 + argon veto
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
  std::vector<float> depth(n, -1), eLAr(n, 0);
  std::vector<short> nFired(n, 0), win[3] = {std::vector<short>(n, 0), std::vector<short>(n, 0), std::vector<short>(n, 0)};
  std::vector<signed char> pv(n, -1);
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it()) // each Ge detector's steps summed per decay
  {
    TString nm = k->GetName();
    if (nm.Length() != 5 || nm[0] != 'V' || !TString(nm(1, 4)).IsDigit()) continue; // a Ge detector: V + 4 digits
    auto t = (TTree *)d->Get(nm);
    std::map<int, float> sum;
    scanTree(t, {"evtid", "edep_in_keV"}, [&](const double *v) { sum[(int)v[0]] += v[1]; });
    for (auto &p : sum)
      if (p.first >= 0 && p.first < n && p.second > cfg::m1_keV) { nFired[p.first]++; r.hits.push_back({p.first, p.second, false, false}); }
  }
  if (auto tl = (TTree *)d->Get("undergroundlar"))
    scanTree(tl, {"evtid", "edep_in_keV"}, [&](const double *v) { if (v[0] >= 0 && v[0] < n) eLAr[(int)v[0]] += v[1]; });
  else printf("WARNING: no stp/undergroundlar in %s - the argon veto keeps everything\n", file.c_str());
  const char *matName[5] = {"EFCu", "OFHC", "SS", "LAr", "outside"};
  scanTree(vtx, {"evtid", "xloc_in_m", "yloc_in_m", "zloc_in_m"}, [&](const double *v) {
    int ev = (int)v[0]; // rows are not in event order with -t: index by evtid
    if (ev < 0 || ev >= n) return;
    double x = v[1] * 1000, y = v[2] * 1000, z = v[3] * 1000;
    const char *m = rt.materialAt(sqrt(x * x + y * y), z);
    int k = 0;
    while (k < 4 && strcmp(m, matName[k])) k++;
    r.nMat[k]++;
    depth[ev] = rt.zTop - z;
    pv[ev] = k < 3 ? k : -1; // EFCu is the mother (0), OFHC 1, SS 2
  });
  const double tot = rt.zTop - rt.zBottom;
  for (auto &h : r.hits)
  {
    h.m1 = nFired[h.ev] == 1;
    h.lar = eLAr[h.ev] <= cfg::lar_keV;
    bool pass[3] = {true, h.m1, h.m1 && h.lar};
    for (int c = 0; c < 3; c++) { r.nHit[c] += pass[c]; if (pass[c] && cfg::inWindow(h.e)) win[c][h.ev]++; }
  }
  for (Long64_t i = 0; i < n; i++) // every decay, as counts
  {
    r.dN[bin(depth[i], tot, cfg::NDB)]++;
    r.dW[bin(depth[i], tot, cfg::NDB)] += win[0][i];
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
      for (int c = 0; c < 3; c++) r.win[c].push_back(win[c][h.ev]);
    }
    h.ev = slot[h.ev];
  }
  delete f;
  return r;
}

//-------------------------------------------------------------------------------
//  5. Combine the runs: files of one isotope are one run split into jobs
static std::vector<std::string> expand(const std::string &pat) // "output/tl208*.root" -> every match, sorted
{
  if (pat.find_first_of("*?") == std::string::npos) return {pat};
  size_t sl = pat.find_last_of('/');
  std::string dir = sl == std::string::npos ? "." : pat.substr(0, sl), base = pat.substr(sl + 1);
  TRegexp re(base.c_str(), kTRUE); // a shell wildcard
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

void background(const char *runs = "output/tl208*.root=Tl208,output/bi214*.root=Bi214", const char *gdml = "KSendcap_l1kGeometry.gdml")
{
  gErrorIgnoreLevel = kWarning; // no "png file has been created" notices in the printout
  RT rt = loadRT(gdml);
  if (!rt.ok) return;
  checkTube(rt);
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
      for (int c = 0; c < 3; c++) { into->win[c].insert(into->win[c].end(), r.win[c].begin(), r.win[c].end()); into->nHit[c] += r.nHit[c]; }
      for (int k = 0; k < 5; k++) into->nMat[k] += r.nMat[k];
      for (int p = 0; p < 3; p++) { into->nPv[p] += r.nPv[p]; for (int b = 0; b < cfg::NZ; b++) into->nZPv[p][b] += r.nZPv[p][b]; }
      for (int b = 0; b < cfg::NDB; b++) { into->dN[b] += r.dN[b]; into->dH[b] += r.dH[b]; into->dW[b] += r.dW[b]; }
      into->nsim += r.nsim; into->nFiles++; into->file += "+" + r.file;
    }
  }
  if (R.empty()) { printf("no usable runs in '%s'\n", runs); return; }

  const double depthTot = rt.zTop - rt.zBottom, L1ks = rt.zTop - rt.seamSS, L2ks = rt.zTop - rt.seamOFHC; // KS seams as depths
  const double UL90 = 2.30, binW = (cfg::eHi - cfg::eLo) / cfg::eBins;
  const char *cutName[3] = {"no cut", "M1", "M1 + argon"};
  double physV[3]; // m^3: mother, OFHC shell, SS shell
  for (int p = 0; p < 3; p++) physV[p] = rt.physVolume(p) * 1e-9;
  for (auto &r : R) for (int p = 0; p < 3; p++) r.mcDensity[p] = r.nPv[p] / physV[p]; // remage fills the mother sparser than its daughters: measure, never assume
  // one simulated decay in section p, where material m sits = rho x A x 1 yr real decays per m^3, over the MC density, per kg of Ge
  auto weight = [&](const Run &r, int m, int p) { return p < 0 || r.mcDensity[p] <= 0 ? 0.0 : cfg::density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[p] * cfg::geMass_kg); };
  auto massOf = [&](double d0, double d1, int m) { return d1 > d0 ? cfg::density[m] * rt.wallVolume(rt.zTop - d1, rt.zTop - d0) * 1e-9 : 0.0; };
  auto ulN = [](long n) { const double t[11] = {2.30, 3.89, 5.32, 6.68, 7.99, 9.27, 10.53, 11.77, 12.99, 14.21, 15.41}; return n <= 10 ? t[n] : n + 1.28 * sqrt((double)n) + 1; }; // Poisson 90% upper limit
  printf("window   : %.0f-%.0f keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's BEW): %.0f keV\n", cfg::winLo, cfg::winHi, cfg::winWidth);
  printf("cuts     : M1 (exactly one detector above %.0f keV), argon veto (UGLAr deposit <= %.0f keV). no detector model: deposited energy\n",
         cfg::m1_keV, cfg::lar_keV);
  printf("activity : chain [uBq/kg] 232Th / 238U: steel %g / %g, Cu %g / %g, EFCu %g / %g; Tl208 is %.2f%% of 232Th, Bi214 %.0f%% of 238U\n\n",
         cfg::chain[0][0], cfg::chain[0][1], cfg::chain[1][0], cfg::chain[1][1], cfg::chain[2][0], cfg::chain[2][1], 100 * cfg::branch[0], 100 * cfg::branch[1]);

  for (auto &r : R)
  {
//-------------------------------------------------------------------------------
//  6. Per run: the decays, the hits, the normalisation
    long w[3] = {0, 0, 0};
    for (int c = 0; c < 3; c++) for (auto k : r.win[c]) w[c] += k;
    printf("=== %s   (%s, %lld decays)\n", r.iso.c_str(), r.file.c_str(), r.nsim);
    printf("[6] decays: EFCu %.1f%%, OFHC %.1f%%, SS %.1f%%; in the argon %ld, outside the wall %ld%s%s\n", 100.0 * r.nMat[0] / r.nsim, 100.0 * r.nMat[1] / r.nsim,
           100.0 * r.nMat[2] / r.nsim, r.nMat[3], r.nMat[4], r.nMat[3] || r.nMat[4] ? "  <- CONFINEMENT BUG" : "", r.nFiles > 1 ? Form("   (merged from %ld files)", r.nFiles) : "");
    printf("    hits: %ld, M1 %ld, M1 + argon %ld; in the window: %ld, %ld, %ld\n", r.nHit[0], r.nHit[1], r.nHit[2], w[0], w[1], w[2]);
    const int asBuilt[3] = {2, 1, 0}; // mother EFCu, OFHC shell Cu, SS shell steel
    const char *pvName[3] = {"mother", "OFHC", "SS"};
    printf("    %-7s %-6s %9s %8s %10s %9s %11s %14s\n", "section", "mat", "V [m^3]", "M [kg]", "A [Bq]", "MC dec.", "MC per m^3", "decays/yr per");
    for (int p = 0; p < 3; p++)
    {
      int m = asBuilt[p];
      double mass = physV[p] * cfg::density[m], A = mass * cfg::activity(r.iso, m);
      printf("    %-7s %-6s %9.5f %8.1f %10.3e %9ld %11.0f %14.3g\n", pvName[p], cfg::mat[m], physV[p], mass, A, r.nPv[p], r.mcDensity[p], r.nPv[p] ? A * cfg::secPerYear / r.nPv[p] : 0.0);
    }
    double dAvg = 0.5 * (r.mcDensity[1] + r.mcDensity[2]);
    printf("    sampling density mother / shells = %.3f   (1.000 would be uniform; the weights use the measured density)\n", dAvg > 0 ? r.mcDensity[0] / dAvg : 0.0);

//-------------------------------------------------------------------------------
//  7. Per run: where the hits come from, and the decays the steel needs
    printf("[7] depth [m]          decays     hits    hits/decay   window\n");
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
//  8. Background index as built, per chain and section:
  printf("[8] background index as built [cts/(keV kg yr)], +- MC statistics; a section with no window hit gets its 90%% limit\n");
  printf("    %-6s %-6s %12s   %-22s %-22s %-22s\n", "chain", "sect", "win hits", cutName[0], cutName[1], cutName[2]);
  double gB[3] = {0, 0, 0}, gV[3] = {0, 0, 0};
  long gN[3] = {0, 0, 0}; // window hits behind each total
  auto total = [](double b, double v, long n, bool goal) -> std::string {
    if (!n) return "no window hit";
    return goal ? Form("%.2e +- %.1e (%.2g x goal)", b, sqrt(v), b / cfg::bgGoal) : Form("%.2e +- %.1e", b, sqrt(v));
  };
  const double edgeKS[4] = {0, L1ks, L2ks, depthTot};
  for (auto &r : R)
  {
    double tB[3] = {0, 0, 0}, tV[3] = {0, 0, 0};
    long tN[3] = {0, 0, 0};
    for (int s = 0; s < 3; s++) // the slabs as built: steel, Cu, EFCu
    {
      double B[3] = {0, 0, 0}, V[3] = {0, 0, 0}, wbar = 0;
      long n[3] = {0, 0, 0}, nsec = 0;
      for (int p = 0; p < 3; p++) // mean weight of a decay in this section, for the limit
        for (int b = bin(edgeKS[s], depthTot, cfg::NZ); b < std::min(cfg::NZ, (int)ceil(edgeKS[s + 1] / depthTot * cfg::NZ)); b++)
        { nsec += r.nZPv[p][b]; wbar += r.nZPv[p][b] * weight(r, s, p) / cfg::winWidth; }
      wbar = nsec ? wbar / nsec : 0;
      for (size_t i = 0; i < r.depth.size(); i++)
      {
        if (slabOf(r.depth[i], L1ks, L2ks) != s) continue;
        double wt = weight(r, s, r.pv[i]) / cfg::winWidth;
        for (int c = 0; c < 3; c++) { B[c] += r.win[c][i] * wt; V[c] += r.win[c][i] * wt * wt; n[c] += r.win[c][i]; }
      }
      std::string cell[3];
      for (int c = 0; c < 3; c++)
      {
        cell[c] = n[c] ? Form("%.2e +- %.1e", B[c], sqrt(V[c])) : Form("< %.2e (90%%)", UL90 * wbar);
        if (n[c]) { tB[c] += B[c]; tV[c] += V[c]; tN[c] += n[c]; }
      }
      printf("    %-6s %-6s %12s   %-22s %-22s %-22s\n", s ? "" : r.iso.c_str(), cfg::mat[s], Form("%ld/%ld/%ld", n[0], n[1], n[2]), cell[0].c_str(), cell[1].c_str(), cell[2].c_str());
    }
    printf("    %-6s %-6s %12s   %-22s %-22s %-22s\n", "", "tube", "", total(tB[0], tV[0], tN[0], false).c_str(), total(tB[1], tV[1], tN[1], false).c_str(),
           total(tB[2], tV[2], tN[2], false).c_str());
    for (int c = 0; c < 3; c++) { gB[c] += tB[c]; gV[c] += tV[c]; gN[c] += tN[c]; }
  }
  printf("    ALL CHAINS   %s: %s   %s: %s   %s: %s\n", cutName[0], total(gB[0], gV[0], gN[0], true).c_str(), cutName[1], total(gB[1], gV[1], gN[1], true).c_str(),
         cutName[2], total(gB[2], gV[2], gN[2], true).c_str());
  printf("    win hits: no cut / M1 / M1 + argon. sections with no window hit are left out of the totals\n");

//-------------------------------------------------------------------------------
//  9. Three-material designs: the MC alone, after all cuts, judged by its 90% upper bound
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
          if (r.win[2][i] && slabOf(r.depth[i], L1, L2) == s) { B += r.win[2][i] * weight(r, s, r.pv[i]) / cfg::winWidth; n += r.win[2][i]; }
      }
      bi += B;
      up += n ? ulN(n) * B / n : UL90 * wmax;
    }
    return up;
  };
  auto passes = [&](double L1, double L2) { double bi; return judge(L1, L2, bi) <= cfg::rtBudget; };
  auto deepestSteel = [&](double L2) { double l1 = -1; for (double L1 = 10; L1 <= L2 - 10 && passes(L1, L2); L1 += 10) l1 = L1; return l1; }; // 1 cm grid, always some Cu
  printf("\n[9] designs: steel to L1, Cu to L2, EFCu below. MC alone, after M1 + argon; passes if its 90%% bound <= the budget %.0e\n", cfg::rtBudget);
  printf("    %-7s %9s %9s %10s %9s %17s %11s %11s %9s\n", "L2 [m]", "EFCu [kg]", "steel to", "steel [kg]", "Cu [kg]", "SS:Cu:EFCu mass", "BI", "90% bound", "x budget");
  std::vector<double> rows = {L2ks};
  for (double L2 = 3000; L2 < depthTot; L2 += 250) rows.push_back(L2);
  std::sort(rows.begin(), rows.end());
  double L2best = -1, L1best = -1;
  for (double L2 : rows)
  {
    double l1 = L2 == L2ks ? L1ks : deepestSteel(L2), bi = 0;
    if (l1 < 0) { printf("    %-7.2f %9.0f %9s   no steel can be shown to pass yet\n", L2 / 1000, massOf(L2, depthTot, 2), "-"); continue; }
    double up = judge(l1, L2, bi), ms = massOf(0, l1, 0), mc = massOf(l1, L2, 1), me = massOf(L2, depthTot, 2), mt = ms + mc + me;
    printf("    %-7.2f %9.0f %7.2f m %10.0f %9.0f %17s %11.2e %11.2e %9.2f%s\n", L2 / 1000, me, l1 / 1000, ms, mc, Form("%.0f : %.0f : %.0f %%", 100 * ms / mt, 100 * mc / mt, 100 * me / mt),
           bi, up, up / cfg::rtBudget, L2 == L2ks ? "   <- KS as built" : "");
    if (L2 != L2ks && up <= cfg::rtBudget && L2 > L2best) { L2best = L2; L1best = l1; }
  }
  if (L2best > 0) printf("    least EFCu on this grid, then most steel: L2 %.2f m, L1 %.2f m\n", L2best / 1000, L1best / 1000);
  else printf("    no design with steel passes on the MC alone: [7] says how many decays its slabs need\n");
  printf("    masses in kg; steel includes the lid, which sits in the top slab\n");

//-------------------------------------------------------------------------------
//  10. Draw: the spectrum as built, and the 90% bound of every design
  TH1D *sp[3];
  for (int c = 0; c < 3; c++) sp[c] = new TH1D(Form("sp%d", c), "RT background in the germanium, as built;energy [keV];cts / (keV kg yr)", cfg::eBins, cfg::eLo, cfg::eHi);
  for (auto &r : R)
    for (auto &h : r.hits)
    {
      double wt = weight(r, slabOf(r.depth[h.ev], L1ks, L2ks), r.pv[h.ev]) / binW;
      sp[0]->Fill(h.e, wt);
      if (h.m1) sp[1]->Fill(h.e, wt);
      if (h.m1 && h.lar) sp[2]->Fill(h.e, wt);
    }
  const int NS = 40;
  auto hMap = new TH2D("hMap", "MC 90% bound after M1 + argon, x budget;L1 steel/Cu seam [m];L2 Cu/EFCu seam [m]", NS, 0, depthTot / 1000, NS, 0, depthTot / 1000);
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
  auto leg = new TLegend(0.6, 0.74, 0.88, 0.88);
  const int col[3] = {kAzure + 2, kOrange + 7, kRed + 1};
  for (int k = 0; k < 3; k++)
  {
    sp[k]->SetLineColor(col[k]);
    sp[k]->SetLineWidth(2);
    sp[k]->Draw(k ? "HIST SAME" : "HIST");
    leg->AddEntry(sp[k], cutName[k], "l");
  }
  sp[0]->SetMinimum(0.1 * std::max(1e-12, sp[2]->GetMinimum(0)));
  leg->Draw();
  c->cd(2)->SetLogz();
  gPad->SetRightMargin(0.15);
  hMap->SetMinimum(std::min(1.0, 0.5 * hMap->GetMinimum(0))); // from the budget itself, so passing designs would stand out
  hMap->GetZaxis()->SetMoreLogLabels();
  hMap->Draw("COLZ");
  auto star = new TMarker(L1ks / 1000, L2ks / 1000, 29);
  star->SetMarkerColor(kRed + 1);
  star->SetMarkerSize(2.2);
  star->Draw();
  std::string tag;
  for (auto &r : R) tag += (tag.empty() ? "" : "_") + r.iso;
  std::string png = outPath(tag + "_background.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s and output/tube.png\n", png.c_str());
}
