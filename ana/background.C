#include "TRegexp.h"
#include <fstream>
#include <sstream>
#include <unordered_map>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  // MAJORANA's background window [Haufe et al., arXiv:2209.10592]: 1950-2350 keV minus 10 keV around Qbb, 2103.5 (Tl208
  // single escape), 2118.5 and 2204.1 keV (Bi214) [Arnquist et al., arXiv:2207.07638]. the BI is its counts / 360 keV
  const double winLo = 1950, winHi = 2350, winGap[4] = {2039.0, 2103.5, 2118.5, 2204.1}, winGapHalf = 5.0;
  bool inWindow(double e)
  {
    if (e < winLo || e > winHi)
      return false;
    for (double g : winGap)
      if (fabs(e - g) <= winGapHalf)
        return false;
    return true;
  }
  const double winWidth = winHi - winLo - 4 * 2 * winGapHalf; // 360 keV
  // the ROI shaded on the spectrum: Qbb +- 2 sigma, sigma from the response's resolution (simprod eresmod: FWHM = sqrt(0.5 + 0.001 E) keV)
  const double Qbb = 2039.0, roiHalf = 2 * sqrt(0.5 + 0.001 * Qbb) / 2.35482;
  // the cuts of the LEGEND-1000 simulation production (legend1000-metadata simprod/config, l1000dsg01: tier/evt, pars/geds/psdcuts)
  const double m1_keV = 25.0;  // a detector fired above this; M1: exactly one did
  const double lar_pe = 4.0;   // the LAr veto rejects a decay whose SiPMs see this many photoelectrons, summed
  const double psd_low = -1.5; // the production's A/E cut: keep a hit whose A/E classifier is above this (low side)
  // ...but the analysis sets it as LEGEND does in data ([5]): at the classifier value that keeps this fraction of the Tl208
  // double-escape peak (DEP, single site like 0vbb) in M1 events, net of the continuum under it
  const double depKeep = 0.90, depE = 1592.5, depHalf = 2.5, sideLo = 5.0, sideHi = 10.0; // keV: peak +-2.5, sidebands 5-10 off it
  const int NC = 4;                                                                       // cut levels: none, M1, M1 + LAr, M1 + LAr + A/E
  const char *cutName[NC] = {"no cut", "M1", "M1 + LAr", "M1 + LAr + A/E"};
  const double geMass_kg = 1000.0; // the array the BI is per kg of
  const double bgGoal = 1e-5;      // LEGEND-1000 goal [cts/(keV kg yr)] at Qbb, for scale
  const double simTol = 0.10;      // a design passes when its BI after all cuts is at most this much above the tube's as built, at 90% CL
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000; // spectrum [keV]
  const int eBins = 200;
  const int NDB = 10; // depth slabs of the hit table
  const int NZ = 625; // ~1 cm depth bins: where each section's decays lie

  const char *mat[3] = {"steel", "Cu", "EFCu"}; // top to bottom; their densities are the GDML's (RT::density)
  // chain activities [uBq/kg], {232Th, 238U}, from Ralph's MaterialMix slides (23 Jun 2026): steel from Bernhard (p. 5),
  // OFHC Cu from the MAJORANA assay paper (p. 15), EFCu from M. Green, CD-1 (p. 5)
  const double chain[3][2] = {{1000, 2500}, {1.1, 1.3}, {0.37, 0.19}};
  const double branch[2] = {0.3594, 1.0};                                                           // decays of the simulated nuclide per chain decay: Tl208 is 35.94% (~36%) of the 232Th chain, the alpha branch of Bi212 (ENSDF); Bi214 100% of the 238U chain
  int iso(const std::string &s) { return s.find("Bi") != std::string::npos ? 1 : 0; }               // 0 Tl208, 1 Bi214
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
      if ((zq < z[i] && zq < z[i + 1]) || (zq > z[i] && zq > z[i + 1]))
        continue;
      double rr = fabs(z[i + 1] - z[i]) < 1e-12 ? std::max(r[i], r[i + 1]) : r[i] + (zq - z[i]) / (z[i + 1] - z[i]) * (r[i + 1] - r[i]);
      best = std::max(best, rr);
    }
    return best;
  }
  bool contains(double rq, double zq) const // even-odd ray cast
  {
    bool in = false;
    for (int i = 0, j = size() - 1; i < size(); j = i++)
      if ((z[i] > zq) != (z[j] > zq) && rq < r[i] + (zq - z[i]) / (z[j] - z[i]) * (r[j] - r[i]))
        in = !in;
    return in;
  }
  // inside or within tol of it (near), at least tol inside (deep). nudged in r AND z: on the lid's shallow cone a float32 z error of 0.2 um is 40 um in r
  bool near(double rq, double zq, double t) const { return contains(rq, zq) || contains(rq + t, zq) || contains(fabs(rq - t), zq) || contains(rq, zq + t) || contains(rq, zq - t); }
  bool deep(double rq, double zq, double t) const { return contains(rq, zq) && contains(rq + t, zq) && contains(fabs(rq - t), zq) && contains(rq, zq + t) && contains(rq, zq - t); }
  double areaAt(double zq) const // cross-section [mm^2]: the crossings pair into rings, so the argon around the lid counts right
  {
    std::vector<double> x;
    for (int i = 0, j = size() - 1; i < size(); j = i++)
      if ((z[i] > zq) != (z[j] > zq))
        x.push_back(r[i] + (zq - z[i]) / (z[j] - z[i]) * (r[j] - r[i]));
    std::sort(x.begin(), x.end());
    double a = 0;
    for (size_t k = 0; k + 1 < x.size(); k += 2)
      a += x[k + 1] * x[k + 1] - x[k] * x[k];
    return M_PI * a;
  }
  double volume(double dz = 0.5) const
  {
    double v = 0;
    for (double zz = zmin; zz < zmax; zz += dz)
      v += areaAt(zz + 0.5 * std::min(dz, zmax - zz)) * std::min(dz, zmax - zz);
    return v;
  }
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
    if (!on)
    {
      on = line.find("<genericPolycone") != std::string::npos && line.find(wanted) != std::string::npos;
      continue;
    }
    if (line.find("</genericPolycone>") != std::string::npos)
      break;
    double rr, zz;
    if (sscanf(line.c_str(), " <rzpoint r=\"%lf\" z=\"%lf\"", &rr, &zz) == 2)
    {
      o.r.push_back(rr);
      o.z.push_back(zz);
    }
  }
  if (o.size())
  {
    o.rmin = *std::min_element(o.r.begin(), o.r.end());
    o.rmax = *std::max_element(o.r.begin(), o.r.end());
    o.zmin = *std::min_element(o.z.begin(), o.z.end());
    o.zmax = *std::max_element(o.z.begin(), o.z.end());
  }
  return o;
}

static std::map<std::string, double> readDensities(const std::string &file) // [kg/m^3] of each logical volume's material, out of the GDML text
{
  std::map<std::string, double> material, volume; // materials come before the volumes that use them
  std::ifstream in(file.c_str());
  std::string mat, vol;
  auto attr = [](const std::string &line, const char *key)
  { size_t a = line.find(key); if (a == std::string::npos) return std::string(); a += strlen(key); return line.substr(a, line.find('"', a) - a); };
  for (std::string line; std::getline(in, line);)
  {
    if (line.find("<material ") != std::string::npos)
      mat = attr(line, "name=\"");
    else if (!mat.empty() && line.find("<D ") != std::string::npos)
    {
      material[mat] = 1000 * atof(attr(line, "value=\"").c_str());
      mat.clear();
    } // no unit: g/cm^3
    else if (line.find("<volume ") != std::string::npos)
      vol = attr(line, "name=\"");
    else if (!vol.empty() && line.find("<materialref ") != std::string::npos)
    {
      volume[vol] = material[attr(line, "ref=\"")];
      vol.clear();
    }
  }
  return volume;
}

struct RT
{
  std::string file;
  Outline wall, argon, ofhcOuter, ofhcInner, ssOuter, ssInner;       // the tube (EFCu mother), the UGLAr inside it, the two shells
  double density[3] = {0, 0, 0};                                     // kg/m^3, the GDML's: steel, Cu, EFCu (the SS shell, the OFHC shell, the mother)
  double zBottom = 0, zTop = 0, seamOFHC = 0, seamSS = 0, zHead = 0; // zHead: the widest point, where the head meets the barrel
  bool ok = false;
  const char *sectionAt(double zq) const { return zq < seamOFHC ? "EFCu" : (zq < seamSS ? "OFHC" : "SS"); }
  const char *materialAt(double rq, double zq, double t = 1e-3) const // 1 um of slack: a float32 vertex lands off the surface it was drawn on
  {
    if (!wall.near(rq, zq, t))
      return "outside";
    if (ofhcOuter.near(rq, zq, t) && !ofhcInner.deep(rq, zq, t))
      return "OFHC";
    if (ssOuter.near(rq, zq, t) && !ssInner.deep(rq, zq, t))
      return "SS";
    return argon.deep(rq, zq, t) ? "LAr" : "EFCu";
  }
  int physAt(double rq, double zq) const
  {
    const char *m = materialAt(rq, zq);
    return !strcmp(m, "OFHC") ? 1 : (!strcmp(m, "SS") ? 2 : (!strcmp(m, "EFCu") ? 0 : -1));
  }                                                              // 0 mother, 1 OFHC, 2 SS
  double wallVolume(double z0, double z1, double dz = 0.5) const // [mm^3]: the tube's solid minus the argon in it
  {
    z0 = std::max(std::min(z0, z1), zBottom);
    z1 = std::min(std::max(z0, z1), zTop);
    double v = 0;
    for (double z = z0; z < z1; z += dz)
    {
      double s = std::min(dz, z1 - z);
      v += (wall.areaAt(z + 0.5 * s) - argon.areaAt(z + 0.5 * s)) * s;
    }
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
  for (int up = 0; up < 6; up++, p = "../" + p)
    if (!gSystem->AccessPathName(p.c_str()))
      return p;
  return "";
}

static RT loadRT(const std::string &gdml)
{
  RT rt;
  rt.file = findFile(gdml);
  if (rt.file.empty())
  {
    printf("ERROR: %s not found\n", gdml.c_str());
    return rt;
  }
  rt.wall = readOutline(rt.file, "reentrance_tube_copper");
  rt.argon = readOutline(rt.file, "liquid_argon_underground");
  rt.ofhcOuter = readOutline(rt.file, "reentrance_tube_layer_copper_ofhc_outer_bound");
  rt.ofhcInner = readOutline(rt.file, "reentrance_tube_layer_copper_ofhc_inner_bound");
  rt.ssOuter = readOutline(rt.file, "reentrance_tube_layer_steel_316L_outer_bound");
  rt.ssInner = readOutline(rt.file, "reentrance_tube_layer_steel_316L_inner_bound");
  if (!rt.wall.size() || !rt.argon.size() || !rt.ofhcOuter.size() || !rt.ofhcInner.size() || !rt.ssOuter.size() || !rt.ssInner.size())
  {
    printf("ERROR: %s has no re-entrant tube\n", rt.file.c_str());
    return rt;
  }
  auto rho = readDensities(rt.file);
  const char *densityOf[3] = {"reentrance_tube_layer_steel_316L", "reentrance_tube_layer_copper_ofhc", "reentrance_tube_copper"};
  for (int k = 0; k < 3; k++)
  {
    double d = rho.count(densityOf[k]) ? rho[densityOf[k]] : 0;
    if (d <= 0)
    {
      printf("ERROR: no density for %s in %s\n", densityOf[k], rt.file.c_str());
      return rt;
    }
    rt.density[k] = d;
  }
  rt.zBottom = rt.wall.zmin;
  rt.zTop = rt.wall.zmax;
  rt.seamOFHC = rt.ofhcOuter.zmin;
  rt.seamSS = rt.ssOuter.zmin;
  rt.zHead = rt.wall.z[std::max_element(rt.wall.r.begin(), rt.wall.r.end()) - rt.wall.r.begin()];
  rt.ok = true;
  return rt;
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

static std::string outPath(const std::string &name)
{
  gSystem->mkdir("output", kTRUE);
  return "output/" + name;
}

//-------------------------------------------------------------------------------
//  3. Check the tube: one clean surface, the endcap flush with the barrel, the shells on the barrel; and the densities read
static void checkTube(const RT &s)
{
  const Outline &w = s.wall;
  const std::vector<double> &r = w.r, &z = w.z;
  const int n = w.size();
  auto turn = [](double ar, double az, double br, double bz, double cr, double cz)
  { double v = (bz - az) * (cr - br) - (br - ar) * (cz - bz); return fabs(v) < 1e-9 ? 0 : (v > 0 ? 1 : 2); };
  int nAxis = 0, nDup = 0, nSpike = 0, nCross = 0;
  for (int i = 0; i < n; i++)
  {
    int i2 = (i + 1) % n;
    if (r[i] == 0)
      nAxis++; // only the bottom tip and the top sit on the axis
    for (int j = 0; j < i; j++)
      if (r[i] == r[j] && z[i] == z[j])
        nDup++; // the same point twice
    if (i > 0 && i < n - 1)
    {
      double ar = r[i] - r[i - 1], az = z[i] - z[i - 1], br = r[i + 1] - r[i], bz = z[i + 1] - z[i], m = sqrt((ar * ar + az * az) * (br * br + bz * bz));
      if (m > 0 && (ar * br + az * bz) / m < -0.99)
        nSpike++; // the outline folds straight back
    }
    for (int j = i + 2; j < n; j++) // segments that are not neighbours must not cross
    {
      int j2 = (j + 1) % n;
      if (i == 0 && j == n - 1)
        continue;
      if (turn(r[i], z[i], r[i2], z[i2], r[j], z[j]) != turn(r[i], z[i], r[i2], z[i2], r[j2], z[j2]) &&
          turn(r[j], z[j], r[j2], z[j2], r[i], z[i]) != turn(r[j], z[j], r[j2], z[j2], r[i2], z[i2]))
        nCross++;
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
  printf("    densities [kg/m^3], the GDML's: steel %.0f, Cu %.0f, EFCu %.0f\n", s.density[0], s.density[1], s.density[2]);
  bool ok = nAxis == 2 && !nDup && !nSpike && !nCross && fabs(step) <= tol && fabs(s.ofhcOuter.rmax - rBarrel) <= tol && fabs(s.ssOuter.rmax - rBarrel) <= tol;
  printf("    tube: %s\n\n", ok ? "PASS" : "FAIL");

  // the tube's shape in one graph: radius against height above its lowest point, each section in its colour, its seams marked
  // as depths from the top
  auto c = new TCanvas("c_rt", "", 700, 900);
  c->SetGrid();
  c->SetMargin(0.17, 0.22, 0.09, 0.07);
  const double hTop = s.zTop - s.zBottom, hLo = 0; // mm
  auto fr = c->DrawFrame(0, hLo, w.rmax * 1.08, hTop * 1.06);
  fr->SetTitle("the re-entrant tube, as built;radius [mm];height above the tube's bottom [mm]");
  fr->GetYaxis()->SetTitleOffset(1.9);
  fr->GetXaxis()->SetNdivisions(505);
  const double cut[4] = {s.zBottom - 1, s.seamOFHC, s.seamSS, s.zTop + 1};
  const char *secName[3] = {"EFCu", "OFHC Cu", "SS 316L"};
  const int secCol[3] = {TColor::GetColor("#c0392b"), TColor::GetColor("#e08e0b"), TColor::GetColor("#5d6d7e")};
  for (int i = 0; i + 1 < n; i++) // each edge of the outline, split at the seams: the barrel takes its shell's colour, the endcap and the lid EFCu's
    for (int k = 0; k < 3; k++)
    {
      double z0 = z[i], z1 = z[i + 1], a = 0, b = 1;
      if (z1 != z0)
      {
        double t0 = (cut[k] - z0) / (z1 - z0), t1 = (cut[k + 1] - z0) / (z1 - z0);
        a = std::max(0.0, std::min(t0, t1));
        b = std::min(1.0, std::max(t0, t1));
      }
      else if (z0 < cut[k] || z0 > cut[k + 1])
        continue;
      if (a >= b)
        continue;
      double pr[2], ph[2];
      for (int e = 0; e < 2; e++)
      {
        double t = e ? b : a;
        pr[e] = r[i] + t * (r[i + 1] - r[i]);
        ph[e] = std::max(hLo, z0 + t * (z1 - z0) - s.zBottom);
      }
      const int m = std::min(pr[0], pr[1]) >= rBarrel - 1 ? k : 0;
      auto g = new TGraph(2, pr, ph);
      g->SetLineColor(secCol[m]);
      g->SetLineWidth(4);
      g->Draw("L");
    }
  for (int k = 0; k < 3; k++)
  {
    const double hMid = 0.5 * (std::max(cut[k] - s.zBottom, 0.0) + std::min(cut[k + 1] - s.zBottom, hTop));
    auto t = new TLatex(w.rmax * 1.10, hMid, secName[k]);
    t->SetTextColor(secCol[k]);
    t->SetTextFont(62);
    t->SetTextSize(0.035);
    t->SetTextAlign(12);
    t->Draw();
  }
  auto lid = new TLatex(w.rmax * 0.5, hTop + 60, "EFCu lid");
  lid->SetTextColor(secCol[0]);
  lid->SetTextFont(62);
  lid->SetTextSize(0.03);
  lid->SetTextAlign(21);
  lid->Draw();
  for (double v : {s.seamOFHC, s.seamSS})
  {
    auto l = new TLine(0, v - s.zBottom, w.rmax * 1.08, v - s.zBottom);
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
    auto t = new TLatex(20, v - s.zBottom + 60, Form("seam, %.2f m from the top", (s.zTop - v) / 1000));
    t->SetTextSize(0.028);
    t->SetTextColor(kGray + 3);
    t->Draw();
  }
  c->SaveAs(outPath("tube.png").c_str());
}

//-------------------------------------------------------------------------------
//  4. One run, reduced to what a reweighting needs:
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
  std::string file, iso;
  Long64_t nsim = 0;
  long nFiles = 1, nMat[5] = {0, 0, 0, 0, 0};                                            // decays in EFCu, OFHC, SS, the argon, outside the wall
  long nPv[3] = {0, 0, 0};                                                               // decays per section: 0 mother (EFCu), 1 OFHC, 2 SS
  long nZPv[3][cfg::NZ] = {{0}};                                                         // the same, per ~1 cm of depth
  double mcDensity[3] = {0, 0, 0};                                                       // simulated decays per m^3 in each section, measured
  long dN[cfg::NDB] = {0}, dH[cfg::NDB] = {0}, dW1[cfg::NDB] = {0}, dWA[cfg::NDB] = {0}; // per depth slab: decays, hits, window hits after M1 and after all cuts
  double S = 0;                                                                          // the share of its M1 window hits that all cuts keep, near the detectors ([9])
  std::vector<float> depth;                                                              // per decay WITH a Ge hit: mm below the top of the tube
  std::vector<signed char> pv;                                                           // its section
  std::vector<short> win[cfg::NC];                                                       // its hits in the window, per cut level
  std::vector<short> winMap[2];                                                          // ...and with the LAr veto read off the map alone: M1 + LAr, + A/E
  std::vector<Hit> hits;
  long nHit[cfg::NC] = {0}, nLarMap = 0; // hits per cut level; M1 hits passing the LAr veto read off the map alone
};

static int slabOf(double dep, double L1, double L2) { return dep < L1 ? 0 : (dep < L2 ? 1 : 2); } // a design's three slabs, from the top
static int bin(double dep, double tot, int n) { return std::min(n - 1, std::max(0, (int)(dep / tot * n))); }

static double gPsdLow = cfg::psd_low; // the A/E cut in use: tuneAoE() sets it before any run is read

static void readHits(TFile *f, Long64_t n, std::vector<Hit> &hits, const std::string &file) // every decay's Ge hits, with the cuts decided
{
  std::vector<float> pe(n, 0), peMap(n, 0);
  std::vector<short> nFired(n, 0);
  std::vector<float> aoe;
  auto tg = (TTree *)f->Get("geds"), tl = (TTree *)f->Get("lar");
  if (!tg || !tl)
  {
    printf("ERROR: no geds or lar tree in %s\n", file.c_str());
    return;
  }
  scanTree(tg, {"evtid", "energy", "aoe_class"}, [&](const double *v) { // one row per detector hit, the response applied
    if (v[0] < 0 || v[0] >= n || v[1] <= cfg::m1_keV)
      return;
    nFired[(int)v[0]]++;
    hits.push_back({(Long64_t)v[0], (float)v[1], false, false, false, false});
    aoe.push_back(v[2]);
  });
  scanTree(tl, {"evtid", "pe", "pe_map"}, [&](const double *v)
           { if (v[0] >= 0 && v[0] < n) { pe[(int)v[0]] += v[1]; peMap[(int)v[0]] += v[2]; } });
  for (size_t i = 0; i < hits.size(); i++)
  {
    Hit &h = hits[i];
    h.m1 = nFired[h.ev] == 1;
    h.lar = pe[h.ev] < cfg::lar_pe;
    h.larMap = peMap[h.ev] < cfg::lar_pe;
    h.psd = aoe[i] > gPsdLow; // an undefined A/E (NaN) fails
  }
}

static Run readRun(const std::string &file, const std::string &iso, const RT &rt)
{
  Run r;
  r.file = file;
  r.iso = iso;
  auto f = TFile::Open(file.c_str());
  auto vtx = (f && !f->IsZombie()) ? (TTree *)f->Get("vtx") : nullptr;
  if (!vtx)
  {
    printf("ERROR: no vtx tree in %s\n", file.c_str());
    return r;
  }
  const Long64_t n = r.nsim = vtx->GetEntries();
  std::vector<float> depth(n, -1);
  std::vector<std::vector<short>> win(cfg::NC, std::vector<short>(n, 0)), winMap(2, std::vector<short>(n, 0));
  std::vector<signed char> pv(n, -1);
  readHits(f, n, r.hits, file);
  const char *matName[5] = {"EFCu", "OFHC", "SS", "LAr", "outside"};
  scanTree(vtx, {"evtid", "xloc", "yloc", "zloc"}, [&](const double *v) { // [m]
    int ev = (int)v[0];                                                   // rows are not in event order with -t: index by evtid
    if (ev < 0 || ev >= n)
      return;
    double x = v[1] * 1000, y = v[2] * 1000, z = v[3] * 1000;
    const char *m = rt.materialAt(sqrt(x * x + y * y), z);
    int k = 0;
    while (k < 4 && strcmp(m, matName[k]))
      k++;
    r.nMat[k]++;
    depth[ev] = rt.zTop - z;
    pv[ev] = k < 3 ? k : -1; // EFCu is the mother (0), OFHC 1, SS 2
  });
  const double tot = rt.zTop - rt.zBottom;
  for (auto &h : r.hits)
  {
    for (int c = 0; c < cfg::NC; c++)
    {
      r.nHit[c] += h.pass(c);
      if (h.pass(c) && cfg::inWindow(h.e))
        win[c][h.ev]++;
    }
    r.nLarMap += h.m1 && h.larMap;
    for (int k = 0; k < 2; k++)
      if (h.passMap(k) && cfg::inWindow(h.e))
        winMap[k][h.ev]++;
  }
  for (Long64_t i = 0; i < n; i++) // every decay, as counts
  {
    r.dN[bin(depth[i], tot, cfg::NDB)]++;
    r.dW1[bin(depth[i], tot, cfg::NDB)] += win[1][i];
    r.dWA[bin(depth[i], tot, cfg::NDB)] += win[cfg::NC - 1][i];
    if (pv[i] >= 0)
    {
      r.nPv[pv[i]]++;
      r.nZPv[pv[i]][bin(depth[i], tot, cfg::NZ)]++;
    }
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
      for (int c = 0; c < cfg::NC; c++)
        r.win[c].push_back(win[c][h.ev]);
      for (int k = 0; k < 2; k++)
        r.winMap[k].push_back(winMap[k][h.ev]);
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

// the A/E cut as LEGEND sets it: the classifier threshold that keeps cfg::depKeep of the DEP of Tl208, in M1 hits of every
// Tl208 run of the tube, the continuum under the peak subtracted with the sidebands. one value for all detectors: they
// share one template and one set of parameters. then the cut on the other lines a 228Th calibration measures it on, and on
// the window's continuum, against the HADES characterization of the L1000 template detector V00000A
static void tuneAoE(const std::vector<std::string> &files)
{
  struct Line
  {
    const char *name;
    double e;
    std::vector<float> peak, side;
    double measured;
  }; // A/E classifier of M1 hits in the peak, the sidebands
  Line line[3] = {{"DEP", cfg::depE, {}, {}, 0.90}, {"SEP", 2103.5, {}, {}, 0.04}, {"FEP", 2614.5, {}, {}, 0.06}};
  std::vector<float> win; // ...and in the window
  for (auto &file : files)
  {
    auto f = TFile::Open(file.c_str());
    auto t = (f && !f->IsZombie()) ? (TTree *)f->Get("geds") : nullptr;
    if (!t)
    {
      delete f;
      continue;
    }
    std::unordered_map<Long64_t, short> nFired;
    scanTree(t, {"evtid", "energy"}, [&](const double *v)
             { if (v[1] > cfg::m1_keV) nFired[(Long64_t)v[0]]++; });
    scanTree(t, {"evtid", "energy", "aoe_class"}, [&](const double *v)
             {
      if (v[1] <= cfg::m1_keV || nFired[(Long64_t)v[0]] != 1) return;
      float a = std::isnan(v[2]) ? -1e30 : v[2]; // an undefined A/E fails
      if (cfg::inWindow(v[1])) win.push_back(a);
      for (auto &l : line)
      {
        double d = fabs(v[1] - l.e);
        if (d < cfg::depHalf) l.peak.push_back(a);
        else if (d > cfg::sideLo && d < cfg::sideHi) l.side.push_back(a);
      } });
    delete f;
  }
  const double k = 2 * cfg::depHalf / (2 * (cfg::sideHi - cfg::sideLo)); // sideband counts per peak width
  auto above = [](std::vector<float> &x, double thr)
  { return (double)(x.end() - std::upper_bound(x.begin(), x.end(), (float)thr)); };
  for (auto &l : line)
  {
    std::sort(l.peak.begin(), l.peak.end());
    std::sort(l.side.begin(), l.side.end());
  }
  std::sort(win.begin(), win.end());
  auto net = [&](Line &l)
  { return l.peak.size() - k * l.side.size(); };
  auto kept = [&](Line &l, double thr)
  { return (above(l.peak, thr) - k * above(l.side, thr)) / net(l); }; // net fraction above thr
  Line &dep = line[0];
  if (net(dep) < 30)
  {
    printf("[5] A/E cut: only %.0f DEP events net, kept at the production's %.2f\n\n", net(dep), cfg::psd_low);
    return;
  }
  double thr = 0;
  while (thr > -50 && kept(dep, thr) < cfg::depKeep)
    thr -= 0.005;
  gPsdLow = thr;
  printf("[5] A/E cut tuned on the Tl208 DEP (%.1f keV) in M1 hits: %zu in +-%.1f keV, %zu in the sidebands %.0f-%.0f keV off it, %.0f net\n",
         cfg::depE, dep.peak.size(), cfg::depHalf, dep.side.size(), cfg::sideLo, cfg::sideHi, net(dep));
  printf("    classifier > %.2f keeps %.0f%% of it (+- %.0f%% statistics); the production's > %.2f kept %.0f%%\n", gPsdLow, 100 * kept(dep, gPsdLow),
         100 * sqrt(cfg::depKeep * (1 - cfg::depKeep) / net(dep)), cfg::psd_low, 100 * kept(dep, cfg::psd_low));
  printf("    kept by it, net of the continuum:");
  for (auto &l : line)
    printf("  %s %.1f %.2f (%.0f net)", l.name, l.e, kept(l, gPsdLow), net(l));
  printf("  window %.2f (%zu)\n    measured on V00000A at HADES:     ", win.empty() ? 0.0 : above(win, gPsdLow) / win.size(), win.size());
  for (auto &l : line)
    printf("  %s %.1f %.2f          ", l.name, l.e, l.measured);
  printf("  window 0.27\n\n");
}

void background(const char *runs = "output/tl208*_hit.root=Tl208,output/bi214*_hit.root=Bi214", const char *gdml = "geom/l1000.gdml")
{
  gErrorIgnoreLevel = kWarning; // no "png file has been created" notices in the printout
  RT rt = loadRT(gdml);
  if (!rt.ok)
    return;
  checkTube(rt);
  std::vector<std::string> tlFiles; // every Tl208 file of the tube, for the A/E cut
  std::stringstream tls(runs);
  for (std::string tok; std::getline(tls, tok, ',');)
  {
    size_t eq = tok.find('=');
    if ((eq == std::string::npos ? "Tl208" : tok.substr(eq + 1)).find("Tl") != std::string::npos)
      for (auto &f : expand(tok.substr(0, eq)))
        tlFiles.push_back(f);
  }
  tuneAoE(tlFiles);
  std::vector<Run> R; // one per isotope
  std::stringstream list(runs);
  for (std::string tok; std::getline(list, tok, ',');) // "a_hit.root=Tl208,b_hit.root=Bi214"
  {
    size_t eq = tok.find('=');
    std::string iso = eq == std::string::npos ? "Tl208" : tok.substr(eq + 1);
    if (iso.find("Tl") == std::string::npos && iso.find("Bi") == std::string::npos)
    {
      printf("ERROR: unknown isotope '%s'\n", iso.c_str());
      continue;
    }
    const std::string pat = tok.substr(0, eq);
    for (auto &f : expand(pat))
    {
      Run r = readRun(f, iso, rt);
      if (!r.nsim)
        continue;
      r.file = pat; // the printout names the pattern, not every file it matched
      Run *into = nullptr;
      for (auto &m : R)
        if (m.iso == iso)
          into = &m;
      if (!into)
      {
        R.push_back(std::move(r));
        continue;
      }
      Long64_t off = into->depth.size(); // slots renumbered, so nothing counts twice
      for (auto h : r.hits)
      {
        h.ev += off;
        into->hits.push_back(h);
      }
      into->depth.insert(into->depth.end(), r.depth.begin(), r.depth.end());
      into->pv.insert(into->pv.end(), r.pv.begin(), r.pv.end());
      for (int c = 0; c < cfg::NC; c++)
      {
        into->win[c].insert(into->win[c].end(), r.win[c].begin(), r.win[c].end());
        into->nHit[c] += r.nHit[c];
      }
      into->nLarMap += r.nLarMap;
      for (int k = 0; k < 2; k++)
        into->winMap[k].insert(into->winMap[k].end(), r.winMap[k].begin(), r.winMap[k].end());
      for (int k = 0; k < 5; k++)
        into->nMat[k] += r.nMat[k];
      for (int p = 0; p < 3; p++)
      {
        into->nPv[p] += r.nPv[p];
        for (int b = 0; b < cfg::NZ; b++)
          into->nZPv[p][b] += r.nZPv[p][b];
      }
      for (int b = 0; b < cfg::NDB; b++)
      {
        into->dN[b] += r.dN[b];
        into->dH[b] += r.dH[b];
        into->dW1[b] += r.dW1[b];
        into->dWA[b] += r.dWA[b];
      }
      into->nsim += r.nsim;
      into->nFiles++;
      if (into->file != pat && into->file.find("+" + pat) == std::string::npos)
        into->file += "+" + pat;
    }
  }
  if (R.empty())
  {
    printf("no usable runs in '%s'\n", runs);
    return;
  }

  const double depthTot = rt.zTop - rt.zBottom, L1built = rt.zTop - rt.seamSS, L2built = rt.zTop - rt.seamOFHC; // the seams as built, as depths
  const double UL90 = 2.30, binW = (cfg::eHi - cfg::eLo) / cfg::eBins;
  const int NC = cfg::NC, ALL = NC - 1;
  const char *const *cutName = cfg::cutName;
  double physV[3]; // m^3: mother, OFHC shell, SS shell
  for (int p = 0; p < 3; p++)
    physV[p] = rt.physVolume(p) * 1e-9;
  for (auto &r : R)
    for (int p = 0; p < 3; p++)
      r.mcDensity[p] = r.nPv[p] / physV[p]; // remage fills the mother sparser than its daughters: measure, never assume
  // one simulated decay in section p, where material m sits = rho x A x 1 yr real decays per m^3, over the MC density, per kg of Ge
  auto weight = [&](const Run &r, int m, int p)
  { return p < 0 || r.mcDensity[p] <= 0 ? 0.0 : rt.density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[p] * cfg::geMass_kg); };
  auto massOf = [&](double d0, double d1, int m)
  { return d1 > d0 ? rt.density[m] * rt.wallVolume(rt.zTop - d1, rt.zTop - d0) * 1e-9 : 0.0; };
  auto ulN = [](long n)
  { const double t[11] = {2.30, 3.89, 5.32, 6.68, 7.99, 9.27, 10.53, 11.77, 12.99, 14.21, 15.41}; return n <= 10 ? t[n] : n + 1.28 * sqrt((double)n) + 1; }; // Poisson 90% upper limit
  printf("window   : %.0f-%.0f keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's BEW): %.0f keV\n", cfg::winLo, cfg::winHi, cfg::winWidth);
  printf("cuts     : M1 (exactly one detector above %.0f keV), LAr veto (the SiPMs see < %.0f photoelectrons), A/E classifier > %.2f (%.0f%% DEP, [5]). energies: the detector response\n",
         cfg::m1_keV, cfg::lar_pe, gPsdLow, 100 * cfg::depKeep);
  printf("activity : chain [uBq/kg] 232Th / 238U: steel %g / %g, Cu %g / %g, EFCu %g / %g; Tl208 is %.2f%% of 232Th, Bi214 %.0f%% of 238U\n\n",
         cfg::chain[0][0], cfg::chain[0][1], cfg::chain[1][0], cfg::chain[1][1], cfg::chain[2][0], cfg::chain[2][1], 100 * cfg::branch[0], 100 * cfg::branch[1]);

  for (auto &r : R)
  {
    //-------------------------------------------------------------------------------
    //  6. Per run: the decays, the hits, the normalisation
    long w[NC] = {0};
    for (int c = 0; c < NC; c++)
      for (auto k : r.win[c])
        w[c] += k;
    printf("=== %s   (%s, %lld decays)\n", r.iso.c_str(), r.file.c_str(), r.nsim);
    printf("[6] decays: EFCu %.1f%%, OFHC %.1f%%, SS %.1f%%; in the argon %ld, outside the wall %ld%s%s\n", 100.0 * r.nMat[0] / r.nsim, 100.0 * r.nMat[1] / r.nsim,
           100.0 * r.nMat[2] / r.nsim, r.nMat[3], r.nMat[4], r.nMat[3] || r.nMat[4] ? "  <- CONFINEMENT BUG" : "", r.nFiles > 1 ? Form("   (merged from %ld files)", r.nFiles) : "");
    printf("    hits: %ld, M1 %ld, + LAr %ld, + A/E %ld; in the window: %ld, %ld, %ld, %ld\n", r.nHit[0], r.nHit[1], r.nHit[2], r.nHit[3], w[0], w[1], w[2], w[3]);
    printf("    LAr veto keeps 1 in %.0f M1 hits; read off the optical map alone, without the argon beyond its edge, 1 in %.0f\n",
           r.nHit[2] ? (double)r.nHit[1] / r.nHit[2] : 0.0, r.nLarMap ? (double)r.nHit[1] / r.nLarMap : 0.0);
    const int asBuilt[3] = {2, 1, 0}; // mother EFCu, OFHC shell Cu, SS shell steel
    const char *pvName[3] = {"mother", "OFHC", "SS"};
    printf("    %-7s %-6s %9s %8s %10s %9s %11s %14s\n", "section", "mat", "V [m^3]", "M [kg]", "A [Bq]", "MC dec.", "MC per m^3", "decays/yr per");
    for (int p = 0; p < 3; p++)
    {
      int m = asBuilt[p];
      double mass = physV[p] * rt.density[m], A = mass * cfg::activity(r.iso, m);
      printf("    %-7s %-6s %9.5f %8.1f %10.3e %9ld %11.0f %14.3g\n", pvName[p], cfg::mat[m], physV[p], mass, A, r.nPv[p], r.mcDensity[p], r.nPv[p] ? A * cfg::secPerYear / r.nPv[p] : 0.0);
    }
    double dAvg = 0.5 * (r.mcDensity[1] + r.mcDensity[2]);
    printf("    sampling density mother / shells = %.3f   (1.000 would be uniform; the weights use the measured density)\n", dAvg > 0 ? r.mcDensity[0] / dAvg : 0.0);

    //-------------------------------------------------------------------------------
    //  7. Per run: where the hits come from, by depth; and the share of the M1 window hits all cuts keep there
    printf("[7] depth [m]          decays   hits/decay   window: M1   all cuts   kept\n");
    for (int b = 0; b < cfg::NDB; b++)
      printf("    %5.2f..%-8.2f %10ld %12.2e %12ld %10ld %6s%s\n", depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000, r.dN[b],
             r.dN[b] ? (double)r.dH[b] / r.dN[b] : 0.0, r.dW1[b], r.dWA[b], r.dW1[b] ? Form("%.4f", (double)r.dWA[b] / r.dW1[b]) : "-",
             b == cfg::NDB - 1 ? "   <- nearest the detectors" : "");
    printf("\n");
  }

  //-------------------------------------------------------------------------------
  //  8. Background index as built, per chain and section:
  printf("[8] background index as built [cts/(keV kg yr)], +- MC statistics; a section with no window hit gets its 90%% limit, 2.30 x its heaviest decay\n");
  printf("    %-6s %-6s %15s   %-22s %-22s %-22s %-22s\n", "chain", "sect", "win hits", cutName[0], cutName[1], cutName[2], cutName[3]);
  double gB[NC] = {0}, gV[NC] = {0};
  long gN[NC] = {0};               // window hits behind each total
  double mB[2] = {0}, mV[2] = {0}; // the same after M1 + LAr and all cuts, the LAr veto off the map alone (steel and empty sections add nothing)
  long mN[2] = {0};
  auto total = [](double b, double v, long n, bool goal) -> std::string
  {
    if (!n)
      return "no window hit";
    return goal ? Form("%.2e +- %.1e (%.2g x goal)", b, sqrt(v), b / cfg::bgGoal) : Form("%.2e +- %.1e", b, sqrt(v));
  };
  const double edgeBuilt[4] = {0, L1built, L2built, depthTot};
  for (auto &r : R)
  {
    double tB[NC] = {0}, tV[NC] = {0};
    long tN[NC] = {0};
    for (int s = 0; s < 3; s++) // the slabs as built: steel, Cu, EFCu
    {
      double B[NC] = {0}, V[NC] = {0}, wmax = 0;
      long n[NC] = {0};
      for (int p = 0; p < 3; p++) // the heaviest decay in the slab, for the limit: runs confined to one section sample a slab unevenly
        for (int b = bin(edgeBuilt[s], depthTot, cfg::NZ); b < std::min(cfg::NZ, (int)ceil(edgeBuilt[s + 1] / depthTot * cfg::NZ)); b++)
          if (r.nZPv[p][b])
          {
            wmax = std::max(wmax, weight(r, s, p) / cfg::winWidth);
            break;
          }
      for (size_t i = 0; i < r.depth.size(); i++)
      {
        if (slabOf(r.depth[i], L1built, L2built) != s)
          continue;
        double wt = weight(r, s, r.pv[i]) / cfg::winWidth;
        for (int c = 0; c < NC; c++)
        {
          B[c] += r.win[c][i] * wt;
          V[c] += r.win[c][i] * wt * wt;
          n[c] += r.win[c][i];
        }
        for (int k = 0; k < 2; k++)
        {
          mB[k] += r.winMap[k][i] * wt;
          mV[k] += r.winMap[k][i] * wt * wt;
          mN[k] += r.winMap[k][i];
        }
      }
      std::string cell[NC];
      for (int c = 0; c < NC; c++)
      {
        cell[c] = n[c] ? Form("%.2e +- %.1e", B[c], sqrt(V[c])) : Form("< %.2e (90%%)", UL90 * wmax);
        if (n[c])
        {
          tB[c] += B[c];
          tV[c] += V[c];
          tN[c] += n[c];
        }
      }
      printf("    %-6s %-6s %15s   %-22s %-22s %-22s %-22s\n", s ? "" : r.iso.c_str(), cfg::mat[s], Form("%ld/%ld/%ld/%ld", n[0], n[1], n[2], n[3]), cell[0].c_str(),
             cell[1].c_str(), cell[2].c_str(), cell[3].c_str());
    }
    printf("    %-6s %-6s %15s   %-22s %-22s %-22s %-22s\n", "", "tube", "", total(tB[0], tV[0], tN[0], false).c_str(), total(tB[1], tV[1], tN[1], false).c_str(),
           total(tB[2], tV[2], tN[2], false).c_str(), total(tB[3], tV[3], tN[3], false).c_str());
    for (int c = 0; c < NC; c++)
    {
      gB[c] += tB[c];
      gV[c] += tV[c];
      gN[c] += tN[c];
    }
  }
  printf("    ALL CHAINS\n");
  for (int c = 0; c < NC; c++)
    printf("      %-16s %s\n", cutName[c], total(gB[c], gV[c], gN[c], true).c_str());
  printf("    ...with the LAr veto read off the optical map alone, no light from the argon beyond its edge:\n");
  for (int k = 0; k < 2; k++)
    printf("      %-16s %s\n", cutName[2 + k], total(mB[k], mV[k], mN[k], true).c_str());
  printf("    win hits: per cut level, as in the columns. sections with no window hit are left out of the totals\n");

  //-------------------------------------------------------------------------------
  //  9. Designs against the tube as built: steel to L1, Cu to L2, EFCu below, for less EFCu (L2 deeper) and more steel (L1 deeper)
  // a design differs from the tube as built only where a slab changed material, so its BI after all cuts is the as-built's plus
  // that change, dB. dB is counted on the changed slabs' M1 window hits, which outnumber those all cuts keep 100-500 to 1, times
  // S, the share of them all cuts keep, measured where the change is ([7], kept): near the detectors, below farEdge, per chain;
  // above it, where the cuts keep more (the optical map gives no light above z 1.375 m) and few hits are, both chains together.
  // its 90% bound adds, per chain, change (Cu or EFCu to steel, EFCu to Cu) and region, UL(n) x S x the heaviest weight change
  // there, with S's own 90% bound above farEdge: a slab with no hit still bounds it
  const double biBuilt = gB[ALL], maxRise = cfg::simTol * biBuilt;
  const int nFar = 7; // the slabs above 4.37 m: above the array
  const double farEdge = depthTot * nFar / cfg::NDB;
  long fA = 0, fM = 0;
  for (auto &r : R)
  {
    long a = 0, m = 0;
    for (int b = 0; b < cfg::NDB; b++)
      (b < nFar ? fA : a) += r.dWA[b], (b < nFar ? fM : m) += r.dW1[b];
    r.S = m ? (double)a / m : 0;
  }
  const double sFar = fM ? (double)fA / fM : 0, sFarUp = fM ? ulN(fA) / fM : 0;
  auto change = [&](double L1, double L2, double &up)
  {
    double dB = 0;
    up = 0;
    for (auto &r : R)
    {
      long n[3][3][2] = {};
      double wmax[3][3][2] = {};  // [as built][design] material (0 steel, 1 Cu, 2 EFCu), [near, far]
      for (int p = 0; p < 3; p++) // the sections sampled at the depths each change covers, ~1 cm at a time
        for (int b = 0; b < cfg::NZ; b++)
        {
          if (!r.nZPv[p][b])
            continue;
          const double d = (b + 0.5) * depthTot / cfg::NZ;
          const int mb = slabOf(d, L1built, L2built), md = slabOf(d, L1, L2), f = d < farEdge;
          if (mb != md)
            wmax[mb][md][f] = std::max(wmax[mb][md][f], (weight(r, md, p) - weight(r, mb, p)) / cfg::winWidth);
        }
      for (size_t i = 0; i < r.depth.size(); i++)
      {
        if (!r.win[1][i])
          continue;
        const int mb = slabOf(r.depth[i], L1built, L2built), md = slabOf(r.depth[i], L1, L2), f = r.depth[i] < farEdge;
        if (mb == md)
          continue;
        dB += (f ? sFar : r.S) * r.win[1][i] * (weight(r, md, r.pv[i]) - weight(r, mb, r.pv[i])) / cfg::winWidth;
        n[mb][md][f] += r.win[1][i];
      }
      for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
          for (int f = 0; f < 2; f++)
            if (wmax[a][b][f] > 0)
              up += (f ? sFarUp : r.S) * ulN(n[a][b][f]) * wmax[a][b][f];
    }
    return dB;
  };
  auto deepestSteel = [&](double L2, bool bound) { // 1 cm at a time below the steel as built, always some Cu
    double l1 = -1, up;
    for (double L1 = L1built; L1 <= L2 - 10; L1 += 10)
    {
      double dB = change(L1, L2, up);
      if ((bound ? up : dB) > maxRise)
        break;
      l1 = L1;
    }
    return l1;
  };
  auto split = [&](double l1, double L2)
  {
    double ms = massOf(0, l1, 0), mc = massOf(l1, L2, 1), me = massOf(L2, depthTot, 2), mt = ms + mc + me;
    return std::string(Form("%4.0f : %3.0f : %3.0f kg = %2.0f : %2.0f : %2.0f %%", ms, mc, me, 100 * ms / mt, 100 * mc / mt, 100 * me / mt));
  };
  printf("\n[9] designs against the tube as built, BI %.2e after all cuts: steel to L1, Cu to L2, EFCu below. a design passes when its BI\n", biBuilt);
  printf("    rises by at most %.0f%%: at 90%% confidence (bound), or on the central value (central). L1 from %.2f m down, 1 cm at a time\n", 100 * cfg::simTol, L1built / 1000);
  printf("    the change counts M1 window hits x S, the share all cuts keep: below %.2f m", farEdge / 1000);
  for (auto &r : R)
    printf(" %s %.2e", r.iso.c_str(), r.S);
  printf("; above it, both chains, %ld of %ld = %.2e (90%% bound %.2e)\n", fA, fM, sFar, sFarUp);
  printf("    %-6s %6s | %-9s %-34s %8s | %-9s %-34s %8s\n", "L2 [m]", "EFCu", "bound: L1", "SS : Cu : EFCu", "rise", "central L1", "SS : Cu : EFCu", "rise");
  std::vector<double> rows = {L2built};
  for (double L2 = 4100; L2 <= 5000; L2 += 100)
    if (L2 > L2built)
      rows.push_back(L2);
  for (double L2 = 5250; L2 < depthTot; L2 += 250)
    rows.push_back(L2);
  double best[2][2] = {{-1, -1}, {-1, -1}}; // [bound, central]: L1, L2 of the least EFCu, then the most steel
  for (double L2 : rows)
  {
    std::string col[2];
    for (int k = 0; k < 2; k++)
    {
      double l1 = deepestSteel(L2, k == 0), up;
      if (l1 < 0)
      {
        col[k] = Form("%-9s %-34s %8s", "-", "rises too much at any L1", "");
        continue;
      }
      double dB = change(l1, L2, up);
      col[k] = Form("%7.2f m %-34s %+7.0f%%", l1 / 1000, split(l1, L2).c_str(), 100 * (k ? dB : up) / biBuilt);
      if (L2 > best[k][1] || (L2 == best[k][1] && l1 > best[k][0]))
      {
        best[k][0] = l1;
        best[k][1] = L2;
      }
    }
    printf("    %-6.2f %3.0f kg | %s | %s%s\n", L2 / 1000, massOf(L2, depthTot, 2), col[0].c_str(), col[1].c_str(), L2 == L2built ? "   <- EFCu as built" : "");
  }
  for (int k = 0; k < 2; k++)
    if (best[k][0] > 0)
      printf("    least EFCu, then most steel, %s: L2 %.2f m, L1 %.2f m: %s\n", k ? "on the central value" : "at 90% confidence", best[k][1] / 1000,
             best[k][0] / 1000, split(best[k][0], best[k][1]).c_str());
  printf("    as built: %s. rise: the 90%% bound on dB (bound) or dB (central), over the BI as built\n", split(L1built, L2built).c_str());
  // the statistics a bound needs: a slab with no M1 window hit bounds steel there at UL(0) x S x its weight change, which falls
  // as more decays are simulated in it
  printf("    a slab with no M1 window hit bounds what steel there adds (x S above %.2f m); to bring that under %.0f%% of the BI as built per chain:\n", farEdge / 1000, 50 * cfg::simTol);
  for (auto &r : R)
    for (int p : {1, 2}) // steel added on the OFHC shell's depths; the steel as built, on the SS shell's
    {
      const double dw = (weight(r, 0, p) - (p == 1 ? weight(r, 1, p) : 0.0)) / cfg::winWidth, b = UL90 * sFar * dw;
      printf("      %-6s %-30s %10ld decays now: bound %5.0f%%; %.1e decays would bring it to %.0f%%\n", r.iso.c_str(), p == 1 ? "steel added in the OFHC shell" : "the steel as built, SS shell",
             r.nPv[p], 100 * b / biBuilt, r.nPv[p] * b / (0.5 * maxRise), 50 * cfg::simTol);
    }

  //-------------------------------------------------------------------------------
  //  10. Draw: the tube's spectrum as built, and every design's BI over the tube's as built
  TH1D *sp[NC];
  for (int c = 0; c < NC; c++)
    sp[c] = new TH1D(Form("sp%d", c), "Re-entrant tube background, as built;energy [keV];cts / (keV kg yr)", cfg::eBins, cfg::eLo, cfg::eHi);
  for (auto &r : R)
    for (auto &h : r.hits)
    {
      double wt = weight(r, slabOf(r.depth[h.ev], L1built, L2built), r.pv[h.ev]) / binW;
      for (int c = 0; c < NC; c++)
        if (h.pass(c))
          sp[c]->Fill(h.e, wt);
    }
  const int NS = 40;
  auto hMap = new TH2D("hMap", "BI after all cuts / BI as built (central);L1: steel down to [m];L2: Cu down to, EFCu below [m]", NS, 0, depthTot / 1000, NS, 0, depthTot / 1000);
  for (int i = 1; i <= NS; i++)
    for (int j = i; j <= NS; j++) // only L1 <= L2 is a design
    {
      double up, dB = change(hMap->GetXaxis()->GetBinCenter(i) * 1000, hMap->GetYaxis()->GetBinCenter(j) * 1000, up);
      hMap->SetBinContent(i, j, std::max(1e-3, (biBuilt + dB) / biBuilt));
    }
  gStyle->SetOptStat(0);
  auto c = new TCanvas("c_bkg", "", 1400, 560);
  c->Divide(2, 1);
  c->cd(1)->SetLogy();
  gPad->SetGrid();
  gPad->SetLeftMargin(0.14);
  gPad->SetTopMargin(0.15);
  auto leg = new TLegend(0.14, 0.86, 0.90, 0.91); // one row, between the title and the frame, clear of the ROI's label
  leg->SetNColumns(NC);
  leg->SetBorderSize(0);
  const int col[NC] = {kAzure + 2, kOrange + 7, kRed + 1, kViolet + 1};
  for (int k = 0; k < NC; k++)
  {
    sp[k]->SetLineColor(col[k]);
    sp[k]->SetLineWidth(2);
    sp[k]->Draw(k ? "HIST SAME" : "HIST");
    leg->AddEntry(sp[k], cutName[k], "l");
  }
  sp[0]->SetMinimum(0.1 * std::max(1e-12, sp[ALL]->GetMinimum(0) > 0 ? sp[ALL]->GetMinimum(0) : sp[2]->GetMinimum(0)));
  gPad->Update();
  { // the ROI over the spectrum
    const double y0 = pow(10, gPad->GetUymin()), y1 = pow(10, gPad->GetUymax());
    const int roiCol = TColor::GetColor("#2e7d32");
    auto b = new TBox(cfg::Qbb - cfg::roiHalf, y0, cfg::Qbb + cfg::roiHalf, y1);
    b->SetFillColorAlpha(roiCol, 0.8);
    b->SetLineWidth(0);
    b->Draw();
    auto t = new TLatex(cfg::Qbb + 15, pow(10, gPad->GetUymin() + 0.9 * (gPad->GetUymax() - gPad->GetUymin())),
                        Form("#splitline{ROI Q_{#beta#beta} #pm 2#sigma}{%.1f-%.1f keV}", cfg::Qbb - cfg::roiHalf, cfg::Qbb + cfg::roiHalf));
    t->SetTextAlign(13);
    t->SetTextColor(roiCol);
    t->SetTextFont(62);
    t->SetTextSize(0.032);
    t->Draw(); // left-aligned: centred text misplaces a subscript
  }
  leg->Draw();
  c->cd(2)->SetLogz();
  gPad->SetRightMargin(0.15);
  hMap->SetMinimum(std::min(0.5, hMap->GetMinimum(0)));
  hMap->GetZaxis()->SetMoreLogLabels();
  hMap->Draw("COLZ");
  for (int k = 0; k < 2; k++) // the deepest steel that passes at each L2: on the central value (solid), at 90% confidence (dashed)
  {
    auto g = new TGraph();
    for (double L2 = L2built; L2 < depthTot; L2 += depthTot / NS)
    {
      double l1 = deepestSteel(L2, k == 1);
      if (l1 >= 0)
        g->AddPoint(l1 / 1000, L2 / 1000);
    }
    g->SetLineColor(kBlack);
    g->SetLineWidth(2);
    g->SetLineStyle(k ? 2 : 1);
    if (g->GetN() > 1)
      g->Draw("L");
  }
  auto star = new TMarker(L1built / 1000, L2built / 1000, 29);
  star->SetMarkerColor(kRed + 1);
  star->SetMarkerSize(2.2);
  star->Draw();
  for (int k = 0; k < 2; k++) // the chosen designs: at 90% confidence (circle), on the central value (square)
    if (best[k][0] > 0)
    {
      auto m = new TMarker(best[k][0] / 1000, best[k][1] / 1000, k ? 25 : 24);
      m->SetMarkerColor(kBlack);
      m->SetMarkerSize(1.8);
      m->SetMarkerStyle(k ? 25 : 24);
      m->Draw();
    }
  std::string png = outPath("background.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s and output/tube.png\n", png.c_str());
}
