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
  // the band shaded on the spectrum: Qbb +- 2 sigma, sigma from the response's resolution (simprod eresmod: FWHM = sqrt(0.5 + 0.001 E) keV)
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
  const double bgGoal = 1e-5;      // LEGEND-1000 goal [cts/(keV kg yr)] at Qbb, for scale
  const double simTol = 0.10;      // a design passes when its BI after all cuts is at most this much above the tube's as built: the BI's own statistical uncertainty
  const double secPerYear = 365.25 * 24 * 3600;
  const double eLo = 1000, eHi = 3000; // spectrum [keV]
  const int eBins = 200;
  const int NDB = 10; // depth slabs of the hit table
  const int NZ = 625; // ~1 cm depth bins: where each section's decays lie

  const char *mat[3] = {"SS", "Cu", "EFCu"}; // top to bottom; their densities are the GDML's (RT::density)
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
  std::vector<std::vector<int>> bins; // per height bin, the edges (i-1, i) that span it: contains() tests only those
  double bz0 = 0, bdz = 0;
  void index(int nb = 1024)
  {
    bins.assign(nb, {});
    bz0 = zmin;
    bdz = (zmax - zmin) / nb;
    for (int i = 0, j = size() - 1; i < size(); j = i++)
    {
      int b0 = std::max(0, (int)((std::min(z[i], z[j]) - bz0) / bdz)), b1 = std::min(nb - 1, (int)((std::max(z[i], z[j]) - bz0) / bdz));
      for (int b = b0; b <= b1; b++)
        bins[b].push_back(i);
    }
  }
  bool contains(double rq, double zq) const // even-odd ray cast, over the edges at that height
  {
    if (zq < zmin || zq >= zmax || bins.empty())
      return false;
    bool in = false;
    for (int i : bins[std::min((int)bins.size() - 1, (int)((zq - bz0) / bdz))])
    {
      int j = i ? i - 1 : size() - 1;
      if ((z[i] > zq) != (z[j] > zq) && rq < r[i] + (zq - z[i]) / (z[j] - z[i]) * (r[j] - r[i]))
        in = !in;
    }
    return in;
  }
  // inside or within tol of it (near), at least tol inside (deep). nudged in r AND z: on a shallow cone a float32 z error of 0.2 um is 40 um in r
  bool near(double rq, double zq, double t) const { return contains(rq, zq) || contains(rq + t, zq) || contains(fabs(rq - t), zq) || contains(rq, zq + t) || contains(rq, zq - t); }
  bool deep(double rq, double zq, double t) const { return contains(rq, zq) && contains(rq + t, zq) && contains(fabs(rq - t), zq) && contains(rq, zq + t) && contains(rq, zq - t); }
  double areaAt(double zq) const // cross-section [mm^2]: the crossings pair into rings, so a hollow solid counts right
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
    o.index();
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

static double readGermanium(const std::string &file, int &n) // [kg]: every germanium detector's solid x its density x its placements, out of the GDML text
{
  std::map<std::string, double> material;                                   // kg/m^3
  std::map<std::string, Outline> solid;                                     // every <genericPolycone>
  std::map<std::string, std::pair<std::string, std::string>> volume;        // logical volume -> its material, its solid
  std::map<std::string, int> placed;                                        // logical volume -> its placements
  std::ifstream in(file.c_str());
  auto attr = [](const std::string &line, const char *key)
  { size_t a = line.find(key); if (a == std::string::npos) return std::string(); a += strlen(key); return line.substr(a, line.find('"', a) - a); };
  std::string mat, vol;
  Outline *cur = nullptr;
  for (std::string line; std::getline(in, line);)
  {
    if (cur)
    {
      double rr, zz;
      if (line.find("</genericPolycone>") != std::string::npos)
      {
        cur->rmin = *std::min_element(cur->r.begin(), cur->r.end());
        cur->rmax = *std::max_element(cur->r.begin(), cur->r.end());
        cur->zmin = *std::min_element(cur->z.begin(), cur->z.end());
        cur->zmax = *std::max_element(cur->z.begin(), cur->z.end());
        cur = nullptr;
      }
      else if (sscanf(line.c_str(), " <rzpoint r=\"%lf\" z=\"%lf\"", &rr, &zz) == 2)
      {
        cur->r.push_back(rr);
        cur->z.push_back(zz);
      }
      continue;
    }
    if (line.find("<genericPolycone") != std::string::npos)
      cur = &solid[attr(line, "name=\"")];
    else if (line.find("<material ") != std::string::npos)
      mat = attr(line, "name=\"");
    else if (!mat.empty() && line.find("<D ") != std::string::npos)
    {
      material[mat] = 1000 * atof(attr(line, "value=\"").c_str()); // no unit: g/cm^3
      mat.clear();
    }
    else if (line.find("<volume ") != std::string::npos)
      vol = attr(line, "name=\"");
    else if (!vol.empty() && line.find("<materialref ") != std::string::npos)
      volume[vol].first = attr(line, "ref=\"");
    else if (!vol.empty() && line.find("<solidref ") != std::string::npos)
      volume[vol].second = attr(line, "ref=\"");
    else if (line.find("<volumeref ") != std::string::npos)
      placed[attr(line, "ref=\"")]++;
  }
  double kg = 0;
  n = 0;
  for (auto &v : volume)
    if (v.second.first.find("Germanium") != std::string::npos && solid.count(v.second.second) && solid[v.second.second].size())
    {
      n += placed[v.first];
      kg += placed[v.first] * material[v.second.first] * solid[v.second.second].volume() * 1e-9;
    }
  return kg;
}

struct RT
{
  std::string file;
  Outline wall, argon, ofhcOuter, ofhcInner, ssOuter, ssInner;       // the tube (EFCu mother), the UGLAr inside it, the two shells
  double density[3] = {0, 0, 0};                                     // kg/m^3, the GDML's: steel, Cu, EFCu (the SS shell, the OFHC shell, the mother)
  double zBottom = 0, zTop = 0, seamOFHC = 0, seamSS = 0, zHead = 0; // zHead: the widest point, where the head meets the barrel
  double geMass = 0;                                                 // kg: the germanium detectors, the GDML's; the BI is per kg of them
  int nGe = 0;
  bool ok = false;
  const char *sectionAt(double zq) const { return zq < seamOFHC ? "EFCu" : (zq < seamSS ? "Cu" : "SS"); }
  const char *materialAt(double rq, double zq, double t = 1e-3) const // 1 um of slack: a float32 vertex lands off the surface it was drawn on
  {
    if (!wall.near(rq, zq, t))
      return "outside";
    if (ofhcOuter.near(rq, zq, t) && !ofhcInner.deep(rq, zq, t))
      return "OFHC";
    if (ssOuter.near(rq, zq, t) && !ssInner.deep(rq, zq, t))
      return "SS";
    if (argon.deep(rq, zq, t))
      return "LAr";
    // the mother's copper, by the seams: EFCu only below the Cu seam; between the seams it is the Cu section's (a 1.4 kg wedge
    // where the EFCu's tapered top meets the Cu shell's); above the steel seam the tube is open
    return zq > seamSS ? "outside" : (zq > seamOFHC ? "OFHC" : "EFCu");
  }
  int physAt(double rq, double zq) const
  {
    const char *m = materialAt(rq, zq);
    return !strcmp(m, "OFHC") ? 1 : (!strcmp(m, "SS") ? 2 : (!strcmp(m, "EFCu") ? 0 : -1));
  }                                                              // 0 mother, 1 OFHC, 2 SS
  // [mm^3]: the tube's wall between two heights: below the steel seam the solid minus the argon in it, above it the steel
  // shell alone: the tube is open at the top, its EFCu ends at the steel seam
  double tubeVolume(double z0, double z1, double dz = 0.5) const
  {
    z0 = std::max(std::min(z0, z1), zBottom);
    z1 = std::min(std::max(z0, z1), zTop);
    double v = 0;
    for (double z = z0; z < z1; z += dz)
    {
      double s = std::min(dz, z1 - z), zm = z + 0.5 * s;
      v += (zm < seamSS ? wall.areaAt(zm) - argon.areaAt(zm) : ssOuter.areaAt(zm) - ssInner.areaAt(zm)) * s;
    }
    return v;
  }
  double physVolume(int p) const // [mm^3]: each section is the tube's wall between its seams: EFCu (0) below the Cu seam, Cu (1), SS (2)
  {
    return p == 1 ? tubeVolume(seamOFHC, seamSS) : (p == 2 ? tubeVolume(seamSS, zTop) : tubeVolume(zBottom, seamOFHC));
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
  rt.geMass = readGermanium(rt.file, rt.nGe);
  if (rt.geMass <= 0)
  {
    printf("ERROR: no germanium detectors in %s\n", rt.file.c_str());
    return rt;
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
  printf("TUBE  read from %s: z %.1f .. %.1f mm (%.3f m), %d outline points, seams EFCu | %.1f | Cu | %.1f | SS\n", s.file.c_str(), s.zBottom, s.zTop,
         (s.zTop - s.zBottom) / 1000, n, s.seamOFHC, s.seamSS);
  printf("    clean surface: on-axis points %d (expect 2), duplicates %d, spikes %d, self-intersections %d (expect 0)\n", nAxis, nDup, nSpike, nCross);
  printf("    endcap: radius %.4f at z %.1f, barrel %.4f: step %+.4f mm%s\n", w.rmax, s.zHead, rBarrel, step, fabs(step) <= tol ? "" : "  <- NOT flush");
  printf("    shells against the barrel: Cu %+.4f, SS %+.4f mm\n", s.ofhcOuter.rmax - rBarrel, s.ssOuter.rmax - rBarrel);
  printf("    wall thickness [mm]:");
  for (double zq : {s.zHead, 0.5 * (s.zHead + s.seamOFHC), s.seamOFHC - 10, zBarrel, 0.5 * (s.seamSS + s.zTop)})
    printf("  %s %.2f (z %.0f)", s.sectionAt(zq), w.radiusAt(zq) - s.argon.radiusAt(zq), zq);
  printf("\n");
  printf("    densities [kg/m^3], the GDML's: SS %.0f, Cu %.0f, EFCu %.0f\n", s.density[0], s.density[1], s.density[2]);
  printf("    germanium, the GDML's: %d detectors, %.1f kg: the BI is per kg of it\n", s.nGe, s.geMass);
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
  for (int i = 0; i + 1 < n; i++) // each edge of the outline, split at the seams: the barrel takes its shell's colour, the endcap EFCu's
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
      if (std::min(pr[0], pr[1]) < rBarrel - 1 && std::min(ph[0], ph[1]) > s.seamSS - s.zBottom)
        continue; // above the steel seam the tube is its shell alone, open at the top
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
  r.nsim = vtx->GetEntries();                                 // the decays
  const Long64_t n = (Long64_t)vtx->GetMaximum("evtid") + 1; // the event numbers, up to it: some may have no vertex, no decay
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
    if (depth[i] < 0) // an event number with no decay
      continue;
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
    printf("A/E CUT  only %.0f DEP events net, kept at the production's %.2f\n\n", net(dep), cfg::psd_low);
    return;
  }
  double thr = 0;
  while (thr > -50 && kept(dep, thr) < cfg::depKeep)
    thr -= 0.005;
  gPsdLow = thr;
  printf("A/E CUT  tuned on the Tl208 DEP (%.1f keV) in M1 hits: %zu in +-%.1f keV, %zu in the sidebands %.0f-%.0f keV off it, %.0f net\n",
         cfg::depE, dep.peak.size(), cfg::depHalf, dep.side.size(), cfg::sideLo, cfg::sideHi, net(dep));
  printf("    classifier > %.2f keeps %.0f%% of it (+- %.0f%% statistics); the production's > %.2f kept %.0f%%\n", gPsdLow, 100 * kept(dep, gPsdLow),
         100 * sqrt(cfg::depKeep * (1 - cfg::depKeep) / net(dep)), cfg::psd_low, 100 * kept(dep, cfg::psd_low));
  printf("    kept by it, net of the continuum:");
  for (auto &l : line)
    printf("  %s %.1f %.2f (%.0f net)", l.name, l.e, kept(l, gPsdLow), net(l));
  printf("  ROI %.2f (%zu)\n    measured on V00000A at HADES:     ", win.empty() ? 0.0 : above(win, gPsdLow) / win.size(), win.size());
  for (auto &l : line)
    printf("  %s %.1f %.2f          ", l.name, l.e, l.measured);
  printf("  ROI 0.27\n\n");
}

void background(const char *runs = "output/output_root/tl208*_hit.root=Tl208,output/output_root/bi214*_hit.root=Bi214", const char *gdml = "geom/l1000geom/l1000.gdml")
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
  const double binW = (cfg::eHi - cfg::eLo) / cfg::eBins;
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
  { return p < 0 || r.mcDensity[p] <= 0 ? 0.0 : rt.density[m] * cfg::activity(r.iso, m) * cfg::secPerYear / (r.mcDensity[p] * rt.geMass); };
  auto massOf = [&](double d0, double d1, int m)
  { return d1 > d0 ? rt.density[m] * rt.tubeVolume(rt.zTop - d1, rt.zTop - d0) * 1e-9 : 0.0; };
  printf("INPUTS  ROI: %.0f-%.0f keV minus 10 keV around 2039, 2103.5, 2118.5, 2204.1 keV (MAJORANA's background region): %.0f keV\n", cfg::winLo, cfg::winHi, cfg::winWidth);
  printf("        cuts: M1 (exactly one detector above %.0f keV), LAr veto (the SiPMs see < %.0f photoelectrons), A/E classifier > %.2f (%.0f%% DEP, A/E CUT)\n",
         cfg::m1_keV, cfg::lar_pe, gPsdLow, 100 * cfg::depKeep);
  printf("        activity: chain [uBq/kg] 232Th / 238U: SS %g / %g, Cu %g / %g, EFCu %g / %g; Tl208 is %.2f%% of 232Th, Bi214 %.0f%% of 238U\n\n",
         cfg::chain[0][0], cfg::chain[0][1], cfg::chain[1][0], cfg::chain[1][1], cfg::chain[2][0], cfg::chain[2][1], 100 * cfg::branch[0], 100 * cfg::branch[1]);

  for (auto &r : R)
  {
    //-------------------------------------------------------------------------------
    //  6. Per run: the decays, the hits, the normalisation
    long w[NC] = {0};
    for (int c = 0; c < NC; c++)
      for (auto k : r.win[c])
        w[c] += k;
    printf("DECAYS  %s  (%s, %lld decays%s)\n", r.iso.c_str(), r.file.c_str(), r.nsim, r.nFiles > 1 ? Form(", %ld files", r.nFiles) : "");
    printf("    where they decayed: EFCu %.1f%%, Cu %.1f%%, SS %.1f%%; in the argon %ld, outside the tube %ld%s\n", 100.0 * r.nMat[0] / r.nsim,
           100.0 * r.nMat[1] / r.nsim, 100.0 * r.nMat[2] / r.nsim, r.nMat[3], r.nMat[4], r.nMat[3] || r.nMat[4] ? "  <- CONFINEMENT BUG" : "");
    printf("    germanium hits, per cut level (each row adds its cut to the ones above), at all energies and inside the %.0f keV ROI:\n", cfg::winWidth);
    printf("    %-16s %14s %14s\n", "cut level", "all energies", "in ROI");
    for (int c = 0; c < NC; c++)
      printf("    %-16s %14ld %14ld\n", cutName[c], r.nHit[c], w[c]);
    printf("    LAr veto keeps 1 in %.0f M1 hits; read off the optical map alone, without the argon beyond its edge, 1 in %.0f\n",
           r.nHit[2] ? (double)r.nHit[1] / r.nHit[2] : 0.0, r.nLarMap ? (double)r.nHit[1] / r.nLarMap : 0.0);
    const int asBuilt[3] = {2, 1, 0}; // mother EFCu, OFHC shell Cu, SS shell steel
    printf("    the normalisation, per section: volume, mass, activity, simulated decays (and per m^3), weight = real decays a year per simulated one\n");
    printf("    %-7s %9s %8s %10s %10s %11s %8s\n", "section", "V [m^3]", "M [kg]", "A [Bq]", "simulated", "per m^3", "weight");
    for (int p = 0; p < 3; p++)
    {
      int m = asBuilt[p];
      double mass = physV[p] * rt.density[m], A = mass * cfg::activity(r.iso, m);
      printf("    %-7s %9.5f %8.1f %10.3e %10ld %11.0f %8.3g\n", cfg::mat[m], physV[p], mass, A, r.nPv[p], r.mcDensity[p], r.nPv[p] ? A * cfg::secPerYear / r.nPv[p] : 0.0);
    }
    double dAvg = 0.5 * (r.mcDensity[1] + r.mcDensity[2]);
    printf("    the runs sampled the EFCu %.3f x as densely as the shells (1 = evenly): each section's weight uses its own density\n", dAvg > 0 ? r.mcDensity[0] / dAvg : 0.0);

    //-------------------------------------------------------------------------------
    //  7. Per run: where the hits come from, by depth; and the share of the M1 window hits all cuts keep there
    printf("DEPTH  %s  by depth below the tube's top: decays, germanium hits per decay, ROI events after M1 and after all cuts,\n", r.iso.c_str());
    printf("       and the cut survival = ROI events after all cuts / after M1\n");
    printf("    depth [m]            decays   hits/decay      ROI: M1   all cuts   survival\n");
    for (int b = 0; b < cfg::NDB; b++)
      printf("    %5.2f..%-8.2f %10ld %12.2e %12ld %10ld %6s%s\n", depthTot * b / cfg::NDB / 1000, depthTot * (b + 1) / cfg::NDB / 1000, r.dN[b],
             r.dN[b] ? (double)r.dH[b] / r.dN[b] : 0.0, r.dW1[b], r.dWA[b], r.dW1[b] ? Form("%.4f", (double)r.dWA[b] / r.dW1[b]) : "-",
             b == cfg::NDB - 1 ? "   <- nearest the detectors" : "");
    printf("\n");
  }

  //-------------------------------------------------------------------------------
  //  8. Background index as built, per chain and section:
  printf("BI AS BUILT  [cts/(keV kg yr)], +- simulation statistics; none: no ROI event in all its simulated decays\n");
  printf("    %-6s %-6s %15s   %-22s %-22s %-22s %-22s\n", "chain", "sect", "ROI hits", cutName[0], cutName[1], cutName[2], cutName[3]);
  double gB[NC] = {0}, gV[NC] = {0};
  long gN[NC] = {0};               // window hits behind each total
  double mB[2] = {0}, mV[2] = {0}; // the same after M1 + LAr and all cuts, the LAr veto off the map alone (steel and empty sections add nothing)
  long mN[2] = {0};
  auto total = [](double b, double v, long n, bool goal) -> std::string
  {
    if (!n)
      return "no ROI hit";
    return goal ? Form("%.2e +- %.1e (%.2g x goal)", b, sqrt(v), b / cfg::bgGoal) : Form("%.2e +- %.1e", b, sqrt(v));
  };
  const double edgeBuilt[4] = {0, L1built, L2built, depthTot};
  for (auto &r : R)
  {
    double tB[NC] = {0}, tV[NC] = {0};
    long tN[NC] = {0};
    for (int s = 0; s < 3; s++) // the slabs as built: steel, Cu, EFCu
    {
      double B[NC] = {0}, V[NC] = {0};
      long n[NC] = {0};
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
        cell[c] = n[c] ? Form("%.2e +- %.1e", B[c], sqrt(V[c])) : std::string("none");
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
  printf("    ROI hits: per cut level, as in the columns. sections with no ROI hit are left out of the totals\n");

  //-------------------------------------------------------------------------------
  //  9. Designs against the tube as built: steel to L1, Cu to L2, EFCu below, for less EFCu (L2 deeper) and more steel (L1 deeper)
  // a design differs from the tube as built only where a slab changed material, so its BI after all cuts is the as-built's plus
  // that change, dB. Too few events survive all cuts away from the detectors to count dB on them, so it is counted on the
  // changed slabs' ROI events after M1, 100-400 times more, times the cut survival (the share of them all cuts keep, see DEPTH),
  // measured where the change is: near the detectors, below farEdge, per chain; above it, both chains together.
  // a design passes when its BI stays within 10% of the as-built's, the BI's own statistical uncertainty
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
  const double sFar = fM ? (double)fA / fM : 0;
  std::vector<std::vector<size_t>> roiM1(R.size()); // per run, the events with an ROI hit after M1: the only ones a change counts
  for (size_t k = 0; k < R.size(); k++)
    for (size_t i = 0; i < R[k].depth.size(); i++)
      if (R[k].win[1][i])
        roiM1[k].push_back(i);
  auto change = [&](double L1, double L2)
  {
    double dB = 0;
    for (size_t k = 0; k < R.size(); k++)
    {
      auto &r = R[k];
      for (size_t i : roiM1[k])
      {
        const int mb = slabOf(r.depth[i], L1built, L2built), md = slabOf(r.depth[i], L1, L2);
        if (mb != md)
          dB += (r.depth[i] < farEdge ? sFar : r.S) * r.win[1][i] * (weight(r, md, r.pv[i]) - weight(r, mb, r.pv[i])) / cfg::winWidth;
      }
    }
    return dB;
  };
  auto deepestSteel = [&](double L2) { // 1 cm at a time below the steel as built, always some Cu
    double l1 = -1;
    for (double L1 = L1built; L1 <= L2 - 10; L1 += 10)
    {
      if (change(L1, L2) > maxRise)
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
  printf("\nDESIGNS  against the tube as built (BI %.2e after all cuts): SS down to L1, Cu down to L2, EFCu below. a design passes when its BI\n", biBuilt);
  printf("    stays within %.0f%% of the as-built's, the BI's own statistical uncertainty. L1 from %.2f m down, 1 cm at a time\n", 100 * cfg::simTol, L1built / 1000);
  printf("    a design's BI change is counted on ROI events after M1 x the cut survival (see DEPTH): below %.2f m", farEdge / 1000);
  for (auto &r : R)
    printf(" %s %.2e", r.iso.c_str(), r.S);
  printf("; above it, both chains, %ld of %ld = %.2e\n", fA, fM, sFar);
  printf("    %-6s %6s | %-10s %-34s %9s\n", "L2 [m]", "EFCu", "deepest L1", "SS : Cu : EFCu", "BI change");
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
      printf("    %-6.2f %3.0f kg | %-10s %s\n", L2 / 1000, massOf(L2, depthTot, 2), "-", "the BI rises more than 10% at any L1");
      continue;
    }
    printf("    %-6.2f %3.0f kg | %8.2f m %-34s %+8.0f%%%s\n", L2 / 1000, massOf(L2, depthTot, 2), l1 / 1000, split(l1, L2).c_str(),
           100 * change(l1, L2) / biBuilt, L2 == L2built ? "   <- EFCu as built" : "");
    if (L2 > best[1] || (L2 == best[1] && l1 > best[0]))
      best[0] = l1, best[1] = L2;
  }
  if (best[0] > 0)
    printf("    least EFCu, then most SS: L2 %.2f m, L1 %.2f m: %s, BI %+.0f%%\n", best[1] / 1000, best[0] / 1000, split(best[0], best[1]).c_str(),
           100 * change(best[0], best[1]) / biBuilt);
  printf("    as built: %s\n", split(L1built, L2built).c_str());

  //-------------------------------------------------------------------------------
  //  10. The tube's spectrum as built, per cut level, for the slide below
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
  gStyle->SetOptStat(0);
  const int col[NC] = {kAzure + 2, kOrange + 7, kRed + 1, kViolet + 1};
  for (int k = 0; k < NC; k++)
    sp[k]->SetLineColor(col[k]);
  sp[0]->SetMinimum(0.1 * std::max(1e-12, sp[ALL]->GetMinimum(0) > 0 ? sp[ALL]->GetMinimum(0) : sp[2]->GetMinimum(0)));

  //-------------------------------------------------------------------------------
  //  11. Two slide figures, one message each: the spectrum with the ROI and the BI marked; the design map near the tube as built
  auto sci = [](double x) { // "9.8 #times 10^{-9}"
    const int e = (int)floor(log10(x));
    return std::string(Form("%.1f #times 10^{%d}", x / pow(10, e), e));
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
    sp[0]->SetTitle("What the as-built tube leaves in the detectors;energy [keV];cts / (keV kg yr)");
    sp[0]->SetMaximum(5 * sp[0]->GetMaximum());
    sp[0]->Draw("HIST");
    gPad->Update();
    const double y0 = pow(10, gPad->GetUymin()), y1 = pow(10, gPad->GetUymax());
    std::vector<double> edges = {cfg::winLo}; // the ROI: 1950-2350 keV, less +-5 keV around each line
    for (double g : cfg::winGap)
      edges.insert(edges.end(), {g - cfg::winGapHalf, g + cfg::winGapHalf});
    edges.push_back(cfg::winHi);
    std::sort(edges.begin(), edges.end());
    const int roiCol = TColor::GetColor("#2e7d32");
    for (size_t i = 0; i + 1 < edges.size(); i += 2)
    {
      auto b = new TBox(edges[i], y0, edges[i + 1], y1);
      b->SetFillColorAlpha(roiCol, 0.15);
      b->SetLineWidth(0);
      b->Draw();
    }
    auto q = new TBox(cfg::Qbb - cfg::roiHalf, y0, cfg::Qbb + cfg::roiHalf, y1);
    q->SetFillColorAlpha(roiCol, 0.9);
    q->SetLineWidth(0);
    q->Draw();
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
    auto goal = new TLine(cfg::eLo, cfg::bgGoal, cfg::eHi, cfg::bgGoal);
    goal->SetLineStyle(2);
    goal->SetLineWidth(2);
    goal->SetLineColor(kGray + 2);
    goal->Draw();
    label(cfg::eHi - 20, 1.6 * cfg::bgGoal, "LEGEND-1000 goal: 10^{-5}", kGray + 2, 0.034, 31);
    auto bi = new TLine(cfg::winLo, biBuilt, cfg::winHi, biBuilt);
    bi->SetLineWidth(4);
    bi->SetLineColor(kBlack);
    bi->Draw();
    label(cfg::winHi + 15, 2.5 * biBuilt, "#splitline{BI = " + sci(biBuilt) + "}{purple, averaged over the ROI}", kBlack, 0.032, 11);
    label(cfg::winLo + 5, pow(10, gPad->GetUymax() - 0.06 * (gPad->GetUymax() - gPad->GetUymin())), "ROI, 360 keV", roiCol, 0.036);
    label(cfg::Qbb + cfg::roiHalf + 4, pow(10, gPad->GetUymax() - 0.13 * (gPad->GetUymax() - gPad->GetUymin())), "Q_{#beta#beta} #pm 2#sigma", roiCol, 0.036);
    label(2614.5 - 12, 0.6 * sp[0]->GetBinContent(sp[0]->FindBin(2614.5)), "Tl-208, 2614 keV", kGray + 3, 0.032, 32);
    cs->SaveAs(outPath("slide_spectrum.png").c_str());
  }
  { // the design map, near the tube as built
    double x0 = 1.6, x1 = 3.2, y0 = 3.9, y1 = 4.8; // m, widened to hold the chosen design
    if (best[0] > 0)
      x1 = std::max(x1, best[0] / 1000 + 0.4), y1 = std::max(y1, best[1] / 1000 + 0.2);
    const int NX = 48, NY = 36;
    auto hz = new TH2D("hZoom", "Cheaper tube designs: BI over the tube's as built;steel (SS) ends at depth L1 [m];EFCu starts at depth L2 [m]", NX, x0, x1, NY, y0, y1);
    for (int i = 1; i <= NX; i++)
      for (int j = 1; j <= NY; j++)
        hz->SetBinContent(i, j, (biBuilt + change(hz->GetXaxis()->GetBinCenter(i) * 1000, hz->GetYaxis()->GetBinCenter(j) * 1000)) / biBuilt);
    auto cd = new TCanvas("c_slide_designs", "", 1100, 800);
    cd->SetMargin(0.10, 0.15, 0.11, 0.08);
    gStyle->SetPalette(kBird);
    hz->SetMinimum(0.9);
    hz->SetMaximum(1.3);
    hz->SetContour(40);
    hz->GetZaxis()->SetTitle("BI / BI as built");
    hz->GetZaxis()->SetTitleOffset(1.3);
    hz->Draw("COLZ");
    double lev[1] = {1 + cfg::simTol};
    auto line = (TH2D *)hz->Clone("hZoomLine");
    line->SetContour(1, lev);
    line->SetLineColor(kBlack);
    line->SetLineWidth(4);
    line->Draw("CONT3 SAME");
    auto leg = new TLegend(0.12, 0.84, 0.62, 0.91);
    leg->SetFillColorAlpha(kWhite, 0.8);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.032);
    leg->AddEntry(line, "BI +10%: the designs left of / below it pass", "l");
    leg->Draw();
    auto star = new TMarker(L1built / 1000, L2built / 1000, 29);
    star->SetMarkerColor(kRed + 1);
    star->SetMarkerSize(3.5);
    star->Draw();
    label(L1built / 1000 + 0.04, L2built / 1000, "as built", kRed + 1, 0.036);
    if (best[0] > 0)
    {
      auto m = new TMarker(best[0] / 1000, best[1] / 1000, 21);
      m->SetMarkerColor(kBlack);
      m->SetMarkerSize(2.4);
      m->Draw();
      const bool right = best[0] / 1000 > 0.5 * (x0 + x1); // label on the side with room
      label(best[0] / 1000 + (right ? -0.04 : 0.04), best[1] / 1000 + 0.06, Form("EFCu from %.2f m, SS to %.2f m", best[1] / 1000, best[0] / 1000), kBlack, 0.034, right ? 32 : 12);
    }
    cd->SaveAs(outPath("slide_designs.png").c_str());
  }
  { // 12. Every design for LEGEND.app's rt (app/): a 2 cm grid through the seams as built, the BI over the as-built
    //     one, with the masses -> output/designs.json
    const double step = 20; // mm
    std::vector<double> g1, g2;
    for (double v = fmod(L1built, step); v <= depthTot; v += step)
      g1.push_back(v);
    for (double v = fmod(L2built, step); v <= depthTot; v += step)
      g2.push_back(v);
    std::ostringstream js;
    auto list = [&](const std::vector<double> &v, double scale, const char *fmt) {
      js << "[";
      for (size_t i = 0; i < v.size(); i++)
        js << (i ? "," : "") << Form(fmt, v[i] * scale);
      js << "]";
    };
    js << "{\"depthTot\":" << depthTot / 1000 << ",\"L1built\":" << L1built / 1000 << ",\"L2built\":" << L2built / 1000
       << ",\"biBuilt\":" << Form("%.4g", biBuilt) << ",\"tol\":" << cfg::simTol << ",\"l1\":";
    list(g1, 1e-3, "%.3f");
    js << ",\"l2\":";
    list(g2, 1e-3, "%.3f");
    std::vector<double> mSS, mCu1, mCu2, mEF;
    for (double v : g1)
      mSS.push_back(massOf(0, v, 0)), mCu1.push_back(massOf(0, v, 1));
    for (double v : g2)
      mCu2.push_back(massOf(0, v, 1)), mEF.push_back(massOf(v, depthTot, 2));
    js << ",\"mSS\":";
    list(mSS, 1, "%.1f");
    js << ",\"mCu1\":";
    list(mCu1, 1, "%.1f");
    js << ",\"mCu2\":";
    list(mCu2, 1, "%.1f");
    js << ",\"mEF\":";
    list(mEF, 1, "%.1f");
    std::ostringstream cen;
    for (size_t i = 0; i < g1.size(); i++)
    {
      cen << (i ? "," : "") << "[";
      for (size_t j = 0; j < g2.size(); j++)
      {
        const char *sep = j ? "," : "";
        if (g1[i] >= g2[j])
        {
          cen << sep << "null";
          continue;
        }
        cen << sep << Form("%.4g", (biBuilt + change(g1[i], g2[j])) / biBuilt);
      }
      cen << "]";
    }
    js << ",\"ratio\":[" << cen.str() << "]}";
    std::ofstream(outPath("designs.json")) << js.str();
  }
  printf("\nwrote output/tube.png, output/slide_spectrum.png, output/slide_designs.png and output/designs.json\n");
}
