//  read the reentrance tube's shape from the GDML in two steps.

#ifndef RT_H // include guard in case two macros have the same header
#define RT_H
#include <algorithm>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include "TROOT.h"
#include "TSystem.h"
#include "TTree.h"
#include "TTreeFormula.h"

//-------------------------------------------------------------------------------
//  1. Find RT Outline:
struct Outline // struct is just a bundle of variables under one name plus functions that work on that bundle. Those functions are called member functions
{
  std::string name;
  std::vector<double> r, z;                      // the corner points, in order
  double rmin = 0, rmax = 0, zmin = 0, zmax = 0; // bounding box
  int size() const { return (int)r.size(); }
  bool empty() const { return r.empty(); }
  // How wide is the shape at height zq?
  double radiusAt(double zq) const // const after a member function promise that the function only reads the struct, never modifies it.
  {
    double best = -1;
    for (int i = 0; i + 1 < size(); i++)
    {
      double z1 = z[i], z2 = z[i + 1];
      if ((zq < z1 && zq < z2) || (zq > z1 && zq > z2))
        continue; // edge misses zq

      double rr;
      if (fabs(z2 - z1) < 1e-12)
      {
        rr = std::max(r[i], r[i + 1]); // horizontal edge: take its wider end
      }
      else
      {
        double frac = (zq - z1) / (z2 - z1); // Linear interpolation: how far along the edge is zq, as a fraction?
        rr = r[i] + frac * (r[i + 1] - r[i]);
      }
      if (rr > best)
        best = rr;
    }
    return best; // Returns -1 if zq is above or below the whole shape.
  }
  // Is the point (rq, zq) inside this outline? using ray casting trick
  bool contains(double rq, double zq) const
  {
    bool inside = false;
    const int n = size();
    for (int i = 0, j = n - 1; i < n; j = i++)
    {
      bool straddles = (z[i] > zq) != (z[j] > zq); // Does this edge straddle our height? If both ends are on the same side of zq, the ray cannot cross it.
      if (!straddles)
        continue;

      double frac = (zq - z[i]) / (z[j] - z[i]); // Where does the edge sit at our height? If it is to our right, the ray (fired to the right) crosses it, so flip inside/outside.
      double rEdge = r[i] + frac * (r[j] - r[i]);
      if (rq < rEdge)
        inside = !inside;
    }
    return inside;
  }
  // tolerant versions: inside or within tol of it (near), at least tol inside (deep). the point is
  // nudged by tol in r AND z: on the lid's shallow cone a float32 z error of 0.2 um is 40 um in r
  bool near(double rq, double zq, double tol) const
  {
    return contains(rq, zq) || contains(rq + tol, zq) || contains(fabs(rq - tol), zq) ||
           contains(rq, zq + tol) || contains(rq, zq - tol);
  }
  bool deep(double rq, double zq, double tol) const
  {
    return contains(rq, zq) && contains(rq + tol, zq) && contains(fabs(rq - tol), zq) &&
           contains(rq, zq + tol) && contains(rq, zq - tol);
  }
  // Cross-section area at height zq [mm^2]. The line z = zq crosses the outline an even number
  // of times; sorted, the crossings pair up into rings [r_in, r_out]. radiusAt only sees the
  // outermost one, which is wrong wherever the solid is a ring rather than a disc: the argon at
  // the very top is closed by a cone, so there it is a ring and the copper lid sits inside it
  double areaAt(double zq) const
  {
    std::vector<double> x;
    for (int i = 0, j = size() - 1; i < size(); j = i++)
      if ((z[i] > zq) != (z[j] > zq))
        x.push_back(r[i] + (zq - z[i]) / (z[j] - z[i]) * (r[j] - r[i]));
    std::sort(x.begin(), x.end());
    double a = 0;
    for (size_t k = 0; k + 1 < x.size(); k += 2)
      a += x[k + 1] * x[k + 1] - x[k] * x[k]; // one ring, pi applied once below
    return M_PI * a;
  }
  // Volume of the whole solid [mm^3]: stack the cross-sections up its height
  double volume(double dz = 0.5) const
  {
    double v = 0;
    for (double zz = zmin; zz < zmax; zz += dz)
    {
      double step = std::min(dz, zmax - zz);
      v += areaAt(zz + 0.5 * step) * step;
    }
    return v;
  }
};

inline Outline rtRead(const std::string &gdmlFile, const char *solid) // Pull reentrancetube solid out of the GDML. tells the compiler it is fine for this function body to appear in several files at once. Header-only code needs it.
{
  Outline out;
  out.name = solid;

  std::ifstream in(gdmlFile.c_str());
  if (!in)
  {
    printf("  ERROR: cannot open %s\n", gdmlFile.c_str());
    return out;
  }
  const std::string wanted = std::string("name=\"") + solid + "\"";
  std::string line;
  bool collecting = false;
  while (std::getline(in, line))
  {
    if (!collecting)
    {
      // std::string::npos is what find() returns when it finds nothing.
      if (line.find("<genericPolycone") != std::string::npos &&
          line.find(wanted) != std::string::npos)
        collecting = true;
      continue;
    }
    if (line.find("</genericPolycone>") != std::string::npos)
      break;
    double rr, zz;
    if (sscanf(line.c_str(), " <rzpoint r=\"%lf\" z=\"%lf\"", &rr, &zz) == 2)
    {
      out.r.push_back(rr);
      out.z.push_back(zz);
    }
  }
  if (!out.empty())
  {
    out.rmin = *std::min_element(out.r.begin(), out.r.end());
    out.rmax = *std::max_element(out.r.begin(), out.r.end());
    out.zmin = *std::min_element(out.z.begin(), out.z.end());
    out.zmax = *std::max_element(out.z.begin(), out.z.end());
  }
  return out;
}

//-------------------------------------------------------------------------------
//  2. Configure RT design:
struct RT
{
  std::string file;
  Outline wall, argon, ofhcOuter, ofhcInner, ssOuter, ssInner;
  double zBottom = 0, zTop = 0;    // full extent of the tube
  double seamOFHC = 0, seamSS = 0; // where one section hands over to the next
  double zHead = 0;                // where the bottom head meets the barrel: the height of the tube's widest point
  bool ok = false;
  const char *sectionAt(double zq) const
  {
    if (zq < seamOFHC)
      return "EFCu";
    if (zq < seamSS)
      return "OFHC";
    return "SS";
  }
  // Volume of the tube wall between two heights [mm^3]: the tube's solid minus the argon inside
  // it, cross-section by cross-section. Uses areaAt, so the lid at the top is counted
  double wallVolume(double z0, double z1, double dz = 0.5) const
  {
    if (z1 < z0)
      std::swap(z0, z1);
    z0 = std::max(z0, zBottom);
    z1 = std::min(z1, zTop);
    double v = 0;
    for (double z = z0; z < z1; z += dz)
    {
      double step = std::min(dz, z1 - z);
      v += (wall.areaAt(z + 0.5 * step) - argon.areaAt(z + 0.5 * step)) * step;
    }
    return v;
  }
  // Volume of each physical volume the source samples [mm^3]: 0 the reentrancetube mother
  // (endcap, lower wall, lid), 1 the ofhc_cu shell, 2 the ss_316l shell. The shells are outer
  // minus inner bound; the mother is whatever wall is left once both shells are taken out
  double physVolume(int pv) const
  {
    double ofhc = ofhcOuter.volume() - ofhcInner.volume(), ss = ssOuter.volume() - ssInner.volume();
    if (pv == 1) return ofhc;
    if (pv == 2) return ss;
    return wallVolume(zBottom, zTop) - ofhc - ss;
  }
  // Which physical volume a point sits in: 0 mother, 1 OFHC, 2 SS, -1 not in the wall at all
  int physAt(double rq, double zq) const
  {
    const char *m = materialAt(rq, zq);
    if (!strcmp(m, "OFHC")) return 1;
    if (!strcmp(m, "SS")) return 2;
    if (!strcmp(m, "EFCu")) return 0;
    return -1;
  }

  const char *materialAt(double rq, double zq, double tol = 1e-3) const // 1 um of slack in r and z: a vertex written as float32 lands off the surface it was generated on, and without slack that reads as a confinement bug
  {
    if (!wall.near(rq, zq, tol))
      return "outside";
    if (ofhcOuter.near(rq, zq, tol) && !ofhcInner.deep(rq, zq, tol))
      return "OFHC";
    if (ssOuter.near(rq, zq, tol) && !ssInner.deep(rq, zq, tol))
      return "SS";

    if (argon.deep(rq, zq, tol))
      return "LAr";
    return "EFCu";
  }
};

inline std::string rtFindFile(const char *name) // run macros anywhere
{
  std::string path = name;
  for (int up = 0; up < 6; up++)
  {
    if (!gSystem->AccessPathName(path.c_str()))
      return path;
    path = "../" + path;
  }
  return "";
}

inline std::string rtResolve(const char *gdmlArg = "") // Which GDML file do we mean?
{
  const char *env = gSystem->Getenv("RT_GDML");
  std::string f;
  if (gdmlArg && *gdmlArg)
    f = rtFindFile(gdmlArg);
  if (f.empty() && env)
    f = rtFindFile(env);
  if (f.empty())
    f = rtFindFile("KSendcap_l1kGeometry.gdml");
  return f;
}

inline RT rtLoad(const char *gdmlArg = "") // Load every outline we need, and measure the seams.
{
  RT rt;
  rt.file = rtResolve(gdmlArg);
  if (rt.file.empty())
  {
    printf("ERROR: no GDML found\n");
    return rt;
  }
  rt.wall = rtRead(rt.file, "reentrancetube");
  rt.argon = rtRead(rt.file, "undergroundlar");
  rt.ofhcOuter = rtRead(rt.file, "ofhc_cu_outer_bound");
  rt.ofhcInner = rtRead(rt.file, "ofhc_cu_inner_bound");
  rt.ssOuter = rtRead(rt.file, "ss_316l_outer_bound");
  rt.ssInner = rtRead(rt.file, "ss_316l_inner_bound");
  if (rt.wall.empty() || rt.ofhcOuter.empty() || rt.ssOuter.empty())
  {
    printf("ERROR: %s does not look like an l1000 geometry\n", rt.file.c_str());
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

inline std::string rtOut(const std::string &name) // every generated file lands in output/, so the source folder stays clean
{
  gSystem->mkdir("output", kTRUE);                       // kTRUE: make parents too, and do not complain if it exists
  return "output/" + name.substr(name.find_last_of("/\\") + 1); // npos+1 == 0, so a bare name passes through unchanged
}

inline bool rtIsGermanium(const TString &n) // a Ge detector's tree is "V" + 4 digits, e.g. V0101
{
  if (n.Length() != 5 || n[0] != 'V')
    return false;
  for (int i = 1; i < 5; i++)
    if (!isdigit(n[i]))
      return false;
  return true;
}

// walk a tree row by row and hand the named columns to f as doubles, whatever type they were stored in:
//   rtScan(t, {"evtid", "edep_in_keV"}, [&](const double *v) { sum[(int)v[0]] += v[1]; });
template <class F> inline void rtScan(TTree *t, const std::vector<const char *> &cols, F f)
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
    f(v.data()); // .data() is a plain pointer to the row's values
  }
  for (auto c : col)
    delete c;
}

inline void rtVerdict(const char *tag, bool ok) // One PASS/FAIL line, printed the same way by every check.
{
  printf("\nRESULT [%s]: %s\n", tag, ok ? "PASS" : "FAIL");
  if (gROOT->IsBatch())
    gSystem->Exit(ok ? 0 : 1);
}
#endif // RT_H
