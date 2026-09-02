//  read the reentrance tube's shape from the GDML in two steps.

#ifndef RT_GEOM_H // include guard in case two macros have the same header
#define RT_GEOM_H
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include "TROOT.h"
#include "TSystem.h"

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
  bool contains(double rq, double zq, double tol = 0) const // tol grows (or with a negative value shrinks) the outline radially
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
      if (rq < rEdge + tol)
        inside = !inside;
    }
    return inside;
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
  bool ok = false;
  const char *sectionAt(double zq) const
  {
    if (zq < seamOFHC)
      return "EFCu";
    if (zq < seamSS)
      return "OFHC";
    return "SS";
  }
  const char *materialAt(double rq, double zq, double tol = 1e-3) const // 1 um of slack: a vertex written as float32 lands a few nm off the surface it was generated on, and without slack that reads as a confinement bug
  {
    if (!wall.contains(rq, zq, tol))
      return "outside";
    if (ofhcOuter.contains(rq, zq, tol) && !ofhcInner.contains(rq, zq, -tol))
      return "OFHC";
    if (ssOuter.contains(rq, zq, tol) && !ssInner.contains(rq, zq, -tol))
      return "SS";

    if (argon.contains(rq, zq, -tol))
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
  rt.ok = true;
  return rt;
}

inline void rtVerdict(const char *tag, bool ok) // One PASS/FAIL line, printed the same way by every check.
{
  printf("\nRESULT [%s]: %s\n", tag, ok ? "PASS" : "FAIL");
  if (gROOT->IsBatch())
    gSystem->Exit(ok ? 0 : 1);
}
#endif // RT_GEOM_H
