//  Analyse bkgrnds: are the RT vertices spread over the WHOLE tube, and in the right material?

#include "rt_geom.h"

//-------------------------------------------------------------------------------
//  1. Read the vertices:
void src_vertices(const char *fn, const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  if (!rt.ok)
    return;
  TFile in(fn);
  TTree *t = (TTree *)in.Get("stp/vtx");
  if (!t)
  {
    printf("no stp/vtx tree in %s\n", fn);
    return;
  }
  const Long64_t n = t->GetEntries();
  if (n == 0)
  {
    printf("stp/vtx is empty\n");
    return;
  }
  TTreeFormula fx("fx", "xloc_in_m", t), fy("fy", "yloc_in_m", t), fz("fz", "zloc_in_m", t); // read type-agnostically: single- or double-precision files both work

  const int NB = 16;
  Long64_t h[NB] = {0};
  Long64_t nEF = 0, nOF = 0, nSS = 0, nLAr = 0, nOut = 0;
  double zmin = 1e9, zmax = -1e9, rmin = 1e9, rmax = -1e9;
  for (Long64_t i = 0; i < n; i++)
  {
    t->GetEntry(i);
    double x = fx.EvalInstance() * 1000.0; // m -> mm, the units of the GDML
    double y = fy.EvalInstance() * 1000.0;
    double z = fz.EvalInstance() * 1000.0;
    double r = sqrt(x * x + y * y);
    zmin = std::min(zmin, z);
    zmax = std::max(zmax, z);
    rmin = std::min(rmin, r);
    rmax = std::max(rmax, r);

    const char *m = rt.materialAt(r, z);
    if (!strcmp(m, "EFCu")) // strcmp returns 0 when the two strings match, so ! means "is equal"
      nEF++;
    else if (!strcmp(m, "OFHC"))
      nOF++;
    else if (!strcmp(m, "SS"))
      nSS++;
    else if (!strcmp(m, "LAr"))
      nLAr++;
    else
      nOut++;
    int b = (int)((z - rt.zBottom) / (rt.zTop - rt.zBottom) * NB);
    h[std::min(NB - 1, std::max(0, b))]++; // clamp so an edge vertex cannot fall outside the array
  }
  printf("file     : %s\n", fn);
  printf("geometry : %s\n", rt.file.c_str());
  printf("vertices : %lld    z %.1f .. %.1f mm    r %.1f .. %.1f mm\n", n, zmin, zmax, rmin, rmax);
  printf("RT extent: %.1f .. %.1f mm    seams %.1f / %.1f\n\n", rt.zBottom, rt.zTop, rt.seamOFHC, rt.seamSS);

  //-------------------------------------------------------------------------------
  //  2. Material split:
  printf("[a] material (exact point-in-solid, not a z-cut):\n"); // a z-cut is wrong at the top: the EFCu lid closing the tube sits ABOVE the SS seam
  printf("      EFCu %8lld (%5.1f%%)\n", nEF, 100.0 * nEF / n);
  printf("      OFHC %8lld (%5.1f%%)\n", nOF, 100.0 * nOF / n);
  printf("      SS   %8lld (%5.1f%%)\n", nSS, 100.0 * nSS / n);
  if (nLAr || nOut)
    printf("      LAr  %8lld    outside %lld    <-- CONFINEMENT BUG\n", nLAr, nOut); // vertices in argon mean the source is not confined to the wall

  //-------------------------------------------------------------------------------
  //  3. Coverage along z:
  printf("\n[b] coverage along z:\n"); // an empty bin is the old failure mode: only one of the three sections got sampled
  Long64_t mx = 1;
  for (int b = 0; b < NB; b++)
    mx = std::max(mx, h[b]);
  int empty = 0;
  for (int b = 0; b < NB; b++)
  {
    double z0 = rt.zBottom + (rt.zTop - rt.zBottom) * b / NB;
    double z1 = rt.zBottom + (rt.zTop - rt.zBottom) * (b + 1) / NB;
    if (h[b] == 0)
      empty++;
    printf("   %8.0f..%8.0f mm | %-40s %lld%s\n", z0, z1,
           TString('#', (int)(40.0 * h[b] / mx)).Data(), h[b], h[b] ? "" : "   <-- EMPTY"); // TString('#', k) is a bar of k hashes
  }
  if (empty)
    printf("\n%d empty bin(s) - the sampler never reached part of the tube.\n", empty);
  rtVerdict("src_vertices", (nLAr == 0) && (nOut == 0) && (empty == 0));
}
