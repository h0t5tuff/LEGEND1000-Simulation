//  are the RT vertices spread over the WHOLE tube, and in the right material?
//    root -l -b -q 'ana/src_vertices.C("output/tl208.root")'

#include "../geom/rt.h"

void src_vertices(const char *fn, const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  if (!rt.ok)
    return;

//-------------------------------------------------------------------------------
//  1. Read the vertices:
  TFile in(fn);
  TTree *t = (TTree *)in.Get("stp/vtx");
  const Long64_t n = t ? t->GetEntries() : 0;
  if (n == 0)
  {
    printf("no vertices in stp/vtx of %s\n", fn);
    return;
  }
  TTreeFormula fx("fx", "xloc_in_m", t), fy("fy", "yloc_in_m", t), fz("fz", "zloc_in_m", t); // works for single- and double-precision files
  const int NB = 16;
  const char *mat[5] = {"EFCu", "OFHC", "SS", "LAr", "outside"};
  Long64_t h[NB] = {0}, nMat[5] = {0};
  double zmin = 1e9, zmax = -1e9, rmin = 1e9, rmax = -1e9, len = rt.zTop - rt.zBottom;
  for (Long64_t i = 0; i < n; i++)
  {
    t->GetEntry(i);
    double x = fx.EvalInstance() * 1000, y = fy.EvalInstance() * 1000, z = fz.EvalInstance() * 1000, r = sqrt(x * x + y * y); // m -> mm, the units of the GDML
    zmin = std::min(zmin, z); zmax = std::max(zmax, z);
    rmin = std::min(rmin, r); rmax = std::max(rmax, r);
    const char *m = rt.materialAt(r, z);
    int k = 0;
    while (k < 4 && strcmp(m, mat[k])) k++; // strcmp is 0 on a match; no match at all leaves k on "outside"
    nMat[k]++;
    h[std::min(NB - 1, std::max(0, (int)((z - rt.zBottom) / len * NB)))]++; // clamp: an edge vertex cannot fall outside the array
  }
  printf("file     : %s\ngeometry : %s\n", fn, rt.file.c_str());
  printf("vertices : %lld    z %.1f .. %.1f mm    r %.1f .. %.1f mm\n", n, zmin, zmax, rmin, rmax);
  printf("RT extent: %.1f .. %.1f mm    seams %.1f / %.1f\n\n", rt.zBottom, rt.zTop, rt.seamOFHC, rt.seamSS);

//-------------------------------------------------------------------------------
//  2. Material split:
  printf("[2] material (exact point-in-solid, not a z-cut):\n"); // a z-cut is wrong at the top: the EFCu lid sits ABOVE the SS seam
  for (int k = 0; k < 3; k++)
    printf("      %-4s %8lld (%5.1f%%)\n", mat[k], nMat[k], 100.0 * nMat[k] / n);
  if (nMat[3] || nMat[4])
    printf("      LAr  %8lld    outside %lld    <-- CONFINEMENT BUG\n", nMat[3], nMat[4]); // the source is not confined to the wall

//-------------------------------------------------------------------------------
//  3. Coverage along z:
  printf("\n[3] coverage along z:\n"); // an empty bin means a whole section never got sampled
  Long64_t mx = *std::max_element(h, h + NB);
  int empty = 0;
  for (int b = 0; b < NB; b++)
  {
    if (h[b] == 0) empty++;
    printf("   %8.0f..%8.0f mm | %-40s %lld%s\n", rt.zBottom + len * b / NB, rt.zBottom + len * (b + 1) / NB,
           TString('#', (int)(40.0 * h[b] / mx)).Data(), h[b], h[b] ? "" : "   <-- EMPTY"); // TString('#', k) is a bar of k hashes
  }
  if (empty)
    printf("\n%d empty bin(s) - the sampler never reached part of the tube.\n", empty);
  rtVerdict("src_vertices", nMat[3] == 0 && nMat[4] == 0 && empty == 0);
}
