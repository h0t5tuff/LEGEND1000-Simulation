//  is the tube one clean surface, is the endcap flush, and how does it differ from the mint geometry
//    root -l -b -q sim/tube.C

#include "rt.h"

static int turn(double ar, double az, double br, double bz, double cr, double cz) // which way a -> b -> c turns: 0 straight, 1 left, 2 right
{
  double v = (bz - az) * (cr - br) - (br - ar) * (cz - bz); // cross product, only the sign matters
  return fabs(v) < 1e-9 ? 0 : (v > 0 ? 1 : 2);
}

static TGraph *silhouette(const Outline &p, int col, int style) // (r,z) outline of one solid
{
  auto g = new TGraph(p.size(), p.r.data(), p.z.data());
  g->SetLineColor(col);
  g->SetLineStyle(style);
  g->SetMarkerColor(col);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  return g;
}

void tube(const char *gdml = "", const char *reference = "l1000.gdml", double tol_mm = 0.05)
{
  RT s = rtLoad(gdml); // one load: every check below reuses these six outlines
  if (!s.ok)
    return rtVerdict("tube", false);
  const Outline &w = s.wall;
  const std::vector<double> &r = w.r, &z = w.z;
  const int n = w.size();
  printf("file : %s\nRT   : z %.1f .. %.1f mm, %d points, seams EFCu | %.1f | OFHC | %.1f | SS\n",
         s.file.c_str(), s.zBottom, s.zTop, n, s.seamOFHC, s.seamSS);

//-------------------------------------------------------------------------------
//  1. The solids that bound the wall:
  printf("\n[1] solids:\n");
  for (const Outline *p : {&s.wall, &s.argon, &s.ofhcOuter, &s.ofhcInner, &s.ssOuter, &s.ssInner})
    printf("      %-22s n=%-4d r %9.4f .. %9.4f   z %9.2f .. %9.2f\n", p->name.c_str(), p->size(), p->rmin, p->rmax, p->zmin, p->zmax);

//-------------------------------------------------------------------------------
//  2. One clean surface:
  int nAxis = 0, nDup = 0, nSpike = 0, nCross = 0;
  for (int i = 0; i < n; i++)
  {
    int i2 = (i + 1) % n; // % n wraps the last point back to the first
    if (r[i] == 0) nAxis++; // only the bottom tip and the top sit on the axis
    for (int j = 0; j < i; j++)
      if (r[i] == r[j] && z[i] == z[j]) nDup++; // the same point twice stalls the outline
    if (i > 0 && i < n - 1)
    {
      double ar = r[i] - r[i - 1], az = z[i] - z[i - 1], br = r[i + 1] - r[i], bz = z[i + 1] - z[i];
      double m = sqrt((ar * ar + az * az) * (br * br + bz * bz));
      if (m > 0 && (ar * br + az * bz) / m < -0.99) nSpike++; // the outline folds straight back on itself
    }
    for (int j = i + 2; j < n; j++) // every pair of segments that are not neighbours
    {
      int j2 = (j + 1) % n;
      if (i == 0 && j == n - 1) continue; // these two share point 0
      if (turn(r[i], z[i], r[i2], z[i2], r[j], z[j]) != turn(r[i], z[i], r[i2], z[i2], r[j2], z[j2]) &&
          turn(r[j], z[j], r[j2], z[j2], r[i], z[i]) != turn(r[j], z[j], r[j2], z[j2], r[i2], z[i2]))
        nCross++; // each segment straddles the other's line: they cross
    }
  }
  bool okSurface = nAxis == 2 && nDup == 0 && nSpike == 0 && nCross == 0;
  printf("\n[2] clean surface: on-axis points %d (expect 2), duplicates %d, spikes %d, self-intersections %d (expect 0)\n",
         nAxis, nDup, nSpike, nCross);

//-------------------------------------------------------------------------------
//  3. Endcap flush with the barrel:
  double zBarrel = 0.5 * (s.seamOFHC + s.seamSS), rBarrel = w.radiusAt(zBarrel); // mid-OFHC, where the tube is certainly a cylinder
  double zEndcap = s.zHead, step = w.rmax - rBarrel; // the endcap radius is the global maximum
  bool okStep = fabs(step) <= tol_mm;
  printf("\n[3] endcap flush: barrel r %.4f (z %.1f), endcap r %.4f (z %.1f), step %+.4f mm (tolerance %.3f)%s\n",
         rBarrel, zBarrel, w.rmax, zEndcap, step, tol_mm, okStep ? "" : "  <- NOT flush");

//-------------------------------------------------------------------------------
//  4. Barrel agrees with both shells:
  bool okShells = true;
  printf("\n[4] shells:\n");
  for (const Outline *p : {&s.ofhcOuter, &s.ssOuter})
  {
    printf("      %-22s rmax %9.4f (barrel %+.4f)\n", p->name.c_str(), p->rmax, p->rmax - rBarrel);
    okShells &= fabs(p->rmax - rBarrel) <= tol_mm;
  }

//-------------------------------------------------------------------------------
//  5. Wall thickness = r(reentrancetube) - r(undergroundlar):
  printf("\n[5] wall thickness:\n      %10s %6s %10s %10s %8s\n", "z [mm]", "sect", "r_out", "r_in", "t [mm]");
  for (double zq : {zEndcap, 0.5 * (zEndcap + s.seamOFHC), s.seamOFHC - 10, zBarrel, s.seamSS - 10, 0.5 * (s.seamSS + s.zTop), s.zTop - 10})
    printf("      %10.1f %6s %10.4f %10.4f %8.4f\n", zq, s.sectionAt(zq), w.radiusAt(zq), s.argon.radiusAt(zq), w.radiusAt(zq) - s.argon.radiusAt(zq));

//-------------------------------------------------------------------------------
//  6. Diff against the mint geometry:
  std::string ref = rtFindFile(reference);
  RT s0;
  if (!ref.empty()) s0 = rtLoad(ref.c_str());
  if (!s0.ok)
    printf("\n[6] reference '%s' not found - no diff\n", reference);
  else
  {
    double rBarrel0 = s0.wall.radiusAt(0.5 * (s0.seamOFHC + s0.seamSS));
    struct { const char *name; double mint, ks; } row[] = {
        {"RT z bottom", s0.zBottom, s.zBottom}, {"RT z top", s0.zTop, s.zTop}, {"RT length", s0.zTop - s0.zBottom, s.zTop - s.zBottom},
        {"endcap rmax", s0.wall.rmax, w.rmax}, {"barrel r", rBarrel0, rBarrel},
        {"EFCu|OFHC seam", s0.seamOFHC, s.seamOFHC}, {"OFHC|SS seam", s0.seamSS, s.seamSS}};
    printf("\n[6] KS vs mint %s:\n      %-16s %11s %11s %10s\n", s0.file.c_str(), "", "mint", "KS", "delta");
    for (auto &q : row)
      printf("      %-16s %11.4f %11.4f %+10.4f\n", q.name, q.mint, q.ks, q.ks - q.mint);
  }

//-------------------------------------------------------------------------------
//  7. Draw: whole tube, endcap, junction
  auto c = new TCanvas("c_rt", "", 1250, 780);
  c->Divide(3, 1);
  struct { const char *title; double z0, z1, r0, r1; } pad[] = {
      {"whole tube", s.zBottom - 150, s.zTop + 150, 0, w.rmax * 1.12},
      {"endcap bottom profile", s.zBottom - 30, zEndcap + 250, 0, w.rmax * 1.06},
      {"endcap / barrel junction", zEndcap - 120, zEndcap + 160, rBarrel - 2.5, w.rmax + 1.0}};
  for (int i = 0; i < 3; i++)
  {
    c->cd(i + 1)->SetGrid();
    gPad->SetLeftMargin(0.17);
    auto fr = gPad->DrawFrame(pad[i].r0, pad[i].z0, pad[i].r1, pad[i].z1); // fixes the axes before anything is drawn
    fr->SetTitle(Form("%s;r [mm];z [mm]", pad[i].title));
    fr->GetYaxis()->SetTitleOffset(1.7);
    fr->GetXaxis()->SetNdivisions(506);
    if (s0.ok) silhouette(s0.wall, kGray + 1, 2)->Draw("L SAME"); // mint, dashed grey
    silhouette(w, kAzure + 2, 1)->Draw(i ? "LP SAME" : "L SAME");   // zoomed pads show the corner points too
    for (double v : {s.seamOFHC, s.seamSS, rBarrel, w.rmax}) // pad 1: the two seams. pad 3: barrel and endcap radius
    {
      bool seam = v == s.seamOFHC || v == s.seamSS;
      if ((seam && i != 0) || (!seam && i != 2)) continue;
      auto l = seam ? new TLine(pad[i].r0, v, pad[i].r1, v) : new TLine(v, pad[i].z0, v, pad[i].z1);
      l->SetLineStyle(seam ? 3 : 2);
      l->SetLineColor(seam ? kGray + 2 : (v == w.rmax ? kRed + 1 : kGreen + 2));
      l->Draw();
    }
    if (i == 2)
    {
      auto t = new TLatex(rBarrel - 2.2, zEndcap + 45, Form("step = %+.4f mm", step));
      t->SetTextColor(kRed + 1);
      t->SetTextSize(0.04);
      t->Draw();
    }
  }
  c->SaveAs(rtOut("tube.png").c_str());
  rtVerdict("tube", okSurface && okStep && okShells);
}
