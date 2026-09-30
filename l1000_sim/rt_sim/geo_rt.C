//  step 2, the meat: is the tube one clean surface, is the endcap flush, and how does it differ from the mint geometry

#include "rt_geom.h"

//-------------------------------------------------------------------------------
//  1. Segment crossing test:
static int rt_orient(double ar, double az, double br, double bz, double cr, double cz) // which way do we turn going a -> b -> c? 0 = straight, 1 = left, 2 = right
{
  double v = (bz - az) * (cr - br) - (br - ar) * (cz - bz); // cross-product sign trick, only the SIGN matters
  if (fabs(v) < 1e-9)
    return 0;
  return (v > 0) ? 1 : 2;
}

static bool rt_crosses(double r1, double z1, double r2, double z2,
                       double r3, double z3, double r4, double z4) // segments p1->p2 and p3->p4 cross iff each pair straddles the other's line
{
  int o1 = rt_orient(r1, z1, r2, z2, r3, z3);
  int o2 = rt_orient(r1, z1, r2, z2, r4, z4);
  int o3 = rt_orient(r3, z3, r4, z4, r1, z1);
  int o4 = rt_orient(r3, z3, r4, z4, r2, z2);
  return (o1 != o2) && (o3 != o4);
}

static void rt_line(const Outline &p) // one row of the solids table
{
  if (p.empty())
  {
    printf("      %-26s (absent)\n", p.name.c_str());
    return;
  }
  printf("      %-26s n=%-4d  r %9.4f .. %9.4f   z %9.2f .. %9.2f\n",
         p.name.c_str(), p.size(), p.rmin, p.rmax, p.zmin, p.zmax);
}

static TGraph *rt_graph(const Outline &p, int col, int style) // (r,z) silhouette, same axes as Bernhard's PDF
{
  auto g = new TGraph(p.size(), p.r.data(), p.z.data());
  g->SetLineColor(col);
  g->SetMarkerColor(col);
  g->SetLineStyle(style);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.3);
  return g;
}

//-------------------------------------------------------------------------------
//  2. Is the outline one clean surface?
void geo_rt(const char *gdml = "", const char *reference = "l1000.gdml", double tol_mm = 0.05)
{
  RT s = rtLoad(gdml); // one load: every check below reuses these six outlines
  if (!s.ok)
  {
    rtVerdict("geo_rt", false);
    return;
  }
  const Outline &w = s.wall;
  const int n = w.size();
  const double *r = w.r.data(); // .data() hands back a plain pointer to the vector's numbers
  const double *z = w.z.data();

  printf("file    : %s\n", s.file.c_str());
  printf("RT      : z %.1f .. %.1f mm   (length %.1f)   %d points\n", s.zBottom, s.zTop, s.zTop - s.zBottom, n);
  printf("seams   : EFCu | %.1f | OFHC | %.1f | SS\n\n", s.seamOFHC, s.seamSS);
  printf("[a] solids that define the RT wall:\n");
  for (const Outline *p : {&s.wall, &s.argon, &s.ofhcOuter, &s.ofhcInner, &s.ssOuter, &s.ssInner})
    rt_line(*p);

  int nAxis = 0; // expect exactly 2: the bottom tip and the top. a third means the outline jumped back to r=0 mid-way
  printf("\n[b] on-axis points (r=0), expect 2:\n");
  for (int i = 0; i < n; i++)
    if (r[i] == 0.0)
    {
      nAxis++;
      printf("      idx %-4d z = %.1f\n", i, z[i]);
    }
  printf("      -> found %d\n", nAxis);

  int nDup = 0; // the same (r,z) twice stalls the outline there. this is the generator artifact that caused the fatal "R/Z segments cross"
  printf("[c] duplicate points, expect 0:\n");
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      if (r[i] == r[j] && z[i] == z[j])
      {
        nDup++;
        printf("      idx %d == idx %d  (r=%.4f z=%.4f)\n", j, i, r[i], z[i]);
      }
  printf("      -> found %d\n", nDup);

  int nSpike = 0; // b->c pointing straight back along a->b, i.e. the outline folds onto itself like a needle
  printf("[d] 180-degree reversals (spikes), expect 0:\n");
  for (int i = 1; i < n - 1; i++)
  {
    double d1r = r[i] - r[i - 1], d1z = z[i] - z[i - 1];
    double d2r = r[i + 1] - r[i], d2z = z[i + 1] - z[i];
    double m1 = sqrt(d1r * d1r + d1z * d1z), m2 = sqrt(d2r * d2r + d2z * d2z);
    if (m1 <= 0 || m2 <= 0)
      continue;
    if ((d1r * d2r + d1z * d2z) / (m1 * m2) < -0.99) // normalised dot product, -1 is a full 180 reversal
    {
      nSpike++;
      printf("      idx %d  (r=%.4f z=%.4f)\n", i, r[i], z[i]);
    }
  }
  printf("      -> found %d\n", nSpike);

  int nCross = 0; // close the loop and test every pair of NON-NEIGHBOURING segments, neighbours legitimately share an end
  printf("[e] self-intersections, expect 0:\n");
  for (int i = 0; i < n; i++)
  {
    int i2 = (i + 1) % n; // % n wraps the last point back to the first
    for (int j = i + 2; j < n; j++)
    {
      int j2 = (j + 1) % n;
      if (i == 0 && j == n - 1)
        continue; // these two touch at point 0
      if (rt_crosses(r[i], z[i], r[i2], z[i2], r[j], z[j], r[j2], z[j2]))
      {
        nCross++;
        if (nCross <= 5)
          printf("      segment %d crosses segment %d\n", i, j);
      }
    }
  }
  printf("      -> found %d\n", nCross);
  bool ok_surface = (nAxis == 2) && (nDup == 0) && (nSpike == 0) && (nCross == 0);

  //-------------------------------------------------------------------------------
  //  3. Is the endcap flush with the wall?
  double z_barrel = 0.5 * (s.seamOFHC + s.seamSS); // read well inside OFHC, where the tube is certainly cylindrical
  double r_barrel = w.radiusAt(z_barrel);
  double r_endcap = w.rmax; // the endcap radius is simply the global maximum
  double step = r_endcap - r_barrel;
  double z_at_max = w.z[0], z_flush = 0;
  for (int i = 0; i < n; i++)
    if (r[i] == w.rmax)
    {
      z_at_max = z[i];
      z_flush = z[i + 1 < n ? i + 1 : i]; // where the profile first settles back onto the barrel
      break;
    }
  printf("\n[f] barrel vs endcap radius:\n");
  printf("      barrel  r = %9.4f mm   (measured at z = %.1f)\n", r_barrel, z_barrel);
  printf("      endcap  r = %9.4f mm   (peak at    z = %.1f)\n", r_endcap, z_at_max);
  printf("      step      = %+9.4f mm   (tolerance %.3f, over z = %.1f .. %.1f)\n",
         step, tol_mm, z_at_max, z_flush);
  bool ok_step = fabs(step) <= tol_mm;
  if (!ok_step)
    printf("      -> the endcap is NOT flush with the barrel: the profile steps by\n"
           "         %.4f mm where the head meets the cylinder.\n",
           fabs(step));

  printf("\n[g] barrel radius vs section shells:\n"); // the barrel must agree with both shell outer bounds
  bool ok_shells = true;
  for (const Outline *p : {&s.ofhcOuter, &s.ssOuter})
  {
    if (p->empty())
      continue;
    printf("      %-20s rmax = %9.4f   (diff %+0.4f)\n", p->name.c_str(), p->rmax, p->rmax - r_barrel);
    ok_shells &= fabs(p->rmax - r_barrel) <= tol_mm;
  }

  printf("\n[h] wall thickness = r(reentrancetube) - r(undergroundlar):\n"); // the shell is bounded outside by the RT and inside by the argon for the whole length, so one subtraction gives the thickness everywhere
  printf("      %10s  %6s  %10s %10s %10s\n", "z [mm]", "sect", "r_out", "r_in", "t [mm]");
  double zs[] = {z_at_max, 0.5 * (z_at_max + s.seamOFHC), s.seamOFHC - 10.0, z_barrel,
                 s.seamSS - 10.0, 0.5 * (s.seamSS + s.zTop), s.zTop - 10.0};
  for (double zq : zs) // range-for: walks every value in zs without an index
  {
    double ro = w.radiusAt(zq), ri = s.argon.radiusAt(zq);
    if (ro < 0 || ri < 0)
      continue;
    printf("      %10.1f  %6s  %10.4f %10.4f %10.4f\n", zq, s.sectionAt(zq), ro, ri, ro - ri);
  }

  //-------------------------------------------------------------------------------
  //  4. Diff against the mint geometry:
  RT s0;
  if (reference && *reference)
  {
    s0 = rtLoad(reference);
    if (!s0.ok)
      printf("\n[i] reference '%s' not found - skipping diff\n", reference);
    else
    {
      printf("\n[i] KS vs mint %s:\n", s0.file.c_str());
      printf("      %-22s %12s %12s %10s\n", "", "mint", "KS", "delta");
      printf("      %-22s %12.2f %12.2f %+10.2f\n", "RT z bottom", s0.zBottom, s.zBottom, s.zBottom - s0.zBottom);
      printf("      %-22s %12.2f %12.2f %+10.2f\n", "RT z top", s0.zTop, s.zTop, s.zTop - s0.zTop);
      printf("      %-22s %12.2f %12.2f %+10.2f\n", "RT length",
             s0.zTop - s0.zBottom, s.zTop - s.zBottom, (s.zTop - s.zBottom) - (s0.zTop - s0.zBottom));
      printf("      %-22s %12.4f %12.4f %+10.4f\n", "endcap rmax", s0.wall.rmax, w.rmax, w.rmax - s0.wall.rmax);
      printf("      %-22s %12.4f %12.4f %+10.4f\n", "barrel r",
             s0.wall.radiusAt(0.5 * (s0.seamOFHC + s0.seamSS)), r_barrel,
             r_barrel - s0.wall.radiusAt(0.5 * (s0.seamOFHC + s0.seamSS)));
      printf("      %-22s %12.2f %12.2f %+10.2f\n", "EFCu|OFHC seam", s0.seamOFHC, s.seamOFHC, s.seamOFHC - s0.seamOFHC);
      printf("      %-22s %12.2f %12.2f %+10.2f\n", "OFHC|SS seam", s0.seamSS, s.seamSS, s.seamSS - s0.seamSS);
    }
  }

  //-------------------------------------------------------------------------------
  //  5. Silhouette: whole tube, endcap, and the junction:
  auto c = new TCanvas("c_rt_verify", "RT verification", 1250, 780);
  c->Divide(3, 1);
  struct
  {
    const char *title;
    double z0, z1, r0, r1;
  } pad[] = {
      {"whole tube", s.zBottom - 150, s.zTop + 150, 0, w.rmax * 1.12},
      {"endcap bottom profile", s.zBottom - 30, z_at_max + 250, 0, w.rmax * 1.06},
      {"endcap / barrel junction", z_at_max - 120, z_at_max + 160, r_barrel - 2.5, w.rmax + 1.0}};

  for (int i = 0; i < 3; i++)
  {
    c->cd(i + 1)->SetGrid();
    gPad->SetLeftMargin(0.17); // room for the z axis title
    gPad->SetBottomMargin(0.12);
    auto fr = gPad->DrawFrame(pad[i].r0, pad[i].z0, pad[i].r1, pad[i].z1); // DrawFrame(xmin,ymin,xmax,ymax) fixes the axes before anything is drawn on them
    fr->SetTitle(Form("%s;r [mm];z [mm]", pad[i].title));
    fr->GetYaxis()->SetTitleOffset(1.7);
    fr->GetXaxis()->SetNdivisions(506); // 5 major ticks, otherwise the zoomed pad crams its labels
    if (s0.ok)
      rt_graph(s0.wall, kGray + 1, 2)->Draw("L SAME"); // mint, dashed grey
    rt_graph(w, kAzure + 2, 1)->Draw(i == 0 ? "L SAME" : "LP SAME");
    if (i == 0)
      for (double zs2 : {s.seamOFHC, s.seamSS}) // mark the section seams
      {
        auto l = new TLine(pad[i].r0, zs2, pad[i].r1, zs2);
        l->SetLineStyle(3);
        l->SetLineColor(kGray + 2);
        l->Draw();
      }
    if (i == 2) // annotate the step that makes this check fail
    {
      for (double rr : {r_barrel, r_endcap})
      {
        auto l = new TLine(rr, pad[i].z0, rr, pad[i].z1);
        l->SetLineStyle(2);
        l->SetLineColor(rr == r_endcap ? kRed + 1 : kGreen + 2);
        l->Draw();
      }
      auto t = new TLatex(r_barrel - 2.2, z_at_max + 45, Form("step = %+.4f mm", step));
      t->SetTextSize(0.04);
      t->SetTextColor(kRed + 1);
      t->Draw();
    }
  }
  if (s0.ok)
  {
    c->cd(1);
    auto leg = new TLegend(0.42, 0.76, 0.90, 0.88);
    leg->SetTextSize(0.035);
    leg->AddEntry(rt_graph(w, kAzure + 2, 1), "KS", "l");
    leg->AddEntry(rt_graph(s0.wall, kGray + 1, 2), "mint l1000", "l");
    leg->Draw();
  }
  c->SaveAs(rtOut("geo_rt.png").c_str());
  rtVerdict("geo_rt", ok_surface && ok_step && ok_shells);
}
