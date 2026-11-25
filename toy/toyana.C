//  The miniature tube study: weight, reweight, scan the seam. Reads remage's two runs of tut/toysim.mac (tut/toygeom.gdml).
//      root -l -b -q tut/toyana.C
//  1. Constants: the GDML's tube and crystal, the true materials and activities (as in ana/ana.C's cfg)
//  2. Read the hits: the energy each gamma leaves in the crystal, and the depth z where that gamma started
//  3. Weight: one simulated gamma stands for w real gammas a year, w = rho A V t / N
//  4. Reweight: a seam at L makes the wall below L EFCu and above it steel; the hits stay, only w changes
//  5. Scan L and pick the new design: the lowest seam whose rise in the analysis window, plus its statistical error,
//     stays within the line (10%)
//  6. Plots: B(E) for both designs, and B(E) for every seam L as a heat map
#include <unordered_map>

namespace cfg
{
  // the two runs: the gammas started in the EFCu section (z 0-2 m) and in the steel section (z 2-4 m)
  const char *inFile[2] = {"tut/toyoutput/toy_efcu.root", "tut/toyoutput/toy_ss.root"};
  // the tube (tut/toygeom.gdml): 4 m tall, inner radius 0.25 m, 1 cm wall; steel above the seam, EFCu below; proposed seam at 2 m
  const double H = 4.0, R = 0.25, T = 0.01, L0 = 2.0;
  const double Vsec = TMath::Pi() * ((R + T) * (R + T) - R * R) * H / 2; // one section's volume [m^3]
  // the crystal: 1 kg of enriched Ge (90% Ge-76), radius 30.61 mm, height 61.22 mm, density 5.5502 g/cm^3
  const double mGe = TMath::Pi() * 3.061 * 3.061 * 6.122 * 5.550157 / 1000; // [kg]
  // materials {EFCu, steel}: density [kg/m^3] (the GDML's) and Th-232 chain activity [Bq/kg] (EFCu: LEGEND-1000 CDR Table XVII;
  // steel: L1000_L.CS_rev2.pdf, Table 2), as in ana/ana.C
  const double rho[2] = {8960, 8000};
  const double th232[2] = {0.37e-6, 1000e-6};
  // the source: only the 2614.5 keV gamma of Tl-208. Tl-208 is 35.94% of the chain, and 99.755% of its decays give that gamma
  const double perChain = 0.3594 * 0.99755;
  const double year = 365.25 * 24 * 3600; // [s]
  const double tol = 0.10;                // the new design may raise the counts a year by at most this: the main study's
                                          // line (cfg::simTol in ana/ana.C), a choice, not a LEGEND requirement
  const double kSig = 1;                  // a seam passes only if its rise + kSig x its statistical error stays within tol
  const double dRes = 0.02;               // how precisely we want the seam [m]: simulate more until the undecided band is this narrow
  const double Qbb = 2039;                // Qbb of Ge-76 [keV]
  const double winHalf = 50;              // analysis window: Qbb +- 50 keV, where the background index B is counted [keV]
  const double roi = 2.5;                 // ROI: the width around Qbb where 0vbb would show [keV]
  const double eMax = 3000, eBin = 50;    // spectrum [keV]; bins are placed so the window is exactly two of them
  const double dL = 0.01, dLmap = 0.05;   // seam step of the scan and of the heat map [m]
}

void toyana()
{
  using namespace cfg;
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);

  //---------------------------------------------------------------------------
  //  2. Read the hits of both runs
  struct Hit
  {
    double E, z;
    int run; // 0: started in the EFCu section, 1: in the steel section
  };
  std::vector<Hit> hit; // one per gamma that reached the crystal
  double N[2] = {0, 0}; // gammas simulated per run
  for (int r = 0; r < 2; r++)
  {
    TFile f(inFile[r]);
    auto ge = (TTree *)f.Get("stp/ge");
    auto vtx = (TTree *)f.Get("stp/vtx");
    if (!ge || !vtx)
    {
      printf("ERROR: no stp/ge or stp/vtx in %s: run tut/toysim.mac first\n", inFile[r]);
      return;
    }
    N[r] = vtx->GetEntries();
    int id;
    float e, z;
    std::unordered_map<int, double> edep; // gamma -> energy in the crystal [keV]
    ge->SetBranchAddress("evtid", &id);
    ge->SetBranchAddress("edep_in_keV", &e);
    for (Long64_t i = 0; i < ge->GetEntries(); i++)
      ge->GetEntry(i), edep[id] += e;
    vtx->SetBranchAddress("evtid", &id);
    vtx->SetBranchAddress("zloc_in_m", &z);
    size_t n0 = hit.size();
    double eTop = 0;
    for (Long64_t i = 0; i < vtx->GetEntries(); i++)
    {
      vtx->GetEntry(i);
      auto it = edep.find(id);
      if (it != edep.end())
        hit.push_back({it->second, z, r}), eTop = std::max(eTop, it->second);
    }
    printf("%s: %.0f gammas of 2614.5 keV, %zu reached the %.3f kg crystal (highest deposit %.1f keV)\n", inFile[r], N[r], hit.size() - n0, mGe, eTop);
  }

  //---------------------------------------------------------------------------
  //  3. Weight: run r threw N[r] gammas uniformly through a section of volume Vsec. If that section were material m, it would
  //     emit rho[m] th232[m] perChain Vsec gammas a second; so one simulated gamma stands for w = that * year / N[r]
  auto wOf = [&](int m, int r) { return rho[m] * th232[m] * perChain * Vsec * year / N[r]; };
  printf("w per simulated gamma: EFCu %.3e (EFCu run), steel %.3e (steel run) real gammas a year\n", wOf(0, 0), wOf(1, 1));

  //---------------------------------------------------------------------------
  //  4. Reweight: a seam at L; the wall below L is EFCu, above it steel. Counts a year and their statistical error
  auto w = [&](const Hit &h, double L) { return wOf(h.z < L ? 0 : 1, h.run); };
  //     Only hits in the analysis window count: that is where the background index is measured
  auto inWin = [&](const Hit &h) { return fabs(h.E - Qbb) <= winHalf; };
  int nWin[2] = {0, 0};
  for (auto &h : hit)
    nWin[h.run] += inWin(h);
  printf("in the analysis window (%.0f-%.0f keV): %d hits from the EFCu run, %d from the steel run\n", Qbb - winHalf, Qbb + winHalf,
         nWin[0], nWin[1]);
  auto counts = [&](double L, double &err) {
    double c = 0, c2 = 0;
    for (auto &h : hit)
      if (inWin(h))
        c += w(h, L), c2 += w(h, L) * w(h, L);
    err = sqrt(c2);
    return c;
  };
  double eBase, base = counts(L0, eBase);
  const double BI = base / (2 * winHalf * mGe);
  printf("as proposed (seam at %.2f m): %.3e +- %.1e counts a year in the window -> B = %.2e counts/(keV kg yr);\n"
         "  in the %.1f keV ROI that is %.2e counts a year\n\n", L0, base, eBase, BI, roi, BI * roi * mGe);

  //---------------------------------------------------------------------------
  //  5. Scan: move the seam down from 2 m. The change r = (counts at L - counts as proposed) / counts as proposed, and its
  //     statistical error. Both designs use the same hits, so only the hits whose material changed carry the change's noise.
  //     With w0 = w(h, L0), wL = w(h, L), d = wL - w0 (sums over hits):
  //       var(change) = sum d^2,  var(as proposed) = sum w0^2,  cov = sum w0 d
  //       var(r) = [var(change) + r^2 var(as proposed) - 2 r cov] / (as proposed)^2 = sum (d - r w0)^2 / (as proposed)^2
  //     The last form splits var(r) into the share of each run, which tells which run to grow.
  //     A seam passes when r + kSig sigma_r <= tol, fails when r - kSig sigma_r > tol; otherwise the statistics can't tell.
  //     The answer lies between the lowest passing seam and the highest failing one; more gammas narrow that band.
  const double massPerM = rho[0] * Vsec / (H / 2); // kg of EFCu per metre of wall
  double share[2]; // var(r) from the hits of each run, at the last call
  auto change = [&](double L, double &sig) {
    double d = 0, vd = 0, cov = 0;
    for (auto &h : hit)
    {
      if (!inWin(h))
        continue;
      double w0 = w(h, L0), dw = w(h, L) - w0;
      d += dw, vd += dw * dw, cov += w0 * dw;
    }
    double r = d / base;
    sig = sqrt(std::max(0.0, vd + r * r * eBase * eBase - 2 * r * cov)) / base;
    share[0] = share[1] = 0;
    for (auto &h : hit)
      if (inWin(h))
      {
        double w0 = w(h, L0), x = (w(h, L) - w0 - r * w0) / base;
        share[h.run] += x * x;
      }
    return r;
  };
  printf(" seam [m]  EFCu [kg]  counts/yr (window)  change +- error    verdict (line %.0f%%, %.0f sigma)\n", 100 * tol, kSig);
  double Lnew = L0, Lfail = -1;
  bool ok = true;
  for (int k = lround(L0 / dL); k >= 0; k--)
  {
    double L = k * dL, err, c = counts(L, err), sig, r = change(L, sig);
    const bool pass = r + kSig * sig <= tol, fail = r - kSig * sig > tol;
    printf("   %4.2f    %6.1f    %.3e          %+7.1f%% +- %4.1f%%    %s\n", L, L * massPerM, c, 100 * r, 100 * sig,
           pass ? "passes" : fail ? "fails" : "can't tell");
    if (ok && pass)
      Lnew = L;
    else
      ok = false;
    if (fail && Lfail < 0)
      Lfail = L;
  }
  double sNew, rNew = change(Lnew, sNew), cNew = base * (1 + rNew);
  printf("\nnew design: seam at %.2f m -> %.1f kg of EFCu replaced by steel (change %+.1f%% +- %.1f%%)\n", Lnew,
         (L0 - Lnew) * massPerM, 100 * rNew, 100 * sNew);
  printf("as proposed is known to +- %.1f%%\n", 100 * eBase / base);
  const double band = Lnew - std::max(Lfail, 0.0);
  printf("the line is crossed between %.2f m (passes) and %.2f m (fails): a band of %.0f cm\n", Lnew, std::max(Lfail, 0.0), 100 * band);
  if (band <= dRes + 1e-9)
    printf("resolved: the band is within the %.0f cm we asked for\n\n", 100 * dRes);
  else
  {
    // which run limits the answer: the shares of the error at the first seam below the new design
    double sig;
    change(Lnew - dL, sig);
    const double f0 = share[0] / (share[0] + share[1]);
    printf("NOT resolved: the band is wider than %.0f cm. At %.2f m the error is +- %.1f%%, %.0f%% of its variance from the EFCu run\n"
           "  and %.0f%% from the steel run: simulate more gammas mostly in the %s run and run again\n\n", 100 * dRes, Lnew - dL, 100 * sig,
           100 * f0, 100 * (1 - f0), f0 > 0.5 ? "EFCu" : "steel");
  }

  //---------------------------------------------------------------------------
  //  6a. B(E) [counts/(keV kg yr)] for the tube as proposed and the new design
  const double e0 = Qbb - winHalf - eBin * ceil((Qbb - winHalf) / eBin); // lowest bin edge: the window's edges fall on bin edges
  const int nE = ceil((eMax - e0) / eBin);
  auto spec = [&](const char *name, double L) {
    auto h = new TH1D(name, ";energy in the crystal [keV];B [counts / (keV kg yr)]", nE, e0, e0 + nE * eBin);
    for (auto &x : hit)
      h->Fill(x.E, w(x, L));
    h->Scale(1 / (eBin * mGe));
    return h;
  };
  auto hP = spec("hP", L0), hN = spec("hN", Lnew);
  auto c1 = new TCanvas("c1", "", 1400, 760);
  c1->SetMargin(0.09, 0.02, 0.10, 0.08);
  c1->SetLogy();
  double top = std::max(hP->GetMaximum(), hN->GetMaximum()), yMin = top;
  for (auto h : {hP, hN})
    for (int i = 1; i <= nE; i++)
      if (h->GetBinContent(i) > 0)
        yMin = std::min(yMin, 0.5 * h->GetBinContent(i));
  const double yMax = 300 * top;
  hP->SetMinimum(yMin);
  hP->SetMaximum(yMax);
  hP->SetFillColor(TColor::GetColor("#2f74c8"));
  hP->SetLineColor(TColor::GetColor("#1d4f91"));
  hP->SetTitle("2614.5 keV gammas (Tl-208) from the tube wall into 1 kg of Ge-76: no argon, no cuts");
  hP->GetYaxis()->SetTitleOffset(1.1);
  hP->GetXaxis()->SetRangeUser(0, eMax);
  hP->Draw("hist");
  hN->SetLineColor(kRed + 1);
  hN->SetLineWidth(2);
  hN->Draw("hist same");
  // the analysis window (light green) with the B as proposed across it, and the ROI (dark green line at Qbb)
  auto win = new TBox(Qbb - winHalf, yMin, Qbb + winHalf, 1.5 * top);
  win->SetFillColorAlpha(TColor::GetColor("#7ccf7a"), 0.45);
  win->Draw();
  auto bi = new TLine(Qbb - winHalf, BI, Qbb + winHalf, BI);
  bi->SetLineWidth(3);
  bi->Draw();

  auto leg = new TLegend(0.12, 0.77, 0.44, 0.91);
  leg->SetTextSize(0.026);
  leg->AddEntry(hP, Form("as proposed (L = %.2f m)", L0), "f");
  leg->AddEntry(hN, Form("new design (L = %.2f m, %+.0f%% in the window)", Lnew, 100 * (cNew / base - 1)), "l");
  leg->AddEntry(win, Form("analysis window: Q_{#beta#beta} #pm %.0f keV", winHalf), "f");
  leg->AddEntry(bi, Form("B in the window, as proposed: %.1e", BI), "l");
  leg->Draw();
  struct Mark
  {
    double e;
    const char *txt;
  };
  for (auto m : {Mark{2614.5, "2614 full energy"}, Mark{2381.7, "2382 Compton edge"}, Mark{2103.5, "2104 SEP"}, Mark{1592.5, "1592 DEP"},
                 Mark{511, "511"}, Mark{Qbb, Form("ROI %.1f keV at Q_{#beta#beta}", roi)}})
  {
    const bool q = m.e == Qbb;
    auto l = new TLine(m.e, yMin, m.e, 1.5 * top);
    l->SetLineStyle(q ? 1 : 2);
    l->SetLineWidth(q ? 3 : 1);
    l->SetLineColor(q ? TColor::GetColor("#1b5e20") : kGray + 2);
    l->Draw();
    auto t = new TLatex(m.e + (q ? 22 : 12), 2 * top, m.txt);
    t->SetTextAngle(90);
    t->SetTextSize(0.026);
    t->SetTextColor(q ? TColor::GetColor("#1b5e20") : kGray + 3);
    t->Draw();
  }
  c1->SaveAs("tut/toyoutput/toy_background.png");

  //---------------------------------------------------------------------------
  //  6b. Heat map: B(E) for every seam L (rows), as proposed and the new design marked
  const int nL = lround(H / dLmap);
  const double eBin2 = 2 * winHalf, e02 = Qbb - winHalf - eBin2 * ceil((Qbb - winHalf) / eBin2); // the window is one bin
  const int nE2 = ceil((eMax - e02) / eBin2);
  auto hm = new TH2D("hm", "B(E) for every seam L [counts / (keV kg yr)];energy in the crystal [keV];seam L: steel above, EFCu below [m]",
                     nE2, e02, e02 + nE2 * eBin2, nL, 0, H);
  for (int j = 1; j <= nL; j++)
  {
    double L = hm->GetYaxis()->GetBinCenter(j);
    for (auto &x : hit)
      hm->Fill(x.E, L, w(x, L));
  }
  hm->Scale(1 / (eBin2 * mGe));
  hm->GetXaxis()->SetRangeUser(0, eMax);
  auto c2 = new TCanvas("c2", "", 1400, 900);
  c2->SetMargin(0.08, 0.13, 0.10, 0.08);
  c2->SetLogz();
  gStyle->SetPalette(kViridis);
  hm->GetZaxis()->SetTitle("B [counts / (keV kg yr)]");
  hm->GetZaxis()->SetTitleOffset(1.3);
  hm->Draw("colz");
  for (double e : {Qbb - winHalf, Qbb + winHalf})
  {
    auto l = new TLine(e, 0, e, H);
    l->SetLineColor(TColor::GetColor("#7ccf7a"));
    l->SetLineWidth(2);
    l->Draw();
  }
  auto tw2 = new TLatex(Qbb, H + 0.04, Form("analysis window #pm%.0f keV", winHalf));
  tw2->SetTextAlign(21);
  tw2->SetTextSize(0.024);
  tw2->SetTextColor(TColor::GetColor("#3d8b3a"));
  tw2->Draw();
  for (auto [L, txt, col] : {std::tuple{L0, "as proposed", kWhite}, std::tuple{Lnew, "new design (#leq +10%)", kRed}})
  {
    auto l = new TLine(0, L, eMax, L); // seam row
    l->SetLineColor(col);
    l->SetLineWidth(3);
    l->SetLineStyle(L == L0 ? 1 : 2);
    l->Draw();
    auto t = new TPaveText(1500, L + (L == L0 ? 0.05 : -0.27), 2500, L + (L == L0 ? 0.27 : -0.05));
    t->AddText(Form("%s: L = %.2f m", txt, L));
    t->SetFillColor(kWhite);
    t->SetBorderSize(1);
    t->SetTextColor(L == L0 ? kBlack : kRed + 1);
    t->SetTextSize(0.028);
    t->Draw();
  }
  c2->SaveAs("tut/toyoutput/toy_heatmap.png");
  printf("wrote tut/toyoutput/toy_background.png and tut/toyoutput/toy_heatmap.png\n");
}
