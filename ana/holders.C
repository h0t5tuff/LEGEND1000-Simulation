//  1. Configuration:
namespace cfg
{
  const double window_ns = 1e4;                                                                                    // one detector event: its steps within 10 us of the first
  const double eMax = 3000, binW = 5;                                                                              // keV
  const char *chain[2] = {"Th-232", "U-238"}, *fill[2] = {"#2f74c8", "#e8862a"}, *ink[2] = {"#1d4f91", "#b35900"}; // spectrum, label
  // the ROI: Qbb +- 2 sigma, sigma from the resolution the tube's response uses (legend1000-metadata simprod l1000dsg01
  // eresmod: FWHM = sqrt(0.5 + 0.001 E) keV, 1.59 keV at Qbb)
  const double Qbb = 2039.0, roiHalf = 2 * sqrt(0.5 + 0.001 * Qbb) / 2.35482;
  struct Line
  {
    double e;
    const char *keV, *name;
    int ch;
    double xLab, yLab;
  }; // label centre [keV], foot [log fraction of the frame]
  const Line lines[] = {{238.6, "238", "Pb-212", 0, 200, 0.90}, {351.9, "352", "Pb-214", 1, 390, 0.79}, {583.2, "583", "Tl-208", 0, 590, 0.90}, {609.3, "609", "Bi-214", 1, 790, 0.79}, {911.2, "911", "Ac-228", 0, 1000, 0.90}, {1120.3, "1120", "Bi-214", 1, 1220, 0.72}, {1764.5, "1764", "Bi-214", 1, 1764, 0.70}, {2447.9, "2448", "Bi-214", 1, 2448, 0.50}, {2614.5, "2614", "Tl-208", 0, 2614, 0.74}};
}

//-------------------------------------------------------------------------------
//  2. One run: its steps, sorted by event, detector and time, then summed into detector events
static Long64_t readRun(const char *file, TH1D *h)
{
  TFile f(file);
  auto t = f.IsZombie() ? nullptr : (TTree *)f.Get("stp/germanium");
  auto v = f.IsZombie() ? nullptr : (TTree *)f.Get("stp/vtx");
  if (!t || !v)
  {
    printf("ERROR: %s has no stp/germanium or stp/vtx\n", file);
    return 0;
  }
  Int_t ev, det;
  Float_t e;
  Double_t tm;
  t->SetBranchStatus("*", 0);
  for (auto b : {"evtid", "det_uid", "edep_in_keV", "time_in_ns"})
    t->SetBranchStatus(b, 1);
  t->SetBranchAddress("evtid", &ev);
  t->SetBranchAddress("det_uid", &det);
  t->SetBranchAddress("edep_in_keV", &e);
  t->SetBranchAddress("time_in_ns", &tm);
  const Long64_t n = t->GetEntries();
  std::vector<Int_t> sev(n), sdet(n);
  std::vector<Float_t> sE(n);
  std::vector<Double_t> sT(n);
  for (Long64_t i = 0; i < n; i++)
  {
    t->GetEntry(i);
    sev[i] = ev;
    sdet[i] = det;
    sE[i] = e;
    sT[i] = tm;
  }
  std::vector<Long64_t> o(n); // the threads write their events interleaved
  std::iota(o.begin(), o.end(), 0);
  std::sort(o.begin(), o.end(), [&](Long64_t a, Long64_t b)
            { return std::tie(sev[a], sdet[a], sT[a]) < std::tie(sev[b], sdet[b], sT[b]); });
  Long64_t nEv = 0, cur = -1;
  double sum = 0, t0 = 0;
  for (Long64_t i : o)
  {
    if (cur >= 0 && (sev[i] != sev[cur] || sdet[i] != sdet[cur] || sT[i] - t0 > cfg::window_ns))
    {
      h->Fill(sum);
      nEv++;
      cur = -1;
    }
    if (cur < 0)
    {
      cur = i;
      t0 = sT[i];
      sum = 0;
    }
    sum += sE[i];
  }
  if (cur >= 0)
  {
    h->Fill(sum);
    nEv++;
  }
  printf("%s: %lld chains, %lld germanium steps -> %lld detector events\n", file, v->GetEntries(), n, nEv);
  return v->GetEntries();
}

//-------------------------------------------------------------------------------
//  3. The spectrum, its lines named
void holders(const char *files = "output/holders_th232.root,output/holders_u238.root", const char *png = "output/holders_spectrum.png")
{
  gErrorIgnoreLevel = kWarning; // no "png file has been created" notice
  gStyle->SetOptStat(0);
  TH1D *h[2];
  Long64_t nCh[2] = {0, 0};
  for (int k = 0; k < 2; k++)
    h[k] = new TH1D(Form("h%d", k), ";energy per detector per event [keV];counts / 5 keV", int(cfg::eMax / cfg::binW), 0, cfg::eMax);
  std::stringstream list(files);
  for (std::string f; std::getline(list, f, ',');) // the chain from the file name
  {
    int k = f.find("th232") != std::string::npos ? 0 : f.find("u238") != std::string::npos ? 1
                                                                                           : -1;
    if (k < 0)
    {
      printf("ERROR: %s: neither th232 nor u238 in its name\n", f.c_str());
      continue;
    }
    nCh[k] += readRun(f.c_str(), h[k]);
  }
  if (!nCh[0] && !nCh[1])
    return;

  auto c = new TCanvas("c", "", 1400, 760);
  c->SetMargin(0.09, 0.02, 0.10, 0.08);
  c->SetLogy();
  const double yMin = 0.5, yMax = 50 * std::max(h[0]->GetMaximum(), h[1]->GetMaximum());
  auto yAt = [&](double frac)
  { return yMin * pow(yMax / yMin, frac); };               // a fraction of the log frame
  auto ndc = [&](double x, double frac, const char *txt) { // text placed in the pad's own coordinates: subscripts misplace on a log axis
    auto t = new TLatex(c->GetLeftMargin() + x / cfg::eMax * (1 - c->GetLeftMargin() - c->GetRightMargin()),
                        c->GetBottomMargin() + frac * (1 - c->GetBottomMargin() - c->GetTopMargin()), txt);
    t->SetNDC();
    t->SetTextFont(62);
    t->SetTextSize(0.026);
    return t;
  };
  h[0]->SetMinimum(yMin);
  h[0]->SetMaximum(yMax);
  h[0]->GetYaxis()->SetTitleOffset(1.1);
  h[0]->Draw("axis");

  auto leg = new TLegend(0.79, 0.82, 0.975, 0.91);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.03);
  for (int k = 0; k < 2; k++)
  {
    if (!nCh[k])
      continue;
    const int col = TColor::GetColor(cfg::fill[k]);
    h[k]->SetFillColorAlpha(col, 0.55);
    h[k]->SetLineColor(col);
    h[k]->Draw("hist same");
    leg->AddEntry(h[k], Form("%s chain (%.0f k)", cfg::chain[k], nCh[k] / 1e3), "f");
  }
  // the ROI over the spectrum, and Qbb
  const int roiCol = TColor::GetColor("#2e7d32");
  auto b = new TBox(cfg::Qbb - cfg::roiHalf, yMin, cfg::Qbb + cfg::roiHalf, yMax);
  b->SetFillColorAlpha(roiCol, 0.8);
  b->SetLineWidth(0);
  b->Draw(); // 2.7 keV: a line at this scale
  auto roiText = ndc(cfg::Qbb + 15, 0.62, Form("#splitline{ROI Q_{#beta#beta} #pm 2#sigma}{%.1f-%.1f keV}", cfg::Qbb - cfg::roiHalf, cfg::Qbb + cfg::roiHalf));
  roiText->SetTextAlign(11);
  roiText->SetTextColor(roiCol);
  roiText->Draw();
  gPad->RedrawAxis();
  leg->Draw();
  auto title = new TLatex(0.5, 0.95, "Radioactivity of the EFCu detector holders: Th-232 and U-238 chains");
  title->SetNDC();
  title->SetTextAlign(22);
  title->SetTextSize(0.04);
  title->Draw();

  for (const auto &L : cfg::lines) // each line, in its chain's colour
  {
    const int ink = TColor::GetColor(cfg::ink[L.ch]);
    const double top = h[L.ch]->GetBinContent(h[L.ch]->FindBin(L.e)), yLab = yAt(L.yLab);
    auto a = new TArrow(L.xLab, yLab / 1.35, L.e, top * 1.25, 0.008, "|>");
    a->SetLineColor(ink);
    a->SetFillColor(ink);
    a->SetLineWidth(2);
    a->Draw();
    auto t = ndc(L.xLab, L.yLab, Form("#splitline{%s keV}{%s (%s)}", L.keV, L.name, cfg::chain[L.ch]));
    t->SetTextAlign(21);
    t->SetTextColor(ink);
    t->Draw();
  }
  c->SaveAs(png);
  printf("wrote %s\n", png);
}
