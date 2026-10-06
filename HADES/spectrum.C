//  what a calibration source looks like in the HPGe, and what that says about the detector
//    root -l -b -q 'spectrum.C("output/ba133.root")'                    FCCD 1.0 mm: peaks, resolution, FCCD scan, output/ba133_spectrum.png
//    root -l -b -q 'spectrum.C("output/ba133.root", 0.8)'               the spectrum at another FCCD
//    root -l -b -q 'spectrum.C("output/ba133.root", 1.0, 1.52, 0.02)'   a measured observable +- its error -> the detector's FCCD
//  args: file, fccd [mm], measured, error, dlf, seed. the source (nuclide, or gamma line and cone) is read from stp/particles

#include <algorithm>
#include <numeric>
#include <vector>

//-------------------------------------------------------------------------------
//  1. Configuration:
namespace cfg
{
  const double resoA = 0.865, resoB = 0.00225; // PLACEHOLDER resolution FWHM(E) = sqrt(A + B E) [keV], ~1.0 keV at 60 keV, ~2.6 keV at 2615 keV: put the detector's calibration here
  const double window_ns = 1e4;                // one detector event = the steps within 10 us of its first: splits the Th-228 chain, whose members decay days apart
  const std::vector<double> scan = {0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0}; // FCCD [mm] scanned
  const double minFit = 300;                   // net counts a peak needs before its width is fitted

  struct Line { const char *name; double e1, e2, emit; bool fit; };                    // a peak from e1 to e2 keV (a doublet if e2 > e1), gammas per decay (0: not one line), fit its width?
  struct Source { const char *name; int Z, A; std::vector<Line> lines; int num, den; }; // the observable: lines[num] / lines[den], or lines[num] per decay if den < 0
  const std::vector<Source> sources = {
      {"Am-241", 95, 241, {{"59.5", 59.5409, 59.5409, 0.3592, true}, {"26.3", 26.3446, 26.3446, 0.0231, false}}, 0, -1},
      {"Ba-133", 56, 133, {{"79.6+81.0", 79.6142, 80.9979, 0, false}, {"276.4", 276.3989, 276.3989, 0.0716, true}, {"302.9", 302.8508, 302.8508, 0.1834, true},
                           {"356.0", 356.0129, 356.0129, 0.6205, true}, {"383.8", 383.8485, 383.8485, 0.0894, true}}, 0, 3},
      {"Co-60", 27, 60, {{"1173.2", 1173.228, 1173.228, 0.9985, true}, {"1332.5", 1332.492, 1332.492, 0.9998, true}, {"2505.7 sum", 2505.720, 2505.720, 0, false}}, 1, -1},
      {"Th-228", 90, 228, {{"238.6+241.0", 238.632, 240.986, 0, false}, {"583.2", 583.187, 583.187, 0.3045, true}, {"727.3", 727.330, 727.330, 0.0667, true},
                           {"860.6", 860.557, 860.557, 0.0448, true}, {"1592.5 DEP", 1592.511, 1592.511, 0, false}, {"1620.5", 1620.50, 1620.50, 0.0147, true},
                           {"2103.5 SEP", 2103.511, 2103.511, 0, false}, {"2614.5", 2614.511, 2614.511, 0.3585, true}}, 7, -1}}; // Th-228 chain in equilibrium, per Th-228 decay
  double sigma(double e) { return sqrt(resoA + resoB * std::max(e, 0.0)) / 2.355; }
}

static double activeness(double d_mm, double fccd, double dlf) // 0 in the dead part (dlf x fccd), linear through the transition, 1 beyond the FCCD
{
  const double dead = fccd * dlf;
  if (d_mm > fccd) return 1.0;
  if (d_mm <= dead) return 0.0;
  return (d_mm - dead) / (fccd - dead);
}

template <class F> static void scanTree(TTree *t, const std::vector<const char *> &cols, F f) // walk a tree row by row, the named columns as doubles whatever their stored type
{
  t->SetBranchStatus("*", 0);
  for (auto c : cols) t->SetBranchStatus(c, 1);
  std::vector<TTreeFormula *> col;
  for (auto c : cols) col.push_back(new TTreeFormula(c, c, t));
  std::vector<double> v(cols.size());
  for (Long64_t i = 0; i < t->GetEntries(); i++)
  {
    t->GetEntry(i);
    for (size_t k = 0; k < col.size(); k++) v[k] = col[k]->EvalInstance();
    f(v.data());
  }
  for (auto c : col) delete c;
}

void spectrum(const char *file, double fccd = 1.0, double meas = -1, double measErr = 0, double dlf = 0.5, int seed = 1)
{
//-------------------------------------------------------------------------------
//  2. The source, from the primaries remage stored:
  auto f = TFile::Open(file);
  auto d = (f && !f->IsZombie()) ? (TDirectory *)f->Get("stp") : nullptr;
  auto tp = d ? (TTree *)d->Get("particles") : nullptr;
  if (!tp) { printf("ERROR: no stp/particles in %s - simulate with run.mac, which stores the primaries\n", file); return; }
  const Long64_t nev = tp->GetEntries();
  int pdg = 0;
  double ekin = 0, cosMin = 1; // the primary's energy [keV] and the widest emission angle from -z, where the crystal is
  scanTree(tp, {"particle", "px_in_MeV", "py_in_MeV", "pz_in_MeV", "ekin_in_MeV"}, [&](const double *v) {
    pdg = (int)v[0];
    ekin = v[4] * 1000;
    const double p = sqrt(v[1] * v[1] + v[2] * v[2] + v[3] * v[3]);
    if (p > 0) cosMin = std::min(cosMin, -v[3] / p);
  });
  cfg::Source src{};
  std::string lineName; // a line-mode source's name: it must outlive src
  double norm = nev;    // what one count is divided by: decays, or for a line the decays that emit it into 4 pi
  TString what;
  if (pdg == 22) // line mode: a gamma aimed into a cone
  {
    const double cone = std::max(1e-12, (1 - cosMin) / 2), deg = acos(cosMin) * 180 / M_PI;
    lineName = Form("%.1f", ekin);
    src = {"gamma", 0, 0, {{lineName.c_str(), ekin, ekin, 0, true}}, 0, -1};
    for (auto &s : cfg::sources)
      for (auto &l : s.lines)
        if (l.emit > 0 && fabs(l.e1 - ekin) < 0.05) { src.name = s.name; src.lines[0] = l; }
    const double emit = src.lines[0].emit;
    norm = nev / cone / (emit > 0 ? emit : 1);
    const TString per = emit > 0 ? TString::Format(" = %.3g decays", norm) : TString(", counts per gamma (no emission probability for this line)");
    what = TString::Format("%s %.4g keV line in a %.2f deg cone: %lld gammas = %.3g in 4 pi%s", src.name, ekin, deg, nev, nev / cone, per.Data());
  }
  else if (pdg > 1000000000) // decay mode: the nuclide itself
  {
    const int Z = pdg / 10000 % 1000, A = pdg / 10 % 1000;
    for (auto &s : cfg::sources)
      if (s.Z == Z && s.A == A) src = s;
    if (!src.name) { printf("ERROR: no lines for Z %d A %d in cfg::sources\n", Z, A); return; }
    what = TString::Format("%s: %lld decays", src.name, nev);
  }
  else { printf("ERROR: primary PDG %d is neither a gamma nor a nucleus\n", pdg); return; }
  const bool ratio = src.den >= 0;
  TString obsName = ratio ? TString::Format("R = (%s) / (%s)", src.lines[src.num].name, src.lines[src.den].name) : TString::Format("%s per decay", src.lines[src.num].name);
  if (pdg == 22 && src.lines[0].emit <= 0) obsName = TString::Format("%s per gamma", src.lines[0].name);

//-------------------------------------------------------------------------------
//  3. Steps -> detector events, deposited and active energy at every FCCD in one pass:
  TTree *tg = nullptr;
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it()) // the germanium table is the one carrying dist_to_surf
  {
    auto c = (TTree *)d->Get(k->GetName());
    if (c && c->GetBranch("dist_to_surf_in_m")) { tg = c; break; }
  }
  if (!tg) { printf("ERROR: no germanium table in %s\n", file); return; }
  std::vector<int> sev;
  std::vector<float> sed, sds;
  std::vector<double> stm;
  scanTree(tg, {"evtid", "edep_in_keV", "dist_to_surf_in_m", "time_in_ns"}, [&](const double *v) {
    sev.push_back((int)v[0]); sed.push_back(v[1]); sds.push_back(v[2] * 1000); stm.push_back(v[3]); // m -> mm
  });
  std::vector<double> F = cfg::scan; // the scan, with the requested FCCD in it
  if (std::find(F.begin(), F.end(), fccd) == F.end()) { F.push_back(fccd); std::sort(F.begin(), F.end()); }
  const int NF = F.size(), k0 = std::find(F.begin(), F.end(), fccd) - F.begin();
  std::vector<size_t> idx(sev.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return sev[a] != sev[b] ? sev[a] < sev[b] : stm[a] < stm[b]; }); // with -t an event's rows need not be together
  std::vector<float> E; // per detector event: deposited, then active at each FCCD
  int cur = -1;
  double t0 = 0;
  for (size_t j : idx)
  {
    if (sev[j] != cur || stm[j] - t0 > cfg::window_ns) { E.insert(E.end(), NF + 1, 0.f); cur = sev[j]; t0 = stm[j]; }
    float *e = &E[E.size() - NF - 1];
    e[0] += sed[j];
    for (int k = 0; k < NF; k++) e[k + 1] += sed[j] * activeness(sds[j], F[k], dlf);
  }
  const size_t nE = E.size() / (NF + 1);
  printf("%s: %s\n  %zu steps -> %zu detector events (%.0f us window)\n", file, what.Data(), sev.size(), nE, cfg::window_ns / 1000);

//-------------------------------------------------------------------------------
//  4. Resolution, then peak counts: a window of +-3 sigma, the background from a side band of half its width on each side:
  gRandom->SetSeed(seed);
  const int NL = src.lines.size();
  std::vector<double> lo(NL), hi(NL);
  for (int l = 0; l < NL; l++) { lo[l] = src.lines[l].e1 - 3 * cfg::sigma(src.lines[l].e1); hi[l] = src.lines[l].e2 + 3 * cfg::sigma(src.lines[l].e2); }
  std::vector<std::vector<double>> nS(NF, std::vector<double>(NL, 0)), nB = nS; // per FCCD, per line: window, side bands
  double emax = 0, sRaw = 0, sAct = 0;
  for (size_t i = 0; i < nE; i++) emax = std::max(emax, (double)E[i * (NF + 1)]);
  emax = std::max(emax, src.lines.back().e2) * 1.03;
  const double bw = emax < 600 ? 0.1 : 0.25;      // fine enough for the width fits
  const int rb = emax < 150 ? 2 : (emax < 600 ? 5 : 4); // display bins: 0.2, 0.5 or 1 keV
  const int nb = ((int)(emax / bw) / rb + 1) * rb;  // a whole number of display bins
  auto hR = new TH1D("hR", "", nb, 0, nb * bw), hA = new TH1D("hA", "", nb, 0, nb * bw), hS = new TH1D("hS", "", nb, 0, nb * bw);
  for (size_t i = 0; i < nE; i++)
  {
    const float *e = &E[i * (NF + 1)];
    const double z = gRandom->Gaus(0, 1); // one draw per event, shared by every FCCD: the scan moves smoothly
    for (int k = 0; k < NF; k++)
    {
      const double es = e[k + 1] + z * cfg::sigma(e[k + 1]);
      for (int l = 0; l < NL; l++)
      {
        const double w = hi[l] - lo[l];
        if (es >= lo[l] && es < hi[l]) nS[k][l]++;
        else if ((es >= lo[l] - w / 2 && es < lo[l]) || (es >= hi[l] && es < hi[l] + w / 2)) nB[k][l]++;
      }
      if (k == k0) { hS->Fill(es); hA->Fill(e[k + 1]); }
    }
    hR->Fill(e[0]);
    sRaw += e[0];
    sAct += e[k0 + 1];
  }
  auto net = [&](int k, int l) { return nS[k][l] - nB[k][l]; };
  auto netErr = [&](int k, int l) { return sqrt(nS[k][l] + nB[k][l]); };

//-------------------------------------------------------------------------------
//  5. The response at the requested FCCD: peaks, efficiencies, widths:
  printf("\nFCCD %.2f mm, DLF %.2f, seed %d: energy lost to the dead layer %.2f %%\n", fccd, dlf, seed, sRaw > 0 ? 100 * (1 - sAct / sRaw) : 0.0);
  printf("  %-12s %-17s %18s %22s   %s\n", "peak [keV]", "window [keV]", "net counts", pdg == 22 && src.lines[0].emit <= 0 ? "per gamma" : "per decay", "FWHM [keV]: fitted (injected)");
  auto gW = new TGraphErrors(); // fitted FWHM^2 against energy: the resolution the run produced
  double eLo = 1e9, eHi = 0;    // the energies it spans
  for (int l = 0; l < NL; l++)
  {
    TString fw = "-";
    const auto &L = src.lines[l];
    if (L.fit && net(k0, l) > cfg::minFit)
    {
      const double s = cfg::sigma(L.e1);
      auto g = new TF1(Form("g%d", l), "gaus(0)+[3]+[4]*0.5*TMath::Erfc((x-[1])/(sqrt(2)*[2]))", L.e1 - 4 * s, L.e1 + 4 * s); // a Gaussian on a step: the level above the peak, plus what scatters below it. chi2: empty bins above a lone line break a likelihood
      const double above = hS->GetBinContent(hS->FindBin(L.e1 + 4 * s)), below = hS->GetBinContent(hS->FindBin(L.e1 - 4 * s));
      g->SetParameters(hS->GetBinContent(hS->FindBin(L.e1)), L.e1, s, above, std::max(0.0, below - above));
      g->SetParLimits(2, 0.3 * s, 3 * s);
      if ((int)hS->Fit(g, "QRNS") != 0) fw = "fit did not converge";
      else
      {
        fw = Form("%.3f +- %.3f (%.3f)", 2.355 * g->GetParameter(2), 2.355 * g->GetParError(2), 2.355 * s);
        gW->SetPoint(gW->GetN(), L.e1, pow(2.355 * g->GetParameter(2), 2));
        gW->SetPointError(gW->GetN() - 1, 0, 2 * pow(2.355, 2) * g->GetParameter(2) * g->GetParError(2));
        eLo = std::min(eLo, L.e1); eHi = std::max(eHi, L.e1);
      }
    }
    printf("  %-12s %7.2f - %-7.2f %10.0f +- %-6.0f %10.3e +- %-8.1e   %s\n", L.name, lo[l], hi[l], net(k0, l), netErr(k0, l), net(k0, l) / norm, netErr(k0, l) / norm, fw.Data());
  }
  if (gW->GetN() >= 2 && eHi > 3 * eLo) // A and B separate only over a wide span
  {
    gW->Fit("pol1", "Q");
    auto p = gW->GetFunction("pol1");
    printf("  resolution from %d peaks, FWHM^2 = A + B E: A %.3f +- %.3f keV^2 (injected %.3f), B %.3e +- %.1e keV (injected %.3e)\n", gW->GetN(),
           p->GetParameter(0), p->GetParError(0), cfg::resoA, p->GetParameter(1), p->GetParError(1), cfg::resoB);
  }

//-------------------------------------------------------------------------------
//  6. The FCCD scan: the observable against the dead layer, and a measurement turned into an FCCD:
  std::vector<double> obs(NF), dobs(NF);
  for (int k = 0; k < NF; k++)
  {
    const double a = net(k, src.num), da = netErr(k, src.num);
    if (ratio)
    {
      const double b = net(k, src.den), db = netErr(k, src.den);
      obs[k] = b > 0 ? a / b : 0;
      dobs[k] = b > 0 && a > 0 ? obs[k] * sqrt(pow(da / a, 2) + pow(db / b, 2)) : 0;
    }
    else { obs[k] = a / norm; dobs[k] = da / norm; }
  }
  printf("\nFCCD scan, DLF %.2f: %s against the dead layer\n  %-10s", dlf, obsName.Data(), "FCCD [mm]");
  if (ratio) printf(" %20s %20s", Form("%s /decay", src.lines[src.num].name), Form("%s /decay", src.lines[src.den].name));
  printf(" %24s\n", ratio ? "R" : obsName.Data());
  for (int k = 0; k < NF; k++)
  {
    printf("  %-10.2f", F[k]);
    if (ratio) printf(" %20.4e %20.4e", net(k, src.num) / norm, net(k, src.den) / norm);
    printf(" %13.5g +- %-8.2g%s\n", obs[k], dobs[k], k == k0 ? "   <- drawn" : "");
  }
  double fMeas = -1, fErr = 0;
  if (meas > 0)
  {
    for (int k = 0; k + 1 < NF; k++)
      if ((obs[k] - meas) * (obs[k + 1] - meas) <= 0 && obs[k] != obs[k + 1]) // the measurement falls in this interval
      {
        const double slope = (obs[k + 1] - obs[k]) / (F[k + 1] - F[k]);
        const double mc = dobs[k] + (dobs[k + 1] - dobs[k]) * (meas - obs[k]) / (obs[k + 1] - obs[k]);
        fMeas = F[k] + (meas - obs[k]) / slope;
        fErr = sqrt(measErr * measErr + mc * mc) / fabs(slope);
        printf("\nmeasured %.5g +- %.2g  ->  FCCD %.3f +- %.3f mm   (measurement %.3f, MC statistics %.3f)\n", meas, measErr, fMeas, fErr,
               measErr / fabs(slope), mc / fabs(slope));
        break;
      }
    if (fMeas < 0) printf("\nmeasured %.5g lies outside the scan (%.5g to %.5g): extend cfg::scan\n", meas, obs.front(), obs.back());
  }

//-------------------------------------------------------------------------------
//  7. Draw: the spectrum at the requested FCCD, and the scan:
  gStyle->SetOptStat(0);
  auto c = new TCanvas("c", "", 1500, 560);
  c->Divide(2, 1);
  auto p1 = c->cd(1);
  p1->SetPad(0, 0, 0.64, 1);
  p1->SetLogy(); p1->SetGrid(); p1->SetLeftMargin(0.1);
  for (auto h : {hR, hA, hS}) h->Rebin(rb);
  hR->SetLineColor(kGray + 2); hA->SetLineColor(kAzure + 2); hS->SetLineColor(kRed + 1);
  hR->SetTitle(Form("%s;energy [keV];counts / %.1f keV", what.Data(), bw * rb));
  hR->GetYaxis()->SetRangeUser(0.5, 3 * hR->GetMaximum());
  hR->Draw("hist"); hA->Draw("hist same"); hS->Draw("hist same");
  auto leg = new TLegend(0.6, 0.74, 0.89, 0.89);
  leg->AddEntry(hR, "deposited", "l");
  leg->AddEntry(hA, Form("after the dead layer (FCCD %.2f mm)", fccd), "l");
  leg->AddEntry(hS, "after the resolution", "l");
  leg->SetBorderSize(0); leg->Draw();
  auto p2 = c->cd(2);
  p2->SetPad(0.64, 0, 1, 1);
  p2->SetGrid(); p2->SetLeftMargin(0.2);
  auto gs = new TGraphErrors(NF, F.data(), obs.data(), nullptr, dobs.data());
  gs->SetTitle(Form("FCCD scan, DLF %.2f;FCCD [mm];%s", dlf, obsName.Data()));
  gs->SetMarkerStyle(20); gs->SetMarkerColor(kAzure + 2); gs->SetLineColor(kAzure + 2);
  gs->Draw("APL");
  gs->GetYaxis()->SetTitleOffset(2.4);
  if (fMeas >= 0)
  {
    auto band = new TBox(F.front(), meas - measErr, F.back(), meas + measErr);
    band->SetFillColorAlpha(kRed, 0.2);
    band->Draw();
    auto m = new TMarker(fMeas, meas, 29);
    m->SetMarkerColor(kRed + 1); m->SetMarkerSize(2);
    m->Draw();
  }
  TString png = TString(file).ReplaceAll(".root", "") + "_spectrum.png";
  c->SaveAs(png);
  printf("wrote %s\n", png.Data());
}
