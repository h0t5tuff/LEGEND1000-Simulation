//  what a germanium hit is made of: a pile of Geant4 steps becomes ONE detector response,
//  and each response traces back to the gamma that caused it
//    root -l -b -q 'ana/hits.C("output/tl208.root")'

#include "../geom/rt.h"

struct Step { int ev, det, pdg, tid, par; float e, x, y, z, dsurf, t; };
const int   pdgOf[5] = {22, 11, -11, 2112, 0}; // the particles that deposit in germanium
const char *pdgName[5] = {"gamma", "e-", "e+", "neutron", "other"};
static int slot(int pdg) { int i = 0; while (i < 4 && pdg != pdgOf[i]) i++; return i; } // no match leaves it on "other"

void hits(const char *fn, double m1_keV = 5.0)
{
  TFile f(fn);
  auto d = f.IsZombie() ? nullptr : (TDirectory *)f.Get("stp");
  if (!d)
  {
    printf("ERROR: no stp/ in %s\n", fn);
    return;
  }

//-------------------------------------------------------------------------------
//  1. Read every step:
  std::vector<Step> steps;
  std::vector<TString> detName;
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it())
  {
    auto t = rtIsGermanium(k->GetName()) ? (TTree *)d->Get(k->GetName()) : nullptr;
    if (!t || !t->GetBranch("edep_in_keV"))
      continue;
    int id = detName.size();
    detName.push_back(k->GetName());
    bool trk = t->GetBranch("trackid"); // only there if the run stored track ids
    rtScan(t, {"evtid", "particle", trk ? "trackid" : "evtid", trk ? "parent_trackid" : "evtid", "edep_in_keV",
               "xloc_in_m", "yloc_in_m", "zloc_in_m", "dist_to_surf_in_m", "time_in_ns"}, [&](const double *v) {
      steps.push_back({(int)v[0], id, (int)v[1], trk ? (int)v[2] : -1, trk ? (int)v[3] : -1, (float)v[4],
                       (float)(v[5] * 1000), (float)(v[6] * 1000), (float)(v[7] * 1000), (float)(v[8] * 1000), (float)v[9]}); // m -> mm
    });
  }
  if (steps.empty())
  {
    printf("no germanium steps in %s\n", fn);
    return;
  }
  std::map<std::pair<int, int>, double> gamma; // the gamma catalogue: (decay, trackid) -> energy at creation [keV]
  if (auto tr = (TTree *)d->Get("tracks"))
    rtScan(tr, {"evtid", "trackid", "particle", "ekin_in_MeV"}, [&](const double *v) { if ((int)v[2] == 22) gamma[{(int)v[0], (int)v[1]}] = v[3] * 1000; });
  std::map<int, double> larE; // what the same decays left in the underground argon
  if (auto tl = (TTree *)d->Get("undergroundlar"))
    rtScan(tl, {"evtid", "edep_in_keV"}, [&](const double *v) { larE[(int)v[0]] += v[1]; });

//-------------------------------------------------------------------------------
//  2. Collapse steps into detector responses:
  struct Resp { double e = 0; int n = 0; };
  std::map<std::pair<int, int>, Resp> resp; // key = (decay, detector), so two detectors in one decay stay separate
  double eBy[5] = {0}, eAll = 0, eTot = 0;
  long nBy[5] = {0}, nAbove = 0;
  for (auto &s : steps)
  {
    auto &r = resp[{s.ev, s.det}];
    r.e += s.e;
    r.n++;
    eBy[slot(s.pdg)] += s.e;
    nBy[slot(s.pdg)]++;
  }
  std::set<int> evts;
  for (auto &p : resp)
    if (p.second.e > m1_keV) { nAbove++; evts.insert(p.first.first); eTot += p.second.e; }
  printf("file   : %s\n", fn);
  printf("\n[2] energy deposition  ->  detector response\n");
  printf("      %-34s %10ld\n", "Geant4 steps in germanium", (long)steps.size());
  printf("      %-34s %10zu\n", "detector responses (evt x det)", resp.size());
  printf("      %-34s %10ld   (above %.0f keV)\n", "responses kept", nAbove, m1_keV);
  printf("      %-34s %10zu\n", "decays that produced one", evts.size());
  printf("      -> on average %.1f steps collapse into one number per detector\n", (double)steps.size() / resp.size());

//-------------------------------------------------------------------------------
//  3. Which particle deposits the energy (gammas carry it in, the electrons they free deposit it: a tight cluster, not a point):
  printf("\n[3] which particle actually deposits the energy:\n");
  printf("      %-9s %12s %14s %10s\n", "particle", "steps", "energy [keV]", "share");
  for (int i = 0; i < 5; i++) eAll += eBy[i];
  for (int i = 0; i < 5; i++)
    if (nBy[i])
      printf("      %-9s %12ld %14.1f %9.1f%%\n", pdgName[i], nBy[i], eBy[i], eAll ? 100 * eBy[i] / eAll : 0.0);

//-------------------------------------------------------------------------------
//  4. One response, step by step:
  std::pair<int, int> show;
  int bestN = 0;
  for (auto &p : resp) // the response above 100 keV with the most steps
    if (p.second.e > 100 && p.second.n > bestN) { bestN = p.second.n; show = p.first; }
  printf("\n[4] one detector response, step by step (decay %d in %s):\n", show.first, detName[show.second].Data());
  printf("      %-9s %12s %12s %12s %12s %12s\n", "particle", "edep [keV]", "x [mm]", "y [mm]", "z [mm]", "depth [mm]");
  double run = 0;
  std::vector<Step> shown;
  for (auto &s : steps)
    if (s.ev == show.first && s.det == show.second)
    {
      run += s.e;
      shown.push_back(s);
      if (shown.size() <= 12)
        printf("      %-9s %12.3f %12.3f %12.3f %12.3f %12.4f\n", pdgName[slot(s.pdg)], s.e, s.x, s.y, s.z, s.dsurf);
    }
  if (shown.size() > 12)
    printf("      ... %zu more steps\n", shown.size() - 12);
  printf("      %-9s %12.3f  <- what the detector reports\n", "SUM", run);

//-------------------------------------------------------------------------------
//  5. Trace each response back to its gamma:
  printf("\n[5] tracing responses back to the gamma that caused them:\n");
  printf("      %-7s %-8s %11s %10s %12s %12s\n", "decay", "detector", "Ge [keV]", "t [ns]", "gamma [keV]", "LAr [keV]");
  long resolved = 0;
  std::map<std::pair<int, int>, std::map<double, double>> byGamma; // response -> gamma energy -> energy it delivered
  std::map<std::pair<int, int>, double> respTime;                  // earliest deposit in that detector
  for (auto &s : steps)
  {
    auto key = std::make_pair(s.ev, s.det);
    auto g = gamma.find({s.ev, s.pdg == 22 ? s.tid : s.par}); // an electron's parent is the gamma that freed it; a gamma is its own
    if (g != gamma.end() && g->second > 0) { resolved++; byGamma[key][g->second] += s.e; }
    if (!respTime.count(key) || s.t < respTime[key]) respTime[key] = s.t;
  }
  int listed = 0;
  for (auto &p : resp)
  {
    if (p.second.e <= m1_keV || listed++ >= 10)
      continue;
    double best = -1, bestE = -1; // credit the response to whichever gamma delivered the most of its energy
    for (auto &g : byGamma[p.first])
      if (g.second > bestE) { bestE = g.second; best = g.first; }
    printf("      %-7d %-8s %11.2f %10.2f %12s %12.2f\n", p.first.first, detName[p.first.second].Data(), p.second.e,
           respTime[p.first], best > 0 ? Form("%.1f", best) : "(secondary)", larE.count(p.first.first) ? larE[p.first.first] : 0.0);
  }
  printf("      %.0f%% of deposits trace to a named gamma; the rest are later generations\n", 100.0 * resolved / steps.size());
  double larTot = 0;
  for (auto &p : larE) larTot += p.second;
  printf("      argon saw %.3e keV over %zu decays, germanium %.3e keV over %zu responses\n", larTot, larE.size(), eTot, resp.size());

//-------------------------------------------------------------------------------
//  6. Draw it:
  gStyle->SetOptStat(0);
  auto hRZ = new TH2D("hRZ", "where the energy lands;r [mm];z [mm]", 90, 0, 900, 90, -1000, 1500);
  auto hXY = new TH2D("hXY", "same hits from above;x [mm];y [mm]", 90, -900, 900, 90, -900, 900);
  auto hStep = new TH1D("hStep", "deposition vs response;energy [keV];counts", 120, 0, 3000);
  auto hResp = new TH1D("hResp", "", 120, 0, 3000);
  auto hN = new TH1D("hN", "steps per detector response;steps;responses", 30, 0.5, 30.5);
  auto hD = new TH1D("hD", "how deep into the crystal;distance to surface [mm];steps", 60, 0, 30);
  for (auto &s : steps)
  {
    hRZ->Fill(sqrt(s.x * s.x + s.y * s.y), s.z, s.e);
    hXY->Fill(s.x, s.y, s.e);
    hStep->Fill(s.e);
    hD->Fill(s.dsurf);
  }
  for (auto &p : resp)
    if (p.second.e > m1_keV) { hResp->Fill(p.second.e); hN->Fill(p.second.n); }
  auto c = new TCanvas("c_hits", "", 1400, 850);
  c->Divide(3, 2);
  struct { TH1 *h; int col; bool log2d; } pad[5] = {{hRZ, 0, true}, {hXY, 0, true}, {hStep, kGray + 2, false}, {hN, kOrange + 7, false}, {hD, kGreen + 2, false}};
  for (int i = 0; i < 5; i++)
  {
    c->cd(i + 1);
    if (pad[i].log2d) { gPad->SetLogz(); gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.15); }
    else { gPad->SetLogy(); gPad->SetGrid(); pad[i].h->SetLineColor(pad[i].col); pad[i].h->SetLineWidth(2); }
    pad[i].h->Draw(pad[i].log2d ? "COLZ" : "HIST");
  }
  c->cd(1);
  if (auto o = (TTree *)d->Get("detector_origins")) // the array itself, for context
  {
    auto g = new TGraph();
    rtScan(o, {"xloc_in_m", "yloc_in_m", "zloc_in_m"}, [&](const double *v) { g->SetPoint(g->GetN(), 1000 * sqrt(v[0] * v[0] + v[1] * v[1]), 1000 * v[2]); });
    g->SetMarkerColor(kGray + 1);
    g->Draw("P SAME");
  }
  c->cd(3);
  hResp->SetLineColor(kAzure + 2);
  hResp->SetLineWidth(2);
  hResp->Draw("HIST SAME");
  auto lg = new TLegend(0.45, 0.75, 0.88, 0.88);
  lg->SetTextSize(0.035);
  lg->AddEntry(hStep, "single Geant4 steps", "l");
  lg->AddEntry(hResp, "summed detector response", "l");
  lg->Draw();
  c->cd(6)->SetGrid(); // the one response from [4], framed on its own extent: a fixed window would crop it or hide it
  double x0 = 1e9, x1 = -1e9, z0 = 1e9, z1 = -1e9;
  for (auto &s : shown)
  {
    x0 = std::min(x0, (double)s.x); x1 = std::max(x1, (double)s.x);
    z0 = std::min(z0, (double)s.z); z1 = std::max(z1, (double)s.z);
  }
  double m = 0.12 * std::max(1.0, std::max(x1 - x0, z1 - z0));
  gPad->DrawFrame(x0 - m, z0 - m, x1 + m, z1 + m)->SetTitle(Form("decay %d in %s: %zu steps -> %.0f keV;x [mm];z [mm]", show.first, detName[show.second].Data(), shown.size(), run));
  auto lg2 = new TLegend(0.60, 0.78, 0.90, 0.90);
  lg2->SetTextSize(0.035);
  for (int isGamma : {1, 0}) // open red = gamma, full blue = electron. marker area follows the energy
  {
    for (auto &s : shown)
      if ((s.pdg == 22) == isGamma)
      {
        auto mk = new TMarker(s.x, s.z, isGamma ? 24 : 20);
        mk->SetMarkerColor(isGamma ? kRed + 1 : kAzure + 2);
        mk->SetMarkerSize(0.6 + 2.2 * sqrt(s.e / std::max(1.0, run)));
        mk->Draw();
      }
    auto key = new TMarker(0, 0, isGamma ? 24 : 20);
    key->SetMarkerColor(isGamma ? kRed + 1 : kAzure + 2);
    lg2->AddEntry(key, isGamma ? "gamma" : "electron", "p");
  }
  lg2->Draw();

  std::string png = rtOut(std::string(TString(fn).ReplaceAll(".root", "").Data()) + "_hits.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s\n", png.c_str());
}
