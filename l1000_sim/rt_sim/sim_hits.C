//  what remage actually gives us: energy deposited in the argon and in the crystals,
//  how a pile of Geant4 steps becomes one detector response, and which gamma caused it
//    root -l -b -q 'sim_hits.C("output/tl208.root")'
//  a gamma entering a crystal Compton-scatters, the electrons it frees deposit in a tight
//  cluster, and the detector reports ONE number. this macro shows both ends of that arrow

#include "rt_geom.h"
#include <vector>

//-------------------------------------------------------------------------------
//  1. Read every step:
struct Step { int ev, det, pdg, tid, par; float e, x, y, z, dsurf, t; };

static const char *pname(int pdg) // the particles that actually deposit in germanium
{
  if (pdg == 22) return "gamma";
  if (pdg == 11) return "e-";
  if (pdg == -11) return "e+";
  if (pdg == 2112) return "neutron";
  return "other";
}

void sim_hits(const char *fn, double m1_keV = 5.0, const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  TFile f(fn);
  auto d = f.IsZombie() ? nullptr : (TDirectory *)f.Get("stp");
  if (!d)
  {
    printf("ERROR: no stp/ in %s\n", fn);
    return;
  }
  std::vector<Step> steps;
  std::vector<TString> detName;
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it())
  {
    TString nm = k->GetName();
    if (nm.Length() != 5 || nm[0] != 'V')
      continue;
    for (int i = 1; i < 5; i++)
      if (!isdigit(nm[i])) { nm = ""; break; }
    if (nm == "")
      continue;
    auto t = (TTree *)d->Get(nm);
    if (!t || !t->GetBranch("edep_in_keV"))
      continue;
    int id = detName.size();
    detName.push_back(nm);
    bool hasTrk = t->GetBranch("trackid");
    TTreeFormula fe("fe", "evtid", t), fp("fp", "particle", t), fd("fd", "edep_in_keV", t),
        fx("fx", "xloc_in_m", t), fy("fy", "yloc_in_m", t), fz("fz", "zloc_in_m", t),
        fs("fs", "dist_to_surf_in_m", t), fm("fm", "time_in_ns", t),
        fi("fi", hasTrk ? "trackid" : "evtid", t), fq("fq", hasTrk ? "parent_trackid" : "evtid", t);
    for (Long64_t j = 0; j < t->GetEntries(); j++)
    {
      t->GetEntry(j);
      steps.push_back({(int)fe.EvalInstance(), id, (int)fp.EvalInstance(),
                       hasTrk ? (int)fi.EvalInstance() : -1, hasTrk ? (int)fq.EvalInstance() : -1,
                       (float)fd.EvalInstance(), (float)(fx.EvalInstance() * 1000),
                       (float)(fy.EvalInstance() * 1000), (float)(fz.EvalInstance() * 1000),
                       (float)(fs.EvalInstance() * 1000), (float)fm.EvalInstance()});
    }
  }
  if (steps.empty())
  {
    printf("no germanium steps in %s\n", fn);
    return;
  }
  // the gamma catalogue: (event, trackid) -> kinetic energy at creation
  std::map<std::pair<int, int>, double> gamma;
  if (auto tr = (TTree *)d->Get("tracks"))
  {
    TTreeFormula e("e", "evtid", tr), i("i", "trackid", tr), p("p", "particle", tr), k("k", "ekin_in_MeV", tr);
    for (Long64_t j = 0; j < tr->GetEntries(); j++)
    {
      tr->GetEntry(j);
      if ((int)p.EvalInstance() == 22)
        gamma[{(int)e.EvalInstance(), (int)i.EvalInstance()}] = k.EvalInstance() * 1000; // -> keV
    }
  }
  // energy the same decays left in the underground argon, per event
  std::map<int, double> larE;
  if (auto tl = (TTree *)d->Get("undergroundlar"))
  {
    TTreeFormula e("e", "evtid", tl), q("q", "edep_in_keV", tl);
    for (Long64_t j = 0; j < tl->GetEntries(); j++)
    {
      tl->GetEntry(j);
      larE[(int)e.EvalInstance()] += q.EvalInstance();
    }
  }
  // a deposit is made by an electron the gamma freed, so one hop up the parentage
  // reaches the gamma; a deposit by the gamma itself is already there
  auto gammaOf = [&](const Step &s) {
    int g = (s.pdg == 22) ? s.tid : s.par;
    auto f = gamma.find({s.ev, g});
    return f == gamma.end() ? -1.0 : f->second;
  };

//-------------------------------------------------------------------------------
//  2. Collapse steps into detector responses:
  // key = event and detector together, so two detectors in one decay stay separate
  struct Resp { double e = 0; int n = 0; double x = 0, y = 0, z = 0; };
  std::map<std::pair<int, int>, Resp> resp;
  double eByPdg[5] = {0}; long nByPdg[5] = {0};
  auto slot = [](int pdg) { return pdg == 22 ? 0 : pdg == 11 ? 1 : pdg == -11 ? 2 : pdg == 2112 ? 3 : 4; };
  for (auto &s : steps)
  {
    auto &r = resp[{s.ev, s.det}];
    r.e += s.e;
    r.n++;
    r.x += s.x * s.e; // energy-weighted centroid: where the response "happened"
    r.y += s.y * s.e;
    r.z += s.z * s.e;
    eByPdg[slot(s.pdg)] += s.e;
    nByPdg[slot(s.pdg)]++;
  }
  std::set<int> evts;
  long nAbove = 0;
  double eTot = 0;
  for (auto &p : resp)
    if (p.second.e > m1_keV) { nAbove++; evts.insert(p.first.first); eTot += p.second.e; }

  printf("file   : %s\n", fn);
  printf("\n[a] energy deposition  ->  detector response\n");
  printf("      %-34s %10ld\n", "Geant4 steps in germanium", (long)steps.size());
  printf("      %-34s %10zu\n", "detector responses (evt x det)", resp.size());
  printf("      %-34s %10ld   (above %.0f keV)\n", "responses kept", nAbove, m1_keV);
  printf("      %-34s %10zu\n", "decays that produced one", evts.size());
  printf("      -> on average %.1f steps collapse into one number per detector\n",
         resp.size() ? (double)steps.size() / resp.size() : 0.0);

  printf("\n[b] which particle actually deposits the energy:\n");
  printf("      %-9s %12s %14s %10s\n", "particle", "steps", "energy [keV]", "share");
  double eAll = 0;
  for (int i = 0; i < 5; i++) eAll += eByPdg[i];
  for (int i = 0; i < 5; i++)
    if (nByPdg[i])
      printf("      %-9s %12ld %14.1f %9.1f%%\n", pname(i == 0 ? 22 : i == 1 ? 11 : i == 2 ? -11 : i == 3 ? 2112 : 0),
             nByPdg[i], eByPdg[i], eAll ? 100 * eByPdg[i] / eAll : 0.0);
  printf("      the gammas carry the energy in, but it is the electrons they free that\n"
         "      deposit it - which is why the response is a tight cluster, not a point\n");

  // walk one real response step by step: the picture in the notebook
  std::pair<int, int> show;
  int bestN = 0;
  for (auto &p : resp)
    if (p.second.e > 100 && p.second.n > bestN) { bestN = p.second.n; show = p.first; }
  printf("\n[c] one detector response, step by step (decay %d in %s):\n", show.first, detName[show.second].Data());
  printf("      %-9s %12s %12s %12s %12s %12s\n", "particle", "edep [keV]", "x [mm]", "y [mm]", "z [mm]", "depth [mm]");
  double run = 0;
  std::vector<Step> shown;
  for (auto &s : steps)
    if (s.ev == show.first && s.det == show.second)
    {
      run += s.e;
      shown.push_back(s);
      if (shown.size() <= 12)
        printf("      %-9s %12.3f %12.3f %12.3f %12.3f %12.4f\n", pname(s.pdg), s.e, s.x, s.y, s.z, s.dsurf);
    }
  if (shown.size() > 12)
    printf("      ... %zu more steps\n", shown.size() - 12);
  printf("      %-9s %12.3f  <- what the detector reports\n", "SUM", run);

  // --- which gamma hit which detector, and when ---
  printf("\n[d] tracing responses back to the gamma that caused them:\n");
  printf("      %-7s %-8s %11s %10s %12s %12s\n", "decay", "detector", "Ge [keV]",
         "t [ns]", "gamma [keV]", "LAr [keV]");
  long resolved = 0, total = 0;
  std::map<std::pair<int, int>, std::map<double, double>> byGamma; // response -> gamma energy -> energy it delivered
  std::map<std::pair<int, int>, double> respTime;                  // earliest deposit in that detector
  for (auto &s : steps)
  {
    auto key = std::make_pair(s.ev, s.det);
    double g = gammaOf(s);
    total++;
    if (g > 0)
    {
      resolved++;
      byGamma[key][g] += s.e;
    }
    if (!respTime.count(key) || s.t < respTime[key])
      respTime[key] = s.t;
  }
  // credit each response to whichever gamma delivered the most of its energy
  std::map<std::pair<int, int>, double> respGamma;
  for (auto &p : byGamma)
  {
    double best = -1, bestE = -1;
    for (auto &g : p.second)
      if (g.second > bestE) { bestE = g.second; best = g.first; }
    respGamma[p.first] = best;
  }
  int listed = 0;
  for (auto &p : resp)
  {
    if (p.second.e <= m1_keV || listed >= 10)
      continue;
    double g = respGamma.count(p.first) ? respGamma[p.first] : -1;
    printf("      %-7d %-8s %11.2f %10.2f %12s %12.2f\n", p.first.first,
           detName[p.first.second].Data(), p.second.e, respTime[p.first],
           g > 0 ? Form("%.1f", g) : "(secondary)",
           larE.count(p.first.first) ? larE[p.first.first] : 0.0);
    listed++;
  }
  printf("      %.0f%% of deposits trace to a named gamma; the rest are later generations\n",
         total ? 100.0 * resolved / total : 0.0);
  double larTot = 0;
  for (auto &p : larE) larTot += p.second;
  printf("      argon saw %.3e keV over %zu decays, germanium %.3e keV over %zu responses\n",
         larTot, larE.size(), eTot, resp.size());

//-------------------------------------------------------------------------------
//  3. Draw it:
  gStyle->SetOptStat(0);
  auto hRZ = new TH2D("hRZ", "where the energy lands;r [mm];z [mm]", 90, 0, 900, 90, -1000, 1500);
  auto hXY = new TH2D("hXY", "same hits from above;x [mm];y [mm]", 90, -900, 900, 90, -900, 900);
  auto hStep = new TH1D("hStep", "", 120, 0, 3000);
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
  c->cd(1)->SetLogz(); gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.15); hRZ->Draw("COLZ");
  if (auto o = (TTree *)d->Get("detector_origins")) // the array itself, for context
  {
    TTreeFormula ox("ox", "xloc_in_m", o), oy("oy", "yloc_in_m", o), oz("oz", "zloc_in_m", o);
    auto g = new TGraph();
    for (Long64_t i = 0; i < o->GetEntries(); i++)
    {
      o->GetEntry(i);
      double x = ox.EvalInstance() * 1000, y = oy.EvalInstance() * 1000;
      g->SetPoint(g->GetN(), sqrt(x * x + y * y), oz.EvalInstance() * 1000);
    }
    g->SetMarkerStyle(1); g->SetMarkerColor(kGray + 1); g->Draw("P SAME");
  }
  c->cd(2)->SetLogz(); gPad->SetRightMargin(0.14); gPad->SetLeftMargin(0.15); hXY->Draw("COLZ");
  c->cd(3)->SetLogy(); gPad->SetGrid();
  hStep->SetTitle("deposition vs response;energy [keV];counts");
  hStep->SetLineColor(kGray + 2); hStep->SetLineWidth(2); hStep->Draw("HIST");
  hResp->SetLineColor(kAzure + 2); hResp->SetLineWidth(2); hResp->Draw("HIST SAME");
  auto lg = new TLegend(0.45, 0.75, 0.88, 0.88); lg->SetTextSize(0.035);
  lg->AddEntry(hStep, "single Geant4 steps", "l");
  lg->AddEntry(hResp, "summed detector response", "l");
  lg->Draw();
  c->cd(4)->SetLogy(); gPad->SetGrid(); hN->SetLineColor(kOrange + 7); hN->SetLineWidth(2); hN->Draw("HIST");
  c->cd(5)->SetLogy(); gPad->SetGrid(); hD->SetLineColor(kGreen + 2); hD->SetLineWidth(2); hD->Draw("HIST");
  c->cd(6); gPad->SetGrid(); // the one response from [c], drawn
  // frame the cluster on its own extent: gammas can scatter tens of mm from where
  // their electrons deposit, so a fixed window either crops them or hides the cluster
  double x0 = 1e9, x1 = -1e9, z0 = 1e9, z1 = -1e9;
  for (auto &s : shown)
  {
    x0 = std::min(x0, (double)s.x); x1 = std::max(x1, (double)s.x);
    z0 = std::min(z0, (double)s.z); z1 = std::max(z1, (double)s.z);
  }
  double pad = 0.12 * std::max(1.0, std::max(x1 - x0, z1 - z0));
  auto fr = gPad->DrawFrame(x0 - pad, z0 - pad, x1 + pad, z1 + pad);
  fr->SetTitle(Form("decay %d in %s: %zu steps -> %.0f keV;x [mm];z [mm]",
                    show.first, detName[show.second].Data(), shown.size(), run));
  for (auto &s : shown) // marker area follows the energy, so the big deposits stand out
  {
    auto m = new TMarker(s.x, s.z, s.pdg == 22 ? 24 : 20);
    m->SetMarkerColor(s.pdg == 22 ? kRed + 1 : kAzure + 2);
    m->SetMarkerSize(0.6 + 2.2 * sqrt(s.e / std::max(1.0, run)));
    m->Draw();
  }
  auto lg2 = new TLegend(0.60, 0.78, 0.90, 0.90); lg2->SetTextSize(0.035);
  auto mg = new TMarker(0, 0, 24); mg->SetMarkerColor(kRed + 1);
  auto me = new TMarker(0, 0, 20); me->SetMarkerColor(kAzure + 2);
  lg2->AddEntry(mg, "gamma", "p"); lg2->AddEntry(me, "electron", "p"); lg2->Draw();

  std::string png = rtOut(std::string(TString(fn).ReplaceAll(".root", "").Data()) + "_hits.png");
  c->SaveAs(png.c_str());
  printf("\nwrote %s\n", png.c_str());
}
