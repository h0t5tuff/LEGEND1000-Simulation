// Analysis: what fraction of decays at each height along the RT wall reaches an HPGe?

#include "rt_geom.h"
#include <set>

//-------------------------------------------------------------------------------
//  1. Spot a germanium tree:
static bool rt_is_germanium_tree(const TString &n) // NtupleUseVolumeName true -> "V0101", false -> "det101" with uid < 5000
{
  if (n.BeginsWith("det"))
  {
    int uid = TString(n(3, n.Length() - 3)).Atoi(); // n(3,len) slices off the "det" prefix
    return uid > 0 && uid < 5000;                   // HPGe < 5000, SiPM >= 5000, LAr >= 10000
  }
  if (n.Length() == 5 && n[0] == 'V')
  {
    for (int i = 1; i < 5; i++)
      if (!isdigit(n[i]))
        return false;
    return true;
  }
  return false; // pen_V####, S####, undergroundlar, ...
}

//-------------------------------------------------------------------------------
//  2. Which events left energy in germanium?
void sim_reach(const char *fn, int nbins = 16, const char *gdml = "")
{
  RT rt = rtLoad(gdml);
  if (!rt.ok)
    return;
  TFile in(fn);
  if (in.IsZombie())
  {
    printf("ERROR: cannot open %s\n", fn);
    return;
  }
  auto d = (TDirectory *)in.Get("stp");
  if (!d)
  {
    printf("ERROR: no stp/ directory in %s\n", fn);
    return;
  }
  std::set<int> hit; // one decay writes many step rows across many detectors, a set counts each event once
  int ndet = 0;
  TIter it(d->GetListOfKeys());
  while (TKey *k = (TKey *)it())
  {
    TString nm = k->GetName();
    if (!rt_is_germanium_tree(nm))
      continue;
    auto t = (TTree *)d->Get(nm);
    if (!t || !t->GetBranch("evtid"))
      continue;
    ndet++;
    int evtid;
    t->SetBranchAddress("evtid", &evtid); // point the branch at our variable, then GetEntry fills it
    for (Long64_t j = 0; j < t->GetEntries(); j++)
    {
      t->GetEntry(j);
      hit.insert(evtid);
    }
  }

  //-------------------------------------------------------------------------------
  //  3. Bin the vertices by height:
  auto vtx = (TTree *)d->Get("vtx");
  if (!vtx)
  {
    printf("ERROR: no stp/vtx tree - the vertex output scheme was off\n");
    return;
  }
  const Long64_t nv = vtx->GetEntries();
  TTreeFormula fx("fx", "xloc_in_m", vtx), fy("fy", "yloc_in_m", vtx),
      fz("fz", "zloc_in_m", vtx), fe("fe", "evtid", vtx);

  std::vector<long> nFired(nbins, 0), nHit(nbins, 0);
  long mEF = 0, mOF = 0, mSS = 0, mBad = 0, hEF = 0, hOF = 0, hSS = 0;
  for (Long64_t i = 0; i < nv; i++)
  {
    vtx->GetEntry(i);
    double x = fx.EvalInstance() * 1000.0, y = fy.EvalInstance() * 1000.0; // m -> mm
    double z = fz.EvalInstance() * 1000.0;
    double r = sqrt(x * x + y * y);
    bool got = hit.count((int)fe.EvalInstance()) > 0; // did this event reach a germanium?

    int b = (int)((z - rt.zBottom) / (rt.zTop - rt.zBottom) * nbins);
    b = std::min(nbins - 1, std::max(0, b));
    nFired[b]++;
    if (got)
      nHit[b]++;
    const char *m = rt.materialAt(r, z);
    if (!strcmp(m, "EFCu"))
    {
      mEF++;
      if (got)
        hEF++;
    }
    else if (!strcmp(m, "OFHC"))
    {
      mOF++;
      if (got)
        hOF++;
    }
    else if (!strcmp(m, "SS"))
    {
      mSS++;
      if (got)
        hSS++;
    }
    else
      mBad++;
  }
  printf("file      : %s\n", fn);
  printf("geometry  : %s\n", rt.file.c_str());
  printf("germanium : %d detector tables\n", ndet);
  printf("events    : %lld fired, %lu reached an HPGe (%.4f %%)\n\n",
         nv, (unsigned long)hit.size(), 100.0 * hit.size() / nv);
  if (mBad)
    printf("WARNING: %ld vertices outside the RT wall - confinement bug\n\n", mBad);

  //-------------------------------------------------------------------------------
  //  4. The radiopurity weighting table:
  printf("[a] by section (uniform-by-volume sampling, so fired counts ~ volume):\n");
  printf("      %-6s %10s %10s %12s\n", "sect", "fired", "reached", "reach [%]");
  struct
  {
    const char *n;
    long f, h;
  } sec[] = {{"EFCu", mEF, hEF}, {"OFHC", mOF, hOF}, {"SS", mSS, hSS}};
  for (auto &q : sec)
    printf("      %-6s %10ld %10ld %12.4f\n", q.n, q.f, q.h, q.f ? 100.0 * q.h / q.f : 0.0);
  printf("\n      Multiply each row by that section's specific activity and mass\n"
         "      to get the rate; do NOT re-split by z (the EFCu top lid sits\n"
         "      above the SS seam, so a z-cut would misassign it).\n");

  //-------------------------------------------------------------------------------
  //  5. Reach vs height:
  auto g = new TGraphErrors();
  printf("\n[b] reach vs height:\n");
  printf("      %10s %10s %10s %10s %12s\n", "z [mm]", "sect", "fired", "reached", "reach [%]");
  for (int b = 0; b < nbins; b++)
  {
    double z0 = rt.zBottom + (rt.zTop - rt.zBottom) * b / nbins;
    double z1 = rt.zBottom + (rt.zTop - rt.zBottom) * (b + 1) / nbins;
    double zc = 0.5 * (z0 + z1);
    if (nFired[b] == 0)
    {
      printf("      %10.0f %10s %10d   (empty)\n", zc, rt.sectionAt(zc), 0);
      continue;
    }
    double p = (double)nHit[b] / nFired[b];
    double ep = sqrt(p * (1 - p) / nFired[b]); // binomial error on a fraction
    int n = g->GetN();
    g->SetPoint(n, zc, 100.0 * (p > 0 ? p : 1e-6)); // zeros pinned to a floor so the log axis can draw them
    g->SetPointError(n, 0.5 * (z1 - z0), 100.0 * ep);
    printf("      %10.0f %10s %10ld %10ld %12.4f\n", zc, rt.sectionAt(zc), nFired[b], nHit[b], 100.0 * p);
  }
  auto c = new TCanvas("c_reach", "", 860, 560);
  c->SetLogy();
  c->SetGrid();
  g->SetTitle("RT wall #rightarrow HPGe;decay height z along the RT wall [mm];events reaching an HPGe [%]");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.1);
  g->SetLineColor(kAzure + 2);
  g->SetMarkerColor(kAzure + 2);
  g->Draw("AP");
  c->Update(); // Update first so GetUymin/GetUymax know the drawn axis range
  for (double zs : {rt.seamOFHC, rt.seamSS})
  {
    auto l = new TLine(zs, c->GetUymin(), zs, c->GetUymax()); // mark the section seams
    l->SetLineStyle(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
  }
  TString png = TString(fn).ReplaceAll(".root", "") + "_reach.png";
  c->SaveAs(png);
  printf("\nwrote %s\n", png.Data());
}
