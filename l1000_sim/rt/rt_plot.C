// reach (% of the fired gammas that hit any Ge) vs source height z.
#include <set>
void rt_plot()
{
  int z[] = {-1000, 0, 500, 1000, 1500, 2000, 3000, 4900}; // same heights as rt_scan.mac
  int nz = sizeof(z) / sizeof(int);
  auto g = new TGraph();

  for (int i = 0; i < nz; i++)
  {
    TFile f(Form("rt_z%d.root", z[i]));
    auto d = (TDirectory *)f.Get("stp");
    std::set<int> ev; // ids of gammas that hit Ge — a set stores each id ONCE (one gamma writes many rows)
    if (d)
    {
      TIter it(d->GetListOfKeys());
      TKey *k;
      while ((k = (TKey *)it()))
      { // loop over every detector tree
        TString n = k->GetName();
        if (!n.BeginsWith("det"))
          continue;
        TString s = n;
        s.Remove(0, 3); // "det101" -> "101"
        int uid = s.Atoi();
        if (uid > 0 && uid < 5000) // HPGe uid < 5000; Sipm uid >=5000, LAr uid>=10000
        {                          // germanium only
          auto t = (TTree *)d->Get(n);
          if (!t)
            continue;
          // each row = one energy-deposit step, tagged with the id of the gamma that made it
          int evtid;
          t->SetBranchAddress("evtid", &evtid);
          for (long j = 0; j < t->GetEntries(); j++)
          {
            t->GetEntry(j);
            ev.insert(evtid);
          } // note that gamma's id (repeats ignored)
        }
      }
    }
    double reach = 100.0 * ev.size() / 100000.0;    // fraction = (distinct gammas that reached Ge) / (100k fired), in %
    g->SetPoint(i, z[i], reach > 0 ? reach : 1e-3); // zeros at the 1-count floor (log axis)
    printf("z=%6d mm   reach=%.3f %%  (%lu of 100000)\n", z[i], reach, (unsigned long)ev.size());
  }

  auto c = new TCanvas("c", "", 720, 500);
  c->SetLogy();
  c->SetGrid();
  g->SetTitle("RT-wall 2614 keV #rightarrow HPGe;source z along RT-wall [mm];fraction reaching HPGe [%]");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.3);
  g->SetLineColor(kRed + 1);
  g->SetMarkerColor(kRed + 1);
  g->Draw("ALP");
  c->SaveAs("rt_reach.png");
  printf("wrote rt_reach.png\n");
}
