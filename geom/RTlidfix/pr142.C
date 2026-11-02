// PR #142: top of the re-entrance tube before vs after, from each l1000.gdml, with the
// remage vertices sampled uniformly in reentrance_tube_copper (confine.mac).
// root -l -b -q pr142.C
TString g;  // the GDML being read
TGraph *rz(TString name) {  // rzpoints of a genericPolycone
  auto p = new TGraph();
  int i = g.Index("<genericPolycone name=\"" + name + "\""), end = g.Index("</genericPolycone>", i);
  while ((i = g.Index("<rzpoint r=\"", i + 1)) > 0 && i < end)
    p->AddPoint(atof(g.Data() + i + 12), atof(g.Data() + g.Index("z=\"", i) + 3));
  return p;
}
double kg(TGraph *p) {  // copper mass of the rz polygon revolved about the axis
  double v = 0, *r = p->GetX(), *z = p->GetY();
  for (int j = 0, n = p->GetN(); j < n; j++) {
    int k = (j + 1) % n;
    v += (r[j] * r[j] + r[j] * r[k] + r[k] * r[k]) * (z[k] - z[j]);
  }
  return fabs(v) * TMath::Pi() / 3 * 8.96e-6;
}
void pr142() {
  auto c = new TCanvas("c", "", 1200, 500);
  c->Divide(2, 1);
  const char *tag[2] = {"before", "after"};
  for (int i = 0; i < 2; i++) {
    std::ifstream f(Form("%s/l1000.gdml", tag[i]));
    g.ReadFile(f);
    auto ar = rz("liquid_argon_underground"), ss = rz("reentrance_tube_layer_steel_316L_outer_bound");
    // copper = its solid minus its daughters: the argon and the OFHC copper and steel layers
    double m = kg(rz("reentrance_tube_copper")) - kg(ar);
    for (TString l : {"copper_ofhc", "steel_316L"})
      m -= kg(rz("reentrance_tube_layer_" + l + "_outer_bound")) - kg(rz("reentrance_tube_layer_" + l + "_inner_bound"));
    // lid = copper between the top face and the argon's top closure (its last 2 points: wall, axis)
    int n = ar->GetN();
    double rw = ar->GetX()[n - 2], top = TMath::MaxElement(ss->GetN(), ss->GetY());
    auto lid = new TGraph(4, new double[4]{rw, 0, 0, rw},
                          new double[4]{ar->GetY()[n - 2], ar->GetY()[n - 1], top, top});
    auto t = TFile::Open(Form("%s/vtx.root", tag[i]))->Get<TTree>("stp/vtx");
    c->cd(i + 1);
    gPad->DrawFrame(0, top - 10, 1050, top + 1,
        Form("%s: copper %.2f kg, %lld of %lld vertices in the top 10 mm;r [mm];z [mm]",
             tag[i], m, t->GetEntries(Form("zloc_in_m > %g", (top - 10) / 1e3)), t->GetEntries()));
    ss->SetFillColor(kGray + 1);
    ar->SetFillColor(kWhite);
    lid->SetFillColor(kOrange);
    ss->Draw("F");
    ar->Draw("F");
    lid->Draw("F");
    t->SetMarkerStyle(20);
    t->SetMarkerSize(0.5);
    t->Draw("zloc_in_m*1e3 : hypot(xloc_in_m, yloc_in_m)*1e3", "", "same");
    gPad->RedrawAxis();
    auto leg = new TLegend(0.12, 0.15, 0.55, 0.3);
    leg->AddEntry(ss, Form("316L steel wall, %.0f mm", ss->GetX()[1] - rw), "f");
    leg->AddEntry(lid, Form("copper lid, %.2f kg", kg(lid)), "f");
    leg->Draw();
  }
  c->SaveAs("pr142.png");
}
