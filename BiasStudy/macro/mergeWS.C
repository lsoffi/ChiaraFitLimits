void mergeWS(TString cat="cat6",TString suffix="cs_btagloose")
{
  TFile *f=TFile::Open(Form("HGG_%s_%s.root",cat.Data(),suffix.Data()),"UPDATE");
  RooWorkspace* w=(RooWorkspace*) f->Get("wspace");
  w->Print();
  TFile *f1=TFile::Open(Form("HGG_%s.root",cat.Data()),"UPDATE");
  RooWorkspace* w1=(RooWorkspace*) f1->Get("wspace");
  RooDataSet* d=(RooDataSet*) w1->data("data_mass");
  d->Print();
  d->SetName("data_mass_nocs");
  f->cd();
  w->import(*d);
  w->Write();
  f->Write();
  f->Close();
}
