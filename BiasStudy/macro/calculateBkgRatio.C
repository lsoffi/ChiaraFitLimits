{
  TFile *f=TFile::Open("HGG_cat7_cs.root");
  RooWorkspace* w=(RooWorkspace*) f->Get("wspace");
  RooExponential* fitFcn=(RooExponential*) w->pdf("pow");
  RooRealVar* mass=(RooRealVar*) w->var("CMS_hgg_mass");
  mass->setRange("sigWindow",130.25-1.5.,130.25+1.5);
  //  mass->setRange("sigWindow",100,110);
  RooAbsReal *intRange = fitFcn->createIntegral(*mass,RooFit::NormSet(*mass),RooFit::Range("sigWindow"));
  std::cout << intRange->getVal() <<std::endl; 
}
