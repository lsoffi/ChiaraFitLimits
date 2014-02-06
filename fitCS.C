{
  TString files[4]=
    {
      "results_gg/tieniCS_8TeV/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_70.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4_cs",
      // "results_gg/tieniCS/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_60.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4_cs",
      "results_gg/tieniCS_cut70/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_70.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4_cs",
      "results_gg/tieniData/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_60.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4",
      "results_gg/tieniData_cut70/events_45.0_25.0_-10000.0_-10000.0_-10000.0_-10000.0_70.0_0.0_0.0_-10000.0_0.0_0.0_0.0_0.0_0.0_300.0_0.0_-1_0_0_1_4",
    }

  gSystem->Load("libRooFit");
  gSystem->Load("libRooFitCore");

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);

  TString names[8]={
    "exp_Norm[massggnewvtx]",
    "2exp_Norm[massggnewvtx]",
    "3exp_Norm[massggnewvtx]",
    "pow_Norm[massggnewvtx]",
    "2pow_Norm[massggnewvtx]",
    "1pol_Norm[massggnewvtx]",
    "2pol_Norm[massggnewvtx]",
    "3pol_Norm[massggnewvtx]"
  }

  TString funName[8]={
    "1exp",
    "2exp",
    "3exp",
    "1pow",
    "2pow",
    "1pol",
    "2pol",
    "3pol"
  }

  int ndof[8]={
    2,
    4,
    6,
    2,
    4,
    2,
    3,
    4
  };

  for (int ifile=0; ifile<1; ifile++) {
      
    cout << "fitting sample " << files[ifile] << endl;

    RooRealVar *massggnewvtx = new RooRealVar("massggnewvtx","#gamma#gamma invmass",100,180,"GeV/c^{2}");
    RooRealVar *ptgg = new RooRealVar("ptgg","#gamma#gamma pt",0.,1000,"GeV/c");
    RooRealVar *ebeb = new RooRealVar("ebeb","ebeb",0.,10.);
    RooRealVar *weight = new RooRealVar("weight","weight",0.,100.);
    RooFormulaVar wF("weightF","event weight","@0*1.",*weight);
    
    RooArgSet totalSet;
    totalSet.add(*massggnewvtx);
    totalSet.add(*ptgg);
    totalSet.add(*ebeb);
    totalSet.add(*weight);
    
    using namespace RooFit;
    
    RooDataSet *datacs = new RooDataSet("datacs","datacs",totalSet);
    datacs=RooDataSet::read(files[ifile]+".txt",totalSet);
    RooRealVar* weightF=(RooRealVar*)datacs->addColumn(wF);
    RooDataSet *wdatacs=new RooDataSet(datacs->GetName(),datacs->GetTitle(),datacs,*datacs->get(),0,weightF->GetName());
    wdatacs->Print();
    
    // chiara: mio 
    // TFile *file = new TFile(files[2]+".root"); 
    // TTree *tree = (TTree*)file->Get("diPhotonEvents"); 
    // RooDataSet *wdatacs = new RooDataSet("wdatacs","wdatacs",tree,totalSet);   // weight->GetName()
    // wdatacs->Print(); 

    TFile *roofitFile = TFile::Open("porcacciamiseria.root","recreate");
    roofitFile->cd();
    wdatacs->Write();
    roofitFile->Close();

    RooArgList cbpars3;
    for(int i=0; i<3; ++i) {
      char name[100];
      sprintf(name,"cb%d_3",i);
      RooRealVar* cb3 = new RooRealVar(name,name,0.1,0.,100.);  
      cbpars3.add( *cb3 );
    }
    
    RooArgList cbpars4;
    for(int i=0; i<4; ++i) {
      char name[100];
      sprintf(name,"cb%d_4",i);
      RooRealVar* cb4 = new RooRealVar(name,name,0.1,0.,100.);
      cbpars4.add( *cb4 );
    }
    
    RooArgList cbpars2;
    for(int i=0; i<2; ++i) {
      char name[100];
      sprintf(name,"cb%d_2",i);
      RooRealVar* cb2 = new RooRealVar(name,name,0.1,0.,100.);
      cbpars2.add( *cb2 );
    }
    
    RooArgList cbpars1;
    for(int i=0; i<1; ++i) {
      char name[100];
      sprintf(name,"cb%d_1",i);
      RooRealVar* cb1 = new RooRealVar(name,name,0.1,0.,100.);
      cbpars1.add( *cb1 );
    }
    
    RooRealVar tau1("tau1","tau1",-0.002,-10.,0.);
    RooRealVar tau2_1("tau2_1","tau2_1",-0.002,-10.,0.);
    RooRealVar tau2_2("tau2_2","tau2_2",-0.002,-10.,0.);
    RooRealVar tau3_1("tau3_1","tau3_1",-0.002,-10.,0.);
    RooRealVar tau3_2("tau3_2","tau3_2",-0.002,-10.,0.);
    RooRealVar tau3_3("tau3_3","tau3_3",-0.002,-10.,0.);
    RooExponential exp1("exp","bkg model exp1", *massggnewvtx, tau1 );
    RooExponential exp2_1("exp2_1","bkg model exp2", *massggnewvtx, tau2_1 );
    RooExponential exp2_2("exp2_2","bkg model exp2", *massggnewvtx, tau2_2 );
    RooExponential exp3_1("exp3_1","bkg model exp3", *massggnewvtx, tau3_1 );
    RooExponential exp3_2("exp3_2","bkg model exp3", *massggnewvtx, tau3_2 );
    RooExponential exp3_3("exp3_3","bkg model exp3", *massggnewvtx, tau3_3 );
    RooRealVar Nexp1("Nexp1","background yield",1.e3,0.,1.e8);
    RooRealVar Nexp2_1("Nexp2_1","background yield",1.e3,0.,1.e8);
    RooRealVar Nexp2_2("Nexp2_2","background yield",1.e3,0.,1.e8);
    RooRealVar Nexp3_1("Nexp3_1","background yield",1.e3,0.,1.e8);
    RooRealVar Nexp3_2("Nexp3_2","background yield",1.e3,0.,1.e8);
    RooRealVar Nexp3_3("Nexp3_3","background yield",1.e3,0.,1.e8);
    RooAddPdf expsum("2exp","2exp",RooArgList(exp2_1,exp2_2),RooArgList(Nexp2_1,Nexp2_2));
    RooAddPdf exp2sum("3exp","3exp",RooArgList(exp3_1,exp3_2,exp3_3),RooArgList(Nexp3_1,Nexp3_2,Nexp3_3));
    
    RooRealVar alpha1("alpha1","alpha1",1.,0.,100.);
    RooRealVar alpha2_1("alpha2_1","alpha2_1",1.,0.,100.);
    RooRealVar alpha2_2("alpha2_2","alpha2_2",1.,0.,100.);
    RooGenericPdf pow1("pow","bkg model pow1", "pow(@0,-@1)",RooArgList(*massggnewvtx,alpha1));
    RooGenericPdf pow2_1("pow2_1","bkg model pow2", "pow(@0,-@1)",RooArgList(*massggnewvtx,alpha2_1));
    RooGenericPdf pow2_2("pow2_2","bkg model pow2", "pow(@0,-@1)",RooArgList(*massggnewvtx,alpha2_2));
    RooRealVar Npow1("Npow1","background yield",1.e3,0.,1.e8);
    RooRealVar Npow2_1("Npow2_1","background yield",1.e3,0.,1.e8);
    RooRealVar Npow2_2("Npow2_2","background yield",1.e3,0.,1.e8);
    RooAddPdf powsum("2pow","2pow",RooArgList(pow2_1,pow2_2),RooArgList(Npow2_1,Npow2_2));
    
    RooBernstein bkg1("0pol","bkg model", *massggnewvtx, cbpars1 ); 
    RooRealVar Nbkg1("1Nbkg","background yield 2",1.e3,0.,1.e8);
    RooBernstein bkg2("1pol","bkg model", *massggnewvtx, cbpars2 ); 
    RooRealVar Nbkg2("2Nbkg","background yield 2",1.e3,0.,1.e8);
    RooBernstein bkg3("2pol","bkg model", *massggnewvtx, cbpars3 ); 
    RooRealVar Nbkg3("3Nbkg","background yield 3",1.e3,0.,1.e8);
    RooBernstein bkg4("3pol","bkg model", *massggnewvtx, cbpars4 ); 
    RooRealVar Nbkg4("3Nbkg","background yield 3",1.e3,0.,1.e8);
    
    RooNLLVar nllexp1("nllexp1","log likelihood var",exp1,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nllexp2("nllexp2","log likelihood var",expsum,*wdatacs,RooFit::Extended(true),RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nllexp3("nllexp3","log likelihood var",exp2sum,*wdatacs,RooFit::Extended(true),RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nllpow1("nllpow1","log likelihood var",pow1,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nllpow2("nllpow2","log likelihood var",powsum,*wdatacs,RooFit::Extended(true),RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nll("nll","log likelihood var",bkg3,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nll2("nll2","log likelihood var",bkg2,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nll1("nll1","log likelihood var",bkg1,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
    RooNLLVar nll4("nll4","log likelihood var",bkg4,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
    
    cout << endl;
    cout << "exp1" << endl;
    RooMinuit mexp1(nllexp1);
    mexp1.migrad();
    RooFitResult* rexp1=mexp1.save();
    RooAbsPdf* pexp1=rexp1->createHessePdf(RooArgSet(tau1));
    cout << endl;
    cout << endl;
    cout << "exp2" << endl;
    RooMinuit mexp2(nllexp2);
    mexp2.migrad();
    RooFitResult* rexp2=mexp2.save();
    RooAbsPdf* pexp2=rexp2->createHessePdf(RooArgSet(tau2_1,tau2_1));
    cout << endl;
    cout << endl;
    cout << "exp3" << endl;
    RooMinuit mexp3(nllexp3);
    mexp3.migrad();
    RooFitResult* rexp3=mexp3.save();
    RooAbsPdf* pexp3=rexp3->createHessePdf(RooArgSet(tau3_1,tau3_2,tau3_3));
    cout << endl;
    cout << endl;
    cout << "mpow1" << endl;
    RooMinuit mpow1(nllpow1);
    mpow1.migrad();
    RooFitResult* rpow1=mpow1.save();
    RooAbsPdf* ppow1=rpow1->createHessePdf(RooArgSet(alpha1));
    cout << endl;
    cout << endl;
    cout << "mpow2" << endl;
    RooMinuit mpow2(nllpow2);
    mpow2.migrad();
    cout << endl;
    cout << endl;
    cout << "nll: 2pol" << endl;
    RooMinuit m(nll);
    m.migrad();
    RooFitResult *rbkg=m.save();
    RooAbsPdf* pbkg3=rbkg->createHessePdf(cbpars3);
    cout << endl;
    cout << endl;
    cout << "nll2: 1pol" << endl;
    RooMinuit m2(nll2);
    m2.migrad();
    RooFitResult *rbkg2=m2.save();
    RooAbsPdf* pbkg2=rbkg2->createHessePdf(cbpars2);
    cout << endl;
    cout << endl;
    cout << "nll1: 0pol" << endl;
    RooMinuit m1(nll1);
    m1.migrad();
    RooFitResult *rbkg1=m1.save();
    RooAbsPdf* pbkg1=rbkg1->createHessePdf(cbpars1);
    cout << endl;
    cout << endl;
    cout << "nll4: 3pol" << endl;
    RooMinuit m4(nll4);
    m4.migrad();
    RooFitResult *rbkg4=m4.save();
    RooAbsPdf* pbkg4=rbkg4->createHessePdf(cbpars4);

    
    float chi2[10];
    RooPlot *p=massggnewvtx->frame(100,180,40);
    wdatacs->plotOn(p,RooFit::DataError(RooAbsData::SumW2));
    exp1.plotOn(p,RooFit::LineColor(kRed));
    std::cout << exp1.GetName() << "," << p->chiSquare(1) << std::endl;
    chi2[0]=p->chiSquare(1);
    //
    expsum.plotOn(p,RooFit::LineColor(kBlue));
    std::cout << expsum.GetName() << "," << p->chiSquare(2) << std::endl;
    chi2[1]=p->chiSquare(2);
    //
    exp2sum.plotOn(p,RooFit::LineColor(8));
    std::cout << exp2sum.GetName() << "," << p->chiSquare(3) << std::endl;
    chi2[2]=p->chiSquare(3);
    //
    pow1.plotOn(p,RooFit::LineColor(kRed),RooFit::LineStyle(kDashed),RooFit::LineWidth(4));
    std::cout << pow1.GetName() << "," << p->chiSquare(1) << std::endl;
    chi2[3]=p->chiSquare(1);
    //
    powsum.plotOn(p,RooFit::LineColor(kBlue),RooFit::LineStyle(kDashed),RooFit::LineWidth(4));
    std::cout << powsum.GetName() << "," << p->chiSquare(2) << std::endl;
    chi2[4]=p->chiSquare(2);
    //
    bkg2.plotOn(p,RooFit::LineColor(kGreen),RooFit::LineStyle(9));
    std::cout << bkg2.GetName() << "," << p->chiSquare(2) << std::endl;
    chi2[5]=p->chiSquare(2);
    // 
    bkg3.plotOn(p,RooFit::LineColor(kYellow),RooFit::LineStyle(9));
    std::cout << bkg3.GetName() << "," << p->chiSquare(3) << std::endl;
    chi2[6]=p->chiSquare(3);
    //
    bkg4.plotOn(p,RooFit::LineColor(kBlack),RooFit::LineStyle(9));
    std::cout << bkg4.GetName() << "," << p->chiSquare(4) << std::endl;
    chi2[7]=p->chiSquare(4);
    //
    //

    p->GetXaxis()->SetTitle("#gamma #gamma mass (GeV/c^{2})");
    p->SetMarkerSize(1.3);
    p->Draw();  
    //
    TLegend* leg=new TLegend(0.6,0.6,0.9,0.89);
    leg->SetFillColor(0);
    leg->SetTextSize(0.03);
    leg->SetBorderSize(0);

    for (int ileg=0;ileg<8;++ileg) {
      TObject *obj = p->findObject(names[ileg].Data());
      if (!obj) {
	Warning("fitCS",Form("Can't find item='%s' in the plot!\n",names[ileg].Data()));
	continue;
      }
      float prob=TMath::Prob(chi2[ileg]*(p->GetNbinsX()-ndof[ileg]),(p->GetNbinsX()-ndof[ileg]));
      leg->AddEntry(obj,Form("%s #chi^{2}:%.3f Prob:%.3f",funName[ileg].Data(),chi2[ileg],prob),"l");
    }
    leg->Draw();
    c1->SaveAs(files[ifile]+"_chiaraFitCS.png");
    


    // saving everything in a workspace
    RooWorkspace* wspace = new RooWorkspace("wspace","workspace");

    // import data
    wspace->import(*wdatacs);

    // import model
    wspace->import(exp1);
    wspace->import(*rexp1);
    wspace->import(*pexp1);
    wspace->import(expsum);
    wspace->import(exp2sum);
    wspace->import(pow1);
    wspace->import(*rpow1);
    wspace->import(*ppow1);
    wspace->import(powsum);
    wspace->import(bkg2);
    wspace->import(*pbkg2);
    wspace->import(*rbkg2);
    wspace->import(bkg1);         //
    wspace->import(*rbkg1);
    wspace->import(*pbkg1);
    wspace->import(bkg3);
    wspace->import(*rbkg);
    wspace->import(*pbkg3);
    wspace->import(bkg4);
    wspace->import(*pbkg4);
    wspace->import(*rbkg4);

    // Print workspace contents
    wspace->Print();
    
    // Save the workspace into a root file
    TFile fileWS(files[ifile]+"_ws.root","RECREATE");
    fileWS.cd();
    wspace->Write();
    fileWS.Close();
    cout << "done!" << endl;
  }
}
  
