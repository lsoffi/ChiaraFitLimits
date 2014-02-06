{
  bool blind=true;
  bool readFromFile=false;
  bool weightedDataset=false;
  bool doSignalHisto=true;
  //  TString suffix="_cs_btagloose";
  TString suffix="";

  TFile *f=TFile::Open("CMS-HGG_ttH_right_2_2.root");
  //TFile *f=TFile::Open("CMS-HGG_drJet.root");

//   TString files[7]=
//     {
//       "events_55.0_25.0_-10000.0_-10000.0_-100000.0_-10000.0_0.0_0.0_-10000.0_0.0_-1_-1_0_4_cs",
//       "events_55.0_25.0_-10000.0_-10000.0_30.0_-100000.0_0.0_0.0_-10000.0_0.0_-1_-1_0_4_cs",
//       "events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_0.0_0.0_-10000.0_0.0_-1_-1_0_4_cs",
//       "events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_0.0_0.0_350.0_0.0_-1_-1_0_4_cs",
//       "events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_0.0_350.0_0.0_-1_-1_0_4_cs",
//       "events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_0.0_-1_-1_0_4_cs",
//       "events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_2.6_-1_-1_0_4_cs",
//     }

   TString files[2]=
     {
       "ttH_lepCS_2jet_btagloose",
       "ttH_hadCS_2jet_btagloose"
     }


  gSystem->Load("libRooFit");
  gSystem->Load("libRooFitCore");

  gROOT->SetStyle("Plain");
  gStyle->SetOptTitle(0);

  TString names[8]=
    {
      "exp_Norm[CMS_hgg_mass]",
      "2exp_Norm[CMS_hgg_mass]",
      "3exp_Norm[CMS_hgg_mass]",
      "pow_Norm[CMS_hgg_mass]",
      "2pow_Norm[CMS_hgg_mass]",
      "2pol_Norm[CMS_hgg_mass]",
      "3pol_Norm[CMS_hgg_mass]",
      "4pol_Norm[CMS_hgg_mass]"
    }

  TString funName[8]=
    {
      "1exp",
      "2exp",
      "3exp",
      "1pow",
      "2pow",
      "2pol",
      "3pol",
      "4pol"
    }

  float chi2[10];
  int ndof[8]=
    {
      2,
      4,
      6,
      2,
      4,
      2,
      3,
      4
    };

  int masses[9]=
    {
      110,
      115,
      120,
      125,
      130,
      135,
      140,
      145,
      150
    };

  RooWorkspace* wdata=(RooWorkspace*)f->Get("cms_hgg_workspace");
  for (int ifile=6; ifile<8; ifile++)
    {
      
      //      RooRealVar *CMS_hgg_mass = new RooRealVar("CMS_hgg_mass","#gamma#gamma invmass",100,180,"GeV/c^{2}");
      RooRealVar *CMS_hgg_mass = (RooRealVar*) wdata->var("CMS_hgg_mass");
//       RooRealVar *ptgg = new RooRealVar("ptgg","#gamma#gamma pt",0.,1000,"GeV/c");
//       RooRealVar *ebeb = new RooRealVar("ebeb","ebeb",0.,10.);
      RooRealVar *weight = new RooRealVar("weight","weight",0.,100.);
      RooFormulaVar wF("weightF","event weight","@0*1.",*weight);
      
       RooArgSet totalSet;
       totalSet.add(*CMS_hgg_mass);
//       totalSet.add(*ptgg);
//       totalSet.add(*ebeb);
       totalSet.add(*weight);
      
      using namespace RooFit;
      RooDataSet *datacs;
      if (readFromFile)
	{
	  datacs = new RooDataSet("datacs","datacs",totalSet);
	  datacs=RooDataSet::read(files[ifile-6]+".txt",totalSet);
	}
      else
	{
	  datacs = (RooDataSet*)wdata->data(Form("data_mass_cat%d",ifile));
	}
      datacs->SetName("data_mass");
      //
      RooDataSet* wdatacs=0;
      if (weightedDataset)
	{
	  RooRealVar* weightF=(RooRealVar*)datacs->addColumn(wF);
	  wdatacs=new RooDataSet(datacs->GetName(),datacs->GetTitle(),datacs,*datacs->get(),0,weightF->GetName());
	}
      else
	{
	  wdatacs=datacs;
	}
      wdatacs->Print();
      
      
      //  RooDataSet *wdatacs = new RooDataSet("wdatacs","wdatacs",diPhotonEvents,totalSet);
      
      //   if (!wdatacs->isWeighted())
      //     {
      //       std::cout << "ERROR DATASET NOT WEIGHTED" << std::endl;
      //       exit(-1);
      //     }
      
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
      RooExponential exp1("exp","bkg model exp1", *CMS_hgg_mass, tau1 );
      RooExponential exp2_1("exp2_1","bkg model exp2", *CMS_hgg_mass, tau2_1 );
      RooExponential exp2_2("exp2_2","bkg model exp2", *CMS_hgg_mass, tau2_2 );
      RooExponential exp3_1("exp3_1","bkg model exp3", *CMS_hgg_mass, tau3_1 );
      RooExponential exp3_2("exp3_2","bkg model exp3", *CMS_hgg_mass, tau3_2 );
      RooExponential exp3_3("exp3_3","bkg model exp3", *CMS_hgg_mass, tau3_3 );
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
      RooGenericPdf pow1("pow","bkg model pow1", "pow(@0,-@1)",RooArgList(*CMS_hgg_mass,alpha1));
      RooGenericPdf pow2_1("pow2_1","bkg model pow2", "pow(@0,-@1)",RooArgList(*CMS_hgg_mass,alpha2_1));
      RooGenericPdf pow2_2("pow2_2","bkg model pow2", "pow(@0,-@1)",RooArgList(*CMS_hgg_mass,alpha2_2));
      RooRealVar Npow1("Npow1","background yield",1.e3,0.,1.e8);
      RooRealVar Npow2_1("Npow2_1","background yield",1.e3,0.,1.e8);
      RooRealVar Npow2_2("Npow2_2","background yield",1.e3,0.,1.e8);
      RooAddPdf powsum("2pow","2pow",RooArgList(pow2_1,pow2_2),RooArgList(Npow2_1,Npow2_2));
      
      RooBernstein bkg1("1pol","bkg model", *CMS_hgg_mass, cbpars1 ); 
      RooRealVar Nbkg1("1Nbkg","background yield 2",1.e3,0.,1.e8);
      RooBernstein bkg2("2pol","bkg model", *CMS_hgg_mass, cbpars2 ); 
      RooRealVar Nbkg2("2Nbkg","background yield 2",1.e3,0.,1.e8);
      RooBernstein bkg3("3pol","bkg model", *CMS_hgg_mass, cbpars3 ); 
      RooRealVar Nbkg3("3Nbkg","background yield 3",1.e3,0.,1.e8);
      RooBernstein bkg4("4pol","bkg model", *CMS_hgg_mass, cbpars4 ); 
      RooRealVar Nbkg4("4Nbkg","background yield 3",1.e3,0.,1.e8);
      
      RooNLLVar nllexp1("nllexp1","log likelihood var",exp1,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nllexp2("nllexp2","log likelihood var",expsum,*wdatacs,RooFit::Extended(true),RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nllexp3("nllexp3","log likelihood var",exp2sum,*wdatacs,RooFit::Extended(true),RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nllpow1("nllpow1","log likelihood var",pow1,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nllpow2("nllpow2","log likelihood var",powsum,*wdatacs,RooFit::Extended(true),RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nll3("nll3","log likelihood var",bkg3,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nll2("nll2","log likelihood var",bkg2,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nll1("nll1","log likelihood var",bkg1,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
      RooNLLVar nll4("nll4","log likelihood var",bkg4,*wdatacs,RooFit::SumW2Error(kTRUE),RooFit::Save());
      
      RooMinuit mexp1(nllexp1);
      mexp1.migrad();
      RooFitResult* rexp1=mexp1.save();
      RooAbsPdf* pexp1=rexp1->createHessePdf(RooArgSet(tau1));
      RooMinuit mexp2(nllexp2);
      mexp2.migrad();
      RooFitResult* rexp2=mexp2.save();
      RooAbsPdf* pexp2=rexp2->createHessePdf(RooArgSet(tau2_1,tau2_1));
      RooMinuit mexp3(nllexp3);
      mexp3.migrad();
      RooFitResult* rexp3=mexp3.save();
      RooAbsPdf* pexp3=rexp3->createHessePdf(RooArgSet(tau3_1,tau3_2,tau3_3));
      RooMinuit mpow1(nllpow1);
      mpow1.migrad();
      RooFitResult* rpow1=mpow1.save();
      RooAbsPdf* ppow1=rpow1->createHessePdf(RooArgSet(alpha1));
      RooMinuit mpow2(nllpow2);
      mpow2.migrad();
      RooMinuit m(nll3);
      m.migrad();
      RooFitResult *rbkg3=m.save();
      RooAbsPdf* pbkg3=rbkg3->createHessePdf(cbpars3);
      RooMinuit m2(nll2);
      m2.migrad();
      RooFitResult *rbkg2=m2.save();
      RooAbsPdf* pbkg2=rbkg2->createHessePdf(cbpars2);
      RooMinuit m1(nll1);
      m1.migrad();
      RooFitResult *rbkg1=m1.save();
      RooAbsPdf* pbkg1=rbkg1->createHessePdf(cbpars1);
      RooMinuit m4(nll4);
      m4.migrad();
      RooFitResult *rbkg4=m4.save();
      RooAbsPdf* pbkg4=rbkg4->createHessePdf(cbpars4);


      RooPlot *p=CMS_hgg_mass->frame(100,180,10);
      wdatacs->plotOn(p,RooFit::Invisible());

      //wdatacs->plotOn(p);
      //  bkg3.Print("t");
      //  model.plotOn(xframe, RooFit::Components(bkg3),RooFit::LineStyle(kDashed)) ;
      exp1.plotOn(p,RooFit::LineColor(kRed));
      std::cout << exp1.GetName() << "," << p->chiSquare(1) << std::endl;
      chi2[0]=p->chiSquare(1);
      expsum.plotOn(p,RooFit::LineColor(kBlue));
      std::cout << expsum.GetName() << "," << p->chiSquare(2) << std::endl;
      chi2[1]=p->chiSquare(2);
      exp2sum.plotOn(p,RooFit::LineColor(8));
      std::cout << exp2sum.GetName() << "," << p->chiSquare(3) << std::endl;
      chi2[2]=p->chiSquare(3);
      pow1.plotOn(p,RooFit::LineColor(kRed),RooFit::LineStyle(kDashed),RooFit::LineWidth(4));
      std::cout << pow1.GetName() << "," << p->chiSquare(1) << std::endl;
      chi2[3]=p->chiSquare(1);
      powsum.plotOn(p,RooFit::LineColor(kBlue),RooFit::LineStyle(kDashed),RooFit::LineWidth(4));
      std::cout << powsum.GetName() << "," << p->chiSquare(2) << std::endl;
      chi2[4]=p->chiSquare(2);
      //   bkg1.plotOn(p,RooFit::LineColor(kYellow));
      //   std::cout << p->chiSquare(1) << std::endl;
      bkg2.plotOn(p,RooFit::LineColor(kGreen),RooFit::LineStyle(9));
      std::cout << bkg2.GetName() << "," << p->chiSquare(2) << std::endl;
      chi2[5]=p->chiSquare(2);
      bkg3.plotOn(p,RooFit::LineColor(kYellow),RooFit::LineStyle(9));
      std::cout << bkg3.GetName() << "," << p->chiSquare(3) << std::endl;
      chi2[6]=p->chiSquare(3);
      bkg4.plotOn(p,RooFit::LineColor(kBlack),RooFit::LineStyle(9));
      std::cout << bkg4.GetName() << "," << p->chiSquare(4) << std::endl;
      chi2[7]=p->chiSquare(4);
      
       if (!blind)
 	wdatacs->plotOn(p,RooFit::DataError(RooAbsData::SumW2));
       else
 	{
 	  CMS_hgg_mass->setRange("unblind_up",140,180);
 	  wdatacs->plotOn(p,RooFit::CutRange("unblind_up"),RooFit::DataError(RooAbsData::SumW2));
	  CMS_hgg_mass->setRange("unblind_down",100,110);
	  wdatacs->plotOn(p,RooFit::CutRange("unblind_down"),RooFit::DataError(RooAbsData::SumW2));
 	}

       p->SetMinimum(0.0001);
       p->SetMaximum(p->GetMaximum()*1.4);
      p->GetXaxis()->SetTitle("#gamma #gamma mass (GeV/c^{2})");
      p->SetMarkerSize(1.3);



       for (int i=0; i<p->numItems(); i++) {
 	TString obj_name=p->nameOf(i); if (obj_name=="") continue;
 	cout << Form("%d. '%s'\n",i,obj_name.Data());
       }

       p->Draw();
//      RooWorkspace w(files[ifile].Data(),"DataCS fits");
      RooWorkspace w("wspace","DataCS fits");
      w.import(*wdatacs);

      w.import(exp1);
      w.import(*rexp1);
      w.import(*pexp1);

      w.import(expsum);
      w.import(exp2sum);

      w.import(pow1);
      w.import(*rpow1);
      w.import(*ppow1);

      w.import(powsum);

      w.import(bkg2);
      w.import(*pbkg2);
      w.import(*rbkg2);

      w.import(*rbkg1);
      w.import(*pbkg1);
      w.import(bkg1);

      w.import(*rbkg3);
      w.import(*pbkg3);
      w.import(bkg3);

      w.import(*pbkg4);
      w.import(bkg4);
      w.import(*rbkg4);

      w.writeToFile(Form("HGG_cat%d%s.root",ifile,suffix.Data()));
      w.Delete();

      TLegend* leg=new TLegend(0.54,0.6,0.84,0.89);
      leg->SetFillColor(0);
      leg->SetTextSize(0.03);
      leg->SetBorderSize(0);
      leg->AddEntry(p->getObject(0),"Data CS","Pl");
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

//       float yhi = 0.90;
//       float xmin= 0.63;
//       float ypass= 0.03;
//       char tmpStr[300];
//       char line[300];
      
       TLatex lat;
       char line[300];
       lat.SetNDC();
       lat.DrawLatex(0.2,0.94,"CMS Preliminary #sqrt{s}=8 TeV L=19.6 fb^{-1}");

//       for (int iline=0;iline<8;++iline)
// 	{
// 	float prob=TMath::Prob(chi2[iline]*(p->GetNbinsX()-ndof[iline]),(p->GetNbinsX()-ndof[iline]));
// 	std::cout << prob << "," << p->GetNbinsX()-ndof[iline] << std::endl;
// 	sprintf(line,"%s #chi ^{2}: %.3f Prob: %.3f",funName[iline].Data(),chi2[iline],prob);
// 	  lat.DrawLatex(xmin,yhi-ypass, line);
// 	  yhi=yhi-ypass;
// 	}
      c1->SaveAs(Form("HGG_cat%d.png",ifile));


      if (doSignalHisto)
	{
	  TFile *out=TFile::Open(Form("histogg_cat%d.root",ifile),"RECREATE");
	  for (int imass=0;imass<9;++imass)
	    {
	      printf("Saving histos for mass %d in cat%d\n",masses[imass],ifile); 
	      TH1F* ggH=(TH1F*)f->Get(Form("th1f_sig_ggh_mass_m%d_cat%d",masses[imass],ifile));
	      TH1F* vbf=(TH1F*)f->Get(Form("th1f_sig_vbf_mass_m%d_cat%d",masses[imass],ifile));
	      TH1F* vh=(TH1F*)f->Get(Form("th1f_sig_wzh_mass_m%d_cat%d",masses[imass],ifile));
	      TH1F* tth=(TH1F*)f->Get(Form("th1f_sig_tth_mass_m%d_cat%d",masses[imass],ifile));
	      out->cd();
	      ggH->Write(Form("var_gluglu_2011_%d",masses[imass]));
	      vbf->Write(Form("var_vbf_2011_%d",masses[imass]));
	      vh->Write(Form("var_wzh_2011_%d",masses[imass]));
	      tth->Write(Form("var_tth_2011_%d",masses[imass]));
	    }
	  out->Write();
	  out->Close();
	}
    }
}
  
