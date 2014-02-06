
using namespace RooFit;
using namespace RooStats ;


Float_t MINmass= 130;
Float_t MAXmass= 1000;

static const Int_t NCAT = 1;  // chiara

std::string filePOSTfix="";
double signalScaler=1.00;




void SetConstantParams(const RooArgSet* params);



TPaveText* get_labelCMS( int legendQuadrant = 0 , std::string year="2012", bool sim=false) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.83;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendQuadrant==2 ) {
    x1 =  0.25;
    y1 = 0.83;
    x2 =  0.42;
    y2 = 0.87;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
    y2 = 0.24;
  } else if( legendQuadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;
  }

  
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
 
    std::string leftText;
   
     
    if (sim)  leftText = "CMS Simulation"; //cwr ->remove 2011
    else {
     leftText = "CMS Preliminary, 19.5 fb^{-1}";
    }
    cmslabel->AddText(leftText.c_str());
    return cmslabel;

}




TPaveText* get_labelSqrt( int legendQuadrant ) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }


  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendQuadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if( legendQuadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->AddText("#sqrt{s} = 8 TeV");
  return label_sqrt;

}







RooArgSet* defineVariables() {

  // define variables of the input ntuple //livia
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",MINmass, MAXmass,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","Reweightings",0,1000,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  //  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  //  RooRealVar* newweight = new RooRealVar("newweight","Reweightings",0., 1000,"");
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9, *evweight, *nvtx);
  
  return ntplVars;
}






void runfits(const Float_t mass=150, Bool_t dobands = false, Int_t perc = 0) {

  //********************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //*******************************************************************//

  TString fileBaseName(TString::Format("HighMass-hgg", mass));    
  TString fileBkgName(TString::Format("HighMass-hgg.inputbkg_8TeV", mass));
  

  TString card_name("HighMass-hgg_models_Bkg_8TeV_test.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
  RooFitResult* fitresults_bern;
  RooFitResult* fitresults_dijet;
  RooFitResult* fitresults_dijetPL;
  RooFitResult* fitresults_dijetEXP;
  RooFitResult* fitresults_dijetEXPOL;
  RooFitResult* fitresults_expol;
  RooFitResult* fitresults_expolPL;
  RooFitResult* fitresults_rookey;
  RooFitResult* fitresults_lau;
  RooFitResult* fitresults_expPAR;
  
  if(mass<600) {
    MINmass=300;
    MAXmass=1000;
    
  }else if(mass>=600 && mass <750){
    MINmass=200;
    MAXmass=1000;

  }else if(mass>=750 && mass <=850){
    MINmass=280;
    MAXmass=1000;

  }

  

  std::cout<<"  MIN MASS: "<<MINmass<<std::endl;

  Double_t MMIN = 130.;
  Double_t MMAX = 1000; 
  w->var("PhotonsMass")->setMin(MMIN);
  w->var("PhotonsMass")->setMax(MMAX);
  

  w->Print("v");
  

 
  cout << endl; cout << "Now AddBkgData" << endl;
  AddBkgData(w, mass);         
  
 
  cout << endl; cout << "Now BkgModelFit" << endl;    
  bool blind = true; 
  bool dobandsHere= false;

  //fitresults_bern = BkgModelFitBernstein(w, dobands, mass, blind);  
  
  //fitresults_expol = BkgModelFitExpolFunc(w, dobands, mass, blind);  
  //fitresults_lau = BkgModelFitLauFunc(w, dobands, mass, blind);  
  fitresults_expPAR = BkgModelFitExpPARFunc(w, dobands, mass, blind);  
  //  fitresults_expolPL = BkgModelFitExpolPLFunc(w, dobands, mass, blind);  
  // fitresults_dijet = BkgModelFitDiJetFunc(w, dobandsHere, mass, blind);     
 
  return;
}

// Data dataset
void AddBkgData(RooWorkspace* w, Float_t mass) {

  // initializations
  Int_t ncat = NCAT;

  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;
  maxMassFit = MAXmass;

  // retrieve the data tree; no common preselection cut applied yet; 
  TString inDir = "";
  TFile dataFile("histograms_CMS-HGG_24072013.root");
  TTree* dataTree = (TTree*) dataFile.Get("Data");

  // Variables
  RooArgSet* ntplVars = defineVariables();

  // common preselection cut
  TString mainCut = TString::Format("PhotonsMass>=(%.1f) && PhotonsMass<=(%.1f)", minMassFit, maxMassFit);   // livia
  
  // Create dataset
  RooDataSet Data("Data","dataset",dataTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "Data, everything: " << endl;
  Data.Print("v");
  cout << "---- nX:  " << Data.sumEntries() << endl;
  cout << endl;

  // split into NCAT categories
  RooDataSet* dataToFit[NCAT];  
  for (int c=0; c<ncat; ++c) {
    int theCat = c+1;

  // 1)  prime 4 cat livia
    if (c==0) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) dataToFit[c] = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));   


    cout << endl; cout << "for category = " << c << endl;
    dataToFit[c]->Print("v");
    cout << "---- nX:  " << dataToFit[c]->sumEntries() << endl;

    w->import(*dataToFit[c],Rename(TString::Format("Data_cat%d",c)));
  }

  cout << "data, no split" << endl;
  // Create full data set without categorization
  RooDataSet* data = (RooDataSet*) Data.reduce(*w->var("PhotonsMass"),mainCut);
  w->import(*data, Rename("Data"));
  cout << endl;
 
  data->Print("v");
  cout << "---- nX:  " << data->sumEntries() << endl; 
}



RooFitResult* BkgModelFitDiJetFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
  Float_t minMassFit, maxMassFit;
  minMassFit = MINmass;//MINmass;
  maxMassFit = MAXmass;
 
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
   
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJet_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_cat%d",c)));
 
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xDiJet_cat%d",c),"","@0/8000.",*w->var("PhotonsMass"));

   
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_DiJet_truth_cat%d",c), "pow(1-@0, @2)/pow(@0, @1+@3*log(@0))", RooArgList(*x, *p1mod, *p2mod,*p3mod));
   

    fitresult[c] = PhotonsMassBkgTmp0->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkgTmp0);
   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit));//,NormRange("bkg range")); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-3;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >850");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJet","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
 
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678, "brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");

    

    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkgTmp0;
      TGraphAsymmErrors *onesigma = new TGraphAsymmErrors();
      TGraphAsymmErrors *twosigma = new TGraphAsymmErrors();
      
      RooRealVar *nlim = new RooRealVar(TString::Format("nlim%d",c),"",0.0,0.0,10.0);
      nlim->removeRange();
      
      RooCurve *nomcurve = dynamic_cast<RooCurve*>(plotPhotonsMassBkg[c]->getObject(1));
      
      for (int i=1; i<(plotPhotonsMassBkg[c]->GetXaxis()->GetNbins()+1); ++i) {
	double lowedge = plotPhotonsMassBkg[c]->GetXaxis()->GetBinLowEdge(i);
	double upedge  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinUpEdge(i);
	double center  = plotPhotonsMassBkg[c]->GetXaxis()->GetBinCenter(i);
	
	double nombkg = nomcurve->interpolate(center);
	nlim->setVal(nombkg);
	PhotonsMass->setRange("errRange",lowedge,upedge);
	RooAbsPdf *epdf = 0;
	epdf = new RooExtendPdf("epdf","",*cpdf,*nlim,"errRange");
	
	RooAbsReal *nll = epdf->createNLL(*(data[c]),Extended());
	RooMinimizer minim(*nll);
	minim.setStrategy(0);
	double clone = 1.0 - 2.0*RooStats::SignificanceToPValue(1.0);
	double cltwo = 1.0 - 2.0*RooStats::SignificanceToPValue(2.0);
	
	minim.migrad();
	minim.minos(*nlim);
	// printf("errlo = %5f, errhi = %5f\n",nlim->getErrorLo(),nlim->getErrorHi());
	
	onesigma->SetPoint(i-1,center,nombkg);
	onesigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	minim.setErrorLevel(0.5*pow(ROOT::Math::normal_quantile(1-0.5*(1-cltwo),1.0), 2)); // the 0.5 is because qmu is -2*NLL
	// eventually if cl = 0.95 this is the usual 1.92!      
	
	minim.migrad();
	minim.minos(*nlim);
	
	twosigma->SetPoint(i-1,center,nombkg);
	twosigma->SetPointError(i-1,0.,0.,-nlim->getErrorLo(),nlim->getErrorHi());
	
	delete nll;
	delete epdf;
      }

      PhotonsMass->setRange("errRange",minMassFit,maxMassFit);
      
      twosigma->SetLineColor(kGreen);
      twosigma->SetFillColor(kGreen);
      twosigma->SetMarkerColor(kGreen);
      twosigma->Draw("L3 SAME");
      
      onesigma->SetLineColor(kYellow);
      onesigma->SetFillColor(kYellow);
      onesigma->SetMarkerColor(kYellow);
      onesigma->Draw("L3 SAME");
      
      plotPhotonsMassBkg[c]->Draw("SAME"); 
    }

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJET_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJET_LOG_M%d.root",c,massI));

  }



  return fitresult;
}





RooFitResult* BkgModelFitLauFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;

  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult_o2[NCAT];
  RooFitResult* fitresult_o4[NCAT];

  RooPlot* plotPhotonsMassBkg[NCAT];

  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
  
    // fit con expol o2  
    RooFormulaVar *p1mod_o2 = new RooFormulaVar(TString::Format("par1Lau_o2_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau1_o2_cat%d",c)));
    RooAbsPdf* PhotonsMassBkg_o2 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Lau_o2_cat%d",c), "(1-@1)*pow(@0,-4.0)+@1*pow(@0,-5.0)", RooArgList(*PhotonsMass, *p1mod_o2));
    
    fitresult_o2[c] = PhotonsMassBkg_o2->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d o2***********************************", mass, c)<<std::endl;
    fitresult_o2[c]->Print("V");


    // fit con expol o4  
    RooFormulaVar *p1mod_o4 = new RooFormulaVar(TString::Format("par1Lau_o4_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau1_o4_cat%d",c)));
    RooFormulaVar *p2mod_o4 = new RooFormulaVar(TString::Format("par2Lau_o4_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau2_o4_cat%d",c)));   
    RooFormulaVar *p3mod_o4 = new RooFormulaVar(TString::Format("par3Lau_o4_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau3_o4_cat%d",c)));   
    RooAbsPdf* PhotonsMassBkg_o4 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Lau_o4_cat%d",c), "(1-@1-@2-@3)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)", RooArgList(*PhotonsMass, *p1mod_o4, *p2mod_o4, *p3mod_o4));
    
    fitresult_o4[c] = PhotonsMassBkg_o4->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d o4***********************************", mass, c)<<std::endl;
    fitresult_o4[c]->Print("V");

    //Perform FTest:

    double minnll_o2 = fitresult_o2[c]->minNll();
    double minnll_o4 = fitresult_o4[c]->minNll();
   
    double chi2 = 2*(minnll_o2-minnll_o4);
    ndof = 2;
    double prob  = TMath::Prob(chi2,ndof);

    std::cout << " F-test Prob(chi2>chi2(data)) == " << prob << std::endl;


    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkg_o2->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    PhotonsMassBkg_o4->plotOn(plotPhotonsMassBkg[c],LineColor(kRed), LineStyle(kDashed),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
  
    
 

    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >850");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: Laurent o1","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_LAU_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_LAU_M%d.pdf",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_LAU_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_LAU_LOG_M%d.pdf",c,massI));

  }

  RooFitResult* r;

  return r;
}



RooFitResult* BkgModelFitExpPARFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;

  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];


  RooPlot* plotPhotonsMassBkg[NCAT];

  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
  
    // fit con expol 
    RooFormulaVar *p1mod= new RooFormulaVar(TString::Format("par1ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_ExpPAR1_cat%d",c)));
    RooFormulaVar *p2mod= new RooFormulaVar(TString::Format("par2ExpPAR_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_ExpPAR2_cat%d",c)));

    RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_Lau_o2_cat%d",c), "exp(-@1*@0)*pow(@0, @2)", RooArgList(*PhotonsMass, *p1mod, *p2mod));
    
    fitresult[c] = PhotonsMassBkg->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");


    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkg->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
   
  
    
 

    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >850");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: ExpPAR","L");
  
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_M%d.pdf",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPAR_LOG_M%d.pdf",c,massI));

  }

  RooFitResult* r;

  return r;
}


 


RooFitResult* BkgModelFitExpolFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;

  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult_o1[NCAT];
  RooFitResult* fitresult_o2[NCAT];
  RooFitResult* fitresult_o3[NCAT];
  RooFitResult* fitresult_o4[NCAT];
  
  RooPlot* plotPhotonsMassBkg[NCAT];

  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
  
    // fit con expol o1  
    RooFormulaVar *p1mod_o1 = new RooFormulaVar(TString::Format("par1Expol_o1_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_o1_cat%d",c)));
    RooFormulaVar *p2mod_o1 = new RooFormulaVar(TString::Format("par2Expol_o1_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_o1_cat%d",c)));   
    RooFormulaVar *x_o1     = new RooFormulaVar(TString::Format("xExpol_o1_cat%d",c),"","@0",*w->var("PhotonsMass"));
    RooAbsPdf* PhotonsMassBkg_o1 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Expol_o1_cat%d",c), "exp(-@0/(@1+@2*@0))", RooArgList(*x_o1, *p1mod_o1, *p2mod_o1));
    
    fitresult_o1[c] = PhotonsMassBkg_o1->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d o1***********************************", mass, c)<<std::endl;
    fitresult_o1[c]->Print("V");


   // fit con expol o2
    RooFormulaVar *p1mod_o2 = new RooFormulaVar(TString::Format("par1Expol_o2_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_o2_cat%d",c)));
    RooFormulaVar *p2mod_o2 = new RooFormulaVar(TString::Format("par2Expol_o2_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_o2_cat%d",c)));   
    RooFormulaVar *p3mod_o2 = new RooFormulaVar(TString::Format("par3Expol_o2_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol3_o2_cat%d",c)));   
    RooFormulaVar *x_o2     = new RooFormulaVar(TString::Format("xExpol_o2_cat%d",c),"","@0",*w->var("PhotonsMass"));
    RooAbsPdf* PhotonsMassBkg_o2 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Expol_o2_cat%d",c), "exp(-(@0+@3*@0*@0)/(@1+@2*@0))", RooArgList(*x_o2, *p1mod_o2, *p2mod_o2,*p3mod_o2));
    
    fitresult_o2[c] = PhotonsMassBkg_o2->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d o2***********************************", mass, c)<<std::endl;
    fitresult_o2[c]->Print("V");
    /*

   // fit con expol o3
    RooFormulaVar *p1mod_o3 = new RooFormulaVar(TString::Format("par1Expol_o3_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_o3_cat%d",c)));
    RooFormulaVar *p2mod_o3 = new RooFormulaVar(TString::Format("par2Expol_o3_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_o3_cat%d",c)));   
    RooFormulaVar *p3mod_o3 = new RooFormulaVar(TString::Format("par3Expol_o3_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol3_o3_cat%d",c)));   
    RooFormulaVar *p4mod_o3 = new RooFormulaVar(TString::Format("par4Expol_o3_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol4_o3_cat%d",c)));   
    RooFormulaVar *x_o3     = new RooFormulaVar(TString::Format("xExpol_o3_cat%d",c),"","@0",*w->var("PhotonsMass"));
    RooAbsPdf* PhotonsMassBkg_o3 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Expol_o3_cat%d",c), "exp(-(@0+@3*@0*@0+@4*@0*@0*@0)/(@1+@2*@0))", RooArgList(*x_o3, *p1mod_o3, *p2mod_o3, *p3mod_o3, *p4mod_o3));
    
    fitresult_o3[c] = PhotonsMassBkg_o3->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d o3***********************************", mass, c)<<std::endl;
    fitresult_o3[c]->Print("V");
    */
    //Perform FTest:

    double minnll_o1 = fitresult_o1[c]->minNll();
    double minnll_o2 = fitresult_o2[c]->minNll();
    //    double minnll_o3 = fitresult_o3[c]->minNll();

    std::cout<<"o1: "<<fitresult_o1[c]->minNll()<<" o2: "<<fitresult_o2[c]->minNll()<<std::endl;//" o3: "<<minnll_o3<<std::endl;
   
    double chi2 = 2*(minnll_o1-minnll_o2);
    ndof = 1;
    double prob  = TMath::Prob(chi2,ndof);

    std::cout << " F-test Prob(chi2>chi2(data)) == " << prob << std::endl;

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkg_o1->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    PhotonsMassBkg_o2->plotOn(plotPhotonsMassBkg[c],LineColor(kRed),LineStyle(kDashed),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    //  PhotonsMassBkg_o3->plotOn(plotPhotonsMassBkg[c],LineColor(kGreen),LineStyle(kDashed),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
 

    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >850");

      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 


   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw();  

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(3),"Parametric Model: Expol o3","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Parametric Model: Expol o2","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: Expol o1","L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw("same");

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    //write down the chi2 of the fit on the
      
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_EXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_EXPOL_M%d.pdf",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_EXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_FTEST_EXPOL_LOG_M%d.pdf",c,massI));

  }

  RooFitResult* r;

  return r;
}


   


void SetConstantParams(const RooArgSet* params) {

  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}



