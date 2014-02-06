
 
 /** \macro H2GGFitter.cc
 *
 * The analysis root trees produced in a simple format 
 *
 *     TFile file(filename,"RECREATE", "X->jj input tree for unbinned maximum-likelihood fit");
 *     TTree* outTree  = new TTree("XTojj","X->jj input tree for unbinned maximum-likelihood fit");
 *     Float_t mass;
 *     Int_t CAT3;
 *     Float_t weight;
 *
 *     outTree->Branch("mass",&mass,"mass/F");
 *     outTree->Branch("weight",&weight,"weight/F");
 *     outTree->Branch("CAT4",&CAT4,"CAT4/I");
 *     {
 *       .............
 *       outTree->Fill();
 *     }
 *
 *     file.Write();
 *     file.Close();
 *     delete outTree;
 *
 * are used as input files. They have to be produced for 
 * data and Monte Carlo signal and background data sets 
 * after all analysis selections to be applied.  *
 */
// Loading:  .L HH2ggbbFitter.cc
// Running:  runfits("hgg120-shapes-combined-Unbinned.root")  
//                

using namespace RooFit;
using namespace RooStats ;

static const Int_t NCAT = 4;  // chiara

std::string filePOSTfix="";
double signalScaler=1.00;

void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*, Float_t);
void MakeRooKeysPDFMCBkg(RooWorkspace*, Float_t);
void SigModelFitGauss(RooWorkspace*, Float_t);
void SigModelFitCBC(RooWorkspace*, Float_t);
RooFitResult*  BkgModelFitBernstein(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitExpo(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitDiJetFunc(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitExpolFunc(RooWorkspace*, Bool_t, Float_t, bool);
void MakePlots(RooWorkspace*, Float_t, RooFitResult*, bool);
void MakeSigWS(RooWorkspace* w, const char* filename);
void MakeBkgWS(RooWorkspace* w, const char* filename);
void SetConstantParams(const RooArgSet* params);
void MakeParameterTrendvsMassGauss();
void MakeParameterTrendvsMassCBC();
void MakePlotMassDataMC(RooWorkspace*, Float_t);
void MakePlotNVTXDataMC(RooWorkspace*, Float_t);



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
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",100, 2000,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* evweight = new RooRealVar("evweight","Reweightings",0,10,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  //  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9,  *evweight, *nvtx);
  
  return ntplVars;
}

void runfits(const Float_t mass=150, Bool_t dobands = false) {

  //******************************************************************//
  //  Running mode  corresponds to the following cases
  //         - full run set:
  //         - create signal and background data sets 
  //         - make and fit signal and background  models 
  //         - write signal and background workspaces in root files
  //         - write data card
  //***********************************ß********************************//

  TString fileBaseName(TString::Format("HighMass-hgg.m%.1f", mass));    
  TString fileBkgName(TString::Format("HighMass-hgg.inputbkg_8TeV", mass));
  
  TString card_name(TString::Format("HighMass-hgg_m%.0f_models_Bkg_8TeV_test.rs", mass));
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();
 
  Double_t MMIN = 100.; //130 //livia   //100  //180   //270
  Double_t MMAX = 800.; //450      //livia //200 // 300  //2000
  w->var("PhotonsMass")->setMin(MMIN);
  w->var("PhotonsMass")->setMax(MMAX);

 
  cout << "Now plot Data vs MC bkg" << endl;
  MakePlotMassDataMC(w, mass);
  //MakePlotNVTXDataMC(w, mass);
  //PlotVtxRecoEff();

  // Make plots for data and fit results
  cout << endl; cout << "Preparing final plots" << endl;
  //MakePlots(w, mass, fitresults_bern, blind);

  //PlotSigShape(w, 0); 
  //PlotSigShape(w, 1); 
  //PlotSigShape(w, 2); 
  //PlotSigShape(w, 3); 

  return;
}

void MakePlots(RooWorkspace* w, Float_t mass, RooFitResult* fitresults, bool blind) {

  Int_t ncat = NCAT;
  
  cout << endl; cout << "Retreive everything:" << endl; 
  w->Print();

  // retrieve data sets from the workspace
  RooDataSet* dataAll   = (RooDataSet*) w->data("Data");
  RooDataSet* signalAll = (RooDataSet*) w->data("SigWeight");
  
  // maximum 9 cat...
  RooDataSet* data[9];  
  RooDataSet* signal[9];
  RooAbsPdf*  PhotonsMassGaussSig[9];
  RooAbsPdf*  PhotonsMassGaussSig_bis[9];
  RooAbsPdf*  PhotonsMassSig[9];
  RooExtendPdf*  PhotonsMassBkg[9];  
  
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
    PhotonsMassGaussSig[c]  = (RooAbsPdf*)  w->pdf(TString::Format("PhotonsMassGaussSig_cat%d",c));
    PhotonsMassGaussSig_bis[c]     = (RooAbsPdf*)  w->pdf(TString::Format("PhotonsMassGaussSig_cat%d_bis",c));
    PhotonsMassSig[c]       = (RooAbsPdf*)  w->pdf("PhotonsMassSig"+TString::Format("_cat%d",c));
    PhotonsMassBkg[c]       = (RooExtendPdf*)  w->pdf(TString::Format("PhotonsMassBkg_cat%d",c));
  }
  
  // retrieve mass observable from the workspace
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
  
  // retrieve pdfs after the fits
  RooAbsPdf* PhotonsMassGaussSigAll = w->pdf("PhotonsMassGaussSig");
  RooAbsPdf* PhotonsMassGaussSigAll_bis    = w->pdf("PhotonsMassGaussSig_bis");
  RooAbsPdf* PhotonsMassSigAll      = w->pdf("PhotonsMassSig");
  RooAbsPdf* PhotonsMassBkgAll      = w->pdf("PhotonsMassBkgAll");

  
 Float_t minMassFit, maxMassFit;
 
 if(mass<200){ //m150
    minMassFit = 100;
    maxMassFit = 200;
  } else if(mass>=200 && mass <500){//m200 and m250
    minMassFit = 180;
    maxMassFit = 500;
  }else{//m300 m400
    minMassFit = 280;
    maxMassFit = 500;
  }


  Float_t MASS(mass);
  int iMass = abs(mass);
  
  Int_t nBinsMass(0.2*mass);  // chiara


  /*  // ****************************SIG ALL
  cout << endl; cout << "Progress plotting: signal" << endl;     

  RooPlot* plotPhotonsMassAll = PhotonsMass->frame(Range(minMassFit,maxMassFit),Bins(nBinsMass));
  signalAll->plotOn(plotPhotonsMassAll);
  
  gStyle->SetOptTitle(0);
  PhotonsMassSigAll->plotOn(plotPhotonsMassAll);
  PhotonsMassSigAll->plotOn(plotPhotonsMassAll,Components("PhotonsMassGaussSig"),LineStyle(kDashed),LineColor(kGreen));
  PhotonsMassSigAll->plotOn(plotPhotonsMassAll,Components("PhotonsMassGaussSig_bis"),LineStyle(kDashed),LineColor(kRed));
  //PhotonsMassSigAll->paramOn(plotPhotonsMassAll, ShowConstants(true), Layout(0.65,0.75,0.85), Format("NEU",AutoPrecision(2)));
  //plotPhotonsMassAll->getAttText()->SetTextSize(0.03);
  
  TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,500,500);
  c1->cd(1);
  plotPhotonsMassAll->Draw();  
  c1->SaveAs("plots/sigmodel_"+TString::Format("%d.png", iMass));
  c1->SaveAs("plots/sigmodel_"+TString::Format("%d.root", iMass));
*/

  // ****************************SIG CAT
  cout << endl; cout << "Progress plotting: signal per categories" << endl;     
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
    
  RooPlot* plotPhotonsMass[9];
  for (int c=0; c<ncat; ++c) {
    plotPhotonsMass[c] = PhotonsMass->frame(Range(mass-30,mass+30),Bins(nBinsMass));
    signal[c]->plotOn(plotPhotonsMass[c],LineColor(kWhite),MarkerColor(kWhite));    
    PhotonsMassSig[c]  ->plotOn(plotPhotonsMass[c]);
    PhotonsMassSig[c]  ->plotOn(plotPhotonsMass[c],Components("PhotonsMassGaussSig"+TString::Format("_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
    PhotonsMassSig[c]  ->plotOn(plotPhotonsMass[c],Components("PhotonsMassGaussSig"+TString::Format("_cat%d_bis",c)),LineStyle(kDashed),LineColor(kRed));
    //PhotonsMassSig[c]  ->paramOn(plotPhotonsMass[c], ShowConstants(true), Layout(0.65,0.75,0.85), Format("NEU",AutoPrecision(2)));
    signal[c]  ->plotOn(plotPhotonsMass[c]);

    //  TCanvas* dummy = new TCanvas("dummy", "dummy",0, 0, 400, 400);
    // TH1F *hist = new TH1F("hist", "hist", 400, minMassFit, maxMassFit);
 
    plotPhotonsMass[c]->SetTitle("");      
    plotPhotonsMass[c]->SetMinimum(0.0);
    plotPhotonsMass[c]->SetMaximum(1.40*plotPhotonsMass[c]->GetMaximum());
    plotPhotonsMass[c]->GetXaxis()->SetTitle("m_{#gamma#gamma} (GeV)");

    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    plotPhotonsMass[c]->Draw();  
    plotPhotonsMass[c]->SetAxisRange(0.000001,plotPhotonsMass[c]->GetMaximum()*1.2,"Y");
    plotPhotonsMass[c]->Draw("SAME");  
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);
    legmc->AddEntry(plotPhotonsMass[c]->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotPhotonsMass[c]->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotPhotonsMass[c]->getObject(3),"Gaussian 2 component","L");
    legmc->AddEntry(plotPhotonsMass[c]->getObject(2),"Gaussian 1 component","L");
    
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    
    TLatex *lat  = new TLatex(mass-20.,0.91*plotPhotonsMass[c]->GetMaximum(),"#scale[1.0]{CMS Preliminary}");
    lat->Draw();
    
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d.root", iMass, c));
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d.png", iMass, c));
    ctmp->SetLogy();
    ctmp->SaveAs("plots/sigmodel_"+TString::Format("%d_cat%d_LOG.png", iMass, c));
  }
 


  // ****************************BKG CAT LOG
  cout << endl; cout << "Progress plotting: background" << endl;     
  TCanvas* c4 = new TCanvas("c4","PhotonsMass Background Categories",0,0,400,400);
  RooPlot* plotPhotonsMassBkg[9];
  for (int c=0; c<ncat; ++c) {
    plotPhotonsMassBkg[c] = PhotonsMass->frame(Range(minMassFit, maxMassFit),nBinsMass);

    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
    PhotonsMassBkg[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
   
     if( blind ) {
      PhotonsMass->setRange("unblind_up",mass+10.,maxMassFit);
      PhotonsMass->setRange("unblind_down",minMassFit,mass-10.);
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_up"));    
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_down"));    
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c], Range("fitrange"));    
      } 
      
    plotPhotonsMassBkg[c]->Draw();  
    gPad->SetLogy(1);
    plotPhotonsMassBkg[c]->SetAxisRange(1.5,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    c4->SaveAs("plots/backgrounds_log_"+TString::Format("cat%d.png", c));
    c4->SaveAs("plots/backgrounds_log_"+TString::Format("cat%d.root", c));
  }
  
  // ****************************BKG LIN
  TCanvas* c5 = new TCanvas("c5","PhotonsMass Background Categories",0,0,400,400);
  RooPlot* plotPhotonsMassBkg[9];
  for (int c=0; c<ncat; ++c) {
    plotPhotonsMassBkg[c] = PhotonsMass->frame(Range(minMassFit, maxMassFit),nBinsMass);
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
    PhotonsMassBkg[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 
   
       if( blind ) {
      PhotonsMass->setRange("unblind_up",mass+10.,maxMassFit);
      PhotonsMass->setRange("unblind_down",minMassFit,mass-10.);
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_up"));    
      data[c]->plotOn(plotPhotonsMassBkg[c], RooFit::CutRange("unblind_down"));    
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c], Range("fitrange"));    
      } 
    

       plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
       plotPhotonsMassBkg[c]->Draw();  
       c5->SaveAs("plots/backgrounds_"+TString::Format("cat%d.png", c));
       c5->SaveAs("plots/backgrounds_"+TString::Format("cat%d.root", c));
  }
}




Double_t effSigma(TH1 *hist) {
  
  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

Double_t effSigma(TH1 *hist) {

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    std::cout << "effsigma: Not a valid histo. nbins = " << nb << std::endl;
    return 0.;
  }

  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    std::cout << "effsigma: Not a valid histo. bwid = " << bwid << std::endl;
    return 0.;
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  if(total < 100.) {
    std::cout << "effsigma: Too few entries " << total << std::endl;
    return 0.;
  }
  Int_t ierr=0;
  Int_t ismin=999;

  Double_t rlim=0.683*total;
  Int_t nrms=rms/(bwid);    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;
    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
    }
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) std::cout << "effsigma: Error of type " << ierr << std::endl;

  return widmin;
}

void R2JJFitter(double mass, std::string postfix="")
{
    filePOSTfix=postfix;
    if(postfix!="")
    {
      // for optimization studies
      MMIN=1000;
      if(mass==1000)
         signalScaler=0.034246;
      if((mass>1000)&&(mass<500))
         signalScaler=0.02469;
      if(mass==2000)
         signalScaler=2.0;
    };
    runfits(mass, 1);
    if(postfix!="")
    {
      // for optimization studies
      MMIN=1000;
      if(mass==1000)
         signalScaler=0.033500;
      if((mass>1000)&&(mass<2000))
         signalScaler=0.02016;
      if(mass==2000)
         signalScaler=2.22222;
    };
    runfits(mass, 0);
    runfits(mass, 2);
}



void MakeParameterTrendvsMassGauss(){

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;
  Float_t masses[5];
  masses[0] = 150.0;
  masses[1] = 200.0;
  masses[2] = 250.0;
  masses[3] = 300.0;
  masses[4] = 400.0.;

  Float_t masses_err[5];
  masses_err[0] = 0.;
  masses_err[1] = 0.;
  masses_err[2] = 0.;
  masses_err[3] = 0.;
  masses_err[4] = 0.;


  TFile* file[5];
  for(imass = 0; imass<5;imass++) file[imass]= new TFile(wsDir+TString::Format("HighMass-hgg.m%.1f_8TeV.inputsig.root",masses[imass]));

  RooWorkspace* ws[5];
  for(imass = 0; imass<5;imass++) ws[imass] = (RooWorkspace*) file[imass]->Get("w_all");

 
  Float_t m0[5];
  Float_t m1[5];
  Float_t sigma0[5];
  Float_t sigma1[5];
  Float_t frac[5];
  Float_t sig[5];
  
  Float_t m0_err[5];
  Float_t m1_err[5];
  Float_t sigma0_err[5];
  Float_t sigma1_err[5];
  Float_t frac_err[5];
  Float_t sig_err[5];

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cat[5];
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
    

  for (int c=0; c<ncat; ++c) { //per ogni categoria

    for(imass = 0; imass<5;imass++){//guardo tutte le masse
      m0[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m0_cat%d",c))->getVal()/masses[imass];
      //  m1[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m1_cat%d",c))->getVal()/masses[imass];
      sigma0[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma0_cat%d",c))->getVal()/masses[imass];
      sigma1[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma1_cat%d",c))->getVal()/masses[imass];
      frac[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_frac_cat%d",c))->getVal();
      //sig[imass]= ws[imass]->var(TString::Format("nsig_cat%d",c))->getVal();
      m0_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m0_cat%d",c))->getError()/masses[imass];
      //  m1_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_m1_cat%d",c))->getError()/masses[imass];
      sigma0_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma0_cat%d",c))->getError()/masses[imass];
      sigma1_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma1_cat%d",c))->getError()/masses[imass];
      frac_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_frac_cat%d",c))->getError();  
      é//sig_err[imass]= ws[imass]->var(TString::Format("nsig_cat%d",c))->getError();   
      std::cout<<masses[imass]<<"    "<<m0[imass]<<"    "<<m0_err[imass]<<std::endl;
  }

 

  
    label_cat[c] = new TPaveText(0.6, 0.8, 0.7, 0.9, "brNDC" );
    label_cat[c]->SetFillColor(kWhite);
    label_cat[c]->SetBorderSize(1.5);
    label_cat[c]->SetTextSize(0.038);
    label_cat[c]->SetTextAlign(11);
    label_cat[c]->SetTextFont(42);
    label_cat[c]->AddText(TString::Format("Cat %d", c));

    TGraphErrors* m0Graph = new TGraphErrors(5,masses, m0, masses_err, m0_err);
    m0Graph->Draw("APE");
    m0Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    m0Graph->GetYaxis()->SetRangeUser(0.8, 1.2);
    m0Graph->GetYaxis()->SetTitle("m0 / m_{#gamma#gamma}");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/m0vsmasses_%d.png", c));

    /*   TGraphErrors* m1Graph = new TGraphErrors(5,masses, m1, masses_err, m1_err);
    m1Graph->Draw("APE");
    m1Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    m1Graph->GetYaxis()->SetTitle("m1 [GeV]");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/m1vsmasses_%d.png", c));*/

    TGraphErrors* sigma0Graph = new TGraphErrors(5,masses, sigma0, masses_err, sigma0_err);
    sigma0Graph->Draw("APE");
    sigma0Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigma0Graph->GetYaxis()->SetTitle("#sigma0 / m_{#gamma#gamma}");
    sigma0Graph->GetYaxis()->SetRangeUser(0., 0.07);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/sigma0vsmasses_%d.png", c));

    TGraphErrors* sigma1Graph = new TGraphErrors(5,masses, sigma1, masses_err, sigma1_err);
    sigma1Graph->Draw("APE");
    sigma1Graph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigma1Graph->GetYaxis()->SetTitle("#sigma1 / m_{#gamma#gamma}");
    sigma1Graph->GetYaxis()->SetRangeUser(0., 0.07);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/sigma1vsmasses_%d.png", c));


    TGraphErrors* fracGraph = new TGraphErrors(5,masses, frac, masses_err, frac_err);
    fracGraph->Draw("APE");
    fracGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    fracGraph->GetYaxis()->SetTitle("%Gauss0 ");
    fracGraph->GetYaxis()->SetRangeUser(0., 2.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/fracvsmasses_%d.png", c));

    /*   TGraphErrors* sigGraph = new TGraphErrors(5,masses, sig, masses_err, sig_err);
    sigGraph->Draw("APE");
    sigGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigGraph->GetYaxis()->SetTitle("Signal yield ");
    sigGraph->GetYaxis()->SetRangeUser(0.001, 100.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/sigYieldvsmasses_%d.png", c));*/




}




}






void MakeParameterTrendvsMassCBC(){
  gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".x /afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/RooCBCrujffPdf.cxx+");
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;
  Float_t masses[5];
  masses[0] = 150.0;
  masses[1] = 200.0;
  masses[2] = 250.0;
  masses[3] = 300.0;
  masses[4] = 400.0.;

  Float_t masses_err[5];
  masses_err[0] = 0.;
  masses_err[1] = 0.;
  masses_err[2] = 0.;
  masses_err[3] = 0.;
  masses_err[4] = 0.;


  TFile* file[5];
  for(imass = 0; imass<5;imass++) file[imass]= new TFile(wsDir+TString::Format("HighMass-hgg.m%.1f_8TeV.inputsig.root",masses[imass]));

  RooWorkspace* ws[5];
  for(imass = 0; imass<5;imass++) ws[imass] = (RooWorkspace*) file[imass]->Get("w_all");

 
  Float_t mean[5];
  Float_t sigma[5];
  Float_t alphaC[5];
  Float_t alphaCB[5];
  Float_t n[5];
  Float_t sig[5];

  Float_t mean_err[5];
  Float_t sigma_err[5];
  Float_t alphaC_err[5];
  Float_t alphaCB_err[5];
  Float_t n_err[5];
  Float_t sig_err[5];
  
  

  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cat[5];
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
    

  for (int c=0; c<ncat; ++c) { //per ogni categoria

    for(imass = 0; imass<5;imass++){//guardo tutte le masse
      mean[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_mean_cat%d",c))->getVal()/masses[imass];
      sigma[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c))->getVal()/masses[imass];
      alphaC[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c))->getVal();
      alphaCB[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c))->getVal();
      n[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_n_cat%d",c))->getVal();
      //  sig[imass]= ws[imass]->var(TString::Format("nsigCBC_cat%d",c))->getVal();
      mean_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_mean_cat%d",c))->getError()/masses[imass];
      sigma_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c))->getError()/masses[imass];
      alphaC_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c))->getError();
      alphaCB_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c))->getError();
      n_err[imass]= ws[imass]->var(TString::Format("PhotonsMass_sig_n_cat%d",c))->getError(); 
      //  sig_err[imass]= ws[imass]->var(TString::Format("nsigCBC_cat%d",c))->getError();    
      std::cout<<masses[imass]<<"    "<<mean[imass]<<"    "<<mean_err[imass]<<std::endl;
  }

 

  
    label_cat[c] = new TPaveText(0.6, 0.8, 0.7, 0.9, "brNDC" );
    label_cat[c]->SetFillColor(kWhite);
    label_cat[c]->SetBorderSize(1.5);
    label_cat[c]->SetTextSize(0.038);
    label_cat[c]->SetTextAlign(11);
    label_cat[c]->SetTextFont(42);
    label_cat[c]->AddText(TString::Format("Cat %d", c));

    TGraphErrors* meanGraph = new TGraphErrors(5,masses, mean, masses_err, mean_err);
    meanGraph->Draw("APE");
    meanGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    meanGraph->GetYaxis()->SetRangeUser(0.8, 1.2);
    meanGraph->GetYaxis()->SetTitle("mean / m_{#gamma#gamma}");
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/meanvsmasses_%d.png", c));

    TGraphErrors* sigmaGraph = new TGraphErrors(5,masses, sigma, masses_err, sigma_err);
    sigmaGraph->Draw("APE");
    sigmaGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigmaGraph->GetYaxis()->SetTitle("sigma / m_{#gamma#gamma}");
    sigmaGraph->GetYaxis()->SetRangeUser(0.,0.1);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/sigmavsmasses_%d.png", c));

    TGraphErrors* alphaCGraph = new TGraphErrors(5,masses, alphaC, masses_err, alphaC_err);
    alphaCGraph->Draw("APE");
    alphaCGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    alphaCGraph->GetYaxis()->SetTitle("#alpha_{C}");
    alphaCGraph->GetYaxis()->SetRangeUser(0., 0.4);
    //Graph->GetYaxis()->SetRangeUser(0., 12.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/alphaCvsmasses_%d.png", c));

    TGraphErrors* alphaCBGraph = new TGraphErrors(5,masses, alphaCB, masses_err, alphaCB_err);
    alphaCBGraph->Draw("APE");
    alphaCBGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    alphaCBGraph->GetYaxis()->SetTitle("#alpha_{CB}");
    alphaCBGraph->GetYaxis()->SetRangeUser(0., 3.);
    //Graph->GetYaxis()->SetRangeUser(0., 20.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/alphaCBvsmasses_%d.png", c));


    TGraphErrors* nGraph = new TGraphErrors(5,masses, n, masses_err, n_err);
    nGraph->Draw("APE");
    nGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    nGraph->GetYaxis()->SetTitle("n ");
    nGraph->GetYaxis()->SetRangeUser(0., 100.);
    //nGraph->GetYaxis()->SetRangeUser(0., 2.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SaveAs(TString::Format("plots/nvsmasses_%d.png", c));

    /*
    TGraphErrors* sigGraph = new TGraphErrors(5,masses, sig, masses_err, sig_err);
    sigGraph->Draw("APE");
    sigGraph->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
    sigGraph->GetYaxis()->SetTitle("Signal yield ");
    sigGraph->GetYaxis()->SetRangeUser(0.001, 100.);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    label_cat[c]->Draw("same");
    c1->SetLogy();
    c1->SaveAs(TString::Format("plots/sigYieldvsmasses_%d.png", c));
    
    */


}




}



void PlotSigShape(RooWorkspace* w, Int_t c) {

  Int_t nmass = 5;
  Int_t masses[5] = {150, 200, 250, 300, 400};

  TString inDir = "";

  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");
  TFile sigFile1("histograms_CMS-HGG_19032013.root");   //ggh prod mode tree livia
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  RooPlot* plot;

  
  
 

  for(int iMass=0;iMass<5;iMass++){ //per ogni massa
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();

  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/ggh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/vbf_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/wzh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/tth_m%d_8TeV", masses[iMass]));
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");

 

  // common preselection cut
  TString mainCut = TString::Format("PhotonsMass>=(%d*0.7) && PhotonsMass<=(%d*1.3)", masses[iMass], masses[iMass]);   // livia

  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
  RooDataSet* signal;

  RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","",TString::Format("@0/%d-1",masses[iMass]),*w->var("PhotonsMass"));

  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);   

  // 1)  prime 4 cat livia
  if (c==0) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
  if (c==1) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
  if (c==2) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
  if (c==3) signal = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));
  

 
  if(iMass==0) plot = massReduced->frame(Range(-0.12, 0.12), Bins(60), Title("Mass Reduced"));

  if(iMass==0) signal->plotOn(plot, DrawOption("L"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kBlack), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==1) signal->plotOn(plot, DrawOption("L"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kPink-8), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==2) signal->plotOn(plot, DrawOption("L"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kBlue+3), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==3) signal->plotOn(plot, DrawOption("L"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kOrange-3), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));
   if(iMass==4) signal->plotOn(plot, DrawOption("L"), MarkerSize(0), MarkerColor(kWhite), LineWidth(2),LineColor(kSpring-8), XErrorSize(0), DataError(RooAbsData::None), Rescale(1./signal->sumEntries()));

  if(iMass==0) plot->Draw();
  else plot->Draw("same");
  }
  plot->SetAxisRange(0.0001,0.4,"Y");
  plot->GetXaxis()->SetTitle("#Delta m ");
  plot->GetXaxis()->SetTitleFont(42);
  plot->GetXaxis()->SetTitleSize(0.05);

  TLegend* legmc = new TLegend(0.62, 0.6, 0.87, 0.89, "", "brNDC");
  legmc->AddEntry(plot->getObject(0),"m_{#gamma#gamma} = 150 GeV","L");
  legmc->AddEntry(plot->getObject(1),"m_{#gamma#gamma} = 200 GeV","L");
  legmc->AddEntry(plot->getObject(2),"m_{#gamma#gamma} = 250 GeV","L");
  legmc->AddEntry(plot->getObject(3),"m_{#gamma#gamma} = 300 GeV","L");
  legmc->AddEntry(plot->getObject(4),"m_{#gamma#gamma} = 400 GeV","L");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw(); 

  TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));
  
  lat->SetTextSize(0.038);
  lat->SetTextAlign(11);
  lat->SetTextFont(42);
  lat->Draw("same");
  lat->SetNDC();


 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");

 
  c1->SaveAs(TString::Format("plots/massesShape_cat%d.png",c));
  c1->SaveAs(TString::Format("plots/massesShape_cat%d.pdf",c));
  c1->SetLogy();
  c1->SaveAs(TString::Format("plots/massesShape_cat%d_LOG.png",c));
  c1->SaveAs(TString::Format("plots/massesShape_cat%d_LOG.pdf",c));
}








void PlotVtxRecoEff() {

  Int_t nmass = 5;
  Int_t masses[5] = {150, 200, 250, 300, 400};
  Double_t massesD[5] = {150., 200., 250., 300., 400.};

  TString inDir = "";

  
  TH1F* massCat0[5];
  TH1F* massCat1[5];
  TH1F* massCat2[5];
  TH1F* massCat3[5];
  TH1F* massCat4[5];

  TH1F* massCat0_vtxOK[5];
  TH1F* massCat1_vtxOK[5];
  TH1F* massCat2_vtxOK[5];
  TH1F* massCat3_vtxOK[5];
  TH1F* massCat4_vtxOK[5];


 
  Double_t effCAT0[5] ;
  Double_t effCAT1[5] ;
  Double_t effCAT2[5] ;
  Double_t effCAT3[5] ;
  Double_t effCAT4[5] ;

  
  for(int iMass=0;iMass<5;iMass++){ //per ogni massa
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();

// common preselection cut
  TString maincut = TString::Format("PhotonsMass>=(%d*0.1) && PhotonsMass<=(%d*1.9)", masses[iMass], masses[iMass]);   // livia
  

  sigTree1->Add(TString::Format("histograms_CMS-HGG_24072013.root/ggh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_24072013.root/vbf_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_24072013.root/wzh_m%d_8TeV", masses[iMass]));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_24072013.root/tth_m%d_8TeV", masses[iMass]));
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");

 
  massCat0[iMass]= new TH1F(TString::Format("massCat0_m%d", masses[iMass]),TString::Format("massCat0_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat1[iMass]= new TH1F(TString::Format("massCat1_m%d", masses[iMass]),TString::Format("massCat1_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat2[iMass]= new TH1F(TString::Format("massCat2_m%d", masses[iMass]),TString::Format("massCat2_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat3[iMass]= new TH1F(TString::Format("massCat3_m%d", masses[iMass]),TString::Format("massCat3_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat4[iMass]= new TH1F(TString::Format("massCat4_m%d", masses[iMass]),TString::Format("massCat4_m%d", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );



  massCat0_vtxOK[iMass]= new TH1F(TString::Format("massCat0_m%d_vtxOK", masses[iMass]),TString::Format("massCat0_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
  
  massCat1_vtxOK[iMass]= new TH1F(TString::Format("massCat1_m%d_vtxOK", masses[iMass]),TString::Format("massCat1_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );

  massCat2_vtxOK[iMass]= new TH1F(TString::Format("massCat2_m%d_vtxOK", masses[iMass]),TString::Format("massCat2_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
  
  massCat3_vtxOK[iMass]= new TH1F(TString::Format("massCat3_m%d_vtxOK", masses[iMass]),TString::Format("massCat3_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
 
  massCat4_vtxOK[iMass]= new TH1F(TString::Format("massCat4_m%d_vtxOK", masses[iMass]),TString::Format("massCat4_m%d_vtxOK", masses[iMass]), 120, masses[iMass]*0.3,masses[iMass]*1.7 );
 
 

  TString cutCat0 = "((abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))";
  TString cutCat1 = "((abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))";
  TString cutCat2 = "((abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))";
  TString cutCat3 = "((abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))";
  TString cutCat4 = "";



  sigTree1->Draw(TString::Format("PhotonsMass>>massCat0_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat0+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat1_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat1+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat2_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat2+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat3_m%d", masses[iMass]),"evweight*("+ maincut+" && "+cutCat3+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat4_m%d", masses[iMass]),"evweight*("+ maincut+")");


  std::cout<<"----------m: "<<masses[iMass]<<" --------"<<std::endl;
  std::cout<<massCat0[iMass]->Integral()<<std::endl;
  std::cout<<massCat1[iMass]->Integral()<<std::endl;
  std::cout<<massCat2[iMass]->Integral()<<std::endl;
  std::cout<<massCat3[iMass]->Integral()<<std::endl;
  std::cout<<massCat4[iMass]->Integral()<<std::endl;

  sigTree1->Draw(TString::Format("PhotonsMass>>massCat0_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat0+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat1_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat1+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat2_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat2+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat3_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1 &&"+cutCat3+")");
  sigTree1->Draw(TString::Format("PhotonsMass>>massCat4_m%d_vtxOK", masses[iMass]),"evweight*("+ maincut+" && abs(gv_z-vtx_z)<1)");

  std::cout<<massCat0_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat1_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat2_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat3_vtxOK[iMass]->Integral()<<std::endl;
  std::cout<<massCat4_vtxOK[iMass]->Integral()<<std::endl;

 

  effCAT0[iMass] = massCat0_vtxOK[iMass]->Integral()/massCat0[iMass]->Integral();
  effCAT1[iMass] = massCat1_vtxOK[iMass]->Integral()/massCat1[iMass]->Integral();
  effCAT2[iMass] = massCat2_vtxOK[iMass]->Integral()/massCat2[iMass]->Integral();
  effCAT3[iMass] = massCat3_vtxOK[iMass]->Integral()/massCat3[iMass]->Integral();
  effCAT4[iMass] = massCat4_vtxOK[iMass]->Integral()/massCat4[iMass]->Integral();

  std::cout<<"------------- eff ------------ "<<std::endl;

  std::cout<<effCAT0[iMass]<<std::endl;
  std::cout<<effCAT1[iMass]<<std::endl;
  std::cout<<effCAT2[iMass]<<std::endl;
  std::cout<<effCAT3[iMass]<<std::endl;
  std::cout<<effCAT4[iMass]<<std::endl;

  }



  Double_t massesErr[5] = {0.,0.,0.,0.,0.};
  Double_t effErr[5] = {0.,0.,0.,0.,0.};
  

  TGraphErrors* effCAT0graph = new TGraphErrors(5, massesD, effCAT0, massesErr, effErr);
  TGraphErrors* effCAT1graph = new TGraphErrors(5, massesD, effCAT1, massesErr, effErr);
  TGraphErrors* effCAT2graph = new TGraphErrors(5, massesD, effCAT2, massesErr, effErr);
  TGraphErrors* effCAT3graph = new TGraphErrors(5, massesD, effCAT3, massesErr, effErr);
  TGraphErrors* effCAT4graph = new TGraphErrors(5, massesD, effCAT4, massesErr, effErr);

  effCAT0graph->SetMarkerColor(kPink-8);
  effCAT1graph->SetMarkerColor(kOrange+1);
  effCAT2graph->SetMarkerColor(kSpring-8);
  effCAT3graph->SetMarkerColor(kBlue-7);
  effCAT4graph->SetMarkerColor(kMagenta-9);
  effCAT0graph->SetLineColor(kPink-8);
  effCAT1graph->SetLineColor(kOrange+1);
  effCAT2graph->SetLineColor(kSpring-8);
  effCAT3graph->SetLineColor(kBlue-7);
  effCAT4graph->SetLineColor(kMagenta-9);

 TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  
  

  TMultiGraph* mg= new TMultiGraph();
  mg->Add(effCAT0graph);
  mg->Add(effCAT1graph);
  mg->Add(effCAT2graph);
  mg->Add(effCAT3graph);
  mg->Add(effCAT4graph);

  mg->Draw("APE");
  mg->GetYaxis()->SetRangeUser(0., 1.);
  mg->GetYaxis()->SetTitle("Vtx MVA Reconstruction efficiency");
  mg->GetXaxis()->SetTitle("m_{#gamma#gamma} [GeV]");
  TLegend* legmc = new TLegend(0.6, 0.2, 0.85, 0.7, "", "brNDC");
  legmc->AddEntry(effCAT0graph, "CAT 0", "PE");
  legmc->AddEntry(effCAT1graph, "CAT 1", "PE");
  legmc->AddEntry(effCAT2graph, "CAT 2", "PE");
  legmc->AddEntry(effCAT3graph, "CAT 3", "PE");
  legmc->AddEntry(effCAT4graph, " All Events", "PE");
  legmc->SetTextSize(0.0206044);  
  legmc->SetTextFont(42);
  legmc->SetBorderSize(0);
  legmc->SetFillStyle(0);
  legmc->Draw("same"); 

  TLatex *lat  = new TLatex(0.65,0.9,"Vertex Recnstruction Efficiency");
  
  lat->SetTextSize(0.038);
  lat->SetTextAlign(11);
  lat->SetTextFont(42);
  // lat->Draw("same");
  lat->SetNDC();


 

  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c1->SetLogy(0);
  c1->SaveAs("plots/VtxRecoEff.png");
  c1->SaveAs("plots/VtxRecoEff.pdf");
 }







void MakePlotMassDataMC(RooWorkspace* w, Float_t mass ) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  RooRealVar* PhotonsMass = w->var("PhotonsMass");
   
  int iMass = abs(mass);      
 // common preselection cut
  TString mainCut = "(PhotonsMass>=130 && PhotonsMass<=1200)";   //130-2000
  RooPlot*  plotPhotonsMassDataMC[NCAT];
  
  //**********DATA***************//
  // create tree
  TFile file("histograms_CMS-HGG_24072013.root");
  TTree* dataTree = (TTree*) file.Get("Data");
   
  //**********G+jets***************//
  // create tree

  TChain* gjTree = new TChain();
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");

  TChain* qcdTree = new TChain();
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf"); 
 
  //**********DIPHOTJET***************//
  // create tree

  TChain* diphotjTree = new TChain();
  diphotjTree->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
   

  //**********DIPHOT***************//
  // create tree

  TChain* diphotTree = new TChain();
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
   

  

  TH1F* h_data[NCAT+1];
  TH1F* h_data_b[NCAT+1];
  TH1F*  h_gj[NCAT+1];
  TH1F*  h_qcd[NCAT+1];
  TH1F*  h_diphot[NCAT+1];
  TH1F*  h_diphotj[NCAT+1];
  TH1F* h_sum[NCAT+1]; 
 
  
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",1);
  Int_t nbin = 87;
  Double_t  min = 130;
  Double_t  max = 1000;

   for (int c=0; c<NCAT+1; ++c) {


     h_data[c]= new TH1F(TString::Format("h_data_cat%d",c), TString::Format("h_data_cat%d",c), nbin, min, max);
     h_data_b[c]= new TH1F(TString::Format("h_data_b_cat%d",c), TString::Format("h_data_b_cat%d",c), nbin, min, max);
     h_gj[c]= new TH1F(TString::Format("h_gj_cat%d",c), TString::Format("h_gj_cat%d",c), nbin, min, max);
     h_qcd[c]= new TH1F(TString::Format("h_qcd_cat%d",c), TString::Format("h_qcd_cat%d",c), nbin, min, max);
     h_diphot[c]= new TH1F(TString::Format("h_diphot_cat%d",c), TString::Format("h_diphot_cat%d",c), nbin, min, max);
     h_diphotj[c]= new TH1F(TString::Format("h_diphotj_cat%d",c), TString::Format("h_diphotj_cat%d",c), nbin, min, max);



    // 1)  prime 4 cat livia
     if (c==0){//&&(PhotonsMass<178 || PhotonsMass >402) 
   
      dataTree->Draw("PhotonsMass>>h_data_cat0", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
      dataTree->Draw("PhotonsMass>>h_data_b_cat0", "("+mainCut+"&&(PhotonsMass<180 || PhotonsMass >850)&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
      gjTree->Draw("PhotonsMass>>h_gj_cat0", "("+mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      gjTree->Draw("PhotonsMass>>h_qcd_cat0", "("+mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotTree->Draw("PhotonsMass>>h_diphot_cat0","("+ mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotjTree->Draw("PhotonsMass>>h_diphotj_cat0","("+ mainCut+"&&(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
     }
     if (c==1){//&&(PhotonsMass<178 || PhotonsMass >402)&&

       dataTree->Draw("PhotonsMass>>h_data_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 )  )*1");
       dataTree->Draw("PhotonsMass>>h_data_b_cat1", "("+mainCut+"&&(PhotonsMass<180 || PhotonsMass >850)&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 )  )*1");
     gjTree->Draw("PhotonsMass>>h_gj_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
     qcdTree->Draw("PhotonsMass>>h_qcd_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
      diphotTree->Draw("PhotonsMass>>h_diphot_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");
     diphotjTree->Draw("PhotonsMass>>h_diphotj_cat1", "("+mainCut+"&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");

       }
     if (c==2){//&&(PhotonsMass<178 || PhotonsMass >402)

     dataTree->Draw("PhotonsMass>>h_data_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
     dataTree->Draw("PhotonsMass>>h_data_b_cat2", "("+mainCut+"&&(PhotonsMass<180 || PhotonsMass >850)&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
     gjTree->Draw("PhotonsMass>>h_gj_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     qcdTree->Draw("PhotonsMass>>h_qcd_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotTree->Draw("PhotonsMass>>h_diphot_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotjTree->Draw("PhotonsMass>>h_diphotj_cat2", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
   }
     if (c==3){//&&(PhotonsMass<178 || PhotonsMass >402)

     dataTree->Draw("PhotonsMass>>h_data_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*1");
     dataTree->Draw("PhotonsMass>>h_data_b_cat3", "("+mainCut+"&&(PhotonsMass<180 || PhotonsMass >850)&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*1");
     gjTree->Draw("PhotonsMass>>h_gj_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
     qcdTree->Draw("PhotonsMass>>h_qcd_cat3", "("+mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
      diphotTree->Draw("PhotonsMass>>h_diphot_cat3","("+ mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ) )*evweight");
      diphotjTree->Draw("PhotonsMass>>h_diphotj_cat3","("+ mainCut+"&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
    
   }
     if (c==4){//&&(PhotonsMass<178 || PhotonsMass >402)

    dataTree->Draw("PhotonsMass>>h_data_cat4","("+ mainCut+")*1");
    dataTree->Draw("PhotonsMass>>h_data_b_cat4","("+ mainCut+"&&(PhotonsMass<180 || PhotonsMass >850))*1");
    gjTree->Draw("PhotonsMass>>h_gj_cat4", "("+mainCut+")*evweight");
    qcdTree->Draw("PhotonsMass>>h_qcd_cat4", "("+mainCut+")*evweight");
    diphotTree->Draw("PhotonsMass>>h_diphot_cat4","("+ mainCut+")*evweight");
    diphotjTree->Draw("PhotonsMass>>h_diphotj_cat4","("+ mainCut+")*evweight");
    

   }

 
     //fisso la normalizzazione di gjets a quella gj+qcd xkè la shape di qcd _pf fa schifo, quindi prenidamoq uella di gj _pf
     Double_t qcdInt = h_qcd[c]->Integral();
     Double_t gJetIntFixed = qcdInt+h_gj[c]->Integral();
     h_gj[c]->Scale(gJetIntFixed/h_gj[c]->Integral());     
     h_gj[c]->SetFillColor(kAzure+8);
     h_gj[c]->Sumw2();
     
     h_diphotj[c]->SetFillColor(kTeal+9);
     h_diphotj[c]->Sumw2();
     
     //sommo i due diphot
     h_diphot[c]->Add(h_diphotj[c]);
     h_diphot[c]->SetFillColor(kSpring+7);
     h_diphot[c]->Sumw2();
     
     
     
     h_sum[c] = (TH1F*) h_gj[c]->Clone();
     h_sum[c]->Add(h_diphot[c]);
     
     h_sum[c]->SetFillColor(kBlack);
     h_sum[c]->SetFillStyle(3004);
     h_sum[c]->SetMarkerSize(0);

     //make kolmogorov test between data and MC
     Double_t CHI2ndf = h_data[c]->Chi2Test(h_sum[c], "UWPCHI2/NDF");
     TPaveText* label = new TPaveText(0.6, 0.67, 0.85, 0.7, "brNDC" );
     label->SetFillColor(kWhite);
     label->SetBorderSize(0.);
     label->SetTextSize(0.038);
     label->SetTextAlign(11);
     label->SetTextFont(42);
     label->AddText(TString::Format("#chi^{2}/NDF: %.3f", CHI2ndf));
     
     THStack hs("hs","hs");
   
   
     hs.Add(h_diphot[c]); 
     hs.Add(h_gj[c]);
   
     std::cout<<"------"<<std::endl;
     std::cout<<h_sum[c]->Integral()<<std::endl;
     std::cout<<h_data[c]->Integral()<<std::endl;
   
 
    ctmp->cd();
    h_data[c]->Sumw2();
    h_data_b[c]->Sumw2();
    h_sum[c]->Sumw2();

    TH1F* h1_ratio1 = (TH1F*)h_data_b[c]->Clone();
    TH1F* h1_ratio1_unblind = (TH1F*)h_data[c]->Clone();
    TH1F* h1_ratio2 = (TH1F*)h_sum[c]->Clone();

    for(int i = 0;i < h1_ratio1->GetNbinsX();i++) std::cout<<" ratio1: "<<h1_ratio1->GetBinContent(i)<<std::endl;//if(h1_ratio1->GetBinContent(i)==0)h1_ratio1->SetBinContent(i, -1);
   
   

    std::cout<<"int 700-850: "<<h_data[c]->Integral(h_data[c]->FindBin(700.),h_data[c]->FindBin(950.) )<<"  int 850-1200: "<<h_data[c]->Integral(h_data[c]->FindBin(950.),h_data[c]->FindBin(1200.) )<<std::endl;



    ctmp->Clear();
    //-------pad 1-------//
    TPad * pad1 = new TPad("pad1", "pad1",0.01,0.13,0.75,1.);  
    
   
    pad1->SetRightMargin(0.1);
    
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();

    h_data_b[c]->SetMarkerSize(0.7);
    h_data_b[c]->Draw("pe");
    
    h_data_b[c]->GetYaxis()->SetTitle(TString::Format("Events/%.2f", (max-min)/nbin));
    //h__data[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    h_data_b[c]->GetYaxis()->SetRangeUser(0.1, h_data_b[c]->GetMaximum()*2);
    hs.Draw("histsame");
    h_data_b[c]->Draw("pesame");
    h_sum[c]->Draw("E2same");
    std::cout<<"--------_> "<<h_data_b[c]->Integral()<<std::endl; 
  
  

    TLegend *leg1;
    if(c!=4)leg1 = new TLegend(0.6075,0.7536441,0.8575,0.9340678, TString::Format("Category %d",c), "brNDC");
    else leg1 = new TLegend(0.6075,0.7536441,0.8575,0.9340678, TString::Format("All Events",c), "brNDC");
    leg1->AddEntry(h_data_b[c],"Data","PE");
    leg1->AddEntry(h_gj[c],"prompt +fake","F");
    //leg1->AddEntry(h_diphotj[c],"DiPhoton + jets", "F");
    leg1->AddEntry(h_diphot[c],"prompt +prompt", "F");
    leg1->AddEntry(h_sum[c], "Bkg Err", "F");
    
   
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");

    label->Draw("same");
 

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    pad1->SetLogy(0);
    pad1->RedrawAxis();

    ctmp->cd();

    //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.001,0.75,0.2);
    pad2->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
    pad2->cd();

    Double_t xmax = h1_ratio1->GetXaxis()->GetXmax();
    Double_t xmin = h1_ratio1->GetXaxis()->GetXmin();
    TLine* line = new TLine(xmin,1.,xmax,1.);
  

    h1_ratio1->SetStats(0);
    
    h1_ratio1->Divide(h_sum[c]);
    h1_ratio1_unblind->Divide(h_sum[c]);
    h1_ratio2->Divide(h_sum[c]);
    h1_ratio1->SetMarkerStyle(20);
    h1_ratio1->SetMarkerSize(1.1);
    //  h1_ratio1->GetXaxis()->SetTitle(xAxis.c_str());
    h1_ratio1->GetYaxis()->SetRangeUser(0.0001, 2.); // cwr zoom
    h1_ratio1->GetYaxis()->SetNdivisions(2,false);
    h1_ratio1->GetYaxis()->SetTitle("Data/Bkg.");
    h1_ratio1->GetYaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    
    h1_ratio1->GetXaxis()->SetTitleSize(0.2);
    h1_ratio1->GetXaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetTitleSize(0.15);
    h1_ratio1->GetYaxis()->SetTitleOffset(0.25);
    h1_ratio1->GetXaxis()->SetTitleOffset(0.5);

    
    
    for(int j = 0;j<=h1_ratio1->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j))  h1_ratio1->SetBinError(j,h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio1->SetBinError(j,0.);
      if(h_sum[c]->GetBinContent(j))  h1_ratio1->SetBinError(j,sqrt(pow(h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j), 2)+ pow(h_data[c]->GetBinError(j)*h_sum[c]->GetBinError(j)/(h_sum[c]->GetBinContent(j)*h_sum[c]->GetBinContent(j)),2)));
      else h1_ratio1->SetBinError(j,0.);
    }
    h1_ratio1->Draw("PEX0");
   
    for(int j = 0;j<=h1_ratio1_unblind->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j))  h1_ratio1_unblind->SetBinError(j,h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio1->SetBinError(j,0.);
      if(h_sum[c]->GetBinContent(j))  h1_ratio1_unblind->SetBinError(j,sqrt(pow(h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j), 2)+ pow(h_data[c]->GetBinError(j)*h_sum[c]->GetBinError(j)/(h_sum[c]->GetBinContent(j)*h_sum[c]->GetBinContent(j)),2)));
      else h1_ratio1_unblind->SetBinError(j,0.);
    }
  
    for(int j = 0;j<=h1_ratio2->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j)) h1_ratio2->SetBinError(j,h_sum[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio2->SetBinError(j,0.);
    }
    h1_ratio2->Draw("E2same");
    
    line->SetLineWidth(1.);
    line->Draw("same");
    

    ctmp->cd();
    //-------pad 3------//
    
   
    TPad * pad3 = new TPad("pad3", "pad3",0.68,0.001,1.,0.2);
    pad3->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad3->SetBottomMargin(0.4);
    pad3->SetRightMargin(0.1);
    pad3->Draw();
    pad3->cd();

    TH1F* h1_res = new TH1F("h1_res", "h1_res",11,-5, 5.);
   
    for(int i = 1;i < h1_ratio1_unblind->GetNbinsX();i++){
      if(h1_ratio1_unblind->GetBinError(i)==0)continue;
      std::cout<<((h1_ratio1_unblind->GetBinContent(i))-1)/h1_ratio1_unblind->GetBinError(i)<<std::endl;
      h1_res->SetFillColor(kBlack);
      h1_res->Fill(((h1_ratio1_unblind->GetBinContent(i))-1)/h1_ratio1_unblind->GetBinError(i));

    }
    std::cout<<h1_res->GetEntries()<<std::endl;
   
    h1_res->GetYaxis()->SetNdivisions(4,false);
    h1_res->GetXaxis()->SetTitle("pull");
    h1_res->GetXaxis()->SetTitleSize(0.15);
    h1_res->GetXaxis()->SetTitleOffset(0.3);
    h1_res->Draw("hbar");


    //-------pad 3------//
    
    ctmp->cd();
    TPad * pad3 = new TPad("pad3", "pad3",0.68,0.2,1.,1);
   
    
    //pad2->SetTopMargin(0.01);
    pad3->SetBottomMargin(0.4);
    pad3->SetRightMargin(0.1);
    pad3->Draw();
    pad3->cd();

    Double_t mean = h1_res->GetMean();
    Double_t sigma = h1_res->GetRMS();

    
    TPaveText* label2 = new TPaveText(0.2, -0.3, 0.65, 0.35, "brNDC" );
    label2->SetFillColor(kWhite);
    label2->SetBorderSize(0.);
    label2->SetTextSize(0.058);
    label2->SetTextAlign(11);
    label2->SetTextFont(42);
    label2->AddText(TString::Format("Mean: %.3f; RMS: %.3f", mean, sigma));


    label2->Draw("same");

    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d.pdf", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d.root", c));
  
    pad1->SetLogy(1);
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d_LOG_ALL_range130-1000.png", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("plots/DATA_MC_MASS_"+TString::Format("cat%d_LOG.root", c));
  }

 
  



}







void MakePlotNVTXDataMC(RooWorkspace* w, Float_t mass ) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  RooRealVar* nvtx = w->var("nvtx");
   
  int iMass = abs(mass);      
 // common preselection cut
  TString mainCut = "";   //130-500
 
  
  //**********DATA***************//
  // create tree
  TFile file("histograms_CMS-HGG_24072013.root");
  TTree* dataTree = (TTree*) file.Get("Data");
   
  //**********G+jets***************//
  // create tree

  TChain* gjTree = new TChain();
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf"); 

  TChain* qcdTree = new TChain();
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf"); 
 
  //**********DIPHOTJET***************//
  // create tree

  TChain* diphotjTree = new TChain();
  diphotjTree->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
   

  //**********DIPHOT***************//
  // create tree

  TChain* diphotTree = new TChain();
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  diphotTree->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
   

  

  TH1F* h_data[NCAT+1];
  TH1F*  h_gj[NCAT+1];
  TH1F*  h_qcd[NCAT+1];
  TH1F*  h_diphot[NCAT+1];
  TH1F*  h_diphotj[NCAT+1];
  TH1F* h_sum[NCAT+1]; 
 
  
  TCanvas* ctmp = new TCanvas("ctmp","nvtx Background Categories",1);
  Int_t nbin = 50;
  Double_t  min = 0;
  Double_t  max = 50;

   for (int c=0; c<NCAT+1; ++c) {


     h_data[c]= new TH1F(TString::Format("h_data_cat%d",c), TString::Format("h_data_cat%d",c), nbin, min, max);
     h_gj[c]= new TH1F(TString::Format("h_gj_cat%d",c), TString::Format("h_gj_cat%d",c), nbin, min, max);
     h_qcd[c]= new TH1F(TString::Format("h_qcd_cat%d",c), TString::Format("h_qcd_cat%d",c), nbin, min, max);
     h_diphot[c]= new TH1F(TString::Format("h_diphot_cat%d",c), TString::Format("h_diphot_cat%d",c), nbin, min, max);
     h_diphotj[c]= new TH1F(TString::Format("h_diphotj_cat%d",c), TString::Format("h_diphotj_cat%d",c), nbin, min, max);



    // 1)  prime 4 cat livia
     if (c==0){//&&(PhotonsMass<178 || PhotonsMass >402) 
   
      dataTree->Draw("nvtx>>h_data_cat0", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
      gjTree->Draw("nvtx>>h_gj_cat0", "("+mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      qcdTree->Draw("nvtx>>h_qcd_cat0", "("+mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotTree->Draw("nvtx>>h_diphot_cat0","("+ mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
      diphotjTree->Draw("nvtx>>h_diphotj_cat0","("+ mainCut+"(abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
     }
     if (c==1){//&&(nvtx<178 || nvtx >402)&&

     dataTree->Draw("nvtx>>h_data_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 )  )*1");
     gjTree->Draw("nvtx>>h_gj_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
     qcdTree->Draw("nvtx>>h_qcd_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ))*evweight ");
      diphotTree->Draw("nvtx>>h_diphot_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");
     diphotjTree->Draw("nvtx>>h_diphotj_cat1", "("+mainCut+" (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) )*evweight");

       }
     if (c==2){//&&(nvtx<178 || nvtx >402)

     dataTree->Draw("nvtx>>h_data_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*1");
     gjTree->Draw("nvtx>>h_gj_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     qcdTree->Draw("nvtx>>h_qcd_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotTree->Draw("nvtx>>h_diphot_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
     diphotjTree->Draw("nvtx>>h_diphotj_cat2", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 ))*evweight");
    
   }
     if (c==3){//&&(nvtx<178 || nvtx >402)

     dataTree->Draw("nvtx>>h_data_cat3", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*1");
     gjTree->Draw("nvtx>>h_gj_cat3", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
     qcdTree->Draw("nvtx>>h_qcd_cat3", "("+mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
      diphotTree->Draw("nvtx>>h_diphot_cat3","("+ mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ) )*evweight");
      diphotjTree->Draw("nvtx>>h_diphotj_cat3","("+ mainCut+" (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ))*evweight ");
    
   }
     if (c==4){//&&(nvtx<178 || nvtx >402)

    dataTree->Draw("nvtx>>h_data_cat4","("+ mainCut+"1)*1");
    gjTree->Draw("nvtx>>h_gj_cat4", "("+mainCut+"1)*evweight");
    qcdTree->Draw("nvtx>>h_qcd_cat4", "("+mainCut+"1)*evweight");
    diphotTree->Draw("nvtx>>h_diphot_cat4","("+ mainCut+"1)*evweight");
    diphotjTree->Draw("nvtx>>h_diphotj_cat4","("+ mainCut+"1)*evweight");
    

   }

 
     //fisso la normalizzazione di gjets a quella gj+qcd xkè la shape di qcd _pf fa schifo, quindi prenidamoq uella di gj _pf
     Double_t qcdInt = h_qcd[c]->Integral();
     Double_t gJetIntFixed = qcdInt+h_gj[c]->Integral();
     h_gj[c]->Scale(gJetIntFixed/h_gj[c]->Integral());   
     h_gj[c]->SetFillColor(kAzure+8);
     h_gj[c]->Sumw2();
   
     h_diphotj[c]->SetFillColor(kTeal+9);
     h_diphotj[c]->Sumw2();

   //sommo i due diphot samples
   h_diphot[c]->Add(h_diphotj[c]);
   h_diphot[c]->SetFillColor(kSpring+7);
   h_diphot[c]->Sumw2();


   
   h_sum[c] = (TH1F*) h_gj[c]->Clone();
   h_sum[c]->Add(h_diphot[c]);
   
   h_sum[c]->SetFillColor(kBlack);
   h_sum[c]->SetFillStyle(3004);
   h_sum[c]->SetMarkerSize(0);
   

   THStack hs("hs","hs");  
   hs.Add(h_diphot[c]); 
   hs.Add(h_gj[c]);
   
   
   std::cout<<h_sum[c]->Integral()<<std::endl;
   std::cout<<h_data[c]->Integral()<<std::endl;
   
 
    ctmp->cd();
      h_data[c]->Sumw2();
    h_sum[c]->Sumw2();
    TH1F* h1_ratio1 = (TH1F*)h_data[c]->Clone();
    TH1F* h1_ratio2 = (TH1F*)h_sum[c]->Clone();
 
    
    //-------pad 1-------//
    TPad * pad1 = new TPad("pad1", "pad1",0.01,0.13,1.,1.);  
    
    //pad1->SetTopMargin(0.1);
    //pad1->SetBottomMargin(0.01);
    pad1->SetRightMargin(0.1);
    
    pad1->SetLogy();
    pad1->Draw();
    pad1->cd();
   
   
    h_data[c]->SetMarkerSize(0.7);
    h_data[c]->Draw("pe");
    
    h_data[c]->GetYaxis()->SetTitle(TString::Format("Events/%.2f", (max-min)/nbin));
    h_data[c]->GetXaxis()->SetTitle("nvtx");
    //  h_data[c]->GetYaxis()->SetRangeUser(0., 28000);
    hs.Draw("histsame");
    h_data[c]->Draw("pesame");
    h_sum[c]->Draw("E2same");
    std::cout<<"--------_> "<<h_data[c]->Integral()<<std::endl; 
  

    TLegend *leg1;
    if(c!=4)leg1 = new TLegend(0.6075,0.7536441,0.8575,0.9340678, TString::Format("Category %d",c), "brNDC");
    else leg1 = new TLegend(0.6075,0.7536441,0.8575,0.9340678, TString::Format("All Events",c), "brNDC");
    leg1->AddEntry(h_data[c],"Data","PE");
    leg1->AddEntry(h_gj[c],"prompt + fake","F");
    leg1->AddEntry(h_diphot[c],"prompt + prompt", "F");
    leg1->AddEntry(h_sum[c], "Bkg Err", "F");
    
   
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
    leg1->Draw("same");
 

    TPaveText* label_cms = get_labelCMS(0, "2012", false);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    pad1->SetLogy(0);
    pad1->RedrawAxis();

    ctmp->cd();

    //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.001,1.,0.2);
    pad2->SetGrid();
    
    //pad2->SetTopMargin(0.01);
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
    pad2->cd();

    Double_t xmax = h1_ratio1->GetXaxis()->GetXmax();
    Double_t xmin = h1_ratio1->GetXaxis()->GetXmin();
    TLine* line = new TLine(xmin,1.,xmax,1.);
  

    h1_ratio1->SetStats(0);
    h1_ratio1->Divide(h_sum[c]);
    h1_ratio1->SetMarkerStyle(20);
    h1_ratio1->SetMarkerSize(1.1);
    //  h1_ratio1->GetXaxis()->SetTitle(xAxis.c_str());
    h1_ratio1->GetYaxis()->SetRangeUser(0., 2.); // cwr zoom
    h1_ratio1->GetYaxis()->SetNdivisions(4,false);
    h1_ratio1->GetYaxis()->SetTitle("Data/Bkg.");
    h1_ratio1->GetYaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitle("nvtx");
    
    h1_ratio1->GetXaxis()->SetTitleSize(0.23);
    h1_ratio1->GetXaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetTitleSize(0.15);
    h1_ratio1->GetYaxis()->SetTitleOffset(0.25);
    h1_ratio1->GetXaxis()->SetTitleOffset(0.5);

    
    
    for(int j = 0;j<=h1_ratio1->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j))  h1_ratio1->SetBinError(j,h_data[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio1->SetBinError(j,0.);
    }
    h1_ratio1->Draw("PEX0");
    
    h1_ratio2->Divide(h_sum[c]);
    for(int j = 0;j<=h1_ratio2->GetNbinsX();j++){
      if(h_sum[c]->GetBinContent(j)) h1_ratio2->SetBinError(j,h_sum[c]->GetBinError(j)/h_sum[c]->GetBinContent(j));
      else h1_ratio2->SetBinError(j,0.);
    }
    h1_ratio2->Draw("E2same");
    
    line->SetLineWidth(1.);
    line->Draw("same");


    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d.pdf", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d.root", c));
  
    pad1->SetLogy(1);
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d_LOG.pdf", c));
    ctmp->SaveAs("plots/DATA_MC_NVTX_"+TString::Format("cat%d_LOG.root", c));
  }

 
  



}




void MakePlotRooKeys(int cat){
  TString wsDir = "BiasStudy/workspaces/";
  TFile* f = new TFile(""+wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV.root");
  RooWorkspace *w = f->Get("w_bias");
  
  RooAbsPdf* r2 = (RooAbsPdf*) *w->pdf(TString::Format("BkgMCKeyPdf_bw2_cat%",cat));
  /*RooKeysPdf* r2 = (RooKeysPdf*) *w->pdf(TString::Format("BkgMCKeyPdf_bw3_cat%",cat));
  RooKeysPdf* r2 = (RooKeysPdfx*) *w->pdf(TString::Format("BkgMCKeyPdf_bw3_cat%",cat));*/
  RooDataSet* mc =(RooDataSet*) *w->data(TString::Format("BkgMCWeight_cat%d",cat));

  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background RooKeysPdf",0,0,500,500);
  RooPlot* plotPhotonsMassBkgMC;
  Int_t nBinsMass(320);
  plotPhotonsMassBkgMC = w->var("PhotonsMass")->frame(130, 1000,nBinsMass);

  r2->plotOn(plotPhotonsMassBkgMC,"L", LineColor(kBlack), LineWidth(2));
  /* r3->plotOn(plotPhotonsMassBkgMC,"L", LineColor(8), LineWidth(2));
  r4->plotOn(plotPhotonsMassBkgMC,"L", LineColor(kOrange), LineWidth(2));*/
  mc->plotOn(plotPhotonsMassBkgMC);
  plotPhotonsMassBkgMC->SetAxisRange(0.001,plotPhotonsMassBkgMC->GetMaximum()*30.,"Y");
    

  TLegend *leg1 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",cat), "brNDC");
  leg1->AddEntry(plotPhotonsMassBkgMC->getObject(0),"Bkg MC","LPE");
  TLegend *leg2 = new TLegend(0.4375,0.7236441,0.85,0.9240678, TString::Format("RooKeysPdf",cat), "brNDC");
  
  leg2->AddEntry(plotPhotonsMassBkgMC->getObject(1),"Bw x 2","L");
  leg2->AddEntry(plotPhotonsMassBkgMC->getObject(2),"Bw x 3","L");
  leg2->AddEntry(plotPhotonsMassBkgMC->getObject(3),"Bw x 4","L");
  
  leg1->SetTextSize(0.035);
  leg1->SetTextFont(42);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  
  leg2->SetTextSize(0.035);
  leg2->SetTextFont(42);
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0);
  
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  ctmp->cd();
  
  plotPhotonsMassBkgMC->Draw();  
  leg1->Draw("same");
  leg2->Draw("same");
  
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  
  ctmp->SetLogy(1);
  
  ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.png", cat));
  ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.root", cat));
}
