
 
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


Int_t MINmass= 130;
Int_t MAXmass= 1000;

static const Int_t NCAT = 4;  // chiara

std::string filePOSTfix="";
double signalScaler=1.00;

void AddSigData(RooWorkspace*, Float_t);
void AddBkgData(RooWorkspace*, Float_t);

void SigModelFitGauss(RooWorkspace*, Float_t);
void SigModelFitCBC(RooWorkspace*, Float_t);

RooHistFunc* getRooHistFunc(int cat, RooRealVar* var);
void SigModelResponseFcnFit(RooWorkspace*);
void SigModelFitConvBW(RooWorkspace*, Float_t, Double_t);


RooFitResult*  BkgModelFitBernstein(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitDiJetFunc(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitLauFunc(RooWorkspace*, Bool_t, Float_t, bool);
RooFitResult*  BkgModelFitExpolFunc(RooWorkspace*, Bool_t, Float_t, bool);



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
RooArgSet* defineVariablesM250() {

  // define variables of the input ntuple //livia
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",130, 1000,"GeV");
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


RooArgSet* defineVariables_newWeight() {

  // define variables of the input ntuple //livia
  RooRealVar* PhotonsMass  = new RooRealVar("PhotonsMass", "M(gg)",MINmass, MAXmass,"GeV");
  RooRealVar* ph1_eta = new RooRealVar("ph1_eta", "eta(g1)",-10,10,"");
  RooRealVar* ph2_eta = new RooRealVar("ph2_eta", "eta(g2)",-10,10,"");
  RooRealVar* ph1_r9 = new RooRealVar("ph1_r9", "R9(g1)",-10,10,"");
  RooRealVar* ph2_r9 = new RooRealVar("ph2_r9", "R9(g2)",-10,10,"");
  RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  RooRealVar* newweight = new RooRealVar("newweight","Reweightings",0., 1000,"");
 
  //  RooRealVar* btagCategory = new RooRealVar("btagCategory","event category",0.9,2.1,"") ;
  
  RooArgSet* ntplVars = new RooArgSet(*PhotonsMass, *ph1_eta, *ph2_eta, *ph1_r9, *ph2_r9,*newweight,   *nvtx);
  
  return ntplVars;
}



TTree *dataset2tree(RooDataSet *dataset, RooArgSet* args, Double_t newweight){

  
  RooArgList argList(*args);

  Double_t variables[50];
  Long64_t nEntries= dataset->numEntries();
  //nEntries=1;
  TTree *tree = new TTree("tree","tree");
  tree->SetDirectory(0);
  TIterator *it1=NULL; 
  it1 = argList.createIterator();
  int index1=0;
  for(RooRealVar *var = (RooRealVar *) it1->Next(); var!=NULL;
      var = (RooRealVar *) it1->Next(),index1++){
    TString name(var->GetName());
    name.ReplaceAll("-","_");
    tree->Branch(name, &(variables[index1]), name+"/D");
  }

  
  tree->Branch("newweight", &newweight, "newweight/D");
  //  tree->Print();

  for(Long64_t jentry=0; jentry<nEntries; jentry++){
    
    TIterator *it=NULL; 
    RooArgList argList1(*(dataset->get(jentry)));
    it = argList1.createIterator();
    //(dataset->get(jentry))->Print();
    int index=0;
    for(RooRealVar *var = (RooRealVar *) it->Next(); var!=NULL;
	var = (RooRealVar *) it->Next(), index++){
      variables[index]=var->getVal();
      //var->Print();
      if(var->GetName() == "evweight") newweight = evweight*2;
    }
   
    delete it;
    tree->Fill();
  }
  tree->ResetBranchAddresses();
  //  tree->Print();
  return tree;
}




void runfits(const Float_t mass=150, Bool_t dobands = false, Int_t perc = 0, Int_t min=130, Int_t max=1000) {

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
  
  //compure roohistfunc with min and max for the fir range 
  RooRealVar* var = new RooRealVar("var","var", 100, 900 );
  TF1* fmax;
  TF1* fmin;
  RooHistFunc* rooFitMin = getRooHistFuncFitMIN(0, var, fmin);
  RooHistFunc* rooFitMax = getRooHistFuncFitMAX(0, var, fmax);
  
  //  makeRooHistPlot(rooFitMin,rooFitMax,fmin,fmin,var);
  var->setVal(mass);
  double newmin = rooFitMin->getVal(*var);
  double newmax = rooFitMax->getVal(*var);
  newmin =200;
  newmax=500;
  MINmass=newmin;
  MAXmass=newmax;
  
  std::cout<<"  MIN MASS: "<<MINmass<<"   MAX MASS: "<<MAXmass<<std::endl;
  
  Double_t MMIN = 130.;
  Double_t MMAX = 1000; 
  w->var("PhotonsMass")->setMin(MMIN);
  w->var("PhotonsMass")->setMax(MMAX);
 
  
  
  // w->Print("v");
  // Add data to the workspace
  cout << endl; cout << "Now AddSigData" << endl;
  if(mass<400 && mass!=350) AddSigData(w, mass);   
  
  // Inside this function you will find a discription our model.
  cout << endl; cout << "Now SigModelFit" << endl;

  //  SigModelFitCBC(w, mass);     
  SigModelResponseFcnFit(w, mass); 

  w->var("PhotonsMass")->setMin(MINmass);
  w->var("PhotonsMass")->setMax(MAXmass);
  

  SigModelFitConvBW(w, mass, perc);      

 
  cout << endl; cout << "Now AddBkgData" << endl;
  AddBkgData(w, mass);         
  
  cout << endl; cout << "Now AddBkgMC" << endl;
  //MakeRooKeysPDFMCBkg(w, mass, true);


  cout << endl; cout << "Now BkgModelFit" << endl;    
  bool blind = true; 
  bool dobandsHere= false;

  //fitresults_bern = BkgModelFitBernstein(w, dobands, mass, blind);  
  
  fitresults_expol = BkgModelFitExpolFunc(w, dobands, mass, blind);
  fitresults_lau = BkgModelFitLauFunc(w, dobands, mass, blind);  
  fitresults_exPAR = BkgModelFitExpPARFunc(w, dobands, mass, blind);  
  //  fitresults_expolPL = BkgModelFitExpolPLFunc(w, dobands, mass, blind);  
  fitresults_dijet = BkgModelFitDiJetFunc(w, dobandsHere, mass, blind);     
  // fitresults_dijetPL = BkgModelFitDiJetPLFunc(w, dobandsHere, mass, blind);     
  //fitresults_dijetEXP = BkgModelFitDiJetEXPFunc(w, dobandsHere, mass, blind);     
  //  fitresults_dijetEXPOL = BkgModelFitDiJetEXPOLFunc(w, dobandsHere, mass, blind);     
  // fitresults_rookey = BkgModelFitRooKeyFunc(w, dobandsHere, mass, blind);     
  
  
   
  // Make statistical treatment. Setup the limit on X production
  cout << endl; cout << "Now make signal workspace" << endl;
  bool isGauss = false;
  // MakeSigWS(w, fileBaseName+"_8TeV",mass, isGauss);     
  //  cout << endl; cout << "Now make background workspace" << endl;
   MakeBkgWS(w, fileBkgName, mass, "Expol");    
   // MakeBkgWS(w, fileBkgName, mass, "ExpolPL");    
   MakeBkgWS(w, fileBkgName, mass, "DiJet"); 
   MakeBkgWS(w, fileBkgName, mass, "Lau"); 
   MakeBkgWS(w, fileBkgName, mass, "ExpPAR"); 
   //MakeBkgWS(w, fileBkgName, mass, "DiJetPL"); 
   //MakeBkgWS(w, fileBkgName, mass, "DiJetEXP"); 
   //MakeBkgWS(w, fileBkgName, mass, "DiJetEXPOL"); 
 
   
   

  std::cout<<"Now prepare WS for Bias study"<<std::endl;
  RooRealVar RooMINmass("RooMINmass", "RooMINmass", MINmass); 
  RooRealVar RooMAXmass("RooMAXmass", "RooMAXmass", MAXmass); 
  w->import(RooMINmass);
  w->import(RooMAXmass);

  //
  MakeWSBiasStudy(w,mass, perc);
  //cout << endl; cout << "Now prepare datacards" << endl;
  //for (int c=0; c<NCAT; c++) MakeDataCard_1Channel(w, fileBaseName, fileBkgName,c);
  
  return;
}

// Signal Data Set
void AddSigData(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  w->import(lumi); 
  
  // Variables
  RooArgSet* ntplVars = defineVariables();
  
  int iMass = abs(mass);      
  TFile sigFile1("histograms_CMS-HGG_19032013.root");   //ggh prod mode tree livia
  
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/ggh_m%d_8TeV", iMass));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/vbf_m%d_8TeV", iMass));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/wzh_m%d_8TeV", iMass));
  sigTree1->Add(TString::Format("histograms_CMS-HGG_19032013.root/tth_m%d_8TeV", iMass));
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");


  // common preselection cut
  TString mainCut = TString::Format("PhotonsMass>=130 && PhotonsMass<=900", mass, mass);   // livia
  
  
  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "sigWeighted" << endl;
  sigWeighted.Print("v");
  cout << "---- nX:  " << sigWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
  RooDataSet* signal[NCAT];
  for (int c=0; c<ncat; ++c) {

    // 0) chiara: 1cat only
    // signal[c] =  (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut);   //chiara, for 1 cat only

   
    // 1)  prime 4 cat livia
    if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut+TString::Format("&& (abs(ph1_eta)>1.566 || abs(ph2_eta)>1.566) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));

    w->import(*signal[c],Rename(TString::Format("SigWeight_cat%d",c)));
    
    cout << "cat " << c << ", signal[c]: " << endl;
    signal[c]->Print("v");
    cout << "---- for category " << c << ", nX for signal[c]:  " << signal[c]->sumEntries() << endl; 
    cout << endl;
  }

  // Create full weighted signal data set without categorization
  RooDataSet* signalAll = (RooDataSet*) sigWeighted.reduce(*w->var("PhotonsMass"),mainCut);
  w->import(*signalAll, Rename("SigWeight"));
  cout << "now signalAll" << endl;
  signalAll->Print("v");
  cout << "---- nX for signalAll:  " << signalAll->sumEntries() << endl; 
  cout << endl;
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


// Fit signal with model gauss pdfs
void SigModelFitGauss(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  RooAbsPdf* PhotonsMassSigGauss[NCAT];
  RooExtendPdf* PhotonsMassSigGaussExt[NCAT];

  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8),maxMassFit(mass*1.2); 


  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));


    RooRealVar* PhotonsMass = w->var("PhotonsMass");  
    PhotonsMass->setUnit("GeV");

    //per fissare la frac
    // ((RooRealVar*) w->var(TString::Format("PhotonsMass_sig_frac_cat%d",c)))->setConstant(true);

    RooFormulaVar m0(TString::Format("m0_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_m0_cat%d",c)));
    RooFormulaVar m1(TString::Format("m1_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_m0_cat%d",c)));
    RooFormulaVar sigma0(TString::Format("sigma0_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_sigma0_cat%d",c)));
    RooFormulaVar sigma1(TString::Format("sigma1_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_sigma1_cat%d",c)));
    RooFormulaVar frac(TString::Format("frac_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_frac_cat%d",c)));

    RooGaussian PhotonsMassGaussSig0(TString::Format("PhotonsMassGaussSig0_cat%d",c), "first gaussian component", *PhotonsMass, m0, sigma0);
    RooGaussian PhotonsMassGaussSig1(TString::Format("PhotonsMassGaussSig1_cat%d",c), "second gaussian component", *PhotonsMass,m0, sigma1);

    

    RooAddPdf* PhotonsMassSigAdd = new RooAddPdf(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c),RooArgList(PhotonsMassGaussSig0, PhotonsMassGaussSig1),frac) ;

    PhotonsMassSigGauss[c] =  (RooAbsPdf*)PhotonsMassSigAdd;
  
    
    //extended pdf
    PhotonsMass->setRange("signal range", minMassFit, maxMassFit);
    RooFormulaVar nsig(TString::Format("nsig_cat%d",c),"","@0",*w->var(TString::Format("nsig_cat%d",c)));
    PhotonsMassSigGaussExt[c] = new RooExtendPdf(TString::Format("PhotonsMassSigExt_cat%d",c),TString::Format("PhotonsMassSigExt_cat%d",c), *PhotonsMassSigGauss[c], nsig, "signal range" );
    

    RooFitResult* fitresults_sig = PhotonsMassSigGauss[c] ->fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"),SumW2Error(kTRUE), Save(kTRUE));
   
    std::cout<<TString::Format("******************************** Signal Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults_sig->Print("V");

    w->import(*PhotonsMassSigGauss[c]);
 
    
    // Plot to verify everything is ok
    RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
    sigToFit[c]->plotOn(plotPhotonsMassAll);
    PhotonsMassSigGauss[c]->plotOn(plotPhotonsMassAll);
    PhotonsMassSigGauss[c]  ->plotOn(plotPhotonsMassAll,Components(TString::Format("PhotonsMassGaussSig0_cat%d",c)),LineStyle(kDashed),LineColor(kGreen));
    PhotonsMassSigGauss[c]  ->plotOn(plotPhotonsMassAll,Components(TString::Format("PhotonsMassGaussSig1_cat%d",c)),LineStyle(kDashed),LineColor(kRed));
  
   
   
    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    c1->cd(1);
    plotPhotonsMassAll->Draw(); 
 

    TLegend *legmc = new TLegend(0.5791457,0.75,0.8291457,0.9340659, TString::Format("Category %d",c), "brNDC");
    legmc->AddEntry(plotPhotonsMassAll->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotPhotonsMassAll->getObject(1),"Parametric Model","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(3),"Gaussian 2 component","L");
    legmc->AddEntry(plotPhotonsMassAll->getObject(2),"Gaussian 1 component","L");
    legmc->SetTextSize(0.0206044);
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();
    
    TPaveText* label_cms = get_labelCMS(0, "2012", true);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");

    int massI(mass);
    c1->SaveAs("preliminaryPlots/prelimSignalGauss"+TString::Format("_M%d_cat%d.png",massI,c));
    c1->SaveAs("preliminaryPlots/prelimSignalGauss"+TString::Format("_M%d_cat%d.root",massI,c));
    c1->SetLogy();
    c1->SaveAs("preliminaryPlots/prelimSignalGauss"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
    c1->SaveAs("preliminaryPlots/prelimSignalGauss"+TString::Format("_M%d_cat%d_LOG.root",massI,c));

    
    // IMPORTANT: fix all pdf parameters to constant
     w->defineSet(TString::Format("SigPdfParam_cat%d",c), RooArgSet(*w->var("PhotonsMass"+TString::Format("_sig_m0_cat%d",c)),
    								   *w->var("PhotonsMass"+TString::Format("_sig_sigma0_cat%d",c)),
								    *w->var("PhotonsMass"+TString::Format("_sig_m1_cat%d",c)),
								    *w->var("PhotonsMass"+TString::Format("_sig_sigma1_cat%d",c)),
		                                                     *w->var("PhotonsMass"+TString::Format("_sig_frac_cat%d",c))));
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
  }
}


//********************//
//Fit function variying the width:


void SigModelResponseFcnFit(RooWorkspace* w, Float_t mass) {


  TFile sigFile1("histograms_CMS-HGG_19032013.root");   
  
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();

  sigTree1->Add("histograms_CMS-HGG_19032013.root/ggh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_19032013.root/vbf_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_19032013.root/wzh_m250_8TeV");
  sigTree1->Add("histograms_CMS-HGG_19032013.root/tth_m250_8TeV");
  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");

  // Variables
  RooArgSet* ntplVars = defineVariablesM250();
  TString mainCut1 = TString::Format("PhotonsMass > 130 && PhotonsMass<1000");   // livia
  RooDataSet sigWeighted("sigWeighted","dataset",sigTree1,*ntplVars,mainCut1,"evweight");
 
  // Variables
  RooRealVar* PhotonsMass = w->var("PhotonsMass");
  
 
  RooRealVar *mH = new RooRealVar("MH", "MH", MINmass, MAXmass);
  mH->setVal(mass);
  mH->setConstant();
  w->import(*mH);

  // w->Print("V");

  RooFormulaVar *massReduced_formula     = new RooFormulaVar("massReduced_formula","","@0/250 -1",*w->var("PhotonsMass"));
  RooRealVar* massReduced = (RooRealVar*) sigWeighted.addColumn(*massReduced_formula);
  massReduced->SetName("massReduced");
  massReduced->SetTitle("massReduced");
  w->import(*massReduced);  
  massReduced->setRange(-0.5, 0.5);
 // common preselection cut
  TString mainCut = TString::Format("massReduced>-0.5 && massReduced <0.5");   // livia
  
  
  RooDataSet* signal[NCAT];
  RooCBShape* ResponseCBpos[NCAT];
  RooCBShape* ResponseCBneg[NCAT];
  RooGaussian* ResponseGauss[NCAT];
  RooAddPdf* ResponseAddGauss[NCAT];
  RooAddPdf* ResponseAdd[NCAT];
  
  TCanvas* c1 = new TCanvas("c1", "c1", 1);
  c1->cd();
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  
  for(int c = 0; c<NCAT; c++){
    
    TLatex *lat  = new TLatex(0.65,0.9,TString::Format("Cat: %d", c));  
    lat->SetTextSize(0.038);
    lat->SetTextAlign(11);
    lat->SetTextFont(42); 
    lat->SetNDC();

   // 1)  prime 4 cat livia
   if (c==0) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==1) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)<1.5 && abs(ph2_eta)<1.5) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
   if (c==2) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
   if (c==3) signal[c] = (RooDataSet*) sigWeighted.reduce(*w->var("massReduced"),mainCut+TString::Format("&& (abs(ph1_eta)>1.5 || abs(ph2_eta)>1.5) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));

  
   //add cb neg +pos

   //cb pos                                                                                                                     
   RooFormulaVar CBpos_mean(TString::Format("ReducedMass_CBpos_sig_mean_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
   RooFormulaVar CBpos_sigma(TString::Format("ReducedMass_CBpos_sig_sigma_cat%d",c), "", "@0", *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)));
   RooFormulaVar CBpos_alphaCB(TString::Format("ReducedMass_CBpos_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)));
   RooFormulaVar CBpos_n(TString::Format("ReducedMass_CBpos_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Npos_cat%d",c)));
     
   //cb neg
   RooFormulaVar CBneg_n(TString::Format("ReducedMass_CBneg_sig_n_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_Nneg_cat%d",c)));
   RooFormulaVar CBneg_alphaCB(TString::Format("ReducedMass_CBneg_sig_alphaCB_cat%d",c),"", "@0", *w->var( TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)));
   
   ResponseCBpos[c] =  new RooCBShape(TString::Format("ResponseCBpos_cat%d",c),TString::Format("ResponseCBpos_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBpos_alphaCB, CBpos_n) ;
   

   ResponseCBneg[c] =  new RooCBShape(TString::Format("ResponseCBneg_cat%d",c),TString::Format("ResponseCBneg_cat%d",c) , *massReduced, CBpos_mean, CBpos_sigma,  CBneg_alphaCB, CBneg_n) ;
   

   
   RooFormulaVar CB_frac(TString::Format("ReducedMass_CBpos_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)));
   w->import(CB_frac);  
   ResponseAdd[c]= new RooAddPdf(TString::Format("ResponseAddPdf_cat%d",c),TString::Format("ResponseAddPdf_cat%d",c) , RooArgList(*ResponseCBneg[c], *ResponseCBpos[c]), CB_frac);
   w->import(*ResponseAdd[c]);
   


   //gauss+cb neg
   RooFormulaVar Gauss_frac(TString::Format("ReducedMass_Gauss_sig_frac_cat%d",c),"", "@0", *w->var(TString::Format("ReducedMass_sig_fracGauss_cat%d",c)));
   w->import(Gauss_frac);
   RooFormulaVar Gauss_sigma(TString::Format("ReducedMass_Gauss_sig_sigma_cat%d",c), "", "@0", *w->var(TString::Format("ReducedMass_sig_sigmaGauss_cat%d",c)));
 
   ResponseGauss[c] = new RooGaussian(TString::Format("ResponseGauss_cat%d",c), TString::Format("ResponseGauss_cat%d",c), *massReduced,CBpos_mean , Gauss_sigma);
   
   ResponseAddGauss[c]= new RooAddPdf(TString::Format("ResponseAddGaussPdf_cat%d",c), TString::Format("ResponseAddGaussPdf_cat%d",c), RooArgList(*ResponseCBneg[c], *ResponseGauss[c]), Gauss_frac);
   //w->import(*ResponseAddGauss[c]);
   
  
   // RooFitResult* fitresults = (RooFitResult* ) ResponseAddGauss[c]->fitTo(*signal[c],SumW2Error(kTRUE),  RooFit::Save(kTRUE));
   RooFitResult* fitresults = (RooFitResult* ) ResponseAdd[c]->fitTo(*signal[c],SumW2Error(kTRUE),Range(-1, 1),  RooFit::Save(kTRUE));
   //  std::cout<<TString::Format("******************************** Signal Fit results CB+Gauss  mass %f cat %d***********************************", mass, c)<<std::endl;
   std::cout<<TString::Format("******************************** Signal Fit results CB+CB  mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults->Print("V");
   

    RooPlot* plotG = massReduced->frame(Range(-1., 1.),Title("Mass Reduced"), Bins(60));
    signal[c]->plotOn(plotG);
   
    // ResponseAddGauss[c]->plotOn(plotG, LineColor(kBlue));
    ResponseAdd[c]->plotOn(plotG, LineColor(kBlue));
    // ResponseAddGauss[c]->plotOn(plotG,Components(TString::Format("ResponseGauss_cat%d",c)), LineColor(kGreen), LineStyle(kDashed));
    // ResponseAddGauss[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBneg_cat%d",c)), LineColor(kGreen), LineStyle(kDashed));
    ResponseAdd[c]->plotOn(plotG,Components(TString::Format("ResponseCBpos_cat%d",c)), LineColor(kRed), LineStyle(kDashed));
  
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*10 );
    plotG->GetXaxis()->SetTitle("#Delta m");
    plotG->GetXaxis()->SetTitleFont(42);
    plotG->GetXaxis()->SetTitleSize(0.05);
    TLegend* legmc = new TLegend(0.6, 0.6, 0.85, 0.89, "", "brNDC");
    legmc->SetTextSize(0.0206044);  
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    
    legmc->AddEntry(plotG->getObject(0),"m_{#gamma#gamma} = 250 GeV","LPE");    
    
    // legmc->AddEntry(plotG->getObject(1),"Sum of CB and Gauss","L");
    legmc->AddEntry(plotG->getObject(1),"Sum of two CB ","L");
    //legmc->AddEntry(plotG->getObject(2),"Gauss","L");    
    legmc->AddEntry(plotG->getObject(2),"CB 1","L");   
    legmc->AddEntry(plotG->getObject(3),"CB 2","L");   
    plotG->Draw();
    
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    
    
    c1->SetLogy();
    
    // c1->SaveAs(TString::Format("plots/responseFcnFitCBGauss_cat%d_LOG.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.pdf",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d_LOG.eps",c)); 
    
    
    plotG->GetYaxis()->SetRangeUser(0.0001,plotG->GetMaximum()*0.12 );
    lat->Draw("same");
    label_cms->Draw("same");
    label_sqrt->Draw("same"); 
    legmc->Draw("same");
    c1->SetLogy(0);  
    //  c1->SaveAs(TString::Format("plots/responseFcnFitCBGauss_cat%d.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d.png",c)); 
    c1->SaveAs(TString::Format("plots/responseFcnFitCBCB_cat%d.pdf",c)); 
    
    
    w->defineSet(TString::Format("ResponseAddPdfParam_cat%d",c),RooArgSet(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									  *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									  *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									  *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									  *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c))));
    SetConstantParams(w->set(TString::Format("ResponseAddPdfParam_cat%d",c)));
    
    
  

  }
  

 

}





// Fit signal with model gauss pdfs
void SigModelFitConvBW(RooWorkspace* w, Float_t mass, Int_t perc) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  
  
  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8);
  Float_t maxMassFit(mass*1.2); 
  RooRealVar* PhotonsMass = w->var("PhotonsMass"); 
  
  std::cout<<"----------------------------------------------------------------------------------------"<<std::endl;
 
  TPaveText* label_cms = get_labelCMS(0, "2012", true);
  TPaveText* label_sqrt = get_labelSqrt(0);
 
  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    
 
    //cb
    RooFormulaVar CBpos_mean_draw(TString::Format("CBpos_mean_draw_cat%d",c),"","@0+@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),*w->var("MH")));
   
    RooFormulaVar CBpos_mean(TString::Format("CBpos_mean_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)));
    RooFormulaVar CBpos_sigma(TString::Format("CBpos_sigma_cat%d",c),"","@0*@1",RooArgList(*w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)),*w->var("MH")) );
   
    RooFormulaVar CBpos_alphaCB(TString::Format("CBpos_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)) );
    RooFormulaVar CBneg_alphaCB(TString::Format("CBneg_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)) );
    RooFormulaVar CBpos_n(TString::Format("CBpos_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)) );
    RooFormulaVar CBneg_n(TString::Format("CBneg_n_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)) );
    RooFormulaVar CBpos_frac(TString::Format("CBpos_frac_cat%d",c),"","@0",*w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)) );
    
    RooCBShape ResCBpos_draw(TString::Format("ResCBpos_draw_cat%d",c),TString::Format("ResCBpos_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg_draw(TString::Format("ResCBneg_draw_cat%d",c),TString::Format("ResCBneg_draw_cat%d",c) , *PhotonsMass, CBpos_mean_draw, CBpos_sigma, CBneg_alphaCB, CBneg_n) ;
    RooCBShape ResCBpos(TString::Format("ResCBpos_cat%d",c),TString::Format("ResCBpos_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBpos_alphaCB, CBpos_n) ;
    RooCBShape ResCBneg(TString::Format("ResCBneg_cat%d",c),TString::Format("ResCBneg_cat%d",c) , *PhotonsMass, CBpos_mean, CBpos_sigma,CBneg_alphaCB, CBneg_n) ;
    PhotonsMass->setBins(40000, "cache");  
    //add CB pos + CB neg
    RooAddPdf ResAddPdf(TString::Format("ResAddPdf_cat%d",c),TString::Format("ResAddPdf_cat%d",c) , RooArgList(ResCBneg, ResCBpos), CBpos_frac);
    RooAddPdf ResAddPdf_draw(TString::Format("ResAddPdf_draw_cat%d",c),TString::Format("ResAddPdf_draw_cat%d",c) , RooArgList(ResCBneg_draw, ResCBpos_draw), CBpos_frac);


   
    //BW

    RooFormulaVar meanBW(TString::Format("meanBW_cat%d",c),"","@0",*w->var("MH"));  
    if(perc==0) width=0.1;
    if(perc==5) width=(w->var("MH")->getVal())*0.05;
    if(perc==10) width=(w->var("MH")->getVal())*0.10;
    if(perc==20) width=(w->var("MH")->getVal())*0.20;
    RooRealVar sigmaBW_var(TString::Format("sigmaBW_var_cat%d",c), TString::Format("sigmaBW_var_cat%d",c), width);
    std::cout<<" width:--------> "<<width<<std::endl;
    sigmaBW_var.setConstant();
    w->import(sigmaBW_var);
    
    
    RooFormulaVar sigmaBW(TString::Format("sigmaBW_cat%d",c),"","@0",*w->var(TString::Format("sigmaBW_var_cat%d",c)));     
    RooBreitWigner SigModelBW(TString::Format("SigModelBW_cat%d",c),TString::Format("SigModelBW_cat%d",c), *PhotonsMass, meanBW, sigmaBW);

  
    //CONV 
    RooFFTConvPdf  ConvolutedRes_CB(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c), *PhotonsMass,SigModelBW, ResAddPdf);
     w->import(ConvolutedRes_CB);
    std::cout<<".............> "<<c<<std::endl;
    RooHistFunc* rooFunc_norm = getRooHistFunc(c, w->var("MH"));
    w->import(*rooFunc_norm);
    //   w->Print("V");

    if(width < 2.&& mass <400 && mass!=350.){ //if i want to plot the fit
      sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
 
      //   RooFitResult* fitresults_CB = (RooFitResult* ) ConvolutedRes_CB.fitTo(*sigToFit[c],Range(MINmass,MAXmass), RooFit::Save(kTRUE));
      std::cout<<TString::Format("******************************** Signal Fit results CB mass %f cat %d***********************************", mass, c)<<std::endl;
      //fitresults_CB->Print("V");
      
      RooPlot* plotOnlyResPdf = w->var("PhotonsMass")->frame(Range(MINmass,MAXmass),Bins(100));
      sigToFit[c]->plotOn(plotOnlyResPdf, LineColor(kRed), LineStyle(kDashed));
      double max = plotOnlyResPdf->GetMaximum();
      
      // Plot to verify everything is ok
      RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(MINmass,MAXmass),Bins(100));
      sigToFit[c]->plotOn(plotPhotonsMassAll);
      SigModelBW.plotOn(plotPhotonsMassAll, LineColor(kGreen), LineStyle(kDashed));
      ConvolutedRes_CB.plotOn(plotPhotonsMassAll, LineColor(kBlue));

      
      TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
      c1->cd(1);
      
      plotPhotonsMassAll->Draw();  
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0001, max*1.2);
      TLegend *legmc = new TLegend(0.5491457,0.75,0.801457,0.9340659, TString::Format("Category %d",c), "brNDC");
      legmc->AddEntry(plotPhotonsMassAll->getObject(0),"Simulation","LPE");
      legmc->AddEntry(plotPhotonsMassAll->getObject(2),"BW #otimes Resolution","L");
      legmc->SetTextSize(0.0206044);
      legmc->SetTextFont(42);
      legmc->SetBorderSize(0);
      legmc->SetFillStyle(0);
      legmc->Draw();
      
      label_cms->Draw("same");
      label_sqrt->Draw("same");
      

      
      int massI(mass);
      c1->SetLogy();
     
      plotPhotonsMassAll->GetXaxis()->SetTitle("m_{#gamma #gamma}[GeV]");
     
      c1->SetLogy(0);
    
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.png",massI, c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d.pdf",massI, c));
      
      c1->SetLogy();
      plotPhotonsMassAll->GetYaxis()->SetRangeUser(0.0001,max*10. );
     
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.pdf",massI,c));
      c1->SaveAs("plots/prelimSignalCBCconvBW"+TString::Format("_M%d_cat%d_LOG.eps",massI,c));

      
  


    }
    
    // IMPORTANT: fix all pdf parameters to constant
   
      w->defineSet(TString::Format("ConvolutedPdfParam_cat%d",c),RooArgSet( *w->var(TString::Format("ReducedMass_sig_sigma_cat%d",c)), 
									    *w->var(TString::Format("ReducedMass_sig_alphaCBpos_cat%d",c)),
									    *w->var(TString::Format("ReducedMass_sig_alphaCBneg_cat%d",c)),
									    *w->var(TString::Format("ReducedMass_sig_Npos_cat%d",c)),
									    *w->var(TString::Format("ReducedMass_sig_Nneg_cat%d",c)),	   
									    *w->var(TString::Format("ReducedMass_sig_frac_cat%d",c)),  
									    *w->var(TString::Format("ReducedMass_sig_mean_cat%d",c)),
									    
									    *w->var(TString::Format("sigmaBW_var_cat%d",c))));
      
      SetConstantParams(w->set(TString::Format("ConvolutedPdfParam_cat%d",c)));
    
      //   w->Print("V");
    
  }

}



//********************//








// Fit signal with model gauss pdfs
void SigModelFitCBC(RooWorkspace* w, Float_t mass) {

  Int_t ncat = NCAT;
  RooDataSet* sigToFit[NCAT];
  RooCBCrujffPdf* PhotonsMassSig[NCAT];
  RooExtendPdf* PhotonsMassSigExt[NCAT];

  Float_t MASS(mass);  
  Float_t minMassFit(mass*0.8),maxMassFit(mass*1.2); 


  // Fit Signal 
  for (int c=0; c<ncat; ++c) {
    cout << "---------- Category = " << c << endl;
    sigToFit[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));


    /* tentativo di fissare m1=m0*/
    RooRealVar* PhotonsMass = w->var("PhotonsMass");  
    PhotonsMass->setUnit("GeV");

   
    RooFormulaVar CBC_mean(TString::Format("CBC_mean_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_mean_cat%d",c)) );
    RooFormulaVar CBC_sigma(TString::Format("CBC_sigma_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_sigma_cat%d",c)) );
    RooFormulaVar CBC_alphaC(TString::Format("CBC_alphaC_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_alphaC_cat%d",c)) );
    RooFormulaVar CBC_alphaCB(TString::Format("CBC_alphaCB_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_alphaCB_cat%d",c)) );
    RooFormulaVar CBC_n(TString::Format("CBC_n_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_sig_n_cat%d",c)) );

    
    PhotonsMassSig[c] =  new RooCBCrujffPdf(TString::Format("PhotonsMassSig_cat%d",c),TString::Format("PhotonsMassSig_cat%d",c) , *PhotonsMass, CBC_mean, CBC_sigma, CBC_alphaC, CBC_alphaCB, CBC_n) ;

    //extended pdf
    PhotonsMass->setRange("signal range", minMassFit, maxMassFit);
    RooFormulaVar nCBC(TString::Format("nCBC_cat%d",c),"","@0",*w->var(TString::Format("nsigCBC_cat%d",c)));
    nCBC.Print();
    PhotonsMassSigExt[c] = new RooExtendPdf(TString::Format("PhotonsMassSigCBCExt_cat%d",c),TString::Format("PhotonsMassSigCBCExt_cat%d",c), *PhotonsMassSig[c], nCBC, "signal range" );
    

    //  w->import(*PhotonsMassSig[c]);
    w->import(*PhotonsMassSig[c]);
  
    //   w->Print();

    RooFitResult* fitresults_sig = PhotonsMassSig[c] ->fitTo(*sigToFit[c],Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));
    std::cout<<TString::Format("******************************** Signal Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresults_sig->Print("V");
    

    // Plot to verify everything is ok
    RooPlot* plotPhotonsMassAll = w->var("PhotonsMass")->frame(Range(minMassFit-20,maxMassFit+20),Bins(100));
    sigToFit[c]->plotOn(plotPhotonsMassAll);
    PhotonsMassSig[c]->plotOn(plotPhotonsMassAll);
     
   
    
    TCanvas* c1 = new TCanvas("c1","PhotonsMass",0,0,800,800);
    c1->cd(1);
    plotPhotonsMassAll->Draw();  
    
    TLegend *legmc = new TLegend(0.5491457,0.75,0.801457,0.9340659, TString::Format("Category %d",c), "brNDC");
    legmc->AddEntry(plotPhotonsMassAll->getObject(5),"Simulation","LPE");
    legmc->AddEntry(plotPhotonsMassAll->getObject(1),"Parametric Model CB Cruiff","L");
    legmc->SetTextSize(0.0206044);
    legmc->SetTextFont(42);
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();

 
    TPaveText* label_cms = get_labelCMS(0, "2012", true);
    TPaveText* label_sqrt = get_labelSqrt(0);
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    int massI(mass);
    c1->SaveAs("preliminaryPlots/prelimSignalCBC"+TString::Format("_M%d_cat%d.png",massI, c));
    c1->SaveAs("preliminaryPlots/prelimSignalCBC"+TString::Format("_M%d_cat%d.root",massI, c));

    c1->SetLogy();
    c1->SaveAs("preliminaryPlots/prelimSignalCBC"+TString::Format("_M%d_cat%d_LOG.png",massI,c));
    c1->SaveAs("preliminaryPlots/prelimSignalCBC"+TString::Format("_M%d_cat%d_LOG.root",massI,c));

    // IMPORTANT: fix all pdf parameters to constant
     w->defineSet(TString::Format("SigPdfParam_cat%d",c), RooArgSet(*w->var("PhotonsMass"+TString::Format("_sig_mean_cat%d",c)),
    								   *w->var("PhotonsMass"+TString::Format("_sig_sigma_cat%d",c)),
								    *w->var("PhotonsMass"+TString::Format("_sig_alphaC_cat%d",c)),
								    *w->var("PhotonsMass"+TString::Format("_sig_alphaCB_cat%d",c)),
                                                                    *w->var("PhotonsMass"+TString::Format("_sig_n_cat%d",c))));
    SetConstantParams(w->set(TString::Format("SigPdfParam_cat%d",c)));
  }
}







RooFitResult* BkgModelFitBernstein(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
  RooBernstein* PhotonsMassBkg_pdf[NCAT];
  RooFitResult* fitresult[NCAT];
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  

  Float_t minMassFit, maxMassFit;
  /* if(mass<200){ //m150
    minMassFit = 100;
    maxMassFit = 200;
  } else if(mass>=200 && mass <300){//m200 and m250
    minMassFit = 180;
    maxMassFit = 300;
    }else{//m300 m400*/
    minMassFit = 120;
    maxMassFit = 480;
    // }
 
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
    RooRealVar x("x", "x", mass, minMassFit, maxMassFit);

    RooRealVar *c0 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c0_cat%d",c));
    RooRealVar *c1 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c1_cat%d",c));
    RooRealVar *c2 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c2_cat%d",c));
    RooRealVar *c3 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c3_cat%d",c));
    RooRealVar *c4 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c4_cat%d",c));
    RooRealVar *c5 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c5_cat%d",c));
    RooRealVar *c6 = w->var(TString::Format("PhotonsMass_bkg_8TeV_c6_cat%d",c));


    PhotonsMassBkg_pdf[c] = new RooBernstein("PhotonsMassBkg"+TString::Format("_cat%d",c),"PhotonsMassBkg"+TString::Format("_cat%d",c), *w->var("PhotonsMass"), RooArgList( *c0, *c1, *c2, *c3, *c4, *c5, *c6)); //*c3


    fitresult[c] = PhotonsMassBkg_pdf[c].fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkg_pdf[c]);
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   
  
    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(200);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(nBinsMass);
    data[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    
    PhotonsMassBkg_pdf[c].plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 

    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = 227;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;

    
    if( blind ) {
      
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >400.4");
      
      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 
      
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
    } 
    

  plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
  plotPhotonsMassBkg[c]->Draw(); 

  TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
  legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
  legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: Bernstein","L");
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  // legdata->SetTextAlign(31);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw();
  
  TPaveText* label_cms = get_labelCMS(0, "2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  
  
    //write down the chi2 of the fit on the


    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678,"brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");
    
    int massI(mass);  
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_BERN_M%d.png",c, massI));
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_BERN_M%d.root",c, massI));

  plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
  ctmp->SetLogy();

  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_BERN_LOG_M%d.png",c, massI));
  ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_BERN_LOG_M%d.root",c, massI));

  dobands = false;
  //********************************************************************************//
  if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkg_pdf[c];
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
  }
  return fitresult;
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

    RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_ExpPAR_truth_cat%d",c), "exp(-@1*@0)*pow(@0, @2)", RooArgList(*PhotonsMass, *p1mod, *p2mod));
    
    fitresult[c] = PhotonsMassBkg->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
    w->import(*PhotonsMassBkg);
 

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkg->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
   
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-2;
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
      
    TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678,"brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");
    
 
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





RooFitResult* BkgModelFitLauFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;

  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];

  RooPlot* plotPhotonsMassBkg[NCAT];

  Float_t minMassFit, maxMassFit;
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;
  

    std::cout<<"---------------->>>>> "<<minMassFit <<"    "<< maxMassFit<<std::endl;
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
   
    // fit con expol o2  
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1Lau_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2Lau_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau2_cat%d",c)));   
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3Lau_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_lau3_cat%d",c))); 

    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xLau_cat%d",c),"","@0",*w->var("PhotonsMass"));

    //    RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_Lau_truth_cat%d",c), "(1-@1)*pow(@0,-5.0)+@1*pow(@0,-6.0)", RooArgList(*x, *p1mod));
     RooAbsPdf* PhotonsMassBkg = new RooGenericPdf(TString::Format("PhotonsMassBkg_Lau_truth_cat%d",c), "(1-@1)*pow(@0,-4.0)+@1*pow(@0,-5.0)", RooArgList(*PhotonsMass, *p1mod));
    
    
    fitresult[c] = PhotonsMassBkg->fitTo(*data[c], Range(minMassFit,maxMassFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   
    w->import(*PhotonsMassBkg);
    std::cout<<TString::Format("******************************** Background Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");



    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background ",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkg->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit)); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-1;
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
        TPaveText* label_chi2 = new TPaveText(0.5524194,0.6419492,0.796371,0.7690678,"brNDC");
    label_chi2->SetFillColor(kWhite);
    label_chi2->SetTextSize(0.035);
    label_chi2->SetTextFont(42);
    label_chi2->SetTextAlign(31); // align right
    label_chi2->AddText(TString::Format("Fit chi square/dof = %.3f", chi2));
    label_chi2->AddText(TString::Format("Chi square Prob = %.3f", prob));
    label_chi2->Draw("same");
  
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_LAU_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_LAU_M%d.pdf",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_LAU_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_LAU_LOG_M%d.pdf",c,massI));

  }

 

  return fitresult;
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






RooFitResult* BkgModelFitDiJetPLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

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
 
    minMassFit = MINmass;
    maxMassFit = MAXmass;

  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
  PhotonsMass->Print("V");

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
  
    // fit a la dijets
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_2_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_2_cat%d",c)));   
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_2_cat%d",c)));
    RooFormulaVar *pow = new RooFormulaVar(TString::Format("powDiJetPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_pow1DiJetPL_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetPL_cat%d",c)));
   

   
    RooAbsPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DIJET_truth_cat%d",c),"pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
    RooAbsPdf* PhotonsMassBkgTmp0PL = new RooGenericPdf(TString::Format("PhotonsMassBkg_PL_truth_cat%d",c), "pow(@0, @1)", RooArgList(*w->var("PhotonsMass"),  *pow));
    
    //rooAdd
    RooAbsPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_DiJetPL_truth_cat%d",c),TString::Format("PhotonsMassBkg_DiJetPL_truth_cat%d",c), RooArgList(*PhotonsMassBkgTmp0DiJet,*PhotonsMassBkgTmp0PL), RooArgList(*pFrac1));

      
    fitresult[c] = PhotonsMassBkgTmpAdd->fitTo(*data[c], RooFit::FitOptions("MHTR"),  Save(kTRUE));
    w->import(*PhotonsMassBkgTmpAdd);
   
    std::cout<<TString::Format("******************************** Background DiJetPL Fit results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    

    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_PL_truth_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_DIJET_truth_cat%d",c)),LineColor(kOrange),LineStyle(kDashed)
); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-5;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass")," PhotonsMass >850");
      TH1F* h_up= new TH1F("h_up", "h_up",nBinsMass, 130, 1000);
      h_up->Sumw2();
      data_up->fillHistogram(h_up, RooArgList(*PhotonsMass));
      TH1F* h_down= new TH1F("h_down", "h_down",nBinsMass, 130, 1000);
      h_down->Sumw2();
      data_down->fillHistogram(h_down, RooArgList(*PhotonsMass));
  
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->GetYaxis()->SetRangeUser(0.1,10000);
    plotPhotonsMassBkg[c]->Draw();  
     if( blind ) {
       h_up->Draw("sameP");
       h_down->Draw("sameP");
     }
    plotPhotonsMassBkg[c]->GetYaxis()->SetRangeUser(0.1,10000);

    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJetPL","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Parametric Model: PL","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(3),"Parametric Model: DiJet","L");
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
    /*  if (dobands) {

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
      }*/

    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETPL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETPL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETPL_M%d.root",c,massI));

    ctmp->SetLogy();
    // plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETPL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETPL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETPL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}






RooFitResult* BkgModelFitDiJetEXPFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooFitResult* fitresult2[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
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
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_3_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_3_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_3_cat%d",c)));
    RooFormulaVar *exp1 = new RooFormulaVar(TString::Format("expDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_exp1DiJetEXP_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXP_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetEXP_cat%d",c)));
   
   
    RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DIJETE_truth_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
  
    RooExponential* PhotonsMassBkgTmp0Exp = new RooExponential(TString::Format("PhotonsMassBkg_EXP_truth_cat%d",c),"", *w->var("PhotonsMass"),  *exp1);
    
   
    RooAddPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_DiJetEXP_truth_cat%d",c),TString::Format("PhotonsMassBkg_DiJetEXP_truth_cat%d",c) , RooArgList(*PhotonsMassBkgTmp0DiJet, *PhotonsMassBkgTmp0Exp), RooArgList(*pFrac1));
    
    fitresult[c] = PhotonsMassBkgTmpAdd->fitTo(*data[c],RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*PhotonsMassBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXP Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_EXP_truth_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_DIJETE_truth_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-5;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
   
    if( blind ) {
 
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass")," PhotonsMass >850");
      TH1F* h_up= new TH1F("h_up", "h_up",nBinsMass, 130, 1000);
      h_up->Sumw2();
      data_up->fillHistogram(h_up, RooArgList(*PhotonsMass));
      TH1F* h_down= new TH1F("h_down", "h_down",nBinsMass, 130, 1000);
      h_down->Sumw2();
      data_down->fillHistogram(h_down, RooArgList(*PhotonsMass));
   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.1,10000,"Y");
    plotPhotonsMassBkg[c]->Draw();  
   if( blind ) {
       h_up->Draw("sameP");
       h_down->Draw("sameP");
     }
  
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJetEXP","L");  
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Parametric Model: EXP","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(3),"Parametric Model: DiJet","L");
  
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
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_M%d.root",c,massI));

    ctmp->SetLogy();
    //  plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXP_LOG_M%d.root",c,massI));

  }



  return fitresult;
}







RooFitResult* BkgModelFitDiJetEXPOLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
 
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;

  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal
  RooDataSet* signal[NCAT];

  RooAbsPdf* PhotonsMassSig[NCAT];
  
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
    
  
    // fit a la dijets
    
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_4_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_4_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_4_cat%d",c)));
    RooFormulaVar *expol1 = new RooFormulaVar(TString::Format("expol1DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_4_cat%d",c)));
    RooFormulaVar *expol2 = new RooFormulaVar(TString::Format("expol2DiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_4_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracDiJetEXPOL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracDiJetEXPOL_cat%d",c)));
   
   
    RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf(TString::Format("PhotonsMassBkg_DiJetEx_truth_cat%d",c), "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod,*p3mod));
    RooGenericPdf* PhotonsMassBkgTmp0Expol = new RooGenericPdf(TString::Format("PhotonsMassBkg_ExpolDiJ_truth_cat%d",c), "exp(-@0/(@1+@2*@0))", RooArgList(*w->var("PhotonsMass"), *expol1, *expol2));
    
   
    RooAddPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_DiJetEXPOL_truth_cat%d",c),TString::Format("PhotonsMassBkg_DiJetEXPOL_truth_cat%d",c) , RooArgList(*PhotonsMassBkgTmp0DiJet, *PhotonsMassBkgTmp0Expol), RooArgList(*pFrac1));
    
    fitresult[c] = PhotonsMassBkgTmpAdd->fitTo(*data[c],RooFit::FitOptions("MHTR"), Save(kTRUE));//RooFit::FitOptions("MHTER"), Range(minMassFit,maxMassFit),    
    w->import(*PhotonsMassBkgTmpAdd);
  
    std::cout<<TString::Format("******************************** Background DiJetEXPOL Fit results mass %f cat %d ***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
   

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(60);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_ExpolDiJ_truth_cat%d",c)),LineColor(kViolet),LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_DiJetEx_truth_cat%d",c)),LineColor(kOrange),LineStyle(kDashed));   

    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-6;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
   

if( blind ) {
 
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 173.5");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass")," PhotonsMass >850");
      TH1F* h_up= new TH1F("h_up", "h_up",nBinsMass, 130, 1000);
      h_up->Sumw2();
      data_up->fillHistogram(h_up, RooArgList(*PhotonsMass));
      TH1F* h_down= new TH1F("h_down", "h_down",nBinsMass, 130, 1000);
      h_down->Sumw2();
      data_down->fillHistogram(h_down, RooArgList(*PhotonsMass));
  

   
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
      } 
       
   
    plotPhotonsMassBkg[c]->GetXaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
    plotPhotonsMassBkg[c]->SetAxisRange(0.1,10000,"Y");
    plotPhotonsMassBkg[c]->Draw();  
    if( blind ) {
       h_up->Draw("sameP");
       h_down->Draw("sameP");
     }
  
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(4),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: DiJetEXPOL","L");  
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Parametric Model: Expol","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(3),"Parametric Model: DiJet","L");
  
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
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    //  plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_DIJETEXPOL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}






RooFitResult* BkgModelFitExpolFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

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
    
  
    // fit con expo pol
   
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_cat%d",c)));
    RooFormulaVar *p3mod = new RooFormulaVar(TString::Format("par3Expol_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol3_cat%d",c)));
 
    //  ((RooRealVar*) w->var(TString::Format("massggnewvtx_bkg_8TeV_slope3_cat%d",c)))->setConstant(true);
    //cout << "---------------- Parameter 3 set to const" << endl;
    //   RooFormulaVar *sqrtS = new RooFormulaVar(TString::Format("sqrtS_cat%d",c),"","@0",*w->var(TString::Format("sqrtS_cat%d",c)));
   
    RooFormulaVar *x     = new RooFormulaVar(TString::Format("xExpol_cat%d",c),"","@0",*w->var("PhotonsMass"));

   
    RooAbsPdf* PhotonsMassBkgTmp0 = new RooGenericPdf(TString::Format("PhotonsMassBkg_Expol_truth_cat%d",c), "exp(-@0/(@1+@2*@0))", RooArgList(*x, *p1mod, *p2mod));
   
    PhotonsMass->setRange("bkg range", MINmass, MAXmass);
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
   
    PhotonsMassBkgTmp0->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range(minMassFit,maxMassFit),NormRange("bkg range")); 
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = nBinsMass-2;
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
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: Expol","L");
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

    
    dobands = false;
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
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}


   



RooFitResult* BkgModelFitExpolPLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

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
    
  
    // fit con expo pol
   
    RooFormulaVar *p1mod = new RooFormulaVar(TString::Format("par1ExpolPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_2_cat%d",c)));
    RooFormulaVar *p2mod = new RooFormulaVar(TString::Format("par2ExpolPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_2_cat%d",c)));
    RooFormulaVar *pow1 = new RooFormulaVar(TString::Format("powExpolPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_pow1ExpolPL_cat%d",c)));
    RooFormulaVar *pFrac1 = new RooFormulaVar(TString::Format("fracExpolPL_cat%d",c),"","@0",*w->var(TString::Format("PhotonsMass_bkg_8TeV_fracExpolPL_cat%d",c)));
      
    RooAbsPdf* PhotonsMassBkgTmp0Expol = new RooGenericPdf(TString::Format("PhotonsMassBkg_EXPOL_truth_cat%d",c), "exp(-@0/8000./(@1+@2*@0/8000.))", RooArgList(*w->var("PhotonsMass"), *p1mod, *p2mod ));
    
    RooAbsPdf* PhotonsMassBkgTmp0PLex = new RooGenericPdf(TString::Format("PhotonsMassBkg_PLex_truth_cat%d",c), "pow(@0,@1)", RooArgList(*w->var("PhotonsMass"),  *pow1));
   
    RooAbsPdf* PhotonsMassBkgTmpAdd = new RooAddPdf(TString::Format("PhotonsMassBkg_ExpolPL_truth_cat%d",c),TString::Format("PhotonsMassBkg_ExpolPL_truth_cat%d",c),RooArgList(*PhotonsMassBkgTmp0Expol, *PhotonsMassBkgTmp0PLex), RooArgList(*pFrac1));
    w->import(*PhotonsMassBkgTmpAdd);
    fitresult[c] = PhotonsMassBkgTmpAdd->fitTo(*data[c],RooFit::FitOptions("MHTR"),   Save(kTRUE));//
   
   

    //   w->Print("V");
    std::cout<<TString::Format("******************************** Background Fit EXPOLPL results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
    //   w->Print("V");

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(40);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_PL_truth_cat%d",c)),LineColor(kViolet), LineStyle(kDashed)); 
    PhotonsMassBkgTmpAdd->plotOn(plotPhotonsMassBkg[c],Components(TString::Format("PhotonsMassBkg_EXPOL_truth_cat%d",c)),LineColor(kOrange), LineStyle(kDashed)); 
    
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = 38;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

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
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: ExpolPL","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Parametric Model: PL","L");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(3),"Parametric Model: Expol","L");
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

    
    dobands = false;
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
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOLPL_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOLPL_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOLPL_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOLPL_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOLPL_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPOLPL_LOG_M%d.root",c,massI));

  }



  return fitresult;
}




RooFitResult* BkgModelFitExpPLFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {
  int nPowLow=1;
  Int_t ncat = NCAT;
  RooDataSet* data[NCAT];
  RooAddPdf* PhotonsMassBkg[NCAT];

  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  Float_t minMassFit, maxMassFit;

    minMassFit = 120;
    maxMassFit = 480;
    
    // Fit data with background pdf for data limit
    RooRealVar* PhotonsMass = w->var("PhotonsMass");  
    PhotonsMass->setUnit("GeV");


    
   
    RooGenericPdf* PhotonsMassBkgPowLow1[NCAT];
    RooGenericPdf* PhotonsMassBkgPowLow2[NCAT];
    RooGenericPdf* PhotonsMassBkgPowLow3[NCAT];
    RooGenericPdf* PhotonsMassBkgPowLow4[NCAT];
  
    for (int c=0; c<ncat; ++c) {
    cout << "---------- category = " << c << endl;
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
    RooRealVar* pow1 = new RooRealVar(TString::Format("PhotonsMass_pow1_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_pow1_cat%d_%dpowlow", c, nPowLow),0.8, 0.,50.);
    RooRealVar* pow2 = new RooRealVar(TString::Format("PhotonsMass_pow2_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_pow2_cat%d_%dpowlow", c, nPowLow),0.8, 0.,50.);
    RooRealVar* pow3 = new RooRealVar(TString::Format("PhotonsMass_pow3_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_pow3_cat%d_%dpowlow", c, nPowLow),0.8, 0.,50.);
    RooRealVar* pow4 = new RooRealVar(TString::Format("PhotonsMass_pow4_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_pow4_cat%d_%dpowlow", c, nPowLow),0.8, 0.,50.);

    RooRealVar* frac12 = new RooRealVar(TString::Format("PhotonsMass_frac12_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_frac12_cat%d_%dpowlow", c, nPowLow),0.1,0.,1.); 
    RooRealVar* frac23 = new RooRealVar(TString::Format("PhotonsMass_frac23_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_frac23_cat%d_%dpowlow", c, nPowLow),0.1,0.,1.); 
    RooRealVar* frac34 = new RooRealVar(TString::Format("PhotonsMass_frac34_cat%d_%dpowlow", c, nPowLow),TString::Format("PhotonsMass_frac34_cat%d_%dpowlow", c, nPowLow),0.1,0.,1.); 
 

   
    
    PhotonsMassBkgPowPlow1[c] = new RooGenericPdf(TString::Format("PhotonsMassBkgPowLow1_ExpPL_truth_cat%d",c),TString::Format("PhotonsMassBkgPowLow_ExpPL_truth_cat%d",c),"pow(@0,-@1)", RooArgList(*PhotonsMass, *pow1));
    PhotonsMassBkgPowLow2[c] = new RooGenericPdf(TString::Format("PhotonsMassBkgPowLow2_ExpPL_truth_cat%d",c,),TString::Format("PhotonsMassBkgPowLow2_ExpPL_truth_cat%d",c,) ,"exp(-@1*@0)", RooArgList(*PhotonsMass, *pow2));
    
    //  if( nPowLow==3)PhotonsMassBkgPowLow3[c] = new RooGenericPdf(TString::Format("PhotonsMassBkgPowLow3_cat%d_%dpowlow",c,nPowLow), TString::Format("PhotonsMassBkgPowLow3_cat%d_%dpowlow",c,nPowLow),"pow(@0,-@1)", RooArgList(*PhotonsMass, *pow3)); 
    // if( nPowLow==4)PhotonsMassBkgPowLow4[c] = new RooGenericPdf(TString::Format("PhotonsMassBkgPowLow4_cat%d_%dpowlow",c,nPowLow), TString::Format("PhotonsMassBkgPowLow4_cat%d_%dpowlow",c,nPowLow),"pow(@0,-@1)", RooArgList(*PhotonsMass, *pow4)); 
    
    
     PhotonsMassBkg[c] = new RooAddPdf(TString::Format("PhotonsMassBkg_ExpPL_truth_cat%d",c),TString::Format("PhotonsMassBkg_ExpPL_truth_cat%d",c),RooArgList(*PhotonsMassBkgPowLow1[c],*PhotonsMassBkgPowLow2[c]), RooArgList(*frac12) );

    // if( nPowLow==3) PhotonsMassBkg[c] = new RooAddPdf(TString::Format("PhotonsMassBkg_PowLow_truth_cat%d_3powlow",c),TString::Format("PhotonsMassBkg_PowLow_truth_cat%d_3powlow",c),RooArgList(*PhotonsMassBkgPowLow1[c],*PhotonsMassBkgPowLow2[c],*PhotonsMassBkgPowLow3[c]), RooArgList(*frac12,*frac23) );
    // if( nPowLow==4) PhotonsMassBkg[c] = new RooAddPdf(TString::Format("PhotonsMassBkg_PowLow_truth_cat%d_4powlow",c),TString::Format("PhotonsMassBkg_PowLow_truth_cat%d_4powlow",c),RooArgList(*PhotonsMassBkgPowLow1[c],*PhotonsMassBkgPowLow2[c],*PhotonsMassBkgPowLow3[c],*PhotonsMassBkgPowLow4[c]), RooArgList(*frac12,*frac23,*frac34) );
    


    fitresult[c] = PhotonsMassBkg[c]->fitTo(*data[c], Strategy(1),Minos(kFALSE), Range(minMassFit,maxMassFit),SumW2Error(kTRUE), Save(kTRUE));
    w->import(*PhotonsMassBkg[c]);
    std::cout<<TString::Format("******************************** Background Fit EXPPLresults mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");

    //   w->Print("v");

    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(40);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(nBinsMass);
    data[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kWhite),MarkerColor(kWhite));    
    
    PhotonsMassBkg[c].plotOn(plotPhotonsMassBkg[c],LineColor(kBlue),Range("fitrange"),NormRange("fitrange")); 

    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof;
    if( nPowLow==2) ndof= 37;
    if( nPowLow==3) ndof= 35;
    if( nPowLow==4) ndof= 33;

    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;

    blind=false;
    if( blind ) {
      
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402.");
      
      data_up->plotOn(plotPhotonsMassBkg[c]);    
      data_down->plotOn(plotPhotonsMassBkg[c]); 
      
    } else {
      data[c]->plotOn(plotPhotonsMassBkg[c]);    
    } 
    

    plotPhotonsMassBkg[c]->SetAxisRange(0.001,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    plotPhotonsMassBkg[c]->Draw(); 
    
    TLegend *legdata = new TLegend(0.3790323,0.7775424,0.6290323,0.9279661, TString::Format("Category %d",c), "brNDC");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(2),"Data","LPE");
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),TString::Format("Parametric Model:ExpPL", nPowLow),"L");
    legdata->SetTextSize(0.035);
    legdata->SetTextFont(42);
    // legdata->SetTextAlign(31);
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw();
    
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
    
    int massI(mass);
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPL_M%d.png",c, massI, nPowLow));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPL_M%d.root",c, massI, nPowLow));
    
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SetLogy();
    
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPL_LOG_M%d.png",c, massI, nPowLow));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_EXPPL_LOG_M%d.root",c, massI, nPowLow));
    
    dobands = false;


    //********************************************************************************//
    if (dobands) {

      RooAbsPdf *cpdf; cpdf = PhotonsMassBkg[c];
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
  }
  return fitresult;
}



void SetConstantParams(const RooArgSet* params) {

  cout << endl; cout << "Entering SetConstantParams" << endl;
  TIterator* iter(params->createIterator());
  for (TObject *a = iter->Next(); a != 0; a = iter->Next()) {
    RooRealVar *rrv = dynamic_cast<RooRealVar *>(a);
    if (rrv) { rrv->setConstant(true); std::cout << " " << rrv->GetName(); }
  }  
}





RooFitResult* BkgModelFitRooKeyFunc(RooWorkspace* w, Bool_t dobands, Float_t mass, bool blind) {

  Int_t ncat = NCAT;
  std::cout<<"isBlind: "<<blind<<std::endl;
  // retrieve pdfs and datasets from workspace to fit with pdf models
  RooDataSet* data[NCAT];
 
  RooFitResult* fitresult[NCAT];;
  RooPlot* plotPhotonsMassBkg[NCAT];

  // dobands and dosignal

  RooAbsPdf* PhotonsMassBkg[NCAT];
  
  Float_t minMassFit, maxMassFit;
 
  minMassFit = MINmass;
  maxMassFit = MAXmass;
  
  // Fit data with background pdf for data limit
  RooRealVar* PhotonsMass = w->var("PhotonsMass");  
  PhotonsMass->setUnit("GeV");
 
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
  

  TFile* f_roo = new TFile("BiasStudy/workspaces/HighMass-hgg.RooKeysPdfMCBkg_8TeV.root");
  RooWorkspace *wRoo = f_roo->Get("w_bias");
  TFile* f_roo123 = new TFile("BiasStudy/workspaces/HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat123.root");
  RooWorkspace *wRoo123 = f_roo123->Get("w_bias");
  
 
  for (int c = 0; c < ncat; ++c) {
    data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    
    if (c==0) PhotonsMassBkg[c] = (RooAbsPdf*) *wRoo->pdf("BkgMCKeyPdf_bw4_cat0");
    else if(c>0) PhotonsMassBkg[c] = (RooAbsPdf*) *wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c));
    PhotonsMassBkg[c]->SetTitle(TString::Format("PhotonsMassBkg_ROOKEY_truth_cat%d",c));
    PhotonsMassBkg[c]->SetName(TString::Format("PhotonsMassBkg_ROOKEY_truth_cat%d",c));
    w->import(*PhotonsMassBkg[c]);
    fitresult[c] = PhotonsMassBkg[c]->fitTo(*data[c],RooFit::FitOptions("MHTR"),RooFit::Range(MINmass, MAXmass),   Save(kTRUE));//
   
   

    //   w->Print("V");
    std::cout<<TString::Format("******************************** Background Fit EXPOLPL results mass %f cat %d***********************************", mass, c)<<std::endl;
    fitresult[c]->Print("V");
    //   w->Print("V");

    //************************************************
    // Plot PhotonsMass background fit results per categories 
    TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
    Int_t nBinsMass(40);
    plotPhotonsMassBkg[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
  
    data[c]->plotOn(plotPhotonsMassBkg[c],RooFit::Invisible());    
   
    PhotonsMassBkg[c]->plotOn(plotPhotonsMassBkg[c],LineColor(kBlue)); 
    
    double chi2 = plotPhotonsMassBkg[c]->chiSquare(3);
    Int_t ndof = 39;
    std::cout<<"------> "<< ndof<<std::endl;
    double prob = TMath::Prob(chi2*ndof, ndof);
    std::cout<<prob<<std::endl;
    blind=false;
    if( blind ) {
     
      RooDataSet* data_down = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass < 178.");
      RooDataSet* data_up = (RooDataSet*) data[c].reduce(*w->var("PhotonsMass"),"PhotonsMass >402");

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
    legdata->AddEntry(plotPhotonsMassBkg[c]->getObject(1),"Parametric Model: RooKey4","L");
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

    
    dobands = false;
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
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_ROOKEY_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_ROOKEY_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_ROOKEY_M%d.root",c,massI));

    ctmp->SetLogy();
    plotPhotonsMassBkg[c]->SetAxisRange(1.3,plotPhotonsMassBkg[c]->GetMaximum()*1.5,"Y");
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_ROOKEY_LOG_M%d.png",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_ROOKEY_LOG_M%d.pdf",c,massI));
    ctmp->SaveAs("preliminaryPlots/prelimBkg"+TString::Format("_cat%d_ROOKEY_LOG_M%d.root",c,massI));

  }



  return fitresult;
}











// Write signal pdfs and datasets into the workspace 
void MakeSigWS(RooWorkspace* w, const char* fileBaseName,float mass, bool isGauss){
  
  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");  

  //********************************//
  // Retrieve P.D.F.s

  for (int c=0; c<ncat; ++c) {
     
      wAll->import(*w->pdf("PhotonsMassSig"+TString::Format("_cat%d",c)));
    
      if(mass<400) wAll->import(*w->data(TString::Format("SigWeight_cat%d",c)));
 
  }
  std::cout << "done with importing signal pdfs" << std::endl;

  // (2) Systematics on energy scale and resolution // chiara: per ora tutte le sistematiche non hanno senso
  // wAll->factory("CMS_hgg_sig_m0_absShift[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_m0_absShift_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat0(massggnewvtx_sig_m0_cat0, CMS_hgg_sig_m0_absShift)");
  // wAll->factory("prod::CMS_hgg_sig_m0_cat1(massggnewvtx_sig_m0_cat1, CMS_hgg_sig_m0_absShift)");

  // (3) Systematics on resolution: create new sigmas
  // wAll->factory("CMS_hgg_sig_sigmaScale[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat0[1,1.0,1.0]");
  // wAll->factory("CMS_hgg_sig_sigmaScale_cat1[1,1.0,1.0]");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat0(massggnewvtx_sig_sigma0_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat0(massggnewvtx_sig_sigma1_cat0, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_sigma_cat1(massggnewvtx_sig_sigma0_cat1, CMS_hgg_sig_sigmaScale)");
  // wAll->factory("prod::CMS_hgg_sig_gsigma_cat1(massggnewvtx_sig_sigma1_cat1, CMS_hgg_sig_sigmaScale)")

  TString filename(wsDir+TString(fileBaseName)+".inputsig_Bias.root");
  wAll->writeToFile(filename);
  cout << "Write signal workspace in: " << filename << " file" << endl;
  
  return;
}

// Write background pdfs and datasets into the workspace 
void MakeBkgWS(RooWorkspace* w, const char* fileBaseName, double mass, std::string fitName) {

  TString wsDir = "workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;  

  RooWorkspace *wAll = new RooWorkspace("w_all","w_all");
  
  //********************************//
  // Retrieve the datasets and PDFs
  RooDataSet* data[NCAT];
  RooAbsPdf* PhotonsMassBkgPdf[NCAT];
 
  for (int c=0; c<ncat; ++c) {
   
   
    PhotonsMassBkgPdf[c] = (RooAbsPdf*) w->pdf(TString::Format("PhotonsMassBkg_%s_truth_cat%d",fitName,c));
 
   data[c] = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
   ((RooRealVar*) data[c]->get()->find("PhotonsMass"))->setBins(320) ;
    RooDataHist* dataBinned = data[c]->binnedClone();
    w->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c)));
    //  w->Print("V");                                    
    wAll->import(*w->pdf("PhotonsMassBkg_"+TString(fitName)+TString::Format("_truth_cat%d",c)));
    //wAll->import(*PhotonsMassBkgPdf[c]);
    wAll->import(*dataBinned, Rename(TString::Format("data_obs_cat%d",c)));
    wAll->import(*w->data(TString::Format("Data_cat%d",c)), Rename(TString::Format("data_unbinned_obs_cat%d",c)));

  }
  std::cout << "done with importing background pdfs" << std::endl;
  

  TString filename;
  filename = (wsDir+TString(fileBaseName)+TString(fitName)+"_truth_Bias.root");
 
  wAll->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;

  std::cout << std::endl; 
  std::cout << "observation:" << std::endl;
  for (int c=0; c<ncat; ++c) {
    std::cout << "cat " << c << ", " << wAll->data(TString::Format("data_obs_cat%d",c))->sumEntries() << endl;
    wAll->data(TString::Format("data_obs_cat%d",c))->Print();
  }
  std::cout << std::endl;
  
  return;
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

// preparing datacards
void MakeDataCard_1Channel(RooWorkspace* w, const char* fileBaseName, const char* fileBkgName, int iChan) {

  TString cardDir = "datacard/"+filePOSTfix;
  Int_t ncat = NCAT;
  TString wsDir   = "workspaces/"+filePOSTfix;

  // **********************
  // Retrieve the datasets
  cout << "Start retrieving dataset" << endl;
  
  RooDataSet* data[9];
  RooDataSet* signal[9];
  for (int c=0; c<ncat; ++c) {
    data[c]   = (RooDataSet*) w->data(TString::Format("Data_cat%d",c));
    signal[c] = (RooDataSet*) w->data(TString::Format("SigWeight_cat%d",c));
  }

  RooRealVar*  lumi = w->var("lumi");

  // *****************************
  // Print Expected event yields
  cout << "======== Expected Events Number =====================" << endl;  
  cout << ".........Measured Data for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events data: " <<  w->data("Data")->sumEntries()  << endl;
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events data cat%d:   ",c) << data[c]->sumEntries()  << endl;
  }
  cout << ".........Expected Signal for L = " << lumi->getVal() << " pb-1 ............................" << endl;  
  cout << "#Events Signal:      " << w->data("SigWeight")->sumEntries()  << endl;
  Float_t siglikeErr[6];
  for (int c = 0; c < ncat; ++c) {
    cout << TString::Format("#Events Signal cat%d: ",c) << signal[c]->sumEntries() << endl;
    siglikeErr[c]=0.6*signal[c]->sumEntries();
  }
  cout << "====================================================" << endl;  


  // *****************************
  // Printdata Data Card int file
  TString filename(cardDir+TString(fileBaseName)+"_"+"_8TeV"+Form("_channel%d.txt",iChan));
  ofstream outFile(filename);

  outFile << "#CMS-HGG HighMass DataCard for Unbinned Limit Setting, " << lumi->getVal() <<  " pb-1 " << endl;
  outFile << "#Run with: combine -d datacardName.txt -U -m *mass* -H ProfileLikelihood -M MarkovChainMC --rMin=0 --rMax=20.0  -b 3000 -i 50000 --optimizeSim=1 --tries 30" << endl;
  outFile << "# Lumi =  " << lumi->getVal() << " pb-1" << endl;
  outFile << "imax *" << endl;
  outFile << "jmax *" << endl;
  outFile << "kmax *" << endl;
  outFile << "---------------" << endl;

  outFile << "shapes data_obs * " << wsDir+TString(fileBkgName)+".root" << " w_all:data_obs_$CHANNEL" << endl;
  outFile << "shapes sig * "      << wsDir+TString(fileBaseName)+"_8TeV"+".inputsig_4Combine.root" << " w_all:PhotonsMassSig_$CHANNEL" << endl;
  outFile << "shapes bkg * "      << wsDir+TString(fileBkgName)+".root" << " w_all:PhotonsMassBkg_$CHANNEL" << endl;

  outFile << "---------------" << endl;
  outFile << Form("bin          cat%d", iChan) << endl;
  outFile <<  "observation   "  <<  Form("%.10lg",data[iChan]->sumEntries()) << endl;
  outFile << "------------------------------" << endl;
  outFile << "bin                 " << Form("cat%d      cat%d", iChan, iChan) << endl;
  outFile << "process                 sig      bkg" << endl;
  outFile << "process                   0        1" << endl;
  // if(signalScaler==1.)
  // signalScaler=1./signal[2]->sumEntries()*20;
  outFile << "rate                   " 
	  << signal[iChan]->sumEntries()*signalScaler << " " << data[iChan]->sumEntries() << endl;
	  outFile << "--------------------------------" << endl;
  outFile << "# signal scaled by " << signalScaler << endl;

  outFile << "lumi_8TeV       lnN  0.950/1.050  - " << endl;
  // outFile << "CMS_VV_eff_g         lnN  0.8/1.20      - # Signal Efficiency" << endl;
  // outFile << "# Parametric shape uncertainties, entered by hand." << endl;
  // outFile << Form("CMS_hgg_sig_m0_absShift    param   1   0.0125   # displacement of the mean w.r.t. nominal in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_sig_sigmaScale     param   1   0.1   # multiplicative correction to sigmas in EB*EX category, good R9",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_cat%d_norm           flatParam  # Normalization uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope2_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // outFile << Form("CMS_hgg_bkg_8TeV_slope3_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  // if (iChan != 2 )  outFile << Form("CMS_hgg_bkg_8TeV_slope1_cat%d         flatParam  # Mean and absolute uncertainty on background slope",iChan) << endl;
  outFile.close();

  cout << "Write data card in: " << filename << " file" << endl;

  return;
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
      if((mass>1000)&&(mass<2000))
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







// Signal Data Set
void MakeRooKeysPDFMCBkg(RooWorkspace* w, Float_t mass, Bool_t isMirror) {


  TString wsDir = "BiasStudy/workspaces/"+filePOSTfix;
  
  RooWorkspace *wBias = new RooWorkspace("w_bias","w_bias");  


  Int_t ncat = NCAT;
  TString inDir = "";

  Float_t MASS(mass);

  // Luminosity:
  Float_t Lum = 19500.0;  
  RooRealVar lumi("lumi","lumi",Lum);
  wBias->import(lumi); 
  // Variables
  RooArgSet* ntplVars = defineVariables();
  RooArgSet* ntplVars_newweight = defineVariables_newWeight();
  int iMass = abs(mass);  
  RooRealVar* PhotonsMass = w->var("PhotonsMass");    
  TFile sigFile1("histograms_CMS-HGG_24072013.root");   //ggh prod mode tree livia
  
  // common preselection cut
  TString mainCut = "PhotonsMass>=(120) && PhotonsMass<=(480)";   // livia
  
  //get sumEntries of QCD 
  TChain* qcdTree=new TChain();
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_30_8TeV_pf");
  qcdTree->Add("histograms_CMS-HGG_24072013.root/qcd_40_8TeV_pf");

 
  RooDataSet qcdMCWeighted("qcdMCWeighted","MC qcd weighted",qcdTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "qcdMCWeighted" << endl;
  qcdMCWeighted.Print("v");
  Double_t qcdInt_ = qcdMCWeighted.sumEntries();
  cout << "---- nX: qcd Int " << qcdInt_ << endl; 


  //get sumEntries of QCD 
  TChain* gjTree=new TChain();
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  gjTree->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");

 

  RooDataSet gjMCWeighted("gjMCWeighted","MC gj weighted",gjTree,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "gjMCWeighted" << endl;
  gjMCWeighted.Print("v");
  Double_t gjInt_ = gjMCWeighted.sumEntries();
  cout << "---- nX: gj Int " << gjInt_ << endl; 



 
  //chain summing up all production modes
  TChain* sigTree1  = new TChain();
  sigTree1->Add("histograms_CMS-HGG_24072013.root/gjet_20_8TeV_pf");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/gjet_40_8TeV_pf");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/diphojet_8TeV");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_25_8TeV");
  sigTree1->Add("histograms_CMS-HGG_24072013.root/dipho_Box_250_8TeV");
  

  sigTree1->SetTitle("sigTree1");
  sigTree1->SetName("sigTree1");
  
 


  
  // Create signal dataset composed with different productions, the weight is already applied in our ntuples
  RooDataSet BkgMCWeighted("BkgMCWeighted","MC BKG weighted",sigTree1,*ntplVars,mainCut,"evweight");
  cout << endl;
  cout << "BkgMCWeighted" << endl;
  BkgMCWeighted.Print("v");
  cout << "---- nX:  " << BkgMCWeighted.sumEntries() << endl; 
  
  // apply a common preselection cut; split in categories
  cout << endl;
 

  RooDataSet* qcdMC[NCAT];
  RooDataSet* gjMC[NCAT];
  TTree* gjAndQcdTree[NCAT];
  RooDataSet* QcdToGjMC[NCAT];

  RooDataSet* BkgMC[NCAT];
  TTree* BkgMCcopyTree[NCAT];
  RooDataSet* BkgMCcopy[NCAT];
  RooKeysPdf* BkgMCKeyPdf[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw2[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw3[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw4[NCAT];
  RooKeysPdf* BkgMCKeyPdf_bw2_noMirr[NCAT];


  RooDerivative* BkgMCKeyPdf_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw3_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw4_D1[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D1_noMirr[NCAT];

  RooDerivative* BkgMCKeyPdf_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw3_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw4_D2[NCAT];
  RooDerivative* BkgMCKeyPdf_bw2_D2_noMirr[NCAT];

  RooPlot* plotPhotonsMassBkgMC[NCAT];
  RooPlot* plotPhotonsMassBkgMC_D1[NCAT];
  RooPlot* plotPhotonsMassBkgMC_D2[NCAT];
  TCanvas* ctmp = new TCanvas("ctmp","PhotonsMass Background Categories",0,0,500,500);
  Int_t nBinsMass(200);
  Double_t  minMassFit = 120;
  Double_t  maxMassFit = 480;

  //  RooArgSet* argset_ = new RooArgSet(*w->var("PhotonsMass"), *w->var("evweight"));

  for (int c=0; c<NCAT; ++c) {

    // 0) chiara: 1cat only
    // signal[c] =  (RooDataSet*) sigWeighted.reduce(*w->var("massggnewvtx"),mainCut);   //chiara, for 1 cat only

   
    // 1)  prime 4 cat livia


  
   
    //reduce QCD dataset
    if (c==0) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) qcdMC[c] = (RooDataSet*) qcdMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));

 
    //reduce gj dataset
    if (c==0) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars, mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) gjMC[c] = (RooDataSet*) gjMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));
 



    gjMC[c]->Print();



    Float_t qcdInt = qcdMC[c]->sumEntries();
    Float_t gjInt = gjMC[c]->sumEntries();
    Float_t gjEntries = gjMC[c]->numEntries();
    std::cout<<"qcd: "<<qcdInt<<" gj: "<<gjInt<<" gjEntries: "<<gjEntries<<std::endl;
    Float_t qcdWeight = qcdInt/gjInt;
 
    
    gjAndQcdTree[c] = (TTree*)dataset2tree(gjMC[c], ntplVars, qcdWeight);

    QcdToGjMC[c] = new RooDataSet(gjMC[c]->GetName(),gjMC[c]->GetTitle(),gjAndQcdTree[c],*ntplVars_newweight, mainCut, "newweight"  );  
    QcdToGjMC[c]->Print("");
  
    
    if (c==0) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==1) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)<1.4442 && abs(ph2_eta)<1.4442) && (ph1_r9<0.94 || ph2_r9<0.94 ) "));
    if (c==2) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 && ph2_r9>0.94 )"));
    if (c==3) BkgMC[c] = (RooDataSet*) BkgMCWeighted.reduce(*ntplVars,mainCut+TString::Format("&& (abs(ph1_eta)>1.4442 || abs(ph2_eta)>1.4442) && (ph1_r9>0.94 || ph2_r9>0.94 ) "));


    BkgMC[c]->Print();
    BkgMCcopyTree[c] = (TTree*) dataset2tree(BkgMC[c], ntplVars, 1.);
 

    BkgMCcopy[c] =  new RooDataSet(BkgMC[c]->GetName(),BkgMC[c]->GetTitle(),BkgMCcopyTree[c] ,*ntplVars_newweight, mainCut, "newweight"  ); 

    BkgMCcopy[c]->append(*QcdToGjMC[c]);
    BkgMCcopy[c]->Print();


    wBias->import(*BkgMCcopy[c],Rename(TString::Format("BkgMCWeight_cat%d",c)));
    
    cout << "cat " << c << ", BkgMC[c]: " << endl;
    BkgMCcopy[c]->Print("v");
    cout << "---- for category " << c << ", nX for [c]:  " << BkgMCcopy[c]->sumEntries() << endl; 
    cout << endl;


    isMirror = true;
    if(isMirror){
    
    //create the rookeyspdf
    BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_truth_cat%d",c),TString::Format("BkgMCKeyPdf_truth_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth );
    BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_truth_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_truth_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,2 );
    BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_truth_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_truth_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,3 );
    BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_truth_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_truth_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::MirrorBoth,4 );

    }else{

      /*    //create the rookeyspdf
      BkgMCKeyPdf[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_cat%d",c),TString::Format("BkgMCKeyPdf_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror );
      BkgMCKeyPdf_bw2[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw2_cat%d",c),TString::Format("BkgMCKeyPdf_bw2_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,2 );
      BkgMCKeyPdf_bw3[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw3_cat%d",c),TString::Format("BkgMCKeyPdf_bw3_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,3 );
      BkgMCKeyPdf_bw4[c] = new RooKeysPdf(TString::Format("BkgMCKeyPdf_bw4_cat%d",c),TString::Format("BkgMCKeyPdf_bw4_cat%d",c),*PhotonsMass,*BkgMCcopy[c],RooKeysPdf::NoMirror,4 );
      */

    }


    wBias->import(*BkgMCKeyPdf[c]);
    wBias->import(*BkgMCKeyPdf_bw2[c]);
    wBias->import(*BkgMCKeyPdf_bw3[c]);
    wBias->import(*BkgMCKeyPdf_bw4[c]);

    plotPhotonsMassBkgMC[c] = PhotonsMass->frame(minMassFit, maxMassFit,nBinsMass);
    BkgMC[c]->plotOn(plotPhotonsMassBkgMC[c],"PE", MarkerColor(kRed), LineColor(kRed), MarkerSize(1.));
    BkgMCKeyPdf[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4[c]->plotOn(plotPhotonsMassBkgMC[c],"L", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC[c]->SetAxisRange(0.001,plotPhotonsMassBkgMC[c]->GetMaximum()*30.,"Y");
    
    /*
    //create the rookeyspdf_D1

    BkgMCKeyPdf_D1[c] =   (RooDerivative*) BkgMCKeyPdf[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw2_D1[c] =  (RooDerivative*)   BkgMCKeyPdf_bw2[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw3_D1[c] =   (RooDerivative*) BkgMCKeyPdf_bw3[c]->derivative(*PhotonsMass, 1, 0.01);
    BkgMCKeyPdf_bw4_D1[c] =   (RooDerivative*) BkgMCKeyPdf_bw4[c]->derivative(*PhotonsMass, 1, 0.01);


    //create the rookeyspdf_D2
    BkgMCKeyPdf_D2[c] =   (RooDerivative*) BkgMCKeyPdf[c]->derivative(*PhotonsMass, 2, 0.01);
    BkgMCKeyPdf_bw2_D2[c] = (RooDerivative*)  BkgMCKeyPdf_bw2[c]->derivative(*PhotonsMass, 2,0.01);
    BkgMCKeyPdf_bw3_D2[c] = (RooDerivative*)   BkgMCKeyPdf_bw3[c]->derivative(*PhotonsMass, 2,0.01);
    BkgMCKeyPdf_bw4_D2[c] =  (RooDerivative*)  BkgMCKeyPdf_bw4[c]->derivative(*PhotonsMass, 2,0.01);


 
    
    plotPhotonsMassBkgMC_D1[c] = PhotonsMass->frame(minMassFit, maxMassFit);
    BkgMCKeyPdf_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4_D1[c]->plotOn(plotPhotonsMassBkgMC_D1[c],"L", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC_D1[c]->SetAxisRange(plotPhotonsMassBkgMC_D1[c]->GetMinimum()*1.3,plotPhotonsMassBkgMC_D1[c]->GetMaximum()*1.3,"Y");   
    
   
    plotPhotonsMassBkgMC_D2[c] = PhotonsMass->frame(minMassFit, maxMassFit);
    BkgMCKeyPdf_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kBlue), LineWidth(2));
    BkgMCKeyPdf_bw2_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kBlack), LineWidth(2));
    BkgMCKeyPdf_bw3_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(8), LineWidth(2));
    BkgMCKeyPdf_bw4_D2[c]->plotOn(plotPhotonsMassBkgMC_D2[c],"", LineColor(kOrange), LineWidth(2));

    plotPhotonsMassBkgMC_D2[c]->SetAxisRange(plotPhotonsMassBkgMC_D2[c]->GetMinimum()*1.3,plotPhotonsMassBkgMC_D2[c]->GetMaximum()*1.3,"Y");
    

    */
 
    TLegend *leg1 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",c), "brNDC");
    leg1->AddEntry(plotPhotonsMassBkgMC[c]->getObject(0),"Bkg MC","LPE");

    TLegend *leg2 = new TLegend(0.4375,0.7236441,0.85,0.9240678, TString::Format("RooKeysPdf",c), "brNDC");

    TLegend *leg3 = new TLegend(0.2175,0.8236441,0.6575,0.9240678, TString::Format("Category %d",c), "brNDC");
    if(isMirror){
    leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(1),"Default Bw","L");
    leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(2),"Bw x 2","L");
    leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(3),"Bw x 3","L");
    leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4","L");
    }else{
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(1),"Default bw NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(2),"Bw x 2 NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(3),"Bw x 3 NoMirr","L");
      leg2->AddEntry(plotPhotonsMassBkgMC[c]->getObject(4),"Bw x 4 NoMirr","L");
    }
    
    leg1->SetTextSize(0.035);
    leg1->SetTextFont(42);
    leg1->SetBorderSize(0);
    leg1->SetFillStyle(0);
   
    leg2->SetTextSize(0.035);
    leg2->SetTextFont(42);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);

    leg3->SetTextSize(0.035);
    leg3->SetTextFont(42);
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
   
   
    TPaveText* label_cms = get_labelCMS(0, "2012", true);
    TPaveText* label_sqrt = get_labelSqrt(0);

  
    ctmp->cd();
    
    plotPhotonsMassBkgMC[c]->Draw();  
    leg1->Draw("same");
    leg2->Draw("same");
    
    label_cms->Draw("same");
    label_sqrt->Draw("same");

   
    ctmp->SetLogy(1);
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_"+TString::Format("cat%d_LOG.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_NoMirr_"+TString::Format("cat%d_LOG.root", c));
  
    }
       ctmp->SetLogy(0);
       /*   //plot D1
    plotPhotonsMassBkgMC_D1[c]->Draw();  
    leg2->SetHeader("RooKeysPDF - 1st Derivative");
    leg2->Draw("same");
    leg3->Draw("same");
    
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_"+TString::Format("cat%d.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D1_NoMirr_"+TString::Format("cat%d_LOG.root", c));
      
    }  
    ctmp->SetLogy(0);
    ctmp->Clear();

    //plot D2
    plotPhotonsMassBkgMC_D2[c]->Draw();  
    leg2->SetHeader("RooKeysPDF - 2nd Derivative");
    leg2->Draw("same");
    leg3->Draw("same");
	
    label_cms->Draw("same");
    label_sqrt->Draw("same");
    if(isMirror){
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.png", c));
    ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_"+TString::Format("cat%d.root", c));
    }else{
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.png", c));
      ctmp->SaveAs("plots/BKG_MC_rookeyspdf_D2_NoMirr_"+TString::Format("cat%d_LOG.root", c));
      
      } */
  }

  

  std::cout << "done with importing MC background datasets and RooKeysPdfs" << std::endl;
  TString filename(wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV.root");
  wBias->writeToFile(filename);
  cout << "Write background workspace in: " << filename << " file" << endl;


  
  return;

}



void MakeWSBiasStudy(RooWorkspace* w, Float_t mass, Int_t perc){

  TString wsDir = "BiasStudy/workspaces/"+filePOSTfix;
  Int_t ncat = NCAT;
  
  RooWorkspace *wAll = new RooWorkspace("wBiasTruth","wBiasTruth");  
  TFile* f_roo = new TFile(""+wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV.root");
  RooWorkspace *wRoo = f_roo->Get("w_bias");
  TFile* f_roo123 = new TFile(""+wsDir+"HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat123.root");
  RooWorkspace *wRoo123 = f_roo123->Get("w_bias");
  
  wAll->import(*w->var("RooMINmass"));
  wAll->import(*w->var("RooMAXmass"));
 
  for (int c=0; c<ncat; ++c) {
  //import Bkg:
  wAll->import(*w->data(TString::Format("data_obs_cat%d",c)), Rename(TString::Format("data_obs_truth_cat%d",c)));
  wAll->import(*w->data(TString::Format("Data_cat%d",c)), Rename(TString::Format("data_unbinned_obs_truth_cat%d",c)));
  wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_DiJet_truth_cat%d",c)));
  //wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_DiJetPL_truth_cat%d",c)));
  //wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_DiJetEXP_truth_cat%d",c)));
  //wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_DiJetEXPOL_truth_cat%d",c)));

  wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_Expol_truth_cat%d",c)));
  wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_Lau_truth_cat%d",c)));
  wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_ExpPAR_truth_cat%d",c)));
  // wAll->import(*w->pdf(TString::Format("PhotonsMassBkg_ExpolPL_truth_cat%d",c)));
 


  //import signal:
  std::cout<<"-----------------> SIGNAL"<<std::endl;
  wAll->import(*w->pdf(TString::Format("PhotonsMassSig_cat%d",c)));
  if(mass<400 && mass!=350)wAll->import(*w->data(TString::Format("SigWeight_cat%d",c)), Rename(TString::Format("SigWeight_truth_cat%d",c)));
  }
  
  //import rookeys
  std::cout<<"-----------------> ROOKEYS"<<std::endl;

  wAll->import(*wRoo->pdf(TString::Format("BkgMCKeyPdf_bw2_cat0",c)));
  wAll->import(*wRoo->pdf(TString::Format("BkgMCKeyPdf_bw3_cat0",c)));
  wAll->import(*wRoo->pdf(TString::Format("BkgMCKeyPdf_bw4_cat0",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw2_cat1",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw3_cat1",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw4_cat1",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw2_cat2",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw3_cat2",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw4_cat2",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw2_cat3",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw3_cat3",c)));
  wAll->import(*wRoo123->pdf(TString::Format("BkgMCKeyPdf_bw4_cat3",c)));
 

  int imass(mass);
  int iperc(perc);
  TString filename(wsDir+TString::Format("HighMass-hgg.m%d.0_w%d_inputBias_truth_8TeV_%d_%d.root", imass,iperc, MINmass, MAXmass));
  wAll->writeToFile(filename);
  cout << "Write gen workspace in: " << filename << " file" << endl;


  return;


}





RooHistFunc* getRooHistFunc(int cat, RooRealVar* var){


  double mass[5] = {150., 200., 250., 300., 400.};
  double c0[5] = {57.2743,73.98, 86.1266,99.099,115.707};
  double c1[5] = {69.246,76.9358,75.4159,72.8432,69.0506 };
  double c2[5] = {27.6973,33.0972,34.1901,37.6118,36.0121 };
  double c3[5] = {54.2063,60.2946,59.8451,60.5352,53.4554 };
  double all[5] = {187.921, 217.376, 226.456,237.32,241.692};


  TH1F* h_all = new TH1F("h_all", "h_all", 6, 130, 450);
  for(int i=0;i<5;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]/19620.<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]/19620.);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]/19620.);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]/19620.);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]/19620.);
    if(cat==4) h_all->SetBinContent(h_all->FindBin(mass[i]),all[i]/19620.);
   
  }
  std::cout<<"--------------------------->"<<cat<<std::endl;

  if(cat==0) h_all->SetBinContent(h_all->FindBin(350),(c0[3]+c0[4])/2/19620.);
  if(cat==1) h_all->SetBinContent(h_all->FindBin(350),(c1[3]+c1[4])/2/19620.);
  if(cat==2) h_all->SetBinContent(h_all->FindBin(350),(c2[3]+c2[4])/2/19620.);
  if(cat==3) h_all->SetBinContent(h_all->FindBin(350),(c3[3]+c3[4])/2/19620.);
  if(cat==4) h_all->SetBinContent(h_all->FindBin(350),(all[3]+all[4])/2/19620.);
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  //  TPaveText* label_cms = get_labelCMS(0, false);
  //TPaveText* label_sqrt = get_labelSqrt(0);  
  h_all->GetYaxis()->SetRangeUser(0.0001, 0.1);
  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Signal Yield");
  h_all->Draw();

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassSig_cat%d_norm",cat),TString::Format("PhotonsMassSig_cat%d_norm",cat), *var,*rooData_all, 3);
  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/signalYield.png");
  c->SaveAs("plots/signalYield.pdf");

  return rooFunc_all;

}



RooHistFunc* getRooHistFuncFitMIN(int cat, RooRealVar* var, TF1* f){


  double mass[8] = {150., 250, 350, 450, 550, 650, 750, 850};
  double c0[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c1[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c2[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
  double c3[8] = {130., 200., 250., 300., 300., 450., 450., 550. };
 
  f = new TF1("f", "[0]+[2]*pow(x, [1])", 100, 900);
  f->SetParameter(0,80 );
  f->SetParameter(1, 1.2 );
  //f->SetParameter(2, 0.8 );
  f->SetParameter(2, 0.14 );
  // f->SetParameter(4, 2 );
  
  TH1F* h_all = new TH1F("h_all", "h_all", 8, 100, 900);
  for(int i=0;i<8;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Minimum");
  h_all->Draw("");
  h_all->Fit("f");
  f->Draw("Lsame");
  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassBkg_cat%d_min",cat),TString::Format("PhotonsMassBkg_cat%d_min",cat), *var,*rooData_all, 3);
  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  //plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/.png");
  c->SaveAs("plots/Bkg_fitMinimum.pdf");



  return rooFunc_all;

}



RooHistFunc* getRooHistFuncFitMAX(int cat, RooRealVar* var, TF1* f){

  
  double mass[8] = {150,250, 350, 450, 550, 650, 750, 850};
  double c0[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c1[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c2[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};
  double c3[8] = {230, 550, 650, 700, 800, 800, 1000, 1000};

 f = new TF1("f", "[0]+[2]*pow(x, [1])", 100, 900);
  f->SetParameter(0,80 );
  f->SetParameter(1, 1.2 );
  //f->SetParameter(2, 0.8 );
  f->SetParameter(2, 0.14 );
  // f->SetParameter(4, 2 );
  

  TH1F* h_all = new TH1F("h_all", "h_all", 8, 100, 900);
  for(int i=0;i<8;i++){
    std::cout<<cat<<std::endl;
    std::cout<<mass[i]<<"  "<<h_all->FindBin(mass[i])<<"   "<<c0[i]<<std::endl;
    if(cat==0) h_all->SetBinContent(h_all->FindBin(mass[i]),c0[i]);
    if(cat==1) h_all->SetBinContent(h_all->FindBin(mass[i]),c1[i]);
    if(cat==2) h_all->SetBinContent(h_all->FindBin(mass[i]),c2[i]);
    if(cat==3) h_all->SetBinContent(h_all->FindBin(mass[i]),c3[i]);
      
  }
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0,"2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);  

  h_all->GetXaxis()->SetTitle("m_{H} [GeV]");
  h_all->GetYaxis()->SetTitle("Fit Range: Maximum");
  h_all->Draw();
  h_all->Fit("f");
  f->Draw("Lsame");
 
  

  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*var,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc(TString::Format("PhotonsMassBkg_cat%d_max",cat),TString::Format("PhotonsMassBkg_cat%d_max",cat), *var,*rooData_all, 7);

  RooPlot* plot = var->frame();
  rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot);
  plot->Draw("same");
  c->SetLogy();

  c->SaveAs("plots/.png");
  c->SaveAs("plots/Bkg_fitMaximum.pdf");

  return rooFunc_all;

}



makeRooHistPlot(RooHistFunc* rooFitMin,RooHistFunc* rooFitMax ,TF1* fmin, TF1* fmax, RooRealVar* var){

  TCanvas* cfit = new TCanvas("cfit", "cfit", 1);
  cfit->cd();
  
  RooPlot* p = var->frame();
  rooFitMin->plotOn(p, LineColor(kAzure+9));
  rooFitMax->plotOn(p, LineColor(kViolet+9));
  p->Draw();
  fmax->SetLineColor(kViolet+9);
  fmax->SetLineStyle(kDashed);
  fmin->SetLineColor(kAzure+9);
  fmin->SetLineStyle(kDashed);
  fmax->Draw("same");
  fmin->Draw("same");
  p->GetYaxis()->SetTitle("m_{#gamma #gamma} [GeV]");
  p->GetXaxis()->SetTitle("m_{H} [GeV]");

  TLegend* leg = new TLegend(0.2, 0.65, 0.45, 0.89, "", "brNDC");
  leg->SetTextSize(0.0206044);  
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(p->getObject(1), "Fit Maximum", "L");
  leg->AddEntry(p->getObject(0), "Fit Minimum", "L");
  leg->Draw("same");
  

  TPaveText* label_cms = get_labelCMS(0, "2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  cfit->SaveAs("plots/bkg_fitRange.png");
  cfit->SaveAs("plots/bkg_fitRange.pdf");
  cfit->SaveAs("~/www/BkgFit/bkg_fitRange.pdf");
  cfit->SaveAs("~/www/BkgFit/bkg_fitRange.png");
 



}
