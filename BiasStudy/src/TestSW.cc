#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>
#include <cassert>

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TTree.h" 

#include "Rtypes.h"

#include "RooMsgService.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooConstVar.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooFit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooNLLVar.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooBernstein.h"
#include "RooMultiVarGaussian.h"
#include "RooChebychev.h"
#include "RooRandom.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooKeysPdf.h"
#include "TSystem.h"
#include "TROOT.h"
#include "RooCBShape.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/RooCBCrujffPdf.h"
#include "/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/RooCBCrujffPdf.cxx"
#include "/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/dict.h"
#include "/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/dict.cc"
#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/include/Riostream.h"

using namespace std;
using namespace RooFit;

#define oldFitFunDefinition
const int nFuncs=18;


string funcNames[nFuncs] = {"1pol","2pol","3pol","4pol","5pol","1exp","2exp","3exp","1pow","2powlow","3powlow","2lau","4lau","6lau","RooKey","DiJet", "Expol","Expol2"};

void checkFunc(string name){
  if (name!="1pol" && name!="2pol" && name!="3pol" && name!="4pol" && name !="5pol" && name!="1exp" && name!="2exp" && name!="3exp" && name!="1pow" && name!="2powlow" && name!="3powlow" && name!="2lau" && name!="4lau" && name!="6lau" && name!="BkgMCKeyPdf"&& name!="DiJet"&& name!="Expol"&& name!="Expol2" && name!= "RooKey") {
    cout << "Invalid function: " << name << endl;
    cout << "Options are: " << endl;
    for (int f=0; f<nFuncs; f++){
      cout << "   " << funcNames[f].c_str() << endl;
    }
    exit(1);
  }
}

void checkMass(int mMC){
  if (mMC!=110 && mMC!=115 && mMC!=120 && mMC!=125  && mMC!=130 && mMC!=135 && mMC!=140 && mMC!=150  && mMC!=200  && mMC!=250  && mMC!=300  && mMC!=400 ) cout << "Invalid mass: " << mMC << endl; 
  assert(mMC==110 || mMC==115 || mMC==120 || mMC==125  || mMC==130 || mMC==135 || mMC==140 || mMC==150 || mMC==200 || mMC==250 || mMC==300 || mMC==400);
}

const int getPar(string name){
//   if (name=="1pol" || name=="1exp" || name=="1pow" || name=="2lau") return 1;
//   else if (name=="2pol") return 2;
//   else if (name=="3pol" || name=="2exp" || name=="exp2_1" || name=="exp2_2" || name=="2pow" || name=="4lau" || name=="pow2_1" || name=="pow2_2") return 3;
//   else if (name=="4pol") return 4;
//   else if (name=="5pol" || name=="3exp" || name=="exp3_1" || name=="exp3_2" || name=="exp3_3" || name=="3pow" || name=="6lau") return 5;
//   else {
//     std::cerr << "[ERROR] getPar not defined" << std::endl;
//     exit(1);
//   }

  if (name=="1exp" || name=="1pow" || name=="2lau") return 1;
  else if (name=="1pol" ||name=="Expol") return 2;
  else if (name=="2pol" || name=="2exp" || name=="exp2_1" || name=="exp2_2" || name=="2pow" || name=="4lau" || name=="pow2_1" || name=="pow2_2" ||name=="DiJet"||name=="Expol2") return 3;
  else if (name=="3pol") return 4;
  else if (name=="4pol" || name=="3exp" || name=="exp3_1" || name=="exp3_2" || name=="exp3_3" || name=="3pow" || name=="6lau") return 5;
  else {
    std::cerr << "[ERROR] getPar for " << name << " not defined" << std::endl;
    exit(1);
  }

}

string getType(string name){
  if (name=="1exp" || name=="2exp" || name=="exp2_1" || name=="exp2_2" || name=="3exp" || name=="exp3_1" || name=="exp3_2" || name=="exp3_3") return "exp";
  else if (name=="1pol" || name=="2pol" || name=="3pol" || name=="4pol" || name=="5pol") return "pol";
  else if (name=="1pow" || name=="2powlow" || name=="3powlow" || name=="pow2_1" || name=="pow2_2") return "powlow";
  else if (name=="2lau" || name=="4lau" || name=="6lau") return "lau";
  else if (name=="DiJet") return "DiJet";
  else if (name=="Expol") return "Expol";
  else if (name=="Expol2") return "Expol2";
  else if (name== "RooKey") return "RooKey";
  else {
    std::cerr << "[ERROR] Type not defined" << std::endl;
    exit(1);
  }
}

int main(int argc, char* argv[]){


  bool help=false;
  bool verbose=false;
  bool doBkgInt=false;
  bool plotGen=false;
  bool saveDataFit=false;
  bool wideRange=false;
  bool sigGen=false;
  bool fixed=false;
  bool onlyB=false;
  float muGen=1;
  int nToys;
  int nJobs=1;
  int jobNumb=0;
  int toyStep=2;
  string genName;
  string fitName;

  int mMC;
  int cat;
  string fileFitName;
  string fileGenName;
  string data_fileName;




  for (int arg=0; arg<argc; arg++) {
    //    std::cout << string(argv[arg]) << std::endl;
    if (string(argv[arg])=="-h" || string(argv[arg])=="--help") help=true;
    if (string(argv[arg])=="-v") verbose=true;
    if (string(argv[arg])=="-b") onlyB=true;
    if (string(argv[arg])=="-bkg") doBkgInt=true;
    if (string(argv[arg])=="-pG") plotGen=true;
    if (string(argv[arg])=="-sDF") saveDataFit=true;
    if (string(argv[arg])=="-gen") genName=string(argv[arg+1]);
    if (string(argv[arg])=="-fit") fitName=string(argv[arg+1]);
    if (string(argv[arg])=="-t") nToys=atoi(argv[arg+1]);
    if (string(argv[arg])=="-p") toyStep=atoi(argv[arg+1]);//la calcola da solo
    if (string(argv[arg])=="-nJ") nJobs=atoi(argv[arg+1]);
    if (string(argv[arg])=="-j") jobNumb=atoi(argv[arg+1]);
    if (string(argv[arg])=="-wide") wideRange=true;


    if (string(argv[arg])=="-f_fit") fileFitName=string(argv[arg+1]);
    if (string(argv[arg])=="-f_gen") fileGenName=string(argv[arg+1]);
    if (string(argv[arg])=="-dataf") data_fileName=string(argv[arg+1]);
    if (string(argv[arg])=="-fixed") fixed=true;
    if (string(argv[arg])=="-sigGen") sigGen=true;
    if (string(argv[arg])=="-muGen")  muGen=atof(argv[arg+1]);
    if (string(argv[arg])=="-mass") mMC=atoi(argv[arg+1]);
    if (string(argv[arg])=="-cat")  cat=atoi(argv[arg+1]);
  }
  if(toyStep<0)toyStep=nToys/1;
  if (!help){
    //    checkFunc(genName);
    //    checkFunc(fitName);
    //    checkMass(mMC);
  }
  if (argc<9 || help){
    cout << "--- Run with following options: ---" << endl;
    cout << "    -t    nToys " << endl;
    cout << "    -p    plotStep " << endl;
    cout << "    -m    to specifiy mass " << endl;

    cout << "    -gen  $i to gen with func $i " << endl;
    cout << "    -fit  $i to fit with func $i " << endl;
    cout << "    -nJ   number of jobs " << endl;
    cout << "    -j    job number " << endl;
    cout << "--- Additional options: ---" << endl;
    cout << "    -v    for diagnostics " << endl;
    cout << "    -bkg  to do bkg int " << endl;
    cout << "    -pG   to plot gen func " << endl;
    cout << "    -sDF  to save data fit " << endl;
    cout << "    -wide to fit from 100-180 " << endl;
    exit(1);
  }


 
  string genType=getType(genName);
  string fitType=getType(fitName);
  
  cout << "--- Running with following options ---" << endl;
  cout << "    nToys:            " << nToys << endl;
  cout << "    plotStep:         " << toyStep << endl;
  cout << "    genFunction:      " << genName << endl;
  cout << "    fitFunction:      " << fitName << endl;
  cout << "    mass:             " << mMC << endl;

  cout << "    nJobs:            " << nJobs << endl;
  cout << "    jobNum:           " << jobNumb << endl;
  cout << "    cat:           " << cat<< endl;
  cout << "    mMC:           " << mMC<< endl;
  if (verbose)     cout << "    print fit results on " << endl;
  if (doBkgInt)    cout << "    bkg integral on " << endl;
  if (plotGen)     cout << "    plot gen function on " << endl;
  if (saveDataFit) cout << "    save data fit on " << endl;
  if (wideRange)   cout << "    wide range fit on " << endl;
  if (sigGen)      cout << "    generate signal + bkg on" << endl;
  if (saveDataFit && nToys>0){
    cout << " ERROR: CANNOT SAVE STARTING PARAMS AND FIT TOYS SIMULTANEOUSLY" << endl;
    cout << " --> either run without -sDF option or without -t or with -t 0 for no toys" << endl;
    exit(1);
  }


  int toysPerJob = nToys/nJobs;
  
  if (!verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  
  //  std::cout<<"//--------------------------------------------------------------------------------//"<<std::endl;
  
  //  std::cout<<"Apro il file di input e prendo il WS che si chiama wBias e salvo la variabile della massa"<<std::endl;
  
 
  
  TFile *genFile = new TFile(fileGenName.c_str());

  std::cout<<"GENWS: "<<fileGenName.c_str()<<std::endl;
  RooWorkspace *genWS;
  
  genWS = (RooWorkspace*)genFile->Get("wBiasTruth");
  if(genWS == NULL){
    std::cerr << "[ERROR] genWS not found in WS" << std::endl;
    exit(1);
  }


  TFile *rooFile;
  if(cat==0 || cat == 1) rooFile = new TFile("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/workspaces/HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat01.root");
  else if(cat==2) rooFile = new TFile("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/workspaces/HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat2.root");
  else if(cat==3) rooFile = new TFile("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/workspaces/HighMass-hgg.RooKeysPdfMCBkg_8TeV_cat3.root");
  RooWorkspace *rooWS;
  
  rooWS = (RooWorkspace*)rooFile->Get("w_bias");
  if(rooWS == NULL){
    std::cerr << "[ERROR] rooWS not found in WS" << std::endl;
    exit(1);
  }




 RooRealVar *mass=(RooRealVar*)genWS->var("PhotonsMass");

 if(mass == NULL){
   std::cerr << "[ERROR] variable mass not found in WS: " << genWS->GetName() << std::endl;
   return 1;
 }

 
 //--------------------------------------------------------------------------------//
 // std::cout<<"Creo cartelle di output"<<std::endl;
 
 
 system("if [ ! -e \"/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults\" ]; then mkdir /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults; fi");
 
 system("mkdir -p /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/SWplots/toys");
 system("mkdir -p /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/SWplots/data");
 system("mkdir -p /afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/SWplots/histos");
 
 
 //--------------------------------------------------------------------------------//
 std::cout<<"Dichiaro le variabili che poi salvero' nel tree e creo il file di output"<<std::endl;
 
 
 RooDataSet *data;
 RooAbsPdf *genFcn;
 RooAbsPdf *genParamatersFcn;
 
 RooAbsPdf *fitFcn;

 //par Expol 
 RooRealVar p1Expol(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_cat%d",cat)),"p1Expol");
 RooRealVar p2Expol(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_cat%d",cat)),"p2Expol");


//par Expol2 
 RooRealVar p1Expol2(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_expol1_2_cat%d",cat)),"p1Expol2");
 RooRealVar p2Expol2(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_expol2_2_cat%d",cat)),"p2Expol2");
 RooRealVar p3Expol2(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_expol3_2_cat%d",cat)),"p3Expol2");

 std::cout<<"dopo expol2"<<std::endl;
 //par DiJet
 RooRealVar p1DiJet(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_slope1_cat%d",cat)),"p1DiJet");
 RooRealVar p2DiJet(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_slope2_cat%d",cat)),"p2DiJet");
 RooRealVar p3DiJet(*genWS->var(TString::Format("PhotonsMass_bkg_8TeV_slope3_cat%d",cat)),"p3DiJet");

  std::cout<<"dopo dijet"<<std::endl;

 double sigSMEvents;
 
 RooRealVar *mu = new RooRealVar("mu","mu",1.,-20,20);
 
 double fwhm;
 double maxMass;
 
 TH1F *bias;
 TH1F *muHist;

 Float_t mu_;
 Float_t muTruth_;
 Float_t sigma_mu_;
 Float_t bkgSig1fwhm_;
 Float_t bkgSig2fwhm_;
 Float_t bkgTrue1fwhm_;
 Float_t bkgTrue2fwhm_;
 Float_t bkgErrSig1fwhm_;
 Float_t bkgErrSig2fwhm_;
 Float_t bkgErrNormSig1fwhm_;
 Float_t bkgErrNormSig2fwhm_;

 Float_t    sigYieldTot_;
 Float_t    sigYieldTotErr_;
 Float_t    sigYield1fwhm_;
 Float_t    sigYield2fwhm_;
 Float_t    sigYield1fwhmErr_;
 Float_t    sigYield2fwhmErr_;
 Float_t    sigYield1fwhmErrNorm_;
 Float_t    sigYield2fwhmErrNorm_;
 Float_t    sigYield1Truefwhm_;
 Float_t    sigYield2Truefwhm_;



 UInt_t iToy_;
 Int_t fitStatus_;
 Int_t cat_;
 Int_t mMC_;
 Float_t chi2_;
 char genName_[30];
 char fitName_[30];
 
 sprintf(genName_, "%s", genName.c_str());
 sprintf(fitName_, "%s", fitName.c_str());
 
 // std::cout << "getName = " << genName_ << std::endl;
 // std::cout << "fitName = " << fitName_ << std::endl;
 
 //std::cout << "[INFO] creating new output file" << std::endl;
 TFile *outFile;
 
 if(!sigGen)  outFile = new TFile(Form("TestSWResults/biasCheck_g%s_f%s_m%d_cat%d_j%d.root",genName.c_str(),fitName.c_str(),mMC,cat,jobNumb),"RECREATE");
 else outFile = new TFile(Form("TestSWResults/biasCheck-sigGen_%g-g%s_f%s_m%d_cat%d_j%d.root",muGen, genName.c_str(),fitName.c_str(),mMC,cat,jobNumb),"RECREATE");
 
 TTree *muTree = new TTree("muTree","muTree");
 muTree->Branch("mass", &mMC, "mass/I");
 muTree->Branch("mu", &mu_, "mu/F");
 muTree->Branch("muTruth", &muTruth_, "muTruth/F");
 muTree->Branch("sigma_mu", &sigma_mu_, "sigma_mu/F");
 muTree->Branch("bkgTrue1fwhm", &bkgTrue1fwhm_, "bkgTrue1fwhm/F");
 muTree->Branch("bkgTrue2fwhm", &bkgTrue2fwhm_, "bkgTrue2fwhm/F");
 muTree->Branch("bkgSig1fwhm", &bkgSig1fwhm_, "bkgSig1fwhm/F");
 muTree->Branch("bkgSig2fwhm", &bkgSig2fwhm_, "bkgSig2fwhm/F");
 muTree->Branch("bkgErrSig1fwhm", &bkgErrSig1fwhm_, "bkgErrSig1fwhm/F");
 muTree->Branch("bkgErrSig2fwhm", &bkgErrSig2fwhm_, "bkgErrSig2fwhm/F");
 muTree->Branch("bkgErrNormSig1fwhm", &bkgErrNormSig1fwhm_, "bkgErrNormSig1fwhm/F");
 muTree->Branch("bkgErrNormSig2fwhm", &bkgErrNormSig2fwhm_, "bkgErrNormSig2fwhm/F");

 
 muTree->Branch("sigYieldTot", &sigYieldTot_, "sigYieldTot/F");
 muTree->Branch("sigYieldTotErr", &sigYieldTotErr_, "sigYieldTotErr/F");
 muTree->Branch("sigYield1fwhm", &sigYield1fwhm_, "sigYield1fwhm/F");
 muTree->Branch("sigYield2fwhm", &sigYield2fwhm_, "sigYield2fwhm/F");
 muTree->Branch("sigYield1fwhmErr", &sigYield1fwhmErr_, "sigYield1fwhmErr/F");
 muTree->Branch("sigYield2fwhmErr", &sigYield2fwhmErr_, "sigYield2fwhmErr/F");
 muTree->Branch("sigYield1fwhmErrNorm", &sigYield1fwhmErrNorm_, "sigYield1fwhmErrNorm/F");
 muTree->Branch("sigYield2fwhmErrNorm", &sigYield2fwhmErrNorm_, "sigYield2fwhmErrNorm/F");
 muTree->Branch("sigYield1Truefwhm", &sigYield1Truefwhm_, "sigYield1Truefwhm/F");
 muTree->Branch("sigYield2Truefwhm", &sigYield2Truefwhm_, "sigYield2Truefwhm/F");



 muTree->Branch("genFun", genName_, "genFun/C");
 muTree->Branch("fitFun", fitName_, "fitFun/C");
 muTree->Branch("iToy", &iToy_, "iToy/i");
 muTree->Branch("cat", &cat_, "cat/i");
 muTree->Branch("mMC", &mMC_, "mMC/i");
 muTree->Branch("chi2", &chi2_, "chi2/F");
 muTree->Branch("fitStatus",&fitStatus_, "fitStatus/I");
 
 
 TLegend *leg = new TLegend(0.65,0.65,0.89,0.89);
 leg->SetLineColor(0);
 leg->SetFillColor(0);
 TH1F *h = new TH1F("h","h",1,0,1);
 h->SetLineColor(kMagenta);
 h->SetLineWidth(3);
 TH1F *h1 = new TH1F("h1","h2",1,0,1);
 h1->SetLineColor(kBlue);
 h1->SetLineWidth(3);
 TH1F *h2 = new TH1F("h2","h2",1,0,1);
 h2->SetLineColor(kRed);
 h2->SetLineStyle(kDashed);
 h2->SetLineWidth(3);
 if (plotGen) leg->AddEntry(h,"Truth","l");
 leg->AddEntry(h2,"bkg part of s+b","l");
 leg->AddEntry(h1,"s+b after fit","l");
 
 
 muHist = new TH1F(Form("muHist_g%s_f%s_m%d",genName.c_str(),fitName.c_str(),mMC),Form("muHist_g%s_f%s_m%d",genName.c_str(),fitName.c_str(),mMC),100,-30,30);
 if (doBkgInt)  bias = new TH1F(Form("b_g%s_f%s_m%d_c%d",genName.c_str(),fitName.c_str(),mMC,cat),Form("b_g%s_f%s_m%d_c%d",genName.c_str(),fitName.c_str(),mMC,cat),100,-100,100);
 
 
 double toData;
 
 std::cout<<"//--------------------------------------------------------------------------------//"<<std::endl;
 // cout << "Dichiaro le funzioni di fit e prendo il dataset dei dati da cui ricavo la normalizzazione per il bkg" << endl;
 
 
 
 cout << "1. Salvo il dataset" << cat << endl;
 data = (RooDataSet*)genWS->data(Form("data_unbinned_obs_truth_cat%d",cat));
 if(data == NULL){
    std::cout << "[WARNING] roodataset from genWS not found" << std::endl;
    exit(0);
  }
 //  cout<<" Le entries del dataset (DATI) sono: "<<data->sumEntries()<<std::endl;
  
  

 //  cout << "2. Salvo la funzione con cui generare il fondo"<< endl;

 //add gen func by livia  
 if(genName == "RooKey") genFcn = (RooKeysPdf*) genWS->pdf(Form("BkgMCKeyPdf_bw2_cat%d",cat));
 
 if(genName =="DiJet") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_DiJet_truth_cat%d",cat));
 if(genName =="Expol") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_Expol_truth_cat%d",cat));
 if(genName =="Expol2") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_Expol2_truth_cat%d",cat));
 if(genName =="2exp") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_Exp_truth_cat%d_2exp",cat));
 if(genName =="2powlow") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_PowLow_truth_cat%d_2powlow",cat));
   
  //se dovesse servire di generare i par iniziale con una multivariata
    genParamatersFcn=0;
    //  if (genName=="2pol") genParamatersFcn = (RooMultiVarGaussian*)genWS->pdf(Form("pdf_nll2"));
  
    if(genFcn==NULL){
      std::cout << "[ERROR] genFcn not defined" << std::endl;
      return 1;
    }

    if (genParamatersFcn!=NULL) genParamatersFcn->Print();


    // ------------------------------------------
    // --- params for fit to gen data ---
    // ------------------------------------------

    //    cout << "3. Salvo la funzione di fit per il bkg" << endl;
  
    if (!fixed)
      {

	//define fit functions



	if(fitName =="Expol")fitFcn = new RooGenericPdf("fitFcn_expol", "exp(-@0/(@1+@2*@0))", RooArgList(*mass, p1Expol, p2Expol));
	if(fitName =="Expol2")fitFcn = new RooGenericPdf("fitFcn_expol2", "exp(-@0*@0/(@1+@2*@0+@3*@0*@0))", RooArgList(*mass, p1Expol2, p2Expol2, p3Expol2));
	if(fitName =="DiJet")fitFcn = new RooGenericPdf("fitFcn_dijet","fitFcn_dijet", "pow((1-(@0/8000.)), @2)/pow((@0/8000.), @1+@3*log(@0/8000.))", RooArgList(*mass, p1DiJet, p2DiJet, p3DiJet));

	if(fitFcn==NULL){
	  std::cout << "[ERROR] fitFcn not defined" << std::endl;
	  return 1;
	}
	

      }
     else
       {
	
	if(fitFcn==NULL){
	  std::cout << "[ERROR] fitFcn not defined" << std::endl;
	  return 1;
	}


	RooArgSet* par=fitFcn->getParameters(*mass);
 
	//	par->Print("v");

	//Setting all parameters fixed
	TIterator* iter = par->createIterator();
	for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
          RooRealVar *v = dynamic_cast<RooRealVar *>(a);
	  std::cout << "Setting " << v->GetName() << " fixed " << std::endl;
	  v->setConstant(true);
	}
	delete iter;
      }


 
    //    cout << "Dichiaro le funzioni di fit per il segnale" << endl;


    // ---- Get mass distributions  and entries ------

   
      cout << "1. Prendo il dataset del segnale e salvo l'integrale" << endl;
      RooDataSet *sigRooData = (RooDataSet*) genWS->data(TString::Format("SigWeight_truth_cat%d",cat));
      sigSMEvents = sigRooData->sumEntries();

  
      /*     cout << "2A. Costruisco la pdf Sig" << endl;
       //CBCrujff      
      RooRealVar aCB(Form("aCB_m%d_cat%d",mMC,cat),Form("aCB_m%d_cat%d",mMC,cat), 1.489, 0., 10.);
      RooRealVar aC(Form("aC_m%d_cat%d",mMC,cat),Form("aC_m%d_cat%d",mMC,cat), 0.127,0., 10.);
      RooRealVar mean(Form("mean_m%d_cat%d",mMC,cat),Form("mean_m%d_cat%d",mMC,cat) , 149.81, 0., 300.);
      RooRealVar n(Form("n_m%d_cat%d",mMC,cat),Form("n_m%d_cat%d",mMC,cat), 7.6, 0., 10.);     
      RooRealVar sigma(Form("sigma_m%d_cat%d",mMC,cat),Form("sigma_m%d_cat%d",mMC,cat), 1.64, 0., 10.);   
      RooCBCrujffPdf* sigMCPdf=  new RooCBCrujffPdf(Form("sigMCPdf_m%d_cat%d",mMC,cat), Form("sigMCPdf_m%d_cat%d",mMC,cat), *mass, mean, sigma,aC, aCB, n); 
      RooArgSet* paramsSig = sigMCPdf->getParameters(*mass);
      paramsSig->readFromFile("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/initParSig.txt");
      paramsSig->writeToStream(std::cout, kFALSE);     
      aCB.setConstant();
      aC.setConstant();
      mean.setConstant();
      n.setConstant();
      sigma.setConstant();

*/
      cout << "2B. Prendo la pdf Sig dal WS" << endl;
      RooFFTConvPdf* sigMCPdf = (RooFFTConvPdf*) genWS->pdf(TString::Format("PhotonsMassSig_cat%d",cat));
	if(sigMCPdf==NULL){
	  std::cout << "[ERROR] sigMCPdf not defined" << std::endl;
	  return 1;
	}
	
 

      cout << "3. Faccio un roodatahist" << endl;
      ((RooRealVar*) sigRooData->get()->find("PhotonsMass"))->setBins(400) ; 
      RooDataHist* sigRooDataHist = sigRooData->binnedClone();
      TH1F* sigDataHist = new TH1F("sigDataHist", "sigDataHist", 400, 130., 450.);
      sigRooDataHist->fillHistogram(sigDataHist, *mass, "" );
     
      int bin1=sigDataHist->FindFirstBinAbove(sigDataHist->GetMaximum()/2.);
      int bin2=sigDataHist->FindLastBinAbove(sigDataHist->GetMaximum()/2.);
      
      fwhm= sigDataHist->GetBinCenter(bin2)- sigDataHist->GetBinCenter(bin1);
      maxMass=sigDataHist->GetBinCenter(sigDataHist->GetMaximumBin());     
      cout << "   Mass: " << mMC << " cat: " << cat << " intSM: " << sigSMEvents << " maxPDF " << maxMass << " fwhmPDF " << fwhm << endl;


     
      
      //    cout << "4. Costruisco la pdf Sig & BKG" << endl;

      Long64_t nBkgExp=8;
      nBkgExp= (Long64_t) data->sumEntries();
 
      RooRealVar bkgYield(Form("bkgYield_m%d_cat%d",mMC,cat),Form("bkgYield_m%d_cat%d",mMC,cat),nBkgExp,0,1e5);
      RooRealVar sigYield(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),muGen*sigSMEvents, -1000., 1e5);
  
      std::cout<<"-----> nBkgExp: "<<nBkgExp<<"   sig: "<<muGen*sigSMEvents<<std::endl;

      RooExtendPdf sigMCPdfE("sigMCPdfE", "sigMCPdfE", *sigMCPdf, sigYield);
      RooExtendPdf fitFcnE("fitFcnE", "fitFcnE", *fitFcn, bkgYield);
      
      // RooAddPdf sigAndBkg(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(fitFcnE,fitFcnE));
      RooAddPdf sigAndBkg(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(fitFcnE,sigMCPdfE));
      
      std::cout<< "after sum"<<std::endl;
  
      
      //plot
      TCanvas *tC = new TCanvas();
      tC->SetLogy(1);
  
      RooPlot *tPlot = mass->frame();
      RooPlot *tPlot2 = mass->frame();
      RooPlot *tPlot3 = mass->frame();
      RooPlot *tPlot4 = mass->frame();

     std::cout<< "before magenta plot "<<std::endl;
      sigRooData->plotOn(tPlot);
      sigMCPdf->plotOn(tPlot,LineColor(kMagenta));
      std::cout<< "after magenta plot "<<std::endl;
      tPlot->Draw();
      tC->Print(Form("SWplots/test_genSig_m%d_cat%d_toy%d.png",mMC, cat,0),"png");
  
      data ->plotOn(tPlot2);
      genFcn->plotOn(tPlot2,LineColor(kRed));
      tPlot2->Draw();
      tC->Print(Form("SWplots/test_gen%s_m%d_cat%d_toy%d.png",genName.c_str(),mMC,cat,0),"png");
     
      data ->plotOn(tPlot3);
      fitFcn->plotOn(tPlot3,LineColor(kBlue));
      tPlot3->Draw();
      RooArgSet* fitPar=fitFcn->getParameters(*mass);
      fitPar->Draw("same");
      tC->Print(Form("SWplots/test_fit%s_m%d_cat%d_toy%d.png",fitName.c_str(),mMC,cat,0),"png");

      
      data ->plotOn(tPlot4);
      sigAndBkg.plotOn(tPlot4,LineColor(kGreen));
      sigAndBkg.plotOn(tPlot4,Components(*(sigAndBkg.pdfList().at(1))),LineColor(kMagenta));
      
      tPlot4->GetYaxis()->SetRangeUser(0.00001, 100000);

      tPlot4->Draw();
    
      tC->Print(Form("SWplots/test_fitSandB%s_m%d_cat%d_toy%d.png", fitName.c_str(),mMC, cat,0),"png");



 

 std::cout<<"//--------------------------------------------------------------------------------//"<<std::endl;
 cout << "------------------------------------------------" << endl;
 cout << "--- Data fitted. Mass distributions obtained ---" << endl;
 cout << "------------------------------------------------" << endl;
 cout << "-------- Generating and fitting toys -----------" << endl;
 
 RooDataSet *genDatUnBin;
 
 RooDataHist *genDat;
 RooRandom::randomGenerator()->SetSeed(0);
 TRandom3 bkgGenRnd;
 RooRealVar sigGenEvents("sigGenEvents","sigGenEvents", 0,1e10);
 RooRealVar bkgGenEvents("bkgGenEvents","bkgGenEvents", 0,1e10);



  for (int itToy=jobNumb*toysPerJob; itToy<(jobNumb+1)*toysPerJob; itToy++){
    cout << "----------------------------------------------" << endl;
    cout << "------------- TOY: " << itToy << "------------" << endl;
    cout << "----------------------------------------------" << endl;

    //	std::cout<<"Genero un num di eventi di bkg pari a qnt sn i dati"<<std::endl;
      
#ifdef NOPOIS
      genDatUnBin = genFcn->generate(*mass,nBkgExp); // genero il bkg
#else
      genDatUnBin = genFcn->generate(*mass,nBkgExp,Extended());
#endif
      
     genDat = new RooDataHist("gen","gen",*mass,*genDatUnBin);

    
     RooDataSet combData(Form("combData_toy%d",itToy),Form("combData_toy%d",itToy),*mass,Import(*genDatUnBin));
     combData.Print();
     bkgGenEvents.setVal(combData.sumEntries());

     mass->setRange("sigWindow",maxMass-fwhm/2.,maxMass+fwhm/2.);
     mass->setRange("sigWindow2",maxMass-fwhm,maxMass+fwhm);
     RooAbsReal *intSigTrueRange= sigMCPdf->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
     RooFormulaVar *normIntSigTrueRange= new RooFormulaVar(Form("normIntSigRange_m%d_cat%d",mMC,cat),Form("normIntSigRange_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(*intSigTrueRange,sigYield));
     RooAbsReal *intSig2TrueRange= sigMCPdf->createIntegral(*mass,NormSet(*mass),Range("sigWindow2"));
     RooFormulaVar *normIntSig2TrueRange= new RooFormulaVar(Form("normIntSig2TrueRange_m%d_cat%d",mMC,cat),Form("normIntSig2TrueRange_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(*intSig2TrueRange,sigYield));
     
     if(sigGen){
       sigGenEvents.setVal(bkgGenRnd.Poisson(muGen*sigSMEvents));
       muTruth_=(sigGenEvents.getVal())/sigSMEvents;
       sigYield1Truefwhm_= normIntSigTrueRange->getVal();
       sigYield2Truefwhm_= normIntSig2TrueRange->getVal();	
       std::cout << "Generating " << sigGenEvents.getVal() << " signal events starting from " << muGen*sigSMEvents << std::endl;
       sigMCPdf->Print();
       RooDataSet *genSim=sigMCPdf->generate(*mass, sigGenEvents.getVal());
       genSim->Print();
       combData.append(*genSim);
       combData.Print();
     }else{
       muTruth_=0;
       sigYield1Truefwhm_=0;
       sigYield2Truefwhm_=0;
     } 
   
     //define fit PDF
     RooAddPdf *SandB;
     SandB = &sigAndBkg;
   

     // ---- fit and save output -----
     
     RooFitResult *res;
     if (onlyB) 
	   {
	     std::cout << "ONLY B FIT" << std::endl;
	     sigYield.setVal(0);
	     mu->setVal(0);
	     sigYield.setConstant(kTRUE);
	     mu->setConstant(kTRUE);
	   }
     
	 //	 res = SandB->fitTo(combData,Save(true), RooFit::Extended()); //
     
     RooNLLVar nll("nll","nll",*SandB,combData, kTRUE) ;	 
     RooMinuit m(nll) ;
     m.optimizeConst(false);
     res = m.fit("mhr");
     
     res->Print("v");
     
     
     cout << "------------------------------------------------" << endl;
     cout << "--- TOY: " << itToy << " Gen: " << genName << " Fit: " << fitName << " mass: " << mMC << endl;
     cout << "------------------------------------------------" << endl;
     res->floatParsFinal().Print("s");

     


     RooAbsReal *intRange = fitFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
     RooAbsReal *intBkgRange = genFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
     RooFormulaVar *normIntRange= new RooFormulaVar(Form("normIntRange_m%d_cat%d",mMC,cat),Form("normIntRange_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(*intRange,bkgYield));
    
     RooAbsReal *int2Range = fitFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow2"));
     RooAbsReal *int2BkgRange = genFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow2"));
     RooFormulaVar *normInt2Range= new RooFormulaVar(Form("normInt2Range_m%d_cat%d",mMC,cat),Form("normInt2Range_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(*int2Range,bkgYield));
     std::cout<<intRange->getVal()<<std::endl;
     std::cout << normIntRange->getVal() << " , "  << normIntRange->getPropagatedError(*res) << " ," << intBkgRange->getVal()* genDatUnBin->sumEntries() << std::endl;
     std::cout << normInt2Range->getVal() << " , "  << normInt2Range->getPropagatedError(*res) << " , " << int2BkgRange->getVal()*genDatUnBin->sumEntries() << std::endl;
     
     if (! onlyB)
       {
	 mu_=((RooRealVar *)(res->floatParsFinal().find(Form("sigYield_m%d_cat%d",mMC,cat))))->getVal()/(sigSMEvents*muGen);
	 sigma_mu_=((RooRealVar *)(res->floatParsFinal().find(Form("sigYield_m%d_cat%d",mMC,cat))))->getError()/(sigSMEvents*muGen)
	   ;
       }
     else
       {
	 mu_=0;
	 sigma_mu_=0;
       }




     RooAbsReal *intSigRange= sigMCPdf->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
     RooFormulaVar *normIntSigRange= new RooFormulaVar(Form("normIntSigRange_m%d_cat%d",mMC,cat),Form("normIntSigRange_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(*intSigRange,sigYield));
     RooAbsReal *intSig2Range= sigMCPdf->createIntegral(*mass,NormSet(*mass),Range("sigWindow2"));
     RooFormulaVar *normIntSig2Range= new RooFormulaVar(Form("normIntSig2Range_m%d_cat%d",mMC,cat),Form("normIntSig2Range_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(*intSig2Range,sigYield));

     sigYieldTot_ = ((RooRealVar *)(res->floatParsFinal().find(Form("sigYield_m%d_cat%d",mMC,cat))))->getVal();
     sigYieldTotErr_=  ((RooRealVar *)(res->floatParsFinal().find(Form("sigYield_m%d_cat%d",mMC,cat))))->getError();
     sigYield1fwhm_= normIntSigRange->getVal();
     sigYield2fwhm_= normIntSig2Range->getVal();
     sigYield1fwhmErr_= normIntSigRange->getPropagatedError(*res);
     sigYield2fwhmErr_= normIntSig2Range->getPropagatedError(*res);
     sigYield1fwhmErrNorm_= normIntSigRange->getPropagatedError(*res)/sigSMEvents ;
     sigYield2fwhmErrNorm_= normIntSig2Range->getPropagatedError(*res)/sigSMEvents ;

     std::cout<<" Signal Yield: "<< sigYieldTot_ <<" +/- "<<sigYieldTotErr_<<std::endl;
     std::cout<<" 1FWHM SIG: "<< sigYield1fwhm_<< " +/- " <<sigYield1fwhmErr_<< " ErrNorm: "<<sigYield1fwhmErrNorm_<<std::endl;
     std::cout<<" 2FWHM SIG: "<< sigYield2fwhm_<< " +/- " <<sigYield2fwhmErr_<< " ErrNorm: "<<sigYield2fwhmErrNorm_<<std::endl;


     bkgSig1fwhm_ = normIntRange->getVal();
     bkgSig2fwhm_ = normInt2Range->getVal();
     bkgTrue1fwhm_ =  intBkgRange->getVal()* genDatUnBin->sumEntries();
     bkgTrue2fwhm_ =  int2BkgRange->getVal()* genDatUnBin->sumEntries();
     bkgErrSig1fwhm_ = normIntRange->getPropagatedError(*res);
     bkgErrSig2fwhm_ = normInt2Range->getPropagatedError(*res);
     bkgErrNormSig1fwhm_ = normIntRange->getPropagatedError(*res)/sigSMEvents ;
     bkgErrNormSig2fwhm_ = normInt2Range->getPropagatedError(*res)/sigSMEvents ;	
     iToy_=itToy;
     cat_ = cat;
     
     
     fitStatus_=res->status();
     
     muTree->Fill();
     delete normIntRange;
     delete intRange;
     delete intBkgRange;
     delete normInt2Range;
     delete int2Range;
     delete int2BkgRange;
     
     outFile->cd();
     toyStep=5;
     
     res->SetName(Form("fitRes_g%s_f%s_m%d_t%d",genName.c_str(),fitName.c_str(),mMC,itToy));
     
     // --- find bkg integral -----
     
     
     if (doBkgInt){

       RooAbsReal *intInRange = fitFcn->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));

       double toBkg = intInRange->getVal()*genDat->numEntries();

       bias->Fill(toData-toBkg);
       
     }
   
     chi2_ = 999.;

     // ---- make plots ----
     if (itToy%toyStep==0 && itToy>=0){
       
       TCanvas *c1 = new TCanvas();
       RooPlot *mFrame = mass->frame(Title(Form("Gen: %s. Fit: %s. Mass %d  cat %d toy %d",genName.c_str(),fitName.c_str(),mMC,cat,itToy)));
       
       mass->setBins(200);
       
       combData.plotOn(mFrame,DataError(RooDataSet::SumW2));//,Binning("coarse"));
       
       if (plotGen) {
	 RooAddPdf genPdf(Form("Sig_Bkg_gen_m%d_cat%d",mMC,cat),Form("Sig_Bkg_gen_m%d_cat%d",mMC,cat),RooArgList(*genFcn,*sigMCPdf),RooArgList(bkgGenEvents,sigGenEvents));
	 if(sigGen)   genPdf.plotOn(mFrame,LineColor(kMagenta)); 
	 else	      genFcn->plotOn(mFrame,LineColor(kMagenta));
       }
       
       SandB->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed),Components(*(SandB->pdfList().at(0))));//,Range(100,160));
       SandB->plotOn(mFrame,LineColor(kGreen),LineStyle(kDashed),Components(*(SandB->pdfList().at(1))));
       SandB->plotOn(mFrame);
       
       
       std::cout<<"-------------------> CHI2: "<<mFrame->chiSquare(3)<<std::endl;
       chi2_=mFrame->chiSquare(3);
       mFrame->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
       mFrame->Draw();
       
       // ---- make legend -----
       TPaveText *text = new TPaveText(0.55,0.7,0.65,0.89,"NDC");
       text->SetLineColor(0);
       text->SetFillColor(0);
       text->AddText(Form("#mu = %.8f", mu_));
       leg->Draw("same");
       text->Draw("same");
       c1->SetLogy(0);
       if(sigGen)    c1->Print(Form("SWplots/toys/fitTogen-sigGen_%g-g%s_f%s_m%d_cat%d_toy%d.png",muGen, genName.c_str(),fitName.c_str(),mMC,cat,itToy),"png");
       else c1->Print(Form("SWplots/toys/fitTogen_g%s_f%s_m%d_cat%d_toy%d.png",genName.c_str(),fitName.c_str(),mMC,cat,itToy),"png");
	   
       c1->SetLogy(1);
       if(sigGen)    c1->Print(Form("SWplots/toys/fitTogen-sigGen_%g-g%s_f%s_m%d_cat%d_toy%d_LOG.png",muGen, genName.c_str(),fitName.c_str(),mMC,cat,itToy),"png");
       else c1->Print(Form("SWplots/toys/fitTogen_g%s_f%s_m%d_cat%d_toy%d_LOG.png",genName.c_str(),fitName.c_str(),mMC,cat,itToy),"png");
       
       
       delete c1;
       delete text;
       delete mFrame;
     }
  }//loop toys


 

  outFile->cd();
  muTree->Write();




  ofstream complete(Form("SWResults/g%s_f%s_j%d.txt",genName.c_str(),fitName.c_str(),jobNumb));
  complete << "Job: " << jobNumb << "/" << nJobs << " gen " << genName << " fit " << fitName << " completed successfully" << endl;
  complete.close();
  cout << "Text file written " << endl;
  outFile->Close();





  /*  delete fitFile;
  delete genFile;
  delete mu;
  delete outFile;
  delete muTree;
  delete leg;
  delete h;
  delete h1;
  delete h2;
  delete muHist;
  delete fitFcn;
  delete genFcn; 
  delete tPlot;
  delete tPlot2;
  delete tPlot3;
  delete tPlot4;
  delete sigMCPdf;
  */


  return 0;
}
