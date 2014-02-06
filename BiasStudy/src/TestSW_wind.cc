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
#include "TPad.h"
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
const int nFuncs=26;


string funcNames[nFuncs] = {"1pol","2pol","3pol","4pol","5pol","1exp","2exp","3exp","1pow","2powlow","3powlow","2lau","4lau","6lau","RooKey2","RooKey3","RooKey4","DiJet","DiJetPL","DiJetEXP","DiJetEXPOL",  "Expol","Exp","ExpolPL","Lau", "ExpPAR"};

void checkFunc(string name){
  if (name!="1pol" && name!="2pol" && name!="3pol" && name!="4pol" && name !="5pol" && name!="1exp" && name!="2exp" && name!="3exp" && name!="1pow" && name!="2powlow" && name!="3powlow" && name!="2lau" && name!="4lau" && name!="6lau" && name!="BkgMCKeyPdf"&& name!="DiJet"&& name!="DiJetPL"&& name!="DiJetEXP"&& name!="DiJetEXPOL"&& name!="Expol"&& name!="ExpolPL"&& name!="Exp"&& name!="ExpPAR"&& name!="Lau" && name!= "RooKey2"&& name!= "RooKey3"&& name!= "RooKey4") {
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

  if (name=="Exp" || name=="1pow" || name=="2lau"|| name=="Lau") return 1;
  else if (name=="1pol" ||name=="Expol") return 2;
  else if (name=="2pol" || name=="2exp" || name=="exp2_1" || name=="exp2_2" || name=="2pow" || name=="4lau" || name=="pow2_1" || name=="pow2_2" ||name=="DiJet"||name=="Expol2") return 3;
  else if (name=="ExpolPL") return 4;
  else if (name=="3pol"||name=="DiJetPL"||name=="DiJetEXP") return 5;
  else if (name=="DiJetEXPOL") return 6;
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

  else if (name=="DiJet") return "DiJet";
  else if (name=="DiJetPL") return "DiJetPL";
  else if (name=="DiJetEXP") return "DiJetEXP";
  else if (name=="DiJetEXPOL") return "DiJetEXPOL";
  else if (name=="Expol") return "Expol";
  else if (name=="ExpolPL") return "ExpolPL";
  else if (name=="Exp") return "Exp";
  else if (name=="Lau") return "Lau";
  else if (name=="ExpPAR") return "ExpPAR";
  else if (name== "RooKey2") return "RooKey2";
  else if (name== "RooKey3") return "RooKey3";
  else if (name== "RooKey4") return "RooKey4";
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
  float nSigGen=1;
  int nToys;
  int nJobs=1;
  int jobNumb=0;
  int toyStep=2;
  string genName;
  string fitName;

  int mMC;
  int width;
  int cat;
  int winSize;
  string fileFitName;
  string fileGenName;
  string data_fileName;

  string suffix;


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
    if (string(argv[arg])=="-nGen")  nSigGen=atof(argv[arg+1]);
    if (string(argv[arg])=="-mass") mMC=atoi(argv[arg+1]);
    if (string(argv[arg])=="-width") width=atoi(argv[arg+1]);
    if (string(argv[arg])=="-cat")  cat=atoi(argv[arg+1]);
    if (string(argv[arg])=="-w") winSize=atoi(argv[arg+1]);
    if (string(argv[arg])=="-su") suffix= string(argv[arg+1]);

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
    cout << "    -w    to specify pm win " << endl;
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

  if(sigGen==true && nSigGen>0){
    if (mMC==150)nSigGen=842;
    if (mMC==250)nSigGen=311;  
    if (mMC==350)nSigGen=180;  
    if (mMC==450)nSigGen=115;      
    if (mMC==550)nSigGen=80;  
    if (mMC==650)nSigGen=56;	  
    if (mMC==750)nSigGen=44;	  
    if (mMC==850)nSigGen=36;   
  }
 
  string genType=getType(genName);
  string fitType=getType(fitName);
  
  cout << "--- Running with following options ---" << endl;
  cout << "    nToys:            " << nToys << endl;
  cout << "    plotStep:         " << toyStep << endl;
  cout << "    genFunction:      " << genName << endl;
  cout << "    fitFunction:      " << fitName << endl;
  cout << "    mass:             " << mMC << endl;
  cout << "    width:             " << width << endl;
  cout << "    winSize:          " << winSize << endl;
  cout << "    nJobs:            " << nJobs << endl;
  cout << "    jobNum:           " << jobNumb << endl;
  cout << "    cat:           " << cat<< endl;
  
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
  
  
 
  //  std::cout<<"Apro il file di input e prendo il WS che si chiama wBias e salvo la variabile della massa"<<std::endl;
  
 
  
  TFile *genFile = new TFile(fileGenName.c_str());

  std::cout<<"GENWS: "<<fileGenName.c_str()<<std::endl;
  RooWorkspace *genWS;
  
  genWS = (RooWorkspace*)genFile->Get("wBiasTruth");
  if(genWS == NULL){
    std::cerr << "[ERROR] genWS not found in WS" << std::endl;
    exit(1);
  }






 RooRealVar *mass=(RooRealVar*)genWS->var("PhotonsMass");
 if(mass == NULL){
   std::cerr << "[ERROR] variable mass not found in WS: " << genWS->GetName() << std::endl;
   return 1;
 }

 std::cout<<"//------------------------------ Range MASS --------------------------------------------------//"<<std::endl;
 
 Float_t mLow =genWS->var("RooMINmass")->getVal();
 Float_t mHigh = genWS->var("RooMAXmass")->getVal();
 
 std::cout<<"mLOW: "<< mLow << "   mHIGH: "<<mHigh<<std::endl;
 mass->setRange("winRange",mLow,mHigh);
 mass->setRange(mLow,mHigh);
  


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
 RooAbsPdf *genFcn2;
 RooAbsPdf *genFcn3;
 RooAbsPdf *genParamatersFcn;
 
 RooAbsPdf *fitFcn;

 //par Exp
 RooRealVar p1Exp("p1Exp", "p1Exp", -0.1, -5., 0.);

 //par Expol 
 RooRealVar p1Expol("p1Expol","p1Expol", 5.7, 0.,50.);
 RooRealVar p2Expol("p2Expol","p2Expol", 0.000005, 0., 10.);

 std::cout<<"after expol"<<std::endl;
 genWS->Print("V");
//par ExpolPL 
 RooRealVar p1ExpolPL("p1ExpolPL","p1ExpolPL",21.7,0., 50.);
 RooRealVar p2ExpolPL("p2ExpolPL","p2ExpolPL",0.05, 0.,10);
 RooRealVar p3ExpolPL("p3ExpolPL","p3ExpolPL", 0.8, 0., 1.);//frac
 RooRealVar p4ExpolPL("p4ExpolPL","p4ExpolPL",-0.8, -100.,0.);

 std::cout<<"after expolPL"<<std::endl;
 //par DiJet
 RooRealVar p1DiJet("p1DiJet","p1DiJet", 5., 0., 20.);
 RooRealVar p2DiJet("p2DiJet","p2DiJet",0.000002, 0., 0.00005);
 RooRealVar p3DiJet("p3DiJet","p3DiJet",0.24, 0., 1.);
 std::cout<<"after dijet"<<std::endl;

 //par DiJetPL
 RooRealVar p1DiJetPL("p1DiJetPL","p1DiJetPL",6.5,1., 20.);
 RooRealVar p2DiJetPL("p2DiJetPL","p2DiJetPL",0.0000002, 0., 50.);
 RooRealVar p3DiJetPL("p3DiJetPL","p3DiJetPL",0.36, 0., 10.);
 RooRealVar p4DiJetPL("p4DiJetPL","p4DiJetPL", 0.8, 0., 1.);//frac
 RooRealVar p5DiJetPL("p5DiJetPL","p5DiJetPL",-0.02, -100., 0.);


 std::cout<<"after dijetPL"<<std::endl;

//par DiJetEXP
 RooRealVar p1DiJetEXP("p1DiJetEXP","p1DiJetEXP",6.52,1., 20);
 RooRealVar p2DiJetEXP("p2DiJetEXP","p2DiJetEXP",0.00000045, 0., 50.);
 RooRealVar p3DiJetEXP("p3DiJetEXP","p3DiJetEXP",0.36, 0., 10);
 RooRealVar p4DiJetEXP("p4DiJetEXP","p4DiJetEXP",0.8, 0., 1.);//frac
 RooRealVar p5DiJetEXP("p5DiJetEXP","p5DiJetEXP",-0.02, -100., 0.);

//par DiJetEXP
 RooRealVar p1DiJetEXPOL("p1DiJetEXPOL","p1DiJetEXPOL",6.52,1., 20);
 RooRealVar p2DiJetEXPOL("p2DiJetEXPOL","p2DiJetEXPOL",0.00000045, 0., 50.);
 RooRealVar p3DiJetEXPOL("p3DiJetEXPOL","p3DiJetEXPOL",0.36, 0., 10);
 RooRealVar p4DiJetEXPOL("p4DiJetEXPOL","p4DiJetEXPOL",0.8, 0., 1.);//frac
 RooRealVar p5DiJetEXPOL("p5DiJetEXPOL","p5DiJetEXPOL",21.7,0., 50.);
 RooRealVar p6DiJetEXPOL("p6DiJetEXPOL","p6DiJetEXPOL",0.05, 0.,10);



 //par Lau
 RooRealVar p1Lau("p1Lau", "p1Lau", 0.99878, 0.8, 1);
 // RooRealVar p2Lau("p2Lau", "p2Lau", 0.0000008, 0., 0.0001);


 //par ExpPAR
 RooRealVar p1ExpPAR("p1ExpPAR", "p1ExpPAR", 0.006, 0.00004, 1);
 RooRealVar p2ExpPAR("p2ExpPAR", "p2ExpPAR", -2., -6., 0.);



  std::cout<<"after dijetEXPOL"<<std::endl;

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
 Float_t    sigYieldTotTruth_;
 Float_t    sigYieldTotErr_;
 

 Int_t iToy_;
 Int_t migradStatus_;
 Int_t fitStatus_;
 Int_t cat_;
 Int_t mMC_;
 Int_t winSize_;
 Float_t chi2_;
 char genName_[30];
 char fitName_[30];

 Float_t fracDiJetEXPOL_;
 Float_t fracDiJetEXP_;
 Float_t fracDiJetPL_;
 Float_t fracExpolPL_;

 Float_t par1Lau_;
 Float_t par2Lau_;
 
 sprintf(genName_, "%s", genName.c_str());
 sprintf(fitName_, "%s", fitName.c_str());
 
 // std::cout << "getName = " << genName_ << std::endl;
 // std::cout << "fitName = " << fitName_ << std::endl;
 
 //std::cout << "[INFO] creating new output file" << std::endl;
 TFile *outFile;
 
 if(!sigGen)  outFile = new TFile(Form("TestSWResults/biasCheck_g%s_f%s_m%d_w%d_cat%d_w%d_j%d_%s.root",genName.c_str(),fitName.c_str(),mMC,width,cat,winSize,jobNumb, suffix.c_str()),"RECREATE");
 else outFile = new TFile(Form("TestSWResults/biasCheck-sigGen_%g-g%s_f%s_m%d_w%d_cat%d_w%d_j%d_%s.root",nSigGen, genName.c_str(),fitName.c_str(),mMC,width,cat,winSize,jobNumb,suffix.c_str()),"RECREATE");
 
 TTree *muTree = new TTree("muTree","muTree");
 muTree->Branch("mass", &mMC, "mass/I");
 muTree->Branch("width", &width, "width/I");
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
 muTree->Branch("sigYieldTotTruth", &sigYieldTotTruth_, "sigYieldTotTruth/F");
 muTree->Branch("sigYieldTotErr", &sigYieldTotErr_, "sigYieldTotErr/F");
 


 muTree->Branch("genFun", genName_, "genFun/C");
 muTree->Branch("fitFun", fitName_, "fitFun/C");
 muTree->Branch("iToy", &iToy_, "iToy/i");
 muTree->Branch("cat", &cat_, "cat/i");
 muTree->Branch("mMC", &mMC_, "mMC/i");
 muTree->Branch("winSize", &winSize_, "winSize/i");
 muTree->Branch("chi2", &chi2_, "chi2/F");
 muTree->Branch("migradStatus",&migradStatus_, "migradStatus/I");
 muTree->Branch("fitStatus",&fitStatus_, "fitStatus/I");
 
 
 muTree->Branch("fracDiJetEXPOL", &fracDiJetEXPOL_, "fracDiJetEXPOL/F");
 muTree->Branch("fracDiJetEXP", &fracDiJetEXP_, "fracDiJetEXP/F");
 muTree->Branch("fracDiJetPL", &fracDiJetPL_, "fracDiJetPL/F");
 muTree->Branch("fracExpolPL", &fracExpolPL_, "fracExpolPL/F");

 muTree->Branch("par1Lau",&par1Lau_, "par1Lau/F");
 muTree->Branch("par2Lau",&par2Lau_, "par2Lau/F");
 
 TLegend *leg = new TLegend(0.65,0.7,0.89,0.89);
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
 TH1F *h3 = new TH1F("h3","h3",1,0,1);
 h3->SetLineColor(kOrange);
 h3->SetLineStyle(kDashed);
 h3->SetLineWidth(3);
 TH1F *h4 = new TH1F("h4","h4",1,0,1);
 h4->SetLineColor(kSpring-9);
 h4->SetLineStyle(kDashed);
 h4->SetLineWidth(3);
 if (plotGen) leg->AddEntry(h,"Truth","l");
 leg->AddEntry(h2,"bkg part of s+b","l");
 leg->AddEntry(h1,"s+b after fit","l");
 leg->AddEntry(h3,"First Bkg Component","l");
 leg->AddEntry(h4,"Second Bkg Component","l");
 
 
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
 if(genName == "RooKey2") genFcn = (RooKeysPdf*) genWS->pdf(Form("BkgMCKeyPdf_bw2_cat%d",cat));
 if(genName == "RooKey3") genFcn = (RooKeysPdf*) genWS->pdf(Form("BkgMCKeyPdf_bw3_cat%d",cat));
 if(genName == "RooKey4") genFcn = (RooKeysPdf*) genWS->pdf(Form("BkgMCKeyPdf_bw4_cat%d",cat));
 
 if(genName =="DiJet") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_DiJet_truth_cat%d",cat));
 if(genName =="DiJetPL") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_DiJetPL_truth_cat%d",cat));
 if(genName =="DiJetEXP") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_DiJetEXP_truth_cat%d",cat));
 if(genName =="DiJetEXPOL") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_DiJetEXPOL_truth_cat%d",cat));
 if(genName =="Expol") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_Expol_truth_cat%d",cat));
 if(genName =="ExpolPL") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_ExpolPL_truth_cat%d",cat));
 if(genName =="Lau") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_Lau_truth_cat%d",cat));
 if(genName =="ExpPAR") genFcn = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_ExpPAR_truth_cat%d",cat));
  
 genFcn2 = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_DiJet_truth_cat%d",cat));
 genFcn3 = (RooGenericPdf*)genWS->pdf(Form("PhotonsMassBkg_Lau_truth_cat%d",cat));
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
	if(fitName == "RooKey4"){
	  fitFcn = (RooKeysPdf*) genWS->pdf(Form("BkgMCKeyPdf_bw4_cat%d",cat));
	  fitFcn->SetName("rookey");
	  fitFcn->SetTitle("rookey");
	}
	if(fitName =="Exp")fitFcn = new RooGenericPdf("fitFcn_exp", "exp((@0*@1))", RooArgList(*mass, p1Exp));
	if(fitName =="Expol")fitFcn = new RooGenericPdf("fitFcn_expol", "exp(-@0/8000./(@1+@2*@0/8000.))", RooArgList(*mass, p1Expol, p2Expol));

	if(fitName =="ExpolPL"){
	  RooAbsPdf* Expol = new RooGenericPdf("EXPOL", "exp(-@0/8000./(@1+@2*@0/8000.))", RooArgList(*mass, p1ExpolPL, p2ExpolPL ));   
	  RooAbsPdf* PLE = new RooGenericPdf("PLE", "pow(@0,@1)", RooArgList(*mass,  p4ExpolPL));
	  fitFcn = new RooAddPdf("fitFcn_expolPL", "fitFcn_expolPL",RooArgList(*Expol, *PLE), RooArgList( p3ExpolPL));
	}
	if(fitName =="DiJet") fitFcn = new RooGenericPdf("fitFcn_dijet","fitFcn_dijet", "pow((1-(@0/8000.)), @2)/pow((@0/8000.), @1+@3*log(@0/8000.))", RooArgList(*mass, p1DiJet, p2DiJet, p3DiJet));

	if(fitName =="DiJetPL"){
	  RooAbsPdf* DiJet = new RooGenericPdf("DIJET", "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*mass, p1DiJetPL, p2DiJetPL, p3DiJetPL ));  
	  RooAbsPdf* PL = new RooGenericPdf("PL", "pow(@0,@1)", RooArgList(*mass,  p5DiJetPL));
	  fitFcn = new RooAddPdf("fitFcn_dijetPL","fitFcn_dijetPL",RooArgList(*DiJet, *PL), RooArgList( p4DiJetPL));
	}
	if(fitName =="DiJetEXP"){
	  RooAbsPdf* DiJetE = new RooGenericPdf("DIJETE", "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*mass, p1DiJetEXP, p2DiJetEXP, p3DiJetEXP )); 
	  RooExponential* EXP = new RooExponential("EXP","EXP", *mass,  p5DiJetEXP); 
	  fitFcn = new RooAddPdf("fitFcn_dijetEXP","fitFcn_dijetEXP",RooArgList(*DiJetE, *EXP), RooArgList( p4DiJetEXP));

	}if(fitName =="DiJetEXPOL"){
	  RooAbsPdf* DiJetEx = new RooGenericPdf("DIJETEx", "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*mass, p1DiJetEXPOL, p2DiJetEXPOL, p3DiJetEXPOL )); 
	  RooAbsPdf* ExpolDiJ = new RooGenericPdf("EXPOLDiJ", "exp(-@0/8000./(@1+@2*@0/8000.))", RooArgList(*mass,p5DiJetEXPOL, p6DiJetEXPOL  ));   
	  fitFcn = new RooAddPdf("fitFcn_dijetEXPOL","fitFcn_dijetEXPOL",RooArgList(*DiJetEx, *ExpolDiJ), RooArgList( p4DiJetEXPOL));

	}
	if(fitName =="Lau")fitFcn = new RooGenericPdf("fitFcn_lau", "(1-@1)*pow(@0,-4.0)+@1*pow(@0,-5.0)", RooArgList(*mass, p1Lau));
	if(fitName =="ExpPAR")fitFcn = new RooGenericPdf("fitFcn_expPAR", "exp(-@1*@0)*pow(@0, @2)", RooArgList(*mass, p1ExpPAR, p2ExpPAR));
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

      cout << "3. Prendo la pdf Sig dal WS" << endl;
      RooFFTConvPdf* sigMCPdf = (RooFFTConvPdf*) genWS->pdf(TString::Format("PhotonsMassSig_cat%d",cat));
	if(sigMCPdf==NULL){
	  std::cout << "[ERROR] sigMCPdf not defined" << std::endl;
	  return 1;
	}

     
      
      cout << "4. Costruisco la pdf Sig & BKG" << endl;

      Long64_t nBkgExp=8;
      nBkgExp= (Long64_t) data->sumEntries(Form("PhotonsMass>%f&&PhotonsMass<%f", mLow, mHigh));

      std::cout<<nBkgExp<<std::endl;

      Float_t nSigExp = nSigGen;
      RooRealVar bkgYield(Form("bkgYield_m%d_cat%d",mMC,cat),Form("bkgYield_m%d_cat%d",mMC,cat),nBkgExp,nBkgExp*0.92, nBkgExp*1.08);//8%
      RooRealVar sigYield(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),nSigExp,-10000., 10000 );//+-100
  
      std::cout<<"-----> nBkgExp: "<<nBkgExp<<"   nSigExp: "<<nSigExp<<std::endl;

      RooExtendPdf sigMCPdfE("sigMCPdfE", "sigMCPdfE", *sigMCPdf, sigYield);
      RooExtendPdf fitFcnE("fitFcnE", "fitFcnE", *fitFcn, bkgYield);
      
      // RooAddPdf sigAndBkg(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(fitFcnE,fitFcnE));
      RooAddPdf sigAndBkg(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(fitFcnE,sigMCPdfE));
      
      std::cout<< "after sum"<<std::endl;
  
      
      //plot
      TCanvas *tC = new TCanvas();
    
  
      RooPlot *tPlot = mass->frame();
      RooPlot *tPlot2 = mass->frame();
      
      RooPlot *tPlot4 = mass->frame();

      
      sigMCPdf->plotOn(tPlot,LineColor(kMagenta));
      tPlot->Draw();
      tC->SetLogy(1);
      tC->Print(Form("SWplots/test_genSig_m%d_w%d_cat%d_toy%d_%s_LOG.png",mMC, width,cat,0, suffix.c_str()),"png");
  
      data ->plotOn(tPlot2);
      genFcn->plotOn(tPlot2,LineColor(kRed));
      tPlot2->Draw();
      tC->Print(Form("SWplots/test_gen%s_m%d_w%d_cat%d_toy%d_%s_LOG.png",genName.c_str(),mMC,width,cat,0, suffix.c_str()),"png");
      tC->SetLogy(0);
      tC->Print(Form("SWplots/test_gen%s_m%d_w%d_cat%d_toy%d_%s.png",genName.c_str(),mMC,width,cat,0, suffix.c_str()),"png");
      tC->Clear();

      //-------pad 1-------//
      TPad * pad1 = new TPad("pad1", "pad1",0.01,0.14,1,1.);  
      pad1->SetRightMargin(0.1);
      pad1->Range(0.01, 0.15, 1., 1);
      pad1->Draw();       
      pad1->cd();
      pad1->SetLogy();
      RooPlot *tPlot3 = mass->frame();
      data ->plotOn(tPlot3, RooFit::Invisible());
      
      fitFcn->plotOn(tPlot3,LineColor(kBlue));
      genFcn->plotOn(tPlot3,LineColor(kRed));
      //      TH1F* h_gen3 = (TH1F*)genFcn3->createHistogram("PhotonsMass", 200);
      TH1F* h_fit = (TH1F*)genFcn2->createHistogram("PhotonsMass", 200);
      TH1F* h_gen = (TH1F*)genFcn->createHistogram("PhotonsMass", 200);


      h_fit->SetLineColor(kBlue);
      h_gen->SetLineColor(kRed);
   

      tPlot3->Draw();
      TLegend* l = new TLegend(0.5, 0.6, 0.8, 0.8);
      //       l->AddEntry(h_gen3, "Lau", "L");
       l->AddEntry(h_fit, "Lau", "L");
       l->AddEntry(h_gen, "Lau", "L");
       l->SetFillColor(kWhite);
       l->SetBorderSize(0);
        l->Draw("same");

      //    RooArgSet* fitPar=fitFcn->getParameters(*mass);
      //    fitPar->Draw("same");

      //-------pad 2------//
      tC->cd();
      TPad * pad2 = new TPad("pad2", "pad2",0.01,0.01,1.,0.13);
       pad2->SetGrid();
      
       pad2->Draw();
       pad2->cd();
      
      
      
       TH1F* h1_ratio=(TH1F*) h_fit->Clone("h1_ratio");
       h1_ratio->SetLineColor(kBlack);
       Double_t xmax = h1_ratio->GetXaxis()->GetXmax();
       Double_t xmin = h1_ratio->GetXaxis()->GetXmin();
       TLine* line = new TLine(xmin,1.,xmax,1.);
       
       
       h1_ratio->SetStats(0);
       h1_ratio->Divide(h_gen);
       h1_ratio->GetYaxis()->SetRangeUser(0.8, 1.2);
       h1_ratio->GetYaxis()->SetNdivisions(2,false);
       line->SetLineColor(kRed);
       h1_ratio->SetLineWidth(2);
       h1_ratio->Draw("hist");
       line->Draw("same");
      




      tC->Print(Form("SWplots/test_fit%s_m%d_w%d_cat%d_toy%d_%s.png",fitName.c_str(),mMC,width,cat,0, suffix.c_str()),"png");
      tC->SaveAs(Form("SWplots/test_fit%s_m%d_w%d_cat%d_toy%d_%s.C",fitName.c_str(),mMC,width,cat,0, suffix.c_str()));
      tC->SetLogy(1);
      tC->Print(Form("SWplots/test_fit%s_m%d_w%d_cat%d_toy%d_%s_LOG.png",fitName.c_str(),mMC,width,cat,0, suffix.c_str()),"png");
      //************//
  
      data ->plotOn(tPlot4);
      sigAndBkg.plotOn(tPlot4,LineColor(kGreen));
      sigAndBkg.plotOn(tPlot4,Components(*(sigAndBkg.pdfList().at(1))),LineColor(kMagenta));
      
      tPlot4->GetYaxis()->SetRangeUser(0.00001, 100000);

      tPlot4->Draw();
    
      tC->Print(Form("SWplots/test_fitSandB%s_m%d_w%d_cat%d_toy%d_%s_LOG.png", fitName.c_str(),mMC,width, cat,0, suffix.c_str()),"png");
      tC->SetLogy(0);
      tC->Print(Form("SWplots/test_fitSandB%s_m%d_w%d_cat%d_toy%d_%s.png", fitName.c_str(),mMC,width, cat,0, suffix.c_str()),"png");


 

 std::cout<<"//--------------------------------------------------------------------------------//"<<std::endl;
 cout << "------------------------------------------------" << endl;
 cout << "--- Data fitted. Mass distributions obtained ---" << endl;
 cout << "------------------------------------------------" << endl;
 cout << "-------- Generating and fitting toys -----------" << endl;
 
 RooDataSet *genDatUnBin;
 
 RooDataHist *genDat;
 RooRandom::randomGenerator()->SetSeed(0);
 TRandom3 bkgGenRnd;
 TRandom3 sigGenRnd;
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

    
     if(sigGen && nSigExp>0){
       sigGenEvents.setVal(nSigExp);
       muTruth_=(sigGenEvents.getVal())/(sigSMEvents*muGen);
       sigYieldTotTruth_=sigGenEvents.getVal();
       std::cout << "Generating " << sigGenEvents.getVal() << " signal events starting from " << nSigExp << std::endl;
       sigMCPdf->Print();
       RooDataSet *genSim=sigMCPdf->generate(*mass, sigGenEvents.getVal());
       genSim->Print();
       combData.append(*genSim);
       combData.Print();
     }else{
       muTruth_=0;
       sigYieldTotTruth_=0;
     } 
   
     //define fit PDF
     RooArgSet* paramsBkg = fitFcn->getParameters(*mass);
     paramsBkg->readFromFile("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/initParBkg.txt");
     paramsBkg->writeToStream(std::cout, kFALSE);

     RooAddPdf *SandB;
     sigYield.setVal(nSigExp);
     bkgYield.setVal(nBkgExp);
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
     m.migrad();
     RooFitResult* r = m.save() ; 
     migradStatus_=r->status();
     m.hesse() ; 
      res = m.save() ;

     // res = m.fit("mhr");
     
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


     sigYieldTot_ = ((RooRealVar *)(res->floatParsFinal().find(Form("sigYield_m%d_cat%d",mMC,cat))))->getVal();
     sigYieldTotErr_=  ((RooRealVar *)(res->floatParsFinal().find(Form("sigYield_m%d_cat%d",mMC,cat))))->getError();
    

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
     
     fracDiJetEXPOL_ = 0;
     fracDiJetEXP_ = 0;
     fracDiJetPL_ = 0;
     fracExpolPL_ = 0;
     par1Lau_ = 0;
     if(fitName=="DiJetEXPOL"){
       fracDiJetEXPOL_ =((RooRealVar *)(res->floatParsFinal().find("p4DiJetEXPOL")))->getVal();
     }else if(fitName=="DiJetEXP"){
       fracDiJetEXP_ =((RooRealVar *)(res->floatParsFinal().find("p4DiJetEXP")))->getVal();
     }else if(fitName=="DiJetPL"){
       fracDiJetPL_ =((RooRealVar *)(res->floatParsFinal().find("p4DiJetPL")))->getVal();
     }else if(fitName=="ExpolPL"){
       fracExpolPL_ =((RooRealVar *)(res->floatParsFinal().find("p3ExpolPL")))->getVal();
     }else if(fitName=="Lau"){
       par1Lau_ = ((RooRealVar *)(res->floatParsFinal().find("p1Lau")))->getVal();
       //   par2Lau_ = ((RooRealVar *)(res->floatParsFinal().find("p2Lau")))->getVal();
     }

     fitStatus_=res->status();
     

     //save chi2 for all toys
     RooPlot *cFrame = mass->frame(Title(Form("Gen: %s. Fit: %s. Mass %d  Width: %d cat %d winSize: %d toy: %d range: %s",genName.c_str(),fitName.c_str(),mMC,width, cat,winSize,itToy, suffix.c_str())));
     combData.plotOn(cFrame,DataError(RooDataSet::SumW2));
     SandB->plotOn(cFrame);
     chi2_ = cFrame->chiSquare(3);


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
   
     //chi2_ = 999.;

     // ---- make plots ----
     if (itToy%toyStep==0 && itToy>=0){
       
       TCanvas *c1 = new TCanvas("c1", "c1",1);
       c1->cd();
       //-------pad 1-------//
       TPad * pad1 = new TPad("pad1", "pad1",0.01,0.25,1,1.);  
       pad1->SetRightMargin(0.1);
       pad1->Draw();       
       pad1->cd();

       RooPlot *mFrame = mass->frame(Title(Form("Gen: %s. Fit: %s. Mass %d  Width: %d cat %d winSize: %d toy: %d range: %s",genName.c_str(),fitName.c_str(),mMC,width, cat,winSize,itToy, suffix.c_str())));
       
       mass->setBins(200);
       
       combData.plotOn(mFrame,DataError(RooDataSet::SumW2));//,Binning("coarse"));
       
       if (plotGen) {
	 RooAddPdf genPdf(Form("Sig_Bkg_gen_m%d_cat%d",mMC,cat),Form("Sig_Bkg_gen_m%d_cat%d",mMC,cat),RooArgList(*genFcn,*sigMCPdf),RooArgList(bkgGenEvents,sigGenEvents));
	 if(sigGen)   genPdf.plotOn(mFrame,LineColor(kMagenta)); 
	 else	      genFcn->plotOn(mFrame,LineColor(kMagenta));
       }
       


       SandB->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed),Components(*(SandB->pdfList().at(0))));
       SandB->plotOn(mFrame,LineColor(kGreen),LineStyle(kDashed),Components(*(SandB->pdfList().at(1))));
       SandB->plotOn(mFrame);
        if(fitName=="DiJetEXPOL"){
	
	 fitFcn->plotOn(mFrame,LineColor(kOrange),Components("DIJETEx"),LineStyle(kDashed));
	 fitFcn->plotOn(mFrame,LineColor(kSpring-9),Components("EXPOLDiJ"),LineStyle(kDashed));

       }else if(fitName=="DiJetEXP"){
	
	  fitFcn->plotOn(mFrame,LineColor(kOrange),Components("DIJETE"),LineStyle(kDashed));	
	  fitFcn->plotOn(mFrame,LineColor(kSpring-9),Components("EXP"),LineStyle(kDashed));

       }else if(fitName=="DiJetPL"){

	 fitFcn->plotOn(mFrame,LineColor(kOrange),Components("DIJET"),LineStyle(kDashed));
	 fitFcn->plotOn(mFrame,LineColor(kSpring-9),Components("PL"),LineStyle(kDashed));

       }else if(fitName=="ExpolPL"){

	 fitFcn->plotOn(mFrame,LineColor(kOrange),Components("EXPOL"),LineStyle(kDashed));
	 fitFcn->plotOn(mFrame,LineColor(kSpring-9),Components("PLE"),LineStyle(kDashed));

       }
      std::cout<<"-------------------> CHI2: "<<mFrame->chiSquare(3)<<std::endl;
      //  chi2_=mFrame->chiSquare(3);
       mFrame->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
       mFrame->GetYaxis()->SetRangeUser(0.001, 100000);
       mFrame->Draw();
       
       // ---- make legend -----
       TPaveText *text = new TPaveText(0.57,0.69,0.65,0.89,"NDC");
       text->SetLineColor(0);
       text->SetFillColor(0);
       //text->AddText(Form("#mu = %.8f", mu_));
       leg->Draw("same");
       //text->Draw("same");
      
 
      

       pad1->SetLogy(0);
       if(sigGen)    c1->Print(Form("SWplots/toys/fitTogen-sigGen_%g-g%s_f%s_m%d_w%d_cat%d_w%d_toy%d_%s.png",nSigGen, genName.c_str(),fitName.c_str(),mMC,width,cat,winSize,itToy, suffix.c_str()),"png");
       else c1->Print(Form("SWplots/toys/fitTogen_g%s_f%s_m%d_w%d_cat%d_toy%d_%s.png",genName.c_str(),fitName.c_str(),mMC,width,cat,itToy, suffix.c_str()),"png");
	   
       pad1->SetLogy(1);
       if(sigGen)    c1->Print(Form("SWplots/toys/fitTogen-sigGen_%g-g%s_f%s_m%d_w%d_cat%d_w%d_toy%d_%s_LOG.png",nSigGen, genName.c_str(),fitName.c_str(),mMC,width,cat,winSize,itToy, suffix.c_str()),"png");
       else c1->Print(Form("SWplots/toys/fitTogen_g%s_f%s_m%d_w%d_cat%d_w%d_toy%d_%s_LOG.png",genName.c_str(),fitName.c_str(),mMC,width,cat,itToy, suffix.c_str()),"png");
    



   
       
       
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
