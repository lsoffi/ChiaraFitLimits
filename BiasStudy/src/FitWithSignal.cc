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

using namespace std;
using namespace RooFit;

#define oldFitFunDefinition
#define PTL45_PTSL25_MET70_2012_FP

const int nFuncs=14;
string funcNames[nFuncs] = {"1pol","2pol","3pol","4pol","5pol","1exp","2exp","3exp","1pow","2pow","3pow","2lau","4lau","6lau"};

void checkFunc(string name){
  if (name!="1pol" && name!="2pol" && name!="3pol" && name!="4pol" && name !="5pol" && name!="1exp" && name!="2exp" && name!="3exp" && name!="1pow" && name!="2pow" && name!="3pow" && name!="2lau" && name!="4lau" && name!="6lau") {
    cout << "Invalid function: " << name << endl;
    cout << "Options are: " << endl;
    for (int f=0; f<nFuncs; f++){
      cout << "   " << funcNames[f].c_str() << endl;
    }
    exit(1);
  }
}

void checkMass(int mL){
  if (mL!=110 && mL!=115 && mL!=120 && mL!=125  && mL!=130 && mL!=135 && mL!=140 && mL!=150) cout << "Invalid mass: " << mL << endl; 
  assert(mL==110 || mL==115 || mL==120 || mL==125  || mL==130 || mL==135 || mL==140 || mL==150);
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
  else if (name=="1pol") return 2;
  else if (name=="2pol" || name=="2exp" || name=="exp2_1" || name=="exp2_2" || name=="2pow" || name=="4lau" || name=="pow2_1" || name=="pow2_2") return 3;
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
  else if (name=="1pow" || name=="2pow" || name=="3pow" || name=="pow2_1" || name=="pow2_2") return "pow";
  else if (name=="2lau" || name=="4lau" || name=="6lau") return "lau";
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
  float muGen=1;
  int nToys;
  int nJobs=1;
  int jobNumb=0;
  int toyStep;
  string genName;
  string fitName;
  int winSize=5;
  int mL=110;
  string fileName;
  string data_fileName;


  std::vector<TString> genFun_vec;
  std::vector<TString> fitFun_vec;
  std::vector<float> mass_vec;
  mass_vec.push_back(110);
  mass_vec.push_back(115);
  mass_vec.push_back(120);
  mass_vec.push_back(130);
  mass_vec.push_back(140);


#ifdef VBF
  std::map<int, double> expEvents_gluglu_map;
  expEvents_gluglu_map[110]=0.752;
  expEvents_gluglu_map[115]=0.740;
  expEvents_gluglu_map[120]=0.761;
  expEvents_gluglu_map[125]=0.790;
  expEvents_gluglu_map[130]=0.751;
  expEvents_gluglu_map[135]=0.697;
  expEvents_gluglu_map[140]=0.570;
  expEvents_gluglu_map[145]=0.550;
  expEvents_gluglu_map[150]=0.412;
  
  
  std::map<int, double> expEvents_vbf_map;
  expEvents_vbf_map[110]=1.762;
  expEvents_vbf_map[115]=1.890;
  expEvents_vbf_map[120]=2.003;
  expEvents_vbf_map[125]=1.957;
  expEvents_vbf_map[130]=1.914;
  expEvents_vbf_map[135]=1.769;
  expEvents_vbf_map[140]=1.588;
  expEvents_vbf_map[145]=1.322;
  expEvents_vbf_map[150]=1.079;
  
  std::map<int, double> expEvents_wzh_map;
  expEvents_wzh_map[110]=0.018;
  expEvents_wzh_map[115]=0.020;
  expEvents_wzh_map[120]=0.019;
  expEvents_wzh_map[125]=0.017;
  expEvents_wzh_map[130]=0.015;
  expEvents_wzh_map[135]=0.013;
  expEvents_wzh_map[140]=0.010;
  expEvents_wzh_map[145]=0.008;
  expEvents_wzh_map[150]=0.007;
#endif


#ifdef PTL45_PTSL25_MET70 
  // chiara: qui ci metto gli eventi che passano la selezione  -----------------------------------
  std::map<int, double> expEvents_gluglu_map;
  expEvents_gluglu_map[110]=0.0484158;         
  expEvents_gluglu_map[115]=0.0554391;
  expEvents_gluglu_map[120]=0.0403891;   
  expEvents_gluglu_map[125]=0.0403891;
  expEvents_gluglu_map[130]=0.0506474;         
  expEvents_gluglu_map[135]=0.0506474;         
  expEvents_gluglu_map[140]=0.0359503;
  expEvents_gluglu_map[145]=0.0359503;
  expEvents_gluglu_map[150]=0.0359503;

  std::map<int, double> expEvents_vbf_map;
  expEvents_vbf_map[110]=0.0154531;
  expEvents_vbf_map[115]=0.0157638;
  expEvents_vbf_map[120]=0.0157943;     
  expEvents_vbf_map[125]=0.0157943;
  expEvents_vbf_map[130]=0.0115509;
  expEvents_vbf_map[135]=0.0115509;
  expEvents_vbf_map[140]=0.0119536;
  expEvents_vbf_map[145]=0.0119536;
  expEvents_vbf_map[150]=0.0119536;

  std::map<int, double> expEvents_wzh_map;
  expEvents_wzh_map[110]=0.38021;
  expEvents_wzh_map[115]=0.336158;
  expEvents_wzh_map[120]=0.330568;      
  expEvents_wzh_map[125]=0.330568;
  expEvents_wzh_map[130]=0.279578;
  expEvents_wzh_map[135]=0.279578;
  expEvents_wzh_map[140]=0.25097;
  expEvents_wzh_map[145]=0.25097;
  expEvents_wzh_map[150]=0.25097;
#endif

#ifdef PTL45_PTSL25_MET70_2012_FP 

  //provisional numbers for 2012 8TeV (not yet final xsec nor selection)
  std::map<int, double> expEvents_gluglu_map;
  expEvents_gluglu_map[110]=0.0;         
  expEvents_gluglu_map[115]=0.0;
  expEvents_gluglu_map[120]=0.0;   
  expEvents_gluglu_map[125]=0.0;
  expEvents_gluglu_map[130]=0.0;         
  expEvents_gluglu_map[135]=0.0;         
  expEvents_gluglu_map[140]=0.0;
  expEvents_gluglu_map[145]=0.0;
  expEvents_gluglu_map[150]=0.0;

  std::map<int, double> expEvents_vbf_map;
  expEvents_vbf_map[110]=0.0;
  expEvents_vbf_map[115]=0.0;
  expEvents_vbf_map[120]=0.0;     
  expEvents_vbf_map[125]=0.0;
  expEvents_vbf_map[130]=0.0;
  expEvents_vbf_map[135]=0.0;
  expEvents_vbf_map[140]=0.0;
  expEvents_vbf_map[145]=0.0;
  expEvents_vbf_map[150]=0.0;

  std::map<int, double> expEvents_wzh_map;
  expEvents_wzh_map[110]=4.08;
  expEvents_wzh_map[115]=3.65;
  expEvents_wzh_map[120]=3.61;      
  expEvents_wzh_map[125]=3.61;
  expEvents_wzh_map[130]=3.0; //provisional
  expEvents_wzh_map[135]=3.0; //provisional
  expEvents_wzh_map[140]=2.61;
  expEvents_wzh_map[145]=2.61;
  expEvents_wzh_map[150]=2.61;
#endif

#ifdef PTL45_PTSL25_MET60 

  std::map<int, double> expEvents_gluglu_map;
  expEvents_gluglu_map[110]=0.109218;         
  expEvents_gluglu_map[115]=0.109345;
  expEvents_gluglu_map[120]=0.112181;   
  expEvents_gluglu_map[125]=0.112181;
  expEvents_gluglu_map[130]=0.123648;         
  expEvents_gluglu_map[135]=0.123648;         
  expEvents_gluglu_map[140]=0.103586;
  expEvents_gluglu_map[145]=0.103586;
  expEvents_gluglu_map[150]=0.103586;

  std::map<int, double> expEvents_vbf_map;
  expEvents_vbf_map[110]=0.0320835;
  expEvents_vbf_map[115]=0.029894;
  expEvents_vbf_map[120]=0.0339233;     
  expEvents_vbf_map[125]=0.0339233;
  expEvents_vbf_map[130]=0.0263439;
  expEvents_vbf_map[135]=0.0263439;
  expEvents_vbf_map[140]=0.0292136;
  expEvents_vbf_map[145]=0.0292136;
  expEvents_vbf_map[150]=0.0292136;

  std::map<int, double> expEvents_wzh_map;
  expEvents_wzh_map[110]=0.480332;
  expEvents_wzh_map[115]=0.429865;
  expEvents_wzh_map[120]=0.416163;      
  expEvents_wzh_map[125]=0.416163;
  expEvents_wzh_map[130]=0.349557;
  expEvents_wzh_map[135]=0.349557;
  expEvents_wzh_map[140]=0.305715;
  expEvents_wzh_map[145]=0.305715;
  expEvents_wzh_map[150]=0.305715;
#endif

  typedef std::map<int,double> expEvents_t;
  std::map<TString, expEvents_t * > expEvents_map;
  expEvents_map["gluglu"]=&expEvents_gluglu_map;
  expEvents_map["vbf"]=&expEvents_vbf_map;
  expEvents_map["wzh"]=&expEvents_wzh_map;

  std::map<int, double> expEvents_tth_map;
  expEvents_tth_map[110]=0.001;
  expEvents_tth_map[115]=0.001;
  expEvents_tth_map[120]=0.001;
  expEvents_tth_map[125]=0.001;
  expEvents_tth_map[130]=0.001;
  expEvents_tth_map[135]=0.001;
  expEvents_tth_map[140]=0.000;
  expEvents_tth_map[145]=0.000;
  expEvents_tth_map[150]=0.000;


  std::map<int,double> xsec_gluglu_map;
  xsec_gluglu_map[100]=24.02;
  xsec_gluglu_map[105]=21.78;
  xsec_gluglu_map[110]=19.84;
  xsec_gluglu_map[115]=18.13;
  xsec_gluglu_map[120]=16.63;
  xsec_gluglu_map[125]=15.31;
  xsec_gluglu_map[130]=14.12;
  xsec_gluglu_map[135]=13.08;
  xsec_gluglu_map[140]=12.13;

  std::map<int,double> xsec_vbf_map;
  xsec_vbf_map[100]=1.546;
  xsec_vbf_map[105]=1.472;
  xsec_vbf_map[110]=1.398;
  xsec_vbf_map[115]=1.332;
  xsec_vbf_map[120]=1.269;
  xsec_vbf_map[125]=1.211;
  xsec_vbf_map[130]=1.154;
  xsec_vbf_map[135]=1.100;
  xsec_vbf_map[140]=1.052;

  std::map<int,double> xsec_wh_map;
  xsec_wh_map[100]=1.186;
  xsec_wh_map[105]=1.018;
  xsec_wh_map[110]=0.8754;
  xsec_wh_map[115]=0.7546;
  xsec_wh_map[120]=0.6561;
  xsec_wh_map[125]=0.5729;
  xsec_wh_map[130]=0.5008;
  xsec_wh_map[135]=0.4390;
  xsec_wh_map[140]=0.3857;

  std::map<int,double> xsec_zh_map;
  xsec_zh_map[100]=0.6313;
  xsec_zh_map[105]=0.5449;
  xsec_zh_map[110]=0.4721;
  xsec_zh_map[115]=0.4107;
  xsec_zh_map[120]=0.3598;
  xsec_zh_map[125]=0.3158;
  xsec_zh_map[130]=0.2778;
  xsec_zh_map[135]=0.2453;
  xsec_zh_map[140]=0.2172;

  std::map<int,double> xsec_wzh_map;

  if(xsec_zh_map.size() != xsec_wh_map.size()) {
    std::cerr << "[ERROR] xset maps for wh and zh have different size!" << std::endl;
    return 1;
  }

  for(std::map<int,double>::const_iterator map_itr= xsec_zh_map.begin(); 
      map_itr != xsec_zh_map.end(); 
      map_itr++){
    // ci si puo' mettere un controllo che wh abbia quel campo
    xsec_wzh_map[map_itr->first]=map_itr->second + xsec_wh_map[map_itr->first];
  }


  for (int arg=0; arg<argc; arg++) {
    if (string(argv[arg])=="-h" || string(argv[arg])=="--help") help=true;
    if (string(argv[arg])=="-v") verbose=true;
    if (string(argv[arg])=="-bkg") doBkgInt=true;
    if (string(argv[arg])=="-pG") plotGen=true;
    if (string(argv[arg])=="-sDF") saveDataFit=true;
    if (string(argv[arg])=="-gen") genName=string(argv[arg+1]);
    if (string(argv[arg])=="-fit") fitName=string(argv[arg+1]);
//     if (string(argv[arg])=="-gen") genFun_vec.push_back(argv[arg+1]);
//     if (string(argv[arg])=="-fit") fitFun_vec.push_back(argv[arg+1]);
    if (string(argv[arg])=="-t") nToys=atoi(argv[arg+1]);
    if (string(argv[arg])=="-p") toyStep=atoi(argv[arg+1]);
    if (string(argv[arg])=="-nJ") nJobs=atoi(argv[arg+1]);
    if (string(argv[arg])=="-j") jobNumb=atoi(argv[arg+1]);
    if (string(argv[arg])=="-wide") wideRange=true;
    if (string(argv[arg])=="-w") winSize=atoi(argv[arg+1]);
    if (string(argv[arg])=="-m") mL=atoi(argv[arg+1]);
    if (string(argv[arg])=="-f") fileName=string(argv[arg+1]);
    if (string(argv[arg])=="-dataf") data_fileName=string(argv[arg+1]);
    if (string(argv[arg])=="-fixed") fixed=true;
    if (string(argv[arg])=="-sigGen") sigGen=true;
    if (string(argv[arg])=="-muGen")  muGen=atof(argv[arg+1]);
    if (string(argv[arg])=="-mass") mass_vec.push_back(atof(argv[arg+1]));
  }
  toyStep=nToys/20;
  if (!help){
    //    checkFunc(genName);
    //    checkFunc(fitName);
    //    checkMass(mL);
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

  const int nGenPar=getPar(genName);
  const int nFitPar=getPar(fitName);
  string genType=getType(genName);
  string fitType=getType(fitName);
  
  cout << "--- Running with following options ---" << endl;
  cout << "    nToys:            " << nToys << endl;
  cout << "    plotStep:         " << toyStep << endl;
  cout << "    genFunction:      " << genName << endl;
  cout << "    fitFunction:      " << fitName << endl;
  cout << "    mass:             " << mL << endl;
  cout << "    winSize:          " << winSize << endl;
  cout << "    nJobs:            " << nJobs << endl;
  cout << "    jobNum:           " << jobNumb << endl;
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

  for(std::map<TString, expEvents_t * >::const_iterator channel_itr=expEvents_map.begin();
      channel_itr!= expEvents_map.end();
      channel_itr++){
    if(!channel_itr->second->count(mL)){
      std::cerr << "[ERROR] mass " << mL << " not found in " << channel_itr->first << std::endl;
      exit(2);
    }
  }

//   if( ! (
// 	 expEvents_gluglu_map.count(mL) && expEvents_vbf_map.count(mL) &&  expEvents_wzh_map.count(mL) && expEvents_tth_map.count(mL)
// 	 )) return 1;

  int toysPerJob = nToys/nJobs;
  
  if (!verbose) RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  // get workspace and mass data and gen pdfs
  TFile *inFile = new TFile(fileName.c_str());
//   RooWorkspace *dataWS = (RooWorkspace*)inFile->Get("cms_hgg_workspace");
  TFile *dataWS = new TFile(data_fileName.c_str()); 
  if(dataWS == NULL){
    std::cerr << "[WARNING] dataWS not found in WS" << std::endl;
  }
  

  TFile *fitFile = new TFile(fileName.c_str());
  RooWorkspace *fitWS;
  //  fitWS= (RooWorkspace*)fitFile->Get("events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_0.0_-1_-1_0_4_cs");
  fitWS= (RooWorkspace*)fitFile->Get("wspace");
  if(fitWS==NULL) fitWS = (RooWorkspace*)fitFile->Get("events_55.0_25.0_-10000.0_-10000.0_30.0_20.0_3.5_2.5_350.0_2.6_-1_-1_0_4_cs");
  if(fitWS == NULL){
    std::cerr << "[ERROR] fitWS not found in WS" << std::endl;
    exit(1);
  }


  RooRealVar *mass = (RooRealVar*)fitWS->var("massggnewvtx");
  if(mass == NULL){
    std::cerr << "[ERROR] variable mass not found in WS: " << fitWS->GetName() << std::endl;
    return 1;
  }

  int mLow = std::max(mL-winSize,100);
  int mHigh = std::min(mL+winSize,180);
  if(wideRange){
    mLow=100;
    mHigh=180;
  }
//   if (winSize==15) {
//     if (mL==110) {
//       mLow = 100;
//       mHigh =130;
//     }
//     else if (mL==150) {
//       mLow = 130;
//       mHigh = 160;
//     }
//   }
//   if (winSize==20){
//     if (mL==110 || mL==115) {
//       mLow = 100;
//       mHigh =140;
//     }
//     else if (mL==150) {
//       mLow = 120;
//       mHigh = 160;
//     }
//   }
  mass->setRange(mLow,mHigh);
  mass->setBins(winSize*8); // ???

  std::cout << "[INFO] Creating directory for output" << std::endl;
  system("if [ ! -e \"TestSWResults\" ]; then mkdir TestSWResults; fi");
  //system("mkdir rm -r plots");
  system("mkdir -p SWplots/toys");
  system("mkdir -p SWplots/data");
  system("mkdir -p SWplots/histos");
  
  // --- declare variables for this code -----
  //  const int nCats=4; // originale
  const int nCats=1; // shervin
  const int nMasses=1;
  const int nWinds=1;
  string winName[nWinds] = {string(Form("%d",winSize))};
  double startPar[nFitPar][nCats];
  double swStartPar[nMasses][nCats][nFitPar];
  RooDataSet *data[nCats];
  RooRealVar *fitPars[5][nCats];
  RooExponential *rooExp[6][nCats];
  RooAbsPdf *genFcn[nCats];
  RooAbsPdf *genParamatersFcn[nCats];
  RooAbsPdf *fitFcn[nCats];
  double sigSMEvents[nMasses][nCats];
  double dataSWEvents[nMasses][nCats];
  RooRealVar *mu = new RooRealVar("mu","mu",0.,-50,50);
  RooHistPdf *sigMCPdf[nMasses][nCats];
  RooRealVar *bkgYield[nMasses][nCats];
  RooFormulaVar *sigYield[nMasses][nCats];
  RooAddPdf *sigAndBkg[nMasses][nCats];
  double fwhm[nMasses][nCats];
  double maxMass[nMasses][nCats];
  
  TH1F *bias[nMasses][nCats];
  TH1F *muHist[nWinds][nMasses];


  Float_t mu_;
  Float_t muTruth_;
  Float_t sigma_mu_;
  Float_t bkgSig1fwhm_;
  Float_t bkgSig2fwhm_;
  Float_t bkgErrSig1fwhm_;
  Float_t bkgErrSig2fwhm_;
  Float_t bkgErrNormSig1fwhm_;
  Float_t bkgErrNormSig2fwhm_;
  UInt_t iToy_;
  Int_t fitStatus_;
  char genName_[30];
  char fitName_[30];
  
  sprintf(genName_, "%s", genName.c_str());
  sprintf(fitName_, "%s", fitName.c_str());
//   TString genName_=genName;
//   TString fitName_=fitName;

  std::cout << "getName = " << genName_ << std::endl;
  std::cout << "fitName = " << fitName_ << std::endl;

  TTree *muTree = new TTree("muTree","muTree");
  muTree->Branch("mass", &mL, "mass/I");
  muTree->Branch("mu", &mu_, "mu/F");
  muTree->Branch("muTruth", &muTruth_, "muTruth/F");
  muTree->Branch("sigma_mu", &sigma_mu_, "sigma_mu/F");
  muTree->Branch("bkgSig1fwhm", &bkgSig1fwhm_, "bkgSig1fwhm/F");
  muTree->Branch("bkgSig2fwhm", &bkgSig2fwhm_, "bkgSig2fwhm/F");
  muTree->Branch("bkgErrSig1fwhm", &bkgErrSig1fwhm_, "bkgErrSig1fwhm/F");
  muTree->Branch("bkgErrSig2fwhm", &bkgErrSig2fwhm_, "bkgErrSig2fwhm/F");
  muTree->Branch("bkgErrNormSig1fwhm", &bkgErrNormSig1fwhm_, "bkgErrNormSig1fwhm/F");
  muTree->Branch("bkgErrNormSig2fwhm", &bkgErrNormSig2fwhm_, "bkgErrNormSig2fwhm/F");
  muTree->Branch("genFun", genName_, "genFun/C");
  muTree->Branch("fitFun", fitName_, "fitFun/C");
  muTree->Branch("iToy", &iToy_, "iToy/i");
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
  const int lowM=mL;
  const int highM=mL;


  std::cout << "[INFO] creating new output file" << std::endl;
  TFile *outFile;

  if(!sigGen)  outFile = new TFile(Form("TestSWResults/biasCheck_g%s_f%s_m%d_w%d_j%d.root",genName.c_str(),fitName.c_str(),mL,winSize,jobNumb),"RECREATE");
  else outFile = new TFile(Form("TestSWResults/biasCheck-sigGen_%g-g%s_f%s_m%d_w%d_j%d.root",muGen, genName.c_str(),fitName.c_str(),mL,winSize,jobNumb),"RECREATE");
  
  int mIt=0;
  for (int mMC=lowM; mMC<=highM; mMC+=5){
    if (mMC==145) continue;
    for (int win=0; win<nWinds; win++) muHist[win][mIt] = new TH1F(Form("muHist_g%s_f%s_m%d_w%s",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str()),Form("muHist_g%s_f%s_m%d_w%s",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str()),100,-30,30);
    if (doBkgInt) for (int cat=0; cat<nCats; cat++) bias[mIt][cat] = new TH1F(Form("b_g%s_f%s_m%d_c%d",genName.c_str(),fitName.c_str(),mMC,cat),Form("b_g%s_f%s_m%d_c%d",genName.c_str(),fitName.c_str(),mMC,cat),100,-100,100);
    mIt++;
  }
  double toData[nMasses][nCats];


  // ------ declare fit variables ------------
  cout << "Declaring variables for fit" << endl;
  for (int cat=0; cat<nCats; cat++){
    
    // ---- get data ----
    RooDataSet *temp = (RooDataSet*)fitWS->data(Form("data_mass_cat%d",cat));
    if(temp == NULL){
      std::cout << "[WARNING] roodataset from fitWS not found" << std::endl;
    }

    mass->setRange("SW",mLow,mHigh);
    if(temp!=NULL)    data[cat] = (RooDataSet*)temp->reduce(CutRange("SW")); 
    cout << "Got data cat " << cat << endl;

    // ---- gen function ----
//     if (genName=="1pol") genFcn[cat] = (RooChebychev*)fitWS->pdf(Form("1pol_cat%d",cat));
//     if (genName=="2pol") genFcn[cat] = (RooChebychev*)fitWS->pdf(Form("2pol_cat%d",cat));
//     if (genName=="3pol") genFcn[cat] = (RooChebychev*)fitWS->pdf(Form("3pol_cat%d",cat));
//     if (genName=="4pol") genFcn[cat] = (RooChebychev*)fitWS->pdf(Form("4pol_cat%d",cat));
//     if (genName=="5pol") genFcn[cat] = (RooChebychev*)fitWS->pdf(Form("5pol_cat%d",cat));

    if (genName=="1pol") genFcn[cat] = (RooBernstein*)fitWS->pdf(Form("1pol"));
    if (genName=="2pol") genFcn[cat] = (RooBernstein*)fitWS->pdf(Form("2pol"));
    if (genName=="3pol") genFcn[cat] = (RooBernstein*)fitWS->pdf(Form("3pol"));
    if (genName=="4pol") genFcn[cat] = (RooBernstein*)fitWS->pdf(Form("4pol"));
    if (genName=="5pol") genFcn[cat] = (RooBernstein*)fitWS->pdf(Form("5pol"));

    if (genName=="1exp")    genFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp"));
    if (genName=="exp2_1") genFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp2_1"));
    if (genName=="exp2_2") genFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp2_2"));
    if (genName=="exp3_1") genFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp3_1"));
    if (genName=="exp3_2") genFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp3_2"));
    if (genName=="exp3_3") genFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp3_3"));

    //     if (genName=="1exp") genFcn[cat] = (RooExponential*)fitWS->pdf(Form("1exp_cat%d",cat));
    //     if (genName=="2exp") genFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("2exp_cat%d",cat));
    //     if (genName=="3exp") genFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("3exp_cat%d",cat));

    if (genName=="1pow") genFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("pow"));
    if (genName=="pow2_1") genFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("pow2_1"));
    if (genName=="pow2_2") genFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("pow2_2"));
    //    if (genName=="3pow") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("3pow_cat%d",cat));

//     if (genName=="1pow") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("1pow_cat%d",cat));
//     if (genName=="2pow") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("2pow_cat%d",cat));
//     if (genName=="3pow") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("3pow_cat%d",cat));
//     if (genName=="2lau") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("2lau_cat%d",cat));
//     if (genName=="4lau") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("4lau_cat%d",cat));
//     if (genName=="6lau") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("6lau_cat%d",cat));

    genParamatersFcn[cat]=0;

      //    if (genName=="1pol") genParamatersFcn[cat] = (RooBernstein*)fitWS->pdf(Form("1pol"));
    if (genName=="2pol") genParamatersFcn[cat] = (RooMultiVarGaussian*)fitWS->pdf(Form("pdf_nll2"));
    if (genName=="3pol") genParamatersFcn[cat] = (RooMultiVarGaussian*)fitWS->pdf(Form("pdf_nll"));
    if (genName=="4pol") genParamatersFcn[cat] = (RooMultiVarGaussian*)fitWS->pdf(Form("pdf_nll4"));
    
//     if (genName=="5pol") genParamatersFcn[cat] = (RooBernstein*)fitWS->pdf(Form("5pol"));

    if (genName=="1exp")    genParamatersFcn[cat] = (RooMultiVarGaussian*)fitWS->pdf(Form("pdf_nllexp1"));
//     if (genName=="exp2_1") genParamatersFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp2_1"));
//     if (genName=="exp2_2") genParamatersFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp2_2"));
//     if (genName=="exp3_1") genParamatersFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp3_1"));
//     if (genName=="exp3_2") genParamatersFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp3_2"));
//     if (genName=="exp3_3") genParamatersFcn[cat] = (RooExponential*)fitWS->pdf(Form("exp3_3"));

    //     if (genName=="1exp") genParamatersFcn[cat] = (RooExponential*)fitWS->pdf(Form("1exp_cat%d",cat));
    //     if (genName=="2exp") genParamatersFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("2exp_cat%d",cat));
    //     if (genName=="3exp") genParamatersFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("3exp_cat%d",cat));

    if (genName=="1pow") genParamatersFcn[cat] = (RooMultiVarGaussian*)fitWS->pdf(Form("pdf_nllpow1"));
//     if (genName=="pow2_1") genParamatersFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("pow2_1"));
//     if (genName=="pow2_2") genParamatersFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("pow2_2"));
    //    if (genName=="3pow") genParamatersFcn[cat] = (RooGenericPdf*)fitWS->pdf(Form("3pow_cat%d",cat));

    if(genFcn[cat]==NULL){
      std::cout << "[ERROR] genFcn not defined" << std::endl;
      return 1;
    }
    genFcn[cat]->Print();
    if (genParamatersFcn[cat]!=NULL)
      genParamatersFcn[cat]->Print();
    // ------------------------------------------
    // --- params for fit to gen data ---
    // ------------------------------------------
    // --- polynomials ---




    if (!fixed)
      {
	if (fitType=="pol"){
	  for (int par=0; par<nFitPar; par++) fitPars[par][cat] = new RooRealVar(Form("polF_p%d_cat%d",par,cat),Form("polF_p%d_cat%d",par,cat),1,0.,100.); 
	  if (fitName=="1pol") fitFcn[cat] = new RooBernstein(Form("1polF_cat%d",cat),Form("1polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat]));
	  if (fitName=="2pol") fitFcn[cat] = new RooBernstein(Form("2polF_cat%d",cat),Form("2polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat]));
	  if (fitName=="3pol") fitFcn[cat] = new RooBernstein(Form("3polF_cat%d",cat),Form("3polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
	  if (fitName=="4pol") fitFcn[cat] = new RooBernstein(Form("4polF_cat%d",cat),Form("4polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat]));
	  if (fitName=="5pol") fitFcn[cat] = new RooBernstein(Form("5polF_cat%d",cat),Form("5polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
	}
	// --- exponentials ---
	if (fitType=="exp"){
	  for (int par=0; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("expF_p%d_cat%d",par,cat),Form("expF_p%d_cat%d",par,cat),-0.1,-1.,0.);
	  for (int par=1; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("expF_p%d_cat%d",par,cat),Form("expF_p%d_cat%d",par,cat),0.1,0.,1.);
	  if (nFitPar>0) rooExp[3][cat] = new RooExponential(Form("rooexpF3_cat%d",cat),Form("rooexpF3_cat%d",cat),*mass,*fitPars[0][cat]);
	  if (nFitPar>2) rooExp[4][cat] = new RooExponential(Form("rooexpF4_cat%d",cat),Form("rooexpF4_cat%d",cat),*mass,*fitPars[2][cat]);
	  if (nFitPar>4) rooExp[5][cat] = new RooExponential(Form("rooexpF5_cat%d",cat),Form("rooexpF5_cat%d",cat),*mass,*fitPars[4][cat]);
	  if (fitName=="1exp") fitFcn[cat] = (RooExponential*)rooExp[3][cat]->clone(Form("1expF_cat%d",cat)); 
	  if (fitName=="2exp" || fitName=="exp2_1" || fitName=="exp2_2") fitFcn[cat] = new RooAddPdf(Form("2expF_cat%d",cat),Form("2expF_cat%d",cat),RooArgList(*rooExp[3][cat],*rooExp[4][cat]),RooArgList(*fitPars[1][cat]));
	  if (fitName=="3exp" || fitName=="exp3_1" || fitName=="exp3_2" || fitName=="exp3_3") fitFcn[cat] = new RooAddPdf(Form("3expF_cat%d",cat),Form("3expF_cat%d",cat),RooArgList(*rooExp[3][cat],*rooExp[4][cat],*rooExp[5][cat]),RooArgList(*fitPars[1][cat],*fitPars[3][cat]));
	}
	// --- power laws ---
	if (fitType=="pow"){
	  for (int par=0; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("powF_p%d_cat%d",par,cat),Form("powF_p%d_cat%d",par,cat),-1.,-50.,0.);
	  for (int par=1; par<nFitPar; par+=2) fitPars[par][cat] = new RooRealVar(Form("powF_p%d_cat%d",par,cat),Form("powF_p%d_cat%d",par,cat),0.8,0.,1.);
	  if (fitName=="1pow") fitFcn[cat] = new RooGenericPdf(Form("1powF_cat%d",cat),Form("1powF_cat%d",cat),"pow(@0,@1)",RooArgList(*mass,*fitPars[0][cat]));
	  if (fitName=="2pow") fitFcn[cat] = new RooGenericPdf(Form("2powF_cat%d",cat),Form("2powF_cat%d",cat),"(1-@2)*pow(@0,@1)+@2*pow(@0,@3)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
	  if (fitName=="pow2_1" || fitName=="pow2_2") fitFcn[cat] = new RooGenericPdf(Form("2powF_cat%d",cat),Form("2powF_cat%d",cat),"(1-@2)*pow(@0,@1)+@2*pow(@0,@3)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
	  if (fitName=="3pow") fitFcn[cat] = new RooGenericPdf(Form("3powF_cat%d",cat),Form("3powF_cat%d",cat),"(1-@2-@4)*pow(@0,@1)+@2*pow(@0,@3)+@4*pow(@0,@5)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
	}
	// --- laurent series ---
	if (fitType=="lau"){
	  for (int par=0; par<nFitPar; par++) fitPars[par][cat] = new RooRealVar(Form("lauF_p%d_cat%d",par,cat),Form("lauF_p%d_cat%d",par,cat),0.5,0.,1.);
	  if (fitName=="2lau") fitFcn[cat] = new RooGenericPdf(Form("2lauF_cat%d",cat),Form("2lauF_cat%d",cat),"(1-@1)*pow(@0,-4.0)+@1*pow(@0,-5.0)",RooArgList(*mass,*fitPars[0][cat]));
	  if (fitName=="4lau") fitFcn[cat] = new RooGenericPdf(Form("4lauF_cat%d",cat),Form("4lauF_cat%d",cat),"(1-@1-@2-@3)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
	  if (fitName=="6lau") fitFcn[cat] = new RooGenericPdf(Form("6lauF_cat%d",cat),Form("6lauF_cat%d",cat),"(1-@1-@2-@3-@4-@5)*pow(@0,-4.0)+@1*pow(@0,-5.0)+@2*pow(@0,-3.0)+@3*pow(@0,-6.0)+@4*pow(@0,-2.0)+@5*pow(@0,-7.0)",RooArgList(*mass,*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
	}
      }
    else
      {
	if (fitName=="1pol") fitFcn[cat] = (RooAbsPdf *) ((RooBernstein*)fitWS->pdf(Form("1pol")))->Clone();
	if (fitName=="2pol") fitFcn[cat] = (RooAbsPdf *) ((RooBernstein*)fitWS->pdf(Form("2pol")))->Clone();
	if (fitName=="3pol") fitFcn[cat] = (RooAbsPdf *) ((RooBernstein*)fitWS->pdf(Form("3pol")))->Clone();
	//     if (fitName=="4pol") fitFcn[cat] = (RooAbsPdf *) ((RooBernstein*)fitWS->pdf(Form("4pol")))->Clone();
	//     if (fitName=="5pol") fitFcn[cat] = (RooAbsPdf *) ((RooBernstein*)fitWS->pdf(Form("5pol")))->Clone();
	
	if (fitName=="1exp")   fitFcn[cat] = (RooAbsPdf *) ((RooExponential*)fitWS->pdf(Form("exp")))->Clone();
	if (fitName=="exp2_1") fitFcn[cat] = (RooAbsPdf *) ((RooExponential*)fitWS->pdf(Form("exp2_1")))->Clone();
	if (fitName=="exp2_2") fitFcn[cat] = (RooAbsPdf *) ((RooExponential*)fitWS->pdf(Form("exp2_2")))->Clone();
	if (fitName=="exp3_1") fitFcn[cat] = (RooAbsPdf *) ((RooExponential*)fitWS->pdf(Form("exp3_1")))->Clone();
	if (fitName=="exp3_2") fitFcn[cat] = (RooAbsPdf *) ((RooExponential*)fitWS->pdf(Form("exp3_2")))->Clone();
	if (fitName=="exp3_3") fitFcn[cat] = (RooAbsPdf *) ((RooExponential*)fitWS->pdf(Form("exp3_3")))->Clone();
	
// 	if(TString(fitName).Contains("pol")){
// 	  for (int par=0; par<nFitPar; par++) fitPars[par][cat] = new RooRealVar(Form("polF_p%d_cat%d",par,cat),Form("polF_p%d_cat%d",par,cat),1,0.,100.); 
// 	  //       if (fitName=="1pol") fitFcn[cat] = new RooBernstein(Form("1polF_cat%d",cat),Form("1polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat]));
// 	  //       if (fitName=="2pol") fitFcn[cat] = new RooBernstein(Form("2polF_cat%d",cat),Form("2polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat]));
// 	  //       if (fitName=="3pol") fitFcn[cat] = new RooBernstein(Form("3polF_cat%d",cat),Form("3polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat]));
// 	  if (fitName=="4pol") fitFcn[cat] = new RooBernstein(Form("4polF_cat%d",cat),Form("4polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat]));
// 	  if (fitName=="5pol") fitFcn[cat] = new RooBernstein(Form("5polF_cat%d",cat),Form("5polF_cat%d",cat),*mass,RooArgList(*fitPars[0][cat],*fitPars[1][cat],*fitPars[2][cat],*fitPars[3][cat],*fitPars[4][cat]));
// 	}
	
	//     if (fitName=="1exp") fitFcn[cat] = (RooExponential*)fitWS->pdf(Form("1exp_cat%d",cat));
	//     if (fitName=="2exp") fitFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("2exp_cat%d",cat));
	//     if (fitName=="3exp") fitFcn[cat] = (RooAddPdf*)fitWS->pdf(Form("3exp_cat%d",cat));
	
	if (fitName=="1pow") fitFcn[cat] =   (RooAbsPdf *) ((RooGenericPdf*)fitWS->pdf(Form("pow")))->Clone();
	if (fitName=="pow2_1") fitFcn[cat] = (RooAbsPdf *) ((RooAddPdf*)fitWS->pdf(Form("pow2_1")))->Clone();
	if (fitName=="pow2_2") fitFcn[cat] = (RooAbsPdf *) ((RooAddPdf*)fitWS->pdf(Form("pow2_2")))->Clone();
	RooArgSet* par=fitFcn[cat]->getParameters(*mass);

	//Setting all parameters fixed
	TIterator* iter = par->createIterator();
	for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
          RooRealVar *v = dynamic_cast<RooRealVar *>(a);
	  std::cout << "Setting " << v->GetName() << " fixed " << std::endl;
	  v->setConstant(true);
	}
	delete iter;
      }




  // ---- Get mass distributions  and entries ------
    cout << "--- Event check: ----- " << endl;
    mIt=0;
    // commentare mIt
    // for(unsigned int mIt=0; mIt < mass_vec.size(); mIt++){
    for (int mMC=lowM; mMC<=highM; mMC+=5){
      //      if (mMC==145) continue;
    // get signal data
      TH1F *gluglu_hist = (TH1F *) dataWS->Get(Form("var_gluglu_2011_%d",mMC));
      TH1F *vbf_hist = (TH1F *) dataWS->Get(Form("var_vbf_2011_%d",mMC));
      TH1F *wzh_hist = (TH1F *) dataWS->Get(Form("var_wzh_2011_%d",mMC));
      
      if(!(gluglu_hist || vbf_hist || wzh_hist)){
	std::cerr << "[ERROR] histogram " << " for m= " << mMC << " not found in " << dataWS->GetName()  << std::endl;
	return 2;
      }
  

      TH1F *sigDataHist = (TH1F *) gluglu_hist->Clone("sigDataHist");
      sigDataHist->Reset();
      sigDataHist->Add(gluglu_hist, xsec_gluglu_map[mMC]);
      sigDataHist->Add(vbf_hist, xsec_vbf_map[mMC]);
      sigDataHist->Add(wzh_hist, xsec_wzh_map[mMC]);

      //We can use the fit in the future...
      int bin1=sigDataHist->FindFirstBinAbove(sigDataHist->GetMaximum()/2.);
      int bin2=sigDataHist->FindLastBinAbove(sigDataHist->GetMaximum()/2.);
      fwhm[mIt][cat]= sigDataHist->GetBinCenter(bin2)- sigDataHist->GetBinCenter(bin1);
      maxMass[mIt][cat]=sigDataHist->GetBinCenter(sigDataHist->GetMaximumBin());

//       RooDataSet *sigData = (RooDataSet*)dataWS->data(Form("var_gluglu_2011_%d",mMC));
//       sigData->append(*((RooDataSet*)dataWS->data(Form("var_vbf_2011_%d",mMC))));
//       sigData->append(*((RooDataSet*)dataWS->data(Form("var_wzh_2011_%d",mMC))));
      //      sigData->append(*((RooDataSet*)dataWS->data(Form("var_tth_mass_m%d_cat%d",mMC,cat))));

    // get expected SM events

      sigSMEvents[mIt][cat] = expEvents_gluglu_map[mMC] + expEvents_vbf_map[mMC] + expEvents_wzh_map[mMC]; 
      std::cerr << "Signal events = " << sigSMEvents[mIt][cat] << "\t" << mIt << "\t" << cat << std::endl;

	// (((TH1F*)inFile->Get(Form("th1f_sig_ggh_mass_m%d_cat%d",mMC,cat)))->Integral())+((TH1F*)inFile->Get(Form("th1f_sig_vbf_mass_m%d_cat%d",mMC,cat)))->Integral()+((TH1F*)inFile->Get(Form("th1f_sig_wzh_mass_m%d_cat%d",mMC,cat)))->Integral()+((TH1F*)inFile->Get(Form("th1f_sig_tth_mass_m%d_cat%d",mMC,cat)))->Integral();
      cout << "   Mass: " << mMC << " cat: " << cat << " intSM: " << sigSMEvents[mIt][cat] << " maxPDF " << maxMass[mIt][cat] << " fwhmPDF " << fwhm[mIt][cat] << endl;

    // construct s+b model
      RooDataHist *sigMC = new RooDataHist(Form("sigMC_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,sigDataHist);
      //      RooDataHist *sigMC = new RooDataHist(Form("sigMC_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,sigData);
      sigMCPdf[mIt][cat] = new RooHistPdf(Form("sigMCPdf_m%d_cat%d",mMC,cat),Form("sigMC_m%d_cat%d",mMC,cat),*mass,*sigMC);

   // ----- def s and b yields and construct s+b model
      bkgYield[mIt][cat] = new RooRealVar(Form("bkgYield_m%d_cat%d",mMC,cat),Form("bkgYield_m%d_cat%d",mMC,cat),1e3,0,1e10);
      sigYield[mIt][cat] = new RooFormulaVar(Form("sigYield_m%d_cat%d",mMC,cat),Form("sigYield_m%d_cat%d",mMC,cat),"@0*@1",RooArgList(RooConst(sigSMEvents[mIt][cat]),*mu));
      sigAndBkg[mIt][cat] = new RooAddPdf(Form("SandB_m%d_cat%d",mMC,cat),Form("SandB_m%d_cat%d",mMC,cat),RooArgList(*fitFcn[cat],*sigMCPdf[mIt][cat]),RooArgList(*bkgYield[mIt][cat],*sigYield[mIt][cat]));

      TCanvas *tC = new TCanvas();
      RooPlot *tPlot = mass->frame();
      sigMCPdf[mIt][0]->plotOn(tPlot,LineColor(kMagenta));
      sigYield[mIt][0]->Print();
      tPlot->Draw();
      tC->Print(Form("SWplots/test_gen%s_cat%d_toy%d.png",genName.c_str(),cat,0),"png");

      mIt++;
    }
  }

  cout << "------------------------------------------------" << endl;
  cout << "--- Data fitted. Mass distributions obtained ---" << endl;
  cout << "------------------------------------------------" << endl;
  cout << "-------- Generating and fitting toys -----------" << endl;

  RooDataSet *genDatUnBin[nCats];

  RooDataHist *genDat[nCats];
  RooRandom::randomGenerator()->SetSeed(0);
  TRandom3 bkgGenRnd;
      RooRealVar sigGenEvents("sigGenEvents","sigGenEvents", 0,1e10);
      RooRealVar bkgGenEvents("bkgGenEvents","bkgGenEvents", 0,1e10);

  for (int itToy=jobNumb*toysPerJob; itToy<(jobNumb+1)*toysPerJob; itToy++){
    cout << "----------------------------------------------" << endl;
    cout << "------------- TOY: " << itToy << "------------" << endl;
    cout << "----------------------------------------------" << endl;
    for (int cat=0; cat<nCats; cat++){

      //      if (genParamatersFcn[cat])
      //	genParamatersFcn[cat]->generate(
#ifdef VBF
      int nBkgExp=120;
#endif
#ifdef PTL45_PTSL25_MET60
      int nBkgExp=61;
#endif
#ifdef PTL45_PTSL25_MET70
      int nBkgExp=26;
#endif
#ifdef PTL45_PTSL25_MET70_2012_FP
      int nBkgExp=127;
#endif
      
#ifdef NOPOIS
      genDatUnBin[cat] = genFcn[cat]->generate(*mass,nBkgExp);
#else
      genDatUnBin[cat] = genFcn[cat]->generate(*mass,nBkgExp,Extended());
#endif
      //      if(data[cat] != NULL) data[cat]->sumEntries();
      //      genDatUnBin[cat] = genFcn[cat]->generate(*mass,bkgGenEvents.getVal(),Extended());

      genDatUnBin[cat]->Print();
      genDat[cat] = new RooDataHist("gen","gen",*mass,*genDatUnBin[cat]);
      
      //      genDat[cat]->Print();
      /*
      TCanvas *tC = new TCanvas();
      RooPlot *tPlot = mass->frame();
      genDat[cat]->plotOn(tPlot,DataError(RooDataSet::SumW2));
      genFcn[cat]->plotOn(tPlot,LineColor(kMagenta));
      tPlot->Draw();
      tC->Print(Form("SWplots/temp/gen%s_cat%d_toy%d.png",genName.c_str(),cat,itToy),"png");
      mu->setVal(0.);
      */
    }
    // --- set up categories and combine data ---- 
    RooCategory category("category","category");
    category.defineType("cat0");
    if(nCats>1){
      category.defineType("cat1");
      if(nCats>2){
	category.defineType("cat2");
	if(nCats>3){
	  category.defineType("cat3");
	}
      }
    }

    //   RooDataHist combData(Form("combData_toy%d",itToy),Form("combData_toy%d",itToy),*mass,Index(category),Import("cat0",*genDat[0]));
       
    //RooDataHist combData(Form("combData_toy%d",itToy),Form("combData_toy%d",itToy),*mass,Index(category),Import("cat0",*genDat[0]),Import("cat1",*genDat[1]),Import("cat2",*genDat[2]),Import("cat3",*genDat[3]));
    

    // ---- construct RooSimultaneous from 4 cats -----
    mIt=0;
    for (int mMC=lowM; mMC<=highM; mMC+=5){
      //      if (mMC==145) continue;
      mass->setRange(mLow,mHigh);
      //mass->setBins(80);
      RooDataSet combData(Form("combData_toy%d",itToy),Form("combData_toy%d",itToy),*mass,Index(category),Import("cat0",*genDatUnBin[0]));
      
      if(sigGen){
	//	combData.Print();
	sigGenEvents.setVal(bkgGenRnd.Poisson(muGen*sigSMEvents[mIt][0]));
	muTruth_=(sigGenEvents.getVal())/sigSMEvents[mIt][0];	
	std::cout << "Generating " << sigGenEvents.getVal() << " signal events starting from " << muGen*sigSMEvents[mIt][0] << std::endl;
	RooDataSet *genSim=sigMCPdf[0][0]->generate(*mass, sigGenEvents.getVal());

	combData.append(*genSim);
	//	combData.Print();
      }else muTruth_=0;
      
      for (int win=0; win<nWinds; win++){
	
        RooSimultaneous simPdf(Form("simPdf%d_toy%d",mMC,itToy),Form("simPdf%d_toy%d",mMC,itToy),category);
        RooSimultaneous simBkgOnlyPdf(Form("simBkgOnlyPdf%d_toy%d",mMC,itToy),Form("simBkgOnlyPdf%d_toy%d",mMC,itToy),category);
        RooAddPdf *SandB[nCats];
        RooAddPdf *B[nCats];
        for (int cat=0; cat<nCats; cat++) {
	  //	  sigAndBkg[mIt][cat]->Print();
	  SandB[cat] = (RooAddPdf*)sigAndBkg[mIt][cat]->Clone(Form("SandB_m%d_w%s_c%d",mMC,winName[win].c_str(),cat));
	  B[cat] = (RooAddPdf*) fitFcn[cat]->Clone(Form("SandB_m%d_w%s_c%d",mMC,winName[win].c_str(),cat));
	  //	  SandB[cat]->Print();
	}
	// da mettere nel ciclo for!
        simPdf.addPdf(*SandB[0],"cat0");
//         simPdf.addPdf(*SandB[1],"cat1");
//         simPdf.addPdf(*SandB[2],"cat2");
//         simPdf.addPdf(*SandB[3],"cat3");

	simBkgOnlyPdf.addPdf(*B[0],"cat0");
	
        // ---- reset starting vals ----
        mu->setVal(0.);

        for (int cat=0; cat<nCats; cat++){
          bkgYield[mIt][cat]->setVal(1000.);
          //for (int p=0; p<nFitPar; p++) fitPars[p][cat]->setVal(swStartPar[mIt][cat][p]);
        }
        // ---- fit and save output -----
	//	combData.Print();
        RooFitResult *res;
	if (verbose) res = simPdf.fitTo(combData,Save(true));
	//        if (verbose) res = simPdf.fitTo(*(genDatUnBin[0]),Save(true));
        else res = simPdf.fitTo(combData,Save(true),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
	res->Print();

        cout << "------------------------------------------------" << endl;
        cout << "--- TOY: " << itToy << " Gen: " << genName << " Fit: " << fitName << " mass: " << mMC << " window: " << winName[win] << endl;
        cout << "------------------------------------------------" << endl;
        res->floatParsFinal().Print("s");

	//Doing a BKG ONLY FIT
//         RooFitResult *res2;
// 	if (verbose) res2 = simBkgOnlyPdf.fitTo(combData,Save(true));
// 	//        if (verbose) res = simPdf.fitTo(*(genDatUnBin[0]),Save(true));
//         else res2 = simBkgOnlyPdf.fitTo(combData,Save(true),PrintLevel(-1),Warnings(false),PrintEvalErrors(-1));
// 	res2->Print();

//         res2->floatParsFinal().Print("s");

	mass->setRange("sigWindow",maxMass[mIt][0]-fwhm[mIt][0]/2.,maxMass[mIt][0]+fwhm[mIt][0]/2.);
	RooAbsReal *intRange = fitFcn[0]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
	RooFormulaVar *normIntRange= new RooFormulaVar(Form("normIntRange_m%d_cat%d",mIt,0),Form("normIntRange_m%d_cat%d",mIt,0),"@0*@1",RooArgList(*intRange,*bkgYield[mIt][0]));
	mass->setRange("sigWindow2",maxMass[mIt][0]-fwhm[mIt][0],maxMass[mIt][0]+fwhm[mIt][0]);
	RooAbsReal *int2Range = fitFcn[0]->createIntegral(*mass,NormSet(*mass),Range("sigWindow2"));
	RooFormulaVar *normInt2Range= new RooFormulaVar(Form("normInt2Range_m%d_cat%d",mIt,0),Form("normInt2Range_m%d_cat%d",mIt,0),"@0*@1",RooArgList(*int2Range,*bkgYield[mIt][0]));

	std::cout << normIntRange->getVal() << " , "  << normIntRange->getPropagatedError(*res) << std::endl;
	std::cout << normInt2Range->getVal() << " , "  << normInt2Range->getPropagatedError(*res) << std::endl;

	mu_=((RooRealVar *)(res->floatParsFinal().find("mu")))->getVal();
	sigma_mu_= ((RooRealVar *)(res->floatParsFinal().find("mu")))->getError();

	bkgSig1fwhm_ = normIntRange->getVal();
	bkgSig2fwhm_ = normInt2Range->getVal();
	bkgErrSig1fwhm_ = normIntRange->getPropagatedError(*res);
	bkgErrSig2fwhm_ = normInt2Range->getPropagatedError(*res);
	bkgErrNormSig1fwhm_ = normIntRange->getPropagatedError(*res)/sigSMEvents[mIt][0] ;
	bkgErrNormSig2fwhm_ = normInt2Range->getPropagatedError(*res)/sigSMEvents[mIt][0] ;	
	iToy_=itToy;

	//	if(res->status()==0)	
	fitStatus_=res->status();

	muTree->Fill();

	delete normInt2Range;
	delete normIntRange;

        outFile->cd();
        res->SetName(Form("fitRes_g%s_f%s_m%d_w%s_t%d",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),itToy));
        //res->Write();
	//SandB[0]->Print();
	//	fitFcn[0]->Print();

        if (itToy>0) muHist[win][mIt]->Fill(mu->getVal()); 
        // --- find bkg integral -----
        for (int cat=0; cat<nCats; cat++){
          if (doBkgInt){
            RooAbsReal *intInRange = fitFcn[cat]->createIntegral(*mass,NormSet(*mass),Range("sigWindow"));
            double toBkg = intInRange->getVal()*genDat[cat]->numEntries();
            bias[mIt][cat]->Fill(toData[mIt][cat]-toBkg);
          }
        // ---- make plots ----
          if (itToy%toyStep==0 && itToy>0){
            TCanvas *c1 = new TCanvas();
            RooPlot *mFrame = mass->frame(Title(Form("Gen: %s. Fit: %s. Mass %d win #pm%s cat %d toy %d",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),cat,itToy)));
            //if (wideRange) mass->setBins(160,"coarse");
            //else mass->setBins(120,"coarse");

	          //mass->setBins(80);

	    combData.plotOn(mFrame,DataError(RooDataSet::SumW2));//,Binning("coarse"));
            //genDat[cat]->plotOn(mFrame,DataError(RooDataSet::SumW2));//,Binning("coarse"));
            if (plotGen) {
	      RooAddPdf genPdf(Form("Sig_Bkg_gen_m%d_cat%d",mMC,cat),Form("Sig_Bkg_gen_m%d_cat%d",mMC,cat),RooArgList(*genFcn[cat],*sigMCPdf[0][cat]),RooArgList(bkgGenEvents,sigGenEvents));
	      if(sigGen)   genPdf.plotOn(mFrame,LineColor(kMagenta)); //genFcn[cat]->plotOn(mFrame,LineColor(kMagenta)); 
	      else	      genFcn[cat]->plotOn(mFrame,LineColor(kMagenta));
            //fitFcn[cat]->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed));
	    }
	    //	    SandB[cat]->plotOn(mFrame,VisualizeError(*res,1,kFALSE),DrawOption("L"),LineWidth(2),LineColor(kOrange),Components(*(SandB[cat]->pdfList().at(0))),LineStyle(kDashed));
            SandB[cat]->plotOn(mFrame,LineColor(kRed),LineStyle(kDashed),Components(*(SandB[cat]->pdfList().at(0))));//,Range(100,160));
            SandB[cat]->plotOn(mFrame);


            mFrame->SetXTitle("m_{#gamma#gamma} (GeV/c^{2})");
            mFrame->Draw();
        // ---- make legend -----
            TPaveText *text = new TPaveText(0.55,0.7,0.65,0.89,"NDC");
            text->SetLineColor(0);
            text->SetFillColor(0);
            text->AddText(Form("#mu = %1.2f", mu->getVal()));
            leg->Draw("same");
            text->Draw("same");
	    if(sigGen)    c1->Print(Form("SWplots/toys/fitTogen-sigGen_%g-g%s_f%s_m%d_w%s_cat%d_toy%d.png",muGen, genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),cat,itToy),"png");
            else c1->Print(Form("SWplots/toys/fitTogen_g%s_f%s_m%d_w%s_cat%d_toy%d.png",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),cat,itToy),"png");
            
            //delete c1;
            //delete h;
            //delete h1;
            //delete h2;
            //delete text;
            //delete leg;
          }
        // --- set starting to vals to that of gen
          //for (int par=0; par<nFitPar; par++) fitPars[par][cat]->setVal(startPar[par][cat]);
          //for (int mIt=0; mIt<nMasses; mIt++) bkgYield[mIt][cat]->setVal(5000);
        }




	//mu->setVal(0.);
      }
      mIt++;
    }
  }

  outFile->cd();
  muTree->Write();

//   mIt=0;

//   for (int mMC=lowM; mMC<=highM; mMC+=5){
//     if (mMC==145) continue;
//     for (int win=0; win<nWinds; win++){
//       muHist[win][mIt]->Write();
//       cout << "m" << mMC << " w" << win << endl;
//       TF1 *gausFit = new TF1("gausFit","gaus",-10,10);
//       muHist[win][mIt]->Fit(gausFit);
//       TCanvas *canv = new TCanvas();
//       muHist[win][mIt]->Draw();
//       canv->Print(Form("SWplots/histos/mu_g%s_f%s_m%d_w%s_j%d.png",genName.c_str(),fitName.c_str(),mMC,winName[win].c_str(),jobNumb),"png");
//       if (doBkgInt){
//         for (int cat=0; cat<nCats; cat++){
//           bias[mIt][cat]->Write();
//           bias[mIt][cat]->Draw();
//           canv->Print(Form("SWplots/histos/biasCheck_g%s_f%s_m%d_cat%d.png",genName.c_str(),fitName.c_str(),mMC,cat),"png");
//           canv->Clear();
//         }
//       }
//       //delete canv;
//     }
//     mIt++;
//   }

  ofstream complete(Form("SWResults/g%s_f%s_j%d.txt",genName.c_str(),fitName.c_str(),jobNumb));
  complete << "Job: " << jobNumb << "/" << nJobs << " gen " << genName << " fit " << fitName << " completed successfully" << endl;
  complete.close();
  cout << "Text file written " << endl;
  outFile->Close();
  inFile->Close();
  fitFile->Close();

  return 0;
}
