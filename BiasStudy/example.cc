#include "RooFit.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooGlobalFunc.h"
#include "RooFit.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooExpolPdf.h"
#include "RooPolynomial.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
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
#include "TROOT.h"
#include "TSystem.h"
#include "RooCBShape.h"
#include "RooMinuit.h"
#include "RooNLLVar.h"
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"

using namespace RooFit;
void test(){

  float mass_ = 250.;
  float minMassFit = 100;
  float maxMassFit = 800; 
  int c = 0;
  RooRealVar* mass = new RooRealVar("mass", "mass",minMassFit, maxMassFit);
 


  //add con dijet

  //par DiJetEXP
  RooRealVar p1DiJetEXP("p1DiJetEXP","p1DiJetEXP",6.52,1., 20);
  RooRealVar p2DiJetEXP("p2DiJetEXP","p2DiJetEXP",0.000001, 0., 100.);
  RooRealVar p3DiJetEXP("p3DiJetEXP","p3DiJetEXP",0.36, 0., 10);

  RooRealVar p4DiJetEXP("p4DiJetEXP","p4DiJetEXP", 0.2,0., 1.);//frac
 
  RooRealVar p5DiJetEXP("p5DiJetEXP","p5DiJetEXP",-0.004, -10000., 0.);
  RooRealVar p6DiJetEXP("p6DiJetEXP","p6DiJetEXP",-0.002, -10000., 0.);
  

  RooFormulaVar *x     = new RooFormulaVar(TString::Format("xDiJetEXP_cat%d",c),"","@0/8000.",*mass);

  //dijet e exp
  RooGenericPdf* PhotonsMassBkgTmp0DiJet = new RooGenericPdf("dijet", "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(*mass, p1DiJetEXP, p2DiJetEXP,p3DiJetEXP)); 

  RooExponential* PhotonsMassBkgTmp0Exp1 = new RooExponential("exp1","", *mass,  p5DiJetEXP);
  RooExponential* PhotonsMassBkgTmp0Exp2 = new RooExponential("exp2","", *mass,  p6DiJetEXP);
  
  //rooadd  
  RooAddPdf* PhotonsMassBkgTmp0 = new RooAddPdf("add", "add" , RooArgList(*PhotonsMassBkgTmp0Exp1, *PhotonsMassBkgTmp0DiJet), RooArgList(p4DiJetEXP));
   
 
  //generate dataset
  RooDataSet* data = (RooDataSet*)PhotonsMassBkgTmp0->generate(*mass, 15000);  
  
  //draw everything
  RooPlot* plot = mass->frame();
  data->plotOn(plot);

  PhotonsMassBkgTmp0->plotOn(plot,RooFit::LineColor(kBlue));
  PhotonsMassBkgTmp0->plotOn(plot,Components("exp1"),RooFit::LineColor(kRed), LineStyle(kDashed));
  PhotonsMassBkgTmp0->plotOn(plot,Components("dijet"),RooFit::LineColor(kGreen), LineStyle(kDashed));
  
  TCanvas* c_test = new TCanvas("ctest", "ctest",1);
  c_test->cd();
  c_test->SetLogy();
  plot->GetYaxis()->SetRangeUser(1., 10000);
  plot->Draw();
  c_test->SaveAs("testOnlyBkgBeforeFit.png");
   
  //fit
  PhotonsMassBkgTmp0->fitTo(*data,RooFit::FitOptions("M"), Save(kTRUE));
  
  RooPlot* plot2 = mass->frame();
  data->plotOn(plot2);

  PhotonsMassBkgTmp0->plotOn(plot2,RooFit::LineColor(kBlue));
  PhotonsMassBkgTmp0->plotOn(plot2,Components("exp1"),RooFit::LineColor(kRed), LineStyle(kDashed));
  PhotonsMassBkgTmp0->plotOn(plot2,Components("dijet"),RooFit::LineColor(kGreen), LineStyle(kDashed));
  
  plot2->Draw();
  c_test->SaveAs("testOnlyBkgAfterFit.png");
  

}


