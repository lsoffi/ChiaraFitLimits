//////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #205
// 
// Options for plotting components of composite p.d.f.s.
//
//
//
// 07/2008 - Wouter Verkerke 
// 
/////////////////////////////////////////////////////////////////////////

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
using namespace RooFit ;


void rf205_compplot()
{
  // S e t u p   c o m p o s i t e    p d f
  // --------------------------------------

  // Declare observable x
  RooRealVar x("x","x",100,800) ;

 
  // Build Chebychev polynomial p.d.f.  
  RooRealVar a0("a0","a0",0.5,0.,1.) ;
  RooRealVar a1("a1","a1",0.2,0.,1.) ;
  RooChebychev bkg1("bkg1","Background 1",x,RooArgSet(a0,a1)) ;

  // Build expontential pdf
  RooRealVar alpha("alpha","alpha",-0.0001, -100., 100.) ;
  RooExponential bkg2("bkg2","Background 2",x,alpha) ;
 // Build expontential pdf
  RooRealVar alpha2("alpha2","alph2",-0.002, -100., 100.) ;
  RooExponential bkg3("bkg3","Background 3",x,alpha2) ;

  
  RooRealVar p1("p1","p1",6.52,1., 20);
  RooRealVar p2("p2","p2",0.0000001, 0., 0.005);
  RooRealVar p3("p3","p3",0.36, 0., 10);
  
  RooGenericPdf DiJet("DiJet", "pow(1-@0/8000., @2)/pow(@0/8000., @1+@3*log(@0/8000.))", RooArgList(x, p1, p2,p3));

  
  // Sum the background components into a composite background p.d.f.
  RooRealVar bkg1frac("bkg1frac","fraction of component 1 in background",0.2,0.,1.) ;
  //  RooAddPdf bkg("bkg","",RooArgList(bkg2,bkg3),bkg1frac) ;
  RooAddPdf bkgD("bkgD","",RooArgList(DiJet,bkg3),bkg1frac) ;
  
  


  // S e t u p   b a s i c   p l o t   w i t h   d a t a   a n d   f u l l   p d f 
  // ------------------------------------------------------------------------------
  
  // Generate a data sample of 1000 events in x from model
  RooDataSet *data = bkgD.generate(x,15000) ;
  bkgD.fitTo(*data,RooFit::FitOptions("MHTER"), Save(kTRUE));// SumW2Error(kTRUE),
  // Plot data and complete PDF overlaid
  RooPlot* xframe  = x.frame(Title("Component plotting of pdf")) ;
  data->plotOn(xframe) ;
  bkgD.plotOn(xframe, LineColor(kBlue)) ;
  bkgD.plotOn(xframe, Components(DiJet), LineColor(kGreen), LineStyle(kDashed)) ;
  bkgD.plotOn(xframe, Components(bkg3), LineColor(kRed), LineStyle(kDashed)) ;

  

  // Draw the frame on the canvas
  TCanvas* c = new TCanvas("rf205_compplot","rf205_compplot",800,400) ;
 
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; xframe->GetYaxis()->SetTitleOffset(1.4) ; xframe->Draw() ;
  


}
