/////////////////////////////////////////////////////////////////////////
//
// 'ADDITION AND CONVOLUTION' RooFit tutorial macro #208
// 
// One-dimensional numeric convolution
// (require ROOT to be compiled with --enable-fftw3)
// 
// pdf = landau(t) (x) gauss(t)
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
#include "RooBreitWigner.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TH1.h"
using namespace RooFit ;



void rf208_convolution()
{
  // S e t u p   c o m p o n e n t   p d f s 
  // ---------------------------------------

  // Construct observable
  RooRealVar t("t","t",-10,30) ;

  // Construct landau(t,ml,sl) ;
  RooRealVar ml("ml","mean bw",5.,-20,20) ;
  RooRealVar sl("sl","sigma bw",1,0.1,10) ;
  RooBreitWigner bw("bw","bw",t,ml,sl) ;
  
  // Construct gauss(t,mg,sg)
  RooRealVar mg("mg","mg",0) ;
  RooRealVar sg("sg","sg",2,0.1,10) ;
  RooGaussian gauss("gauss","gauss",t,mg,sg) ;


  // C o n s t r u c t   c o n v o l u t i o n   p d f 
  // ---------------------------------------

  // Set #bins to be used for FFT sampling to 10000
  t.setBins(10000,"cache") ; 

  // Construct landau (x) gauss
  RooFFTConvPdf lxg("lxg","bw (X) gauss",t,bw,gauss) ;



  // S a m p l e ,   f i t   a n d   p l o t   c o n v o l u t e d   p d f 
  // ----------------------------------------------------------------------

  // Sample 1000 events in x from gxlx
  RooDataSet* data = lxg.generate(t,10000) ;

  // Fit gxlx to data
  lxg.fitTo(*data) ;

  // Plot data, landau pdf, landau (X) gauss pdf
  RooPlot* frame = t.frame(Title("landau (x) gauss convolution")) ;
  data->plotOn(frame) ;
  lxg.plotOn(frame) ;
  bw.plotOn(frame,LineStyle(kDashed)) ;


  // Draw frame on canvas
  new TCanvas("rf208_convolution","rf208_convolution",600,600) ;
  gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.4) ; frame->Draw() ;

  //add a variable to the dataset
  RooFormulaVar *r_formula     = new RooFormulaVar("r_formula","","@0",t);
  RooRealVar* r = (RooRealVar*) data->addColumn(*r_formula);
  r->SetName("r");
  r->SetTitle("r");

  RooDataSet* data_r =(RooDataSet*) data->reduce(*r, "");
  r->setRange("sigrange",-10.,30.);
  RooPlot* r_frame = r->frame(Range("sigRange"),Title(" r (x) gauss convolution")) ;
  data_r->plotOn(r_frame, MarkerColor(kRed));
  r_frame->GetXaxis()->SetRangeUser(-10., 30.);
  r_frame->Draw() ;
}



