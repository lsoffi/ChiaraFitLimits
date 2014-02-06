/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/include/Riostream.h" 

#include "RooCBCrujffPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooCBCrujffPdf) 

 RooCBCrujffPdf::RooCBCrujffPdf(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _m0,
                        RooAbsReal& _sigma,
                        RooAbsReal& _alpha,
                        RooAbsReal& _alphaCB,
                        RooAbsReal& _nCB) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
    m0(" m0"," m0",this,_m0),
    sigma(" sigma"," sigma",this,_sigma),
    alpha(" alpha"," alpha",this,_alpha),
    alphaCB(" alphaCB"," alphaCB",this,_alphaCB),
    nCB(" nCB"," nCB",this,_nCB)
 { 
 } 


 RooCBCrujffPdf::RooCBCrujffPdf(const RooCBCrujffPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
    m0(" m0",this,other. m0),
    sigma(" sigma",this,other. sigma),
    alpha(" alpha",this,other. alpha),
    alphaCB(" alphaCB",this,other. alphaCB),
    nCB(" nCB",this,other. nCB)
 { 
 } 



 Double_t RooCBCrujffPdf::evaluate() const 
 { 
   double dx = (m-m0) ;
   double val;
   if(dx<0) { // Cruijff                                                                                                                                                   
     double f = 2*sigma*sigma + alpha*dx*dx ;
     val = exp(-dx*dx/f);
   } else { // CryBall                                                                                                                                                     
     if(dx/sigma < alphaCB) {
       val = exp(-dx*dx/2./sigma/sigma);
     } else {
       double A = pow(nCB/fabs(alphaCB),nCB)*exp(-pow(alphaCB,2.)/2.);
       double B = nCB/fabs(alphaCB)-fabs(alphaCB);
       val = A*pow(B+dx/sigma,-nCB);
     }
   }
   return val;
 } 



