/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/root/5.32.00-cms5/include/Riostream.h" 

#include "RooExpolPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooExpolPdf) 

 RooExpolPdf::RooExpolPdf(const char *name, const char *title, 
			  RooAbsReal& _m,
			  RooAbsReal& _p1,
			  RooAbsReal& _p2) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
   p1("p1","p1",this,_p1),
   p2("p2","p2",this,_p2)
   
 { 
 } 


 RooExpolPdf::RooExpolPdf(const RooExpolPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
   p1("p1",this,other.p1),
   p2("p2",this,other.p2) 
 { 
 } 



 Double_t RooExpolPdf::evaluate() const 
 { 
   double val;
   val = exp(-m/(p1+p2*m));
   return val;
 } 



