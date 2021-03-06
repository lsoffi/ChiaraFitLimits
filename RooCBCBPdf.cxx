/***************************************************************************** 
 * Project: RooFit                                                           * 
 *                                                                           * 
 * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/ 

// Your description goes here... 

#include "Riostream.h" 

#include "RooCBCBPdf.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <math.h> 
#include "TMath.h" 

ClassImp(RooCBCBPdf) 

 RooCBCBPdf::RooCBCBPdf(const char *name, const char *title, 
                        RooAbsReal& _m,
                        RooAbsReal& _m0,
                        RooAbsReal& _sigma,
                        RooAbsReal& _alpha1,
                        RooAbsReal& _alpha2,
                        RooAbsReal& _n1,
                        RooAbsReal& _n2) :
   RooAbsPdf(name,title), 
   m("m","m",this,_m),
    m0(" m0"," m0",this,_m0),
    sigma(" sigma"," sigma",this,_sigma),
    alpha1(" alpha1"," alpha1",this,_alpha1),
    alpha2(" alpha2"," alpha2",this,_alpha2),
    n1(" n1"," n1",this,_n1),
    n2(" n2"," n2",this,_n2)
 { 
 } 


 RooCBCBPdf::RooCBCBPdf(const RooCBCBPdf& other, const char* name) :  
   RooAbsPdf(other,name), 
   m("m",this,other.m),
    m0(" m0",this,other.m0),
    sigma(" sigma",this,other.sigma),
    alpha1(" alpha1",this,other.alpha1),
    alpha2(" alpha2",this,other.alpha2),
    n1(" n1",this,other.n1),
    n2(" n2",this,other.n2)
 { 
 } 



 Double_t RooCBCBPdf::evaluate() const 
 { 
   double dx = (m-m0) ;
   double val;
   if(dx>0) {          
      if(dx/sigma < alpha1) {
       val = exp(-dx*dx/2./sigma/sigma);
       } else {
       double A = pow(n1/fabs(alpha1),n1)*exp(-pow(alpha1,2.)/2.);
       double B = n1/fabs(alpha1)-fabs(alpha1);
       val = A*pow(B+dx/sigma,-n1);
       }
   }else{          
     if(dx/sigma > alpha2) {
       val = exp(-dx*dx/2./sigma/sigma);
     } else {
       double A = pow(n2/fabs(alpha2),n2)*exp(-pow(alpha2,2.)/2.);
       double B = n2/fabs(alpha2)-fabs(alpha2);
       val = A*pow(B+dx/sigma,-n2);
       }
   }
   return val;
 } 



