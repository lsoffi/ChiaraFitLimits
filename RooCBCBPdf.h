/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOCBCBPDF
#define ROOCBCBPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class RooCBCBPdf : public RooAbsPdf {
public:
  RooCBCBPdf() {} ; 
  RooCBCBPdf(const char *name, const char *title,
	      RooAbsReal& _m,
	      RooAbsReal& _m0,
	      RooAbsReal& _sigma,
	      RooAbsReal& _alpha1,
	      RooAbsReal& _alpha2,
	      RooAbsReal& _n1,
	      RooAbsReal& _n2);
  RooCBCBPdf(const RooCBCBPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooCBCBPdf(*this,newname); }
  inline virtual ~RooCBCBPdf() { }

protected:

  RooRealProxy m ;
  RooRealProxy  m0 ;
  RooRealProxy  sigma ;
  RooRealProxy  alpha1 ;
  RooRealProxy  alpha2 ;
  RooRealProxy  n1 ;
  RooRealProxy  n2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooCBCBPdf,1) // Your description goes here...
};
 
#endif
