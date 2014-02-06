/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOEXPOLPDF
#define ROOEXPOLPDF
#include "RooFit.h"
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooExpolPdf : public RooAbsPdf {
public:
  RooExpolPdf() {} ; 
  RooExpolPdf(const char *name, const char *title,
	      RooAbsReal& _m,
	      RooAbsReal& _p1,
	      RooAbsReal& _p2);
  RooExpolPdf(const RooExpolPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooExpolPdf(*this,newname); }
  inline virtual ~RooExpolPdf() { }

protected:

  RooRealProxy m;
  RooRealProxy p1;
  RooRealProxy p2;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooExpolPdf,1) // Your description goes here...
};
 
#endif
