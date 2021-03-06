/*****************************************************************************
 * Project: RooFit                                                           *
 * Package: RooFitModels                                                     *
 *    File: $Id$
 * Authors:                                                                  *
 *   DMS, Diego Martinez Santos, Nikhef, Diego.Martinez.Santos@cern.ch       *
 *                                                                           *
 * Copyright (c) 2013, Nikhef. All rights reserved.                          *
 *                                                                           *
 * Redistribution and use in source and binary forms,                        *
 * with or without modification, are permitted according to the terms        *
 * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             *
 *****************************************************************************/

#ifndef V0HHFITTER_ROOHYPATIA2_H
#define V0HHFITTER_ROOHYPATIA2_H 1

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooAbsReal;
 
class RooHypatia2 : public RooAbsPdf {
public:
  RooHypatia2() {} ; 
  RooHypatia2(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _l,
	      RooAbsReal& _zeta,
	      RooAbsReal& _fb,
	      RooAbsReal& _sigma,
	      RooAbsReal& _mu,
	      RooAbsReal& _a,
	      RooAbsReal& _n,
	      RooAbsReal& _a2,
	      RooAbsReal& _n2);
  RooHypatia2(const RooHypatia2& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new RooHypatia2(*this,newname); }
  inline virtual ~RooHypatia2() { }

protected:

  RooRealProxy x ;
  RooRealProxy l ;
  RooRealProxy zeta ;
  RooRealProxy fb ;
  RooRealProxy sigma ;
  RooRealProxy mu ;
  RooRealProxy a ;
  RooRealProxy n ;
  RooRealProxy a2 ;
  RooRealProxy n2 ;
  
  Double_t evaluate() const ;

private:

  ClassDef(RooHypatia2,1) 
};
 
#endif // V0HHFITTER_ROOHYPATIA2_H
