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

#include "TMath.h"
#include "Math/SpecFunc.h"

#include "RooAbsReal.h" 
#include "gsl/gsl_sf_bessel.h"
#include <math.h> 
#include "V0hhFitter/RooHypatia2.h" 
#include "Math/SpecFuncMathMore.h"

const Double_t sq2pi = TMath::Sqrt(2.*TMath::ACos(-1.));
const Double_t sq2pi_inv = 1./sq2pi;
const Double_t logsq2pi = TMath::Log(sq2pi);
const Double_t log_de_2 = TMath::Log(2.);

Double_t low_x_BK2(Double_t nu,Double_t x){
  return TMath::Gamma(nu)*TMath::Power(2.,nu-1.)*TMath::Power(x,-nu);
}


Double_t low_x_LnBK2(Double_t nu, Double_t x){
  return TMath::Log(TMath::Gamma(nu)) + (nu-1.)*log_de_2 - nu * TMath::Log(x);
}

Double_t BK2(Double_t ni, Double_t x) {
  Double_t nu = TMath::Abs(ni);
  if ( x < 1.e-06 && nu > 0.) return low_x_BK2(nu,x);
  if ( x < 1.e-04 && nu > 0. && nu < 55.) return low_x_BK2(nu,x);
  if ( x < 0.1 && nu >= 55.) return low_x_BK2(nu,x);

  return ROOT::Math::cyl_bessel_k(nu, x);
}

Double_t LnBK2(double ni, double x) {
  Double_t nu = TMath::Abs(ni);
  if ( x < 1.e-06 && nu > 0.) return low_x_LnBK2(nu,x);
  if ( x < 1.e-04 && nu > 0. && nu < 55.) return low_x_LnBK2(nu,x);
  if ( x < 0.1 && nu >= 55.) return low_x_LnBK2(nu,x);

  return TMath::Log(ROOT::Math::cyl_bessel_k(nu, x));
}


Double_t LogEval2(Double_t d, Double_t l, Double_t alpha, Double_t beta, Double_t delta) {

  Double_t gamma = alpha;
  Double_t dg = delta*gamma;
  Double_t thing = delta*delta + d*d;
  Double_t logno = l*TMath::Log(gamma/delta) - logsq2pi -LnBK2(l, dg);
  
  return TMath::Exp(logno + beta*d +(0.5-l)*(TMath::Log(alpha)-0.5*TMath::Log(thing)) + LnBK2(l-0.5,alpha*TMath::Sqrt(thing))); 
}


Double_t diff_eval2(Double_t d, Double_t l, Double_t alpha, Double_t beta, Double_t delta){

  Double_t gamma = alpha;
  Double_t dg = delta*gamma;
  Double_t thing = delta*delta + d*d;
  Double_t sqthing = TMath::Sqrt(thing);
  Double_t alphasq = alpha*sqthing;
  Double_t no = TMath::Power(gamma/delta,l)/BK2(l,dg)*sq2pi_inv;
  Double_t ns1 = 0.5-l;
  
  return no*TMath::Power(alpha, ns1)*TMath::Power(thing, l/2. - 1.25)*(-d*alphasq*(BK2(l - 1.5, alphasq) + BK2(l + 0.5, alphasq)) + (2.*(beta*thing + d*l) - d)*BK2(ns1, alphasq))*TMath::Exp(beta*d)/2.;
}

RooHypatia2::RooHypatia2(const char *name, const char *title, 
                        RooAbsReal& _x,
                        RooAbsReal& _l,
                        RooAbsReal& _zeta,
                        RooAbsReal& _fb,
                        RooAbsReal& _sigma,
                        RooAbsReal& _mu,
                        RooAbsReal& _a,
                        RooAbsReal& _n,
                        RooAbsReal& _a2,
                        RooAbsReal& _n2) :
   RooAbsPdf(name,title), 
   x("x","x",this,_x),
   l("l","l",this,_l),
   zeta("zeta","zeta",this,_zeta),
   fb("fb","fb",this,_fb),
   sigma("sigma","sigma",this,_sigma),
   mu("mu","mu",this,_mu),
   a("a","a",this,_a),
   n("n","n",this,_n),
   a2("a2","a2",this,_a2),
   n2("n2","n2",this,_n2)
{ 
} 


RooHypatia2::RooHypatia2(const RooHypatia2& other, const char* name) :  
   RooAbsPdf(other,name), 
   x("x",this,other.x),
   l("l",this,other.l),
   zeta("zeta",this,other.zeta),
   fb("fb",this,other.fb),
   sigma("sigma",this,other.sigma),
   mu("mu",this,other.mu),
   a("a",this,other.a),
   n("n",this,other.n),
   a2("a2",this,other.a2),
   n2("n2",this,other.n2)
{ 
} 



Double_t RooHypatia2::evaluate() const 
{ 
   Double_t d = x-mu;
   Double_t cons0 = TMath::Sqrt(zeta);
   Double_t alpha, beta, delta,  cons1, phi, A, B, k1, k2;
   Double_t asigma = a*sigma;
   Double_t a2sigma = a2*sigma;
   Double_t out = 0.;
   if (zeta!= 0.) {
     phi = BK2(l+1.,zeta)/BK2(l,zeta); 
     cons1 = sigma/TMath::Sqrt(phi);
     alpha  = cons0/cons1;
     beta = fb;
     delta = cons0*cons1;
     
     if (d < -asigma){
       
       k1 = LogEval2(-asigma,l,alpha,beta,delta);
       k2 = diff_eval2(-asigma,l,alpha,beta,delta); 
       B = -asigma + n*k1/k2;
       A = k1*TMath::Power(B+asigma,n);
       out = A*TMath::Power(B-d,-n);
     }
     else if (d>a2sigma) {
       k1 = LogEval2(a2sigma,l,alpha,beta,delta);
       k2 = diff_eval2(a2sigma,l,alpha,beta,delta);
       
       B = -a2sigma - n2*k1/k2;
       
       A = k1*TMath::Power(B+a2sigma,n2);
       
       out =  A*TMath::Power(B+d,-n2);
       
     }
     else {
       out = LogEval2(d,l,alpha,beta,delta);
     }
     


   }
   else if (l < 0.) {
     beta = fb;
     cons1 = -2.*l;
     delta = sigma;
     if (d < -asigma ) {
       cons1 = TMath::Exp(-beta*asigma);
       phi = 1. + a*a;
       k1 = cons1*TMath::Power(phi,l-0.5);
       k2 = beta*k1- cons1*(l-0.5)*TMath::Power(phi,l-1.5)*2*a/delta;
       B = -asigma + n*k1/k2;
       A = k1*TMath::Power(B+asigma,n);
       out = A*TMath::Power(B-d,-n);
     }
     else if (d > a2sigma) {
       cons1 = TMath::Exp(beta*a2sigma);
       phi = 1. + a2*a2;
       k1 = cons1*TMath::Power(phi,l-0.5);
       k2 = beta*k1+ cons1*(l-0.5)*TMath::Power(phi,l-1.5)*2.*a2/delta;
       B = -a2sigma - n2*k1/k2;
       A = k1*TMath::Power(B+a2sigma,n2);
       out =  A*TMath::Power(B+d,-n2);
       
     }
     else { out = TMath::Exp(beta*d)*TMath::Power(1. + d*d/(delta*delta),l-0.5);}
   }
   else {
   }
   return out;
}
