#ifndef DoSWeight_Sanity_H
#define DoSWeight_Sanity_H

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
//#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TPaveLabel.h"
#include "TStopwatch.h"

Double_t DoSWeight_Sanity(Int_t run = 1, Int_t trackType = 3,
                          Bool_t logFlag = false);
#endif
