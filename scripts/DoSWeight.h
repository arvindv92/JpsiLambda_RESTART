#ifndef DoSWeight_H
#define DoSWeight_H
#endif
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
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

void DoSWeight(Int_t run = 1, Int_t trackType = 3, Bool_t logFlag = false, Bool_t zeroFlag = false);

void AddModel(RooWorkspace* ws = nullptr, Int_t lowRange = 5200, Int_t highRange = 6000, Int_t nEntries = 0);

void AddData(RooWorkspace* ws = nullptr, Int_t run = 1, TTree *treeIn = nullptr);

void DosPlot(RooWorkspace* ws = nullptr, Int_t run = 1, const char* type = "LL", TTree *treeOut = nullptr, TTree *treeOut_training = nullptr, Bool_t zeroFlag = false);
