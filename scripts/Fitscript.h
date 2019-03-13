#ifndef Fitscript_H
#define Fitscript_H
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLegend.h"
#include "TSystem.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
//#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooKeysPdf.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooDataHist.h"
#include "TError.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TRolke.h>
using namespace RooFit;
using namespace std;

std::vector <Float_t> Fitscript( Int_t run = 2, Int_t finalBDTConf = 1, Int_t isoConf = 2,
                                 const char *isoVersion = "v0",
                                 Int_t mylow = 4700, Int_t constraintflag = 1,
                                 Float_t bdtCut_nonZero = 0.535, Float_t bdtCut_Zero = 0.415);

#endif
