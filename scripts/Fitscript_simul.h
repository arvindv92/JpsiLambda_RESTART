#ifndef Fitscript_simul_H
#define Fitscript_simul_H
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
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooOneSidedProfileLL.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/TestStatSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/TestStatistic.h"
#include "RooHypatia2.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TRolke.h>
#include "TROOT.h"
#include "TIterator.h"
#include "TLatex.h"

void Fitscript_simul(const char *option = "best", Int_t myLow = 4700, Int_t myHigh = 6000, Int_t Lst1405_rwtype = 0, Float_t bdtCut = 0.0);

// void Fitscript_simul(Int_t finalBDTConf_nonZero = 1, Int_t finalBDTConf_Zero = 1,
//                      Int_t isoConf = 2, const char *isoVersion = "v0",
//                      Float_t bdtCut_nonZero = 0.535, Float_t bdtCut_Zero = 0.415);

#endif
