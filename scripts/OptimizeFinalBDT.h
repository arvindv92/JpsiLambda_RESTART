#ifndef OptimizeFinalBDT_H
#define OptimizeFinalBDT_H

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TList.h"
#include "TIterator.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include <iostream>
#include <vector>
using namespace std;

std::vector <Double_t> OptimizeFinalBDT(Int_t run = 1, const char* isoVersion = "v1",
                                        Int_t isoConf = 1, Int_t bdtConf = 1,
                                        Bool_t isoFlag = true,
                                        Bool_t logFlag = false,
                                        const char* FOM = "Sig",
                                        const char *part = "lambda");

#endif
