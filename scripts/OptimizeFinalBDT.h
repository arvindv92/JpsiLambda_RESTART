#ifndef OptimizeFinalBDT_H
#define OptimizeFinalBDT_H
#endif

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TList.h"
#include "TIterator.h"
#include <iostream>
#include <vector>
using namespace std;

std::vector <Double_t> OptimizeFinalBDT(Int_t run = 1, Int_t trackType = 3,
		      const char* isoVersion = "v1", Int_t isoConf = 1, Int_t bdtConf = 1,
		      Bool_t isoFlag = true, Bool_t logFlag = false, Bool_t newFlag = false,TString FOM = "sig");
