#ifndef TrainFinalBDT_H
#define TrainFinalBDT_H

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TROOT.h"
//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/CrossValidation.h"
#endif

using namespace std;
void TrainFinalBDT(Int_t run = 1, const char* isoVersion = "v1", Int_t isoConf = 1,
                   Bool_t isoFlag = true, Int_t bdtConf = 1, Bool_t logFlag = false,
                   Bool_t simFlag = false);

#endif
