#ifndef TrainIsolation_H
#define TrainIsolation_H
#endif
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
#endif

using namespace std;

void TrainIsolation(Int_t run = 1,Int_t trackType = 3, const char* isoVersion = "v1", Bool_t logFlag = false);
