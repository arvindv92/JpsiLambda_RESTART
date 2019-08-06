#ifndef ApplyFinalBDT_H
#define ApplyFinalBDT_H

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
using namespace std;

void ApplyFinalBDT(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0,
                   Int_t trackType = 3, const char* isoVersion = "v1", Int_t isoConf = 1,
                   Int_t bdtConf = 1, Int_t flag = 1, Bool_t isoFlag = true,
                   Bool_t zeroFlag = false, Bool_t logFlag = false,
                   Bool_t simFlag = false);

#endif
