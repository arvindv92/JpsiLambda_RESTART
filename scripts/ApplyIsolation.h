#ifndef ApplyIsolation_H
#define ApplyIsolation_H

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

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;
using namespace std;

void ApplyIsolation(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0,
                    Int_t trackType = 3, Int_t flag = 1, const char* isoVersion = "v1",
                    Int_t isoConf = 1, Bool_t logFlag = false);
#endif
