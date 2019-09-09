#ifndef CutOutKs_H
#define CutOutKs_H

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
using namespace std;

void CutOutKs(Int_t run = 1, Int_t year = 2011, Bool_t isData = true, Int_t mcType = 0, Bool_t logFlag = false);

#endif
