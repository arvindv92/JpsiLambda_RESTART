#ifndef CutFinalBDT_H
#define CutFinalBDT_H
#endif
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>

using namespace std;
void CutFinalBDT(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0, Int_t trackType = 3, const char* isoVersion = "v1", Int_t isoConf = 1, Int_t bdtConf = 1, Float_t bdtCut = -1.0, Float_t bdtCut_ZeroTracks = -1.0, Bool_t isoFlag = true, Bool_t newFlag = false, Bool_t logFlag = false);
