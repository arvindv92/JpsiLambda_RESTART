#ifndef Trigger_H
#define Trigger_H
#endif

#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "CollateFiles.h"
#include <iostream>
#include <fstream>

using namespace std;

void Trigger(Int_t run = 1, Int_t year = 2011, Bool_t isData = true, Int_t mcType = 0, Bool_t testing = false, Bool_t loose = true, Bool_t logFlag = false);
