#ifndef Trigger_JpsiXi_H
#define Trigger_JpsiXi_H
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "CollateFiles_JpsiXi.h"
#include <iostream>
#include <fstream>

using namespace std;

void Trigger_JpsiXi(Int_t run = 1, Int_t year = 2011, Bool_t isData = true,
                    Bool_t testing = false, Bool_t logFlag = false);

#endif
