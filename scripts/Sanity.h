#ifndef Sanity_H
#define Sanity_H
#endif

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
using namespace std;

void Sanity(Int_t run = 1, Int_t year = 2011, Bool_t isData = true, Int_t mcType = 0, Bool_t logFlag = false);
