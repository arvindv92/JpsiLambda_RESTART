#ifndef Cuts_JpsiXi_H
#define Cuts_JpsiXi_H
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

void Cuts_JpsiXi(Int_t run = 1, Int_t year = 2011, Bool_t isData = true,
                 Bool_t logFlag = false);
#endif
