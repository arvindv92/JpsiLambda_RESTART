#ifndef CollateFiles_H
#define CollateFiles_H
#include <TChain.h>
#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TCollection.h>
#include <TFileCollection.h>
void CollateFiles(Int_t run = 1, Bool_t isData = true, Int_t mcType = 1, TChain** h1 = NULL, TChain** h2 = NULL, Bool_t testing = false, Bool_t loose = true);

#endif
