#ifndef CollateFiles_H
#define CollateFiles_H

#include <TChain.h>
#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TCollection.h>
#include <TFileCollection.h>
void CollateFiles(Int_t run = 1, Int_t year = 2011, Bool_t isData = true,
                  Int_t mcType = 1, TChain** h1 = nullptr, TChain** h2 = nullptr,
                  Bool_t testing = false, Bool_t logFlag = true);
#endif
