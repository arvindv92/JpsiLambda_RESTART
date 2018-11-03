#ifndef collateFiles_H
#define collateFiles_H
#include <TChain.h>
#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TSystem.h>
#include <TCollection.h>
#include <TFileCollection.h>
void collateFiles(Int_t run = 1, Int_t isData = 1, Int_t mcType = 1, TChain** h1 = NULL, TChain** h2 = NULL);

#endif
