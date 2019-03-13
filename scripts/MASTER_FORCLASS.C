#include "Trigger.h"
#include "Sanity.h"
#include "CutOutKs.h"
#include "DoSWeight.h"
#include "TrainIsolation.h"
#include "ApplyIsolation.h"
#include "TrainFinalBDT.h"
#include "ApplyFinalBDT.h"
#include "OptimizeFinalBDT.h"
#include "CutFinalBDT.h"

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include <iostream>
#include <vector>

using namespace std;

void MASTER()
{
  TStopwatch sw;
  sw.Start();//starts timing your script

  gSystem->Exec("date");//prints out date and time. Very useful.
  
  gROOT->ProcessLine(".L CollateFiles.C+");//load and compile scripts. 
  gROOT->ProcessLine(".L Trigger.C+");
  gROOT->ProcessLine(".L Sanity.C+");
  gROOT->ProcessLine(".L CutOutKs.C+");

  gROOT->ProcessLine(".L DoSWeight.C+");
  gROOT->ProcessLine(".L TrainIsolation.C+");
  gROOT->ProcessLine(".L ApplyIsolation.C+");
  gROOT->ProcessLine(".L TrainFinalBDT.C+");
  gROOT->ProcessLine(".L ApplyFinalBDT.C+");
  gROOT->ProcessLine(".L OptimizeFinalBDT.C+");
  //  gROOT->ProcessLine(".L CutFinalBDT.C+");

  gSystem->Load("DoSWeight_C.so");//load dynamic libraries
  gSystem->Load("TrainIsolation_C.so");
  gSystem->Load("ApplyIsolation_C.so");
  gSystem->Load("TrainFinalBDT_C.so");
  gSystem->Load("ApplyFinalBDT_C.so");
  gSystem->Load("OptimizeFinalBDT_C.so");
  //  gSystem->Load("CutFinalBDT_C.so");

  Int_t run              = 2;// Run 1 or Run 2?
  Int_t mcType           = 0;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
  Int_t trackType        = 3;// 3 for LL, 5 for DD.
  Int_t isoConf          = 1;// config for isolation BDT. Currently 1 or 2 supported
  Int_t finalBDTconf     = 1;// config for final BDT. Currently 1 or 2 supported

  Bool_t testing         = false;// when true, analysis will only run over a subset of data
  Bool_t loose           = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
  Bool_t isData          = true;// Data or MC?
  Bool_t isoFlag         = true;// when true, isolation will be used in final BDT.
  Bool_t logFlag         = true;// set to false only while testing.

  Float_t bdtCut            = 0.0;
  Float_t bdtCut_ZeroTracks = 0.0;

  Double_t myChi2_nonZero = 0.0, myChi2_Zero = 0.0;
  const char* isoVersion  = ""; // which version of isolation BDT? changes input variables used in isolation. v0,1 supported.

  TString FOM                = "Punzi";//Sig for S/sqrt(S+B), Punzi for Punzi FOM with a = 3
  Int_t runArray[2]          = {1,2};
  Int_t isoConfArray[2]      = {1,2};
  Int_t finalBDTconfArray[2] = {1,2};
  Int_t mcTypeArray[5]       = {1,2,3,4,5};
  TString isoVersionArray[2] = {"v0","v1"};

  std::vector <Double_t> bdtCuts;

  for(Int_t h = 0; h <= 4; h++) //Loop over different MC Types
    {
      mcType = mcTypeArray[h];
      cout<<"$$$$$$$$$$$ Processing MC Type "<<mcType<< " $$$$$$$$$$$$$$"<<endl;
      for(Int_t i = 1; i<=2; i++)//Loop over different runs
	{
	  run = i;
	  //Trigger Cut
	  cout<<"***Trigger***"<<endl;
	  Trigger(run, year, isData, mcType, testing, loose, logFlag);
	  //Sanity Cuts
	  cout<<"***Sanity***"<<endl;
	  Sanity(run, year, isData, mcType, logFlag);
	  //Cut Out Ks0
	  cout<<"***CutOutKs***"<<endl;
	  CutOutKs(run, year, isData, mcType, trackType, logFlag);
	}
      //sWeight data
      if(isData)
	{
	  //sWeight nonZeroTracks data
	  cout<<"****DoSWeight run "<<run<<" nonZeroTracks"<<endl;
	  
	  myChi2_nonZero = DoSWeight(run, trackType, logFlag, false);

	  if(myChi2_nonZero > 2.0)
	    {
	      cout<<"####Fit failed for nonZeroTracks sWeighting. Exiting####"<<endl;
	      exit(1);
	    }
	  //sWeight ZeroTracks data
	  cout<<"****DoSWeight run "<<run<<" ZeroTracks"<<endl;
	  
	  myChi2_Zero = DoSWeight(run, trackType, logFlag, true);
	  
	  if(myChi2_Zero > 2.0)
	    {
	      cout<<"####Fit failed for ZeroTracks sWeighting. Exiting####"<<endl;
	      exit(1);
	    }
	}
    }
  sw.Stop();
  cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
