#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>

using namespace std;

void MASTER_FIXED()
{
  Bool_t testing         = false;// when true, analysis will only run over a subset of data
  Bool_t loose           = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
  Int_t run              = 1;// Run 1 or Run 2?
  Bool_t isData          = false;// Data or MC?
  Int_t mcType           = 1;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
  Int_t trackType        = 3;// 3 for LL, 5 for DD.
  Bool_t logFlag         = true;// set to false only while testing.
  Int_t year             = 2011;

  Int_t runArray[2] = {1,2};
  Int_t mcTypeArray[3] = {1,2,3};
  Bool_t isDataArray[2] = {true, false};

  for(Int_t i = 0; i<=2; i++)//Loop over different MC types
    {
      for(Int_t j = 0; j<=1; j++)//Loop over runs
	{
	  mcType = mcTypeArray[i];
	  run = runArray[j];

	  cout<<"*******PROCESSING RUN "<<run<<" MC TYPE "<<mcType<<"********"<<endl;
	  //MC
	  //Trigger Cut
	  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	  gROOT->ProcessLine(".L Trigger.C+");
	  gROOT->ProcessLine(Form("Trigger(%d, %d, %d, %d, %d, %d, %d)",run, year, isData, mcType, testing, loose, logFlag));
	  gROOT->Reset();

	  //Sanity Cuts
	  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	  gROOT->ProcessLine(".L Sanity.C+");
	  gROOT->ProcessLine(Form("Sanity(%d, %d, %d, %d, %d)", run, year, isData, mcType, logFlag));
	  gROOT->Reset();

	  //Cut Out Ks0
	  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	  gROOT->ProcessLine(".L CutOutKs.C+");
	  gROOT->ProcessLine(Form("CutOutKs(%d, %d, %d, %d, %d, %d)", run, year, isData, mcType, trackType, logFlag));
	  gROOT->Reset();

	  //sWeight data
	  // if(isData)
	  // {
	  //      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	  //      gROOT->ProcessLine(".L DosPlot.C++");
	  //      gROOT->ProcessLine(Form("DosPlot(%d, %d, %d)", run, trackType, logFlag));
	  //      gROOT->Reset();
	  // }
	}
    }

}
