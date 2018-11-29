#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>

using namespace std;

void MASTER_FIXED()
{
	Bool_t testing         = false; // when true, analysis will only run over a subset of data
	Bool_t loose           = true;  // when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run              = 1;     // Run 1 or Run 2?
	Bool_t isData          = true;  // Data or MC?
	Int_t mcType           = 0;     // 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;     // 3 for LL, 5 for DD.
	Bool_t logFlag         = true;  // set to false only while testing.

	Int_t runArray[2] = {1,2};
	Bool_t isDataArray[2] = {true, false};

	for(Int_t i = 0; i<=1; i++)
	{
		run = runArray[i];
		//Data
		//Trigger Cut
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L Trigger.C++");
		gROOT->ProcessLine(TString::Format("Trigger(%d, %d, %d, %d, %d, %d)",run, isData, mcType, testing, loose, logFlag).Data());
		gROOT->Reset();

		//Sanity Cuts
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L Sanity.C++");
		gROOT->ProcessLine(TString::Format("Sanity(%d, %d, %d, %d)", run, isData, mcType, logFlag).Data());
		gROOT->Reset();

		//Cut Out Ks0
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L CutOutKs.C++");
		gROOT->ProcessLine(TString::Format("CutOutKs(%d, %d, %d, %d, %d)", run, isData, mcType, trackType, logFlag).Data());
		gROOT->Reset();

		//sWeight data
		if(isData)
		{
			gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			gROOT->ProcessLine(".L DosPlot.C++");
			gROOT->ProcessLine(TString::Format("DosPlot(%d, %d, %d)", run, trackType, logFlag).Data());
			gROOT->Reset();
		}
	}

}
