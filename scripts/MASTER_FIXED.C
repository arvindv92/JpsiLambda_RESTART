#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>

using namespace std;

void MASTER_FIXED(Int_t year)
{
	Bool_t testing         = false; // when true, analysis will only run over a subset of data
	Bool_t loose           = true;  // when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run              = 2;     // Run 1 or Run 2?
	Bool_t isData          = true;  // Data or MC?
	Int_t mcType           = 0;     // 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;     // 3 for LL, 5 for DD.
	Bool_t logFlag         = true;  // set to false only while testing.

	Int_t runArray[2] = {1,2};
	Bool_t isDataArray[2] = {true, false};

	for(Int_t i = 1; i<=1; i++)
	{
		run = runArray[i];
		//Data
		//Trigger Cut
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L Trigger.C++");
		gROOT->ProcessLine(Form("Trigger(%d, %d, %d, %d, %d, %d, %d)",run, year, isData, mcType, testing, loose, logFlag));
		gROOT->Reset();

		//Sanity Cuts
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L Sanity.C++");
		gROOT->ProcessLine(Form("Sanity(%d, %d, %d, %d, %d)", run, year, isData, mcType, logFlag));
		gROOT->Reset();

		//Cut Out Ks0
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L CutOutKs.C++");
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
