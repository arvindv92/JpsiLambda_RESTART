#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include "Trigger.h"
#include "Sanity.h"
#include "CutOutKs.h"
#include <iostream>

using namespace std;

void MASTER_FIXED()
{
	TStopwatch sw;
	sw.Start();

	gSystem->Exec("date");

	gROOT->ProcessLine(".L CollateFiles.C+");
	gROOT->ProcessLine(".L Trigger.C+");
	gROOT->ProcessLine(".L Sanity.C+");
	gROOT->ProcessLine(".L CutOutKs.C+");

	gSystem->Load("CollateFiles_C.so");
	gSystem->Load("Trigger_C.so");
	gSystem->Load("Sanity_C.so");
	gSystem->Load("CutOutKs_C.so");

	Bool_t testing         = false;// when true, analysis will only run over a subset of data
	Bool_t loose           = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run              = 2;// Run 1 or Run 2?
	Bool_t isData          = true;// Data or MC?
	Int_t mcType           = 0;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;// 3 for LL, 5 for DD.
	Bool_t logFlag         = true;// set to false only while testing.
	Int_t year             = 2018;

	Int_t runArray[2] = {1,2};
	Int_t yearArray[3] = {2015,2016,2017};
	// Int_t mcTypeArray[2] = {4,5};
	// Bool_t isDataArray[2] = {true, false};

	// for(Int_t i = 0; i<=2; i++)//Loop over different years
	//   {
	//     year = yearArray[i];
	//     for(Int_t j = 1; j<=1; j++)//Loop over runs
	//      {
	// mcType = mcTypeArray[i];
	//  run = runArray[j];

	cout<<"*******PROCESSING RUN "<<run<<" DATA YEAR "<<year<<"********"<<endl;
	//MC
	//Trigger Cut
	cout<<"***Trigger***"<<endl;
	Trigger(run, year, isData, mcType, testing, loose, logFlag);
	//Sanity Cuts
	cout<<"***Sanity***"<<endl;
	Sanity(run, year, isData, mcType, logFlag);
	//Cut Out Ks0
	cout<<"***CutOutKs***"<<endl;
	CutOutKs(run, year, isData, mcType, trackType, logFlag);
	//sWeight data
	// if(isData)
	// {
	//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	//      gROOT->ProcessLine(".L DosPlot.C++");
	//      gROOT->ProcessLine(Form("DosPlot(%d, %d, %d)", run, trackType, logFlag));
	//      gROOT->Reset();
	// }
	//	}
	//    }
	sw.Stop();
	cout << "==> MASTER_FIXED is done! Huzzah!: "; sw.Print();
}
