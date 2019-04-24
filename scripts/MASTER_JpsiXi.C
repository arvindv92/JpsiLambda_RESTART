#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include "Trigger_JpsiXi.h"
#include "Cuts_JpsiXi.h"
#include <iostream>

using namespace std;

void MASTER_JpsiXi()
{
	TStopwatch sw;
	sw.Start();

	gSystem->Exec("date");

	gROOT->ProcessLine(".L CollateFiles_JpsiXi.C+");
	gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
	gROOT->ProcessLine(".L Cuts_JpsiXi.C+");

	gSystem->Load("CollateFiles_JpsiXi_C.so");
	gSystem->Load("Trigger_JpsiXi_C.so");
	gSystem->Load("Cuts_JpsiXi_C.so");

	Bool_t testing         = false;// when true, analysis will only run over a subset of data
	Bool_t loose           = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run              = 1;// Run 1 or Run 2?
	Bool_t isData          = false;// Data or MC?
	Bool_t logFlag         = true;// set to false only while testing.
	Int_t year            = 2011;
	Int_t runArray[2]          = {1,2};

	/*  //  Data- Run 1
	   cout<<"Processing Run 1"<<endl;
	   cout<<"Processing 2011"<<endl;
	   //  Trigger Cut
	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
	   gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
	                        run, year, isData, testing, loose, logFlag));
	   gROOT->Reset();

	   //Sanity Cuts
	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
	   gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
	   gROOT->Reset();


	   cout<<"Processing 2012"<<endl;
	   year = 2012;
	   //  Trigger Cut
	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
	   gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
	                        run, year, isData, testing, loose, logFlag));
	   gROOT->Reset();

	   //Sanity Cuts
	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
	   gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
	   gROOT->Reset();


	 */
	//  Simulation

	for(Int_t i = 1; i<=1; i++)
	{
		run = runArray[i];
		cout<<"Processing Run "<<run<<endl;

		//  Trigger Cut
		cout<<"***Trigger***"<<endl;
		Trigger_JpsiXi(run, year, isData, testing, loose, logFlag);

		//Sanity Cuts
		cout<<"***Sanity***"<<endl;
		Cuts_JpsiXi(run, year, isData, logFlag);
	}
	//********************************************************************



	sw.Stop();
	cout << "==> MASTER_JpsiXi is done! Huzzah!: "; sw.Print();
}
