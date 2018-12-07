/********************************
   Author : Aravindhan V.
 *********************************/
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>

using namespace std;
void massCut(Int_t run = 1, Int_t trackType = 3, Int_t low = 5200, Int_t high = 6000, TString inFileName = "", TString outFileName = "", Bool_t logFlag = false)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   high and low refer to the mass cut you want to make
 */
{
	cout<<"******Starting massCut*********"<<endl;

	TStopwatch sw;
	sw.Start();

	Double_t Lb_DTF_M_JpsiLConstr = 0.;
	Int_t entries_init = 0, entries_final = 0;
	TFile *fileIn(0), *fileOut(0);
	TTree *treeIn(0), *treeOut(0);
	const char *logFileName = "", *masscut = "";

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

	masscut = (TString::Format("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d",low,high)).Data();

	if(trackType == 3)
	{
		cout<<"Processing LL run"<<run<<endl;
		logFileName = "masscut_LL_log.txt";
	}
	else if(trackType == 5)
	{
		cout<<"Processing DD run"<<run<<endl;
		logFileName = "masscut_DD_log.txt";
	}

	if(logFlag)
	{
		gROOT->ProcessLine((TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName)).Data());
	}

	fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/%s",run,inFileName.Data()),"READ");
	fileOut = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/%s",run,outFileName.Data()),"RECREATE");
	treeIn = (TTree*)fileIn->Get("MyTuple");

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	cout<<"I am making the following mass cuts. Sit tight"<<endl;
	cout<<masscut<<endl;

	treeIn->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&Lb_DTF_M_JpsiLConstr);
	treeOut = (TTree*)treeIn->CloneTree(0);

	for (Int_t i = 0; i < entries_init; i++)
	{
		treeIn->GetEntry(i);

		if(Lb_DTF_M_JpsiLConstr > low && Lb_DTF_M_JpsiLConstr < high)
		{
			treeOut->Fill();
		}
	}

	entries_final = treeOut->GetEntries();
	cout<<"Outgoing entries = "<<entries_final<<endl;

	if(logFlag) gROOT->ProcessLine(".>");

	fileOut->Write();
	fileOut->Close();
	fileIn->Close();

	sw.Stop();
	cout << "==> End of massCut! Cheers!: "; sw.Print();
}
