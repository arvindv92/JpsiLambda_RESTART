/********************************
   Author : Aravindhan V.
 *********************************/
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TROOT.h>
#include <iostream>

using namespace std;
void massCut(Int_t run = 1, Int_t trackType = 3, Int_t low = 5200, Int_t high = 6000)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   high and low refer to the mass cut you want to make
 */
{
	Bool_t logFlag = 0;
	const char *logFileName = "", *type = "", *masscut = "";
	Double_t Lb_DTF_M_JpsiLConstr = 0.;
	Int_t entries_init = 0, entries_final = 0;

	masscut = (TString::Format("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d",low,high)).Data();

	if(trackType == 3)
	{
		cout<<"Processing LL"<<endl;
		logFileName = "masscut_LL_log.txt";
		type = "LL";
	}
	else if(trackType == 5)
	{
		cout<<"Processing DD"<<endl;
		logFileName = "masscut_DD_log.txt";
		type = "DD";
	}
	TFile *filein(0), *fileout(0);
	TTree *treein(0), *treeout(0);

	if(logFlag == 1)
	{
		gROOT->ProcessLine((TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName)).Data());
	}

	filein = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type),"READ");
	fileout = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_masscut_%s.root",run,type),"RECREATE");
	treein = (TTree*)filein->Get("MyTuple");

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<filein->GetName()<<endl;
	cout<<"Output file = "<<fileout->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treein->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	cout<<"I am making the following mass cuts. Sit tight"<<endl;
	cout<<masscut<<endl;

	treein->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&Lb_DTF_M_JpsiLConstr);
	treeout = (TTree*)treein->CloneTree(0);

	for (Int_t i = 0; i < entries_init; i++)
	{
		treein->GetEntry(i);

		if(Lb_DTF_M_JpsiLConstr > low && Lb_DTF_M_JpsiLConstr < high)
		{
			treeout->Fill();
		}
	}

	entries_final = treeout->GetEntries();
	cout<<"Outgoing entries = "<<entries_final<<endl;

	if(logFlag == 1) gROOT->ProcessLine(".>");

	fileout->Write();
	fileout->Close();
	filein->Close();
}
