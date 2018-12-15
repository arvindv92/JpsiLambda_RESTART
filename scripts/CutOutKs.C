/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply cuts to remove the Ks0 -> Lambda0 mis-ID reflection.
 *********************************/
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
using namespace std;

void CutOutKs(Int_t run = 1, Int_t year = 2011, Bool_t isData = true, Int_t mcType = 0, Int_t trackType = 3, Bool_t logFlag = false)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{
	cout<<"***********Starting CutOutKs***********"<<endl;

	TStopwatch sw;
	sw.Start();

	TFile *fileIn = nullptr, *fileOut = nullptr, *fileOut_nonZeroTracks = nullptr, *fileOut_ZeroTracks = nullptr;;
	TTree *treeIn = nullptr, *treeOut = nullptr, *treeOut_nonZeroTracks = nullptr, *treeOut_ZeroTracks = nullptr;

	const char *L_dmcut = "", *pidcut = "", *massWindow = "", *type = "";
	Int_t entries_init = 0, entries_final = 0, entries_final_nonZeroTracks = 0, entries_final_ZeroTracks = 0, entries_gen = 0, nTracks = 0;
	Float_t eff_excl = 0.0, eff_excl_err = 0.0, eff_incl = 0.0, eff_incl_err = 0.0;
	Double_t L_dm = 0.0, p_PIDp = 0.0, Lb_DTF_L_WMpipi_JpsiConstr = 0.0;
	Bool_t genFlag = false;
	TString logFileName = "";

	fstream genFile;

	massWindow  = "Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520";
	type        = (trackType == 3) ? ("LL") : ("DD");
	logFileName = (trackType == 3) ? TString::Format("cutoutks_%d_LL_log.txt",year) : TString::Format("cutoutks_%d_DD_log.txt",year);

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.
	cout<<"WD = "<<gSystem->pwd()<<endl;

	//Setup logging, input and output files.
	if(!isData)  // MC
	{
		switch(mcType)
		{
		case 1:         //JpsiLambda
			if(logFlag)
			{
				gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName.Data()).Data());
			}
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data()))
			{
				genFile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_sanity_%s.root",run,type),"READ");
			treeIn = (TTree*)fileIn->Get("MyTuple");
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type),"RECREATE");
			treeOut = (TTree*)treeIn->CloneTree(0);
			fileOut_nonZeroTracks = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_nonZeroTracks.root",run,type),"RECREATE");
			treeOut_nonZeroTracks = (TTree*)treeIn->CloneTree(0);
			fileOut_ZeroTracks = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_ZeroTracks.root",run,type),"RECREATE");
			treeOut_ZeroTracks = (TTree*)treeIn->CloneTree(0);
			break;

		case 2:         //JpsiSigma
			if(logFlag)
			{
				gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName.Data()).Data());
			}
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data()))
			{
				genFile.open((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_sanity_%s.root",run,type),"READ");
			treeIn = (TTree*)fileIn->Get("MyTuple");
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_cutoutks_%s.root",run,type),"RECREATE");
			treeOut = (TTree*)treeIn->CloneTree(0);
			fileOut_nonZeroTracks = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_cutoutks_%s_nonZeroTracks.root",run,type),"RECREATE");
			treeOut_nonZeroTracks = (TTree*)treeIn->CloneTree(0);
			fileOut_ZeroTracks = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_cutoutks_%s_ZeroTracks.root",run,type),"RECREATE");
			treeOut_ZeroTracks = (TTree*)treeIn->CloneTree(0);
			break;

		case 3:         //JpsiXi
			if(logFlag)
			{
				gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiXi/run%d/%s",run,logFileName.Data()).Data());
			}
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data()))
			{
				genFile.open((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_sanity_%s.root",run,type),"READ");
			treeIn = (TTree*)fileIn->Get("MyTuple");
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cutoutks_%s.root",run,type),"RECREATE");
			treeOut = (TTree*)treeIn->CloneTree(0);
			fileOut_nonZeroTracks = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cutoutks_%s_nonZeroTracks.root",run,type),"RECREATE");
			treeOut_nonZeroTracks = (TTree*)treeIn->CloneTree(0);
			fileOut_ZeroTracks = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cutoutks_%s_ZeroTracks.root",run,type),"RECREATE");
			treeOut_ZeroTracks = (TTree*)treeIn->CloneTree(0);
			break;
		}
	}
	else   // Data
	{
		if(logFlag)
		{
			gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName.Data()).Data());
		}
		fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_%s_%d.root",run,type,year),"READ");
		treeIn = (TTree*)fileIn->Get("MyTuple");
		fileOut = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_%d.root",run,type,year),"RECREATE");
		treeOut = (TTree*)treeIn->CloneTree(0);
		fileOut_nonZeroTracks = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_%d_nonZeroTracks.root",run,type,year),"RECREATE");
		treeOut_nonZeroTracks = (TTree*)treeIn->CloneTree(0);
		fileOut_ZeroTracks = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_%d_ZeroTracks.root",run,type,year),"RECREATE");
		treeOut_ZeroTracks = (TTree*)treeIn->CloneTree(0);
	}
	//end setup of input, output, logging
	cout<<"******************************************"<<endl;
	cout<<"==> Starting CutOutKs: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Processing Run "<<run<<" YEAR "<<year<<" "<<type<<((isData) ? (" Data") : (" MC type "))<<mcType<<endl;
	cout<<"******************************************"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"nonZeroTracks Output file = "<<fileOut_nonZeroTracks->GetName()<<endl;
	cout<<"ZeroTracks Output file = "<<fileOut_ZeroTracks->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("L_dm",&L_dm);
	treeIn->SetBranchAddress("p_PIDp",&p_PIDp);
	treeIn->SetBranchAddress("Lb_DTF_L_WMpipi_JpsiConstr",&Lb_DTF_L_WMpipi_JpsiConstr);
	treeIn->SetBranchAddress("Added_n_Particles",&nTracks);

	if(trackType == 3)
	{
		L_dmcut = "L_dm < 7.5";
		pidcut = "p_PIDp > 10";
		cout<<"Making LL file"<<endl;
	}
	else if(trackType == 5)
	{
		L_dmcut = "L_dm < 10";
		pidcut = "p_PIDp > 15";
		cout<<"Making DD file"<<endl;
	}

	cout<<"I am making the following cutoutks cuts only on events within 480-520 MeV in Lb_DTF_L_WMpipi_JpsiConstr. Events outside this window aren't cut on. Sit tight"<<endl;
	cout<<L_dmcut<<endl;
	cout<<pidcut<<endl;

	for (Int_t i = 0; i < entries_init; i++)
	{
		if(i % 100000 == 0) cout<<i<<endl;
		treeIn->GetEntry(i);
		if(Lb_DTF_L_WMpipi_JpsiConstr < 480 || Lb_DTF_L_WMpipi_JpsiConstr > 520)
		{
			if(nTracks > 0) treeOut_nonZeroTracks->Fill();
			else if(nTracks == 0) treeOut_ZeroTracks->Fill();
			treeOut->Fill();
			continue;
		}
		else
		{
			if(trackType == 3)
			{
				if (Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520)
				{
					if (p_PIDp > 10)
					{
						if(L_dm < 7.5)
						{
							if(nTracks > 0) treeOut_nonZeroTracks->Fill();
							else if(nTracks == 0) treeOut_ZeroTracks->Fill();
							treeOut->Fill();
						}
					}
				}
			}
			if(trackType == 5)
			{
				if (Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520)
				{
					if (p_PIDp > 15)
					{
						if(L_dm < 10)
						{
							if(nTracks > 0) treeOut_nonZeroTracks->Fill();
							if(nTracks == 0) treeOut_ZeroTracks->Fill();
							treeOut->Fill();
						}
					}
				}
			}
		}
	}

	entries_final = treeOut->GetEntries();
	entries_final_nonZeroTracks = treeOut_nonZeroTracks->GetEntries();
	entries_final_ZeroTracks = treeOut_ZeroTracks->GetEntries();

	cout<<"Outgoing Entries = "<<entries_final<<endl;
	cout<<"Outgoing Entries with nonZeroTracks = "<<entries_final_nonZeroTracks<<endl;
	cout<<"Outgoing Entries with ZeroTracks = "<<entries_final_ZeroTracks<<endl;

	if(!isData)//Calculate inclusive and exclusive efficiencies for MC
	{
		if(entries_init != 0)
		{
			eff_excl = (Float_t)entries_final*100/entries_init;
			eff_excl_err = sqrt(eff_excl*(100.0-eff_excl)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<type<<" CutOutKs cuts made with exclusive efficiency = "<<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genFile>>entries_gen;//NEEDS TO BE TESTED.
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl = (Float_t)entries_final*100/entries_gen;
				eff_incl_err = sqrt(eff_incl*(100.0-eff_incl)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<type<<" CutOutKs cuts made with inclusive efficiency = "<<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;

			genFile.close();
		}
	}

	fileOut->Write();
	fileOut->Close();
	fileOut_nonZeroTracks->Write();
	fileOut_nonZeroTracks->Close();
	fileOut_ZeroTracks->Write();
	fileOut_ZeroTracks->Close();
	fileIn->Close();

	sw.Stop();
	cout << "==> CutOutKs is done! Death to Ks0!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
}
