#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

using namespace std;

void cutoutks(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0, Int_t trackType = 3)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{
	TFile *filein(0), *fileout(0);
	TTree *treein(0), *treeout(0);

	TCut L_dmcut, pidcut,masswindow, totalcut;

	Int_t entries_init = 0, entries_final = 0, entries_gen = 0;
	Float_t eff_excl = 0.0, eff_excl_err = 0.0, eff_incl = 0.0, eff_incl_err = 0.0;
	Double_t L_dm = 0.0, p_PIDp = 0.0, Lb_DTF_L_WMpipi_JpsiConstr = 0.0;
	Bool_t logFlag = 0, genFlag = 0;
	//TEST CHANGE TO TEST GITHUB
	fstream genfile;
	masswindow = "Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520";

	const char* type = (trackType == 3) ? ("LL") : ("DD");
	const char* logFileName;

	if(trackType == 3)
	{
		cout<<"Processing LL"<<endl;
		logFileName = "cutoutks_LL_log.txt";
	}
	if(trackType == 5)
	{
		cout<<"Processing DD"<<endl;
		logFileName = "cutoutks_DD_log.txt";
	}

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

	switch(isData)//Setup logging, input and output files.
	{
	case 0: // MC
		switch(mcType)
		{
		case 1: //JpsiLambda
			switch(run)
			{
			case 1:
				if(logFlag == 1)
				{
					gROOT->ProcessLine(strcat(".> logs/mc/JpsiLambda/run1/",logFileName));//What is ProcessLineFast()?
				}
				if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run1/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiLambda/run1/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_sanity_%s.root",type),"READ");
				fileout = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_cutoutks_%s.root",type),"RECREATE");
				break;

			case 2:
				if(logFlag == 1)
				{
					gROOT->ProcessLine(strcat(".> logs/mc/JpsiLambda/run2/cutoutks_LL_log.txt",logFileName));
				}
				if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run2/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiLambda/run2/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run2/jpsilambd_sanity_%s.root",type),"READ");
				fileout = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_cutoutks_%s.root",type),"RECREATE");
				break;
			}
			break;

		case 2: //JpsiSigma
			switch(run)
			{
			case 1:
				if(logFlag == 1)
				{
					gROOT->ProcessLine(strcat(".> logs/mc/JpsiSigma/run1/cutoutks_LL_log.txt",logFileName));
				}
				if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run1/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiSigma/run1/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_%s.root",type),"READ");
				fileout = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_cutoutks_%s.root",type),"RECREATE");
				break;

			case 2:
				if(logFlag == 1)
				{
					gROOT->ProcessLine(strcat(".> logs/mc/JpsiSigma/run2/cutoutks_LL_log.txt",logFileName));
				}
				if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run2/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiSigma/run2/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_%s.root",type),"READ");
				fileout = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_cutoutks_%s.root",type),"RECREATE");
				break;
			}
			break;

		case 3: //JpsiXi
			switch(run)
			{
			case 1:
				if(logFlag == 1)
				{
					gROOT->ProcessLine(strcat(".> logs/mc/JpsiXi/run1/cutoutks_LL_log.txt",logFileName));
				}
				if(!gSystem->AccessPathName("logs/mc/JpsiXi/run1/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiXi/run1/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_%s.root",type),"READ");
				fileout = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run1/jpsixi_cutoutks_%s.root",type),"RECREATE");
				break;

			case 2:
				if(logFlag == 1)
				{
					gROOT->ProcessLine(strcat(".> logs/mc/JpsiXi/run2/cutoutks_LL_log.txt",logFileName));
				}
				if(!gSystem->AccessPathName("logs/mc/JpsiXi/run2/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiXi/run2/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_%s.root",type),"READ");
				fileout = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run2/jpsixi_cutoutks_%s.root",type),"RECREATE");
				break;
			}
			break;
		}
		treein = (TTree*)filein->Get("MyTuple");
		break;

	case 1: // Data
		switch(run)
		{
		case 1:
			if(logFlag == 1)
			{
				gROOT->ProcessLine(strcat(".> logs/data/JpsiLambda/run1/cutoutks_LL_log.txt",logFileName));
			}
			filein = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_%s.root",type),"READ");
			fileout = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_%s.root",type),"RECREATE");
			break;

		case 2:
			if(logFlag == 1)
			{
				gROOT->ProcessLine(strcat(".> logs/data/JpsiLambda/run2/cutoutks_LL_log.txt",logFileName));
			}
			filein = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_%s.root",type),"READ");
			fileout = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_%s.root",type),"RECREATE");
			break;
		}
		treein = (TTree*)filein->Get("MyTuple");
		break;
	}
	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<filein->GetName()<<endl;
	cout<<"Output file = "<<fileout->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treein->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treein->SetBranchAddress("L_dm",&L_dm);
	treein->SetBranchAddress("p_PIDp",&p_PIDp);
	treein->SetBranchAddress("Lb_DTF_L_WMpipi_JpsiConstr",&Lb_DTF_L_WMpipi_JpsiConstr);

	treeout = (TTree*)treein->CloneTree(0);

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

	totalcut = (L_dmcut && pidcut && masswindow) || (!masswindow);

	cout<<"I am making the following cutoutks cuts. Sit tight"<<endl;
	totalcut.Print();

	for (Int_t i = 0; i < entries_init; i++)
	{
		if(i % 100000 == 0) cout<<i<<endl;
		treein->GetEntry(i);
		if(trackType == 3)
		{
			if((L_dm < 7.5 && p_PIDp > 10 && Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520) || !(Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520))
			{
				treeout->Fill();
			}
		}
		else if(trackType == 5)
		{
			if((L_dm < 10 && p_PIDp > 15 && Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520) || !(Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520))
			{
				treeout->Fill();
			}
		}
	}

	entries_final = treeout->GetEntries();
	cout<<"Outgoing Entries = "<<entries_final<<endl;

	if(isData==0)
	{
		if(entries_init != 0)
		{
			eff_excl = (Float_t)entries_final*100/entries_init;
			eff_excl_err = sqrt( eff_excl*(100.0-eff_excl)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<type<<" Sanity cuts made with exclusive efficiency = "<<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genfile>>entries_gen;//NEEDS TO BE TESTED.
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl = (Float_t)entries_final*100/entries_gen;
				eff_incl_err = sqrt(eff_incl*(100.0-eff_incl)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<type<<" Sanity cuts made with inclusive efficiency = "<<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
	}

	if(logFlag) gROOT->ProcessLine(".>");

	fileout->Write();
	fileout->Close();
	filein->Close();
}
