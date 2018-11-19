/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply loose sanity cuts to triggered data/MC, and then separate out into LL and DD files.
 *********************************/
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
using namespace std;

void Sanity(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
 */
/*
   NB: BKGCAT applies to MC, but not data. Read  Steve's analysis and think about exactly what TM condition is to be applied
   Lb_DTF_CTAU_L -> Lb_TAU
   DTFchi2 cut is now 0-50
   PID cuts are applied to reject junk events with PID's of 0 or -1000 for Run1 and reject -1000 for Run2
 */
{
	TStopwatch sw;
	sw.Start();

	TFile *fileIn = nullptr, *fileOut_LL = nullptr, *fileOut_DD = nullptr;
	TTree *treeIn = nullptr, *treeOut_LL = nullptr, *treeOut_DD = nullptr;
	Int_t entries_init = 0, entries_final_LL = 0, entries_final_DD, entries_gen = 0, Lb_BKGCAT = 0, p_TRACK_Type = 0;
	Float_t eff_excl_LL = 0.0, eff_excl_LL_err = 0.0,eff_incl_LL = 0.0,eff_incl_LL_err = 0.0;
	Float_t eff_excl_DD = 0.0, eff_excl_DD_err = 0.0,eff_incl_DD = 0.0,eff_incl_DD_err = 0.0;
	Float_t Lb_ConsLb_chi2 = 0.0, Lb_ConsLb_nDOF = 0.0;
	Double_t Lb_TAU = 0.0, Lb_ETA = 0.0;
	Double_t pi_PIDp = 0.0, pi_PIDK = 0.0, p_PIDp = 0.0, p_PIDK = 0.0, L_TAU = 0.0;
	Bool_t logFlag = false, genFlag = false;
	const char* logFileName = "sanity_log.txt";
	TCut dtfCut, lifetimeCut, pidCut, etaCut, tmCut;

	fstream genFile;//contains no of generated candidates in MC.

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.
	cout<<"WD = "<<gSystem->pwd()<<endl;

	if(!isData)  // MC
	{
		switch(mcType)
		{
		case 1:                 //JpsiLambda
			if(logFlag) gROOT->ProcessLine((TString::Format(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName)).Data());
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data()))
			{
				genFile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
			treeIn = (TTree*)fileIn->Get("MyTuple");
			fileOut_LL = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_sanity_LL.root",run),"RECREATE");
			treeOut_LL = (TTree*)treeIn->CloneTree(0);
			fileOut_DD = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_sanity_DD.root",run),"RECREATE");
			treeOut_DD = (TTree*)treeIn->CloneTree(0);
			break;

		case 2: //JpsiSigma
			if(logFlag) gROOT->ProcessLine((TString::Format(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName)).Data());
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data()))
			{
				genFile.open((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_triggered.root",run));
			treeIn = (TTree*)fileIn->Get("MyTuple");
			fileOut_LL = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_sanity_LL.root",run),"RECREATE");
			treeOut_LL = (TTree*)treeIn->CloneTree(0);
			fileOut_DD = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_sanity_DD.root",run),"RECREATE");
			treeOut_DD = (TTree*)treeIn->CloneTree(0);
			break;

		case 3: //JpsiXi
			if(logFlag) gROOT->ProcessLine((TString::Format(".> logs/mc/JpsiXi/run%d/%s",run,logFileName)).Data());
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data()))
			{
				genFile.open((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_triggered.root",run));
			treeIn = (TTree*)fileIn->Get("MyTuple");
			fileOut_LL = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_sanity_LL.root",run),"RECREATE");
			treeOut_LL = (TTree*)treeIn->CloneTree(0);
			fileOut_DD = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_sanity_DD.root",run),"RECREATE");
			treeOut_DD = (TTree*)treeIn->CloneTree(0);
			break;
		}
		treeIn->SetBranchAddress("Lb_BKGCAT",&Lb_BKGCAT);
	} // end MC block
	else //Data
	{
		if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%ss",run,logFileName).Data());
		fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		fileOut_LL = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL.root",run),"RECREATE");
		treeOut_LL = (TTree*)treeIn->CloneTree(0);
		fileOut_DD = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_DD.root",run),"RECREATE");
		treeOut_DD = (TTree*)treeIn->CloneTree(0);
	}//end Data block

	cout<<"Starting Sanity"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"LL Output file = "<<fileOut_LL->GetName()<<endl;
	cout<<"DD Output file = "<<fileOut_DD->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("Lb_TAU",&Lb_TAU);
	treeIn->SetBranchAddress("L_TAU",&L_TAU);
	treeIn->SetBranchAddress("Lb_ConsLb_chi2",&Lb_ConsLb_chi2);
	treeIn->SetBranchAddress("Lb_ConsLb_nDOF",&Lb_ConsLb_nDOF);
	treeIn->SetBranchAddress("Lb_ETA",&Lb_ETA);
	treeIn->SetBranchAddress("pi_PIDp",&pi_PIDp);
	treeIn->SetBranchAddress("pi_PIDK",&pi_PIDK);
	treeIn->SetBranchAddress("p_PIDp",&p_PIDp);
	treeIn->SetBranchAddress("p_PIDK",&p_PIDK);
	treeIn->SetBranchAddress("p_TRACK_Type",&p_TRACK_Type);

	lifetimeCut = "(Lb_TAU > 0)&&(L_TAU > 0)";
	dtfCut = "(Lb_ConsLb_chi2/Lb_ConsLb_nDOF > 0 && Lb_ConsLb_chi2/Lb_ConsLb_nDOF < 50)";
	etaCut = "(Lb_ETA > 2 && Lb_ETA < 6)";

	if(run == 1)
	{
		pidCut = "(pi_PIDp != 0 && pi_PIDp != -1000 && pi_PIDK != 0 && p_PIDp !=0 && p_PIDp!= -1000 && p_PIDK !=0)";//This also ensures that pi_PIDK != -1000 && p_PIDK != -1000
	}
	else if(run == 2)
	{
		pidCut = "(pi_PIDp != -1000 && p_PIDp!= -1000)";
	}
	if(isData)
	{
		if(mcType == 1 || mcType == 2)
		{
			tmCut  = "(Lb_BKGCAT == 0 || Lb_BKGCAT == 50)";
		}
		if(mcType == 3)
		{
			tmCut = "(Lb_BKGCAT == 40)";
		}
	}

	cout<<"I am making the following sanity cuts, and then separating in LL and DD files. Sit tight"<<endl;

	lifetimeCut.Print();
	dtfCut.Print();
	pidCut.Print();
	etaCut.Print();
	if(!isData) tmCut.Print();

	for(Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn->GetEntry(i);
		if(Lb_TAU > 0 && L_TAU > 0) //require lifetimes to be positive, universal
		{
			if(Lb_ETA > 2 && Lb_ETA < 6) //acceptance cut, universal
			{
				if((Lb_ConsLb_chi2/Lb_ConsLb_nDOF) > 0 && (Lb_ConsLb_chi2/Lb_ConsLb_nDOF) < 50)// DTF chi2 cut, universal
				{
					if((run == 1 && pi_PIDp != 0 && pi_PIDp != -1000 && pi_PIDK != 0 && p_PIDp !=0 && p_PIDp!= -1000 && p_PIDK !=0)|| (run == 2 && p_PIDp != -1000 && pi_PIDp != -1000))
					{
						if(isData == 1 || (isData == 0 && mcType == 3 && Lb_BKGCAT == 40) || (isData == 0 && (mcType == 1 || mcType == 2) && (Lb_BKGCAT == 0 || Lb_BKGCAT == 50)))//TM, needs more thought.
						{
							if(p_TRACK_Type == 3)
							{
								treeOut_LL->Fill();
							}
							else
							{
								treeOut_DD->Fill();
							}
						}
					}
				}
			}
		}
	}

	entries_final_LL = treeOut_LL->GetEntries();
	entries_final_DD = treeOut_DD->GetEntries();
	cout<<"Outgoing LL entries = "<<entries_final_LL<<endl;
	cout<<"Outgoing DD entries = "<<entries_final_DD<<endl;

	if(isData==0)//Efficiency calculation for MC.
	{
		if(entries_init != 0)//Exclusive efficiency calculation
		{
			eff_excl_LL = (Float_t)entries_final_LL*100/entries_init;
			eff_excl_LL_err = sqrt( eff_excl_LL*(100.0-eff_excl_LL)/entries_init);

			eff_excl_DD = (Float_t)entries_final_DD*100/entries_init;
			eff_excl_DD_err = sqrt( eff_excl_DD*(100.0-eff_excl_DD)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<"LL Sanity cuts made with exclusive efficiency = "<<eff_excl_LL<<"% +/- " <<eff_excl_LL_err<<" %"<<endl;
		cout<<"DD Sanity cuts made with exclusive efficiency = "<<eff_excl_DD<<"% +/- " <<eff_excl_DD_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genFile>>entries_gen;//Inclusive efficiency calculation.
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl_LL = (Float_t)entries_final_LL*100/entries_gen;
				eff_incl_LL_err = sqrt(eff_incl_LL*(100.0-eff_incl_LL)/entries_gen);

				eff_incl_DD = (Float_t)entries_final_DD*100/entries_gen;
				eff_incl_DD_err = sqrt(eff_incl_DD*(100.0-eff_incl_DD)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<"LL Sanity cuts made with inclusive efficiency = "<<eff_incl_LL<<"% +/- " <<eff_incl_LL_err<<" %"<<endl;
			cout<<"DD Sanity cuts made with inclusive efficiency = "<<eff_incl_DD<<"% +/- " <<eff_incl_DD_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
	}

	fileOut_LL->Write();
	fileOut_LL->Close();

	fileOut_DD->Write();
	fileOut_DD->Close();

	sw.Stop();
	cout << "--- End of Sanity: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
	//if(logFlag) gSystem->Exec("cat sanity_log.txt");
}
