/********************************
   Author : Aravindhan V.
 *********************************/
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>
using namespace std;

void sanity(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
 */

/*
   NB: BKGCAT applied to MC, but not data. Read  Steve's analysis and think about exactly what TM condition is to be applied
   Lb_DTF_CTAU_L -> Lb_TAU
   DTFchi2 cut is now 0-50
   PID cuts are applied on Run1 but not Run2 right now. See if some similar PID cuts can be applied on Run2.
 */
{
	//if(logFlag) gROOT->ProcessLine(".> sanity_log.txt");//This should be 0 only while testing.

	TFile *filein(0), *fileout_LL(0), *fileout_DD(0);
	TTree *treein(0), *treeout_LL(0), *treeout_DD(0);

	Int_t entries_init = 0, entries_final_LL = 0, entries_final_DD, entries_gen = 0, Lb_BKGCAT = 0, p_TRACK_Type = 0;
	Float_t eff_excl_LL = 0.0, eff_excl_LL_err = 0.0,eff_incl_LL = 0.0,eff_incl_LL_err = 0.0;
	Float_t eff_excl_DD = 0.0, eff_excl_DD_err = 0.0,eff_incl_DD = 0.0,eff_incl_DD_err = 0.0;
	Double_t Lb_TAU = 0.0, Lb_ConsLb_chi2 = 0.0, Lb_ConsLb_nDOF = 0.0, Lb_ETA = 0.0;
	Double_t pi_PIDp = 0.0, pi_PIDK = 0.0, p_PIDp = 0.0, p_PIDK = 0.0, MyL_TAU = 0.0;
	Bool_t logFlag = false, genFlag = false;
	const char* logFileName = "sanity_log.txt";
	fstream genfile;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

	switch(isData)
	{
	case 0: // MC
		switch(mcType)
		{
		case 1: //JpsiLambda
			if(logFlag) gROOT->ProcessLine((TString::Format(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName)).Data());
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data()))
			{
				genfile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
			fileout_LL = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_sanity_LL.root",run),"RECREATE");
			fileout_DD = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_sanity_DD.root",run),"RECREATE");
			/*switch(run)
			   {
			   case 1:
			        if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run1/sanity_log.txt");
			        if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run1/gen_log.txt"))
			        {
			                genfile.open("logs/mc/JpsiLambda/run1/gen_log.txt");
			                genFlag = true;
			        }
			        filein = TFile::Open("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_triggered.root");
			        fileout_LL = new TFile("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_sanity_LL.root","RECREATE");
			        fileout_DD = new TFile("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_sanity_DD.root","RECREATE");
			        break;

			   case 2:
			        if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run2/sanity_log.txt");
			        if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run2/gen_log.txt"))
			        {
			                genfile.open("logs/mc/JpsiLambda/run2/gen_log.txt");
			                genFlag = true;
			        }
			        filein = TFile::Open("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_triggered.root");
			        fileout_LL = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root","RECREATE");
			        fileout_DD = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity_DD.root","RECREATE");
			        break;
			   }
			 */
			break;

		case 2: //JpsiSigma
			if(logFlag) gROOT->ProcessLine((TString::Format(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName)).Data());
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data()))
			{
				genfile.open((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_triggered.root",run));
			fileout_LL = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_sanity_LL.root",run),"RECREATE");
			fileout_DD = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_sanity_DD.root",run),"RECREATE");
			/*	switch(run)
			        {
			        case 1:
			                if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run1/sanity_log.txt");
			                if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run1/gen_log.txt"))
			                {
			                        genfile.open("logs/mc/JpsiSigma/run1/gen_log.txt");
			                        genFlag = true;
			                }
			                filein = TFile::Open("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_triggered.root");
			                fileout_LL = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_LL.root","RECREATE");
			                fileout_DD = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_DD.root","RECREATE");
			                break;

			        case 2:
			                if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run2/sanity_log.txt");
			                if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run2/gen_log.txt"))
			                {
			                        genfile.open("logs/mc/JpsiSigma/run2/gen_log.txt");
			                        genFlag = true;
			                }
			                filein = TFile::Open("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_triggered.root");
			                fileout_LL = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_LL.root","RECREATE");
			                fileout_DD = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_DD.root","RECREATE");
			                break;
			        }
			 */
			break;

		case 3: //JpsiXi
			if(logFlag) gROOT->ProcessLine((TString::Format(".> logs/mc/JpsiXi/run%d/%s",run,logFileName)).Data());
			if(!gSystem->AccessPathName((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data()))
			{
				genfile.open((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data());
				genFlag = true;
			}
			filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_triggered.root",run));
			fileout_LL = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_sanity_LL.root",run),"RECREATE");
			fileout_DD = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_sanity_DD.root",run),"RECREATE");
			/*	switch(run)
			        {
			        case 1:
			                if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run1/sanity_log.txt");
			                if(!gSystem->AccessPathName("logs/mc/JpsiXi/run1/gen_log.txt"))
			                {
			                        genfile.open("logs/mc/JpsiXi/run1/gen_log.txt");
			                        genFlag = true;
			                }
			                filein = TFile::Open("rootFiles/mcFiles/JpsiXi/run1/jpsixi_triggered.root");
			                fileout_LL = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_LL.root","RECREATE");
			                fileout_DD = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_DD.root","RECREATE");
			                break;

			        case 2:
			                if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run2/sanity_log.txt");
			                if(!gSystem->AccessPathName("logs/mc/JpsiXi/run2/gen_log.txt"))
			                {
			                        genfile.open("logs/mc/JpsiXi/run2/gen_log.txt");
			                        genFlag = true;
			                }
			                filein = TFile::Open("rootFiles/mcFiles/JpsiXi/run2/jpsixi_triggered.root");
			                fileout_LL = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_LL.root","RECREATE");
			                fileout_DD = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_DD.root","RECREATE");
			                break;
			        }
			 */
			break;
		}
		treein = (TTree*)filein->Get("MyTuple");
		treein->SetBranchAddress("Lb_BKGCAT",&Lb_BKGCAT);
		break;

	case 1: // Data
		if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%ss",run,logFileName).Data());
		filein = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
		fileout_LL = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL.root",run),"RECREATE");
		fileout_DD = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_DD.root",run),"RECREATE");
		/*switch(run)
		   {
		   case 1:
		        if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run1/sanity_log.txt");
		        filein = TFile::Open("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_triggered.root");
		        fileout_LL = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_LL.root","RECREATE");
		        fileout_DD = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_DD.root","RECREATE");
		        break;

		   case 2:
		        if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run2/sanity_log.txt");
		        filein = TFile::Open("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_triggered.root");
		        fileout_LL = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root","RECREATE");
		        fileout_DD = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_DD.root","RECREATE");
		        break;
		   }
		 */
		treein = (TTree*)filein->Get("MyTuple");
		break;
	}
	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<filein->GetName()<<endl;
	cout<<"LL Output file = "<<fileout_LL->GetName()<<endl;
	cout<<"DD Output file = "<<fileout_DD->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treein->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treein->SetBranchAddress("Lb_TAU",&Lb_TAU);
	treein->SetBranchAddress("L_TAU",&MyL_TAU);
	treein->SetBranchAddress("Lb_ConsLb_chi2",&Lb_ConsLb_chi2);
	treein->SetBranchAddress("Lb_ConsLb_nDOF",&Lb_ConsLb_nDOF);
	treein->SetBranchAddress("Lb_ETA",&Lb_ETA);
	treein->SetBranchAddress("pi_PIDp",&pi_PIDp);
	treein->SetBranchAddress("pi_PIDK",&pi_PIDK);
	treein->SetBranchAddress("p_PIDp",&p_PIDp);
	treein->SetBranchAddress("p_PIDK",&p_PIDK);
	treein->SetBranchAddress("p_TRACK_Type",&p_TRACK_Type);

	treeout_LL = (TTree*)treein->CloneTree(0);
	treeout_DD = (TTree*)treein->CloneTree(0);

	TCut dtfcut, lifetimecut, pidcut, etacut, tmcut;

	lifetimecut = "(Lb_TAU > 0)&&(L_TAU > 0)";
	dtfcut = "(Lb_ConsLb_chi2/Lb_ConsLb_nDOF > 0 && Lb_ConsLb_chi2/Lb_ConsLb_nDOF < 50)";
	etacut = "(Lb_ETA > 2 && Lb_ETA < 6)";

	if(run == 1)
	{
		pidcut = "(pi_PIDp != 0 && pi_PIDp != -1000 && pi_PIDK != 0 && p_PIDp !=0 && p_PIDp!= -1000 && p_PIDK !=0)";//This also ensures that pi_PIDK != -1000 && p_PIDK != -1000
	}
	else if(run == 2)
	{
		pidcut = "(pi_PIDp != -1000 && p_PIDp!= -1000)";
	}
	if(isData == 0)
	{
		if(mcType == 1 || mcType == 2)
		{
			tmcut  = "(Lb_BKGCAT == 0 || Lb_BKGCAT == 50)";
		}
		if(mcType == 3)
		{
			tmcut = "(Lb_BKGCAT == 40)";
		}
	}

	cout<<"I am making the following sanity cuts, and then separating in LL and DD files. Sit tight"<<endl;

	lifetimecut.Print();
	dtfcut.Print();
	pidcut.Print();
	etacut.Print();
	if(!isData) tmcut.Print();

	for(Int_t i = 0; i < entries_init; i++)
	{
		treein->GetEntry(i);
		if(Lb_TAU > 0 && MyL_TAU > 0) //require lifetimes to be positive, universal
		{
			if(Lb_ETA > 2 && Lb_ETA < 6) //acceptance cut, universal
			{
				if(Lb_ConsLb_chi2/Lb_ConsLb_nDOF > 0 && Lb_ConsLb_chi2/Lb_ConsLb_nDOF < 50)// DTF chi2 cut, universal
				{
					if((run == 1 && pi_PIDp != 0 && pi_PIDp != -1000 && pi_PIDK != 0 && p_PIDp !=0 && p_PIDp!= -1000 && p_PIDK !=0)|| (p_PIDp != -1000 && pi_PIDp != -1000 && run == 2))
					{
						if(isData == 1 || (isData == 0 && mcType == 3 && Lb_BKGCAT == 40) || (isData == 0 && (mcType == 1 || mcType == 2) && (Lb_BKGCAT == 0 || Lb_BKGCAT == 50)))//TM, needs more thought.
						{
							if(p_TRACK_Type == 3)
							{
								treeout_LL->Fill();
							}
							else
							{
								treeout_DD->Fill();
							}

						}
					}
				}
			}
		}
	}

	entries_final_LL = treeout_LL->GetEntries();
	entries_final_DD = treeout_DD->GetEntries();
	cout<<"Outgoing LL entries = "<<entries_final_LL<<endl;
	cout<<"Outgoing DD entries = "<<entries_final_DD<<endl;

	if(isData==0)
	{
		if(entries_init != 0)
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
			genfile>>entries_gen;//NEEDS TO BE TESTED.
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

	fileout_LL->cd();//This writing machanism needs to be tested
	treeout_LL->Write();
	fileout_LL->Close();

	fileout_DD->cd();
	treeout_DD->Write();
	fileout_DD->Close();

	if(logFlag) gROOT->ProcessLine(".>");
	//if(logFlag) gSystem->Exec("cat sanity_log.txt");
}
