/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply cuts to triggered Jpsi Xi data/MC
   Split by year is only done for data, to allow parallel processing of different years.
   For MC, all years for a given run are processed together.
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

void Cuts_JpsiXi(Int_t run = 1, Int_t year = 2011, Bool_t isData = true, Bool_t logFlag = false)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
 */
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	//Set up logging
	if(isData && logFlag)
		gROOT->ProcessLine(Form(".> logs/data/JpsiXi/run%d/cuts_JpsiXi_%d_log.txt",run,year));
	else if(!isData && logFlag)
		gROOT->ProcessLine((Form(".> logs/mc/JpsiXi/run%d/cuts_JpsiXi_%d_log.txt",run,year)));

	TStopwatch sw;
	sw.Start();

	cout<<"******************************************"<<endl;
	cout<<"==> Starting Cuts_JpsiXi: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	Int_t entries_init  = 0, entries_final_LL   = 0, entries_gen = 0;
	Int_t Xib_BKGCAT    = 0, p_TRACK_Type       = 0;
	Int_t Xib_TRUEID    = 0, Xib_ENDVERTEX_NDOF = 0;
	Float_t eff_excl_LL = 0., eff_excl_LL_err   = 0.;
	Float_t eff_incl_LL = 0., eff_incl_LL_err   = 0.;

	Double_t L_TAU              = 0., L_M    = 0., L_FD_OWNPV = 0.;
	Double_t L_ENDVERTEX_ZERR   = 0., p_PT   = 0., pi_PT      = 0.;
	Double_t Xib_ENDVERTEX_CHI2 = 0., Xib_PT = 0.;

	Double_t Xib_IPCHI2_OWNPV    = 0., Xib_TAU                 = 0.;
	Double_t BachPi_IPCHI2_OWNPV = 0., Xi_M                    = 0.;
	Double_t Xi_DIRA_ORIVX       = 0., Xi_TAU                  = 0.;
	Double_t Xib_ETA             = 0., Xib_DTF_M_JpsiXiLConstr = 0.;

	Bool_t genFlag = false;

	TFile *fileIn = nullptr, *fileOut_LL = nullptr;
	TTree *treeIn = nullptr, *treeOut_LL = nullptr;

	fstream genFile;//contains no of generated candidates in MC.

	// Set up input, output
	if(!isData) // MC
	{
		if(!gSystem->AccessPathName((Form("logs/mc/JpsiXi/run%d/gen_log.txt",run))))
		{
			genFile.open((Form("logs/mc/JpsiXi/run%d/gen_log.txt",run)));
			genFlag = true;
		}
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_triggered.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		fileOut_LL = new TFile(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cut_LL.root",run),"RECREATE");
		treeOut_LL = (TTree*)treeIn->CloneTree(0);

		treeIn->SetBranchAddress("Xib_BKGCAT",&Xib_BKGCAT);
		treeIn->SetBranchAddress("Xib_TRUEID",&Xib_TRUEID);
	} // end MC block
	else //Data
	{
		fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiXi/run%d/jpsixi_triggered_%d.root",run,year));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		fileOut_LL = new TFile(Form("rootFiles/dataFiles/JpsiXi/run%d/jpsixi_cut_LL_%d.root",run,year),"RECREATE");
		treeOut_LL = (TTree*)treeIn->CloneTree(0);
	}//end Data block
	 //end setup of input, output

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"LL Output file = "<<fileOut_LL->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("L_TAU",&L_TAU);
	treeIn->SetBranchAddress("L_M",&L_M);
	treeIn->SetBranchAddress("L_FD_OWNPV",&L_FD_OWNPV);
	treeIn->SetBranchAddress("L_ENDVERTEX_ZERR",&L_ENDVERTEX_ZERR);
	treeIn->SetBranchAddress("p_PT",&p_PT);
	treeIn->SetBranchAddress("pi_PT",&pi_PT);
	treeIn->SetBranchAddress("Xib_ENDVERTEX_CHI2",&Xib_ENDVERTEX_CHI2);
	treeIn->SetBranchAddress("Xib_ENDVERTEX_NDOF",&Xib_ENDVERTEX_NDOF);
	treeIn->SetBranchAddress("Xib_IPCHI2_OWNPV",&Xib_IPCHI2_OWNPV);
	treeIn->SetBranchAddress("Xib_TAU",&Xib_TAU);
	treeIn->SetBranchAddress("BachPi_IPCHI2_OWNPV",&BachPi_IPCHI2_OWNPV);
	treeIn->SetBranchAddress("Xi_M",&Xi_M);
	treeIn->SetBranchAddress("Xi_DIRA_ORIVX",&Xi_DIRA_ORIVX);
	treeIn->SetBranchAddress("Xi_TAU",&Xi_TAU);
	treeIn->SetBranchAddress("Xib_ETA",&Xib_ETA);
	treeIn->SetBranchAddress("Xib_DTF_M_JpsiXiLConstr",&Xib_DTF_M_JpsiXiLConstr);
	treeIn->SetBranchAddress("Xib_PT",&Xib_PT);

	for(Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn->GetEntry(i);
		if(Xib_TAU > 0.0002 && Xi_TAU > 0 && L_TAU > 0) //require lifetimes to be positive, universal
		{
			if(abs(L_M-1115.7) < 8 && abs(Xi_M - L_M -206.2) < 10)
			{
				if((L_FD_OWNPV/L_ENDVERTEX_ZERR) > 5 && BachPi_IPCHI2_OWNPV > 16 && Xib_IPCHI2_OWNPV < 12)
				{
					if((Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF < 10) && Xi_DIRA_ORIVX > 0.995)
					{
						if(p_PT > 500 && pi_PT > 100)
						{
							if(Xib_ETA > 2 && Xib_ETA < 6 && Xib_PT < 20000) //acceptance cut, universal
							{
								if(isData == 1 || (isData == 0 && (abs(Xib_TRUEID) == 5132 || (Xib_BKGCAT == 60 && (abs(Xib_DTF_M_JpsiXiLConstr - 5795) < 35))) ))//TM, needs more thought.
								{
									treeOut_LL->Fill();
								}
							}
						}
					}
				}
			}
		}
	}

	entries_final_LL = treeOut_LL->GetEntries();
	cout<<"Outgoing LL entries = "<<entries_final_LL<<endl;

	if(isData==0)//Efficiency calculation for MC.
	{
		if(entries_init != 0)//Exclusive efficiency calculation
		{
			eff_excl_LL     = (Float_t)entries_final_LL*100/entries_init;
			eff_excl_LL_err = sqrt( eff_excl_LL*(100.0-eff_excl_LL)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<"LL cuts made with exclusive efficiency = "
		    <<eff_excl_LL<<"% +/- " <<eff_excl_LL_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genFile>>entries_gen;//Inclusive efficiency calculation.
			cout<<"Original generated number = "<<entries_gen<<endl;

			if(entries_gen != 0)
			{
				eff_incl_LL     = (Float_t)entries_final_LL*100/entries_gen;
				eff_incl_LL_err = sqrt(eff_incl_LL*(100.0-eff_incl_LL)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<"LL cuts made with inclusive efficiency = "
			    <<eff_incl_LL<<"% +/- " <<eff_incl_LL_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
	}

	fileOut_LL->Write();
	fileOut_LL->Close();

	sw.Stop();
	cout << "==> Cuts_JpsiXi is done! Mazel Tov!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
}
