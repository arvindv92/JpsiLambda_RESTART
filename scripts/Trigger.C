/********************************
   Author : Aravindhan V.
         The purpose of this script is to apply trigger cuts on data/MC coming out of DaVinci.
 *********************************/
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "CollateFiles.h"
#include <iostream>
#include <fstream>

using namespace std;
/* DOCUMENTATION
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   testing = true to run only over a subset of data
   loose = true to run over data from loose stripping line. Only LL for loose stripping line
 */
void Trigger(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0, Bool_t testing = false, Bool_t loose = true)
{
	gROOT->ProcessLine(".L CollateFiles.C++");
	TStopwatch sw;
	sw.Start();

	Int_t entries_init = 0, entries_final = 0, entries_gen = 0;
	Float_t eff_excl = 0., eff_excl_err = 0., eff_incl = 0., eff_incl_err = 0.;
	Bool_t hlt1DiMuonHighMass = false, hlt1TrackMuon = false, hlt1TrackAllL0 = false, hlt2DiMuonDetached = false;
	Bool_t logFlag = false; //This should be 0 only while testing.
	Bool_t collateFlag = true; //If you don't want to re-collate MC, set this to zero. For example, if only the trigger condition changes.
	TFile *fileOut = nullptr, *fileIn = nullptr;
	TTree *treeIn = nullptr, *treeIn_gen = nullptr, *treeOut = nullptr, *myTree = nullptr;
	const char *triggerCut = "(Lb_Hlt1DiMuonHighMassDecision_TOS==1||Lb_Hlt1TrackMuonDecision_TOS==1||Lb_Hlt1TrackAllL0Decision_TOS==1)&&(Lb_Hlt2DiMuonDetachedJPsiDecision_TOS==1)";
	ofstream genFile;//The number of generated MC entries in every MC case will be written out to this file, so that it can be accessed for calculating exclusive efficiencies later

	//gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	// Set up input, output, logging
	if(isData)
	{
		TH1D *lumiHist = nullptr, *lumiErrHist = nullptr;
		Double_t lumi = 0., lumiErr = 0., lumiMean = 0., lumiErrMean = 0.;
		Int_t lumiEntries = 0, lumiErrEntries = 0;
		TChain *h1 = new TChain("Lb2JpsiLTree/MyTuple");
		TChain *h2 = new TChain("GetIntegratedLuminosity/LumiTuple");

		if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/trigger_log.txt",run).Data());

		CollateFiles(run, isData, mcType, &h1, &h2, testing, loose);//CollateFiles will cd to the massdump folder and then back to JpsiLambda_RESTART
		fileOut = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run),"RECREATE");

		// if(run == 1)
		// {
		//      if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run1/trigger_log.txt");
		//
		//      CollateFiles(run, isData, mcType, &h1, &h2, testing, loose);//CollateFiles will cd to the massdump folder and then back to JpsiLambda_RESTART
		//      fileOut = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_triggered.root","RECREATE");
		// }
		// else if(run == 2)
		// {
		//      if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run2/trigger_log.txt");
		//
		//      CollateFiles(run, isData, mcType, &h1, &h2, testing, loose);
		//      fileOut = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_triggered.root","RECREATE");
		// }

		h2->Draw("IntegratedLuminosity>>lumiHist","","goff");
		h2->Draw("IntegratedLuminosityErr>>lumiErrHist","","goff");

		lumiHist    = (TH1D*)gDirectory->Get("lumiHist");
		lumiErrHist = (TH1D*)gDirectory->Get("lumiErrHist");

		lumiMean       = lumiHist->GetMean();
		lumiEntries    = lumiHist->GetEntries();
		lumiErrMean    = lumiErrHist->GetMean();
		lumiErrEntries = lumiErrHist->GetEntries();

		lumi = lumiMean*lumiEntries;
		lumiErr = lumiErrMean*lumiErrEntries;

		cout<<"Processing "<<lumi<<" +/- "<<lumiErr<< "luminosity"<<endl;

		myTree = (TTree*)h1;
	}//end Data block

	else //MC
	{
		if(mcType == 1)//Jpsi Lambda
		{
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiLambda/run%d/trigger_log.txt",run).Data());
			genFile.open(TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run).Data());

			cout<<"PROCESSING MC for Run "<<run<<" Jpsi Lambda"<<endl;

			if(collateFlag) CollateFiles(run,isData,mcType);

			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda.root",run));
			treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");

			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run),"RECREATE");
			// if(run == 1)
			// {
			//      if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run1/trigger_log.txt");
			//      genFile.open("logs/mc/JpsiLambda/run1/gen_log.txt");
			//
			//      cout<<"PROCESSING MC for Run 1 Jpsi Lambda"<<endl;
			//
			//      if(collateFlag) CollateFiles(run,isData,mcType);
			//
			//      fileIn = TFile::Open("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda.root");
			//      treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			//      treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
			//
			//      fileOut = new TFile("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_triggered.root","RECREATE");
			// }
			// if(run == 2)
			// {
			//      if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run2/trigger_log.txt");
			//      genFile.open("logs/mc/JpsiLambda/run2/gen_log.txt");
			//
			//      cout<<"PROCESSING MC for Run 2 Jpsi Lambda"<<endl;
			//
			//      if(collateFlag) CollateFiles(run,isData,mcType);
			//
			//      fileIn = TFile::Open("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda.root");
			//      treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			//      treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
			//
			//      fileOut = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_triggered.root","RECREATE");
			// }
		}
		if(mcType == 2)//Jpsi Sigma
		{
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiSigma/run%d/trigger_log.txt",run).Data());
			genFile.open(TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run).Data());

			cout<<"PROCESSING MC for Run "<<run<<" Jpsi Lambda"<<endl;

			if(collateFlag) CollateFiles(run,isData,mcType);

			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma.root",run));
			treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");

			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_triggered.root",run),"RECREATE");
			// if(run == 1)
			// {
			//      if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run1/trigger_log.txt");
			//      genFile.open("logs/mc/JpsiSigma/run1/gen_log.txt");
			//
			//      cout<<"PROCESSING MC for Run 1 JpsiSigma"<<endl;
			//
			//      if(collateFlag) CollateFiles(run,isData,mcType);
			//
			//      fileIn = TFile::Open("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma.root");
			//      treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			//      treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
			//
			//      fileOut = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_triggered.root","RECREATE");
			// }
			// if(run == 2)
			// {
			//      if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run2/trigger_log.txt");
			//      genFile.open("logs/mc/JpsiSigma/run2/gen_log.txt");
			//
			//      cout<<"PROCESSING MC for Run 2 JpsiSigma"<<endl;
			//
			//      if(collateFlag) CollateFiles(run,isData,mcType);
			//
			//      fileIn = TFile::Open("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma.root");
			//      treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			//      treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
			//
			//      fileOut = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_triggered.root","RECREATE");
			// }
		}
		if(mcType == 3)//Jpsi Xi
		{
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiXi/run%d/trigger_log.txt",run).Data());
			genFile.open(TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run).Data());

			cout<<"PROCESSING MC for Run "<<run<<" Jpsi Lambda"<<endl;

			if(collateFlag) CollateFiles(run,isData,mcType);

			fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi.root",run));
			treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");

			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_triggered.root",run),"RECREATE");
			// if(run == 1)
			// {
			//      if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run1/trigger_log.txt");
			//      genFile.open("logs/mc/JpsiXi/run1/gen_log.txt");
			//
			//      cout<<"PROCESSING MC for Run 1 JpsiXi"<<endl;
			//
			//      if(collateFlag) CollateFiles(run,isData,mcType);
			//
			//      fileIn = TFile::Open("rootFiles/mcFiles/JpsiXi/run1/jpsixi.root");
			//      treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			//      treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
			//
			//      fileOut = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_triggered.root","RECREATE");
			// }
			// if(run == 2)
			// {
			//      if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run2/trigger_log.txt");
			//      genFile.open("logs/mc/JpsiXi/run2/gen_log.txt");
			//
			//      cout<<"PROCESSING MC for Run 2 JpsiXi"<<endl;
			//
			//      if(collateFlag) CollateFiles(run,isData,mcType);
			//
			//      fileIn = TFile::Open("rootFiles/mcFiles/JpsiXi/run2/jpsixi.root");
			//      treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
			//      treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
			//
			//      fileOut = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_triggered.root","RECREATE");
			// }
		}
		entries_gen = treeIn_gen->GetEntries();
		genFile<<entries_gen<<endl;
		genFile.close();
		entries_init = treeIn->GetEntries();
		myTree = treeIn;
	}//end MC block
	 //end setup of input, output, logging

	cout<<"******************************************"<<endl;
	cout<<"==> Starting Trigger: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	cout<<"******************************************"<<endl;
	if(!isData) cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = myTree->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeOut = (TTree*)myTree->CloneTree(0);

	myTree->SetBranchAddress("Lb_Hlt1DiMuonHighMassDecision_TOS",&hlt1DiMuonHighMass);
	myTree->SetBranchAddress("Lb_Hlt1TrackMuonDecision_TOS",&hlt1TrackMuon);
	myTree->SetBranchAddress("Lb_Hlt1TrackAllL0Decision_TOS",&hlt1TrackAllL0);
	myTree->SetBranchAddress("Lb_Hlt2DiMuonDetachedJPsiDecision_TOS",&hlt2DiMuonDetached);

	cout<<"I am making the following trigger cuts. Sit tight"<<endl;
	cout<<triggerCut<<endl;

	for(Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0) cout<<i<<endl;

		myTree->GetEntry(i);
		if(hlt2DiMuonDetached)
		{
			if(hlt1DiMuonHighMass||hlt1TrackMuon||hlt1TrackAllL0)
			{
				treeOut->Fill();
			}
		}
	}

	entries_final = treeOut->GetEntries();
	cout<<"Outgoing entries = "<<entries_final<<endl;

	if(!isData)//Calculate exclusive and inclusive efficiencies for MC
	{
		if(entries_init != 0)
		{
			eff_excl = (Float_t)entries_final*100/entries_init;
			eff_excl_err = sqrt( eff_excl*(100.0-eff_excl)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<"Trigger cut made with exclusive efficiency = "<<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(entries_gen != 0)
		{
			eff_incl = (Float_t)entries_final*100/entries_gen;
			eff_incl_err = sqrt( eff_incl*(100.0-eff_incl)/entries_gen);
		}
		cout<<"******************************************"<<endl;
		cout<<"Trigger cut made with inclusive efficiency = "<<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;
	}

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	sw.Stop();
	cout << "==> Trigger is done! Huzzah!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
	//if(logFlag) gSystem->Exec("cat trigger_log.txt");
}
