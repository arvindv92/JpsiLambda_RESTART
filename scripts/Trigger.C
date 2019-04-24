/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply trigger cuts on data/MC coming out of DaVinci.
 *********************************/
#include "Trigger.h"

void Trigger(Int_t run, Int_t year, Bool_t isData, Int_t mcType, Bool_t testing,
             Bool_t loose, Bool_t logFlag)
/* DOCUMENTATION
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   mcType = 0 when running over data.
   When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   MCType = 4 for Bu_JpsiX, 5 for Bd_JpsiX
   testing = true to run only over a subset of data
   loose = true to run over data from loose stripping line. Only LL for loose stripping line
 */
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	const char *folder = "", *part = "";
	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Bu_JpsiX";
		part = "bu_jpsix";
		break;
	}
	case 5:
	{
		folder = "Bd_JpsiX";
		part = "bd_jpsix";
		break;
	}
	case 6:
	{
		folder = "Lst1405";
		part = "lst1405_BONN";
		break;
	}
	case 7:
	{
		folder = "Lst1520";
		part = "lst1520";
		break;
	}
	case 8:
	{
		folder = "Lst1600";
		part = "lst1600";
		break;
	}
	case 9:
	{
		folder = "chiC1";
		part = "chic1";
		break;
	}
	}
	//Set up logging
	if(isData && logFlag)
	{
		//gROOT->ProcessLine(Form(".> logs/data/JpsiLambda/run%d/trigger_%d_log.txt",run,year));
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/Trigger_%d_log.txt",
		                             run,year),"w");
	}
	else if(!isData && logFlag)
	{
		//gROOT->ProcessLine(Form(".> logs/mc/JpsiLambda/%s/run%d/trigger_log.txt",folder,run));
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Trigger_log.txt",
		                             folder,run),"w");
	}
	// else if(!isData && logFlag && mcType == 2)
	// {
	//      gROOT->ProcessLine(Form(".> logs/mc/JpsiLambda/JpsiSigma/run%d/trigger_log.txt",run));
	// }
	// else if(!isData && logFlag && mcType == 3)
	// {
	//      gROOT->ProcessLine(Form(".> logs/mc/JpsiLambda/JpsiXi/run%d/trigger_log.txt",run));
	// }

	TStopwatch sw;
	sw.Start();

	cout<<"******************************************"<<endl;
	cout<<"==> Starting Trigger: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	gSystem->Exec("date");
	cout<<"******************************************"<<endl;

	//gROOT->ProcessLine(".L scripts/CollateFiles.C++");//Chains/hadds input files

	Int_t entries_init = 0, entries_final = 0, entries_gen = 0;
	Float_t eff_excl   = 0., eff_excl_err = 0.;
	Float_t eff_incl   = 0., eff_incl_err = 0.;

	Bool_t hlt1DiMuonHighMass = false, hlt1TrackMuon      = false;
	Bool_t hlt1TrackAllL0     = false, hlt2DiMuonDetached = false;
	Bool_t collateFlag        = false;//NOTE collateFlag is FALSE!!!
	/*f you don't want to re-collate MC, set this to false.
	   For example, if only the trigger condition changes.*/

	TFile *fileOut = nullptr, *fileIn     = nullptr;
	TFile *fileIn_gen = nullptr;
	TTree *treeIn  = nullptr, *treeIn_gen = nullptr;
	TTree *treeOut = nullptr, *myTree     = nullptr;
	const char *triggerCut = "(Lb_Hlt1DiMuonHighMassDecision_TOS==1||"
	                         "Lb_Hlt1TrackMuonDecision_TOS==1||"
	                         "Lb_Hlt1TrackAllL0Decision_TOS==1)"
	                         "&&(Lb_Hlt2DiMuonDetachedJPsiDecision_TOS==1)";

	/*The number of generated MC entries in every MC case will be written out to
	   this file, so that it can be accessed for calculating exclusive efficiencies
	   later*/
	ofstream genFile;
	// Set up input, output
	if(isData)
	{
		Double_t lumi     = 0., lumiErr = 0.;
		Double_t lumiMean = 0., lumiErrMean = 0.;
		Int_t lumiEntries = 0, lumiErrEntries = 0;

		TH1D *lumiHist = nullptr, *lumiErrHist = nullptr;
		TChain *h1 = new TChain("Lb2JpsiLTree/MyTuple");
		TChain *h2 = new TChain("GetIntegratedLuminosity/LumiTuple");

		//CollateFiles will cd to the massdump folder and then back to JpsiLambda_RESTART
		CollateFiles(run, year, isData, mcType, &h1, &h2, testing, loose, logFlag);

		fileOut = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered_%d.root",
		                         run,year),"RECREATE");

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
		genFile.open(Form("logs/mc/JpsiLambda/%s/run%d/gen_log.txt",folder,run));

		if(collateFlag)
		{
			CollateFiles(run, year, isData, mcType);
		}
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s.root",
		                          folder,run,part));
		fileIn_gen = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/lst1405.root",
		                              folder,run));
		treeIn_gen = (TTree*)fileIn_gen->Get("MCTuple/MCDecayTree");
		// treeIn     = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
		treeIn     = (TTree*)fileIn->Get("MyTuple");//TEMP


		fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_triggered.root",
		                         folder,run,part),"RECREATE");

		entries_gen = treeIn_gen->GetEntries();
		genFile<<entries_gen<<endl;
		genFile.close();

		entries_init = treeIn->GetEntries();
		myTree = treeIn;
	}//end MC block
	 //end setup of input, output

	cout<<"******************************************"<<endl;
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

			cout<<"******************************************"<<endl;
			cout<<"Trigger cut made with exclusive efficiency = "<<
			        eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
		if(entries_gen != 0)
		{
			eff_incl = (Float_t)entries_final*100/entries_gen;
			eff_incl_err = sqrt( eff_incl*(100.0-eff_incl)/entries_gen);

			cout<<"******************************************"<<endl;
			cout<<"Trigger cut made with inclusive efficiency = "<<
			        eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
	}

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	sw.Stop();
	cout << "==> Trigger is done! Huzzah!: "; sw.Print();

	//if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);
	//if(logFlag) gSystem->Exec("cat trigger_log.txt");
}
