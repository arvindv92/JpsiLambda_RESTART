/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply trigger cuts on Xib -> J/psi Xi data/MC coming out of DaVinci.
 *********************************/
#include "Trigger_JpsiXi.h"
/* DOCUMENTATION
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   testing = true to run only over a subset of data
 */
void Trigger_JpsiXi(Int_t run, Int_t year, Bool_t isData, Bool_t testing, Bool_t logFlag)
{
	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");

	//Set up logging
	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiXi/run%d/Trigger_%d_log.txt",
		                             run,year),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiXi/run%d/Trigger_log.txt",
		                             run),"w");
	}

	TStopwatch sw;
	sw.Start();

	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting JpsiXi Trigger: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	gSystem->Exec("date");
	cout<<"******************************************"<<endl;

//	gROOT->ProcessLine(".L scripts/CollateFiles_JpsiXi.C++");//Chains/hadds input files

	Int_t entries_gen      = 0, entries_init_TM = 0;
	Int_t entries_init     = 0, entries_final   = 0;
	Int_t entries_final_TM = 0;

	Float_t eff_excl   = 0., eff_excl_err = 0.;
	Float_t eff_incl   = 0., eff_incl_err = 0.;

	Bool_t hlt1TrackAllL0     = false, hlt2DiMuonDetached = false;
	Bool_t hlt1DiMuonHighMass = false, hlt1TrackMuon      = false;

	TFile *fileOut = nullptr, *fileIn     = nullptr;
	TTree *treeIn  = nullptr, *treeIn_gen = nullptr;
	TTree *treeOut = nullptr, *myTree     = nullptr;

	const char *triggerCut = "(Xib_Hlt1DiMuonHighMassDecision_TOS==1||"
	                         "Xib_Hlt1TrackMuonDecision_TOS==1||"
	                         "Xib_Hlt1TrackAllL0Decision_TOS==1)&&"
	                         "(Xib_Hlt2DiMuonDetachedJPsiDecision_TOS==1)";

	ofstream genFile;//contains no of generated candidates in MC.

	// Set up input, output
	if(isData)
	{
		Double_t lumi     = 0., lumiErr       = 0.;
		Double_t lumiMean = 0., lumiErrMean   = 0.;
		Int_t lumiEntries = 0, lumiErrEntries = 0;

		TH1D *lumiHist = nullptr, *lumiErrHist = nullptr;

		TChain *h1 = new TChain("Xib2JpsiXiTree/MyTuple");
		TChain *h2 = new TChain("GetIntegratedLuminosity/LumiTuple");

		//CollateFiles will cd to the massdump folder and then back to JpsiLambda_RESTART
		// CollateFiles_JpsiXi(run, year, isData, &h1, &h2, testing, logFlag);
		gSystem->Exec(".L CollateFiles_JpsiXi.C");
		gSystem->Exec(Form("CollateFiles_JpsiXi(%d,%d,%d,%p,%p,%d,%d)",run,year,isData,(void*)&h1,(void*)&h2,testing,logFlag));

		fileOut = new TFile(Form("rootFiles/dataFiles/JpsiXi/run%d/jpsixi_triggered_%d.root",
		                         run,year),"RECREATE");

		h2->Draw("IntegratedLuminosity>>lumiHist","","goff");
		h2->Draw("IntegratedLuminosityErr>>lumiErrHist","","goff");

		lumiHist    = (TH1D*)gDirectory->Get("lumiHist");
		lumiErrHist = (TH1D*)gDirectory->Get("lumiErrHist");

		lumiMean       = lumiHist->GetMean();
		lumiEntries    = lumiHist->GetEntries();
		lumiErrMean    = lumiErrHist->GetMean();
		lumiErrEntries = lumiErrHist->GetEntries();

		lumi    = lumiMean*lumiEntries;
		lumiErr = lumiErrMean*lumiErrEntries;

		cout<<"Processing "<<lumi<<" +/- "<<lumiErr<< "luminosity"<<endl;

		myTree = (TTree*)h1;
	}//end Data block
	else //MC
	{
		genFile.open(Form("logs/mc/JpsiXi/run%d/gen_log.txt",run));

		cout<<"PROCESSING MC for Run "<<run<<" Jpsi Xi"<<endl;

		// if(collateFlag) CollateFiles_JpsiXi(run, year, isData);

		fileIn     = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi.root",
		                              run));
		if(!fileIn)
		{
			cout<<"Unable to open MC input file for Trigger ! Exiting!"<<endl;
			exit(1);
		}
		treeIn_gen = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
		treeIn     = (TTree*)fileIn->Get("Xib2JpsiXiTree/MyTuple");

		fileOut = new TFile(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_triggered.root",
		                         run),"RECREATE");

		entries_gen = treeIn_gen->GetEntries();
		genFile<<entries_gen<<endl;
		genFile.close();

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

	myTree->SetBranchAddress("Xib_Hlt1DiMuonHighMassDecision_TOS",&hlt1DiMuonHighMass);
	myTree->SetBranchAddress("Xib_Hlt1TrackMuonDecision_TOS",&hlt1TrackMuon);
	myTree->SetBranchAddress("Xib_Hlt1TrackAllL0Decision_TOS",&hlt1TrackAllL0);
	myTree->SetBranchAddress("Xib_Hlt2DiMuonDetachedJPsiDecision_TOS",&hlt2DiMuonDetached);

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

	if(!isData)
	{
		entries_init_TM = myTree->GetEntries("abs(Xib_TRUEID) == 5132 || (Xib_BKGCAT == 60 && (fabs(Xib_DTF_M_JpsiXiLConstr - 5795) < 35))");
		entries_final_TM = treeOut->GetEntries("abs(Xib_TRUEID) == 5132 || (Xib_BKGCAT == 60 && (fabs(Xib_DTF_M_JpsiXiLConstr - 5795) < 35))");
	}

	if(!isData)//Calculate exclusive and inclusive efficiencies for MC
	{
		if(entries_init != 0)
		{
			eff_excl = (Float_t)entries_final_TM*100/entries_init_TM;
			eff_excl_err = sqrt( eff_excl*(100.0-eff_excl)/entries_init_TM);
		}
		cout<<"******************************************"<<endl;
		cout<<"Trigger cut made with exclusive efficiency = "<<eff_excl
		    <<"% +/- " <<eff_excl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(entries_gen != 0)
		{
			eff_incl = (Float_t)entries_final_TM*100/entries_gen;
			eff_incl_err = sqrt( eff_incl*(100.0-eff_incl)/entries_gen);
		}
		cout<<"******************************************"<<endl;
		cout<<"Trigger cut made with inclusive efficiency = "<<eff_incl
		    <<"% +/- " <<eff_incl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;
	}

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	sw.Stop();
	cout << "==> Trigger is done! Huzzah!: "; sw.Print();

	if(logFlag) gSystem->RedirectOutput(0);
	//if(logFlag) gSystem->Exec("cat trigger_log.txt");
}
