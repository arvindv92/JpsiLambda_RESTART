#include "MVmodel.C"
#include "BONN.C"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TEntryList.h"
#include "iostream"
#include <vector>
#include <algorithm>

using namespace std;

//SOMETHING IS GOING WRONG. LOT OF THE WEIGHTS ARE ZERO.
void reweight_Lst1405(Int_t run = 1)
{
	TFile *fileIn        = nullptr, *fileIn_gen = nullptr;// *mcfile_MV       = nullptr, *mcfile_BONN = nullptr;
	TFile *genFile_MV    = nullptr, *genFile_BONN    = nullptr;
	TTree *treeIn        = nullptr, *treeIn_gen      = nullptr;
//TTree *treeout_MV    = nullptr, *treeout_BONN    = nullptr;
	TTree *genTreeout_MV = nullptr, *genTreeout_BONN = nullptr;

	Double_t Lst_PE       = 0.0, Lst_PX         = 0.0, Lst_PY     = 0.0, Lst_PZ   = 0.0;
	Double_t Jpsi_PE      = 0.0, Jpsi_PX        = 0.0, Jpsi_PY    = 0.0, Jpsi_PZ  = 0.0;
	Double_t Lambda_PE    = 0.0, Lambda_PX      = 0.0, Lambda_PY  = 0.0,Lambda_PZ = 0.0;
	Double_t mass_gen     = 0.0, MVweight       = 0.0, BONNweight = 0.0;
	Double_t MVweight_gen = 0.0, BONNweight_gen = 0.0;

	Int_t nentries = 0, nentries_gen = 0, massbin_MV = 0, massbin_BONN = 0;
	Int_t low_MV   = 0, high_MV      = 0, low_BONN   = 0, high_BONN    = 0;//limits for Lst mass in which the models work
	Int_t nbins_MV = 0, nbins_BONN   = 0, bkgcat     = 0;

	ULong64_t evtno_rec = 0, evtno_gen = 0;
	UInt_t runno_rec    = 0, runno_gen = 0;

	vector<ULong64_t> evtno;
	vector<UInt_t> runno;
	vector<Double_t>lstmass_gen;

	low_MV    = 1300;//these values are rough estimates. Look at paper and get precise limits
	high_MV   = 2523;
	low_BONN  = 1328;
	high_BONN = 2519;

	Int_t binwidth = 1;

	nbins_MV   = (Int_t)(high_MV - low_MV)/binwidth;
	nbins_BONN = (Int_t)(high_BONN - low_BONN)/binwidth;

	TH1D *theory_MVmodel   = new TH1D("theory_MVmodel","",nbins_MV,low_MV,high_MV);
	TH1D *theory_BONNmodel = new TH1D("theory_BONNmodel","",nbins_BONN,low_BONN,high_BONN);

	fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/lst1405_cutoutks_LL_ZeroTracks.root",run),"UPDATE");
	treeIn     = (TTree*)fileIn->Get("MyTuple");

	fileIn_gen = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/lst1405.root",run),"READ");
	treeIn_gen = (TTree*)fileIn_gen->Get("MCTuple/MCDecayTree");

	// genFile_MV = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/RW/lst1405_gen_MV.root",run),"RECREATE");
	// genTreeout_MV = new TTree("MCDecayTree","");

	// genFile_BONN = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/RW/lst1405_gen_BONN.root",run),"RECREATE");
	// genTreeout_BONN = new TTree("MCDecayTree","");

	cout<<"Done copying Trees"<<endl;

	treeIn_gen->Draw(Form("sqrt(pow(Lambda_1405_0_TRUEP_E,2) - pow(Lambda_1405_0_TRUEP_X,2) - pow(Lambda_1405_0_TRUEP_Y,2) - pow(Lambda_1405_0_TRUEP_Z,2) )>>hmv(%d,%d,%d)",nbins_MV,low_MV,high_MV),"","goff");
	treeIn_gen->Draw(Form("sqrt(pow(Lambda_1405_0_TRUEP_E,2) - pow(Lambda_1405_0_TRUEP_X,2) - pow(Lambda_1405_0_TRUEP_Y,2) - pow(Lambda_1405_0_TRUEP_Z,2) )>>hbonn(%d,%d,%d)",nbins_BONN,low_BONN,high_BONN),"","goff");

	TH1D *hmv = (TH1D*)gDirectory->Get("hmv");
	hmv->Scale(1.0/hmv->Integral());

	TH1D *hbonn = (TH1D*)gDirectory->Get("hbonn");
	hbonn->Scale(1.0/hbonn->Integral());

	for(Int_t i=0; i<nbins_MV; i++)
	{
		theory_MVmodel->Fill((i+low_MV)*1.0,(Double_t)MVmodel(i+low_MV));
	}
	for(Int_t i=0; i<nbins_BONN; i++)
	{
		theory_BONNmodel->Fill((i+low_BONN)*1.0,(Double_t)invariantmassdistribution_sig0pi0((i+low_BONN)*1.0/1000));//takes input in GeV
	}

	theory_MVmodel->Scale(1.0/theory_MVmodel->Integral());
	theory_BONNmodel->Scale(1.0/theory_BONNmodel->Integral());

	theory_MVmodel->SetLineColor(kMagenta+2);
	hmv->SetLineColor(kRed);
	hbonn->SetLineColor(kRed);

	TH1D *wts_MV = (TH1D*)theory_MVmodel->Clone("wts_MV");
	TH1D *wts_BONN = (TH1D*)theory_BONNmodel->Clone("wts_BONN");

	wts_MV->Divide(hmv);
	wts_BONN->Divide(hbonn);

	nentries = treeIn->GetEntries();
	nentries_gen = treeIn_gen->GetEntries();

	treeIn_gen->SetBranchAddress("Lambda_1405_0_TRUEP_E",&Lst_PE);
	treeIn_gen->SetBranchAddress("Lambda_1405_0_TRUEP_X",&Lst_PX);
	treeIn_gen->SetBranchAddress("Lambda_1405_0_TRUEP_Y",&Lst_PY);
	treeIn_gen->SetBranchAddress("Lambda_1405_0_TRUEP_Z",&Lst_PZ);

	treeIn_gen->SetBranchAddress("J_psi_1S_TRUEP_E",&Jpsi_PE);
	treeIn_gen->SetBranchAddress("J_psi_1S_TRUEP_X",&Jpsi_PX);
	treeIn_gen->SetBranchAddress("J_psi_1S_TRUEP_Y",&Jpsi_PY);
	treeIn_gen->SetBranchAddress("J_psi_1S_TRUEP_Z",&Jpsi_PZ);

	treeIn_gen->SetBranchAddress("Lambda0_TRUEP_E",&Lambda_PE);
	treeIn_gen->SetBranchAddress("Lambda0_TRUEP_X",&Lambda_PX);
	treeIn_gen->SetBranchAddress("Lambda0_TRUEP_Y",&Lambda_PY);
	treeIn_gen->SetBranchAddress("Lambda0_TRUEP_Z",&Lambda_PZ);

	treeIn_gen->SetBranchAddress("eventNumber",&evtno_gen);
	treeIn_gen->SetBranchAddress("runNumber",&runno_gen);

	treeIn->SetBranchAddress("eventNumber",&evtno_rec);
	treeIn->SetBranchAddress("runNumber",&runno_rec);
	treeIn->SetBranchAddress("Lb_BKGCAT",&bkgcat);

	TBranch *MVwtbranch = treeIn->Branch("MVweight",&MVweight,"MVweight/D");
	TBranch *BONNwtbranch = treeIn->Branch("BONNweight",&BONNweight,"BONNweight/D");

	// TBranch *MVwtbranch_gen = genTreeout_MV->Branch("MVweight",&MVweight_gen,"MVweight/D");
	// TBranch *BONNwtbranch_gen = genTreeout_BONN->Branch("BONNweight",&BONNweight_gen,"BONNweight/D");
	Int_t ctr = 0;

	//First loop over reco tree
	cout<<"First loop over reco tree"<<endl;
	for(Int_t i=0; i<nentries; i++)
	{
		if(i%10000 == 0)
			cout<<i<<endl;
		treeIn->GetEntry(i);
		if(bkgcat==50)
		{
			evtno.push_back(evtno_rec);
			runno.push_back(runno_rec);
		}
	}

	cout<<"Length of evtno vector = "<<evtno.size()<<" and length of runno vector = "<<runno.size()<<endl;
	if(evtno.size()!=runno.size())
	{
		cout<<"Vector lengths dont match!! Exiting!"<<endl;
		exit(1);
	}

	int len = evtno.size();
	lstmass_gen.assign(len,0.0);

	//Loop over generated tree
	cout<<"Looping over generated tree which has "<<nentries_gen<<" entries"<<endl;
	Int_t ctr_notfound = 0;
	int index_rno;
	int index_eno;
	for(Int_t i=0; i<nentries_gen; i++)
	{
		if(i%100000 == 0)
			cout<<i<<endl;
		treeIn_gen->GetEntry(i);

		//*****Reweight the generator MC also*****
		// mass_gen       = sqrt(Lst_PE*Lst_PE - Lst_PX*Lst_PX - Lst_PY*Lst_PY - Lst_PZ*Lst_PZ);
		// massbin_MV     = wts_MV->FindBin(mass_gen);
		// massbin_BONN   = wts_BONN->FindBin(mass_gen);
		// MVweight_gen   = wts_MV->GetBinContent(massbin_MV);
		// BONNweight_gen = wts_BONN->GetBinContent(massbin_BONN);
		// genTreeout_MV->Fill();
		// genTreeout_BONN->Fill();
		//****************************************

		//The purpose of the loop below is to fill the generated Lst mass corresponding to a reconstructed event
		for(Int_t j=0; j<len; j++)
		{
			if((runno_gen == runno[j]) && (evtno_gen == evtno[j]))
			{
				// cout<<"Lst_PE = "<<Lst_PE<<"\tLst_PX = "<<Lst_PX<<"\tLst_PY = "<<Lst_PY<<"\tLst_PZ = "<<Lst_PZ<<endl;
				lstmass_gen[j] = sqrt(Lst_PE*Lst_PE - Lst_PX*Lst_PX - Lst_PY*Lst_PY - Lst_PZ*Lst_PZ);
			}
		}
	}

	cout<<"Done looping over generated tree."<<endl;

	// cout<<"ctr_notfound = "<<ctr_notfound<<endl;
	cout<<"Length of evtno vector = "<<evtno.size()<<" and length of lstmass vector = "<<lstmass_gen.size()<<endl;
	//How many entries of lstmass_gen are zero at this stage?
	Int_t zeroct = 0;
	Int_t newct = 0;

	for(Int_t i=0; i<len; i++)
	{
		if(lstmass_gen[i] == 0.0)
		{
			zeroct++;
		}
	}
	cout<<zeroct<<" entries of lstmass_gen are 0.0"<<endl;
	cout<<"Looping over reco tree again"<<endl;

	ctr = 0;
	for(Int_t i=0; i<nentries; i++)
	{
		//Hopefully the entries in the vector should be in the same order as the
		// bkgcat=50 events in the reco tree.
		treeIn->GetEntry(i);
		if(i%10000 == 0)
			cout<<i<<endl;
		if(bkgcat!=50)
		{
			MVweight   = 1.0;
			BONNweight = 1.0;
			MVwtbranch->Fill();
			BONNwtbranch->Fill();
			continue;
		}
		else
		{
			mass_gen = lstmass_gen[ctr];
			massbin_MV = wts_MV->FindBin(mass_gen);
			massbin_BONN = wts_BONN->FindBin(mass_gen);
			MVweight = wts_MV->GetBinContent(massbin_MV);
			BONNweight = wts_BONN->GetBinContent(massbin_BONN);
			MVwtbranch->Fill();
			BONNwtbranch->Fill();
			ctr++;
		}

	}
	cout<<"ctr = "<<ctr<<endl;

	fileIn->cd();
	treeIn->Write("",TObject::kOverwrite);
	fileIn->Close();

	// cout<<"Writing out reweighted generator"<<endl;
	// genFile_MV->cd();
	// genTreeout_MV->Write();
	// genFile_MV->Close();

	// genFile_BONN->cd();
	// genTreeout_BONN->Write();
	// genFile_BONN->Close();
}
