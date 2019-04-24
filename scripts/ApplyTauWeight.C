#include "TFile.h"
#include "TTree.h"
#include "iostream"

using namespace std;

void ApplyTauWeight(Int_t run = 1, Int_t mcType = 0, Bool_t isGen = false)
//mcType = 1 for JpsiLambda MC
//mcType = 2 for JpsiSigma MC
//isGen = true when you are RW'ing Generator Info
{
	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;
	const char *folder = "", *part = "";
	Double_t Lb_TAU = 0.0;
	Int_t nEntries = 0;
	Float_t wt_tau = 0.0;

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
	}
	fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/%s.root",folder,run,part));

	if(!isGen)
	{
		treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
		treeIn->SetBranchAddress("Lb_TAU",&Lb_TAU);
		fileOut = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_rec.root",folder,run),"RECREATE");
		treeOut = new TTree("MyTuple","MyTuple");
	}
	else
	{
		treeIn = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
		treeIn->SetBranchAddress("Lambda_b0_TRUETAU",&Lb_TAU);
		fileOut = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_gen.root",folder,run),"RECREATE");
		treeOut = new TTree("MyTuple","MyTuple");
	}

	treeOut->Branch("wt_tau",&wt_tau,"wt_tau/F");

	nEntries = treeIn->GetEntries();

	for(Int_t i=0; i<nEntries; i++)
	{
		if(isGen)
		{
			if(i%100000==0)
			{
				cout<<i<<endl;
			}
		}
		if(!isGen)
		{
			if(i%10000==0)
			{
				cout<<i<<endl;
			}
		}

		treeIn->GetEntry(i);
		if(run == 1)
		{
			wt_tau = exp(-1000*Lb_TAU/1.470)/exp(-1000*Lb_TAU/1.425);
		}
		else if(run == 2)
		{
			wt_tau = exp(-1000*Lb_TAU/1.470)/exp(-1000*Lb_TAU/1.450);
		}
		if(Lb_TAU > 1)
		{
			wt_tau = 1.0;
		}
		treeOut->Fill();
	}
	fileOut->Write();
	fileOut->Close();


}