#include "TFile.h"
#include "TTree.h"
#include "iostream"

using namespace std;

void ApplyTauWeight(Int_t run = 1, Int_t mcType = 0, Bool_t isGen = false)
//mcType = 1 for JpsiLambda MC
//mcType = 2 for JpsiSigma MC
//mcType = 6 for JpsiLst(1405) MC
//mcType = 7 for JpsiLst(1520) MC
//mcType = 8 for JpsiLst(1600) MC
//mcType = 9 for chiC1 Lambda MC
//isGen = true when you are RW'ing Generator Info
{
	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;
	const char *folder = "", *part = "";
	Double_t Lb_TAU = 0.0;
	Int_t nEntries = 0;
	Float_t wt_tau = 0.0;

	cout<<"********Processing Run "<<run<<" MC Type "<<mcType<<" isGen "<<isGen<<endl;
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
	case 6:
	{
		folder = "Lst1405";
		part = "lst1405";
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
	if(!isGen)
		fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_pidgen.root",folder,run,part),"UPDATE");
	else
		fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/%s.root",folder,run,part));

	if(!isGen)
	{
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->SetBranchAddress("Lb_TAU",&Lb_TAU);
		// fileOut = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_rec.root",folder,run),"RECREATE");
		// treeOut = new TTree("MyTuple","MyTuple");
		fileOut = fileIn;
		treeOut = treeIn;
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

	if(!isGen)
	{
		fileOut->cd();
		treeOut->Write("",TObject::kOverwrite);
	}
	else
	{
		fileOut->Write();
		fileOut->Close();
	}

}
