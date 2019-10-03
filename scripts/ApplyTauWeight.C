#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TSystem.h"
#include <iostream>

using namespace std;

void ApplyTauWeight(Int_t run = 1, Int_t mcType = 0, Bool_t isGen = false, Bool_t logFlag = true)
/*
    mcType = 1 for Lb -> J/psi Lambda MC
    mcType = 2 for Lb -> J/psi Sigma MC        (reco'd JpsiLambda)
    mcType = 4 for Lb -> J/psi Lambda*(1405)MC (reco'd JpsiLambda)
    mcType = 5 for Lb -> J/psi Lambda*(1520)MC (reco'd JpsiLambda)
    mcType = 6 for Lb -> J/psi Lambda*(1600)MC (reco'd JpsiLambda)

    isGen = true to reweight generator MC, false to reweight reconstructed MC
 */
{
	TStopwatch sw;
	sw.Start();

	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;

	const char *folder = "", *part = "";

	Double_t Lb_TAU = 0.0;

	Int_t nEntries = 0;

	Float_t wt_tau = 0.0, wt_tau_plus = 0.0, wt_tau_minus = 0.0;

	cout<<"******** Processing Run "<<run<<" MC Type "<<mcType<<" isGen "<<isGen<<endl;

	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part   = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part   = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part   = "jpsisigma";
		break;
	}
	case 4:
	{
		folder = "Lst1405";
		part   = "lst1405";
		break;
	}
	case 5:
	{
		folder = "Lst1520";
		part   = "lst1520";
		break;
	}
	case 6:
	{
		folder = "Lst1600";
		part   = "lst1600";
		break;
	}
	default:
	{
		cout<<"$$$ MC Type doesn't match any of the allowed cases. Exiting! $$$"<<endl;
		exit(1);
	}
	}

	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
	if(logFlag) gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Alias.txt",folder,run),"w");

	if(!isGen)
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_pidgen.root",folder,run,part),"UPDATE");
	else
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s.root",folder,run,part),"READ");

	if(!isGen)
	{
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->SetBranchAddress("Lb_TAU",&Lb_TAU);
		// fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_rec.root",folder,run),"RECREATE");
		// treeOut = new TTree("MyTuple","MyTuple");
		fileOut = fileIn;
		treeOut = treeIn;
	}
	else
	{
		treeIn = (TTree*)fileIn->Get("MCTuple/MCDecayTree");
		treeIn->SetBranchAddress("Lambda_b0_TRUETAU",&Lb_TAU);
		fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_gen.root",folder,run),"RECREATE");
		treeOut = new TTree("MyTuple","MyTuple");
	}

	TBranch *wtBranch       = treeOut->Branch("wt_tau",&wt_tau,"wt_tau/F");
	TBranch *wtBranch_plus  = treeOut->Branch("wt_tau_plus",&wt_tau_plus,"wt_tau_plus/F");
	TBranch *wtBranch_minus = treeOut->Branch("wt_tau_minus",&wt_tau_minus,"wt_tau_minus/F");

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
		//World average for Lb lifetime is (1.470 +/- 0.010) ps  (0.7% uncertainty).
		if(run == 1)
		{
			wt_tau       = exp(-1000*Lb_TAU/1.470)/exp(-1000*Lb_TAU/1.425);
			wt_tau_plus  = exp(-1000*Lb_TAU/1.480)/exp(-1000*Lb_TAU/1.425);
			wt_tau_minus = exp(-1000*Lb_TAU/1.460)/exp(-1000*Lb_TAU/1.425);
		}
		else if(run == 2)
		{
			wt_tau       = exp(-1000*Lb_TAU/1.470)/exp(-1000*Lb_TAU/1.450);
			wt_tau_plus  = exp(-1000*Lb_TAU/1.480)/exp(-1000*Lb_TAU/1.450);
			wt_tau_minus = exp(-1000*Lb_TAU/1.460)/exp(-1000*Lb_TAU/1.450);
		}
		if(Lb_TAU > 1)
		{
			wt_tau = 1.0;
		}
		if(!isGen)
		{
			wtBranch->Fill();
			wtBranch_plus->Fill();
			wtBranch_minus->Fill();
		}
		else
			treeOut->Fill();
	}

	if(!isGen)
	{
		fileOut->cd();
		treeOut->Write("",TObject::kOverwrite);
		fileOut->Close();
	}
	else
	{
		fileOut->cd();
		treeOut->Write();
		fileOut->Close();
	}

	sw.Stop();
	cout << "==> ApplyTauWeights is done! Huzzah!: "; sw.Print();
	if(logFlag) gSystem->RedirectOutput(0);
}
