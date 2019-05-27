#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TBranch.h"
#include "iostream"

using namespace std;

void ApplyTrackingCorr(Int_t run = 1)
{
	TFile *wtFile = nullptr, *fileIn = nullptr;
	TTree *treeIn = nullptr;
	TH2D *trackEff = nullptr;

	Double_t BachPi_P = 0.0, BachPi_ETA = 0.0;
	Double_t wt = 1.0, wtErr = 0.0;
	Int_t nEntries = 0, bin = 0;

	TBranch *wtBranch = nullptr, *wtErrBranch = nullptr;

	if(run == 1)
	{
		wtFile = TFile::Open("ratio2012S20.root");
	}
	if(run == 2)
	{
		wtFile = TFile::Open("Ratio_Long_P-ETA_2016_25ns.root");
	}
	trackEff = (TH2D*)wtFile->Get("Ratio");

	fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cut_LL.root",run),"UPDATE");
	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->SetBranchAddress("BachPi_P",&BachPi_P);
	treeIn->SetBranchAddress("BachPi_ETA",&BachPi_ETA);

	nEntries = treeIn->GetEntries();

	wtBranch    = treeIn->Branch("wt_tracking",&wt,"wt_tracking/D");
	wtErrBranch = treeIn->Branch("wtErr_tracking",&wtErr,"wtErr_tracking/D");

	for(Int_t i=0; i<nEntries; i++)
	{
		treeIn->GetEntry(i);
		if(i%10000 == 0)
			cout<<i<<endl;

		if(BachPi_ETA > 1.9 && BachPi_ETA < 4.9 && BachPi_P > 5000 && BachPi_P < 200000)
		{
			if(run==1) bin = trackEff->FindBin(0.001*BachPi_P,BachPi_ETA); // p in GeV
			if(run==2) bin = trackEff->FindBin(BachPi_P,BachPi_ETA); // p in MeV
			wt    = trackEff->GetBinContent(bin);
			wtErr = trackEff->GetBinError(bin);
		}
		else
		{
			wt = 1.0;
			wtErr = 0.05;
		}
		wtBranch->Fill();
		wtErrBranch->Fill();
	}

	fileIn->cd();
	treeIn->Write("",TObject::kOverwrite);
	fileIn->Close();
	wtFile->Close();
}
