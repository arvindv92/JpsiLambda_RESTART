#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include <iostream>

using namespace std;

#define Open TFile::Open

void PolarizationEffs(Int_t run = 1)
{
	TFile *fileIn_nonZero_JpsiLambda = 0, *fileIn_nonZero_JpsiSigma = 0;
	TTree *treeIn_nonZero_JpsiLambda = 0, *treeIn_nonZero_JpsiSigma = 0;

	TFile *fileIn_Zero_JpsiLambda = 0, *fileIn_Zero_JpsiSigma = 0;
	TTree *treeIn_Zero_JpsiLambda = 0, *treeIn_Zero_JpsiSigma = 0;

	TFile *fileIn_gen_JpsiLambda = 0, *fileIn_gen_JpsiSigma = 0;
	TTree *treeIn_gen_JpsiLambda = 0, *treeIn_gen_JpsiSigma = 0;

	Int_t isoConf = -1;
	if(run == 1)
	{
		isoConf = 2;
	}
	else if(run == 2)
	{
		isoConf = 1;
	}

	Float_t bdtCut_nonZero = -1.0, bdtCut_Zero = -1.0;

	if(run == 1)
	{
		bdtCut_nonZero = 0.475;
		bdtCut_Zero    = 0.365;
	}
	else if(run == 2)
	{
		bdtCut_nonZero = 0.555;
		bdtCut_Zero    = 0.495;
	}

	Double_t Lb_PX_gen_JpsiLambda = 0.0;
	Double_t Lb_PY_gen_JpsiLambda = 0.0;
	Double_t Lb_PZ_gen_JpsiLambda = 0.0;
	Double_t Lb_PE_gen_JpsiLambda = 0.0;

	Double_t L_PX_gen_JpsiLambda = 0.0;
	Double_t L_PY_gen_JpsiLambda = 0.0;
	Double_t L_PZ_gen_JpsiLambda = 0.0;
	Double_t L_PE_gen_JpsiLambda = 0.0;

	Double_t gbWt_gen_JpsiLambda = 0.0;
	Float_t tauWt_gen_JpsiLambda = 0.0;

	Int_t nEntries_gen_JpsiLambda = 0;

	Double_t Lb_PX_rec_JpsiLambda = 0.0;
	Double_t Lb_PY_rec_JpsiLambda = 0.0;
	Double_t Lb_PZ_rec_JpsiLambda = 0.0;
	Double_t Lb_PE_rec_JpsiLambda = 0.0;

	Double_t L_PX_rec_JpsiLambda = 0.0;
	Double_t L_PY_rec_JpsiLambda = 0.0;
	Double_t L_PZ_rec_JpsiLambda = 0.0;
	Double_t L_PE_rec_JpsiLambda = 0.0;

	Double_t gbWt_rec_JpsiLambda = 0.0;
	Float_t tauWt_rec_JpsiLambda = 0.0;

	Int_t nEntries_rec_JpsiLambda = 0;

	Double_t Lb_PX_gen_JpsiSigma = 0.0;
	Double_t Lb_PY_gen_JpsiSigma = 0.0;
	Double_t Lb_PZ_gen_JpsiSigma = 0.0;
	Double_t Lb_PE_gen_JpsiSigma = 0.0;

	Double_t L_PX_gen_JpsiSigma = 0.0;
	Double_t L_PY_gen_JpsiSigma = 0.0;
	Double_t L_PZ_gen_JpsiSigma = 0.0;
	Double_t L_PE_gen_JpsiSigma = 0.0;

	Double_t gbWt_gen_JpsiSigma = 0.0;
	Float_t tauWt_gen_JpsiSigma = 0.0;

	Int_t nEntries_gen_JpsiSigma = 0;

	Double_t Lb_PX_rec_JpsiSigma = 0.0;
	Double_t Lb_PY_rec_JpsiSigma = 0.0;
	Double_t Lb_PZ_rec_JpsiSigma = 0.0;
	Double_t Lb_PE_rec_JpsiSigma = 0.0;

	Double_t L_PX_rec_JpsiSigma = 0.0;
	Double_t L_PY_rec_JpsiSigma = 0.0;
	Double_t L_PZ_rec_JpsiSigma = 0.0;
	Double_t L_PE_rec_JpsiSigma = 0.0;

	Double_t gbWt_rec_JpsiSigma = 0.0;
	Float_t tauWt_rec_JpsiSigma = 0.0;

	Int_t nEntries_rec_JpsiSigma = 0;

	Double_t BDT_JpsiLambda = 0.0;
	Double_t BDT_JpsiSigma = 0.0;

	TLorentzVector Lb_gen_JpsiLambda, Lb_rec_JpsiLambda, L_gen_JpsiLambda, L_rec_JpsiLambda;
	TLorentzVector Lb_gen_JpsiSigma, Lb_rec_JpsiSigma, L_gen_JpsiSigma, L_rec_JpsiSigma;

	TVector3 n;//perpendicular to production plane of Lb
	TVector3 beam;
	beam.SetXYZ(0.0,0.0,1.0);

	Double_t theta_JpsiLambda_gen = 0.0, theta_JpsiLambda_rec = 0.0;
	Double_t theta_JpsiSigma_gen = 0.0, theta_JpsiSigma_rec = 0.0;

	TH1D *Lb_genHist_JpsiLambda = new TH1D("Lb_genHist_JpsiLambda","Lb_genHist_JpsiLambda",100,-1,1);
	TH1D *Lb_genHist_wt_JpsiLambda = new TH1D("Lb_genHist_wt_JpsiLambda","Lb_genHist_wt_JpsiLambda",100,-1,1);

	TH1D *Lb_genHist_JpsiSigma = new TH1D("Lb_genHist_JpsiSigma","Lb_genHist_JpsiSigma",100,-1,1);
	TH1D *Lb_genHist_wt_JpsiSigma = new TH1D("Lb_genHist_wt_JpsiSigma","Lb_genHist_wt_JpsiSigma",100,-1,1);

	TH1D *Lb_recHist_JpsiLambda = new TH1D("Lb_recHist_JpsiLambda","Lb_recHist_JpsiLambda",100,-1,1);
	TH1D *Lb_recHist_wt_JpsiLambda = new TH1D("Lb_recHist_wt_JpsiLambda","Lb_recHist_wt_JpsiLambda",100,-1,1);

	TH1D *Lb_recHist_JpsiSigma = new TH1D("Lb_recHist_JpsiSigma","Lb_recHist_JpsiSigma",100,-1,1);
	TH1D *Lb_recHist_wt_JpsiSigma = new TH1D("Lb_recHist_wt_JpsiSigma","Lb_recHist_wt_JpsiSigma",100,-1,1);

	fileIn_gen_JpsiLambda = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda.root",run));
	treeIn_gen_JpsiLambda = (TTree*)fileIn_gen_JpsiLambda->Get("MCTuple/MCDecayTree");
	treeIn_gen_JpsiLambda->AddFriend("MyTuple",Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/RW/tauWeights_gen.root",run));
	if(run == 1)
	{
		treeIn_gen_JpsiLambda->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/RW/gbWeights_gen_new.root");
	}
	else if(run == 2)
	{
		treeIn_gen_JpsiLambda->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/RW/gbWeights_gen.root");
	}
	nEntries_gen_JpsiLambda = treeIn_gen_JpsiLambda->GetEntries();

	fileIn_gen_JpsiSigma = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma.root",run));
	treeIn_gen_JpsiSigma = (TTree*)fileIn_gen_JpsiSigma->Get("MCTuple/MCDecayTree");
	treeIn_gen_JpsiSigma->AddFriend("MyTuple",Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/RW/tauWeights_gen.root",run));
	if(run == 1)
	{
		treeIn_gen_JpsiSigma->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run1/RW/gbWeights_gen_new.root");
	}
	else if(run == 2)
	{
		treeIn_gen_JpsiSigma->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run2/RW/gbWeights_gen.root");
	}
	nEntries_gen_JpsiSigma = treeIn_gen_JpsiSigma->GetEntries();

	fileIn_nonZero_JpsiLambda = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root",run));
	treeIn_nonZero_JpsiLambda = (TTree*)fileIn_nonZero_JpsiLambda->Get("MyTuple");
	treeIn_nonZero_JpsiLambda->AddFriend("MyTuple",Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_LL_FinalBDT2_iso%d_v0_noPID.root",run,isoConf));

	fileIn_Zero_JpsiLambda = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_cutoutks_LL_ZeroTracks_noPID.root",run));
	treeIn_Zero_JpsiLambda = (TTree*)fileIn_Zero_JpsiLambda->Get("MyTuple");
	treeIn_Zero_JpsiLambda->AddFriend("MyTuple",Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_zeroTracksLL_FinalBDT2_noPID.root",run));

	fileIn_nonZero_JpsiSigma = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_cutoutks_LL_nonZeroTracks_noPID.root",run));
	treeIn_nonZero_JpsiSigma = (TTree*)fileIn_nonZero_JpsiSigma->Get("MyTuple");
	treeIn_nonZero_JpsiSigma->AddFriend("MyTuple",Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_LL_FinalBDT2_iso%d_v0_noPID.root",run,isoConf));

	fileIn_Zero_JpsiSigma = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_cutoutks_LL_ZeroTracks_noPID.root",run));
	treeIn_Zero_JpsiSigma = (TTree*)fileIn_Zero_JpsiSigma->Get("MyTuple");
	treeIn_Zero_JpsiSigma->AddFriend("MyTuple",Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_zeroTracksLL_FinalBDT2_noPID.root",run));

	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda_b0_TRUEP_E",&Lb_PE_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda_b0_TRUEP_X",&Lb_PX_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda_b0_TRUEP_Y",&Lb_PY_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda_b0_TRUEP_Z",&Lb_PZ_gen_JpsiLambda);

	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_E",&L_PE_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_X",&L_PX_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_Y",&L_PY_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_Z",&L_PZ_gen_JpsiLambda);

	treeIn_gen_JpsiLambda->SetBranchAddress("wt_tau",&tauWt_gen_JpsiLambda);
	if(run == 1)
	{
		treeIn_gen_JpsiLambda->SetBranchAddress("gb_wts_new",&gbWt_gen_JpsiLambda);
	}
	else if(run == 2)
	{
		treeIn_gen_JpsiLambda->SetBranchAddress("gb_wts",&gbWt_gen_JpsiLambda);
	}

	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda_b0_TRUEP_E",&Lb_PE_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda_b0_TRUEP_X",&Lb_PX_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda_b0_TRUEP_Y",&Lb_PY_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda_b0_TRUEP_Z",&Lb_PZ_gen_JpsiSigma);

	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_E",&L_PE_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_X",&L_PX_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_Y",&L_PY_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_Z",&L_PZ_gen_JpsiSigma);

	treeIn_gen_JpsiSigma->SetBranchAddress("wt_tau",&tauWt_gen_JpsiSigma);
	if(run == 1)
	{
		treeIn_gen_JpsiSigma->SetBranchAddress("gb_wts_new",&gbWt_gen_JpsiSigma);
	}
	else if(run == 2)
	{
		treeIn_gen_JpsiSigma->SetBranchAddress("gb_wts",&gbWt_gen_JpsiSigma);
	}

	TFile *tempFile1 = new TFile("tempFile1.root","RECREATE");

	TTree *treeIn_nonZero_Lambda_cut = (TTree*)treeIn_nonZero_JpsiLambda->CopyTree(Form("BDT2 > %f",bdtCut_nonZero));
	TTree *treeIn_Zero_Lambda_cut    = (TTree*)treeIn_Zero_JpsiLambda->CopyTree(Form("BDT2 > %f",bdtCut_Zero));

	TList *list_Lambda = new TList;
	list_Lambda->Add(treeIn_nonZero_Lambda_cut);
	list_Lambda->Add(treeIn_Zero_Lambda_cut);

	TTree *combTree_Lambda = TTree::MergeTrees(list_Lambda);
	combTree_Lambda->SetName("combTree_Lambda");

	nEntries_rec_JpsiLambda = combTree_Lambda->GetEntries();

	TFile *tempFile2 = new TFile("tempFile2.root","RECREATE");

	TTree *treeIn_nonZero_Sigma_cut = (TTree*)treeIn_nonZero_JpsiSigma->CopyTree(Form("BDT2 > %f",bdtCut_nonZero));
	TTree *treeIn_Zero_Sigma_cut    = (TTree*)treeIn_Zero_JpsiSigma->CopyTree(Form("BDT2 > %f",bdtCut_Zero));

	TList *list_Sigma = new TList;
	list_Sigma->Add(treeIn_nonZero_Sigma_cut);
	list_Sigma->Add(treeIn_Zero_Sigma_cut);

	TTree *combTree_Sigma = TTree::MergeTrees(list_Sigma);
	combTree_Sigma->SetName("combTree_Sigma");

	nEntries_rec_JpsiSigma = combTree_Sigma->GetEntries();

	combTree_Lambda->SetBranchAddress("Lb_TRUEP_E",&Lb_PE_rec_JpsiLambda);
	combTree_Lambda->SetBranchAddress("Lb_TRUEP_X",&Lb_PX_rec_JpsiLambda);
	combTree_Lambda->SetBranchAddress("Lb_TRUEP_Y",&Lb_PY_rec_JpsiLambda);
	combTree_Lambda->SetBranchAddress("Lb_TRUEP_Z",&Lb_PZ_rec_JpsiLambda);

	combTree_Lambda->SetBranchAddress("L_TRUEP_E",&L_PE_rec_JpsiLambda);
	combTree_Lambda->SetBranchAddress("L_TRUEP_X",&L_PX_rec_JpsiLambda);
	combTree_Lambda->SetBranchAddress("L_TRUEP_Y",&L_PY_rec_JpsiLambda);
	combTree_Lambda->SetBranchAddress("L_TRUEP_Z",&L_PZ_rec_JpsiLambda);

	combTree_Lambda->SetBranchAddress("wt_tau",&tauWt_rec_JpsiLambda);
	if(run == 1)
	{
		combTree_Lambda->SetBranchAddress("gb_wts_new",&gbWt_rec_JpsiLambda);
	}
	else if(run == 2)
	{
		combTree_Lambda->SetBranchAddress("gb_wts",&gbWt_rec_JpsiLambda);
	}

	combTree_Sigma->SetBranchAddress("Lb_TRUEP_E",&Lb_PE_rec_JpsiSigma);
	combTree_Sigma->SetBranchAddress("Lb_TRUEP_X",&Lb_PX_rec_JpsiSigma);
	combTree_Sigma->SetBranchAddress("Lb_TRUEP_Y",&Lb_PY_rec_JpsiSigma);
	combTree_Sigma->SetBranchAddress("Lb_TRUEP_Z",&Lb_PZ_rec_JpsiSigma);

	combTree_Sigma->SetBranchAddress("L_TRUEP_E",&L_PE_rec_JpsiSigma);
	combTree_Sigma->SetBranchAddress("L_TRUEP_X",&L_PX_rec_JpsiSigma);
	combTree_Sigma->SetBranchAddress("L_TRUEP_Y",&L_PY_rec_JpsiSigma);
	combTree_Sigma->SetBranchAddress("L_TRUEP_Z",&L_PZ_rec_JpsiSigma);

	combTree_Sigma->SetBranchAddress("wt_tau",&tauWt_rec_JpsiSigma);
	if(run == 1)
	{
		combTree_Sigma->SetBranchAddress("gb_wts_new",&gbWt_rec_JpsiSigma);
	}
	else if(run == 2)
	{
		combTree_Sigma->SetBranchAddress("gb_wts",&gbWt_rec_JpsiSigma);
	}
	cout<<"*****JpsiLambda generated MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_gen_JpsiLambda; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn_gen_JpsiLambda->GetEntry(i);

		Lb_gen_JpsiLambda.SetPxPyPzE(Lb_PX_gen_JpsiLambda,Lb_PY_gen_JpsiLambda,Lb_PZ_gen_JpsiLambda,Lb_PE_gen_JpsiLambda);
		L_gen_JpsiLambda.SetPxPyPzE(L_PX_gen_JpsiLambda,L_PY_gen_JpsiLambda,L_PZ_gen_JpsiLambda,L_PE_gen_JpsiLambda);

		//Boost Lambda momentum to Lb rest frame
		L_gen_JpsiLambda.Boost(-1*Lb_PX_gen_JpsiLambda/Lb_PE_gen_JpsiLambda,-1*Lb_PY_gen_JpsiLambda/Lb_PE_gen_JpsiLambda,-1*Lb_PZ_gen_JpsiLambda/Lb_PE_gen_JpsiLambda);

		n = beam.Cross(Lb_gen_JpsiLambda.Vect());
		theta_JpsiLambda_gen  = (L_gen_JpsiLambda.Vect()).Angle(n);

		Lb_genHist_JpsiLambda->Fill(cos(theta_JpsiLambda_gen));
		Lb_genHist_wt_JpsiLambda->Fill(cos(theta_JpsiLambda_gen),tauWt_gen_JpsiLambda*gbWt_gen_JpsiLambda);
	}
	cout<<"*****JpsiLambda reconstructed MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_gen_JpsiSigma; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn_gen_JpsiSigma->GetEntry(i);

		Lb_gen_JpsiSigma.SetPxPyPzE(Lb_PX_gen_JpsiSigma,Lb_PY_gen_JpsiSigma,Lb_PZ_gen_JpsiSigma,Lb_PE_gen_JpsiSigma);
		L_gen_JpsiSigma.SetPxPyPzE(L_PX_gen_JpsiSigma,L_PY_gen_JpsiSigma,L_PZ_gen_JpsiSigma,L_PE_gen_JpsiSigma);

		//Boost Lambda momentum to Lb rest frame
		L_gen_JpsiSigma.Boost(-1*Lb_PX_gen_JpsiSigma/Lb_PE_gen_JpsiSigma,-1*Lb_PY_gen_JpsiSigma/Lb_PE_gen_JpsiSigma,-1*Lb_PZ_gen_JpsiSigma/Lb_PE_gen_JpsiSigma);

		n = beam.Cross(Lb_gen_JpsiSigma.Vect());
		theta_JpsiSigma_gen  = (L_gen_JpsiSigma.Vect()).Angle(n);

		Lb_genHist_JpsiSigma->Fill(cos(theta_JpsiSigma_gen));
		Lb_genHist_wt_JpsiSigma->Fill(cos(theta_JpsiSigma_gen),tauWt_gen_JpsiSigma*gbWt_gen_JpsiSigma);
	}
	cout<<"*****JpsiSigma generated MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_rec_JpsiLambda; i++)
	{
		if(i%1000 == 0)
		{
			cout<<i<<endl;
		}
		combTree_Lambda->GetEntry(i);

		Lb_rec_JpsiLambda.SetPxPyPzE(Lb_PX_rec_JpsiLambda,Lb_PY_rec_JpsiLambda,Lb_PZ_rec_JpsiLambda,Lb_PE_rec_JpsiLambda);
		L_rec_JpsiLambda.SetPxPyPzE(L_PX_rec_JpsiLambda,L_PY_rec_JpsiLambda,L_PZ_rec_JpsiLambda,L_PE_rec_JpsiLambda);

		//Boost Lambda momentum to Lb rest frame
		L_rec_JpsiLambda.Boost(-1*Lb_PX_rec_JpsiLambda/Lb_PE_rec_JpsiLambda,-1*Lb_PY_rec_JpsiLambda/Lb_PE_rec_JpsiLambda,-1*Lb_PZ_rec_JpsiLambda/Lb_PE_rec_JpsiLambda);

		n = beam.Cross(Lb_rec_JpsiLambda.Vect());
		theta_JpsiLambda_rec  = (L_rec_JpsiLambda.Vect()).Angle(n);

		Lb_recHist_JpsiLambda->Fill(cos(theta_JpsiLambda_rec));
		Lb_recHist_wt_JpsiLambda->Fill(cos(theta_JpsiLambda_rec),tauWt_rec_JpsiLambda*gbWt_rec_JpsiLambda);
	}

	cout<<"*****JpsiSigma reconstructed MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_rec_JpsiSigma; i++)
	{
		if(i%1000 == 0)
		{
			cout<<i<<endl;
		}
		combTree_Sigma->GetEntry(i);

		Lb_rec_JpsiSigma.SetPxPyPzE(Lb_PX_rec_JpsiSigma,Lb_PY_rec_JpsiSigma,Lb_PZ_rec_JpsiSigma,Lb_PE_rec_JpsiSigma);
		L_rec_JpsiSigma.SetPxPyPzE(L_PX_rec_JpsiSigma,L_PY_rec_JpsiSigma,L_PZ_rec_JpsiSigma,L_PE_rec_JpsiSigma);

		//Boost Lambda momentum to Lb rest frame
		L_rec_JpsiSigma.Boost(-1*Lb_PX_rec_JpsiSigma/Lb_PE_rec_JpsiSigma,-1*Lb_PY_rec_JpsiSigma/Lb_PE_rec_JpsiSigma,-1*Lb_PZ_rec_JpsiSigma/Lb_PE_rec_JpsiSigma);

		n = beam.Cross(Lb_rec_JpsiSigma.Vect());
		theta_JpsiSigma_rec  = (L_rec_JpsiSigma.Vect()).Angle(n);

		Lb_recHist_JpsiSigma->Fill(cos(theta_JpsiSigma_rec));
		Lb_recHist_wt_JpsiSigma->Fill(cos(theta_JpsiSigma_rec),tauWt_rec_JpsiSigma*gbWt_rec_JpsiSigma);
	}
	Lb_genHist_JpsiLambda->Sumw2();
	Lb_genHist_wt_JpsiLambda->Sumw2();
	Lb_genHist_JpsiSigma->Sumw2();
	Lb_genHist_wt_JpsiSigma->Sumw2();
	Lb_recHist_JpsiLambda->Sumw2();
	Lb_recHist_wt_JpsiLambda->Sumw2();
	Lb_recHist_JpsiSigma->Sumw2();
	Lb_recHist_wt_JpsiSigma->Sumw2();

	TH1D *eff_JpsiLambda = (TH1D*)Lb_recHist_JpsiLambda->Clone();
	eff_JpsiLambda->SetNameTitle("eff_JpsiLambda","eff_JpsiLambda");
	eff_JpsiLambda->Divide(Lb_genHist_JpsiLambda);

	TH1D *eff_JpsiLambda_wt = (TH1D*)Lb_recHist_wt_JpsiLambda->Clone();
	eff_JpsiLambda_wt->SetNameTitle("eff_JpsiLambda_wt","eff_JpsiLambda_wt");
	eff_JpsiLambda_wt->Divide(Lb_genHist_wt_JpsiLambda);

	TH1D *eff_JpsiSigma = (TH1D*)Lb_recHist_JpsiSigma->Clone();
	eff_JpsiSigma->SetNameTitle("eff_JpsiSigma","eff_JpsiSigma");
	eff_JpsiSigma->Divide(Lb_genHist_JpsiSigma);

	TH1D *eff_JpsiSigma_wt = (TH1D*)Lb_recHist_wt_JpsiSigma->Clone();
	eff_JpsiSigma_wt->SetNameTitle("eff_JpsiSigma_wt","eff_JpsiSigma_wt");
	eff_JpsiSigma_wt->Divide(Lb_genHist_wt_JpsiSigma);

	TH1D *ratio = (TH1D*)eff_JpsiLambda->Clone();
	ratio->SetNameTitle("ratio","ratio");
	ratio->Divide(eff_JpsiSigma);

	TH1D *ratio_wt = (TH1D*)eff_JpsiLambda_wt->Clone();
	ratio_wt->SetNameTitle("ratio_wt","ratio_wt");
	ratio_wt->Divide(eff_JpsiSigma_wt);

	Lb_genHist_JpsiLambda->Draw();
	new TCanvas();
	Lb_genHist_wt_JpsiLambda->Draw();
	new TCanvas();
	Lb_genHist_JpsiSigma->Draw();
	new TCanvas();
	Lb_genHist_wt_JpsiSigma->Draw();
	new TCanvas();
	Lb_recHist_JpsiLambda->Draw();
	new TCanvas();
	Lb_recHist_wt_JpsiLambda->Draw();
	new TCanvas();
	Lb_recHist_JpsiSigma->Draw();
	new TCanvas();
	Lb_recHist_wt_JpsiSigma->Draw();
	new TCanvas();

	eff_JpsiLambda->Draw();
	new TCanvas();
	eff_JpsiLambda_wt->Draw();
	new TCanvas();

	eff_JpsiSigma->Draw();
	new TCanvas();
	eff_JpsiSigma_wt->Draw();
	new TCanvas();

	ratio->Draw();
	new TCanvas();
	ratio_wt->Draw();
}
