#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLine.h"
#include "TF1.h"
#include <iostream>

using namespace std;
using namespace TMath;

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

	Double_t Jpsi_PX_gen_JpsiLambda = 0.0;
	Double_t Jpsi_PY_gen_JpsiLambda = 0.0;
	Double_t Jpsi_PZ_gen_JpsiLambda = 0.0;
	Double_t Jpsi_PE_gen_JpsiLambda = 0.0;

	Double_t L_PX_gen_JpsiLambda = 0.0;
	Double_t L_PY_gen_JpsiLambda = 0.0;
	Double_t L_PZ_gen_JpsiLambda = 0.0;
	Double_t L_PE_gen_JpsiLambda = 0.0;

	Double_t p_PX_gen_JpsiLambda = 0.0;
	Double_t p_PY_gen_JpsiLambda = 0.0;
	Double_t p_PZ_gen_JpsiLambda = 0.0;
	Double_t p_PE_gen_JpsiLambda = 0.0;

	Double_t mu_PX_gen_JpsiLambda = 0.0;
	Double_t mu_PY_gen_JpsiLambda = 0.0;
	Double_t mu_PZ_gen_JpsiLambda = 0.0;
	Double_t mu_PE_gen_JpsiLambda = 0.0;

	Double_t gbWt_gen_JpsiLambda = 0.0;
	Float_t tauWt_gen_JpsiLambda = 0.0;

	Int_t nEntries_gen_JpsiLambda = 0;

	Double_t Lb_PX_rec_JpsiLambda = 0.0;
	Double_t Lb_PY_rec_JpsiLambda = 0.0;
	Double_t Lb_PZ_rec_JpsiLambda = 0.0;
	Double_t Lb_PE_rec_JpsiLambda = 0.0;

	Double_t Jpsi_PX_rec_JpsiLambda = 0.0;
	Double_t Jpsi_PY_rec_JpsiLambda = 0.0;
	Double_t Jpsi_PZ_rec_JpsiLambda = 0.0;
	Double_t Jpsi_PE_rec_JpsiLambda = 0.0;

	Double_t L_PX_rec_JpsiLambda = 0.0;
	Double_t L_PY_rec_JpsiLambda = 0.0;
	Double_t L_PZ_rec_JpsiLambda = 0.0;
	Double_t L_PE_rec_JpsiLambda = 0.0;

	Double_t p_PX_rec_JpsiLambda = 0.0;
	Double_t p_PY_rec_JpsiLambda = 0.0;
	Double_t p_PZ_rec_JpsiLambda = 0.0;
	Double_t p_PE_rec_JpsiLambda = 0.0;

	Double_t mu_PX_rec_JpsiLambda = 0.0;
	Double_t mu_PY_rec_JpsiLambda = 0.0;
	Double_t mu_PZ_rec_JpsiLambda = 0.0;
	Double_t mu_PE_rec_JpsiLambda = 0.0;

	Double_t gbWt_rec_JpsiLambda = 0.0;
	Float_t tauWt_rec_JpsiLambda = 0.0;

	Int_t nEntries_rec_JpsiLambda = 0;

	Double_t Lb_PX_gen_JpsiSigma = 0.0;
	Double_t Lb_PY_gen_JpsiSigma = 0.0;
	Double_t Lb_PZ_gen_JpsiSigma = 0.0;
	Double_t Lb_PE_gen_JpsiSigma = 0.0;

	Double_t Jpsi_PX_gen_JpsiSigma = 0.0;
	Double_t Jpsi_PY_gen_JpsiSigma = 0.0;
	Double_t Jpsi_PZ_gen_JpsiSigma = 0.0;
	Double_t Jpsi_PE_gen_JpsiSigma = 0.0;

	Double_t L_PX_gen_JpsiSigma = 0.0;
	Double_t L_PY_gen_JpsiSigma = 0.0;
	Double_t L_PZ_gen_JpsiSigma = 0.0;
	Double_t L_PE_gen_JpsiSigma = 0.0;

	Double_t p_PX_gen_JpsiSigma = 0.0;
	Double_t p_PY_gen_JpsiSigma = 0.0;
	Double_t p_PZ_gen_JpsiSigma = 0.0;
	Double_t p_PE_gen_JpsiSigma = 0.0;

	Double_t mu_PX_gen_JpsiSigma = 0.0;
	Double_t mu_PY_gen_JpsiSigma = 0.0;
	Double_t mu_PZ_gen_JpsiSigma = 0.0;
	Double_t mu_PE_gen_JpsiSigma = 0.0;

	Double_t gbWt_gen_JpsiSigma = 0.0;
	Float_t tauWt_gen_JpsiSigma = 0.0;

	Int_t nEntries_gen_JpsiSigma = 0;

	Double_t Lb_PX_rec_JpsiSigma = 0.0;
	Double_t Lb_PY_rec_JpsiSigma = 0.0;
	Double_t Lb_PZ_rec_JpsiSigma = 0.0;
	Double_t Lb_PE_rec_JpsiSigma = 0.0;

	Double_t Jpsi_PX_rec_JpsiSigma = 0.0;
	Double_t Jpsi_PY_rec_JpsiSigma = 0.0;
	Double_t Jpsi_PZ_rec_JpsiSigma = 0.0;
	Double_t Jpsi_PE_rec_JpsiSigma = 0.0;

	Double_t L_PX_rec_JpsiSigma = 0.0;
	Double_t L_PY_rec_JpsiSigma = 0.0;
	Double_t L_PZ_rec_JpsiSigma = 0.0;
	Double_t L_PE_rec_JpsiSigma = 0.0;

	Double_t p_PX_rec_JpsiSigma = 0.0;
	Double_t p_PY_rec_JpsiSigma = 0.0;
	Double_t p_PZ_rec_JpsiSigma = 0.0;
	Double_t p_PE_rec_JpsiSigma = 0.0;

	Double_t mu_PX_rec_JpsiSigma = 0.0;
	Double_t mu_PY_rec_JpsiSigma = 0.0;
	Double_t mu_PZ_rec_JpsiSigma = 0.0;
	Double_t mu_PE_rec_JpsiSigma = 0.0;

	Double_t gbWt_rec_JpsiSigma = 0.0;
	Float_t tauWt_rec_JpsiSigma = 0.0;

	Int_t nEntries_rec_JpsiSigma = 0;

	Double_t chi2 = 0.0;

	TString title_y    = "#frac{#epsilon_{#Lambda_b #rightarrow J/#psi #Lambda}}{#epsilon_{#Lambda_b #rightarrow J/#psi #Sigma}}";
	TString title_y_wt = "Weighted #frac{#epsilon_{#Lambda_b #rightarrow J/#psi #Lambda}}{#epsilon_{#Lambda_b #rightarrow J/#psi #Sigma}}";

	TLorentzVector Lb_gen_JpsiLambda, Lb_rec_JpsiLambda, L_gen_JpsiLambda, L_rec_JpsiLambda;
	TLorentzVector Jpsi_gen_JpsiLambda, Jpsi_rec_JpsiLambda, p_gen_JpsiLambda, p_rec_JpsiLambda;
	TLorentzVector mu_gen_JpsiLambda, mu_rec_JpsiLambda;

	TLorentzVector Lb_gen_JpsiSigma, Lb_rec_JpsiSigma, L_gen_JpsiSigma, L_rec_JpsiSigma;
	TLorentzVector Jpsi_gen_JpsiSigma, Jpsi_rec_JpsiSigma, p_gen_JpsiSigma, p_rec_JpsiSigma;
	TLorentzVector mu_gen_JpsiSigma, mu_rec_JpsiSigma;

	TVector3 n;//p_Lb X p_beam
	TVector3 z1;//parallel to Lambda momentum in the Lb rest frame
	TVector3 y1;//n X (Lambda momentum in the Lb rest frame)
	TVector3 z2;//parallel to Jpsi momentum in the Lb rest frame
	TVector3 y2;//n X (Lambda momentum in the Lb rest frame)

	TVector3 z1_prime, z2_prime;

	TVector3 beam;
	beam.SetXYZ(0.0,0.0,1.0); //beam is defined to be along z axis

	Double_t theta, theta1, phi1, theta2, phi2;

	TH1D *theta_gen_JpsiLambda = new TH1D("theta_gen_JpsiLambda","theta_gen_JpsiLambda",50,-1,1);
	TH1D *theta_gen_wt_JpsiLambda = new TH1D("theta_gen_wt_JpsiLambda","theta_gen_wt_JpsiLambda",50,-1,1);

	TH1D *theta_gen_JpsiSigma = new TH1D("theta_gen_JpsiSigma","theta_gen_JpsiSigma",50,-1,1);
	TH1D *theta_gen_wt_JpsiSigma = new TH1D("theta_gen_wt_JpsiSigma","theta_gen_wt_JpsiSigma",50,-1,1);

	TH1D *theta_rec_JpsiLambda = new TH1D("theta_rec_JpsiLambda","theta_rec_JpsiLambda",50,-1,1);
	TH1D *theta_rec_wt_JpsiLambda = new TH1D("theta_rec_wt_JpsiLambda","theta_rec_wt_JpsiLambda",50,-1,1);

	TH1D *theta_rec_JpsiSigma = new TH1D("theta_rec_JpsiSigma","theta_rec_JpsiSigma",50,-1,1);
	TH1D *theta_rec_wt_JpsiSigma = new TH1D("theta_rec_wt_JpsiSigma","theta_rec_wt_JpsiSigma",50,-1,1);

	TH1D *theta1_gen_JpsiLambda = new TH1D("theta1_gen_JpsiLambda","theta1_gen_JpsiLambda",50,-1,1);
	TH1D *theta1_gen_wt_JpsiLambda = new TH1D("theta1_gen_wt_JpsiLambda","theta1_gen_wt_JpsiLambda",50,-1,1);

	TH1D *theta1_gen_JpsiSigma = new TH1D("theta1_gen_JpsiSigma","theta1_gen_JpsiSigma",50,-1,1);
	TH1D *theta1_gen_wt_JpsiSigma = new TH1D("theta1_gen_wt_JpsiSigma","theta1_gen_wt_JpsiSigma",50,-1,1);

	TH1D *theta1_rec_JpsiLambda = new TH1D("theta1_rec_JpsiLambda","theta1_rec_JpsiLambda",50,-1,1);
	TH1D *theta1_rec_wt_JpsiLambda = new TH1D("theta1_rec_wt_JpsiLambda","theta1_rec_wt_JpsiLambda",50,-1,1);

	TH1D *theta1_rec_JpsiSigma = new TH1D("theta1_rec_JpsiSigma","theta1_rec_JpsiSigma",50,-1,1);
	TH1D *theta1_rec_wt_JpsiSigma = new TH1D("theta1_rec_wt_JpsiSigma","theta1_rec_wt_JpsiSigma",50,-1,1);

	TH1D *theta2_gen_JpsiLambda = new TH1D("theta2_gen_JpsiLambda","theta2_gen_JpsiLambda",50,-1,1);
	TH1D *theta2_gen_wt_JpsiLambda = new TH1D("theta2_gen_wt_JpsiLambda","theta2_gen_wt_JpsiLambda",50,-1,1);

	TH1D *theta2_gen_JpsiSigma = new TH1D("theta2_gen_JpsiSigma","theta2_gen_JpsiSigma",50,-1,1);
	TH1D *theta2_gen_wt_JpsiSigma = new TH1D("theta2_gen_wt_JpsiSigma","theta2_gen_wt_JpsiSigma",50,-1,1);

	TH1D *theta2_rec_JpsiLambda = new TH1D("theta2_rec_JpsiLambda","theta2_rec_JpsiLambda",50,-1,1);
	TH1D *theta2_rec_wt_JpsiLambda = new TH1D("theta2_rec_wt_JpsiLambda","theta2_rec_wt_JpsiLambda",50,-1,1);

	TH1D *theta2_rec_JpsiSigma = new TH1D("theta2_rec_JpsiSigma","theta2_rec_JpsiSigma",50,-1,1);
	TH1D *theta2_rec_wt_JpsiSigma = new TH1D("theta2_rec_wt_JpsiSigma","theta2_rec_wt_JpsiSigma",50,-1,1);

	TH1D *phi1_gen_JpsiLambda = new TH1D("phi1_gen_JpsiLambda","phi1_gen_JpsiLambda",50,0,1);
	TH1D *phi1_gen_wt_JpsiLambda = new TH1D("phi1_gen_wt_JpsiLambda","phi1_gen_wt_JpsiLambda",50,0,1);

	TH1D *phi1_gen_JpsiSigma = new TH1D("phi1_gen_JpsiSigma","phi1_gen_JpsiSigma",50,0,1);
	TH1D *phi1_gen_wt_JpsiSigma = new TH1D("phi1_gen_wt_JpsiSigma","phi1_gen_wt_JpsiSigma",50,0,1);

	TH1D *phi1_rec_JpsiLambda = new TH1D("phi1_rec_JpsiLambda","phi1_rec_JpsiLambda",50,0,1);
	TH1D *phi1_rec_wt_JpsiLambda = new TH1D("phi1_rec_wt_JpsiLambda","phi1_rec_wt_JpsiLambda",50,0,1);

	TH1D *phi1_rec_JpsiSigma = new TH1D("phi1_rec_JpsiSigma","phi1_rec_JpsiSigma",50,0,1);
	TH1D *phi1_rec_wt_JpsiSigma = new TH1D("phi1_rec_wt_JpsiSigma","phi1_rec_wt_JpsiSigma",50,0,1);

	TH1D *phi2_gen_JpsiLambda = new TH1D("phi2_gen_JpsiLambda","phi2_gen_JpsiLambda",50,0,1);
	TH1D *phi2_gen_wt_JpsiLambda = new TH1D("phi2_gen_wt_JpsiLambda","phi2_gen_wt_JpsiLambda",50,0,1);

	TH1D *phi2_gen_JpsiSigma = new TH1D("phi2_gen_JpsiSigma","phi2_gen_JpsiSigma",50,0,1);
	TH1D *phi2_gen_wt_JpsiSigma = new TH1D("phi2_gen_wt_JpsiSigma","phi2_gen_wt_JpsiSigma",50,0,1);

	TH1D *phi2_rec_JpsiLambda = new TH1D("phi2_rec_JpsiLambda","phi2_rec_JpsiLambda",50,0,1);
	TH1D *phi2_rec_wt_JpsiLambda = new TH1D("phi2_rec_wt_JpsiLambda","phi2_rec_wt_JpsiLambda",50,0,1);

	TH1D *phi2_rec_JpsiSigma = new TH1D("phi2_rec_JpsiSigma","phi2_rec_JpsiSigma",50,0,1);
	TH1D *phi2_rec_wt_JpsiSigma = new TH1D("phi2_rec_wt_JpsiSigma","phi2_rec_wt_JpsiSigma",50,0,1);

	TLine *line1 = new TLine(-1,1,1,1);
	TLine *line2 = new TLine(0,1,1,1);

	TF1 *unif = new TF1("unif","1",-1,1);
	TF1 *unif1 = new TF1("unif1","1",0,1);

	line1->SetLineColor(kRed);
	line2->SetLineColor(kRed);

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

	treeIn_gen_JpsiLambda->SetBranchAddress("J_psi_1S_TRUEP_E",&Jpsi_PE_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("J_psi_1S_TRUEP_X",&Jpsi_PX_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("J_psi_1S_TRUEP_Y",&Jpsi_PY_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("J_psi_1S_TRUEP_Z",&Jpsi_PZ_gen_JpsiLambda);

	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_E",&L_PE_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_X",&L_PX_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_Y",&L_PY_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("Lambda0_TRUEP_Z",&L_PZ_gen_JpsiLambda);

	treeIn_gen_JpsiLambda->SetBranchAddress("pplus_TRUEP_E",&p_PE_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("pplus_TRUEP_X",&p_PX_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("pplus_TRUEP_Y",&p_PY_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("pplus_TRUEP_Z",&p_PZ_gen_JpsiLambda);

	treeIn_gen_JpsiLambda->SetBranchAddress("muplus_TRUEP_E",&mu_PE_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("muplus_TRUEP_X",&mu_PX_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("muplus_TRUEP_Y",&mu_PY_gen_JpsiLambda);
	treeIn_gen_JpsiLambda->SetBranchAddress("muplus_TRUEP_Z",&mu_PZ_gen_JpsiLambda);

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

	treeIn_gen_JpsiSigma->SetBranchAddress("J_psi_1S_TRUEP_E",&Jpsi_PE_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("J_psi_1S_TRUEP_X",&Jpsi_PX_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("J_psi_1S_TRUEP_Y",&Jpsi_PY_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("J_psi_1S_TRUEP_Z",&Jpsi_PZ_gen_JpsiSigma);

	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_E",&L_PE_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_X",&L_PX_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_Y",&L_PY_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("Lambda0_TRUEP_Z",&L_PZ_gen_JpsiSigma);

	treeIn_gen_JpsiSigma->SetBranchAddress("pplus_TRUEP_E",&p_PE_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("pplus_TRUEP_X",&p_PX_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("pplus_TRUEP_Y",&p_PY_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("pplus_TRUEP_Z",&p_PZ_gen_JpsiSigma);

	treeIn_gen_JpsiSigma->SetBranchAddress("muplus_TRUEP_E",&mu_PE_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("muplus_TRUEP_X",&mu_PX_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("muplus_TRUEP_Y",&mu_PY_gen_JpsiSigma);
	treeIn_gen_JpsiSigma->SetBranchAddress("muplus_TRUEP_Z",&mu_PZ_gen_JpsiSigma);

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

	TTree *treeIn_rec_JpsiLambda = TTree::MergeTrees(list_Lambda);
	treeIn_rec_JpsiLambda->SetName("treeIn_rec_JpsiLambda");

	nEntries_rec_JpsiLambda = treeIn_rec_JpsiLambda->GetEntries();

	TFile *tempFile2 = new TFile("tempFile2.root","RECREATE");

	TTree *treeIn_nonZero_Sigma_cut = (TTree*)treeIn_nonZero_JpsiSigma->CopyTree(Form("BDT2 > %f",bdtCut_nonZero));
	TTree *treeIn_Zero_Sigma_cut    = (TTree*)treeIn_Zero_JpsiSigma->CopyTree(Form("BDT2 > %f",bdtCut_Zero));

	TList *list_Sigma = new TList;
	list_Sigma->Add(treeIn_nonZero_Sigma_cut);
	list_Sigma->Add(treeIn_Zero_Sigma_cut);

	TTree *treeIn_rec_JpsiSigma = TTree::MergeTrees(list_Sigma);
	treeIn_rec_JpsiSigma->SetName("treeIn_rec_JpsiSigma");

	nEntries_rec_JpsiSigma = treeIn_rec_JpsiSigma->GetEntries();

	treeIn_rec_JpsiLambda->SetBranchAddress("Lb_TRUEP_E",&Lb_PE_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("Lb_TRUEP_X",&Lb_PX_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("Lb_TRUEP_Y",&Lb_PY_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("Lb_TRUEP_Z",&Lb_PZ_rec_JpsiLambda);

	treeIn_rec_JpsiLambda->SetBranchAddress("Jpsi_TRUEP_E",&Jpsi_PE_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("Jpsi_TRUEP_X",&Jpsi_PX_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("Jpsi_TRUEP_Y",&Jpsi_PY_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("Jpsi_TRUEP_Z",&Jpsi_PZ_rec_JpsiLambda);

	treeIn_rec_JpsiLambda->SetBranchAddress("L_TRUEP_E",&L_PE_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("L_TRUEP_X",&L_PX_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("L_TRUEP_Y",&L_PY_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("L_TRUEP_Z",&L_PZ_rec_JpsiLambda);

	treeIn_rec_JpsiLambda->SetBranchAddress("p_TRUEP_E",&p_PE_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("p_TRUEP_X",&p_PX_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("p_TRUEP_Y",&p_PY_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("p_TRUEP_Z",&p_PZ_rec_JpsiLambda);

	treeIn_rec_JpsiLambda->SetBranchAddress("muplus_TRUEP_E",&mu_PE_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("muplus_TRUEP_X",&mu_PX_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("muplus_TRUEP_Y",&mu_PY_rec_JpsiLambda);
	treeIn_rec_JpsiLambda->SetBranchAddress("muplus_TRUEP_Z",&mu_PZ_rec_JpsiLambda);

	treeIn_rec_JpsiLambda->SetBranchAddress("wt_tau",&tauWt_rec_JpsiLambda);
	if(run == 1)
	{
		treeIn_rec_JpsiLambda->SetBranchAddress("gb_wts_new",&gbWt_rec_JpsiLambda);
	}
	else if(run == 2)
	{
		treeIn_rec_JpsiLambda->SetBranchAddress("gb_wts",&gbWt_rec_JpsiLambda);
	}

	treeIn_rec_JpsiSigma->SetBranchAddress("Lb_TRUEP_E",&Lb_PE_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("Lb_TRUEP_X",&Lb_PX_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("Lb_TRUEP_Y",&Lb_PY_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("Lb_TRUEP_Z",&Lb_PZ_rec_JpsiSigma);

	treeIn_rec_JpsiSigma->SetBranchAddress("Jpsi_TRUEP_E",&Jpsi_PE_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("Jpsi_TRUEP_X",&Jpsi_PX_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("Jpsi_TRUEP_Y",&Jpsi_PY_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("Jpsi_TRUEP_Z",&Jpsi_PZ_rec_JpsiSigma);

	treeIn_rec_JpsiSigma->SetBranchAddress("L_TRUEP_E",&L_PE_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("L_TRUEP_X",&L_PX_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("L_TRUEP_Y",&L_PY_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("L_TRUEP_Z",&L_PZ_rec_JpsiSigma);

	treeIn_rec_JpsiSigma->SetBranchAddress("p_TRUEP_E",&p_PE_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("p_TRUEP_X",&p_PX_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("p_TRUEP_Y",&p_PY_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("p_TRUEP_Z",&p_PZ_rec_JpsiSigma);

	treeIn_rec_JpsiSigma->SetBranchAddress("muplus_TRUEP_E",&mu_PE_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("muplus_TRUEP_X",&mu_PX_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("muplus_TRUEP_Y",&mu_PY_rec_JpsiSigma);
	treeIn_rec_JpsiSigma->SetBranchAddress("muplus_TRUEP_Z",&mu_PZ_rec_JpsiSigma);

	treeIn_rec_JpsiSigma->SetBranchAddress("wt_tau",&tauWt_rec_JpsiSigma);
	if(run == 1)
	{
		treeIn_rec_JpsiSigma->SetBranchAddress("gb_wts_new",&gbWt_rec_JpsiSigma);
	}
	else if(run == 2)
	{
		treeIn_rec_JpsiSigma->SetBranchAddress("gb_wts",&gbWt_rec_JpsiSigma);
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
		Jpsi_gen_JpsiLambda.SetPxPyPzE(Jpsi_PX_gen_JpsiLambda,Jpsi_PY_gen_JpsiLambda,Jpsi_PZ_gen_JpsiLambda,Jpsi_PE_gen_JpsiLambda);
		L_gen_JpsiLambda.SetPxPyPzE(L_PX_gen_JpsiLambda,L_PY_gen_JpsiLambda,L_PZ_gen_JpsiLambda,L_PE_gen_JpsiLambda);
		p_gen_JpsiLambda.SetPxPyPzE(p_PX_gen_JpsiLambda,p_PY_gen_JpsiLambda,p_PZ_gen_JpsiLambda,p_PE_gen_JpsiLambda);
		mu_gen_JpsiLambda.SetPxPyPzE(mu_PX_gen_JpsiLambda,mu_PY_gen_JpsiLambda,mu_PZ_gen_JpsiLambda,mu_PE_gen_JpsiLambda);

		//Boost Lambda momentum to Lb rest frame
		L_gen_JpsiLambda.Boost(-1*Lb_PX_gen_JpsiLambda/Lb_PE_gen_JpsiLambda,-1*Lb_PY_gen_JpsiLambda/Lb_PE_gen_JpsiLambda,-1*Lb_PZ_gen_JpsiLambda/Lb_PE_gen_JpsiLambda);
		//Boost J/psi momentum to Lb rest frame
		Jpsi_gen_JpsiLambda.Boost(-1*Lb_PX_gen_JpsiLambda/Lb_PE_gen_JpsiLambda,-1*Lb_PY_gen_JpsiLambda/Lb_PE_gen_JpsiLambda,-1*Lb_PZ_gen_JpsiLambda/Lb_PE_gen_JpsiLambda);

		//Boost muon momentum to J/psi rest frame
		mu_gen_JpsiLambda.Boost(-1*Jpsi_PX_gen_JpsiLambda/Jpsi_PE_gen_JpsiLambda,-1*Jpsi_PY_gen_JpsiLambda/Jpsi_PE_gen_JpsiLambda,-1*Jpsi_PZ_gen_JpsiLambda/Jpsi_PE_gen_JpsiLambda);
		//Boost proton momentum to Lambda rest frame
		p_gen_JpsiLambda.Boost(-1*L_PX_gen_JpsiLambda/L_PE_gen_JpsiLambda,-1*L_PY_gen_JpsiLambda/L_PE_gen_JpsiLambda,-1*L_PZ_gen_JpsiLambda/L_PE_gen_JpsiLambda);

		n = (beam.Cross(Lb_gen_JpsiLambda.Vect())).Unit();
		z1 = (L_gen_JpsiLambda.Vect()).Unit();
		z2 = (Jpsi_gen_JpsiLambda.Vect()).Unit();
		y1 = (n.Cross(z1));
		y2 = (n.Cross(z2));

		theta  = (L_gen_JpsiLambda.Vect()).Angle(n);
		theta1 = (p_gen_JpsiLambda.Vect()).Angle(z1);
		theta2 = (mu_gen_JpsiLambda.Vect()).Angle(z2);

		z1_prime = z1;
		z1_prime *= ((p_gen_JpsiLambda.Vect()).Dot(z1));

		z2_prime = z2;
		z2_prime *= ((mu_gen_JpsiLambda.Vect()).Dot(z2));

		phi1   = PiOver2() - (p_gen_JpsiLambda.Vect() - z1_prime).Angle(y1);
		phi2   = PiOver2() - (mu_gen_JpsiLambda.Vect() - z2_prime).Angle(y2);

		theta_gen_JpsiLambda->Fill(cos(theta));
		theta_gen_wt_JpsiLambda->Fill(cos(theta),tauWt_gen_JpsiLambda*gbWt_gen_JpsiLambda);

		theta1_gen_JpsiLambda->Fill(cos(theta1));
		theta1_gen_wt_JpsiLambda->Fill(cos(theta1),tauWt_gen_JpsiLambda*gbWt_gen_JpsiLambda);

		theta2_gen_JpsiLambda->Fill(cos(theta2));
		theta2_gen_wt_JpsiLambda->Fill(cos(theta2),tauWt_gen_JpsiLambda*gbWt_gen_JpsiLambda);

		phi1_gen_JpsiLambda->Fill(cos(phi1));
		phi1_gen_wt_JpsiLambda->Fill(cos(phi1),tauWt_gen_JpsiLambda*gbWt_gen_JpsiLambda);

		phi2_gen_JpsiLambda->Fill(cos(phi2));
		phi2_gen_wt_JpsiLambda->Fill(cos(phi2),tauWt_gen_JpsiLambda*gbWt_gen_JpsiLambda);
	}
	cout<<"*****JpsiSigma generated MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_gen_JpsiSigma; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn_gen_JpsiSigma->GetEntry(i);

		Lb_gen_JpsiSigma.SetPxPyPzE(Lb_PX_gen_JpsiSigma,Lb_PY_gen_JpsiSigma,Lb_PZ_gen_JpsiSigma,Lb_PE_gen_JpsiSigma);
		Jpsi_gen_JpsiSigma.SetPxPyPzE(Jpsi_PX_gen_JpsiSigma,Jpsi_PY_gen_JpsiSigma,Jpsi_PZ_gen_JpsiSigma,Jpsi_PE_gen_JpsiSigma);
		L_gen_JpsiSigma.SetPxPyPzE(L_PX_gen_JpsiSigma,L_PY_gen_JpsiSigma,L_PZ_gen_JpsiSigma,L_PE_gen_JpsiSigma);
		p_gen_JpsiSigma.SetPxPyPzE(p_PX_gen_JpsiSigma,p_PY_gen_JpsiSigma,p_PZ_gen_JpsiSigma,p_PE_gen_JpsiSigma);
		mu_gen_JpsiSigma.SetPxPyPzE(mu_PX_gen_JpsiSigma,mu_PY_gen_JpsiSigma,mu_PZ_gen_JpsiSigma,mu_PE_gen_JpsiSigma);

		//Boost Lambda momentum to Lb rest frame
		L_gen_JpsiSigma.Boost(-1*Lb_PX_gen_JpsiSigma/Lb_PE_gen_JpsiSigma,-1*Lb_PY_gen_JpsiSigma/Lb_PE_gen_JpsiSigma,-1*Lb_PZ_gen_JpsiSigma/Lb_PE_gen_JpsiSigma);
		//Boost J/psi momentum to Lb rest frame
		Jpsi_gen_JpsiSigma.Boost(-1*Lb_PX_gen_JpsiSigma/Lb_PE_gen_JpsiSigma,-1*Lb_PY_gen_JpsiSigma/Lb_PE_gen_JpsiSigma,-1*Lb_PZ_gen_JpsiSigma/Lb_PE_gen_JpsiSigma);

		//Boost muon momentum to J/psi rest frame
		mu_gen_JpsiSigma.Boost(-1*Jpsi_PX_gen_JpsiSigma/Jpsi_PE_gen_JpsiSigma,-1*Jpsi_PY_gen_JpsiSigma/Jpsi_PE_gen_JpsiSigma,-1*Jpsi_PZ_gen_JpsiSigma/Jpsi_PE_gen_JpsiSigma);
		//Boost proton momentum to Lambda rest frame
		p_gen_JpsiSigma.Boost(-1*L_PX_gen_JpsiSigma/L_PE_gen_JpsiSigma,-1*L_PY_gen_JpsiSigma/L_PE_gen_JpsiSigma,-1*L_PZ_gen_JpsiSigma/L_PE_gen_JpsiSigma);

		n = (beam.Cross(Lb_gen_JpsiSigma.Vect())).Unit();
		z1 = (L_gen_JpsiSigma.Vect()).Unit();
		z2 = (Jpsi_gen_JpsiSigma.Vect()).Unit();
		y1 = (n.Cross(z1));
		y2 = (n.Cross(z2));

		theta  = (L_gen_JpsiSigma.Vect()).Angle(n);
		theta1 = (p_gen_JpsiSigma.Vect()).Angle(z1);
		theta2 = (mu_gen_JpsiSigma.Vect()).Angle(z2);

		z1_prime = z1;
		z1_prime *= ((p_gen_JpsiSigma.Vect()).Dot(z1));

		z2_prime = z2;
		z2_prime *= ((mu_gen_JpsiSigma.Vect()).Dot(z2));

		phi1   = PiOver2() - (p_gen_JpsiSigma.Vect() - z1_prime).Angle(y1);
		phi2   = PiOver2() - (mu_gen_JpsiSigma.Vect() - z2_prime).Angle(y2);

		theta_gen_JpsiSigma->Fill(cos(theta));
		theta_gen_wt_JpsiSigma->Fill(cos(theta),tauWt_gen_JpsiSigma*gbWt_gen_JpsiSigma);

		theta1_gen_JpsiSigma->Fill(cos(theta1));
		theta1_gen_wt_JpsiSigma->Fill(cos(theta1),tauWt_gen_JpsiSigma*gbWt_gen_JpsiSigma);

		theta2_gen_JpsiSigma->Fill(cos(theta2));
		theta2_gen_wt_JpsiSigma->Fill(cos(theta2),tauWt_gen_JpsiSigma*gbWt_gen_JpsiSigma);

		phi1_gen_JpsiSigma->Fill(cos(phi1));
		phi1_gen_wt_JpsiSigma->Fill(cos(phi1),tauWt_gen_JpsiSigma*gbWt_gen_JpsiSigma);

		phi2_gen_JpsiSigma->Fill(cos(phi2));
		phi2_gen_wt_JpsiSigma->Fill(cos(phi2),tauWt_gen_JpsiSigma*gbWt_gen_JpsiSigma);
	}
	cout<<"*****JpsiLambda reconstructed MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_rec_JpsiLambda; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn_rec_JpsiLambda->GetEntry(i);

		Lb_rec_JpsiLambda.SetPxPyPzE(Lb_PX_rec_JpsiLambda,Lb_PY_rec_JpsiLambda,Lb_PZ_rec_JpsiLambda,Lb_PE_rec_JpsiLambda);
		Jpsi_rec_JpsiLambda.SetPxPyPzE(Jpsi_PX_rec_JpsiLambda,Jpsi_PY_rec_JpsiLambda,Jpsi_PZ_rec_JpsiLambda,Jpsi_PE_rec_JpsiLambda);
		L_rec_JpsiLambda.SetPxPyPzE(L_PX_rec_JpsiLambda,L_PY_rec_JpsiLambda,L_PZ_rec_JpsiLambda,L_PE_rec_JpsiLambda);
		p_rec_JpsiLambda.SetPxPyPzE(p_PX_rec_JpsiLambda,p_PY_rec_JpsiLambda,p_PZ_rec_JpsiLambda,p_PE_rec_JpsiLambda);
		mu_rec_JpsiLambda.SetPxPyPzE(mu_PX_rec_JpsiLambda,mu_PY_rec_JpsiLambda,mu_PZ_rec_JpsiLambda,mu_PE_rec_JpsiLambda);

		//Boost Lambda momentum to Lb rest frame
		L_rec_JpsiLambda.Boost(-1*Lb_PX_rec_JpsiLambda/Lb_PE_rec_JpsiLambda,-1*Lb_PY_rec_JpsiLambda/Lb_PE_rec_JpsiLambda,-1*Lb_PZ_rec_JpsiLambda/Lb_PE_rec_JpsiLambda);
		//Boost J/psi momentum to Lb rest frame
		Jpsi_rec_JpsiLambda.Boost(-1*Lb_PX_rec_JpsiLambda/Lb_PE_rec_JpsiLambda,-1*Lb_PY_rec_JpsiLambda/Lb_PE_rec_JpsiLambda,-1*Lb_PZ_rec_JpsiLambda/Lb_PE_rec_JpsiLambda);

		//Boost muon momentum to J/psi rest frame
		mu_rec_JpsiLambda.Boost(-1*Jpsi_PX_rec_JpsiLambda/Jpsi_PE_rec_JpsiLambda,-1*Jpsi_PY_rec_JpsiLambda/Jpsi_PE_rec_JpsiLambda,-1*Jpsi_PZ_rec_JpsiLambda/Jpsi_PE_rec_JpsiLambda);
		//Boost proton momentum to Lambda rest frame
		p_rec_JpsiLambda.Boost(-1*L_PX_rec_JpsiLambda/L_PE_rec_JpsiLambda,-1*L_PY_rec_JpsiLambda/L_PE_rec_JpsiLambda,-1*L_PZ_rec_JpsiLambda/L_PE_rec_JpsiLambda);

		n = (beam.Cross(Lb_rec_JpsiLambda.Vect())).Unit();
		z1 = (L_rec_JpsiLambda.Vect()).Unit();
		z2 = (Jpsi_rec_JpsiLambda.Vect()).Unit();
		y1 = (n.Cross(z1));
		y2 = (n.Cross(z2));

		theta  = (L_rec_JpsiLambda.Vect()).Angle(n);
		theta1 = (p_rec_JpsiLambda.Vect()).Angle(z1);
		theta2 = (mu_rec_JpsiLambda.Vect()).Angle(z2);

		z1_prime = z1;
		z1_prime *= ((p_rec_JpsiLambda.Vect()).Dot(z1));

		z2_prime = z2;
		z2_prime *= ((mu_rec_JpsiLambda.Vect()).Dot(z2));

		phi1   = PiOver2() - (p_rec_JpsiLambda.Vect() - z1_prime).Angle(y1);
		phi2   = PiOver2() - (mu_rec_JpsiLambda.Vect() - z2_prime).Angle(y2);

		theta_rec_JpsiLambda->Fill(cos(theta));
		theta_rec_wt_JpsiLambda->Fill(cos(theta),tauWt_rec_JpsiLambda*gbWt_rec_JpsiLambda);

		theta1_rec_JpsiLambda->Fill(cos(theta1));
		theta1_rec_wt_JpsiLambda->Fill(cos(theta1),tauWt_rec_JpsiLambda*gbWt_rec_JpsiLambda);

		theta2_rec_JpsiLambda->Fill(cos(theta2));
		theta2_rec_wt_JpsiLambda->Fill(cos(theta2),tauWt_rec_JpsiLambda*gbWt_rec_JpsiLambda);

		phi1_rec_JpsiLambda->Fill(cos(phi1));
		phi1_rec_wt_JpsiLambda->Fill(cos(phi1),tauWt_rec_JpsiLambda*gbWt_rec_JpsiLambda);

		phi2_rec_JpsiLambda->Fill(cos(phi2));
		phi2_rec_wt_JpsiLambda->Fill(cos(phi2),tauWt_rec_JpsiLambda*gbWt_rec_JpsiLambda);
	}

	cout<<"*****JpsiSigma reconstructed MC*****"<<endl;
	for(Int_t i = 0; i<nEntries_rec_JpsiSigma; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn_rec_JpsiSigma->GetEntry(i);

		Lb_rec_JpsiSigma.SetPxPyPzE(Lb_PX_rec_JpsiSigma,Lb_PY_rec_JpsiSigma,Lb_PZ_rec_JpsiSigma,Lb_PE_rec_JpsiSigma);
		Jpsi_rec_JpsiSigma.SetPxPyPzE(Jpsi_PX_rec_JpsiSigma,Jpsi_PY_rec_JpsiSigma,Jpsi_PZ_rec_JpsiSigma,Jpsi_PE_rec_JpsiSigma);
		L_rec_JpsiSigma.SetPxPyPzE(L_PX_rec_JpsiSigma,L_PY_rec_JpsiSigma,L_PZ_rec_JpsiSigma,L_PE_rec_JpsiSigma);
		p_rec_JpsiSigma.SetPxPyPzE(p_PX_rec_JpsiSigma,p_PY_rec_JpsiSigma,p_PZ_rec_JpsiSigma,p_PE_rec_JpsiSigma);
		mu_rec_JpsiSigma.SetPxPyPzE(mu_PX_rec_JpsiSigma,mu_PY_rec_JpsiSigma,mu_PZ_rec_JpsiSigma,mu_PE_rec_JpsiSigma);

		//Boost Lambda momentum to Lb rest frame
		L_rec_JpsiSigma.Boost(-1*Lb_PX_rec_JpsiSigma/Lb_PE_rec_JpsiSigma,-1*Lb_PY_rec_JpsiSigma/Lb_PE_rec_JpsiSigma,-1*Lb_PZ_rec_JpsiSigma/Lb_PE_rec_JpsiSigma);
		//Boost J/psi momentum to Lb rest frame
		Jpsi_rec_JpsiSigma.Boost(-1*Lb_PX_rec_JpsiSigma/Lb_PE_rec_JpsiSigma,-1*Lb_PY_rec_JpsiSigma/Lb_PE_rec_JpsiSigma,-1*Lb_PZ_rec_JpsiSigma/Lb_PE_rec_JpsiSigma);

		//Boost muon momentum to J/psi rest frame
		mu_rec_JpsiSigma.Boost(-1*Jpsi_PX_rec_JpsiSigma/Jpsi_PE_rec_JpsiSigma,-1*Jpsi_PY_rec_JpsiSigma/Jpsi_PE_rec_JpsiSigma,-1*Jpsi_PZ_rec_JpsiSigma/Jpsi_PE_rec_JpsiSigma);
		//Boost proton momentum to Lambda rest frame
		p_rec_JpsiSigma.Boost(-1*L_PX_rec_JpsiSigma/L_PE_rec_JpsiSigma,-1*L_PY_rec_JpsiSigma/L_PE_rec_JpsiSigma,-1*L_PZ_rec_JpsiSigma/L_PE_rec_JpsiSigma);

		n = (beam.Cross(Lb_rec_JpsiSigma.Vect())).Unit();
		z1 = (L_rec_JpsiSigma.Vect()).Unit();
		z2 = (Jpsi_rec_JpsiSigma.Vect()).Unit();
		y1 = (n.Cross(z1));
		y2 = (n.Cross(z2));

		theta  = (L_rec_JpsiSigma.Vect()).Angle(n);
		theta1 = (p_rec_JpsiSigma.Vect()).Angle(z1);
		theta2 = (mu_rec_JpsiSigma.Vect()).Angle(z2);

		z1_prime = z1;
		z1_prime *= ((p_rec_JpsiSigma.Vect()).Dot(z1));

		z2_prime = z2;
		z2_prime *= ((mu_rec_JpsiSigma.Vect()).Dot(z2));

		phi1   = PiOver2() - (p_rec_JpsiSigma.Vect() - z1_prime).Angle(y1);
		phi2   = PiOver2() - (mu_rec_JpsiSigma.Vect() - z2_prime).Angle(y2);

		theta_rec_JpsiSigma->Fill(cos(theta));
		theta_rec_wt_JpsiSigma->Fill(cos(theta),tauWt_rec_JpsiSigma*gbWt_rec_JpsiSigma);

		theta1_rec_JpsiSigma->Fill(cos(theta1));
		theta1_rec_wt_JpsiSigma->Fill(cos(theta1),tauWt_rec_JpsiSigma*gbWt_rec_JpsiSigma);

		theta2_rec_JpsiSigma->Fill(cos(theta2));
		theta2_rec_wt_JpsiSigma->Fill(cos(theta2),tauWt_rec_JpsiSigma*gbWt_rec_JpsiSigma);

		phi1_rec_JpsiSigma->Fill(cos(phi1));
		phi1_rec_wt_JpsiSigma->Fill(cos(phi1),tauWt_rec_JpsiSigma*gbWt_rec_JpsiSigma);

		phi2_rec_JpsiSigma->Fill(cos(phi2));
		phi2_rec_wt_JpsiSigma->Fill(cos(phi2),tauWt_rec_JpsiSigma*gbWt_rec_JpsiSigma);
	}
	theta_gen_JpsiLambda->Sumw2();
	theta_gen_wt_JpsiLambda->Sumw2();
	theta_gen_JpsiSigma->Sumw2();
	theta_gen_wt_JpsiSigma->Sumw2();
	theta_rec_JpsiLambda->Sumw2();
	theta_rec_wt_JpsiLambda->Sumw2();
	theta_rec_JpsiSigma->Sumw2();
	theta_rec_wt_JpsiSigma->Sumw2();

	theta1_gen_JpsiLambda->Sumw2();
	theta1_gen_wt_JpsiLambda->Sumw2();
	theta1_gen_JpsiSigma->Sumw2();
	theta1_gen_wt_JpsiSigma->Sumw2();
	theta1_rec_JpsiLambda->Sumw2();
	theta1_rec_wt_JpsiLambda->Sumw2();
	theta1_rec_JpsiSigma->Sumw2();
	theta1_rec_wt_JpsiSigma->Sumw2();

	theta2_gen_JpsiLambda->Sumw2();
	theta2_gen_wt_JpsiLambda->Sumw2();
	theta2_gen_JpsiSigma->Sumw2();
	theta2_gen_wt_JpsiSigma->Sumw2();
	theta2_rec_JpsiLambda->Sumw2();
	theta2_rec_wt_JpsiLambda->Sumw2();
	theta2_rec_JpsiSigma->Sumw2();
	theta2_rec_wt_JpsiSigma->Sumw2();

	phi1_gen_JpsiLambda->Sumw2();
	phi1_gen_wt_JpsiLambda->Sumw2();
	phi1_gen_JpsiSigma->Sumw2();
	phi1_gen_wt_JpsiSigma->Sumw2();
	phi1_rec_JpsiLambda->Sumw2();
	phi1_rec_wt_JpsiLambda->Sumw2();
	phi1_rec_JpsiSigma->Sumw2();
	phi1_rec_wt_JpsiSigma->Sumw2();

	phi2_gen_JpsiLambda->Sumw2();
	phi2_gen_wt_JpsiLambda->Sumw2();
	phi2_gen_JpsiSigma->Sumw2();
	phi2_gen_wt_JpsiSigma->Sumw2();
	phi2_rec_JpsiLambda->Sumw2();
	phi2_rec_wt_JpsiLambda->Sumw2();
	phi2_rec_JpsiSigma->Sumw2();
	phi2_rec_wt_JpsiSigma->Sumw2();

	TH1D *eff_theta_JpsiLambda = (TH1D*)theta_rec_JpsiLambda->Clone();
	eff_theta_JpsiLambda->SetNameTitle("eff_theta_JpsiLambda","eff_theta_JpsiLambda");
	eff_theta_JpsiLambda->Divide(theta_gen_JpsiLambda);

	TH1D *eff_theta_JpsiLambda_wt = (TH1D*)theta_rec_wt_JpsiLambda->Clone();
	eff_theta_JpsiLambda_wt->SetNameTitle("eff_theta_JpsiLambda_wt","eff_theta_JpsiLambda_wt");
	eff_theta_JpsiLambda_wt->Divide(theta_gen_wt_JpsiLambda);

	TH1D *eff_theta_JpsiSigma = (TH1D*)theta_rec_JpsiSigma->Clone();
	eff_theta_JpsiSigma->SetNameTitle("eff_theta_JpsiSigma","eff_theta_JpsiSigma");
	eff_theta_JpsiSigma->Divide(theta_gen_JpsiSigma);

	TH1D *eff_theta_JpsiSigma_wt = (TH1D*)theta_rec_wt_JpsiSigma->Clone();
	eff_theta_JpsiSigma_wt->SetNameTitle("eff_theta_JpsiSigma_wt","eff_theta_JpsiSigma_wt");
	eff_theta_JpsiSigma_wt->Divide(theta_gen_wt_JpsiSigma);

	TH1D *ratio_theta = (TH1D*)eff_theta_JpsiLambda->Clone();
	ratio_theta->SetNameTitle("ratio_theta","ratio_theta");
	ratio_theta->GetYaxis()->SetTitle(title_y);
	ratio_theta->GetXaxis()->SetTitle("cos(#theta)");
	ratio_theta->Divide(eff_theta_JpsiLambda,eff_theta_JpsiSigma,(1.0/eff_theta_JpsiLambda->Integral()),(1.0/eff_theta_JpsiSigma->Integral()));
	chi2 = ratio_theta->Chisquare(unif);
	cout<<"ratio_theta chi2/ndf = "<<chi2/50<<endl;

	TH1D *ratio_theta_wt = (TH1D*)eff_theta_JpsiLambda_wt->Clone();
	ratio_theta_wt->SetNameTitle("ratio_theta_wt","ratio_theta_wt");
	ratio_theta_wt->GetYaxis()->SetTitle(title_y_wt);
	ratio_theta_wt->GetXaxis()->SetTitle("cos(#theta)");
	ratio_theta_wt->Divide(eff_theta_JpsiLambda_wt,eff_theta_JpsiSigma_wt,(1.0/eff_theta_JpsiLambda_wt->Integral()),(1.0/eff_theta_JpsiSigma_wt->Integral()));
	chi2 = ratio_theta_wt->Chisquare(unif);
	cout<<"ratio_theta_wt chi2/ndf = "<<chi2/50<<endl;

	TH1D *eff_theta1_JpsiLambda = (TH1D*)theta1_rec_JpsiLambda->Clone();
	eff_theta1_JpsiLambda->SetNameTitle("eff_theta1_JpsiLambda","eff_theta1_JpsiLambda");
	eff_theta1_JpsiLambda->Divide(theta1_gen_JpsiLambda);

	TH1D *eff_theta1_JpsiLambda_wt = (TH1D*)theta1_rec_wt_JpsiLambda->Clone();
	eff_theta1_JpsiLambda_wt->SetNameTitle("eff_theta1_JpsiLambda_wt","eff_theta1_JpsiLambda_wt");
	eff_theta1_JpsiLambda_wt->Divide(theta1_gen_wt_JpsiLambda);

	TH1D *eff_theta1_JpsiSigma = (TH1D*)theta1_rec_JpsiSigma->Clone();
	eff_theta1_JpsiSigma->SetNameTitle("eff_theta1_JpsiSigma","eff_theta1_JpsiSigma");
	eff_theta1_JpsiSigma->Divide(theta1_gen_JpsiSigma);

	TH1D *eff_theta1_JpsiSigma_wt = (TH1D*)theta1_rec_wt_JpsiSigma->Clone();
	eff_theta1_JpsiSigma_wt->SetNameTitle("eff_theta1_JpsiSigma_wt","eff_theta1_JpsiSigma_wt");
	eff_theta1_JpsiSigma_wt->Divide(theta1_gen_wt_JpsiSigma);

	TH1D *ratio_theta1 = (TH1D*)eff_theta1_JpsiLambda->Clone();
	ratio_theta1->SetNameTitle("ratio_theta1","ratio_theta1");
	ratio_theta1->GetYaxis()->SetTitle(title_y);
	ratio_theta1->GetXaxis()->SetTitle("cos(#theta_{1})");
	ratio_theta1->Divide(eff_theta1_JpsiLambda,eff_theta1_JpsiSigma,(1.0/eff_theta1_JpsiLambda->Integral()),(1.0/eff_theta1_JpsiSigma->Integral()));
	chi2 = ratio_theta1->Chisquare(unif);
	cout<<"ratio_theta1 chi2/ndf = "<<chi2/50<<endl;

	TH1D *ratio_theta1_wt = (TH1D*)eff_theta1_JpsiLambda_wt->Clone();
	ratio_theta1_wt->SetNameTitle("ratio_theta1_wt","ratio_theta1_wt");
	ratio_theta1_wt->GetYaxis()->SetTitle(title_y_wt);
	ratio_theta1_wt->GetXaxis()->SetTitle("cos(#theta_{1})");
	ratio_theta1_wt->Divide(eff_theta1_JpsiLambda_wt,eff_theta1_JpsiSigma_wt,(1.0/eff_theta1_JpsiLambda_wt->Integral()),(1.0/eff_theta1_JpsiSigma_wt->Integral()));
	chi2 = ratio_theta1_wt->Chisquare(unif);
	cout<<"ratio_theta1_wt chi2/ndf = "<<chi2/50<<endl;

	TH1D *eff_theta2_JpsiLambda = (TH1D*)theta2_rec_JpsiLambda->Clone();
	eff_theta2_JpsiLambda->SetNameTitle("eff_theta2_JpsiLambda","eff_theta2_JpsiLambda");
	eff_theta2_JpsiLambda->Divide(theta2_gen_JpsiLambda);

	TH1D *eff_theta2_JpsiLambda_wt = (TH1D*)theta2_rec_wt_JpsiLambda->Clone();
	eff_theta2_JpsiLambda_wt->SetNameTitle("eff_theta2_JpsiLambda_wt","eff_theta2_JpsiLambda_wt");
	eff_theta2_JpsiLambda_wt->Divide(theta2_gen_wt_JpsiLambda);

	TH1D *eff_theta2_JpsiSigma = (TH1D*)theta2_rec_JpsiSigma->Clone();
	eff_theta2_JpsiSigma->SetNameTitle("eff_theta2_JpsiSigma","eff_theta2_JpsiSigma");
	eff_theta2_JpsiSigma->Divide(theta2_gen_JpsiSigma);

	TH1D *eff_theta2_JpsiSigma_wt = (TH1D*)theta2_rec_wt_JpsiSigma->Clone();
	eff_theta2_JpsiSigma_wt->SetNameTitle("eff_theta2_JpsiSigma_wt","eff_theta2_JpsiSigma_wt");
	eff_theta2_JpsiSigma_wt->Divide(theta2_gen_wt_JpsiSigma);

	TH1D *ratio_theta2 = (TH1D*)eff_theta2_JpsiLambda->Clone();
	ratio_theta2->SetNameTitle("ratio_theta2","ratio_theta2");
	ratio_theta2->GetYaxis()->SetTitle(title_y);
	ratio_theta2->GetXaxis()->SetTitle("cos(#theta_{2})");
	ratio_theta2->Divide(eff_theta2_JpsiLambda,eff_theta2_JpsiSigma,(1.0/eff_theta2_JpsiLambda->Integral()),(1.0/eff_theta2_JpsiSigma->Integral()));
	chi2 = ratio_theta2->Chisquare(unif);
	cout<<"ratio_theta2 chi2/ndf = "<<chi2/50<<endl;

	TH1D *ratio_theta2_wt = (TH1D*)eff_theta2_JpsiLambda_wt->Clone();
	ratio_theta2_wt->SetNameTitle("ratio_theta2_wt","ratio_theta2_wt");
	ratio_theta2_wt->GetYaxis()->SetTitle(title_y_wt);
	ratio_theta2_wt->GetXaxis()->SetTitle("cos(#theta_{2})");
	ratio_theta2_wt->Divide(eff_theta2_JpsiLambda_wt,eff_theta2_JpsiSigma_wt,(1.0/eff_theta2_JpsiLambda_wt->Integral()),(1.0/eff_theta2_JpsiSigma_wt->Integral()));
	chi2 = ratio_theta2_wt->Chisquare(unif);
	cout<<"ratio_theta2_wt chi2/ndf = "<<chi2/50<<endl;

	TH1D *eff_phi1_JpsiLambda = (TH1D*)phi1_rec_JpsiLambda->Clone();
	eff_phi1_JpsiLambda->SetNameTitle("eff_phi1_JpsiLambda","eff_phi1_JpsiLambda");
	eff_phi1_JpsiLambda->Divide(phi1_gen_JpsiLambda);

	TH1D *eff_phi1_JpsiLambda_wt = (TH1D*)phi1_rec_wt_JpsiLambda->Clone();
	eff_phi1_JpsiLambda_wt->SetNameTitle("eff_phi1_JpsiLambda_wt","eff_phi1_JpsiLambda_wt");
	eff_phi1_JpsiLambda_wt->Divide(phi1_gen_wt_JpsiLambda);

	TH1D *eff_phi1_JpsiSigma = (TH1D*)phi1_rec_JpsiSigma->Clone();
	eff_phi1_JpsiSigma->SetNameTitle("eff_phi1_JpsiSigma","eff_phi1_JpsiSigma");
	eff_phi1_JpsiSigma->Divide(phi1_gen_JpsiSigma);

	TH1D *eff_phi1_JpsiSigma_wt = (TH1D*)phi1_rec_wt_JpsiSigma->Clone();
	eff_phi1_JpsiSigma_wt->SetNameTitle("eff_phi1_JpsiSigma_wt","eff_phi1_JpsiSigma_wt");
	eff_phi1_JpsiSigma_wt->Divide(phi1_gen_wt_JpsiSigma);

	TH1D *ratio_phi1 = (TH1D*)eff_phi1_JpsiLambda->Clone();
	ratio_phi1->SetNameTitle("ratio_phi1","ratio_phi1");
	ratio_phi1->GetYaxis()->SetTitle(title_y);
	ratio_phi1->GetXaxis()->SetTitle("cos(#phi_{1})");
	ratio_phi1->Divide(eff_phi1_JpsiLambda,eff_phi1_JpsiSigma,(1.0/eff_phi1_JpsiLambda->Integral()),(1.0/eff_phi1_JpsiSigma->Integral()));
	chi2 = ratio_phi1->Chisquare(unif);
	cout<<"ratio_phi1 chi2/ndf = "<<chi2/50<<endl;

	TH1D *ratio_phi1_wt = (TH1D*)eff_phi1_JpsiLambda_wt->Clone();
	ratio_phi1_wt->SetNameTitle("ratio_phi1_wt","ratio_phi1_wt");
	ratio_phi1_wt->GetYaxis()->SetTitle(title_y_wt);
	ratio_phi1_wt->GetXaxis()->SetTitle("cos(#phi_{1})");
	ratio_phi1_wt->Divide(eff_phi1_JpsiLambda_wt,eff_phi1_JpsiSigma_wt,(1.0/eff_phi1_JpsiLambda_wt->Integral()),(1.0/eff_phi1_JpsiSigma_wt->Integral()));
	chi2 = ratio_phi1_wt->Chisquare(unif);
	cout<<"ratio_phi1_wt chi2/ndf = "<<chi2/50<<endl;

	TH1D *eff_phi2_JpsiLambda = (TH1D*)phi2_rec_JpsiLambda->Clone();
	eff_phi2_JpsiLambda->SetNameTitle("eff_phi2_JpsiLambda","eff_phi2_JpsiLambda");
	eff_phi2_JpsiLambda->Divide(phi2_gen_JpsiLambda);

	TH1D *eff_phi2_JpsiLambda_wt = (TH1D*)phi2_rec_wt_JpsiLambda->Clone();
	eff_phi2_JpsiLambda_wt->SetNameTitle("eff_phi2_JpsiLambda_wt","eff_phi2_JpsiLambda_wt");
	eff_phi2_JpsiLambda_wt->Divide(phi2_gen_wt_JpsiLambda);

	TH1D *eff_phi2_JpsiSigma = (TH1D*)phi2_rec_JpsiSigma->Clone();
	eff_phi2_JpsiSigma->SetNameTitle("eff_phi2_JpsiSigma","eff_phi2_JpsiSigma");
	eff_phi2_JpsiSigma->Divide(phi2_gen_JpsiSigma);

	TH1D *eff_phi2_JpsiSigma_wt = (TH1D*)phi2_rec_wt_JpsiSigma->Clone();
	eff_phi2_JpsiSigma_wt->SetNameTitle("eff_phi2_JpsiSigma_wt","eff_phi2_JpsiSigma_wt");
	eff_phi2_JpsiSigma_wt->Divide(phi2_gen_wt_JpsiSigma);

	TH1D *ratio_phi2 = (TH1D*)eff_phi2_JpsiLambda->Clone();
	ratio_phi2->SetNameTitle("ratio_phi2","ratio_phi2");
	ratio_phi2->GetYaxis()->SetTitle(title_y);
	ratio_phi2->GetXaxis()->SetTitle("cos(#phi_{2})");
	ratio_phi2->Divide(eff_phi2_JpsiLambda,eff_phi2_JpsiSigma,(1.0/eff_phi2_JpsiLambda->Integral()),(1.0/eff_phi2_JpsiSigma->Integral()));
	chi2 = ratio_phi2->Chisquare(unif);
	cout<<"ratio_phi2 chi2/ndf = "<<chi2/50<<endl;

	TH1D *ratio_phi2_wt = (TH1D*)eff_phi2_JpsiLambda_wt->Clone();
	ratio_phi2_wt->SetNameTitle("ratio_phi2_wt","ratio_phi2_wt");
	ratio_phi2_wt->GetYaxis()->SetTitle(title_y_wt);
	ratio_phi2_wt->GetXaxis()->SetTitle("cos(#phi_{2})");
	ratio_phi2_wt->Divide(eff_phi2_JpsiLambda_wt,eff_phi2_JpsiSigma_wt,(1.0/eff_phi2_JpsiLambda_wt->Integral()),(1.0/eff_phi2_JpsiSigma_wt->Integral()));
	chi2 = ratio_phi2_wt->Chisquare(unif);
	cout<<"ratio_phi2_wt chi2/ndf = "<<chi2/50<<endl;

	// theta_gen_JpsiLambda->Draw();
	// new TCanvas();
	// theta_gen_wt_JpsiLambda->Draw();
	// new TCanvas();
	// theta_gen_JpsiSigma->Draw();
	// new TCanvas();
	// theta_gen_wt_JpsiSigma->Draw();
	// new TCanvas();
	// theta_rec_JpsiLambda->Draw();
	// new TCanvas();
	// theta_rec_wt_JpsiLambda->Draw();
	// new TCanvas();
	// theta_rec_JpsiSigma->Draw();
	// new TCanvas();
	// theta_rec_wt_JpsiSigma->Draw();
	// new TCanvas();
	//
	// eff_JpsiLambda->Draw();
	// new TCanvas();
	// eff_JpsiLambda_wt->Draw();
	// new TCanvas();
	//
	// eff_JpsiSigma->Draw();
	// new TCanvas();
	// eff_JpsiSigma_wt->Draw();
	// new TCanvas();

	TCanvas *c1 = new TCanvas();
	ratio_theta->Draw();
	c1->Update();
	line1->Draw();

	TCanvas *c2 = new TCanvas();
	ratio_theta_wt->Draw();
	c2->Update();
	line1->Draw();

	TCanvas *c3 = new TCanvas();
	ratio_theta1->Draw();
	c3->Update();
	line1->Draw();

	TCanvas *c4 = new TCanvas();
	ratio_theta1_wt->Draw();
	c4->Update();
	line1->Draw();

	TCanvas *c5 = new TCanvas();
	ratio_theta2->Draw();
	c5->Update();
	line1->Draw();

	TCanvas *c6 = new TCanvas();
	ratio_theta2_wt->Draw();
	c6->Update();
	line1->Draw();

	TCanvas *c7 = new TCanvas();
	ratio_phi1->Draw();
	c7->Update();
	line2->Draw();

	TCanvas *c8 = new TCanvas();
	ratio_phi1_wt->Draw();
	c8->Update();
	line2->Draw();

	TCanvas *c9 = new TCanvas();
	ratio_phi2->Draw();
	c9->Update();
	line2->Draw();

	TCanvas *c10 = new TCanvas();
	ratio_phi2_wt->Draw();
	c10->Update();
	line2->Draw();
}
