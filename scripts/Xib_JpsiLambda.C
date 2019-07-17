#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <iostream>
#include <fstream>

using namespace std;

#define Open TFile::Open
void Xib_JpsiLambda()
{
	//Combined BF result
	Float_t BF_Ratio_Comb     = 0.0;
	Float_t StatErr_BF_Ratio_Comb = 0.0;
	Float_t SystErr_BF_Ratio_Comb = 0.0;
	Float_t Err_BF_Ratio_Comb = 0.0;

	//Combined BF result - with MC RW
	Float_t BF_Ratio_Comb_wt     = 0.0;
	Float_t StatErr_BF_Ratio_Comb_wt = 0.0;
	Float_t SystErr_BF_Ratio_Comb_wt = 0.0;
	Float_t Err_BF_Ratio_Comb_wt = 0.0;

	//Run 1 BF Result
	Float_t BF_Ratio_Run1     = 0.0;
	Float_t StatErr_BF_Ratio_Run1 = 0.0;
	Float_t SystErr_BF_Ratio_Run1 = 0.0;
	Float_t Err_BF_Ratio_Run1 = 0.0;

	//Run 2 BF Result
	Float_t BF_Ratio_Run2     = 0.0;
	Float_t StatErr_BF_Ratio_Run2 = 0.0;
	Float_t SystErr_BF_Ratio_Run2 = 0.0;
	Float_t Err_BF_Ratio_Run2 = 0.0;

	//Weighted Run 1 BF Result
	Float_t BF_Ratio_Run1_wt     = 0.0;
	Float_t StatErr_BF_Ratio_Run1_wt = 0.0;
	Float_t SystErr_BF_Ratio_Run1_wt = 0.0;
	Float_t Err_BF_Ratio_Run1_wt = 0.0;

	//Weighted Run 2 BF Result
	Float_t BF_Ratio_Run2_wt     = 0.0;
	Float_t StatErr_BF_Ratio_Run2_wt = 0.0;
	Float_t SystErr_BF_Ratio_Run2_wt = 0.0;
	Float_t Err_BF_Ratio_Run2_wt = 0.0;

	//Fit yield for Xib0 -> J/psi Lambda
	Float_t N_Xib0_JpsiLambda_Run1 = 6.75;
	Float_t N_Xib0_JpsiLambda_Run2 = 14.36;

	Float_t StatErr_N_Xib0_JpsiLambda_Run1 = 2.96;//44% error
	Float_t StatErr_N_Xib0_JpsiLambda_Run2 = 4.58;//32% error

	Float_t SystErr_N_Xib0_JpsiLambda_Run1 = 0.03;
	Float_t SystErr_N_Xib0_JpsiLambda_Run2 = 0.26;

	//Total MC Eff for Xib0 -> J/psi Lambda
	Float_t Eff_Xib0_JpsiLambda_Run1 = 0.0;
	Float_t Eff_Xib0_JpsiLambda_Run2 = 0.0;

	Float_t Err_Eff_Xib0_JpsiLambda_Run1 = 0.0;
	Float_t Err_Eff_Xib0_JpsiLambda_Run2 = 0.0;

	//Generator MC Eff for Xib0 -> J/psi Lambda
	Float_t GenEff_Xib0_JpsiLambda_Run1 = 0.0;
	Float_t GenEff_Xib0_JpsiLambda_Run2 = 0.0;

	Float_t Err_GenEff_Xib0_JpsiLambda_Run1 = 0.0;
	Float_t Err_GenEff_Xib0_JpsiLambda_Run2 = 0.0;

	//Reconstructed MC Eff for Xib0 -> J/psi Lambda
	Float_t RecEff_Xib0_JpsiLambda_Run1 = 0.0;
	Float_t RecEff_Xib0_JpsiLambda_Run2 = 0.0;

	Float_t Err_RecEff_Xib0_JpsiLambda_Run1 = 0.0;
	Float_t Err_RecEff_Xib0_JpsiLambda_Run2 = 0.0;

	//Weighted Total MC Eff for Xib0 -> J/psi Lambda
	Float_t Eff_Xib0_JpsiLambda_Run1_wt = 0.0;
	Float_t Eff_Xib0_JpsiLambda_Run2_wt = 0.0;

	Float_t Err_Eff_Xib0_JpsiLambda_Run1_wt = 0.0;
	Float_t Err_Eff_Xib0_JpsiLambda_Run2_wt = 0.0;

	//Weighted Reconstruction MC Eff for Xib0 -> J/psi Lambda
	Float_t RecEff_Xib0_JpsiLambda_Run1_wt = 0.0;
	Float_t RecEff_Xib0_JpsiLambda_Run2_wt = 0.0;

	Float_t Err_RecEff_Xib0_JpsiLambda_Run1_wt = 0.0;
	Float_t Err_RecEff_Xib0_JpsiLambda_Run2_wt = 0.0;

	//Reconstruction MC Eff for Xib- -> J/psi Xi-
	Float_t RecEff_Xibm_JpsiXim_Run1 = 0.0631/100;
	Float_t RecEff_Xibm_JpsiXim_Run2 = 0.0650/100;

	Float_t Err_RecEff_Xibm_JpsiXim_Run1 = 0.0014/100;
	Float_t Err_RecEff_Xibm_JpsiXim_Run2 = 0.0014/100;

	//Weighted Reconstruction MC Eff for Xib- -> J/psi Xi-
	Float_t RecEff_Xibm_JpsiXim_Run1_wt = 0.0683/100;
	Float_t RecEff_Xibm_JpsiXim_Run2_wt = 0.0681/100;

	Float_t Err_RecEff_Xibm_JpsiXim_Run1_wt = 0.0017/100;
	Float_t Err_RecEff_Xibm_JpsiXim_Run2_wt = 0.0015/100;

	//Generator MC Eff for Xib- -> J/psi Xi-
	Float_t GenEff_Xibm_JpsiXim_Run1 = 0.17267;
	Float_t GenEff_Xibm_JpsiXim_Run2 = 0.18019;

	Float_t Err_GenEff_Xibm_JpsiXim_Run1 = 0.00056;
	Float_t Err_GenEff_Xibm_JpsiXim_Run2 = 0.00043;

	Float_t Eff_Xibm_JpsiXim_Run1 = GenEff_Xibm_JpsiXim_Run1 * RecEff_Xibm_JpsiXim_Run1;
	Float_t Eff_Xibm_JpsiXim_Run2 = GenEff_Xibm_JpsiXim_Run2 * RecEff_Xibm_JpsiXim_Run2;

	Float_t Err_Eff_Xibm_JpsiXim_Run1 = Eff_Xibm_JpsiXim_Run1*sqrt(pow(Err_GenEff_Xibm_JpsiXim_Run1/GenEff_Xibm_JpsiXim_Run1,2)
	                                                               + pow(Err_RecEff_Xibm_JpsiXim_Run1/RecEff_Xibm_JpsiXim_Run1,2));
	Float_t Err_Eff_Xibm_JpsiXim_Run2 = Eff_Xibm_JpsiXim_Run2*sqrt(pow(Err_GenEff_Xibm_JpsiXim_Run2/GenEff_Xibm_JpsiXim_Run2,2)
	                                                               + pow(Err_RecEff_Xibm_JpsiXim_Run2/RecEff_Xibm_JpsiXim_Run2,2));

	cout<<"Run 1 Unweighted Xib- -> Jpsi Xi- Eff = "<<Eff_Xibm_JpsiXim_Run1*100<<"% +/- "<<Err_Eff_Xibm_JpsiXim_Run1*100<<endl;
	cout<<"Run 2 Unweighted Xib- -> Jpsi Xi- Eff = "<<Eff_Xibm_JpsiXim_Run2*100<<"% +/- "<<Err_Eff_Xibm_JpsiXim_Run2*100<<endl;

	Float_t Eff_Xibm_JpsiXim_Run1_wt = GenEff_Xibm_JpsiXim_Run1 * RecEff_Xibm_JpsiXim_Run1_wt;
	Float_t Eff_Xibm_JpsiXim_Run2_wt = GenEff_Xibm_JpsiXim_Run2 * RecEff_Xibm_JpsiXim_Run2_wt;

	Float_t Err_Eff_Xibm_JpsiXim_Run1_wt = Eff_Xibm_JpsiXim_Run1_wt*sqrt(pow(Err_GenEff_Xibm_JpsiXim_Run1/GenEff_Xibm_JpsiXim_Run1,2)
	                                                                     + pow(Err_RecEff_Xibm_JpsiXim_Run1_wt/RecEff_Xibm_JpsiXim_Run1_wt,2));
	Float_t Err_Eff_Xibm_JpsiXim_Run2_wt = Eff_Xibm_JpsiXim_Run2_wt*sqrt(pow(Err_GenEff_Xibm_JpsiXim_Run2/GenEff_Xibm_JpsiXim_Run2,2)
	                                                                     + pow(Err_RecEff_Xibm_JpsiXim_Run2_wt/RecEff_Xibm_JpsiXim_Run2_wt,2));

	cout<<"Run 1 Weighted Xib- -> Jpsi Xi- Eff = "<<Eff_Xibm_JpsiXim_Run1_wt*100<<"% +/- "<<Err_Eff_Xibm_JpsiXim_Run1_wt*100<<endl;
	cout<<"Run 2 Weighted Xib- -> Jpsi Xi- Eff = "<<Eff_Xibm_JpsiXim_Run2_wt*100<<"% +/- "<<Err_Eff_Xibm_JpsiXim_Run2_wt*100<<endl;

	Float_t B_Xi0_LambdaPi0 = 99.524/100;
	Float_t Err_B_Xi0_LambdaPi0 = 0.012/100;

	Float_t N_Xibm_JpsiXim_Run1 = 36.361;//Yield of fully reco'd Xib- -> J/psi Xi-, comes from fit to data
	Float_t StatErr_N_Xibm_JpsiXim_Run1 = 7.757;//Fit error
	Float_t SystErr_N_Xibm_JpsiXim_Run1 = N_Xibm_JpsiXim_Run1*0.09;//Syst error coming from model choice
	Float_t Err_N_Xibm_JpsiXim_Run1 = sqrt(pow(StatErr_N_Xibm_JpsiXim_Run1,2) + pow(SystErr_N_Xibm_JpsiXim_Run1,2));

	Float_t N_Xibm_JpsiXim_Run2 = 173.281;
	Float_t StatErr_N_Xibm_JpsiXim_Run2 = 16.421;
	Float_t SystErr_N_Xibm_JpsiXim_Run2 = N_Xibm_JpsiXim_Run2*0.06;
	Float_t Err_N_Xibm_JpsiXim_Run2 = sqrt(pow(StatErr_N_Xibm_JpsiXim_Run2,2) + pow(SystErr_N_Xibm_JpsiXim_Run2,2));

	Float_t Scale_Factor = 1.087;
	Float_t Err_Scale_Factor = 0.036;

	Float_t N_Xib0_JpsiXi0_Run1 = N_Xibm_JpsiXim_Run1/Scale_Factor;
	Float_t StatErr_N_Xib0_JpsiXi0_Run1 = N_Xib0_JpsiXi0_Run1*sqrt(pow(StatErr_N_Xibm_JpsiXim_Run1/N_Xibm_JpsiXim_Run1,2));
	Float_t SystErr_N_Xib0_JpsiXi0_Run1 = N_Xib0_JpsiXi0_Run1*sqrt(pow(SystErr_N_Xibm_JpsiXim_Run1/N_Xibm_JpsiXim_Run1,2)
	                                                               + pow(Err_Scale_Factor/Scale_Factor,2));

	Float_t Err_N_Xib0_JpsiXi0_Run1 = sqrt(pow(StatErr_N_Xib0_JpsiXi0_Run1,2) + pow(SystErr_N_Xib0_JpsiXi0_Run1,2));

	Float_t N_Xib0_JpsiXi0_Run2 = N_Xibm_JpsiXim_Run2/Scale_Factor;
	Float_t StatErr_N_Xib0_JpsiXi0_Run2 = N_Xib0_JpsiXi0_Run2*sqrt(pow(StatErr_N_Xibm_JpsiXim_Run2/N_Xibm_JpsiXim_Run2,2));
	Float_t SystErr_N_Xib0_JpsiXi0_Run2 = N_Xib0_JpsiXi0_Run2*sqrt(pow(SystErr_N_Xibm_JpsiXim_Run2/N_Xibm_JpsiXim_Run2,2)
	                                                               + pow(Err_Scale_Factor/Scale_Factor,2));

	Float_t Err_N_Xib0_JpsiXi0_Run2 = sqrt(pow(StatErr_N_Xib0_JpsiXi0_Run2,2) + pow(SystErr_N_Xib0_JpsiXi0_Run2,2));

	Int_t nGen_Run1 = 0, nGen_Run2 = 0;
	Int_t nGen_Run1_wt = 0, nGen_Run2_wt = 0;

	Float_t trackingUnc = 0.05; //5%
	Float_t xiVtxUnc                = 0.014;//1.4%

	//Get efficiencies
	TFile *fileIn_nonZero_Run1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/Xib0/run1/xib0_cutoutks_LL_nonZeroTracks_noPID.root","READ");
	TTree *treeIn_nonZero_Run1 = (TTree*)fileIn_nonZero_Run1->Get("MyTuple");
	treeIn_nonZero_Run1->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/Xib0/run1/xib0_LL_FinalBDT2_iso2_v0_noPID.root");

	TFile *fileIn_Zero_Run1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/Xib0/run1/xib0_cutoutks_LL_ZeroTracks_noPID.root","READ");
	TTree *treeIn_Zero_Run1 = (TTree*)fileIn_Zero_Run1->Get("MyTuple");
	treeIn_Zero_Run1->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/Xib0/run1/xib0_zeroTracksLL_FinalBDT2_noPID.root");

	treeIn_nonZero_Run1->Draw("gb_wts_new>>wt_Run1_nonZero","BDT2 > 0.475","goff");
	treeIn_Zero_Run1->Draw("gb_wts_new>>wt_Run1_Zero","BDT2 > 0.365","goff");

	TH1F *wt_Run1_nonZero = (TH1F*)gDirectory->Get("wt_Run1_nonZero");
	TH1F *wt_Run1_Zero = (TH1F*)gDirectory->Get("wt_Run1_Zero");

	fstream genFile_Xib0_Run1;
	genFile_Xib0_Run1.open("../logs/mc/JpsiLambda/Xib0/run1/gen_log.txt");

	genFile_Xib0_Run1>>nGen_Run1; //Get number of generated events

	TFile *genWtsFile_Run1 = nullptr;
	genWtsFile_Run1 = Open("../rootFiles/mcFiles/JpsiLambda/Xib0/run1/RW/gbWeights_gen_new.root");

	TTree *genWtsTree_Run1 = (TTree*)genWtsFile_Run1->Get("MyTuple");

	genWtsTree_Run1->Draw("gb_wts_new>>genWt_Run1","","goff");

	TH1F *genWt_Run1 = (TH1F*)gDirectory->Get("genWt_Run1");
	nGen_Run1_wt = genWt_Run1->GetEntries()*genWt_Run1->GetMean();

	fstream genEffFile_Run1;
	genEffFile_Run1.open("../logs/mc/JpsiLambda/Xib0/run1/Generator_Effs_Combined.txt");

	genEffFile_Run1>>GenEff_Xib0_JpsiLambda_Run1;    //Get generator efficiency
	genEffFile_Run1>>Err_GenEff_Xib0_JpsiLambda_Run1;//and error on above

	cout<<"Run 1 Generator Effs = "<<GenEff_Xib0_JpsiLambda_Run1*100
	    <<" % +/- "<<Err_GenEff_Xib0_JpsiLambda_Run1*100<<" %"<<endl;

	Float_t num_Run1_wt = (wt_Run1_nonZero->GetMean()*wt_Run1_nonZero->GetEntries()) +
	                      (wt_Run1_Zero->GetMean()*wt_Run1_Zero->GetEntries());

	Float_t num_Run1 = treeIn_nonZero_Run1->GetEntries("BDT2 > 0.475")
	                   + treeIn_Zero_Run1->GetEntries("BDT2 > 0.365");

	RecEff_Xib0_JpsiLambda_Run1     = num_Run1*1.0/nGen_Run1;//Calc. reco eff.
	Err_RecEff_Xib0_JpsiLambda_Run1 = sqrt(RecEff_Xib0_JpsiLambda_Run1*(1-RecEff_Xib0_JpsiLambda_Run1)/nGen_Run1); //statistical error on recon. eff.

	RecEff_Xib0_JpsiLambda_Run1_wt     = num_Run1_wt*1.0/nGen_Run1_wt; //Calc. weighted reco eff.
	Err_RecEff_Xib0_JpsiLambda_Run1_wt = sqrt(RecEff_Xib0_JpsiLambda_Run1_wt*(1-RecEff_Xib0_JpsiLambda_Run1_wt)/nGen_Run1_wt); //statistical error on weighted recon. eff.

	cout<<"Run 1 UNWEIGHTED Lambda Recons. Effs = "<<RecEff_Xib0_JpsiLambda_Run1*100
	    <<" % +/- "<<Err_RecEff_Xib0_JpsiLambda_Run1*100<<" %"<<endl;

	cout<<"Run 1 WEIGHTED Lambda Recons. Effs = "<<RecEff_Xib0_JpsiLambda_Run1_wt*100
	    <<" % +/- "<<Err_RecEff_Xib0_JpsiLambda_Run1_wt*100<<" %"<<endl;

	Eff_Xib0_JpsiLambda_Run1     = RecEff_Xib0_JpsiLambda_Run1*GenEff_Xib0_JpsiLambda_Run1; // Calc. total eff.
	Err_Eff_Xib0_JpsiLambda_Run1 = Eff_Xib0_JpsiLambda_Run1*sqrt(pow((Err_GenEff_Xib0_JpsiLambda_Run1/GenEff_Xib0_JpsiLambda_Run1),2) +
	                                                             pow((Err_RecEff_Xib0_JpsiLambda_Run1/RecEff_Xib0_JpsiLambda_Run1),2));                                                                                                                                                                                                                                                                   // and stat error on tot. eff.

	Eff_Xib0_JpsiLambda_Run1_wt     = RecEff_Xib0_JpsiLambda_Run1_wt*GenEff_Xib0_JpsiLambda_Run1; // Calc. total eff.
	Err_Eff_Xib0_JpsiLambda_Run1_wt = Eff_Xib0_JpsiLambda_Run1_wt*sqrt(pow((Err_GenEff_Xib0_JpsiLambda_Run1/GenEff_Xib0_JpsiLambda_Run1),2) +
	                                                                   pow((Err_RecEff_Xib0_JpsiLambda_Run1_wt/RecEff_Xib0_JpsiLambda_Run1_wt),2)); // and stat error on tot. eff.
	cout<<"************************************************"<<endl;
	cout<<"Run 1 UNWEIGHTED Eff = "<<Eff_Xib0_JpsiLambda_Run1*100
	    <<" % +/- "<<Err_Eff_Xib0_JpsiLambda_Run1*100<<" %"<<endl;
	cout<<"Run 1 WEIGHTED Eff = "<<Eff_Xib0_JpsiLambda_Run1_wt*100
	    <<" % +/- "<<Err_Eff_Xib0_JpsiLambda_Run1_wt*100<<" %"<<endl;
	cout<<"************************************************"<<endl;

	BF_Ratio_Run1    = (N_Xib0_JpsiLambda_Run1/Eff_Xib0_JpsiLambda_Run1) * B_Xi0_LambdaPi0 / (N_Xib0_JpsiXi0_Run1/Eff_Xibm_JpsiXim_Run1);
	BF_Ratio_Run1_wt = (N_Xib0_JpsiLambda_Run1/Eff_Xib0_JpsiLambda_Run1_wt) * B_Xi0_LambdaPi0 / (N_Xib0_JpsiXi0_Run1/Eff_Xibm_JpsiXim_Run1_wt);

	SystErr_BF_Ratio_Run1 = BF_Ratio_Run1 * sqrt(pow(SystErr_N_Xib0_JpsiLambda_Run1/N_Xib0_JpsiLambda_Run1,2)
	                                             + pow(SystErr_N_Xib0_JpsiXi0_Run1/N_Xib0_JpsiXi0_Run1,2)
	                                             + pow(Err_Eff_Xib0_JpsiLambda_Run1/Eff_Xib0_JpsiLambda_Run1,2)
	                                             + pow(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0,2)
	                                             + pow(Err_Eff_Xibm_JpsiXim_Run1/Eff_Xibm_JpsiXim_Run1,2)
	                                             + pow(trackingUnc,2)
	                                             + pow(xiVtxUnc,2));

	StatErr_BF_Ratio_Run1 = BF_Ratio_Run1 * sqrt(pow(StatErr_N_Xib0_JpsiLambda_Run1/N_Xib0_JpsiLambda_Run1,2)
	                                             + pow(StatErr_N_Xib0_JpsiXi0_Run1/N_Xib0_JpsiXi0_Run1,2));

	Err_BF_Ratio_Run1 = sqrt(pow(StatErr_BF_Ratio_Run1,2) + pow(SystErr_BF_Ratio_Run1,2));

	SystErr_BF_Ratio_Run1_wt = BF_Ratio_Run1_wt * sqrt(pow(SystErr_N_Xib0_JpsiLambda_Run1/N_Xib0_JpsiLambda_Run1,2)
	                                                   + pow(SystErr_N_Xib0_JpsiXi0_Run1/N_Xib0_JpsiXi0_Run1,2)
	                                                   + pow(Err_Eff_Xib0_JpsiLambda_Run1_wt/Eff_Xib0_JpsiLambda_Run1_wt,2)
	                                                   + pow(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0,2)
	                                                   + pow(Err_Eff_Xibm_JpsiXim_Run1_wt/Eff_Xibm_JpsiXim_Run1_wt,2)
	                                                   + pow(trackingUnc,2)
	                                                   + pow(xiVtxUnc,2));

	StatErr_BF_Ratio_Run1_wt = BF_Ratio_Run1_wt * sqrt(pow(StatErr_N_Xib0_JpsiLambda_Run1/N_Xib0_JpsiLambda_Run1,2)
	                                                   + pow(StatErr_N_Xib0_JpsiXi0_Run1/N_Xib0_JpsiXi0_Run1,2));

	Err_BF_Ratio_Run1_wt = sqrt(pow(StatErr_BF_Ratio_Run1_wt,2) + pow(SystErr_BF_Ratio_Run1_wt,2));

	//*****Temporarily assume Run2 effs for Xib0->J/psi Lambda is same as Run1 eff, until I get Run 2 MC******
	Eff_Xib0_JpsiLambda_Run2        = Eff_Xib0_JpsiLambda_Run1;
	Eff_Xib0_JpsiLambda_Run2_wt     = Eff_Xib0_JpsiLambda_Run1_wt;
	Err_Eff_Xib0_JpsiLambda_Run2    = Err_Eff_Xib0_JpsiLambda_Run1;
	Err_Eff_Xib0_JpsiLambda_Run2_wt = Err_Eff_Xib0_JpsiLambda_Run1_wt;
	//*************************

	BF_Ratio_Run2    = (N_Xib0_JpsiLambda_Run2/Eff_Xib0_JpsiLambda_Run2) * B_Xi0_LambdaPi0 / (N_Xib0_JpsiXi0_Run2/Eff_Xibm_JpsiXim_Run2);
	BF_Ratio_Run2_wt = (N_Xib0_JpsiLambda_Run2/Eff_Xib0_JpsiLambda_Run2_wt) * B_Xi0_LambdaPi0 / (N_Xib0_JpsiXi0_Run2/Eff_Xibm_JpsiXim_Run2_wt);

	SystErr_BF_Ratio_Run2 = BF_Ratio_Run2 * sqrt(pow(SystErr_N_Xib0_JpsiLambda_Run2/N_Xib0_JpsiLambda_Run2,2)
	                                             + pow(SystErr_N_Xib0_JpsiXi0_Run2/N_Xib0_JpsiXi0_Run2,2)
	                                             + pow(Err_Eff_Xib0_JpsiLambda_Run2/Eff_Xib0_JpsiLambda_Run2,2)
	                                             + pow(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0,2)
	                                             + pow(Err_Eff_Xibm_JpsiXim_Run2/Eff_Xibm_JpsiXim_Run2,2)
	                                             + pow(trackingUnc,2)
	                                             + pow(xiVtxUnc,2));

	StatErr_BF_Ratio_Run2 = BF_Ratio_Run2 * sqrt(pow(StatErr_N_Xib0_JpsiLambda_Run2/N_Xib0_JpsiLambda_Run2,2)
	                                             + pow(StatErr_N_Xib0_JpsiXi0_Run2/N_Xib0_JpsiXi0_Run2,2));

	Err_BF_Ratio_Run2 = sqrt(pow(StatErr_BF_Ratio_Run2,2) + pow(SystErr_BF_Ratio_Run2,2));

	SystErr_BF_Ratio_Run2_wt = BF_Ratio_Run2_wt * sqrt(pow(SystErr_N_Xib0_JpsiLambda_Run2/N_Xib0_JpsiLambda_Run2,2)
	                                                   + pow(SystErr_N_Xib0_JpsiXi0_Run2/N_Xib0_JpsiXi0_Run2,2)
	                                                   + pow(Err_Eff_Xib0_JpsiLambda_Run2_wt/Eff_Xib0_JpsiLambda_Run2_wt,2)
	                                                   + pow(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0,2)
	                                                   + pow(Err_Eff_Xibm_JpsiXim_Run2_wt/Eff_Xibm_JpsiXim_Run2_wt,2)
	                                                   + pow(trackingUnc,2)
	                                                   + pow(xiVtxUnc,2));

	StatErr_BF_Ratio_Run2_wt = BF_Ratio_Run2_wt * sqrt(pow(StatErr_N_Xib0_JpsiLambda_Run2/N_Xib0_JpsiLambda_Run2,2)
	                                                   + pow(StatErr_N_Xib0_JpsiXi0_Run2/N_Xib0_JpsiXi0_Run2,2));

	Err_BF_Ratio_Run2_wt = sqrt(pow(StatErr_BF_Ratio_Run2_wt,2) + pow(SystErr_BF_Ratio_Run2_wt,2));

	cout<<"******RUN 1 RESULT***********"<<endl;
	// cout<<"Ncorr Xib0-> J/psi Lambda = "<<(N_Xib0_JpsiLambda_Run1/Eff_Xib0_JpsiLambda_Run1)<<endl;
	// cout<<"Ncorr Xib0-> J/psi Xi = "<<(N_Xib0_JpsiXi0_Run1/Eff_Xibm_JpsiXim_Run1)<<endl;

	cout<<"UNWEIGHTED RATIO = "<<BF_Ratio_Run1<<" +/- "<<StatErr_BF_Ratio_Run1<<" +/- "<<SystErr_BF_Ratio_Run1<<endl;
	cout<<"WEIGHTED RATIO   = " <<BF_Ratio_Run1_wt<<" +/- "<<StatErr_BF_Ratio_Run1_wt<<" +/- "<<SystErr_BF_Ratio_Run1_wt<<endl;

	cout<<"******RUN 2 RESULT***********"<<endl;
	// cout<<"Ncorr Xib0-> J/psi Lambda = "<<(N_Xib0_JpsiLambda_Run2/Eff_Xib0_JpsiLambda_Run2)<<endl;
	// cout<<"Ncorr Xib0-> J/psi Xi = "<<(N_Xib0_JpsiXi0_Run2/Eff_Xibm_JpsiXim_Run2)<<endl;

	cout<<"UNWEIGHTED RATIO = "<<BF_Ratio_Run2<<" +/- "<<StatErr_BF_Ratio_Run2<<" +/- "<<SystErr_BF_Ratio_Run2<<endl;
	cout<<"WEIGHTED RATIO   = " <<BF_Ratio_Run2_wt<<" +/- "<<StatErr_BF_Ratio_Run2_wt<<" +/- "<<SystErr_BF_Ratio_Run2_wt<<endl;

	//Combined Result
	Err_BF_Ratio_Comb    = (Err_BF_Ratio_Run1*Err_BF_Ratio_Run2)/sqrt(pow(Err_BF_Ratio_Run1,2)+pow(Err_BF_Ratio_Run2,2));
	Err_BF_Ratio_Comb_wt = (Err_BF_Ratio_Run1_wt*Err_BF_Ratio_Run2_wt)/sqrt(pow(Err_BF_Ratio_Run1_wt,2)+pow(Err_BF_Ratio_Run2_wt,2));

	StatErr_BF_Ratio_Comb    = (StatErr_BF_Ratio_Run1*StatErr_BF_Ratio_Run2)/sqrt(pow(StatErr_BF_Ratio_Run1,2)+pow(StatErr_BF_Ratio_Run2,2));
	StatErr_BF_Ratio_Comb_wt = (StatErr_BF_Ratio_Run1_wt*StatErr_BF_Ratio_Run2_wt)/sqrt(pow(StatErr_BF_Ratio_Run1_wt,2)+pow(StatErr_BF_Ratio_Run2_wt,2));

	SystErr_BF_Ratio_Comb    = sqrt( pow(Err_BF_Ratio_Comb,2) - pow(StatErr_BF_Ratio_Comb,2));
	SystErr_BF_Ratio_Comb_wt = sqrt( pow(Err_BF_Ratio_Comb_wt,2) - pow(StatErr_BF_Ratio_Comb_wt,2));

	BF_Ratio_Comb    = ((BF_Ratio_Run1/pow(Err_BF_Ratio_Run1,2)) + (BF_Ratio_Run2/pow(Err_BF_Ratio_Run2,2))) * pow(Err_BF_Ratio_Comb,2);
	BF_Ratio_Comb_wt = ((BF_Ratio_Run1_wt/pow(Err_BF_Ratio_Run1_wt,2)) + (BF_Ratio_Run2_wt/pow(Err_BF_Ratio_Run2_wt,2))) * pow(Err_BF_Ratio_Comb_wt,2);

	cout<<"******Combined Result************"<<endl;
	cout<<"UNWEIGHTED RATIO = ("<<BF_Ratio_Comb*100<<" +/- "<<StatErr_BF_Ratio_Comb*100<<" +/- "<<SystErr_BF_Ratio_Comb*100<<")*10^{-2}"<<endl;
	cout<<"WEIGHTED RATIO = ("<<BF_Ratio_Comb_wt*100<<" +/- "<<StatErr_BF_Ratio_Comb_wt*100<<" +/- "<<SystErr_BF_Ratio_Comb_wt*100<<")*10^{-2}"<<endl;

	cout<<"*******Systematics Breakup Run 1**********"<<endl;
	cout<<"TOTAL                :"<<(SystErr_BF_Ratio_Run1/BF_Ratio_Run1)*100<<endl;
	cout<<"N_Xib0_JpsiLambda    :"<<(SystErr_N_Xib0_JpsiLambda_Run1/N_Xib0_JpsiLambda_Run1)*100<<endl;
	cout<<"N_Xibm_JpsiXim 	    :"<<(SystErr_N_Xibm_JpsiXim_Run1/N_Xibm_JpsiXim_Run1)*100<<endl;
	cout<<"Scale_Factor 	    :"<<(Err_Scale_Factor/Scale_Factor)*100<<endl;
	cout<<"Eff_Xib0_JpsiLambda  :"<<(Err_Eff_Xib0_JpsiLambda_Run1/Eff_Xib0_JpsiLambda_Run1)*100<<endl;
	cout<<"Eff_Xibm_JpsiXim     :"<<(Err_Eff_Xibm_JpsiXim_Run1/Eff_Xibm_JpsiXim_Run1)*100<<endl;
	cout<<"B_Xi0_LambdaPi0      :"<<(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0)*100<<endl;
	cout<<"Tracking Bach.pi-    :"<<trackingUnc*100<<endl;
	cout<<"Xi Vtx. Unc.         :"<<xiVtxUnc*100<<endl;

	cout<<"*******Systematics Breakup Run 2**********"<<endl;
	cout<<"TOTAL                :"<<(SystErr_BF_Ratio_Run2/BF_Ratio_Run2)*100<<endl;
	cout<<"N_Xib0_JpsiLambda    :"<<(SystErr_N_Xib0_JpsiLambda_Run2/N_Xib0_JpsiLambda_Run2)*100<<endl;
	cout<<"N_Xibm_JpsiXim 	    :"<<(SystErr_N_Xibm_JpsiXim_Run2/N_Xibm_JpsiXim_Run2)*100<<endl;
	cout<<"Scale_Factor 	    :"<<(Err_Scale_Factor/Scale_Factor)*100<<endl;
	cout<<"Eff_Xib0_JpsiLambda  :"<<(Err_Eff_Xib0_JpsiLambda_Run2/Eff_Xib0_JpsiLambda_Run2)*100<<endl;
	cout<<"Eff_Xibm_JpsiXim     :"<<(Err_Eff_Xibm_JpsiXim_Run2/Eff_Xibm_JpsiXim_Run2)*100<<endl;
	cout<<"B_Xi0_LambdaPi0      :"<<(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0)*100<<endl;
	cout<<"Tracking Bach.pi-    :"<<trackingUnc*100<<endl;
	cout<<"Xi Vtx. Unc.         :"<<xiVtxUnc*100<<endl;

	cout<<"*******Systematics Breakup Run 1 Weighted **********"<<endl;
	cout<<"TOTAL                :"<<(SystErr_BF_Ratio_Run1_wt/BF_Ratio_Run1_wt)*100<<endl;
	cout<<"N_Xib0_JpsiLambda    :"<<(SystErr_N_Xib0_JpsiLambda_Run1/N_Xib0_JpsiLambda_Run1)*100<<endl;
	cout<<"N_Xibm_JpsiXim 	    :"<<(SystErr_N_Xibm_JpsiXim_Run1/N_Xibm_JpsiXim_Run1)*100<<endl;
	cout<<"Scale_Factor 	    :"<<(Err_Scale_Factor/Scale_Factor)*100<<endl;
	cout<<"Eff_Xib0_JpsiLambda  :"<<(Err_Eff_Xib0_JpsiLambda_Run1_wt/Eff_Xib0_JpsiLambda_Run1_wt)*100<<endl;
	cout<<"Eff_Xibm_JpsiXim     :"<<(Err_Eff_Xibm_JpsiXim_Run1_wt/Eff_Xibm_JpsiXim_Run1_wt)*100<<endl;
	cout<<"B_Xi0_LambdaPi0      :"<<(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0)*100<<endl;
	cout<<"Tracking Bach.pi-    :"<<trackingUnc*100<<endl;
	cout<<"Xi Vtx. Unc.         :"<<xiVtxUnc*100<<endl;

	cout<<"*******Systematics Breakup Run 2 Weighted **********"<<endl;
	cout<<"TOTAL                :"<<(SystErr_BF_Ratio_Run2_wt/BF_Ratio_Run2_wt)*100<<endl;
	cout<<"N_Xib0_JpsiLambda    :"<<(SystErr_N_Xib0_JpsiLambda_Run2/N_Xib0_JpsiLambda_Run2)*100<<endl;
	cout<<"N_Xibm_JpsiXim 	    :"<<(SystErr_N_Xibm_JpsiXim_Run2/N_Xibm_JpsiXim_Run2)*100<<endl;
	cout<<"Scale_Factor 	    :"<<(Err_Scale_Factor/Scale_Factor)*100<<endl;
	cout<<"Eff_Xib0_JpsiLambda  :"<<(Err_Eff_Xib0_JpsiLambda_Run2_wt/Eff_Xib0_JpsiLambda_Run2_wt)*100<<endl;
	cout<<"Eff_Xibm_JpsiXim     :"<<(Err_Eff_Xibm_JpsiXim_Run2_wt/Eff_Xibm_JpsiXim_Run2_wt)*100<<endl;
	cout<<"B_Xi0_LambdaPi0      :"<<(Err_B_Xi0_LambdaPi0/B_Xi0_LambdaPi0)*100<<endl;
	cout<<"Tracking Bach.pi-    :"<<trackingUnc*100<<endl;
	cout<<"Xi Vtx. Unc.         :"<<xiVtxUnc*100<<endl;

}
