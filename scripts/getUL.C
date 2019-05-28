#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLegend.h"
#include "TSystem.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
//#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooKeysPdf.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "TAxis.h"
#include "TH1.h"
#include "RooDataHist.h"
#include "TError.h"
#include "RooFitResult.h"
#include "RooHistPdf.h"
#include "RooDataHist.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooOneSidedProfileLL.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterPlot.h"
#include "RooStats/TestStatSampler.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/TestStatistic.h"
#include "RooHypatia2.h"
#include "RooMinimizer.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <TRolke.h>
#include "TROOT.h"
#include "TIterator.h"
#include "TLatex.h"

using namespace std;
#define Open TFile::Open

using namespace RooFit;
using namespace RooStats;

void getUL(Int_t logFlag, const char *option, Int_t config, Int_t fitType)
//fitType = 0 for nomial fit. Hypatia2 signal + Exponential bkg
//fitType = 1 for alternate fit. Double Gaussian signal + Exponential bkg
//fitType = 2 for alternate fit. Hypatia2 signal + First order order cheby background.
{
	if(logFlag)
		gSystem->RedirectOutput(Form("../logs/data/JpsiLambda/UpperLimit/config%d_tight.txt",config),"w");

	gSystem->Load("RooHypatia2_cpp.so"); //Load library for Hypatia shape
	gROOT->ProcessLine(".x lhcbStyle.C");

	Float_t bdtCut_nonZero[2] = {0.0,0.0};
	Float_t bdtCut_Zero[2]    = {0.0,0.0};

	Int_t bdtConf_nonZero[2]  = {0,0};
	Int_t bdtConf_Zero[2]     = {0,0};
	Int_t isoConf[2]          = {0,0};

	const char *isoVersion[2] = {"",""};

	Float_t xibnorm_LL[2]         = {0.0,0.0};
	Float_t xibnorm_LL_staterr[2] = {0.0,0.0};
	Float_t xibnorm_LL_systerr[2] = {0.0,0.0};
	Float_t xibnorm_LL_err[2]     = {0.0,0.0};

	Float_t xibnorm_LL_wt[2]         = {0.0,0.0};
	Float_t xibnorm_LL_staterr_wt[2] = {0.0,0.0};
	Float_t xibnorm_LL_systerr_wt[2] = {0.0,0.0};
	Float_t xibnorm_LL_err_wt[2]     = {0.0,0.0};

	// Lb->J/psi L MC eff & errs
	Int_t nGen_Lambda[2]          = {0,0}; // Generated yield
	Float_t nGen_Lambda_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_Lambda_gen[2]     = {0.0,0.0}; // Generator Eff.
	Float_t eff_Lambda_gen_err[2] = {0.0,0.0}; // Generator Eff. Stat. Err.

	Float_t eff_Lambda_rec[2]     = {0.0,0.0}; // Reco. Eff.
	Float_t eff_Lambda_rec_err[2] = {0.0,0.0}; // Reco. Eff. Stat. Err.
	Float_t eff_Lambda[2]         = {0.0,0.0}; // Overall Eff.
	Float_t eff_Lambda_systerr[2] = {0.0,0.0}; // Overall Eff. Stat. Err.

	Float_t eff_Lambda_rec_wt[2]     = {0.0,0.0}; // Reco. Eff.
	Float_t eff_Lambda_rec_err_wt[2] = {0.0,0.0}; // Reco. Eff. Stat. Err.
	Float_t eff_Lambda_wt[2]         = {0.0,0.0}; // Overall Eff.
	Float_t eff_Lambda_systerr_wt[2] = {0.0,0.0}; // Overall Eff. Stat. Err.

	// Lb->J/psi Sigma MC eff & errs
	Int_t nGen_Sigma[2]          = {0,0};
	Float_t nGen_Sigma_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_Sigma_gen[2]     = {0.0,0.0};
	Float_t eff_Sigma_gen_err[2] = {0.0,0.0};

	Float_t eff_Sigma_rec[2]     = {0.0,0.0};
	Float_t eff_Sigma_rec_err[2] = {0.0,0.0};
	Float_t eff_Sigma[2]         = {0.0,0.0};
	Float_t eff_Sigma_systerr[2] = {0.0,0.0};

	Float_t eff_Sigma_rec_wt[2]     = {0.0,0.0};
	Float_t eff_Sigma_rec_err_wt[2] = {0.0,0.0};
	Float_t eff_Sigma_wt[2]         = {0.0,0.0};
	Float_t eff_Sigma_systerr_wt[2] = {0.0,0.0};

	// eff(Lb -> J/psi Sigma) / eff(Lb -> J/psi Lambda)
	Float_t eff_ratio[2]          = {0.0,0.0};
	Float_t eff_ratio_staterr[2]  = {0.0,0.0};
	Float_t eff_ratio_systerr[2]  = {0.0,0.0};
	Float_t eff_ratio_err[2]      = {0.0,0.0};

	Float_t eff_ratio_wt[2]          = {0.0,0.0};
	Float_t eff_ratio_staterr_wt[2]  = {0.0,0.0};
	Float_t eff_ratio_systerr_wt[2]  = {0.0,0.0};
	Float_t eff_ratio_err_wt[2]      = {0.0,0.0};

	TFile *fileIn_nonZero[2];
	TFile *fileIn_Zero[2];
	TTree *treeIn_nonZero[2];
	TTree *treeIn_Zero[2];

	Int_t sigWindow_low  = 5365;
	Int_t sigWindow_high = 5600;

	Int_t bkgWindow_low  = 5700;
	Int_t bkgWindow_high = 5935;

	Int_t low = 5000;
	Int_t high = 6000;
	Int_t binwidth = 4;
	Int_t nbins    = (Int_t)(high-low)/binwidth;
	TH1D *myhist[2];
	RooDataHist *ds[2];

	Int_t Nobs[2]       = {0,0};
	Int_t Ncomb[2]      = {0,0};
	Float_t Nsig[2]     = {0.0,0.0};

	Float_t Nlb[2]              = {0.0,0.0};
	Float_t Nsig_tot[2]         = {0.0,0.0};
	Float_t Nsig_tot_STATERR[2] = {0.0,0.0};
	Float_t Nsig_tot_SYSTERR[2] = {0.0,0.0};

	Float_t Nlb_tot[2]     = {0.0,0.0};
	Float_t Nlb_tot_STATERR[2] = {0.0,0.0};
	Float_t Nxib[2]        = {0.0,0.0};
	Float_t xib_STATERR[2]  = {0.0,0.0};
	Float_t xib_SYSTERR[2]  = {0.0,0.0};

	Float_t lb_STATERR[2]   = {0.0,0.0};
	Float_t nsig_STATERR[2] = {0.0,0.0};
	Float_t nsig_SYSTERR[2] = {0.0,0.0};
	Int_t nentries[2]      = {0,0};

	RooAbsReal* xibInt[2];
	Double_t xibINT[2]= {0.0,0.0};

	RooAbsReal* lbInt[2];
	Double_t lbINT[2]= {0.0,0.0};

	RooAbsReal* sigmaInt[2];
	Double_t sigmaINT[2]= {0.0,0.0};

	//************POI****************************************************
	Float_t R[2] = {0.0,0.0};
	Float_t R_STATERR[2] = {0.0,0.0};
	Float_t R_SYSTERR[2] = {0.0,0.0};
	Float_t R_ERR[2] = {0.0,0.0};

	Float_t R_comb     = 0.0;
	Float_t R_comb_err = 0.0;
	//*******************************************************************
	if(!strncmp(option,"switch",5))
	{
		switch(config)
		{
		case 1:
		{
			isoVersion[0]   = "v0";
			isoConf[0]      = 1;
			isoVersion[1]   = "v0";
			isoConf[1]      = 1;
			bdtConf_nonZero[0] = 1;
			bdtConf_nonZero[1] = 1;
			bdtConf_Zero[0] = 1;
			bdtConf_Zero[1] = 1;
			break;
		}
		case 2:
		{
			isoVersion[0]   = "v0";
			isoConf[0]      = 1;
			isoVersion[1]   = "v0";
			isoConf[1]      = 1;
			bdtConf_nonZero[0] = 2;
			bdtConf_nonZero[1] = 2;
			bdtConf_Zero[0] = 2;
			bdtConf_Zero[1] = 2;
			break;
		}
		case 3:
		{
			isoVersion[0]   = "v0";
			isoConf[0]      = 2;
			isoVersion[1]   = "v0";
			isoConf[1]      = 2;
			bdtConf_nonZero[0] = 1;
			bdtConf_nonZero[1] = 1;
			bdtConf_Zero[0] = 1;
			bdtConf_Zero[1] = 1;
			break;
		}
		case 4:
		{
			isoVersion[0]   = "v0";
			isoConf[0]      = 2;
			isoVersion[1]   = "v0";
			isoConf[1]      = 2;
			bdtConf_nonZero[0] = 2;
			bdtConf_nonZero[1] = 2;
			bdtConf_Zero[0] = 2;
			bdtConf_Zero[1] = 2;
			break;
		}
		case 5:
		{
			isoVersion[0]   = "v1";
			isoConf[0]      = 1;
			isoVersion[1]   = "v1";
			isoConf[1]      = 1;
			bdtConf_nonZero[0] = 1;
			bdtConf_nonZero[1] = 1;
			bdtConf_Zero[0] = 1;
			bdtConf_Zero[1] = 1;
			break;
		}
		case 6:
		{
			isoVersion[0]   = "v1";
			isoConf[0]      = 1;
			isoVersion[1]   = "v1";
			isoConf[1]      = 1;
			bdtConf_nonZero[0] = 2;
			bdtConf_nonZero[1] = 2;
			bdtConf_Zero[0] = 2;
			bdtConf_Zero[1] = 2;
			break;
		}
		case 7:
		{
			isoVersion[0]   = "v1";
			isoConf[0]      = 2;
			isoVersion[1]   = "v1";
			isoConf[1]      = 2;
			bdtConf_nonZero[0] = 1;
			bdtConf_nonZero[1] = 1;
			bdtConf_Zero[0] = 1;
			bdtConf_Zero[1] = 1;
			break;
		}
		case 8:
		{
			isoVersion[0]   = "v1";
			isoConf[0]      = 2;
			isoVersion[1]   = "v1";
			isoConf[1]      = 2;
			bdtConf_nonZero[0] = 2;
			bdtConf_nonZero[1] = 2;
			bdtConf_Zero[0] = 2;
			bdtConf_Zero[1] = 2;
			break;
		}
		}
		bdtCut_nonZero[0] = 0.475;//0.375 - 0.1;
		bdtCut_nonZero[1] = 0.555;//0.535 - 0.1;

		bdtCut_Zero[0] = 0.365;//0.285;
		bdtCut_Zero[1] = 0.455;//0.415;
	}
	else if(!strncmp(option,"best",4)) //Set parameters for best fit
	{
		isoVersion[0] = "v0";//"v1";
		isoVersion[1] = "v0";//"v0";

		isoConf[0] = 2;//1;
		isoConf[1] = 1;//2;

		bdtConf_nonZero[0] = 2;//2;
		bdtConf_nonZero[1] = 2;//1;

		bdtConf_Zero[0] = 2;//1;
		bdtConf_Zero[1] = 2;//1;

		bdtCut_nonZero[0] = 0.475;//0.375 - 0.1;
		bdtCut_nonZero[1] = 0.555;//0.535 - 0.1;

		bdtCut_Zero[0] = 0.365;//0.285;
		bdtCut_Zero[1] = 0.495;//0.455;
	}

	//******Systematics that are set by hand now*************************
	Float_t xib_syst             = 0.05;
	Float_t eff_ratio_syst       = 0.02;

	// ************************Master Workspace**************************
	RooWorkspace w("w");

	//*********MASTER VARIABLE*******************************************
	w.factory(Form("Lb_DTF_M_JpsiLConstr[%d,%d]",4000,7000));

	RooRealVar *myVar = w.var("Lb_DTF_M_JpsiLConstr");

	myVar->setRange("signal_window",sigWindow_low,sigWindow_high);
	//********Get data*******************************
	const char *dataPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda";
	for(Int_t run=1; run<=2; run++)
	{
		Int_t i = run-1;
		fileIn_nonZero[i] = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root",
		                              dataPath,run),"READ");
		treeIn_nonZero[i] = (TTree*)fileIn_nonZero[i]->Get("MyTuple");
		fileIn_Zero[i] = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_ZeroTracks_noPID.root",
		                           dataPath,run),"READ");
		treeIn_Zero[i] = (TTree*)fileIn_Zero[i]->Get("MyTuple");
		treeIn_nonZero[i]->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                            dataPath,run,bdtConf_nonZero[i],
		                                            isoConf[i],isoVersion[i]));
		treeIn_Zero[i]->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_zeroTracksLL_FinalBDT%d_noPID.root",
		                                         dataPath,run,bdtConf_Zero[i]));

		Nobs[i] = treeIn_nonZero[i]->GetEntries(Form("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d && BDT%d > %f",sigWindow_low,sigWindow_high,bdtConf_nonZero[i],bdtCut_nonZero[i]))
		          + treeIn_Zero[i]->GetEntries(Form("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d && BDT%d > %f",sigWindow_low,sigWindow_high,bdtConf_Zero[i],bdtCut_Zero[i]));

		Ncomb[i] = treeIn_nonZero[i]->GetEntries(Form("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d && BDT%d > %f",bkgWindow_low,bkgWindow_high,bdtConf_nonZero[i],bdtCut_nonZero[i]))
		           + treeIn_Zero[i]->GetEntries(Form("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d && BDT%d > %f",bkgWindow_low,bkgWindow_high,bdtConf_Zero[i],bdtCut_Zero[i]));

		treeIn_nonZero[i]->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_nonzero%d(%d,%d,%d)",run,nbins,low,high),
		                        Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");

		treeIn_Zero[i]->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_zero%d(%d,%d,%d)",run,nbins,low,high),
		                     Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");

		TH1D *myhist_nonzero = (TH1D*)gDirectory->Get(Form("myhist_nonzero%d",run));
		TH1D *myhist_zero = (TH1D*)gDirectory->Get(Form("myhist_zero%d",run));

		myhist[i] = new TH1D(Form("myhist%d",run),"",nbins,low,high);

		myhist[i]->Add(myhist_zero,myhist_nonzero);

		nentries[i] = myhist[i]->Integral();
		cout<<"Run "<<run<<" nentries = "<<nentries[i]<<endl;

		ds[i] = new RooDataHist(Form("ds%d",run),Form("ds%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),myhist[i]);

		cout<<"Done making RooDataHist"<<endl;
		(ds[i])->Print();

		w.import(*(ds[i]));
	}
	//************************************************

	//*****Run script to get Xib normalization***************************
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	cout<<"Get Xib Normalization"<<endl;
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		gSystem->Exec(Form("python -c \'from GetXibNorm import GetNorm;"
		                   " GetNorm(%d, \"%s\", %d, %d, %d, %f, %f, 0)\'",
		                   run, isoVersion[i], isoConf[i], bdtConf_nonZero[i],
		                   bdtConf_Zero[i], bdtCut_nonZero[i], bdtCut_Zero[i]));

		ifstream infile(Form("../logs/mc/JpsiXi/run%d/xibNorm_log.txt",run));

		infile>>xibnorm_LL[i];         // Get normalization
		infile>>xibnorm_LL_staterr[i]; // Get statistical error on norm. This comes from fit error
		infile>>xibnorm_LL_systerr[i]; // Get systematic error on norm.

		infile>>xibnorm_LL_wt[i];         // Get weighted normalization
		infile>>xibnorm_LL_staterr_wt[i]; // Get weighted statistical error on norm. This comes from the fit error
		infile>>xibnorm_LL_systerr_wt[i]; // Get weighted systematic error on norm.

		// xibnorm_LL_systerr[i] = xibnorm_LL[i]*xib_syst;
		xibnorm_LL_err[i]     = sqrt(pow(xibnorm_LL_staterr[i],2) + pow(xibnorm_LL_systerr[i],2)); //Combining stat and syst in quadrature

		// xibnorm_LL_systerr_wt[i] = xibnorm_LL_wt[i]*xib_syst;
		xibnorm_LL_err_wt[i]     = sqrt(pow(xibnorm_LL_staterr_wt[i],2) + pow(xibnorm_LL_systerr_wt[i],2)); //Combining stat and syst in quadrature

		cout<<"************************************************"<<endl;
		cout<<"The UNWEIGHTED LL Xib normalization for Run "<<run
		    <<" is "<<xibnorm_LL[i]<<" +/- "<<xibnorm_LL_staterr[i]
		    <<" +/- "<<xibnorm_LL_systerr[i]<<endl;

		cout<<"The WEIGHTED LL Xib normalization for Run "<<run
		    <<" is "<<xibnorm_LL_wt[i]<<" +/- "<<xibnorm_LL_staterr_wt[i]
		    <<" +/- "<<xibnorm_LL_systerr_wt[i]<<endl;
		cout<<"************************************************"<<endl;
	}
	//************************************************

	//****Get J/psi Lambda efficiencies from MC*******
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	cout<<"Get Lb -> J/psi Lambda efficiency"<<endl;
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	const char* lambdaMCPath = "../rootFiles/mcFiles/JpsiLambda/JpsiLambda";

	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_Lambda = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root",
		                                           lambdaMCPath,run));
		TTree *mcTreeIn_nonZero_Lambda = (TTree*)mcFileIn_nonZero_Lambda->Get("MyTuple");

		TFile *mcFileIn_Zero_Lambda    = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_ZeroTracks_noPID.root",
		                                           lambdaMCPath,run));
		TTree *mcTreeIn_Zero_Lambda    = (TTree*)mcFileIn_Zero_Lambda->Get("MyTuple");

		mcTreeIn_nonZero_Lambda->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                                  lambdaMCPath,run,bdtConf_nonZero[i],
		                                                  isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_Lambda->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_zeroTracksLL_FinalBDT%d_noPID.root",
		                                               lambdaMCPath,run,bdtConf_Zero[i]));

		fstream genFile_Lambda;
		genFile_Lambda.open((Form("../logs/mc/JpsiLambda/JpsiLambda/run%d/gen_log.txt",run)));

		genFile_Lambda>>nGen_Lambda[i]; //Get number of generated events

		TFile *genWtsFile_Lambda = nullptr;
		if(run == 1)
		{
			genWtsFile_Lambda = Open(Form("%s/run%d/RW/gbWeights_gen_new.root",
			                              lambdaMCPath,run));
		}
		else if(run == 2)
		{
			genWtsFile_Lambda = Open(Form("%s/run%d/RW/gbWeights_gen.root",
			                              lambdaMCPath,run));
		}
		TTree *genWtsTree_Lambda = (TTree*)genWtsFile_Lambda->Get("MyTuple");

		genWtsTree_Lambda->AddFriend("MyTuple",Form("%s/run%d/RW/tauWeights_gen.root",
		                                            lambdaMCPath,run));
		if(run == 1)
		{
			genWtsTree_Lambda->Draw("gb_wts_new*wt_tau>>genWt_Lambda","","goff");
		}
		else if(run == 2)
		{
			genWtsTree_Lambda->Draw("gb_wts*wt_tau>>genWt_Lambda","","goff");
		}

		TH1F *genWt_Lambda = (TH1F*)gDirectory->Get("genWt_Lambda");
		nGen_Lambda_wt[i] = genWt_Lambda->GetEntries()*genWt_Lambda->GetMean();

		fstream genEffFile_Lambda;
		genEffFile_Lambda.open(Form("../logs/mc/JpsiLambda/JpsiLambda/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_Lambda>>eff_Lambda_gen[i];    //Get generator efficiency
		genEffFile_Lambda>>eff_Lambda_gen_err[i];//and error on above

		cout<<"Run "<<run<<" Lambda Generator Effs = "<<eff_Lambda_gen[i]*100
		    <<" % +/- "<<eff_Lambda_gen_err[i]*100<<" %"<<endl;

		if(run == 1)
		{
			mcTreeIn_nonZero_Lambda->Draw("gb_wts_new*wt_tau>>wt_lambda_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_Lambda->Draw("gb_wts_new*wt_tau>>wt_lambda_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}
		else if(run == 2)
		{
			mcTreeIn_nonZero_Lambda->Draw("gb_wts*wt_tau>>wt_lambda_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_Lambda->Draw("gb_wts*wt_tau>>wt_lambda_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}

		TH1F *wt_lambda_nonZero = (TH1F*)gDirectory->Get("wt_lambda_nonZero");
		TH1F *wt_lambda_Zero = (TH1F*)gDirectory->Get("wt_lambda_Zero");

		Int_t num_Lambda_wt = (wt_lambda_nonZero->GetMean()*wt_lambda_nonZero->GetEntries()) +
		                      (wt_lambda_Zero->GetMean()*wt_lambda_Zero->GetEntries());

		Int_t num_Lambda = mcTreeIn_nonZero_Lambda->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]))
		                   + mcTreeIn_Zero_Lambda->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i])); //NOTE NO TM HERE

		eff_Lambda_rec[i]     = num_Lambda*1.0/nGen_Lambda[i]; //Calc. reco eff.
		eff_Lambda_rec_err[i] = sqrt(eff_Lambda_rec[i]*(1-eff_Lambda_rec[i])/nGen_Lambda[i]); //statistical error on recon. eff.

		eff_Lambda_rec_wt[i]     = num_Lambda_wt*1.0/nGen_Lambda_wt[i]; //Calc. weighted reco eff.
		eff_Lambda_rec_err_wt[i] = sqrt(eff_Lambda_rec_wt[i]*(1-eff_Lambda_rec_wt[i])/nGen_Lambda_wt[i]); //statistical error on weighted recon. eff.

		cout<<"Run "<<run<<" UNWEIGHTED Lambda Recons. Effs = "<<eff_Lambda_rec[i]*100
		    <<" % +/- "<<eff_Lambda_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED Lambda Recons. Effs = "<<eff_Lambda_rec_wt[i]*100
		    <<" % +/- "<<eff_Lambda_rec_err_wt[i]*100<<" %"<<endl;

		eff_Lambda[i]     = eff_Lambda_rec[i]*eff_Lambda_gen[i]; // Calc. total eff.
		eff_Lambda_systerr[i] = eff_Lambda[i]*sqrt(pow((eff_Lambda_gen_err[i]/eff_Lambda_gen[i]),2) +
		                                           pow((eff_Lambda_rec_err[i]/eff_Lambda_rec[i]),2));                                                                                                                               // and stat error on tot. eff.

		eff_Lambda_wt[i]     = eff_Lambda_rec_wt[i]*eff_Lambda_gen[i];                                                                                                                                                                          // Calc. total eff.
		eff_Lambda_systerr_wt[i] = eff_Lambda_wt[i]*sqrt(pow((eff_Lambda_gen_err[i]/eff_Lambda_gen[i]),2) +
		                                                 pow((eff_Lambda_rec_err_wt[i]/eff_Lambda_rec_wt[i]),2));                                                                                                                                                                                                                                                                                                                                                                                                                                                 // and stat error on tot. eff.

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi Lambda Eff = "<<eff_Lambda[i]*100
		    <<" % +/- "<<eff_Lambda_systerr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi Lambda Eff = "<<eff_Lambda_wt[i]*100
		    <<" % +/- "<<eff_Lambda_systerr_wt[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;
	}
	//*******************************************************************

	RooHistPdf* SIG[2];
	RooKeysPdf* SIG_KEYS[2];
	RooDataSet* ds_sig[2];
	RooDataSet* ds_sig_wt[2];

	//****Get J/psi Sigma efficiencies and shape from MC*****************
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	cout<<"Get Lb -> J/psi Sigma efficiency and shape"<<endl;
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;

	const char* sigmaPath = "/data1/avenkate/JpsiLambda_RESTART/"
	                        "rootFiles/mcFiles/JpsiLambda/JpsiSigma";

	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_Sigma = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_nonZeroTracks_noPID.root",
		                                          sigmaPath,run));
		TTree *mcTreeIn_nonZero_Sigma = (TTree*)mcFileIn_nonZero_Sigma->Get("MyTuple");

		TFile *mcFileIn_Zero_Sigma    = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_ZeroTracks_noPID.root",
		                                          sigmaPath,run));
		TTree *mcTreeIn_Zero_Sigma    = (TTree*)mcFileIn_Zero_Sigma->Get("MyTuple");

		mcTreeIn_nonZero_Sigma->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                                 sigmaPath,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_Sigma->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_zeroTracksLL_FinalBDT%d_noPID.root",
		                                              sigmaPath,run,bdtConf_Zero[i]));

		//***************Efficiency***************************************
		fstream genFile_Sigma;
		genFile_Sigma.open((Form("../logs/mc/JpsiLambda/JpsiSigma/run%d/gen_log.txt",run)));

		genFile_Sigma>>nGen_Sigma[i]; // Get number of generated events

		TFile *genWtsFile_Sigma = nullptr;

		if(run == 1)
		{
			genWtsFile_Sigma = Open(Form("%s/run%d/RW/gbWeights_gen_new.root",
			                             sigmaPath,run));
		}
		else if(run == 2)
		{
			genWtsFile_Sigma = Open(Form("%s/run%d/RW/gbWeights_gen.root",
			                             sigmaPath,run));
		}
		TTree *genWtsTree_Sigma = (TTree*)genWtsFile_Sigma->Get("MyTuple");

		genWtsTree_Sigma->AddFriend("MyTuple",Form("%s/run%d/RW/tauWeights_gen.root",
		                                           sigmaPath,run));

		if(run == 1)
		{
			genWtsTree_Sigma->Draw("gb_wts_new*wt_tau>>genWt_Sigma","","goff");
		}
		else if(run == 2)
		{
			genWtsTree_Sigma->Draw("gb_wts*wt_tau>>genWt_Sigma","","goff");
		}
		TH1F *genWt_Sigma = (TH1F*)gDirectory->Get("genWt_Sigma");

		nGen_Sigma_wt[i] = genWt_Sigma->GetEntries()*genWt_Sigma->GetMean();

		fstream genEffFile_Sigma;
		genEffFile_Sigma.open(Form("../logs/mc/JpsiLambda/JpsiSigma/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_Sigma>>eff_Sigma_gen[i]; // Get generator efficiency
		genEffFile_Sigma>>eff_Sigma_gen_err[i]; // and error on above

		cout<<"Run "<<run<<" Sigma Generator Effs = "<<eff_Sigma_gen[i]*100
		    <<" % +/- "<<eff_Sigma_gen_err[i]*100<<" %"<<endl;

		if(run == 1)
		{
			mcTreeIn_nonZero_Sigma->Draw("gb_wts_new*wt_tau>>wt_Sigma_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_Sigma->Draw("gb_wts_new*wt_tau>>wt_Sigma_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}
		else if(run == 2)
		{
			mcTreeIn_nonZero_Sigma->Draw("gb_wts*wt_tau>>wt_Sigma_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_Sigma->Draw("gb_wts*wt_tau>>wt_Sigma_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}
		TH1F *wt_Sigma_nonZero = (TH1F*)gDirectory->Get("wt_Sigma_nonZero");
		TH1F *wt_Sigma_Zero = (TH1F*)gDirectory->Get("wt_Sigma_Zero");

		Int_t num_Sigma_wt = (wt_Sigma_nonZero->GetMean()*wt_Sigma_nonZero->GetEntries()) +
		                     (wt_Sigma_Zero->GetMean()*wt_Sigma_Zero->GetEntries());

		Int_t num_Sigma = mcTreeIn_nonZero_Sigma->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i])) +
		                  mcTreeIn_Zero_Sigma->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));

		eff_Sigma_rec[i]     = num_Sigma*1.0/nGen_Sigma[i]; // Calc. reco. eff.
		eff_Sigma_rec_err[i] = sqrt(eff_Sigma_rec[i]*(1-eff_Sigma_rec[i])/nGen_Sigma[i]); //stat error on recon. eff.

		eff_Sigma_rec_wt[i]     = num_Sigma_wt*1.0/nGen_Sigma_wt[i]; //Calc. weighted reco eff.
		eff_Sigma_rec_err_wt[i] = sqrt(eff_Sigma_rec_wt[i]*(1-eff_Sigma_rec_wt[i])/nGen_Sigma_wt[i]); //statistical error on weighted recon. eff.

		cout<<"Run "<<run<<" Sigma Recons. Effs = "
		    <<eff_Sigma_rec[i]*100<<" % +/- "<<eff_Sigma_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED Sigma Recons. Effs = "<<eff_Sigma_rec_wt[i]*100
		    <<" % +/- "<<eff_Sigma_rec_err_wt[i]*100<<" %"<<endl;

		eff_Sigma[i] = eff_Sigma_gen[i] * eff_Sigma_rec[i]; // Calc overall eff.
		eff_Sigma_systerr[i] = eff_Sigma[i]*sqrt(pow((eff_Sigma_gen_err[i]/eff_Sigma_gen[i]),2) +
		                                         pow((eff_Sigma_rec_err[i]/eff_Sigma_rec[i]),2));                                                                                                                         // and stat. error on above

		eff_Sigma_wt[i]     = eff_Sigma_rec_wt[i]*eff_Sigma_gen[i];                                                                                                                                                                          // Calc. total eff.
		eff_Sigma_systerr_wt[i] = eff_Sigma_wt[i]*sqrt(pow((eff_Sigma_gen_err[i]/eff_Sigma_gen[i]),2) +
		                                               pow((eff_Sigma_rec_err_wt[i]/eff_Sigma_rec_wt[i]),2));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             // and stat error on tot. eff.

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi Sigma Eff = "<<eff_Sigma[i]*100
		    <<" % +/- "<<eff_Sigma_systerr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi Sigma Eff = "<<eff_Sigma_wt[i]*100
		    <<" % +/- "<<eff_Sigma_systerr_wt[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;

		eff_ratio[i]         = eff_Sigma[i]/eff_Lambda[i]; // Calc eff ratio.
		eff_ratio_staterr[i] = eff_ratio[i]*sqrt(pow((eff_Sigma_systerr[i]/eff_Sigma[i]),2)+pow((eff_Lambda_systerr[i]/eff_Lambda[i]),2)); // stat err on ratio
		eff_ratio_systerr[i] = eff_ratio[i]*eff_ratio_syst;
		eff_ratio_err[i]     = sqrt(pow(eff_ratio_staterr[i],2) + pow(eff_ratio_systerr[i],2));//combine in quadrature

		eff_ratio_wt[i]         = eff_Sigma_wt[i]/eff_Lambda_wt[i]; // Calc eff ratio.
		eff_ratio_staterr_wt[i] = eff_ratio_wt[i]*sqrt(pow((eff_Sigma_systerr_wt[i]/eff_Sigma_wt[i]),2)+pow((eff_Lambda_systerr_wt[i]/eff_Lambda_wt[i]),2)); // stat err on ratio
		eff_ratio_systerr_wt[i] = eff_ratio_wt[i]*eff_ratio_syst;
		eff_ratio_err_wt[i]     = sqrt(pow(eff_ratio_staterr_wt[i],2) + pow(eff_ratio_systerr_wt[i],2));//combine in quadrature

		cout<<"***************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Sigma/Lambda efficiency ratio = "<<eff_ratio[i]<<" +/- "<<eff_ratio_err[i]<<endl;
		cout<<"Run "<<run<<" WEIGHTED Sigma/Lambda efficiency ratio   = "<<eff_ratio_wt[i]<<" +/- "<<eff_ratio_err_wt[i]<<endl;
		cout<<"***************************************"<<endl;

		//******************************************************************

		//****************Shape*********************************************

		mcTreeIn_Zero_Sigma->SetBranchStatus("*",0);
		mcTreeIn_Zero_Sigma->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_Zero_Sigma->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
		mcTreeIn_Zero_Sigma->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_Zero_Sigma->SetBranchStatus("gb_wts",1);
		mcTreeIn_Zero_Sigma->SetBranchStatus("wt_tau",1);
		if(run == 1)
			mcTreeIn_Zero_Sigma->SetBranchStatus("gb_wts_new",1);

		mcTreeIn_nonZero_Sigma->SetBranchStatus("*",0);
		mcTreeIn_nonZero_Sigma->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_nonZero_Sigma->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		mcTreeIn_nonZero_Sigma->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_nonZero_Sigma->SetBranchStatus("gb_wts",1);
		mcTreeIn_nonZero_Sigma->SetBranchStatus("wt_tau",1);
		if(run == 1)
			mcTreeIn_nonZero_Sigma->SetBranchStatus("gb_wts_new",1);
		TFile *tempFile = new TFile("tempFile_sig.root","RECREATE");

		TTree* mcTreeIn_Zero_Sigma_cut    = (TTree*)mcTreeIn_Zero_Sigma->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));//Not TRUTH MATCHING HERE!
		TTree* mcTreeIn_nonZero_Sigma_cut = (TTree*)mcTreeIn_nonZero_Sigma->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));//Not TRUTH MATCHING HERE!

		TList *list_sig = new TList;
		list_sig->Add(mcTreeIn_Zero_Sigma_cut);
		list_sig->Add(mcTreeIn_nonZero_Sigma_cut);

		TTree *combTree_sig = TTree::MergeTrees(list_sig);
		combTree_sig->SetName("combTree_sig");

		RooRealVar *gbWtVar = nullptr;

		if(run == 1)
		{
			gbWtVar  = new RooRealVar("gb_wts_new","gb Weight Var",-100.,100.);
		}
		else if(run == 2)
		{
			gbWtVar  = new RooRealVar("gb_wts","gb Weight Var",-100.,100.);
		}

		RooRealVar *tauWtVar = new RooRealVar("wt_tau","tau Weight Var",-100.,100.);

		ds_sig[i] = new RooDataSet("ds_sig","ds_sig",combTree_sig,RooArgSet(*myVar,*gbWtVar,*tauWtVar));
		ds_sig[i]->Print();

		RooFormulaVar *totWt = new RooFormulaVar("totWt","@0*@1",RooArgList(*gbWtVar,*tauWtVar));
		RooRealVar *totWt_var = (RooRealVar*)ds_sig[i]->addColumn(*totWt);

		ds_sig_wt[i] = new RooDataSet("ds_sig_wt","ds_sig_wt",RooArgSet(*myVar,*totWt_var),Import(*(ds_sig[i])),WeightVar(*totWt_var));
		ds_sig_wt[i]->Print();

		SIG_KEYS[i] = new RooKeysPdf(Form("SIG%d",run),Form("SIG%d",run),*myVar,*(ds_sig_wt[i]),RooKeysPdf::MirrorBoth,1);

		RooPlot *framesigma = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		framesigma->SetTitle("J/#psi #Sigma");
		// ds_sigma->plotOn(framesigma,Name("sigmadata"));
		ds_sig[i]->plotOn(framesigma,Name("sigmadata_nowt"),LineColor(kGreen));
		ds_sig_wt[i]->plotOn(framesigma,Name("sigmadata"),LineColor(kBlack));
		// sigmashape.plotOn(framesigma,Name("sigmafit"),LineColor(kBlue));
		(*(SIG_KEYS[i])).plotOn(framesigma,Name("sigmafitsmooth"),LineColor(kRed),LineStyle(kDashed));

		TCanvas *csigma = new TCanvas(Form("JpsiSigma%d",run),Form("JpsiSigma%d",run));
		framesigma->Draw();
		w.import(*(SIG_KEYS[i]));

		cout<<"Done importing Jpsi Sigma shape"<<endl;

		sigmaInt[i] = SIG_KEYS[i]->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
		sigmaINT[i] = sigmaInt[i]->getValV();

		// cout<<"Run "<<run<<" sigma fraction = "<<sigmaINT[i]<<endl;
	}
	//********************************************************************


	RooDataSet* ds_xi[2];
	RooDataSet* ds_xi_wt[2];
	RooKeysPdf* XIB_KEYS[2];
	//******************Get shape from Xib background********************
	// RooRealVar xibmass("Lb_DTF_M_JpsiLConstr","xibmass",5200.,5740.);
	const char* xibPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/JpsiXi";

	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		TFile *filein_xi_nonZero = Open(Form("%s/run%d/jpsixi_cutoutks_LL_nonZeroTracks_noPID.root",xibPath,run));
		TTree *treein_xi_nonZero = (TTree*)filein_xi_nonZero->Get("MyTuple");

		TFile *filein_xi_Zero = Open(Form("%s/run%d/jpsixi_cutoutks_LL_ZeroTracks_noPID.root",xibPath,run));
		TTree *treein_xi_Zero = (TTree*)filein_xi_Zero->Get("MyTuple");

		treein_xi_nonZero->AddFriend("MyTuple",Form("%s/run%d/jpsixi_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                            xibPath,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		treein_xi_Zero->AddFriend("MyTuple",Form("%s/run%d/jpsixi_zeroTracksLL_FinalBDT%d_noPID.root",
		                                         xibPath,run,bdtConf_Zero[i]));
		treein_xi_Zero->SetBranchStatus("*",0);
		treein_xi_Zero->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		treein_xi_Zero->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
		treein_xi_Zero->SetBranchStatus("Lb_BKGCAT",1);
		treein_xi_Zero->SetBranchStatus("gb_wts",1);
		if(run == 1)
		{
			treein_xi_Zero->SetBranchStatus("gb_wts_new",1);
		}

		treein_xi_nonZero->SetBranchStatus("*",0);
		treein_xi_nonZero->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		treein_xi_nonZero->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		treein_xi_nonZero->SetBranchStatus("Lb_BKGCAT",1);
		treein_xi_nonZero->SetBranchStatus("gb_wts",1);
		if(run == 1)
		{
			treein_xi_Zero->SetBranchStatus("gb_wts_new",1);
		}

		TFile *tempFile = new TFile("tempFile.root","RECREATE");

		TTree* treein_xi_Zero_cut    = (TTree*)treein_xi_Zero->CopyTree(Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));
		TTree* treein_xi_nonZero_cut = (TTree*)treein_xi_nonZero->CopyTree(Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));

		TList *list = new TList;
		list->Add(treein_xi_Zero_cut);
		list->Add(treein_xi_nonZero_cut);

		TTree *combTree = TTree::MergeTrees(list);
		combTree->SetName("combTree");

		RooRealVar *gbWtVar = nullptr;

		if(run == 1)
		{
			gbWtVar  = new RooRealVar("gb_wts","gb Weight Var",-100.,100.);
		}
		else if(run == 2)
		{
			gbWtVar  = new RooRealVar("gb_wts","gb Weight Var",-100.,100.);
		}
		ds_xi[i] = new RooDataSet("ds_xi","ds_xi",combTree,RooArgSet(*myVar,*gbWtVar));
		ds_xi[i]->Print();

		ds_xi_wt[i] = new RooDataSet("ds_xi_wt","ds_xi_wt",RooArgSet(*myVar,*gbWtVar),Import(*(ds_xi[i])),WeightVar(*gbWtVar));
		ds_xi_wt[i]->Print();
		// treein_xi_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_nonZero%d(%d,%d,%d)",run,nbins,low,high),
		//                         Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");//TRUTH MATCHING HERE
		// treein_xi_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_Zero%d(%d,%d,%d)",run,nbins,low,high),
		//                      Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");//TRUTH MATCHING HERE
		//
		// TH1D *hxib_nonZero = (TH1D*)gDirectory->Get(Form("hxib_nonZero%d",run));
		// TH1D *hxib_Zero    = (TH1D*)gDirectory->Get(Form("hxib_Zero%d",run));
		// TH1D *hxib         = new TH1D(Form("hxib%d",run),"",nbins,low,high);
		// hxib->Add(hxib_nonZero,hxib_Zero);
		//
		// TH1D *hxib_smooth  = (TH1D*)hxib->Clone(Form("hxib_smooth%d",run));
		// hxib_smooth->Smooth(40);//TODO TWEAK THIS
		//
		// RooDataHist *ds_xib        = new RooDataHist(Form("ds_xib%d",run),Form("ds_xib%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),hxib);
		// RooDataHist *ds_xib_smooth = new RooDataHist(Form("ds_xib_smooth%d",run),Form("ds_xib_smooth%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),hxib_smooth);
		//
		// ds_xib_smooth->Print();

		XIB_KEYS[i] = new RooKeysPdf(Form("XIB%d",run),Form("XIB%d",run),*myVar,*(ds_xi_wt[i]),RooKeysPdf::NoMirror);
		// XIB[i] = new RooHistPdf(Form("XIB%d",run),Form("XIB%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_xib_smooth,0);

		RooPlot *framexib = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		framexib->SetTitle("J/#psi #Xi");
		ds_xi[i]->plotOn(framexib,Name("xibdata_nowt"),LineColor(kGreen));
		ds_xi_wt[i]->plotOn(framexib,Name("xibdata"),LineColor(kBlack));
		// ds_xib->plotOn(framexib,Name("xibdata"));
		// (*(XIB[i])).plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));
		(*(XIB_KEYS[i])).plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));

		TCanvas *cxib = new TCanvas(Form("JpsiXi%d",run),Form("JpsiXi%d",run));
		framexib->Draw();

		// tempFile->Close();
		// w.import(*(XIB[i]));
		w.import(*(XIB_KEYS[i]));

		cout<<"Done importing Xib shape"<<endl;

		xibInt[i] = XIB_KEYS[i]->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
		xibINT[i] = xibInt[i]->getValV();

		Nxib[i]     = xibnorm_LL_wt[i]*xibINT[i];
		xib_STATERR[i] = xibnorm_LL_staterr_wt[i]*xibINT[i];
		xib_SYSTERR[i] = xibnorm_LL_systerr_wt[i]*xibINT[i];

		cout<<"************************************************"<<endl;
		cout<<"The WEIGHTED Xib normalization inside signal region for Run "<<run
		    <<" is "<<Nxib[i]<<" +/- "<<xib_STATERR[i]
		    <<" +/- "<<xib_SYSTERR[i]<<endl;
		cout<<"************************************************"<<endl;

	}
	//*********Hypatia signal shape for Lambda_b0************************

	w.factory("RooHypatia2::Lb_Run1(Lb_DTF_M_JpsiLConstr,lambda_Run1[-2.0,-4.0,0.0],0,0,"
	          "sigma_Run1[10.,1.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[1.7,1.0,3.0],"
	          "2 ,a2_Run1[3.0,1.0,4.0], 2)");

	w.factory("RooHypatia2::Lb_Run2(Lb_DTF_M_JpsiLConstr,lambda_Run2[-2.5,-4.0,0.0],0,0,"
	          "sigma_Run2[10.,1.,20.], mean_Run2[5619.6,5619,5621], a1_Run2[1.5,1.0,3.0],"
	          "2 ,a2_Run2[1.5,1.0,3.0], 2)");

	w.factory("Gaussian::Lb_Run1_Gaus1(Lb_DTF_M_JpsiLConstr, mean_Run1_Gaus[5619.6,5619,5621],"
	          "sigma_Run1_Gaus1[5.,1.,20.])");
	w.factory("Gaussian::Lb_Run1_Gaus2(Lb_DTF_M_JpsiLConstr, mean_Run1_Gaus,"
	          "sigma_Run1_Gaus2[10.,1.,20.])");
	w.factory("Gaussian::Lb_Run2_Gaus1(Lb_DTF_M_JpsiLConstr, mean_Run2_Gaus[5619.6,5619,5621],"
	          "sigma_Run2_Gaus1[5.,1.,20.])");
	w.factory("Gaussian::Lb_Run2_Gaus2(Lb_DTF_M_JpsiLConstr, mean_Run2_Gaus,"
	          "sigma_Run2_Gaus2[10.,1.,20.])");

	cout<<"Done defining J/psi Lambda Hypatia shapes"<<endl;
	//*******************************************************************

	//*********Continuum backgruond****************
	w.factory("Exponential::Bkg_Run1(Lb_DTF_M_JpsiLConstr,expr('-1*pow(10,slope_Run1)',slope_Run1[-6.0,-11.0,-2.0]))");
	w.factory("Exponential::Bkg_Run2(Lb_DTF_M_JpsiLConstr,expr('-1*pow(10,slope_Run2)',slope_Run2[-6.0,-11.0,-2.0]))");

	w.var("slope_Run1")->setError(0.25);
	w.var("slope_Run2")->setError(0.25);

	w.factory("Chebychev::Bkg_Run1_Cheby(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0]})");
	w.factory("Chebychev::Bkg_Run2_Cheby(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-0.1,0.1]})");

	cout<<"Done defining cont. bkg shapes"<<endl;
	//*******************************************************************

	w.factory("nLb_Run1[5600,2000,7000]");
	w.factory("nLb_Run2[20000,10000,22000]");

	w.factory(Form("nBkg_Run1[200,1,%d]",nentries[0]));
	w.factory(Form("nBkg_Run2[1000,1,%d]",nentries[1]));

	w.factory("SUM:model1(nLb_Run1*Lb_Run1 , nBkg_Run1*Bkg_Run1)");
	w.factory("SUM:model2(nLb_Run2*Lb_Run2 , nBkg_Run2*Bkg_Run2)");

	w.factory("SUM:model1_Cheby(nLb_Run1*Lb_Run1 , nBkg_Run1*Bkg_Run1_Cheby)");
	w.factory("SUM:model2_Cheby(nLb_Run2*Lb_Run2 , nBkg_Run2*Bkg_Run2_Cheby)");

	w.factory("SUM:Lb_Run1_Gaus(f_Run1[0.5,0.0,1.0]*Lb_Run1_Gaus1 , Lb_Run1_Gaus2)");
	w.factory("SUM:Lb_Run2_Gaus(f_Run2[0.5,0.0,1.0]*Lb_Run2_Gaus1 , Lb_Run2_Gaus2)");

	w.factory("SUM:model1_Gaus(nLb_Run1*Lb_Run1_Gaus , nBkg_Run1*Bkg_Run1)");
	w.factory("SUM:model2_Gaus(nLb_Run2*Lb_Run2_Gaus , nBkg_Run2*Bkg_Run2)");

	RooAbsPdf* model1 = w.pdf("model1"); // get the model
	RooAbsPdf* model2 = w.pdf("model2"); // get the model

	RooAbsPdf* model1_Cheby = w.pdf("model1_Cheby"); // get the model
	RooAbsPdf* model2_Cheby = w.pdf("model2_Cheby"); // get the model

	RooAbsPdf* model1_Gaus = w.pdf("model1_Gaus"); // get the model
	RooAbsPdf* model2_Gaus = w.pdf("model2_Gaus"); // get the model

	//***********************MAKE COMBINED DATASET************************************
	RooCategory sample("sample","sample");
	sample.defineType("run1");
	sample.defineType("run2");

	RooAbsData *combData;
	combData = new RooDataHist("combData","combined data",*(w.var("Lb_DTF_M_JpsiLConstr")),
	                           Index(sample),Import("run1",*myhist[0]),
	                           Import("run2",*myhist[1]));
	combData->Print();
	RooSimultaneous simPdf("simPdf","simultaneous pdf",sample);
	simPdf.addPdf(*model1,"run1");
	simPdf.addPdf(*model2,"run2");
	w.import(*combData);
	w.import(simPdf);

	RooSimultaneous simPdf_Cheby("simPdf_Cheby","simultaneous pdf",sample);
	simPdf_Cheby.addPdf(*model1_Cheby,"run1");
	simPdf_Cheby.addPdf(*model2_Cheby,"run2");
	w.import(simPdf_Cheby);

	RooSimultaneous simPdf_Gaus("simPdf_Gaus","simultaneous Gaussian pdf",sample);
	simPdf_Gaus.addPdf(*model1_Gaus,"run1");
	simPdf_Gaus.addPdf(*model2_Gaus,"run2");
	w.import(simPdf_Gaus);

	//************************DO THE FIT***********************
	w.var("Lb_DTF_M_JpsiLConstr")->setRange("ref",5500,5800);
	RooSimultaneous *fitPdf = nullptr;
	const char *sigPdf_Run1 = "", *bkgPdf_Run1 = "";
	const char *sigPdf_Run2 = "", *bkgPdf_Run2 = "";

	if(fitType == 0)
	{
		fitPdf = &simPdf;
		sigPdf_Run1 = "Lb_Run1";
		bkgPdf_Run1 = "Bkg_Run1";
		sigPdf_Run2 = "Lb_Run2";
		bkgPdf_Run2 = "Bkg_Run2";
	}
	else if(fitType == 1)
	{
		fitPdf = &simPdf_Gaus;
		sigPdf_Run1 = "Lb_Run1_Gaus";
		bkgPdf_Run1 = "Bkg_Run1";
		sigPdf_Run2 = "Lb_Run2_Gaus";
		bkgPdf_Run2 = "Bkg_Run2";
	}
	else if(fitType == 2)
	{
		fitPdf = &simPdf_Cheby;
		sigPdf_Run1 = "Lb_Run1";
		bkgPdf_Run1 = "Bkg_Run1_Cheby";
		sigPdf_Run2 = "Lb_Run2";
		bkgPdf_Run2 = "Bkg_Run2_Cheby";
	}
	RooFitResult *res = fitPdf->fitTo(*combData,Extended(), Save(), Hesse(false), Strategy(1), SumCoefRange("ref"),Range("ref"));
	//*******************************************************************

	//*********************PLOTTING STUFF*********************************************

	TLatex *myLatex = new TLatex();
	myLatex->SetTextFont(42);
	myLatex->SetTextColor(1);
	myLatex->SetTextAlign(12);
	myLatex->SetNDC(kTRUE);

	TCanvas* c_run1 = new TCanvas("Run1","Run1", 1200, 800);

	RooPlot *frame_run1 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	frame_run1->SetTitle("Run1 Fit");
	frame_run1->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame_run1->GetYaxis()->SetTitle("Candidates/(4 MeV/#it{c}^{2})");

	combData->plotOn(frame_run1,Name("data_Run1"),Cut("sample==sample::run1"),DataError(RooAbsData::Poisson));
	fitPdf->plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Name("fit_Run1"));
	fitPdf->plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf(sigPdf_Run1))),Name("lb_Run1"),LineColor(kMagenta+2));
	// fitPdf->plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb1_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	// fitPdf->plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb2_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	fitPdf->plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf(bkgPdf_Run1))),LineColor(kRed),Name("bkg_Run1"));

	frame_run1->GetYaxis()->SetRangeUser(0,20);
	frame_run1->GetXaxis()->SetRangeUser(5300,5900);
	RooArgSet *allpar_run1 = fitPdf->getParameters(*(ds[0]));
	RooArgSet *floatpar_run1 = (RooArgSet*)allpar_run1->selectByAttrib("Constant",kFALSE);
	floatpar_run1->Print();
	int floatpars_run1 = (floatpar_run1->selectByAttrib("Constant",kFALSE))->getSize() - 1;//-1 because sample also gets included in this list
	cout<<"run1 float pars = "<<floatpars_run1<<endl;
	Double_t chi2_run1 = frame_run1->chiSquare("fit_Run1","data_Run1",floatpars_run1);
	cout<<"chi square2/dof = "<<chi2_run1<<endl;

	Int_t fit_ndof_run1 = nbins - floatpars_run1;
	cout<<"chi square2 = "<<chi2_run1*fit_ndof_run1<<endl;

	///////////
	// TPad *pad1 = new TPad("pad1","pad1",0.0,0.2,1.0,1.0);
	// TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.2);
	//
	// pad1->SetGridx();
	// pad1->SetGridy();
	// pad2->SetGridx();
	// pad2->SetGridy();
	//
	// pad1->SetBottomMargin(0.04);
	// pad2->SetTopMargin(0);
	// pad2->SetBottomMargin(0.3);
	// pad2->SetBorderMode(0);
	// pad1->SetBorderMode(0);
	// c_run1->SetBorderMode(0);
	// pad2->Draw();
	// pad1->Draw();
	// pad1->cd();
	// gPad->SetTopMargin(0.06);
	// pad1->Update();

	frame_run1->Draw();

	// TLatex l_run1;
	// l_run1.SetTextSize(0.025);
	// l_run1.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run1));

	// c_run1->Modified();

	auto legend_run1 = new TLegend(0.65,0.7,0.85,0.9);
	legend_run1->SetTextSize(0.045);
	legend_run1->AddEntry("data_Run1","Data","lp");
	legend_run1->AddEntry("fit_Run1","Total Fit","l");
	legend_run1->AddEntry("lb_Run1","J/#psi #Lambda shape","l");
	legend_run1->AddEntry("bkg_Run1","Comb. Bkg. shape","l");
	legend_run1->Draw("same");

	myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");

	c_run1->Update();

	// // Pull distribution
	// RooPlot *frame_run1x2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	// // RooPlot *framex2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), low,5800,100);
	// RooHist* hpull_run1 = frame_run1->pullHist("data_Run1","fit_Run1");
	// frame_run1x2->addPlotable(hpull_run1,"P");
	// hpull_run1->SetLineColor(kBlack);
	// hpull_run1->SetMarkerColor(kBlack);
	// frame_run1x2->SetTitle(0);
	// frame_run1x2->GetYaxis()->SetTitle("Pull");
	// frame_run1x2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	// frame_run1x2->GetYaxis()->SetTitleSize(0.15);
	// frame_run1x2->GetYaxis()->SetLabelSize(0.15);
	// frame_run1x2->GetXaxis()->SetTitleSize(0.15);
	// frame_run1x2->GetXaxis()->SetLabelSize(0.15);
	// frame_run1x2->GetYaxis()->CenterTitle();
	// frame_run1x2->GetYaxis()->SetTitleOffset(0.25);
	// frame_run1x2->GetXaxis()->SetTitleOffset(0.75);
	// frame_run1x2->GetYaxis()->SetNdivisions(505);
	// frame_run1x2->GetYaxis()->SetRangeUser(-5.0,5.0);
	// pad2->cd();
	// frame_run1x2->Draw();

	// c_run1->cd();
	// pad1->cd();

	// cout<<"Run1 Pull Mean Y = "<<hpull_run1->GetMean(2)<<endl;
	// cout<<"Run1 Pull RMS Y = "<<hpull_run1->GetRMS(2)<<endl;

	TCanvas* c_run2 = new TCanvas("Run2","Run2", 1200, 800);

	RooPlot *frame_run2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);

	frame_run2->SetTitle("Run2 Fit");
	frame_run2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame_run2->GetYaxis()->SetTitle("Candidates/(4 MeV/#it{c}^{2})");

	combData->plotOn(frame_run2,Name("data_Run2"),Cut("sample==sample::run2"),DataError(RooAbsData::Poisson));
	fitPdf->plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Name("fit_Run2"));
	fitPdf->plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf(sigPdf_Run2))),Name("lb_Run2"),LineColor(kMagenta+2));
	// fitPdf->plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb1_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	// fitPdf->plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb2_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	fitPdf->plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf(bkgPdf_Run2))),LineColor(kRed),Name("bkg_Run2"));

	frame_run2->GetYaxis()->SetRangeUser(0,60);
	frame_run2->GetXaxis()->SetRangeUser(5300,5900);
	// Double_t chiSquare1 = frame_run2->chiSquare("fit_run2","data_run2");
	// cout<<"chi square1/dof = "<<chiSquare1<<endl;
	RooArgSet *floatpar_run2 = fitPdf->getParameters(*(ds[1]));
	floatpar_run2->Print();
	int floatpars_run2 = (floatpar_run2->selectByAttrib("Constant",kFALSE))->getSize() -1;//-1 because sample also gets included in this list
	cout<<"run2 float pars = "<<floatpars_run2<<endl;
	Double_t chi2_run2 = frame_run2->chiSquare("fit_Run2","data_Run2",floatpars_run2);
	cout<<"chi square2/dof = "<<chi2_run2<<endl;

	Int_t fit_ndof_run2 = nbins - floatpars_run2;
	cout<<"chi square2 = "<<chi2_run2*fit_ndof_run2<<endl;

	///////////
	// TPad *pad3 = new TPad("pad3","pad3",0.0,0.2,1.0,1.0);
	// TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,0.2);
	//
	// pad3->SetGridx();
	// pad3->SetGridy();
	// pad4->SetGridx();
	// pad4->SetGridy();
	//
	// pad3->SetBottomMargin(0.04);
	// pad4->SetTopMargin(0);
	// pad4->SetBottomMargin(0.3);
	// pad4->SetBorderMode(0);
	// pad3->SetBorderMode(0);
	// c_run1->SetBorderMode(0);
	// pad4->Draw();
	// pad3->Draw();
	// pad3->cd();
	// gPad->SetTopMargin(0.06);
	// pad3->Update();

	frame_run2->Draw();
	// c_run2->cd();

	auto legend_run2 = new TLegend(0.65,0.7,0.85,0.9);
	legend_run2->SetTextSize(0.045);
	legend_run2->AddEntry("data_Run2","Data","lp");
	legend_run2->AddEntry("fit_Run2","Total Fit","l");
	legend_run2->AddEntry("lb_Run2","J/#psi #Lambda shape","l");
	legend_run2->AddEntry("bkg_Run2","Comb. Bkg. shape","l");
	legend_run2->Draw("same");

	myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");

	// TLatex l_run2;
	// l_run2.SetTextSize(0.025);
	// l_run2.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run2));

	c_run2->Update();

	// c_run1->SaveAs("../plots/ANA/finalFit_run1_zoomY.pdf");
	// c_run2->SaveAs("../plots/ANA/finalFit_run2_zoomY.pdf");


	// // Pull distribution
	// RooPlot *frame_run2x2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	// // RooPlot *framex2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), low,5800,100);
	// RooHist* hpull_run2 = frame_run2->pullHist("data_Run2","fit_Run2");
	// frame_run2x2->addPlotable(hpull_run2,"P");
	// hpull_run2->SetLineColor(kBlack);
	// hpull_run2->SetMarkerColor(kBlack);
	// frame_run2x2->SetTitle(0);
	// frame_run2x2->GetYaxis()->SetTitle("Pull");
	//
	// frame_run2x2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	// frame_run2x2->GetYaxis()->SetTitleSize(0.15);
	// frame_run2x2->GetYaxis()->SetLabelSize(0.15);
	// frame_run2x2->GetXaxis()->SetTitleSize(0.15);
	// frame_run2x2->GetXaxis()->SetLabelSize(0.15);
	// frame_run2x2->GetYaxis()->CenterTitle();
	// frame_run2x2->GetYaxis()->SetTitleOffset(0.25);
	// frame_run2x2->GetXaxis()->SetTitleOffset(0.75);
	// frame_run2x2->GetYaxis()->SetNdivisions(505);
	// frame_run2x2->GetYaxis()->SetRangeUser(-5.0,5.0);
	// pad4->cd();
	// frame_run2x2->Draw();
	//
	// c_run2->cd();

	Double_t chi2_global = chi2_run1*fit_ndof_run1 + chi2_run2*fit_ndof_run2;
	Int_t ndof_global = fit_ndof_run1 + fit_ndof_run2;

	Double_t chi2_ndof_global = chi2_global/ndof_global;

	cout<<"****************************"<<endl;
	cout<<"Global Fit chi2/dof = "<<chi2_ndof_global<<endl;
	cout<<"****************************"<<endl;

	for(Int_t i=0; i<=1; i++)
	{
		const char *myvar = Form("nLb_Run%d",i+1);
		const char *mystr = nullptr;

		if(fitType == 0 || fitType == 2)
		{
			mystr = Form("Lb_Run%d",i+1);
		}
		else if(fitType == 1)
		{
			mystr = Form("Lb_Run%d_Gaus",i+1);
		}
		lbInt[i]       = w.pdf(mystr)->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
		lbINT[i]       = lbInt[i]->getValV();
		Nlb[i]         = lbINT[i]*(w.var(myvar)->getVal());
		lb_STATERR[i]   = (lbInt[i]->getValV())*(w.var(myvar)->getError());
		Nlb_tot[i]     = w.var(myvar)->getVal();
		Nlb_tot_STATERR[i] = w.var(myvar)->getError();

		Nsig[i]        = Nobs[i] - Nlb[i] - Nxib[i] - Ncomb[i];
		nsig_STATERR[i] = sqrt( Nobs[i] + Ncomb[i] + pow(lb_STATERR[i],2) + pow(xib_STATERR[i],2));
		nsig_SYSTERR[i] = sqrt( pow(xib_SYSTERR[i],2));
		Nsig_tot[i]     = Nsig[i]/sigmaINT[i];
		Nsig_tot_STATERR[i] = nsig_STATERR[i]/sigmaINT[i];
		Nsig_tot_SYSTERR[i] = nsig_SYSTERR[i]/sigmaINT[i];

		R[i] = (Nsig_tot[i]/Nlb_tot[i])/eff_ratio_wt[i];
		R_STATERR[i] = R[i]*sqrt( pow(Nsig_tot_STATERR[i]/Nsig_tot[i],2) +
		                          pow(Nlb_tot_STATERR[i]/Nlb_tot[i],2));

		R_SYSTERR[i] = R[i]*sqrt( pow(Nsig_tot_SYSTERR[i]/Nsig_tot[i], 2) +
		                          pow(eff_ratio_err_wt[i]/eff_ratio_wt[i],2) );
		R_ERR[i] = sqrt( pow(R_STATERR[i],2) + pow(R_SYSTERR[i],2) );

		cout<<"*********RUN "<<i+1<<"***********"<<endl;
		cout<<"Nobs    = "<<Nobs[i]<<" +/- "<<sqrt((float)Nobs[i])<<endl;
		cout<<"Ncomb   = "<<Ncomb[i]<<" +/- "<<sqrt((float)Ncomb[i])<<endl;
		cout<<"Nxib    = "<<Nxib[i]<<" +/- "<<xib_STATERR[i]<<" +/- "<<xib_SYSTERR[i]<<endl;
		cout<<"NLb     = "<<Nlb[i]<<" +/- "<<lb_STATERR[i]<<endl;
		cout<<"_____________________________________"<<endl;
		cout<<"Nsig    = "<<Nsig[i]<<" +/- "<<nsig_STATERR[i]<<" +/- "<<nsig_SYSTERR[i]<<endl;
		cout<<"Frac    = "<<sigmaINT[i]<<endl;
		cout<<"NLb_tot = "<<Nlb_tot[i]<<" +/- "<<Nlb_tot_STATERR[i]<<endl;
		cout<<"eff_JpsiLambda = "<<eff_Lambda_wt[i]*100<<" % +/- "
		    <<eff_Lambda_systerr_wt[i]*100<<" %"<<endl;
		cout<<"eff_JpsiSigma = "<<eff_Sigma_wt[i]*100<<" % +/- "
		    <<eff_Sigma_systerr_wt[i]*100<<" %"<<endl;
		cout<<"Eff Rat = "<<eff_ratio_wt[i]<<" +/- "<<eff_ratio_err_wt[i]<<endl;
		cout<<"_____________________________________"<<endl;
		cout<<"R       = "<<R[i]<<" +/- "<<R_STATERR[i]<<" +/- "<<R_SYSTERR[i]<<endl;
		cout<<"90% CL Upper Limit = "<<R[i]+(1.28*R_ERR[i])<<endl;
		cout<<"*********************************"<<endl;

		cout<<endl;
	}

	R_comb_err = 1/sqrt(pow(1/R_ERR[0],2) + pow(1/R_ERR[1],2));
	R_comb = (R[0]/pow(R_ERR[0],2) + R[1]/pow(R_ERR[1],2)) * pow(R_comb_err,2);
	cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	cout<<"COMBINED R = "<<R_comb<<" +/- "<<R_comb_err<<endl;
	cout<<"COMBINED 90% CL Upper Limit = "<<R_comb+(1.28*R_comb_err)<<endl;
	cout<<"COMBINED 95% CL Upper Limit = "<<R_comb+(1.65*R_comb_err)<<endl;
	cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"<<endl;
	if(logFlag)
		gSystem->RedirectOutput(0);
}
