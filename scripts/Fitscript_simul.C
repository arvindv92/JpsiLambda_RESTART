#include "Fitscript_simul.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

#define Open TFile::Open

void Fitscript_simul(const char *option, Int_t myLow, Int_t myHigh, Int_t Lst1405_rwtype, Int_t bkgType, Float_t bdtCut, const char* fileName)
//option = "best" to fit to data with best BDT cuts
//myLow and myHigh define the fit range
//rwType 0 is no RW, 1 is MV RW, 2 is BONN RW
//bkgType = 0 for Exponential. 1 for 2nd order Chebychev. 2 for 3rd order Chebychev
{
	Bool_t calcUL   = true;
	Bool_t isBinned = true; //set to false if you want unbinned ML fit.

	// Fit params
	// If final fit is a binned fit ,it will use binwidth for binning (MeV)
	// If it is an unbinned it, it will just use the binning for visualization
	// Int_t myLow      = myLow, myHigh = myHigh; //Define range in which fit is performed
	Int_t binwidth = 4;
	Int_t nbins    = (Int_t)(myHigh-myLow)/binwidth;

	Int_t sigWindow_low  = 5365;
	Int_t sigWindow_high = 5600;


	// gSystem->RedirectOutput(Form("../logs/data/JpsiLambda/UpperLimit/Fit_HypatiaSig_ExpoBkg_%d_%d_%dMeVBins.txt",
	//                              myLow,myHigh,binwidth),"w");

	// gSystem->RedirectOutput("tempLog.txt","a");

	gSystem->Exec("date");
	gSystem->Load("RooHypatia2_cpp.so"); //Load library for Hypatia shape

	// Tail params for Lb Hypatia
	// Float_t a1[2]  = {1.98,2.68};
	// Float_t a2[2]  = {1.89,3.00};
	// Float_t lam[2] = {-2.32,-2.01};
	//
	// Float_t a1_err[2]  = {0.25,0.22};
	// Float_t a2_err[2]  = {0.29,0.33};
	// Float_t lam_err[2] = {0.19,0.07};

	Float_t bdtCut_nonZero[2] = {0.0,0.0};
	Float_t bdtCut_Zero[2]    = {0.0,0.0};

	Int_t bdtConf_nonZero[2]  = {0,0};
	Int_t bdtConf_Zero[2]     = {0,0};
	Int_t isoConf[2]          = {0,0};

	const char *isoVersion[2] = {"",""};

	if(!strncmp(option,"best",4)) //Set parameters for best fit
	{
		isoVersion[0] = "v0";//"v1";
		isoVersion[1] = "v0";//"v0";

		isoConf[0] = 2;//1;
		isoConf[1] = 1;//2;

		bdtConf_nonZero[0] = 2;//2;
		bdtConf_nonZero[1] = 2;//1;

		bdtConf_Zero[0] = 2;//1;
		bdtConf_Zero[1] = 2;//1;

		bdtCut_nonZero[0] = 0.475-bdtCut;
		bdtCut_nonZero[1] = 0.555-bdtCut;

		bdtCut_Zero[0] = 0.365;
		bdtCut_Zero[1] = 0.495;
	}

	// Xib normalization & errs

	Float_t trackingErr = 0.05; // rel error. on Xib Norm includes 4.5% tracking unc, and 2% uncertainty for material effects from hadronic interactions
	Float_t xiVtxUnc    = 0.014; // rel error. on Xib Norm Unc from vertexing the Xi-. Comes from Steves ANA

	Float_t XibNorm[2]         = {0.0,0.0}; // Unweighted normalization inside signal window
	Float_t XibNorm_StatErr[2] = {0.0,0.0}; // Abs. error
	Float_t XibNorm_SystErr[2] = {0.0,0.0}; // Abs. error

	Float_t XibNorm_wt[2]                = {0.0,0.0}; //Weighted normalization inside signal window
	Float_t XibNorm_wt_StatErr[2]        = {0.0,0.0}; // Abs. error
	Float_t XibNorm_wt_SystErr[2]        = {0.0,0.0}; // Abs. error

	// Float_t XibNorm[2]         = {0.0,0.0};
	// Float_t XibNorm_StatErr[2] = {0.0,0.0};
	// Float_t XibNorm_systerr[2] = {0.0,0.0};
	// Float_t XibNorm_err[2]     = {0.0,0.0};
	//
	// Float_t XibNorm_wt[2]         = {0.0,0.0};
	// Float_t XibNorm_wt_StatErr[2] = {0.0,0.0};
	// Float_t XibNorm_systerr_wt[2] = {0.0,0.0};
	// Float_t XibNorm_err_wt[2]     = {0.0,0.0};
	//
	// Lb->J/psi L MC eff & errs
	Int_t nGen_Lambda[2]          = {0,0}; // Generated yield
	Float_t nGen_Lambda_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_Lambda_gen[2]     = {0.0,0.0}; // Generator Eff.
	Float_t eff_Lambda_gen_err[2] = {0.0,0.0}; // Generator Eff. Stat. Err.

	Float_t eff_Lambda_rec[2]     = {0.0,0.0}; // Reco. Eff.
	Float_t eff_Lambda_rec_err[2] = {0.0,0.0}; // Reco. Eff. Stat. Err.
	// Float_t eff_Lambda[2]         = {0.0,0.0}; // Overall Eff.
	// Float_t eff_Lambda_staterr[2] = {0.0,0.0}; // Overall Eff. Stat. Err.

	Float_t eff_Lambda_rec_wt[2]     = {0.0,0.0}; // Reco. Eff.
	Float_t eff_Lambda_rec_err_wt[2] = {0.0,0.0}; // Reco. Eff. Stat. Err.
	// Float_t eff_Lambda_wt[2]         = {0.0,0.0}; // Overall Eff.
	// Float_t eff_Lambda_staterr_wt[2] = {0.0,0.0}; // Overall Eff. Stat. Err.

	// Lb->J/psi 1405 MC eff & errs
	Int_t nGen_1405[2]          = {0,0};
	Float_t nGen_1405_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_1405_gen[2]     = {0.0,0.0};
	Float_t eff_1405_gen_err[2] = {0.0,0.0};

	Float_t eff_1405_rec[2]      = {0.0,0.0};
	Float_t eff_1405_rec_err[2]  = {0.0,0.0};
	Float_t eff_1405[2]          = {0.0,0.0};
	Float_t eff_1405_SystErr[2]  = {0.0,0.0};

	Float_t eff_1405_wt_rec[2]      = {0.0,0.0};
	Float_t eff_1405_wt_rec_err[2]  = {0.0,0.0};
	Float_t eff_1405_wt[2]          = {0.0,0.0};
	Float_t eff_1405_wt_SystErr[2]  = {0.0,0.0};

	// Lb->J/psi 1520 MC eff & errs
	Int_t nGen_1520[2]          = {0,0};
	Float_t nGen_1520_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_1520_gen[2]     = {0.0,0.0};
	Float_t eff_1520_gen_err[2] = {0.0,0.0};

	Float_t eff_1520_rec[2]      = {0.0,0.0};
	Float_t eff_1520_rec_err[2]  = {0.0,0.0};
	Float_t eff_1520[2]          = {0.0,0.0};
	Float_t eff_1520_SystErr[2]  = {0.0,0.0};

	Float_t eff_1520_rec_wt[2]      = {0.0,0.0};
	Float_t eff_1520_wt_rec_err[2]  = {0.0,0.0};
	Float_t eff_1520_wt[2]          = {0.0,0.0};
	Float_t eff_1520_wt_SystErr[2]  = {0.0,0.0};

	// Lb->J/psi 1600 MC eff & errs
	Int_t nGen_1600[2]          = {0,0};
	Float_t nGen_1600_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_1600_gen[2]     = {0.0,0.0};
	Float_t eff_1600_gen_err[2] = {0.0,0.0};

	Float_t eff_1600_rec[2]      = {0.0,0.0};
	Float_t eff_1600_rec_err[2]  = {0.0,0.0};
	Float_t eff_1600[2]          = {0.0,0.0};
	Float_t eff_1600_SystErr[2]  = {0.0,0.0};

	Float_t eff_1600_rec_wt[2]      = {0.0,0.0};
	Float_t eff_1600_wt_rec_err[2]  = {0.0,0.0};
	Float_t eff_1600_wt[2]          = {0.0,0.0};
	Float_t eff_1600_wt_SystErr[2]  = {0.0,0.0};

	// // Lb->chiC1 Lambda MC eff & errs
	// Int_t nGen_chic1[2]          = {0,0};
	// Float_t nGen_chic1_wt[2]          = {0,0};
	// Float_t nGen_chiC1_wt[2]     = {0,0}; // Generated yield weighted
	// Float_t eff_chic1_gen[2]     = {0.0,0.0};
	// Float_t eff_chic1_gen_err[2] = {0.0,0.0};
	//
	// Float_t eff_chic1_rec[2]      = {0.0,0.0};
	// Float_t eff_chic1_rec_err[2]  = {0.0,0.0};
	// Float_t eff_chic1[2]          = {0.0,0.0};
	// Float_t eff_chic1_staterr[2]  = {0.0,0.0};
	//
	// Float_t eff_chic1_rec_wt[2]      = {0.0,0.0};
	// Float_t eff_chic1_rec_err_wt[2]  = {0.0,0.0};
	// Float_t eff_chic1_wt[2]          = {0.0,0.0};
	// Float_t eff_chic1_staterr_wt[2]  = {0.0,0.0};

	// Lb->J/psi Sigma MC eff & errs
	Int_t nGen_Sigma[2]          = {0,0};
	Float_t nGen_Sigma_wt[2]     = {0,0}; // Generated yield weighted
	Float_t eff_Sigma_gen[2]     = {0.0,0.0};
	Float_t eff_Sigma_gen_err[2] = {0.0,0.0};

	Float_t eff_Sigma_rec[2]     = {0.0,0.0};
	Float_t eff_Sigma_rec_err[2] = {0.0,0.0};
	// Float_t eff_Sigma[2]         = {0.0,0.0};
	// Float_t eff_Sigma_staterr[2] = {0.0,0.0};

	Float_t eff_Sigma_rec_wt[2]     = {0.0,0.0};
	Float_t eff_Sigma_rec_err_wt[2] = {0.0,0.0};
	// Float_t eff_Sigma_wt[2]         = {0.0,0.0};
	// Float_t eff_Sigma_staterr_wt[2] = {0.0,0.0};

	// eff(Lb -> J/psi Sigma) / eff(Lb -> J/psi Lambda)
	// Float_t eff_ratio[2]          = {0.0,0.0};
	// Float_t eff_ratio_StatErr[2]  = {0.0,0.0};
	// Float_t eff_ratio_SystErr[2]  = {0.0,0.0};
	// Float_t eff_ratio_Err[2]      = {0.0,0.0};
	//
	// Float_t eff_ratio_wt[2]          = {0.0,0.0};
	// Float_t eff_ratio_wt_StatErr[2]  = {0.0,0.0};
	// Float_t eff_ratio_wt_SystErr[2]  = {0.0,0.0};
	// Float_t eff_ratio_Err_wt[2]      = {0.0,0.0};

	// eff(Lb -> J/psi LAMBDA(1405)) / eff(Lb -> J/psi Lambda)
	Float_t eff_ratio_1405[2]          = {0.0,0.0};
	Float_t eff_ratio_StatErr_1405[2]  = {0.0,0.0};
	Float_t eff_ratio_SystErr_1405[2]  = {0.0,0.0};
	Float_t eff_ratio_Err_1405[2]      = {0.0,0.0};

	Float_t eff_ratio_wt_1405[2]          = {0.0,0.0};
	Float_t eff_ratio_wt_StatErr_1405[2]  = {0.0,0.0};
	Float_t eff_ratio_wt_SystErr_1405[2]  = {0.0,0.0};
	Float_t eff_ratio_Err_1405_wt[2]      = {0.0,0.0};

	// eff(Lb -> J/psi LAMBDA(1520)) / eff(Lb -> J/psi Lambda)
	Float_t eff_ratio_1520[2]          = {0.0,0.0};
	Float_t eff_ratio_StatErr_1520[2]  = {0.0,0.0};
	Float_t eff_ratio_SystErr_1520[2]  = {0.0,0.0};
	Float_t eff_ratio_Err_1520[2]      = {0.0,0.0};

	Float_t eff_ratio_wt_1520[2]          = {0.0,0.0};
	Float_t eff_ratio_wt_StatErr_1520[2]  = {0.0,0.0};
	Float_t eff_ratio_wt_SystErr_1520[2]  = {0.0,0.0};
	Float_t eff_ratio_Err_1520_wt[2]      = {0.0,0.0};

	// eff(Lb -> J/psi LAMBDA(1600)) / eff(Lb -> J/psi Lambda)
	Float_t eff_ratio_1600[2]          = {0.0,0.0};
	Float_t eff_ratio_StatErr_1600[2]  = {0.0,0.0};
	Float_t eff_ratio_SystErr_1600[2]  = {0.0,0.0};
	Float_t eff_ratio_Err_1600[2]      = {0.0,0.0};

	Float_t eff_ratio_wt_1600[2]          = {0.0,0.0};
	Float_t eff_ratio_wt_StatErr_1600[2]  = {0.0,0.0};
	Float_t eff_ratio_wt_SystErr_1600[2]  = {0.0,0.0};
	Float_t eff_ratio_Err_1600_wt[2]      = {0.0,0.0};

	// // eff(Lb -> chiC1 Lambda) / eff(Lb -> J/psi Lambda)
	// Float_t eff_ratio_chic1[2]          = {0.0,0.0};
	// Float_t eff_ratio_StatErr_chic1[2]  = {0.0,0.0};
	// Float_t eff_ratio_SystErr_chic1[2]  = {0.0,0.0};
	// Float_t eff_ratio_Err_chic1[2]      = {0.0,0.0};
	//
	// Float_t eff_ratio_chic1_wt[2]          = {0.0,0.0};
	// Float_t eff_ratio_StatErr_chic1_wt[2]  = {0.0,0.0};
	// Float_t eff_ratio_SystErr_chic1_wt[2]  = {0.0,0.0};
	// Float_t eff_ratio_Err_chic1_wt[2]      = {0.0,0.0};

	Float_t N_JpsiLambda[2]         = {0.0,0.0};
	Float_t N_JpsiLambda_StatErr[2] = {0.0,0.0}; // abs. error. = fit error on yield
	Float_t N_JpsiLambda_RelSyst[2] = {0.9,0.8}; // relative % unc. comes from choice of fit model
	Float_t N_JpsiLambda_SystErr[2] = {0.0,0.0}; // abs. error

	Float_t N_JpsiSigma[2] = {0.0,0.0};
	Float_t N_JpsiSigma_StatErr[2] = {0.0,0.0};
	Float_t N_JpsiSigma_SystErr[2] = {0.0,0.0};

	Float_t N_JpsiSigma_wt[2] = {0.0,0.0};
	Float_t N_JpsiSigma_wt_StatErr[2] = {0.0,0.0};
	Float_t N_JpsiSigma_wt_SystErr[2] = {0.0,0.0};

	Float_t window_JpsiLambda[2] = {0.0,0.0};
	Float_t window_JpsiSigma[2]  = {0.0,0.0};
	Float_t window_JpsiXi[2]     = {0.0,0.0};
	Float_t fitwindow_JpsiXi[2]  = {0.0,0.0};

	Float_t N_JpsiLambda_window[2] = {0.0,0.0};
	Float_t N_JpsiLambda_window_StatErr[2] = {0.0,0.0};
	Float_t N_JpsiLambda_window_SystErr[2] = {0.0,0.0};

	Float_t Nobs[2]       = {0.0,0.0}; // no. of events observed in data inside signal window
	Float_t Nobs_StatErr[2] = {0.0,0.0};

	Float_t Ncomb[2]      = {0.0,0.0}; // no. of events counted in sideband window above peak in data
	Float_t Ncomb_StatErr[2] = {0.0,0.0};

	Float_t eff_JpsiLambda[2]         = {0.0,0.0}; //% Unweighted efficiency
	Float_t eff_JpsiLambda_SystErr[2] = {0.0,0.0}; //% Binomial uncertainty

	Float_t eff_JpsiSigma[2]         = {0.0,0.0}; //% Unweighted efficiency
	Float_t eff_JpsiSigma_SystErr[2] = {0.0,0.0}; //% Binomial uncertainty

	Float_t eff_JpsiLambda_wt[2]         = {0.0,0.0}; //% Weighted efficiency
	Float_t eff_JpsiLambda_wt_SystErr[2] = {0.0,0.0}; //% Binomial uncertainty

	Float_t eff_JpsiSigma_wt[2]         = {0.0,0.0}; //% Weighted efficiency
	Float_t eff_JpsiSigma_wt_SystErr[2] = {0.0,0.0}; //% Binomial uncertainty

	Float_t eff_ratio[2]         = {0.0,0.0}; // Ratio of unweighted efficiencies
	Float_t eff_ratio_StatErr[2] = {0.0,0.0}; // Absolute error. Will remain zero
	Float_t eff_ratio_SystErr[2] = {0.0,0.0}; // Comes from binomial uncertainties

	Float_t eff_ratio_wt[2]         = {0.0,0.0}; // Ratio of Weighted efficiencies
	Float_t eff_ratio_wt_StatErr[2] = {0.0,0.0}; // Absolute error. Will remain zero
	Float_t eff_ratio_wt_SystErr[2] = {0.0,0.0}; // Comes from binomial uncertainties


	// ***************Systematics put in by hand ************************
	// Float_t xib_syst             = 0.05;
	// Float_t eff_ratio_syst       = 0.02;
	// Float_t eff_ratio_syst_1405  = 0.02;
	// Float_t eff_ratio_syst_1520  = 0.02;
	// Float_t eff_ratio_syst_1600  = 0.02;
	// Float_t eff_ratio_syst_chic1 = 0.02;
	// ***************************Flags**********************************
	// Flags controlling shapes
	Int_t lst1405flag    = 1;
	Int_t lst1520flag    = 1;
	Int_t lst1600flag    = 1;
	Int_t lst1810flag    = 0;
	// Int_t chic1flag      = 0;
	Int_t xibflag        = 1;
	Int_t sigmaflag      = 1;

	// ************************Master Workspace**************************
	RooWorkspace w("w");

	//*********MASTER VARIABLE*******************************************
	w.factory(Form("Lb_DTF_M_JpsiLConstr[%d,%d]",myLow,myHigh));

	RooRealVar *myVar = w.var("Lb_DTF_M_JpsiLConstr");

	myVar->setRange("signal_window",sigWindow_low,sigWindow_high); //this define the signal window
	myVar->setRange("fit_window",5500,5800); // this defines the fit window for the lambda_b fit

	w.var("Lb_DTF_M_JpsiLConstr")->setBins((Int_t)(myHigh-myLow)/binwidth);
	//*******************************************************************

	//*********CONTROL VERBOSITY*****************************************
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
	//*******************************************************************

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

		infile>>XibNorm[i];         // Get normalization
		infile>>XibNorm_StatErr[i]; // Get statistical error on norm. This is the fit error
		infile>>XibNorm_SystErr[i]; // Get systematic error on norm.

		infile>>XibNorm_wt[i];         // Get weighted normalization
		infile>>XibNorm_wt_StatErr[i]; // Get weighted statistical error on norm. This is the fit error
		infile>>XibNorm_wt_SystErr[i]; // Get weighted systematic error on norm.
		//
		// XibNorm_systerr[i] = XibNorm[i]*xib_syst;
		// XibNorm_err[i]     = sqrt(pow(XibNorm_StatErr[i],2) + pow(XibNorm_systerr[i],2)); //Combining stat and syst in quadrature
		//
		// XibNorm_systerr_wt[i] = XibNorm_wt[i]*xib_syst;
		// XibNorm_err_wt[i]     = sqrt(pow(XibNorm_wt_StatErr[i],2) + pow(XibNorm_systerr_wt[i],2)); //Combining stat and syst in quadrature
		//
		// XibNorm[i]         = XibNorm[i]*2; //ACCOUNT FOR XIB0
		// XibNorm_err[i]     = XibNorm_err[i]*1.414; //ACCOUNT FOR XIB0
		//
		// XibNorm_wt[i]         = XibNorm_wt[i]*2; //ACCOUNT FOR XIB0
		// XibNorm_err_wt[i]     = XibNorm_err_wt[i]*1.414; //ACCOUNT FOR XIB0

		cout<<"************************************************"<<endl;
		cout<<"The UNWEIGHTED LL Xib normalization for Run "<<run
		    <<" is "<<XibNorm[i]<<" +/- "<<XibNorm_StatErr[i]
		    <<" +/- "<<XibNorm_SystErr[i]<<endl;

		cout<<"The WEIGHTED LL Xib normalization for Run "<<run
		    <<" is "<<XibNorm_wt[i]<<" +/- "<<XibNorm_wt_StatErr[i]
		    <<" +/- "<<XibNorm_wt_SystErr[i]<<endl;
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
		                   + mcTreeIn_Zero_Lambda->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));                                                       //NOTE NO TM HERE

		eff_Lambda_rec[i]     = num_Lambda*1.0/nGen_Lambda[i]; //Calc. reco eff.
		eff_Lambda_rec_err[i] = sqrt(eff_Lambda_rec[i]*(1-eff_Lambda_rec[i])/nGen_Lambda[i]); //statistical error on recon. eff.

		eff_Lambda_rec_wt[i]     = num_Lambda_wt*1.0/nGen_Lambda_wt[i]; //Calc. weighted reco eff.
		eff_Lambda_rec_err_wt[i] = sqrt(eff_Lambda_rec_wt[i]*(1-eff_Lambda_rec_wt[i])/nGen_Lambda_wt[i]); //statistical error on weighted recon. eff.

		cout<<"Run "<<run<<" UNWEIGHTED Lambda Recons. Effs = "<<eff_Lambda_rec[i]*100
		    <<" % +/- "<<eff_Lambda_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED Lambda Recons. Effs = "<<eff_Lambda_rec_wt[i]*100
		    <<" % +/- "<<eff_Lambda_rec_err_wt[i]*100<<" %"<<endl;

		eff_JpsiLambda[i]     = eff_Lambda_rec[i]*eff_Lambda_gen[i]; // Calc. total eff.
		eff_JpsiLambda_SystErr[i] = eff_JpsiLambda[i]*sqrt(pow((eff_Lambda_gen_err[i]/eff_Lambda_gen[i]),2) +
		                                                   pow((eff_Lambda_rec_err[i]/eff_Lambda_rec[i]),2));                                                                                                                                                                                                                                                                             // and stat error on tot. eff.

		eff_JpsiLambda_wt[i]     = eff_Lambda_rec_wt[i]*eff_Lambda_gen[i];                                                                                                                                                                          // Calc. total eff.
		eff_JpsiLambda_wt_SystErr[i] = eff_JpsiLambda_wt[i]*sqrt(pow((eff_Lambda_gen_err[i]/eff_Lambda_gen[i]),2) +
		                                                         pow((eff_Lambda_rec_err_wt[i]/eff_Lambda_rec_wt[i]),2));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 // and stat error on tot. eff.

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi Lambda Eff = "<<eff_JpsiLambda[i]*100
		    <<" % +/- "<<eff_JpsiLambda_SystErr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi Lambda Eff = "<<eff_JpsiLambda_wt[i]*100
		    <<" % +/- "<<eff_JpsiLambda_wt_SystErr[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;
	}
	//*******************************************************************

	RooHistPdf* SIG[2];
	RooKeysPdf* SIG_KEYS[2];
	RooDataSet* ds_sig[2];
	RooDataSet* ds_sig_wt[2];
	RooAbsReal* sigmaInt[2];

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

		genFile_Sigma>>nGen_Sigma[i];                                                                                                 // Get number of generated events

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

		genEffFile_Sigma>>eff_Sigma_gen[i];                                                                                                 // Get generator efficiency
		genEffFile_Sigma>>eff_Sigma_gen_err[i];                                                                                                 // and error on above

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

		eff_Sigma_rec[i]     = num_Sigma*1.0/nGen_Sigma[i];                                                                                                 // Calc. reco. eff.
		eff_Sigma_rec_err[i] = sqrt(eff_Sigma_rec[i]*(1-eff_Sigma_rec[i])/nGen_Sigma[i]);                                                                                                 //stat error on recon. eff.

		eff_Sigma_rec_wt[i]     = num_Sigma_wt*1.0/nGen_Sigma_wt[i];                                                                                                 //Calc. weighted reco eff.
		eff_Sigma_rec_err_wt[i] = sqrt(eff_Sigma_rec_wt[i]*(1-eff_Sigma_rec_wt[i])/nGen_Sigma_wt[i]);                                                                                                 //statistical error on weighted recon. eff.

		cout<<"Run "<<run<<" Sigma Recons. Effs = "
		    <<eff_Sigma_rec[i]*100<<" % +/- "<<eff_Sigma_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED Sigma Recons. Effs = "<<eff_Sigma_rec_wt[i]*100
		    <<" % +/- "<<eff_Sigma_rec_err_wt[i]*100<<" %"<<endl;

		eff_JpsiSigma[i] = eff_Sigma_gen[i] * eff_Sigma_rec[i];                                                                                                 // Calc overall eff.
		eff_JpsiSigma_SystErr[i] = eff_JpsiSigma[i]*sqrt(pow((eff_Sigma_gen_err[i]/eff_Sigma_gen[i]),2) +
		                                                 pow((eff_Sigma_rec_err[i]/eff_Sigma_rec[i]),2));                                                                                                                                                                                                                 // and stat. error on above

		eff_JpsiSigma_wt[i]     = eff_Sigma_rec_wt[i]*eff_Sigma_gen[i];                                                                                                                                                                                                                                                                          // Calc. total eff.
		eff_JpsiSigma_wt_SystErr[i] = eff_JpsiSigma_wt[i]*sqrt(pow((eff_Sigma_gen_err[i]/eff_Sigma_gen[i]),2) +
		                                                       pow((eff_Sigma_rec_err_wt[i]/eff_Sigma_rec_wt[i]),2));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     // and stat error on tot. eff.

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi Sigma Eff = "<<eff_JpsiSigma[i]*100
		    <<" % +/- "<<eff_JpsiSigma_SystErr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi Sigma Eff = "<<eff_JpsiSigma_wt[i]*100
		    <<" % +/- "<<eff_JpsiSigma_wt_SystErr[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;

		eff_ratio[i]         = eff_JpsiLambda[i]/eff_JpsiSigma[i];                                                                                                 // Calc eff ratio.
		eff_ratio_SystErr[i] = eff_ratio[i]*sqrt(pow((eff_JpsiSigma_SystErr[i]/eff_JpsiSigma[i]),2)+pow((eff_JpsiLambda_SystErr[i]/eff_JpsiLambda[i]),2));                                                                                                 // stat err on ratio
		eff_ratio_StatErr[i] = 0;

		eff_ratio_wt[i]         = eff_JpsiLambda_wt[i]/eff_JpsiSigma_wt[i];                                                                                                 // Calc eff ratio.
		eff_ratio_wt_SystErr[i] = eff_ratio_wt[i]*sqrt(pow((eff_JpsiSigma_wt_SystErr[i]/eff_JpsiSigma_wt[i]),2)+pow((eff_JpsiLambda_wt_SystErr[i]/eff_JpsiLambda_wt[i]),2));                                                                                                 // stat err on ratio
		eff_ratio_wt_StatErr[i] = 0;

		cout<<"***************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Lambda/Sigma efficiency ratio = "<<eff_ratio[i]<<" +/- "<<eff_ratio_SystErr[i]<<endl;
		cout<<"Run "<<run<<" WEIGHTED Lambda/Sigma efficiency ratio   = "<<eff_ratio_wt[i]<<" +/- "<<eff_ratio_wt_SystErr[i]<<endl;
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

		TTree* mcTreeIn_Zero_Sigma_cut    = (TTree*)mcTreeIn_Zero_Sigma->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));                                                                                                //Not TRUTH MATCHING HERE!
		TTree* mcTreeIn_nonZero_Sigma_cut = (TTree*)mcTreeIn_nonZero_Sigma->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));                                                                                                //Not TRUTH MATCHING HERE!

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

		// RooPlot *framesigma = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		// framesigma->SetTitle("J/#psi #Sigma");
		// // ds_sigma->plotOn(framesigma,Name("sigmadata"));
		// ds_sig[i]->plotOn(framesigma,Name("sigmadata_nowt"),LineColor(kGreen));
		// ds_sig_wt[i]->plotOn(framesigma,Name("sigmadata"),LineColor(kBlack));
		// // sigmashape.plotOn(framesigma,Name("sigmafit"),LineColor(kBlue));
		// (*(SIG_KEYS[i])).plotOn(framesigma,Name("sigmafitsmooth"),LineColor(kRed),LineStyle(kDashed));
		//
		// TCanvas *csigma = new TCanvas(Form("JpsiSigma%d",run),Form("JpsiSigma%d",run));
		// framesigma->Draw();

		w.import(*(SIG_KEYS[i]));

		cout<<"Done importing Jpsi Sigma shape"<<endl;

		sigmaInt[i] = SIG_KEYS[i]->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
		window_JpsiSigma[i] = sigmaInt[i]->getValV();

		// cout<<"Run "<<run<<" sigma fraction = "<<sigmaINT[i]<<endl;
	}
	//********************************************************************

	RooDataSet* ds_1405[2];
	RooKeysPdf* KEYS_1405[2];

	//****Get J/psi Lst(1405) efficiencies and shape from MC*******
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	cout<<"Get Lb -> J/psi Lst(1405) efficiency and shape"<<endl;
	cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<endl;
	const char* Lst1405Path = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/Lst1405";
	const char* rwSuffix    = "";
	const char* rwType      = "";
	if(Lst1405_rwtype==0)
	{
		nGen_1405[0] = 3902198;
		nGen_1405[1] = 3460130;
		rwSuffix     = "";
	}
	else if(Lst1405_rwtype==1)//MV RW
	{
		nGen_1405[0] = 3890603;
		nGen_1405[1] = 3450028;
		rwSuffix     = "_MV";
		rwType       = "MV";
	}
	else if(Lst1405_rwtype==2)//BONN RW
	{
		nGen_1405[0] = 3901943;
		nGen_1405[1] = 3460080;
		rwSuffix     = "_BONN";
		rwType       = "BONN";
	}

	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_1405 = Open(Form("%s/run%d/lst1405_cutoutks_LL_nonZeroTracks_noPID.root",
		                                         Lst1405Path,run));
		TTree *mcTreeIn_nonZero_1405 = (TTree*)mcFileIn_nonZero_1405->Get("MyTuple");

		TFile *mcFileIn_Zero_1405    = Open(Form("%s/run%d/lst1405_cutoutks_LL_ZeroTracks_noPID.root",
		                                         Lst1405Path,run));
		TTree *mcTreeIn_Zero_1405    = (TTree*)mcFileIn_Zero_1405->Get("MyTuple");

		mcTreeIn_nonZero_1405->AddFriend("MyTuple",Form("%s/run%d/lst1405_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                                Lst1405Path,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_1405->AddFriend("MyTuple",Form("%s/run%d/lst1405_zeroTracksLL_FinalBDT%d_noPID.root",
		                                             Lst1405Path,run,bdtConf_Zero[i]));

		TFile *genWtsFile_1405 = nullptr;
		TTree *genWtsTree_1405 = nullptr;
		if(run == 1)
		{
			genWtsFile_1405 = Open(Form("%s/run%d/RW/gbWeights_gen_new.root",
			                            Lst1405Path,run));
		}
		else if(run == 2)
		{
			genWtsFile_1405 = Open(Form("%s/run%d/RW/gbWeights_gen.root",
			                            Lst1405Path,run));
		}
		genWtsTree_1405 = (TTree*)genWtsFile_1405->Get("MyTuple");

		genWtsTree_1405->AddFriend("MyTuple",Form("%s/run%d/RW/tauWeights_gen.root",
		                                          Lst1405Path,run));

		if(run == 1)
		{
			genWtsTree_1405->Draw("gb_wts_new*wt_tau>>genWt_1405","","goff");
		}
		else if(run == 2)
		{
			genWtsTree_1405->Draw("gb_wts*wt_tau>>genWt_1405","","goff");
		}

		TH1F *genWt_1405 = (TH1F*)gDirectory->Get("genWt_1405");
		nGen_1405_wt[i]  = genWt_1405->GetEntries()*genWt_1405->GetMean();

		fstream genEffFile_1405;
		genEffFile_1405.open(Form("../logs/mc/JpsiLambda/Lst1405/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_1405>>eff_1405_gen[i]; // Get generator efficiency
		genEffFile_1405>>eff_1405_gen_err[i]; // and error on above

		cout<<"Run "<<run<<" Lst1405 Generator Effs = "<<eff_1405_gen[i]*100
		    <<" % +/- "<<eff_1405_gen_err[i]*100<<" %"<<endl;

		Int_t num_1405 = 0, num_1405_wt;
		if(Lst1405_rwtype==0)
		{
			num_1405 = mcTreeIn_nonZero_1405->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i])) +
			           mcTreeIn_Zero_1405->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));      //NOTE NO TM HERE

			if(run == 1)
			{
				mcTreeIn_nonZero_1405->Draw("gb_wts_new*wt_tau>>wt_1405_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
				mcTreeIn_Zero_1405->Draw("gb_wts_new*wt_tau>>wt_1405_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
			}
			else if(run == 2)
			{
				mcTreeIn_nonZero_1405->Draw("gb_wts*wt_tau>>wt_1405_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
				mcTreeIn_Zero_1405->Draw("gb_wts*wt_tau>>wt_1405_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
			}
			TH1F *wt_1405_nonZero = (TH1F*)gDirectory->Get("wt_1405_nonZero");
			TH1F *wt_1405_Zero    = (TH1F*)gDirectory->Get("wt_1405_Zero");

			num_1405_wt = (wt_1405_nonZero->GetMean()*wt_1405_nonZero->GetEntries()) +
			              (wt_1405_Zero->GetMean()*wt_1405_Zero->GetEntries());

		}
		else //HAVE TO DEAL WITH THIS
		{
			mcTreeIn_nonZero_1405->Draw(Form("%sweight>>hrw0",rwType),Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_1405->Draw(Form("%sweight>>hrw1",rwType),Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
			TH1D* hrw0 = (TH1D*)gDirectory->Get("hrw0");
			TH1D* hrw1 = (TH1D*)gDirectory->Get("hrw1");
			num_1405   = (hrw0->GetMean()*hrw0->GetEntries()) +(hrw1->GetMean()*hrw1->GetEntries());

			mcTreeIn_nonZero_1405->Draw(Form("%sweight*gb_wts*wt_tau>>wt_1405_nonZero",rwType),Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_1405->Draw(Form("%sweight*gb_wts*wt_tau>>wt_1405_Zero",rwType),Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");

			TH1F *wt_1405_nonZero = (TH1F*)gDirectory->Get("wt_1405_nonZero");
			TH1F *wt_1405_Zero    = (TH1F*)gDirectory->Get("wt_1405_Zero");

			num_1405_wt = (wt_1405_nonZero->GetMean()*wt_1405_nonZero->GetEntries()) +
			              (wt_1405_Zero->GetMean()*wt_1405_Zero->GetEntries());

		}

		eff_1405_rec[i]     = num_1405*1.0/nGen_1405[i]; // Calc. reco. eff.
		eff_1405_rec_err[i] = sqrt(eff_1405_rec[i]*(1-eff_1405_rec[i])/nGen_1405[i]); //stat error on recon. eff.

		eff_1405_wt_rec[i]     = num_1405_wt*1.0/nGen_1405_wt[i]; //Calc. weighted reco eff.
		eff_1405_wt_rec_err[i] = sqrt(eff_1405_wt_rec[i]*(1-eff_1405_wt_rec[i])/nGen_1405_wt[i]); //statistical error on weighted recon. eff.

		cout<<"Run "<<run<<" UNWEIGHTED /\(1405) Recons. Effs = "<<eff_1405_rec[i]*100
		    <<" % +/- "<<eff_1405_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED /\(1405) Recons. Effs = "<<eff_1405_wt_rec[i]*100
		    <<" % +/- "<<eff_1405_wt_rec_err[i]*100<<" %"<<endl;

		eff_1405[i]         = eff_1405_rec[i]*eff_1405_gen[i];
		eff_1405_SystErr[i] = eff_1405[i]*sqrt(pow((eff_1405_gen_err[i]/eff_1405_gen[i]),2) +
		                                       pow((eff_1405_rec_err[i]/eff_1405_rec[i]),2));

		eff_1405_wt[i]         = eff_1405_wt_rec[i]*eff_1405_gen[i];
		eff_1405_wt_SystErr[i] = eff_1405_wt[i]*sqrt(pow((eff_1405_gen_err[i]/eff_1405_gen[i]),2) +
		                                             pow((eff_1405_wt_rec_err[i]/eff_1405_wt_rec[i]),2));

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi /\(1405) Eff = "<<eff_1405[i]*100
		    <<" % +/- "<<eff_1405_SystErr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi /\(1405) Eff = "<<eff_1405_wt[i]*100
		    <<" % +/- "<<eff_1405_wt_SystErr[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;

		eff_ratio_1405[i]            = eff_1405[i]/eff_JpsiLambda[i];
		eff_ratio_SystErr_1405[i]    = eff_ratio_1405[i]*sqrt(pow((eff_1405_SystErr[i]/eff_1405[i]),2)+pow((eff_JpsiLambda_SystErr[i]/eff_JpsiLambda[i]),2)); // stat err on ratio
		eff_ratio_StatErr_1405[i]    = 0.0;
		eff_ratio_Err_1405[i]        = sqrt(pow(eff_ratio_StatErr_1405[i],2) + pow(eff_ratio_SystErr_1405[i],2));//combine in quadrature

		eff_ratio_wt_1405[i]            = eff_1405_wt[i]/eff_JpsiLambda_wt[i];
		eff_ratio_wt_SystErr_1405[i]    = eff_ratio_wt_1405[i]*sqrt(pow((eff_1405_wt_SystErr[i]/eff_1405_wt[i]),2)+pow((eff_JpsiLambda_wt_SystErr[i]/eff_JpsiLambda_wt[i]),2)); // stat err on ratio
		eff_ratio_wt_StatErr_1405[i]    = 0.0;
		eff_ratio_Err_1405_wt[i]        = sqrt(pow(eff_ratio_wt_StatErr_1405[i],2) + pow(eff_ratio_wt_SystErr_1405[i],2));//combine in quadrature

		cout<<"***************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED /\(1405)/Lambda efficiency ratio = "<<eff_ratio_1405[i]<<" +/- "<<eff_ratio_Err_1405[i] <<endl;
		cout<<"Run "<<run<<" WEIGHTED /\(1405)/Lambda efficiency ratio   = "<<eff_ratio_wt_1405[i]<<" +/- "<<eff_ratio_Err_1405_wt[i]<<endl;
		cout<<"***************************************"<<endl;

		//*******************************************************************

		//****************Shape*********************************************
		cout<<"***************************************"<<endl;
		cout<<"Get J/psi Lambda(1405) shape from MC"<<endl;
		cout<<"***************************************"<<endl;

		mcTreeIn_Zero_1405->SetBranchStatus("*",0);
		mcTreeIn_Zero_1405->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_Zero_1405->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
		mcTreeIn_Zero_1405->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_Zero_1405->SetBranchStatus("gb_wts",1);
		mcTreeIn_Zero_1405->SetBranchStatus("wt_tau",1);

		mcTreeIn_nonZero_1405->SetBranchStatus("*",0);
		mcTreeIn_nonZero_1405->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_nonZero_1405->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		mcTreeIn_nonZero_1405->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_nonZero_1405->SetBranchStatus("gb_wts",1);
		mcTreeIn_nonZero_1405->SetBranchStatus("wt_tau",1);

		if(run == 1)
		{
			mcTreeIn_Zero_1405->SetBranchStatus("gb_wts_new",1);
			mcTreeIn_nonZero_1405->SetBranchStatus("gb_wts_new",1);
		}

		TFile *tempFile_1405 = new TFile("tempFile_1405.root","RECREATE");

		TTree* mcTreeIn_Zero_1405_cut    = (TTree*)mcTreeIn_Zero_1405->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));//NOTE NO TRUTH MATCHING HERE!
		TTree* mcTreeIn_nonZero_1405_cut = (TTree*)mcTreeIn_nonZero_1405->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));//NOTE NO TRUTH MATCHING HERE!

		TList *list_1405 = new TList;
		list_1405->Add(mcTreeIn_Zero_1405_cut);
		list_1405->Add(mcTreeIn_nonZero_1405_cut);

		TTree *combTree_1405 = TTree::MergeTrees(list_1405);
		combTree_1405->SetName("combTree_1405");

		if(Lst1405_rwtype!=0)
		{
			if(run == 1)
			{
				ds_1405[i] = new RooDataSet("ds_1405","ds_1405",combTree_1405,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,Form("%sweight*gb_wts_new*wt_tau",rwType));
			}
			else if(run == 2)
			{
				ds_1405[i] = new RooDataSet("ds_1405","ds_1405",combTree_1405,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,Form("%sweight*gb_wts*wt_tau",rwType));
			}
		}
		else
		{
			if(run == 1)
			{
				ds_1405[i] = new RooDataSet("ds_1405","ds_1405",combTree_1405,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,"gb_wts_new*wt_tau");
			}
			else if(run == 2)
			{
				ds_1405[i] = new RooDataSet("ds_1405","ds_1405",combTree_1405,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,"gb_wts*wt_tau");
			}
		}
		ds_1405[i]->Print();

		KEYS_1405[i] = new RooKeysPdf(Form("LST1405_Run%d",run),Form("LST1405_Run%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*(ds_1405[i]),RooKeysPdf::MirrorBoth,1);

		// mcTreeIn_nonZero_1405->Draw(Form("Lb_DTF_M_JpsiLConstr>>h1405_nonZero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                            Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");        //Not TRUTH MATCHING HERE!
		//
		// mcTreeIn_Zero_1405->Draw(Form("Lb_DTF_M_JpsiLConstr>>h1405_Zero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                         Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");        //Not TRUTH MATCHING HERE!
		//
		// TH1D *h1405_nonZero = (TH1D*)gDirectory->Get(Form("h1405_nonZero%d",run));
		// TH1D *h1405_Zero    = (TH1D*)gDirectory->Get(Form("h1405_Zero%d",run));
		// TH1D *h1405         = new TH1D(Form("h1405%d",run),"",nbins,myLow,myHigh);
		//
		// h1405->Add(h1405_nonZero,h1405_Zero);
		// TH1D *h1405_smooth  = (TH1D*)h1405->Clone(Form("h1405_smooth%d",run));
		// h1405_smooth->Smooth(2);
		//
		// RooDataHist *ds_1405        = new RooDataHist(Form("ds_1405%d",run),Form("ds_1405%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),h1405);
		// RooDataHist *ds_1405_smooth = new RooDataHist(Form("ds_1405_smooth%d",run),Form("ds_1405_smooth%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),h1405_smooth);
		//
		// ds_1405->Print();
		//
		// RooHistPdf 1405shape(Form("1405shape%d",run),Form("1405shape%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_1405,0);
		// SIG[i] = new RooHistPdf(Form("SIG%d",run),Form("SIG%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_1405_smooth,0);

		// RooPlot *frame1405 = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		// frame1405->SetTitle("J/#psi #Lambda(1405)");
		// ds_1405[i]->plotOn(frame1405,Name("1405data"));
		// // 1405shape.plotOn(frame1405,Name("1405fit"),LineColor(kBlue));
		// (*(KEYS_1405[i])).plotOn(frame1405,Name("1405fitsmooth"),LineColor(kRed),LineStyle(kDashed));
		//
		// TCanvas *c1405 = new TCanvas(Form("JpsiLst(1405)%d",run),Form("JpsiLst(1405)%d",run));
		// frame1405->Draw();

		w.import(*(KEYS_1405[i]));

		cout<<"Done importing Jpsi Lst(1405) shape"<<endl;
	}
	//*******************************************************************

	RooDataSet* ds_1520[2];
	RooKeysPdf* KEYS_1520[2];

	//****Get J/psi Lst(1520) efficiencies and shape from MC*******
	const char* Lst1520Path = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/Lst1520";

	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_1520 = Open(Form("%s/run%d/lst1520_cutoutks_LL_nonZeroTracks_noPID.root",
		                                         Lst1520Path,run));
		TTree *mcTreeIn_nonZero_1520 = (TTree*)mcFileIn_nonZero_1520->Get("MyTuple");

		TFile *mcFileIn_Zero_1520    = Open(Form("%s/run%d/lst1520_cutoutks_LL_ZeroTracks_noPID.root",
		                                         Lst1520Path,run));
		TTree *mcTreeIn_Zero_1520    = (TTree*)mcFileIn_Zero_1520->Get("MyTuple");

		mcTreeIn_nonZero_1520->AddFriend("MyTuple",Form("%s/run%d/lst1520_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                                Lst1520Path,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_1520->AddFriend("MyTuple",Form("%s/run%d/lst1520_zeroTracksLL_FinalBDT%d_noPID.root",
		                                             Lst1520Path,run,bdtConf_Zero[i]));
		fstream genFile_1520;
		genFile_1520.open((Form("../logs/mc/JpsiLambda/Lst1520/run%d/gen_log.txt",run)));

		genFile_1520>>nGen_1520[i]; // Get number of generated events

		TFile *genWtsFile_1520 = Open(Form("%s/run%d/RW/gbWeights_gen.root",
		                                   Lst1520Path,run));
		TTree *genWtsTree_1520 = (TTree*)genWtsFile_1520->Get("MyTuple");

		genWtsTree_1520->AddFriend("MyTuple",Form("%s/run%d/RW/tauWeights_gen.root",
		                                          Lst1520Path,run));

		genWtsTree_1520->Draw("gb_wts*wt_tau>>genWt_1520","","goff");

		TH1F *genWt_1520 = (TH1F*)gDirectory->Get("genWt_1520");
		nGen_1520_wt[i] = genWt_1520->GetEntries()*genWt_1520->GetMean();

		fstream genEffFile_1520;
		genEffFile_1520.open(Form("../logs/mc/JpsiLambda/Lst1520/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_1520>>eff_1520_gen[i]; // Get generator efficiency
		genEffFile_1520>>eff_1520_gen_err[i]; // and error on above

		cout<<"Run "<<run<<" Lst1520 Generator Effs = "<<eff_1520_gen[i]*100
		    <<" % +/- "<<eff_1520_gen_err[i]*100<<" %"<<endl;

		Int_t num_1520 = mcTreeIn_nonZero_1520->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i])) +
		                 mcTreeIn_Zero_1520->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));    //NOTE NO TM HERE

		if(run == 1)
		{
			mcTreeIn_nonZero_1520->Draw("gb_wts_new*wt_tau>>wt_1520_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_1520->Draw("gb_wts_new*wt_tau>>wt_1520_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}
		else if(run == 2)
		{
			mcTreeIn_nonZero_1520->Draw("gb_wts*wt_tau>>wt_1520_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_1520->Draw("gb_wts*wt_tau>>wt_1520_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");

		}
		TH1F *wt_1520_nonZero = (TH1F*)gDirectory->Get("wt_1520_nonZero");
		TH1F *wt_1520_Zero    = (TH1F*)gDirectory->Get("wt_1520_Zero");

		Int_t num_1520_wt = (wt_1520_nonZero->GetMean()*wt_1520_nonZero->GetEntries()) +
		                    (wt_1520_Zero->GetMean()*wt_1520_Zero->GetEntries());

		eff_1520_rec[i]     = num_1520*1.0/nGen_1520[i]; // Calc. reco. eff.
		eff_1520_rec_err[i] = sqrt(eff_1520_rec[i]*(1-eff_1520_rec[i])/nGen_1520[i]); //syst error on recon. eff.

		eff_1520_rec_wt[i]     = num_1520_wt*1.0/nGen_1520_wt[i]; //Calc. weighted reco eff.
		eff_1520_wt_rec_err[i] = sqrt(eff_1520_rec_wt[i]*(1-eff_1520_rec_wt[i])/nGen_1520_wt[i]); //syst error on weighted recon. eff.

		cout<<"Run "<<run<<" UNWEIGHTED /\(1520) Recons. Effs = "<<eff_1520_rec[i]*100
		    <<" % +/- "<<eff_1520_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED /\(1520) Recons. Effs = "<<eff_1520_rec_wt[i]*100
		    <<" % +/- "<<eff_1520_wt_rec_err[i]*100<<" %"<<endl;

		eff_1520[i]         = eff_1520_rec[i]*eff_1520_gen[i];
		eff_1520_SystErr[i] = eff_1520[i]*sqrt(pow((eff_1520_gen_err[i]/eff_1520_gen[i]),2) +
		                                       pow((eff_1520_rec_err[i]/eff_1520_rec[i]),2));

		eff_1520_wt[i]         = eff_1520_rec_wt[i]*eff_1520_gen[i];
		eff_1520_wt_SystErr[i] = eff_1520_wt[i]*sqrt(pow((eff_1520_gen_err[i]/eff_1520_gen[i]),2) +
		                                             pow((eff_1520_wt_rec_err[i]/eff_1520_rec_wt[i]),2));

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi /\(1520) Eff = "<<eff_1520[i]*100
		    <<" % +/- "<<eff_1520_SystErr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi /\(1520) Eff = "<<eff_1520_wt[i]*100
		    <<" % +/- "<<eff_1520_wt_SystErr[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;

		eff_ratio_1520[i]            = eff_1520[i]/eff_JpsiLambda[i];
		eff_ratio_SystErr_1520[i]    = eff_ratio_1520[i]*sqrt(pow((eff_1520_SystErr[i]/eff_1520[i]),2)+pow((eff_JpsiLambda_SystErr[i]/eff_JpsiLambda[i]),2)); // stat err on ratio
		eff_ratio_StatErr_1520[i]    = 0.0;
		eff_ratio_Err_1520[i]        = sqrt(pow(eff_ratio_StatErr_1520[i],2) + pow(eff_ratio_SystErr_1520[i],2));//combine in quadrature

		eff_ratio_wt_1520[i]            = eff_1520_wt[i]/eff_JpsiLambda_wt[i];
		eff_ratio_wt_SystErr_1520[i]    = eff_ratio_wt_1520[i]*sqrt(pow((eff_1520_wt_SystErr[i]/eff_1520_wt[i]),2)+pow((eff_JpsiLambda_wt_SystErr[i]/eff_JpsiLambda_wt[i]),2)); // stat err on ratio
		eff_ratio_wt_StatErr_1520[i]    = 0.0;
		eff_ratio_Err_1520_wt[i]        = sqrt(pow(eff_ratio_wt_StatErr_1520[i],2) + pow(eff_ratio_wt_SystErr_1520[i],2));//combine in quadrature

		cout<<"***************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED /\(1520)/Lambda efficiency ratio = "<<eff_ratio_1520[i]<<" +/- "<<eff_ratio_Err_1520[i] <<endl;
		cout<<"Run "<<run<<" WEIGHTED /\(1520)/Lambda efficiency ratio   = "<<eff_ratio_wt_1520[i]<<" +/- "<<eff_ratio_Err_1520_wt[i]<<endl;
		cout<<"***************************************"<<endl;

		//*******************************************************************

		//****************Shape*********************************************
		cout<<"***************************************"<<endl;
		cout<<"Get J/psi Lambda(1520) shape from MC"<<endl;
		cout<<"***************************************"<<endl;

		mcTreeIn_Zero_1520->SetBranchStatus("*",0);
		mcTreeIn_Zero_1520->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_Zero_1520->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
		mcTreeIn_Zero_1520->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_Zero_1520->SetBranchStatus("gb_wts",1);
		mcTreeIn_Zero_1520->SetBranchStatus("wt_tau",1);

		mcTreeIn_nonZero_1520->SetBranchStatus("*",0);
		mcTreeIn_nonZero_1520->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_nonZero_1520->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		mcTreeIn_nonZero_1520->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_nonZero_1520->SetBranchStatus("gb_wts",1);
		mcTreeIn_nonZero_1520->SetBranchStatus("wt_tau",1);

		if(run == 1)
		{
			mcTreeIn_Zero_1520->SetBranchStatus("gb_wts_new",1);
			mcTreeIn_nonZero_1520->SetBranchStatus("gb_wts_new",1);
		}

		TFile *tempFile_1520 = new TFile("tempFile_1520.root","RECREATE");

		TTree* mcTreeIn_Zero_1520_cut    = (TTree*)mcTreeIn_Zero_1520->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));//NOTE NO TRUTH MATCHING HERE!
		TTree* mcTreeIn_nonZero_1520_cut = (TTree*)mcTreeIn_nonZero_1520->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));//NOTE NO TRUTH MATCHING HERE!

		TList *list_1520 = new TList;
		list_1520->Add(mcTreeIn_Zero_1520_cut);
		list_1520->Add(mcTreeIn_nonZero_1520_cut);

		TTree *combTree_1520 = TTree::MergeTrees(list_1520);
		combTree_1520->SetName("combTree_1520");

		ds_1520[i] = new RooDataSet("ds_1520","ds_1520",combTree_1520,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,"gb_wts*wt_tau");
		ds_1520[i]->Print();

		KEYS_1520[i] = new RooKeysPdf(Form("LST1520_Run%d",run),Form("LST1520_Run%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*(ds_1520[i]),RooKeysPdf::MirrorBoth,1);

		// mcTreeIn_nonZero_1520->Draw(Form("Lb_DTF_M_JpsiLConstr>>h1520_nonZero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                            Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");        //Not TRUTH MATCHING HERE!
		//
		// mcTreeIn_Zero_1520->Draw(Form("Lb_DTF_M_JpsiLConstr>>h1520_Zero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                         Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");        //Not TRUTH MATCHING HERE!
		//
		// TH1D *h1520_nonZero = (TH1D*)gDirectory->Get(Form("h1520_nonZero%d",run));
		// TH1D *h1520_Zero    = (TH1D*)gDirectory->Get(Form("h1520_Zero%d",run));
		// TH1D *h1520         = new TH1D(Form("h1520%d",run),"",nbins,myLow,myHigh);
		//
		// h1520->Add(h1520_nonZero,h1520_Zero);
		// TH1D *h1520_smooth  = (TH1D*)h1520->Clone(Form("h1520_smooth%d",run));
		// h1520_smooth->Smooth(2);
		//
		// RooDataHist *ds_1520        = new RooDataHist(Form("ds_1520%d",run),Form("ds_1520%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),h1520);
		// RooDataHist *ds_1520_smooth = new RooDataHist(Form("ds_1520_smooth%d",run),Form("ds_1520_smooth%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),h1520_smooth);
		//
		// ds_1520->Print();
		//
		// RooHistPdf 1520shape(Form("1520shape%d",run),Form("1520shape%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_1520,0);
		// SIG[i] = new RooHistPdf(Form("SIG%d",run),Form("SIG%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_1520_smooth,0);

		// RooPlot *frame1520 = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		// frame1520->SetTitle("J/#psi #Lambda(1520)");
		// ds_1520[i]->plotOn(frame1520,Name("1520data"));
		// (*(KEYS_1520[i])).plotOn(frame1520,Name("1520fitsmooth"),LineColor(kRed),LineStyle(kDashed));

		// TCanvas *c1520 = new TCanvas(Form("JpsiLst(1520)%d",run),Form("JpsiLst(1520)%d",run));
		// frame1520->Draw();

		w.import(*(KEYS_1520[i]));

		cout<<"Done importing Jpsi Lst(1520) shape"<<endl;
	}
	//*******************************************************************

	RooDataSet* ds_1600[2];
	RooKeysPdf* KEYS_1600[2];

	//****Get J/psi Lst(1600) efficiencies and shape from MC*******
	const char* Lst1600Path = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/Lst1600";

	for(Int_t run = 1; run<=2; run++)
	{
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_1600 = Open(Form("%s/run%d/lst1600_cutoutks_LL_nonZeroTracks_noPID.root",
		                                         Lst1600Path,run));
		TTree *mcTreeIn_nonZero_1600 = (TTree*)mcFileIn_nonZero_1600->Get("MyTuple");

		TFile *mcFileIn_Zero_1600    = Open(Form("%s/run%d/lst1600_cutoutks_LL_ZeroTracks_noPID.root",
		                                         Lst1600Path,run));
		TTree *mcTreeIn_Zero_1600    = (TTree*)mcFileIn_Zero_1600->Get("MyTuple");

		mcTreeIn_nonZero_1600->AddFriend("MyTuple",Form("%s/run%d/lst1600_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                                Lst1600Path,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_1600->AddFriend("MyTuple",Form("%s/run%d/lst1600_zeroTracksLL_FinalBDT%d_noPID.root",
		                                             Lst1600Path,run,bdtConf_Zero[i]));
		fstream genFile_1600;
		genFile_1600.open((Form("../logs/mc/JpsiLambda/Lst1600/run%d/gen_log.txt",run)));

		genFile_1600>>nGen_1600[i]; // Get number of generated events

		TFile *genWtsFile_1600 = Open(Form("%s/run%d/RW/gbWeights_gen.root",
		                                   Lst1600Path,run));
		TTree *genWtsTree_1600 = (TTree*)genWtsFile_1600->Get("MyTuple");

		genWtsTree_1600->AddFriend("MyTuple",Form("%s/run%d/RW/tauWeights_gen.root",
		                                          Lst1600Path,run));

		genWtsTree_1600->Draw("gb_wts*wt_tau>>genWt_1600","","goff");

		TH1F *genWt_1600 = (TH1F*)gDirectory->Get("genWt_1600");
		nGen_1600_wt[i] = genWt_1600->GetEntries()*genWt_1600->GetMean();

		fstream genEffFile_1600;
		genEffFile_1600.open(Form("../logs/mc/JpsiLambda/Lst1600/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_1600>>eff_1600_gen[i]; // Get generator efficiency
		genEffFile_1600>>eff_1600_gen_err[i]; // and error on above

		cout<<"Run "<<run<<" Lst1600 Generator Effs = "<<eff_1600_gen[i]*100
		    <<" % +/- "<<eff_1600_gen_err[i]*100<<" %"<<endl;

		Int_t num_1600 = mcTreeIn_nonZero_1600->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i])) +
		                 mcTreeIn_Zero_1600->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));    //NOTE NO TM HERE

		if(run == 1)
		{
			mcTreeIn_nonZero_1600->Draw("gb_wts_new*wt_tau>>wt_1600_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_1600->Draw("gb_wts_new*wt_tau>>wt_1600_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}
		else if(run == 2)
		{
			mcTreeIn_nonZero_1600->Draw("gb_wts*wt_tau>>wt_1600_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_1600->Draw("gb_wts*wt_tau>>wt_1600_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");

		}
		TH1F *wt_1600_nonZero = (TH1F*)gDirectory->Get("wt_1600_nonZero");
		TH1F *wt_1600_Zero    = (TH1F*)gDirectory->Get("wt_1600_Zero");

		Int_t num_1600_wt = (wt_1600_nonZero->GetMean()*wt_1600_nonZero->GetEntries()) +
		                    (wt_1600_Zero->GetMean()*wt_1600_Zero->GetEntries());

		eff_1600_rec[i]     = num_1600*1.0/nGen_1600[i]; // Calc. reco. eff.
		eff_1600_rec_err[i] = sqrt(eff_1600_rec[i]*(1-eff_1600_rec[i])/nGen_1600[i]); //syst error on recon. eff.

		eff_1600_rec_wt[i]     = num_1600_wt*1.0/nGen_1600_wt[i]; //Calc. weighted reco eff.
		eff_1600_wt_rec_err[i] = sqrt(eff_1600_rec_wt[i]*(1-eff_1600_rec_wt[i])/nGen_1600_wt[i]); //syst error on weighted recon. eff.

		cout<<"Run "<<run<<" UNWEIGHTED /\(1600) Recons. Effs = "<<eff_1600_rec[i]*100
		    <<" % +/- "<<eff_1600_rec_err[i]*100<<" %"<<endl;

		cout<<"Run "<<run<<" WEIGHTED /\(1600) Recons. Effs = "<<eff_1600_rec_wt[i]*100
		    <<" % +/- "<<eff_1600_wt_rec_err[i]*100<<" %"<<endl;

		eff_1600[i]         = eff_1600_rec[i]*eff_1600_gen[i];
		eff_1600_SystErr[i] = eff_1600[i]*sqrt(pow((eff_1600_gen_err[i]/eff_1600_gen[i]),2) +
		                                       pow((eff_1600_rec_err[i]/eff_1600_rec[i]),2));

		eff_1600_wt[i]         = eff_1600_rec_wt[i]*eff_1600_gen[i];
		eff_1600_wt_SystErr[i] = eff_1600_wt[i]*sqrt(pow((eff_1600_gen_err[i]/eff_1600_gen[i]),2) +
		                                             pow((eff_1600_wt_rec_err[i]/eff_1600_rec_wt[i]),2));

		cout<<"************************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED Jpsi /\(1600) Eff = "<<eff_1600[i]*100
		    <<" % +/- "<<eff_1600_SystErr[i]*100<<" %"<<endl;
		cout<<"Run "<<run<<" WEIGHTED Jpsi /\(1600) Eff = "<<eff_1600_wt[i]*100
		    <<" % +/- "<<eff_1600_wt_SystErr[i]*100<<" %"<<endl;
		cout<<"************************************************"<<endl;

		eff_ratio_1600[i]            = eff_1600[i]/eff_JpsiLambda[i];
		eff_ratio_SystErr_1600[i]    = eff_ratio_1600[i]*sqrt(pow((eff_1600_SystErr[i]/eff_1600[i]),2)+pow((eff_JpsiLambda_SystErr[i]/eff_JpsiLambda[i]),2)); // stat err on ratio
		eff_ratio_StatErr_1600[i]    = 0.0;
		eff_ratio_Err_1600[i]        = sqrt(pow(eff_ratio_StatErr_1600[i],2) + pow(eff_ratio_SystErr_1600[i],2));//combine in quadrature

		eff_ratio_wt_1600[i]            = eff_1600_wt[i]/eff_JpsiLambda_wt[i];
		eff_ratio_wt_SystErr_1600[i]    = eff_ratio_wt_1600[i]*sqrt(pow((eff_1600_wt_SystErr[i]/eff_1600_wt[i]),2)+pow((eff_JpsiLambda_wt_SystErr[i]/eff_JpsiLambda_wt[i]),2)); // stat err on ratio
		eff_ratio_wt_StatErr_1600[i]    = 0.0;
		eff_ratio_Err_1600_wt[i]        = sqrt(pow(eff_ratio_wt_StatErr_1600[i],2) + pow(eff_ratio_wt_SystErr_1600[i],2));//combine in quadrature

		cout<<"***************************************"<<endl;
		cout<<"Run "<<run<<" UNWEIGHTED /\(1600)/Lambda efficiency ratio = "<<eff_ratio_1600[i]<<" +/- "<<eff_ratio_Err_1600[i] <<endl;
		cout<<"Run "<<run<<" WEIGHTED /\(1600)/Lambda efficiency ratio   = "<<eff_ratio_wt_1600[i]<<" +/- "<<eff_ratio_Err_1600_wt[i]<<endl;
		cout<<"***************************************"<<endl;

		//*******************************************************************

		//****************Shape*********************************************
		cout<<"***************************************"<<endl;
		cout<<"Get J/psi Lambda(1600) shape from MC"<<endl;
		cout<<"***************************************"<<endl;

		mcTreeIn_Zero_1600->SetBranchStatus("*",0);
		mcTreeIn_Zero_1600->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_Zero_1600->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
		mcTreeIn_Zero_1600->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_Zero_1600->SetBranchStatus("gb_wts",1);
		mcTreeIn_Zero_1600->SetBranchStatus("wt_tau",1);

		mcTreeIn_nonZero_1600->SetBranchStatus("*",0);
		mcTreeIn_nonZero_1600->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_nonZero_1600->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		mcTreeIn_nonZero_1600->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_nonZero_1600->SetBranchStatus("gb_wts",1);
		mcTreeIn_nonZero_1600->SetBranchStatus("wt_tau",1);

		if(run == 1)
		{
			mcTreeIn_Zero_1600->SetBranchStatus("gb_wts_new",1);
			mcTreeIn_nonZero_1600->SetBranchStatus("gb_wts_new",1);
		}

		TFile *tempFile_1600 = new TFile("tempFile_1600.root","RECREATE");

		TTree* mcTreeIn_Zero_1600_cut    = (TTree*)mcTreeIn_Zero_1600->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));//NOTE NO TRUTH MATCHING HERE!
		TTree* mcTreeIn_nonZero_1600_cut = (TTree*)mcTreeIn_nonZero_1600->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));//NOTE NO TRUTH MATCHING HERE!

		TList *list_1600 = new TList;
		list_1600->Add(mcTreeIn_Zero_1600_cut);
		list_1600->Add(mcTreeIn_nonZero_1600_cut);

		TTree *combTree_1600 = TTree::MergeTrees(list_1600);
		combTree_1600->SetName("combTree_1600");

		ds_1600[i] = new RooDataSet("ds_1600","ds_1600",combTree_1600,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,"gb_wts*wt_tau");
		ds_1600[i]->Print();

		KEYS_1600[i] = new RooKeysPdf(Form("LST1600_Run%d",run),Form("LST1600_Run%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*(ds_1600[i]),RooKeysPdf::MirrorBoth,1);

		// mcTreeIn_nonZero_1600->Draw(Form("Lb_DTF_M_JpsiLConstr>>h1600_nonZero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                            Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");        //Not TRUTH MATCHING HERE!
		//
		// mcTreeIn_Zero_1600->Draw(Form("Lb_DTF_M_JpsiLConstr>>h1600_Zero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                         Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");        //Not TRUTH MATCHING HERE!
		//
		// TH1D *h1600_nonZero = (TH1D*)gDirectory->Get(Form("h1600_nonZero%d",run));
		// TH1D *h1600_Zero    = (TH1D*)gDirectory->Get(Form("h1600_Zero%d",run));
		// TH1D *h1600         = new TH1D(Form("h1600%d",run),"",nbins,myLow,myHigh);
		//
		// h1600->Add(h1600_nonZero,h1600_Zero);
		// TH1D *h1600_smooth  = (TH1D*)h1600->Clone(Form("h1600_smooth%d",run));
		// h1600_smooth->Smooth(2);
		//
		// RooDataHist *ds_1600        = new RooDataHist(Form("ds_1600%d",run),Form("ds_1600%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),h1600);
		// RooDataHist *ds_1600_smooth = new RooDataHist(Form("ds_1600_smooth%d",run),Form("ds_1600_smooth%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),h1600_smooth);
		//
		// ds_1600->Print();
		//
		// RooHistPdf 1600shape(Form("1600shape%d",run),Form("1600shape%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_1600,0);
		// SIG[i] = new RooHistPdf(Form("SIG%d",run),Form("SIG%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_1600_smooth,0);

		// RooPlot *frame1600 = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		// frame1600->SetTitle("J/#psi #Lambda(1600)");
		// ds_1600[i]->plotOn(frame1600,Name("1600data"));
		// (*(KEYS_1600[i])).plotOn(frame1600,Name("1600fitsmooth"),LineColor(kRed),LineStyle(kDashed));

		// TCanvas *c1600 = new TCanvas(Form("JpsiLst(1600)%d",run),Form("JpsiLst(1600)%d",run));
		// frame1600->Draw();

		w.import(*(KEYS_1600[i]));

		cout<<"Done importing Jpsi Lst(1600) shape"<<endl;
	}
	//*******************************************************************
	//
	// RooDataSet* ds_chic1[2];
	// RooKeysPdf* KEYS_chic1[2];
	//
	// //****Get chiC1 Lambda efficiencies and shape from MC*******
	// const char* chiC1Path = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/chiC1";
	//
	// for(Int_t run = 1; run<=2; run++)
	// {
	//      Int_t i = run-1;
	//      TFile *mcFileIn_nonZero_chic1 = Open(Form("%s/run%d/chic1_cutoutks_LL_nonZeroTracks_noPID.root",
	//                                                chiC1Path,run));
	//      TTree *mcTreeIn_nonZero_chic1 = (TTree*)mcFileIn_nonZero_chic1->Get("MyTuple");
	//
	//      TFile *mcFileIn_Zero_chic1    = Open(Form("%s/run%d/chic1_cutoutks_LL_ZeroTracks_noPID.root",
	//                                                chiC1Path,run));
	//      TTree *mcTreeIn_Zero_chic1    = (TTree*)mcFileIn_Zero_chic1->Get("MyTuple");
	//
	//      mcTreeIn_nonZero_chic1->AddFriend("MyTuple",Form("%s/run%d/chic1_LL_FinalBDT%d_iso%d_%s_noPID.root",
	//                                                       chiC1Path,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
	//      mcTreeIn_Zero_chic1->AddFriend("MyTuple",Form("%s/run%d/chic1_zeroTracksLL_FinalBDT%d_noPID.root",
	//                                                    chiC1Path,run,bdtConf_Zero[i]));
	//      fstream genFile_chic1;
	//      genFile_chic1.open((Form("../logs/mc/JpsiLambda/chiC1/run%d/gen_log.txt",run)));
	//
	//      genFile_chic1>>nGen_chic1[i]; // Get number of generated events
	//
	//      TFile *genWtsFile_chic1 = Open(Form("%s/run%d/RW/gbWeights_gen.root",
	//                                          chiC1Path,run));
	//      TTree *genWtsTree_chic1 = (TTree*)genWtsFile_chic1->Get("MyTuple");
	//
	//      genWtsTree_chic1->AddFriend("MyTuple",Form("%s/run%d/RW/tauWeights_gen.root",
	//                                                 chiC1Path,run));
	//
	//      genWtsTree_chic1->Draw("gb_wts*wt_tau>>genWt_chic1","","goff");
	//
	//      TH1F *genWt_chic1 = (TH1F*)gDirectory->Get("genWt_chic1");
	//      nGen_chic1_wt[i] = genWt_chic1->GetEntries()*genWt_chic1->GetMean();
	//
	//      //***** DONT HAVE GENERATOR EFFS FOR chic1 YET
	//      // fstream genEffFile_chic1;
	//      // genEffFile_chic1.open(Form("../logs/mc/JpsiLambda/chiC1/run%d/Generator_Effs_Combined.txt",run));
	//      //
	//      // genEffFile_chic1>>eff_chic1_gen[i]; // Get generator efficiency
	//      // genEffFile_chic1>>eff_chic1_gen_err[i]; // and error on above
	//      //
	//      // cout<<"Run "<<run<<" chic1 Generator Effs = "<<eff_chic1_gen[i]*100
	//      //     <<" % +/- "<<eff_chic1_gen_err[i]*100<<" %"<<endl;
	//
	//      //FOR NOW, EFF RATIO IS JUST RATIO OF RECO EFFS. ONCE I GET GENERATOR EFFS, I'LL PUT THEM IN
	//
	//      Int_t num_chic1 = mcTreeIn_nonZero_chic1->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i])) +
	//                        mcTreeIn_Zero_chic1->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));                                               //NOTE NO TM HERE
	//
	//      mcTreeIn_nonZero_chic1->Draw("gb_wts*wt_tau>>wt_chic1_nonZero",Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]));
	//      mcTreeIn_Zero_chic1->Draw("gb_wts*wt_tau>>wt_chic1_Zero",Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));
	//
	//      TH1F *wt_chic1_nonZero = (TH1F*)gDirectory->Get("wt_chic1_nonZero");
	//      TH1F *wt_chic1_Zero    = (TH1F*)gDirectory->Get("wt_chic1_Zero");
	//
	//      Int_t num_chic1_wt = (wt_chic1_nonZero->GetMean()*wt_chic1_nonZero->GetEntries()) +
	//                           (wt_chic1_Zero->GetMean()*wt_chic1_Zero->GetEntries());
	//
	//      eff_chic1_rec[i]     = num_chic1*1.0/nGen_chic1[i]; // Calc. reco. eff.
	//      eff_chic1_rec_err[i] = sqrt(eff_chic1_rec[i]*(1-eff_chic1_rec[i])/nGen_chic1[i]); //stat error on recon. eff.
	//
	//      eff_chic1_rec_wt[i]     = num_chic1_wt*1.0/nGen_chic1_wt[i]; //Calc. weighted reco eff.
	//      eff_chic1_rec_err_wt[i] = sqrt(eff_chic1_rec_wt[i]*(1-eff_chic1_rec_wt[i])/nGen_chic1_wt[i]); //statistical error on weighted recon. eff.
	//
	//      cout<<"Run "<<run<<" UNWEIGHTED chic1 Recons. Effs = "<<eff_chic1_rec[i]*100
	//          <<" % +/- "<<eff_chic1_rec_err[i]*100<<" %"<<endl;
	//
	//      cout<<"Run "<<run<<" WEIGHTED chic1 Recons. Effs = "<<eff_chic1_rec_wt[i]*100
	//          <<" % +/- "<<eff_chic1_rec_err_wt[i]*100<<" %"<<endl;
	//
	//      eff_chic1[i] = eff_chic1_rec[i]; //TEMPORARY
	//      eff_chic1_staterr[i] = eff_chic1_rec_err[i]; //TEMPORARY
	//
	//      eff_chic1_wt[i]         = eff_chic1_rec_wt[i]; //TEMPORARY
	//      eff_chic1_staterr_wt[i] = eff_chic1_rec_err_wt[i]; //TEMPORARY
	//
	//      // eff_chic1[i] = eff_chic1_gen[i] * eff_chic1_rec[i]; // Calc overall eff.
	//      // eff_chic1_staterr[i] = eff_chic1[i]*sqrt(pow((eff_chic1_gen_err[i]/eff_chic1_gen[i]),2) +
	//      //                                        pow((eff_chic1_rec_err[i]/eff_chic1_rec[i]),2));   // and stat. error on above
	//
	//      cout<<"************************************************"<<endl;
	//      cout<<"Run "<<run<<" UNWEIGHTED Jpsi chic1 Eff = "<<eff_chic1[i]*100
	//          <<" % +/- "<<eff_chic1_staterr[i]*100<<" %"<<endl;
	//      cout<<"Run "<<run<<" WEIGHTED Jpsi chic1 Eff = "<<eff_chic1_wt[i]*100
	//          <<" % +/- "<<eff_chic1_staterr_wt[i]*100<<" %"<<endl;
	//      cout<<"************************************************"<<endl;
	//
	//      // eff_ratio[i]         = eff_chic1[i]/eff_JpsiLambda[i]; // Calc eff ratio.
	//
	//      eff_ratio_chic1[i] = eff_chic1[i]/eff_Lambda_rec[i]; //TEMPORARY - FOR NOW EFF RATIO IS JUST RATIO OF REC EFFS.
	//      // eff_ratio_StatErr_chic1[i] = eff_ratio_chic1[i]*sqrt(pow((eff_chic1_staterr[i]/eff_chic1[i]),2)+pow((eff_JpsiLambda_SystErr[i]/eff_JpsiLambda[i]),2)); // stat err on ratio
	//      eff_ratio_StatErr_chic1[i] = eff_ratio_chic1[i]*sqrt(pow((eff_chic1_staterr[i]/eff_chic1[i]),2)+pow((eff_Lambda_rec_err[i]/eff_Lambda_rec[i]),2)); // TEMPORARY
	//      eff_ratio_SystErr_chic1[i] = eff_ratio_chic1[i]*eff_ratio_syst_chic1;
	//      eff_ratio_Err_chic1[i] = sqrt(pow(eff_ratio_StatErr_chic1[i],2) + pow(eff_ratio_SystErr_chic1[i],2));//combine in quadrature
	//
	//      eff_ratio_chic1_wt[i]            = eff_chic1_wt[i]/eff_Lambda_rec_wt[i]; //TEMPORARY - FOR NOW EFF RATIO IS JUST RATIO OF REC EFFS.
	//      // eff_ratio_StatErr_chic1_wt[i] = eff_ratio_chic1_wt[i]*sqrt(pow((eff_chic1_staterr_wt[i]/eff_chic1_wt[i]),2)+pow((eff_JpsiLambda_wt_SystErr[i]/eff_JpsiLambda_wt[i]),2)); // stat err on ratio
	//      eff_ratio_StatErr_chic1_wt[i]    = eff_ratio_chic1_wt[i]*sqrt(pow((eff_chic1_staterr_wt[i]/eff_chic1_wt[i]),2)+pow((eff_Lambda_rec_err_wt[i]/eff_Lambda_rec_wt[i]),2)); // TEMPORARY
	//      eff_ratio_SystErr_chic1_wt[i]    = eff_ratio_chic1_wt[i]*eff_ratio_syst_chic1;
	//      eff_ratio_Err_chic1_wt[i]        = sqrt(pow(eff_ratio_StatErr_chic1_wt[i],2) + pow(eff_ratio_SystErr_chic1_wt[i],2));//combine in quadrature
	//
	//      cout<<"***************************************"<<endl;
	//      cout<<"Run "<<run<<" UNWEIGHTED chic1/Lambda efficiency ratio = "<<eff_ratio[i]<<" +/- "<<eff_ratio_Err[i]<<endl;
	//      cout<<"Run "<<run<<" WEIGHTED chic1/Lambda efficiency ratio   = "<<eff_ratio_wt[i]<<" +/- "<<eff_ratio_Err_wt[i]<<endl;
	//      cout<<"***************************************"<<endl;
	//      //*******************************************************************
	//
	//      //****************Shape*********************************************
	//      cout<<"***************************************"<<endl;
	//      cout<<"Get J/psi Lambda(chic1) shape from MC"<<endl;
	//      cout<<"***************************************"<<endl;
	//
	//      mcTreeIn_Zero_chic1->SetBranchStatus("*",0);
	//      mcTreeIn_Zero_chic1->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	//      mcTreeIn_Zero_chic1->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
	//      mcTreeIn_Zero_chic1->SetBranchStatus("Lb_BKGCAT",1);
	//      mcTreeIn_Zero_chic1->SetBranchStatus("gb_wts",1);
	//      mcTreeIn_Zero_chic1->SetBranchStatus("wt_tau",1);
	//
	//      mcTreeIn_nonZero_chic1->SetBranchStatus("*",0);
	//      mcTreeIn_nonZero_chic1->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	//      mcTreeIn_nonZero_chic1->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
	//      mcTreeIn_nonZero_chic1->SetBranchStatus("Lb_BKGCAT",1);
	//      mcTreeIn_nonZero_chic1->SetBranchStatus("gb_wts",1);
	//      mcTreeIn_nonZero_chic1->SetBranchStatus("wt_tau",1);
	//
	//      TFile *tempFile_chic1 = new TFile("tempFile_chic1.root","RECREATE");
	//
	//      TTree* mcTreeIn_Zero_chic1_cut    = (TTree*)mcTreeIn_Zero_chic1->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));//NOTE NO TRUTH MATCHING HERE!
	//      TTree* mcTreeIn_nonZero_chic1_cut = (TTree*)mcTreeIn_nonZero_chic1->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));//NOTE NO TRUTH MATCHING HERE!
	//
	//      TList *list_chic1 = new TList;
	//      list_chic1->Add(mcTreeIn_Zero_chic1_cut);
	//      list_chic1->Add(mcTreeIn_nonZero_chic1_cut);
	//
	//      TTree *combTree_chic1 = TTree::MergeTrees(list_chic1);
	//      combTree_chic1->SetName("combTree_chic1");
	//
	//      ds_chic1[i] = new RooDataSet("ds_chic1","ds_chic1",combTree_chic1,RooArgSet(*(w.var("Lb_DTF_M_JpsiLConstr"))),0,"gb_wts*wt_tau");
	//      ds_chic1[i]->Print();
	//
	//      KEYS_chic1[i] = new RooKeysPdf(Form("chic1_Run%d",run),Form("chic1_Run%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*(ds_chic1[i]),RooKeysPdf::MirrorBoth,1);
	//
	//      // mcTreeIn_nonZero_chic1->Draw(Form("Lb_DTF_M_JpsiLConstr>>hchic1_nonZero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
	//      //                            Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");        //Not TRUTH MATCHING HERE!
	//      //
	//      // mcTreeIn_Zero_chic1->Draw(Form("Lb_DTF_M_JpsiLConstr>>hchic1_Zero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
	//      //                         Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");        //Not TRUTH MATCHING HERE!
	//      //
	//      // TH1D *hchic1_nonZero = (TH1D*)gDirectory->Get(Form("hchic1_nonZero%d",run));
	//      // TH1D *hchic1_Zero    = (TH1D*)gDirectory->Get(Form("hchic1_Zero%d",run));
	//      // TH1D *hchic1         = new TH1D(Form("hchic1%d",run),"",nbins,myLow,myHigh);
	//      //
	//      // hchic1->Add(hchic1_nonZero,hchic1_Zero);
	//      // TH1D *hchic1_smooth  = (TH1D*)hchic1->Clone(Form("hchic1_smooth%d",run));
	//      // hchic1_smooth->Smooth(2);
	//      //
	//      // RooDataHist *ds_chic1        = new RooDataHist(Form("ds_chic1%d",run),Form("ds_chic1%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),hchic1);
	//      // RooDataHist *ds_chic1_smooth = new RooDataHist(Form("ds_chic1_smooth%d",run),Form("ds_chic1_smooth%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),hchic1_smooth);
	//      //
	//      // ds_chic1->Print();
	//      //
	//      // RooHistPdf chic1shape(Form("chic1shape%d",run),Form("chic1shape%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_chic1,0);
	//      // SIG[i] = new RooHistPdf(Form("SIG%d",run),Form("SIG%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_chic1_smooth,0);
	//
	//      RooPlot *framechic1 = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
	//      framechic1->SetTitle("#chi_{C1} #Lambda");
	//      // ds_chic1->plotOn(framechic1,Name("chic1data"));
	//      ds_chic1[i]->plotOn(framechic1,Name("chiC1data"));
	//      // chic1shape.plotOn(framechic1,Name("chic1fit"),LineColor(kBlue));
	//      (*(KEYS_chic1[i])).plotOn(framechic1,Name("chiC1fitsmooth"),LineColor(kRed),LineStyle(kDashed));
	//
	//      TCanvas *cchic1 = new TCanvas(Form("chic1_Run%d",run),Form("chic1_Run%d",run));
	//      framechic1->Draw();
	//      w.import(*(KEYS_chic1[i]));
	//
	//      cout<<"Done importing chiC1 Lambda shape"<<endl;
	// }
	//*******************************************************************

	RooDataSet* ds_xi[2];
	RooDataSet* ds_xi_wt[2];
	RooKeysPdf* XIB_KEYS[2];

	Double_t xibCentral[2], xibErr[2], xibLow[2], xibHigh[2];
	xibLow[0] = 0;
	xibLow[1] = 0;

	xibHigh[0] = 200;
	xibHigh[1] = 400;

	//******************Get shape from Xib background********************
	// RooRealVar xibmass("Lb_DTF_M_JpsiLConstr","xibmass",5200.,5740.);
	const char* xibPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/JpsiXi";
	RooAbsReal* xibInt[2];

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
		// treein_xi_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_nonZero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                         Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");//TRUTH MATCHING HERE
		// treein_xi_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_Zero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
		//                      Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");//TRUTH MATCHING HERE
		//
		// TH1D *hxib_nonZero = (TH1D*)gDirectory->Get(Form("hxib_nonZero%d",run));
		// TH1D *hxib_Zero    = (TH1D*)gDirectory->Get(Form("hxib_Zero%d",run));
		// TH1D *hxib         = new TH1D(Form("hxib%d",run),"",nbins,myLow,myHigh);
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

		// RooPlot *framexib = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		// framexib->SetTitle("J/#psi #Xi");
		// ds_xi[i]->plotOn(framexib,Name("xibdata_nowt"),LineColor(kGreen));
		// ds_xi_wt[i]->plotOn(framexib,Name("xibdata"),LineColor(kBlack));
		// // (*(XIB[i])).plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));
		// (*(XIB_KEYS[i])).plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));
		//
		// TCanvas *cxib = new TCanvas(Form("JpsiXi%d",run),Form("JpsiXi%d",run));
		// framexib->Draw();

		// w.import(*(XIB[i]));
		w.import(*(XIB_KEYS[i]));

		cout<<"Done importing Xib shape"<<endl;

		xibInt[i] = XIB_KEYS[i]->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
		window_JpsiXi[i] = xibInt[i]->getValV();

		xibCentral[i] = XibNorm_wt[i];
		xibErr[i]     = sqrt(pow(XibNorm_wt_StatErr[i],2)+pow(XibNorm_wt_SystErr[i],2));

		cout<<"************************************************"<<endl;
		cout<<"The UNWEIGHTED Xib normalization inside signal region for Run "<<run
		    <<" is "<<XibNorm[i]<<" +/- "<<XibNorm_StatErr[i]
		    <<" +/- "<<XibNorm_SystErr[i]<<endl;
		cout<<"The WEIGHTED Xib normalization inside signal region for Run "<<run
		    <<" is "<<XibNorm_wt[i]<<" +/- "<<XibNorm_wt_StatErr[i]
		    <<" +/- "<<XibNorm_wt_SystErr[i]<<endl;
		cout<<"************************************************"<<endl;

	}
	//*******************************************************************

	//*********Double Crystal Ball signal shape for Lambda_b0************
	// w.factory("RooCBShape::Lb1_Run1(Lb_DTF_M_JpsiLConstr,mean_Run1[5619.6,5619,5621],"
	//           "sigma_Run1[10.,0.,20.], alpha1_Run1[1.028, 0.8,1.3], 10.0)" );
	// w.factory("RooCBShape::Lb2_Run1(Lb_DTF_M_JpsiLConstr,mean_Run1,sigma_Run1,alpha2_Run1[-1.097,-1.4,-0.7], 10.0)");
	// w.factory("SUM::Lb_Run1(0.5*Lb1_Run1 , 0.5*Lb2_Run1)");
	//
	// w.factory("RooCBShape::Lb1_Run2(Lb_DTF_M_JpsiLConstr,mean_Run2[5619.6,5619,5621],"
	//           "sigma_Run2[10.,0.,20.], alpha1_Run2[1.028, 0.8,1.3], 10.0)" );
	// w.factory("RooCBShape::Lb2_Run2(Lb_DTF_M_JpsiLConstr,mean_Run2,sigma_Run2,alpha2_Run2[-1.097,-1.4,-0.7], 10.0)");
	// w.factory("SUM::Lb_Run2(0.5*Lb1_Run2 , 0.5*Lb2_Run2)");
	//*******************************************************************

	//*********Hypatia signal shape for Lambda_b0************************

	// w.factory("RooHypatia2::Lb_Run1(Lb_DTF_M_JpsiLConstr,lambda_Run1[-2.0,-4.0,0.0],0,0,"
	//           "sigma_Run1[10.,1.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[1.7,1.0,3.0],"
	//           "2 ,a2_Run1[3.0,1.0,4.0], 2)");
	//
	// w.factory("RooHypatia2::Lb_Run2(Lb_DTF_M_JpsiLConstr,lambda_Run2[-2.5,-4.0,0.0],0,0,"
	//           "sigma_Run2[10.,1.,20.], mean_Run2[5619.6,5619,5621], a1_Run2[1.5,1.0,3.0],"
	//           "2 ,a2_Run2[1.5,1.0,3.0], 2)");

	w.factory("RooHypatia2::Lb_Run1(Lb_DTF_M_JpsiLConstr,lambda_Run1[-2.0,-4.0,0.0],0,0,"
	          "sigma_Run1[10.,1.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[15.0],"
	          "2 ,a2_Run1[15.0], 2)");

	w.factory("RooHypatia2::Lb_Run2(Lb_DTF_M_JpsiLConstr,lambda_Run2[-2.5,-4.0,0.0],0,0,"
	          "sigma_Run2[10.,1.,20.], mean_Run2[5619.6,5619,5621], a1_Run2[15.0],"
	          "2 ,a2_Run2[15.0], 2)");

	cout<<"Done defining J/psi Lambda Hypatia shapes"<<endl;
	//*******************************************************************

	//*********Continuum backgruond****************

	if(bkgType == 0)
	{
		cout<<"*****UING EXPONENTIAL BKG SHAPE*****"<<endl;
		w.factory("Exponential::Bkg_Run1(Lb_DTF_M_JpsiLConstr,expr('-1*pow(10,slope_Run1)',slope_Run1[-6.0,-11.0,-2.0]))");
		w.factory("Exponential::Bkg_Run2(Lb_DTF_M_JpsiLConstr,expr('-1*pow(10,slope_Run2)',slope_Run2[-6.0,-11.0,-2.0]))");

		w.var("slope_Run1")->setError(0.25);
		w.var("slope_Run2")->setError(0.25);
	}
	else if(bkgType == 1)
	{
		cout<<"*****UING 2nd ORDER CHEBYCHEV BKG SHAPE*****"<<endl;
		w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0], c1_Run1[0.0,-1.0,1.0]})");
		w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-2.0,2.0], c1_Run2[0.0,-1.0,1.0]})");
	}
	else if(bkgType == 2)
	{
		cout<<"*****UING 3rd ORDER CHEBYCHEV BKG SHAPE*****"<<endl;
		w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0], c1_Run1[0.0,-1.0,1.0], c2_Run1[0.0,-1.0,1.0]})");
		w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-2.0,2.0], c1_Run2[0.0,-1.0,1.0], c2_Run2[0.0,-1.0,1.0]})");
	}
	else if(bkgType == 3)
	{
		cout<<"*****UING 4th ORDER CHEBYCHEV BKG SHAPE*****"<<endl;
		w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0], c1_Run1[0.0,-1.0,1.0], c2_Run1[0.0,-1.0,1.0], c3_Run1[0.0,-1.0,1.0]})");
		w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-2.0,2.0], c1_Run2[0.0,-1.0,1.0], c2_Run2[0.0,-1.0,1.0], c3_Run2[0.0,-1.0,1.0]})");
	}

	cout<<"Done defining cont. bkg exp. shapes"<<endl;

	// w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[-0.5,-2.0,2.0], c1_Run1[0.5,-1.0,1.0], c2_Run1[0.,-1.0,1.0]})");
	// w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[-0.5,-2.0,2.0], c1_Run2[-0.5,-1.0,1.0], c2_Run2[0.,-1.0,1.0]})");

	//*******************************************************************

	//*********Gaussian Lump for misc. Lambda*'s ************************
	w.factory("Gaussian::lstLump_Run1(Lb_DTF_M_JpsiLConstr,miscLstMean_Run1[5009.,4990.,5040.],"
	          "miscLstSigma_Run1[60.,30.,90.])");

	w.factory("Gaussian::lstLump_Run2(Lb_DTF_M_JpsiLConstr,miscLstMean_Run2[5009.,4990.,5040.],"
	          "miscLstSigma_Run2[60.,30.,90.])");
	cout<<"Done defining Misc Lambda* Gaussian shapes"<<endl;
	//*******************************************************************

	RooDataHist *ds[2];
	RooDataSet *ds_unb[2];
	Int_t nentries[2];
	TH1D *myhist[2];

	//********Data for Fit to simulation*********************************
	TH1D *mcHist[2];

//	mcHist[0] = new TH1D("","",nbins,myLow,myHigh);
//	mcHist[1] = new TH1D("","",nbins,myLow,myHigh);

	RooDataHist *mc_ds[2];
	RooDataSet *ds_sim[2], *ds_sim_wt[2];
	Int_t mcNentries[2];

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

		mcTreeIn_Zero_Lambda->SetBranchStatus("*",0);
		mcTreeIn_Zero_Lambda->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_Zero_Lambda->SetBranchStatus(Form("BDT%d",bdtConf_Zero[i]),1);
		mcTreeIn_Zero_Lambda->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_Zero_Lambda->SetBranchStatus("gb_wts",1);
		mcTreeIn_Zero_Lambda->SetBranchStatus("wt_tau",1);
		if(run == 1)
			mcTreeIn_Zero_Lambda->SetBranchStatus("gb_wts_new",1);

		mcTreeIn_nonZero_Lambda->SetBranchStatus("*",0);
		mcTreeIn_nonZero_Lambda->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		mcTreeIn_nonZero_Lambda->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		mcTreeIn_nonZero_Lambda->SetBranchStatus("Lb_BKGCAT",1);
		mcTreeIn_nonZero_Lambda->SetBranchStatus("gb_wts",1);
		mcTreeIn_nonZero_Lambda->SetBranchStatus("wt_tau",1);
		if(run == 1)
			mcTreeIn_nonZero_Lambda->SetBranchStatus("gb_wts_new",1);

		TFile *tempFile = new TFile("tempFile_sim.root","RECREATE");

		TTree* mcTreeIn_Zero_Lambda_cut    = (TTree*)mcTreeIn_Zero_Lambda->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));                                                                                        //Not TRUTH MATCHING HERE!
		TTree* mcTreeIn_nonZero_Lambda_cut = (TTree*)mcTreeIn_nonZero_Lambda->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));                                                                                        //Not TRUTH MATCHING HERE!

		TList *list_sim = new TList;
		list_sim->Add(mcTreeIn_Zero_Lambda_cut);
		list_sim->Add(mcTreeIn_nonZero_Lambda_cut);

		TTree *combTree_sim = TTree::MergeTrees(list_sim);
		combTree_sim->SetName("combTree_sim");

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

		ds_sim[i] = new RooDataSet(Form("ds_sim%d",run),Form("ds_sim%d",run),combTree_sim,RooArgSet(*myVar,*gbWtVar,*tauWtVar));
		ds_sim[i]->Print();

		RooFormulaVar *totWt = new RooFormulaVar("totWt","@0*@1",RooArgList(*gbWtVar,*tauWtVar));
		RooRealVar *totWt_var = (RooRealVar*)ds_sim[i]->addColumn(*totWt);

		ds_sim_wt[i] = new RooDataSet(Form("ds_sim_wt%d",run),Form("ds_sim_wt%d",run),RooArgSet(*myVar,*totWt_var),Import(*(ds_sim[i])),WeightVar(*totWt_var));
		ds_sim_wt[i]->Print();

		if(run == 1)
		{
			mcTreeIn_nonZero_Lambda->Draw(Form("Lb_DTF_M_JpsiLConstr>>wt_Lambda_nonZero(%d,%d,%d)",nbins*2,myLow,myHigh),Form("(BDT%d > %f)*gb_wts_new*wt_tau", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_Lambda->Draw(Form("Lb_DTF_M_JpsiLConstr>>wt_Lambda_Zero(%d,%d,%d)",nbins*2,myLow,myHigh),Form("(BDT%d > %f)*gb_wts_new*wt_tau", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}
		else if(run == 2)
		{
			mcTreeIn_nonZero_Lambda->Draw(Form("Lb_DTF_M_JpsiLConstr>>wt_Lambda_nonZero(%d,%d,%d)",nbins*2,myLow,myHigh),Form("(BDT%d > %f)*gb_wts*wt_tau", bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
			mcTreeIn_Zero_Lambda->Draw(Form("Lb_DTF_M_JpsiLConstr>>wt_Lambda_Zero(%d,%d,%d)",nbins*2,myLow,myHigh),Form("(BDT%d > %f)*gb_wts*wt_tau", bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		}

		mcHist[i] = new TH1D(Form("mcHist%d",run),"",nbins*2,myLow,myHigh);
		TH1D *wt_Lambda_nonZero = (TH1D*)gDirectory->Get("wt_Lambda_nonZero");
		TH1D *wt_Lambda_Zero = (TH1D*)gDirectory->Get("wt_Lambda_Zero");

		(mcHist[i])->Sumw2();
		(mcHist[i])->Add(wt_Lambda_nonZero,wt_Lambda_Zero);

		mcNentries[i] = mcHist[i]->Integral();
		cout<<"mcNentries = "<<mcNentries[i]<<endl;

		// mc_ds[i] = new RooDataHist(Form("mc_ds%d",run),Form("mc_ds%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),mcHist[i]);
		// //  RooDataSet ds("ds","ds",treein,Lb_DTF_M_JpsiLConstr);
		// cout<<"Done making MC RooDataHist"<<endl;
		// (mc_ds[i])->Print();

		// w.import(*(mc_ds[i]));
		w.import(*(ds_sim_wt[i]));
	}

	w.factory("Exponential::mcbkg_run1(Lb_DTF_M_JpsiLConstr,tau_run1[-0.0007,-0.01,-0.0000001])");
	w.factory("Exponential::mcbkg_run2(Lb_DTF_M_JpsiLConstr,tau_run2[-0.0007,-0.01,-0.0000001])");

	w.factory(Form("mcnsig_run1[1,%d]",mcNentries[0]));
	w.factory(Form("mcnsig_run2[1,%d]",mcNentries[1]));

	w.factory(Form("mcnbkg_run1[1,%d]",mcNentries[0]));
	w.factory(Form("mcnbkg_run2[1,%d]",mcNentries[1]));

	w.factory("SUM:mcFit_Run1(mcnsig_run1*Lb_Run1, mcnbkg_run1*mcbkg_run1)");
	w.factory("SUM:mcFit_Run2(mcnsig_run2*Lb_Run2, mcnbkg_run2*mcbkg_run2)");

// (w.pdf("mcFit_Run1"))->fitTo(*(mc_ds[0]),Range(5500,5740));
// (w.pdf("mcFit_Run2"))->fitTo(*(mc_ds[1]),Range(5500,5740));

	(w.pdf("mcFit_Run1"))->fitTo(*(ds_sim_wt[0]),Range(5500,5740),Strategy(2));
	(w.pdf("mcFit_Run2"))->fitTo(*(ds_sim_wt[1]),Range(5500,5740),Strategy(2));

	TCanvas *sim_Run1      = new TCanvas("sim_Run1","sim_Run1",1200,800);
	RooPlot *simframe_Run1 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),5500,5740,120);

//	(mc_ds[0])->plotOn(simframe_Run1,Name("simdata_Run1"),DataError(RooAbsData::Poisson));
	(ds_sim_wt[0])->plotOn(simframe_Run1,Name("simdata_Run1"),DataError(RooAbsData::Poisson));
	(w.pdf("mcFit_Run1"))->plotOn(simframe_Run1,Name("simfit_Run1"));
	simframe_Run1->Draw();

	TCanvas *sim_Run2      = new TCanvas("sim_Run2","sim_Run2",1200,800);
	RooPlot *simframe_Run2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),5500,5740,120);

//	(mc_ds[1])->plotOn(simframe_Run2,Name("simdata_Run2"),DataError(RooAbsData::Poisson));
	(ds_sim_wt[1])->plotOn(simframe_Run2,Name("simdata_Run2"),DataError(RooAbsData::Poisson));
	(w.pdf("mcFit_Run2"))->plotOn(simframe_Run2,Name("simfit_Run2"));
	simframe_Run2->Draw();

	auto mcvars1 = w.pdf("Lb_Run1")->getVariables();
// auto mcp1 = (RooRealVar*)mcvars1->find("a1_Run1");
// mcp1->setConstant(kTRUE);
// auto mcp2 = (RooRealVar*)mcvars1->find("a2_Run1");
// mcp2->setConstant(kTRUE);
	auto mcp3 = (RooRealVar*)mcvars1->find("lambda_Run1");
	mcp3->setConstant(kTRUE);
	delete mcvars1;

	auto mcvars2 = w.pdf("Lb_Run2")->getVariables();
// auto mcp4 = (RooRealVar*)mcvars2->find("a1_Run2");
// mcp4->setConstant(kTRUE);
// auto mcp5 = (RooRealVar*)mcvars2->find("a2_Run2");
// mcp5->setConstant(kTRUE);
	auto mcp6 = (RooRealVar*)mcvars2->find("lambda_Run2");
	mcp6->setConstant(kTRUE);
	delete mcvars2;
	//*******************************************************************

	//*********Input Data************************************************
	const char *dataPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda";

	for(Int_t run = 1; run<=2; run++) {
		cout<<"Importing Run "<<run<<" data"<<endl;
		Int_t i = run-1;

		TFile *filein_nonZero = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root",
		                                  dataPath,run),"READ");
		TTree *treein_nonZero = (TTree*)filein_nonZero->Get("MyTuple");

		TFile *filein_Zero = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_ZeroTracks_noPID.root",
		                               dataPath,run),"READ");
		TTree *treein_Zero = (TTree*)filein_Zero->Get("MyTuple");

		treein_nonZero->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_LL_FinalBDT%d_iso%d_%s_noPID.root",
		                                         dataPath,run,bdtConf_nonZero[i],
		                                         isoConf[i],isoVersion[i]));
		treein_Zero->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_zeroTracksLL_FinalBDT%d_noPID.root",
		                                      dataPath,run,bdtConf_Zero[i]));

		if(!isBinned)
		{
			TFile *tempFile_data = new TFile("tempFile_data.root","RECREATE");

			TTree* treein_Zero_cut    = (TTree*)treein_Zero->CopyTree(Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));                                                                                                                                                //Not TRUTH MATCHING HERE!
			TTree* treein_nonZero_cut = (TTree*)treein_nonZero->CopyTree(Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));                                                                                                                                                //Not TRUTH MATCHING HERE!

			TList *list_data = new TList;
			list_data->Add(treein_Zero_cut);
			list_data->Add(treein_nonZero_cut);

			TTree *combTree_data = TTree::MergeTrees(list_data);
			combTree_data->SetName("combTree_data");
			nentries[i] = combTree_data->GetEntries();
			cout<<"Run "<<run<<" nentries = "<<nentries[i]<<endl;

			ds_unb[i] = new RooDataSet(Form("ds%d",run),Form("ds%d",run),combTree_data,*(w.var("Lb_DTF_M_JpsiLConstr")));
			(ds_unb[i])->Print();
			w.import(*(ds_unb[i]));
		}
		else
		{
			treein_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_nonzero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
			                     Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");

			treein_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_zero%d(%d,%d,%d)",run,nbins,myLow,myHigh),
			                  Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");

			TH1D *myhist_nonzero = (TH1D*)gDirectory->Get(Form("myhist_nonzero%d",run));
			TH1D *myhist_zero = (TH1D*)gDirectory->Get(Form("myhist_zero%d",run));

			myhist[i] = new TH1D(Form("myhist%d",run),"",nbins,myLow,myHigh);

			myhist[i]->Sumw2();
			myhist[i]->Add(myhist_zero,myhist_nonzero);

			nentries[i] = myhist[i]->Integral();
			cout<<"Run "<<run<<" nentries = "<<nentries[i]<<endl;

			ds[i] = new RooDataHist(Form("ds%d",run),Form("ds%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),myhist[i]);

			cout<<"Done making RooDataHist"<<endl;
			(ds[i])->Print();

			w.import(*(ds[i]));
		}

		cout<<"Done importing Run "<<run<<" data"<<endl;
	}
	//*******************************************************************

	//************************MAKE COMBINED MODEL************************

	w.factory("R[0,-5000,5000]"); // R*10^5 is the parameter of interest  This is shared b/w Run1 and Run2
	// NB: R is allowed to fluctuate negative in the fit. I don't completely understand why.
	// R*10^5 =  [N(Jpsi Sigma)/N(Jpsi Lambda)] * [eff(Jpsi Lambda)/eff(Jpsi Sigma)]
	// w.var("R")->setError(0.1);

	w.factory("R_1405[0.1,0.0,1.0]"); //R_1405 = N_corr(Lb -> Jpsi Lambda(1405)) / N_corr(Lb -> Jpsi Lambda)
	w.factory("R_1520[0.2,0.0,1.0]"); //R_1520 = N_corr(Lb -> Jpsi Lambda(1520)) / N_corr(Lb -> Jpsi Lambda)
	w.factory("R_1600[0.2,0.0,1.0]"); //R_1600 = N_corr(Lb -> Jpsi Lambda(1600)) / N_corr(Lb -> Jpsi Lambda)
	w.factory("R_chic1[0.1,0.0,0.5]"); //R_chic1 = N_corr(Lb -> chic1 Lambda / N_corr(Lb -> Jpsi Lambda)

	// w.factory(Form("sigmaEff1[%f,0,1]",eff_Sigma[0]));
	// w.factory(Form("lambdaEff1[%f,0,1]",eff_Lambda[0]));
	// w.factory(Form("sigmaEff2[%f,0,1]",eff_Sigma[1]));
	// w.factory(Form("lambdaEff2[%f,0,1]",eff_Lambda[1]));
	//
	// w.factory(Form("Gaussian::sigmaEff_constraint1(gsigmaEff1[%f,0,1],sigmaEff1,%f)",eff_Sigma[0],eff_Sigma_err[0]));
	// w.factory(Form("Gaussian::lambdaEff_constraint1(glambdaEff1[%f,0,1],lambdaEff1,%f)",eff_Lambda[0],eff_Lambda_err[0]));
	//
	// w.factory(Form("Gaussian::sigmaEff_constraint2(gsigmaEff2[%f,0,1],sigmaEff2,%f)",eff_Sigma[1],eff_Sigma_err[1]));
	// w.factory(Form("Gaussian::lambdaEff_constraint2(glambdaEff2[%f,0,1],lambdaEff2,%f)",eff_Lambda[1],eff_Lambda_err[1]));
	//
	// w.var("gsigmaEff1")->setConstant();
	// w.var("glambdaEff1")->setConstant();
	//
	// w.var("gsigmaEff2")->setConstant();
	// w.var("glambdaEff2")->setConstant();


	//**************** Sigma/Lb Efficiency Ratio*************************
	w.factory(Form("eff_ratio1[%f,0.,5.]",eff_ratio_wt[0])); //eff(Lb -> Jpsi Lambda)/eff(Lb -> Jpsi Sigma)
	w.factory(Form("eff_ratio2[%f,0.,5.]",eff_ratio_wt[1]));

	w.factory(Form("Gaussian::eff_ratio_constraint1(geff_ratio1[%f,0,5],eff_ratio1,%f)",eff_ratio_wt[0],eff_ratio_wt_SystErr[0]));
	w.factory(Form("Gaussian::eff_ratio_constraint2(geff_ratio2[%f,0,5],eff_ratio2,%f)",eff_ratio_wt[1],eff_ratio_wt_SystErr[1]));

	w.var("geff_ratio1")->setConstant();
	w.var("geff_ratio2")->setConstant();
	//*******************************************************************

	//**************** 1405/Lb Efficiency Ratio**************************
	w.factory(Form("eff_ratio_1405_1[%f,0.,5.]",eff_ratio_wt_1405[0])); //eff(Lb -> Jpsi Lambda(1405))/eff(Lb -> Jpsi Lambda)
	w.factory(Form("eff_ratio_1405_2[%f,0.,5.]",eff_ratio_wt_1405[1]));

	w.factory(Form("Gaussian::eff_ratio_1405_constraint1(geff_ratio_1405_1[%f,0,5],eff_ratio_1405_1,%f)",eff_ratio_wt_1405[0],eff_ratio_Err_1405_wt[0]));
	w.factory(Form("Gaussian::eff_ratio_1405_constraint2(geff_ratio_1405_2[%f,0,5],eff_ratio_1405_2,%f)",eff_ratio_wt_1405[1],eff_ratio_Err_1405_wt[1]));

	w.var("geff_ratio_1405_1")->setConstant();
	w.var("geff_ratio_1405_2")->setConstant();
	//*******************************************************************

	//**************** 1520/Lb Efficiency Ratio**************************
	w.factory(Form("eff_ratio_1520_1[%f,0.,5.]",eff_ratio_wt_1520[0])); //eff(Lb -> Jpsi Lambda(1520))/eff(Lb -> Jpsi Lambda)
	w.factory(Form("eff_ratio_1520_2[%f,0.,5.]",eff_ratio_wt_1520[1]));

	w.factory(Form("Gaussian::eff_ratio_1520_constraint1(geff_ratio_1520_1[%f,0,5],eff_ratio_1520_1,%f)",eff_ratio_wt_1520[0],eff_ratio_Err_1520_wt[0]));
	w.factory(Form("Gaussian::eff_ratio_1520_constraint2(geff_ratio_1520_2[%f,0,5],eff_ratio_1520_2,%f)",eff_ratio_wt_1520[1],eff_ratio_Err_1520_wt[1]));

	w.var("geff_ratio_1520_1")->setConstant();
	w.var("geff_ratio_1520_2")->setConstant();
	//*******************************************************************

	//**************** 1600/Lb Efficiency Ratio**************************
	w.factory(Form("eff_ratio_1600_1[%f,0.,5.]",eff_ratio_wt_1600[0])); //eff(Lb -> Jpsi Lambda(1600))/eff(Lb -> Jpsi Lambda)
	w.factory(Form("eff_ratio_1600_2[%f,0.,5.]",eff_ratio_wt_1600[1]));

	w.factory(Form("Gaussian::eff_ratio_1600_constraint1(geff_ratio_1600_1[%f,0,5],eff_ratio_1600_1,%f)",eff_ratio_wt_1600[0],eff_ratio_Err_1600_wt[0]));
	w.factory(Form("Gaussian::eff_ratio_1600_constraint2(geff_ratio_1600_2[%f,0,5],eff_ratio_1600_2,%f)",eff_ratio_wt_1600[1],eff_ratio_Err_1600_wt[1]));

	w.var("geff_ratio_1600_1")->setConstant();
	w.var("geff_ratio_1600_2")->setConstant();
	//*******************************************************************

	//**************** chic1/Lb Efficiency Ratio**************************
	// w.factory(Form("eff_ratio_chic1_1[%f,0.,1.]",eff_ratio_chic1_wt[0])); //eff(Lb -> Jpsi Lambda(chic1))/eff(Lb -> Jpsi Lambda)
	// w.factory(Form("eff_ratio_chic1_2[%f,0.,1.]",eff_ratio_chic1_wt[1]));
	//
	// w.factory(Form("Gaussian::eff_ratio_chic1_constraint1(geff_ratio_chic1_1[%f,0,1],eff_ratio_chic1_1,%f)",eff_ratio_chic1_wt[0],eff_ratio_Err_chic1_wt[0]));
	// w.factory(Form("Gaussian::eff_ratio_chic1_constraint2(geff_ratio_chic1_2[%f,0,1],eff_ratio_chic1_2,%f)",eff_ratio_chic1_wt[1],eff_ratio_Err_chic1_wt[1]));
	//
	// w.var("geff_ratio_chic1_1")->setConstant();
	// w.var("geff_ratio_chic1_2")->setConstant();
	//*******************************************************************

	w.factory("nLb_Run1[5000,1000,8000]");
	w.factory("nLb_Run2[17000,1000,24000]");

	// w.factory("n1405_Run1[500,1,5000]");
	// w.factory("n1405_Run2[500,1,5000]");
	// w.factory("n1520_Run1[3000,1,5000]");
	// w.factory("n1520_Run2[3000,1,5000]");
	// w.factory("n1600_Run1[3000,1,10000]");
	// w.factory("n1600_Run2[3000,1,10000]");
	w.factory("nMiscLst_Run1[1500,1,4000]");
	w.factory("nMiscLst_Run2[4000,1,9000]");
	w.factory(Form("nBkg_Run1[1,%d]",nentries[0]));
	w.factory(Form("nBkg_Run2[1,%d]",nentries[1]));

	//****************Xib Bkg Yield**************************************

	//What should the limits on nXib be?

	w.factory(Form("nXib1[%f,%f,%f]",xibCentral[0],xibLow[0],xibHigh[0]));
	w.factory(Form("nXib2[%f,%f,%f]",xibCentral[1],xibLow[1],xibHigh[1]));

	w.factory(Form("Gaussian::nXib_constraint1(gnXib1[%f,%f,%f],nXib1,%f)",xibCentral[0],xibLow[0],xibHigh[0],xibErr[0]));
	w.factory(Form("Gaussian::nXib_constraint2(gnXib2[%f,%f,%f],nXib2,%f)",xibCentral[1],xibLow[1],xibHigh[1],xibErr[1]));

	w.var("gnXib1")->setConstant();
	w.var("gnXib2")->setConstant();
	//*******************************************************************
	// w.factory("expr::nSigma1('pow(10,-5)*R*nLb_Run1*sigmaEff1/(lambdaEff1)',R,nLb_Run1,sigmaEff1,lambdaEff1)");
	// w.factory("expr::nSigma2('pow(10,-5)*R*nLb_Run2*sigmaEff2/(lambdaEff2)',R,nLb_Run2,sigmaEff2,lambdaEff2)");

	//*****************Jpsi Sigma yield**********************************
	w.factory("expr::nSigma1('pow(10,-5)*R*nLb_Run1/eff_ratio1',R,nLb_Run1,eff_ratio1)");
	w.factory("expr::nSigma2('pow(10,-5)*R*nLb_Run2/eff_ratio2',R,nLb_Run2,eff_ratio2)");
	//*******************************************************************

	//*****************Jpsi Lambda(1405) yield***************************
	w.factory("expr::n1405_Run1('R_1405*nLb_Run1*eff_ratio_1405_1',R_1405,nLb_Run1,eff_ratio_1405_1)");
	w.factory("expr::n1405_Run2('R_1405*nLb_Run2*eff_ratio_1405_2',R_1405,nLb_Run2,eff_ratio_1405_2)");
	//*******************************************************************

	//*****************Jpsi Lambda(1520) yield***************************
	w.factory("expr::n1520_Run1('R_1520*nLb_Run1*eff_ratio_1520_1',R_1520,nLb_Run1,eff_ratio_1520_1)");
	w.factory("expr::n1520_Run2('R_1520*nLb_Run2*eff_ratio_1520_2',R_1520,nLb_Run2,eff_ratio_1520_2)");
	//*******************************************************************

	//*****************Jpsi Lambda(1600) yield***************************
	w.factory("expr::n1600_Run1('R_1600*nLb_Run1*eff_ratio_1600_1',R_1600,nLb_Run1,eff_ratio_1600_1)");
	w.factory("expr::n1600_Run2('R_1600*nLb_Run2*eff_ratio_1600_2',R_1600,nLb_Run2,eff_ratio_1600_2)");
	//*******************************************************************

	//*****************chic1 Lambda yield***************************
	// w.factory("expr::nchic1_Run1('R_chic1*nLb_Run1*eff_ratio_chic1_1',R_chic1,nLb_Run1,eff_ratio_chic1_1)");
	// w.factory("expr::nchic1_Run2('R_chic1*nLb_Run2*eff_ratio_chic1_2',R_chic1,nLb_Run2,eff_ratio_chic1_2)");
	//*******************************************************************
	w.factory("SUM:model1(nSigma1*SIG1 , nLb_Run1*Lb_Run1 , nXib1*XIB1 , n1405_Run1*LST1405_Run1 ,"
	          " n1520_Run1*LST1520_Run1 , n1600_Run1*LST1600_Run1 , nMiscLst_Run1*lstLump_Run1 , nBkg_Run1*Bkg_Run1)");
	w.factory("SUM:model2(nSigma2*SIG2 , nLb_Run2*Lb_Run2 , nXib2*XIB2 , n1405_Run2*LST1405_Run2 ,"
	          " n1520_Run2*LST1520_Run2 , n1600_Run2*LST1600_Run2 , nMiscLst_Run2*lstLump_Run2 , nBkg_Run2*Bkg_Run2)");

	if(!lst1405flag)
	{
		w.var("R_1405")->setVal(0);
		w.var("R_1405")->setConstant();
	}
	if(!lst1520flag)
	{
		w.var("R_1520")->setVal(0);
		w.var("R_1520")->setConstant();
	}
	if(!lst1600flag)
	{
		w.var("R_1600")->setVal(0);
		w.var("R_1600")->setConstant();
	}
	// if(!chic1flag)
	// {
	//      w.var("R_chic1")->setVal(0);
	//      w.var("R_chic1")->setConstant();
	// }

	// w.factory("PROD::model_const1(model1,sigmaEff_constraint1,lambdaEff_constraint1,nXib_constraint1)");//Multiply model by constraint terms
	// w.factory("PROD::model_const2(model2,sigmaEff_constraint2,lambdaEff_constraint2,nXib_constraint2)");//Multiply model by constraint terms

	w.defineSet("poi","R"); //parameters of interest

	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,lambda_Run1,a1_Run1,a2_Run1,nBkg_Run1,slope_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,lambda_Run2,a1_Run2,a2_Run2,nBkg_Run2,slope_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff1,lambdaEff1,nXib1,n1405_Run2,n1520_Run2");

	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,lambda_Run1,a1_Run1,"
	//             "a2_Run1,nBkg_Run1,slope_Run1,nMiscLst_Run1,miscLstMean_Run1,"
	//             "miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,"
	//             "n1520_Run1");
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,lambda_Run2,a1_Run2,"
	//             "a2_Run2,nBkg_Run2,slope_Run2,nMiscLst_Run2,miscLstMean_Run2,"
	//             "miscLstSigma_Run2,sigmaEff2,lambdaEff2,nXib2,n1405_Run2,"
	//             "n1520_Run2");// define set of nuisance parameters


	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,"
	//             "nBkg_Run1,slope_Run1,nMiscLst_Run1,miscLstMean_Run1,"
	//             "miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,"
	//             "n1520_Run1");
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,"
	//             "nBkg_Run2,slope_Run2,nMiscLst_Run2,miscLstMean_Run2,"
	//             "miscLstSigma_Run2,sigmaEff2,lambdaEff2,nXib2,n1405_Run2,"
	//             "n1520_Run2");// define set of nuisance parameters

	w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,"
	            "lambda_Run1,a1_Run1,a2_Run1,"
	            "nBkg_Run1,nMiscLst_Run1,miscLstMean_Run1,"
	            "miscLstSigma_Run1,eff_ratio1,nXib1,eff_ratio_1405_1,"
	            "eff_ratio_1520_1,eff_ratio_1600_1");

	w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,"
	            "lambda_Run2,a1_Run2,a2_Run2,"
	            "nBkg_Run2,nMiscLst_Run2,miscLstMean_Run2,"
	            "miscLstSigma_Run2,eff_ratio2,nXib2,eff_ratio_1405_2,"
	            "eff_ratio_1520_2,eff_ratio_1600_2");

	if(bkgType == 0)
	{
		w.extendSet("nuisParams","slope_Run1,slope_Run2");
	}
	else if(bkgType == 1)
	{
		w.extendSet("nuisParams","c0_Run1,c0_Run2,c1_Run1,c1_Run2");
	}
	else if(bkgType == 2)
	{
		w.extendSet("nuisParams","c0_Run1,c0_Run2,c1_Run1,c1_Run2,c2_Run1,c2_Run2");
	}
	else if(bkgType == 3)
	{
		w.extendSet("nuisParams","c0_Run1,c0_Run2,c1_Run1,c1_Run2,"
		            "c2_Run1,c2_Run2,c3_Run1,c3_Run2");
	}
	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,alpha1_Run1,alpha2_Run1,nBkg_Run1,slope_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,alpha1_Run2,alpha2_Run2,nBkg_Run2,slope_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff2,lambdaEff2,nXib2,n1405_Run2,n1520_Run2");

	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,alpha1_Run1,alpha2_Run1,nBkg_Run1,a0_Run1,a1_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,alpha1_Run2,alpha2_Run2,nBkg_Run2,a0_Run2,a1_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff2,lambdaEff2,nXib2,n1405_Run2,n1520_Run2");

	if(strncmp(option,"mcFit",5))//Define global obs if not MC fit
	{
		// w.defineSet("globObs","gsigmaEff1,glambdaEff1,gnXib1,gsigmaEff2,"
		//             "glambdaEff2,gnXib2"); //define set of global observables
		w.defineSet("globObs","geff_ratio1,gnXib1,geff_ratio2,"
		            "gnXib2,geff_ratio_1405_1,geff_ratio_1405_2,"
		            "geff_ratio_1520_1,geff_ratio_1520_2,"
		            "geff_ratio_1600_1,geff_ratio_1600_2");   //define set of global observables
	}
	//*******************************************************************

	//First fit background shapes to data above peak
	w.pdf("Bkg_Run1")->fitTo(*(ds[0]),Range(5900,myHigh));
	w.pdf("Bkg_Run2")->fitTo(*(ds[1]),Range(5900,myHigh));

	TCanvas *SB_Run1      = new TCanvas("SB_Run1","SB_Run1",1200,800);
	RooPlot *SBframe_Run1 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),5900,myHigh,(myHigh-5900)/4);

	(ds[0])->plotOn(SBframe_Run1,Name("SBdata_Run1"),DataError(RooAbsData::Poisson));
	(w.pdf("Bkg_Run1"))->plotOn(SBframe_Run1,Name("SBfit_Run1"));
	SBframe_Run1->Draw();

	TCanvas *SB_Run2      = new TCanvas("SB_Run2","SB_Run2",1200,800);
	RooPlot *SBframe_Run2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),5900,myHigh,(myHigh-5900)/4);

	(ds[1])->plotOn(SBframe_Run2,Name("SBdata_Run2"),DataError(RooAbsData::Poisson));
	(w.pdf("Bkg_Run2"))->plotOn(SBframe_Run2,Name("SBfit_Run2"));
	SBframe_Run2->Draw();

	if(bkgType == 0)
	{
		w.factory(Form("Gaussian::slope_Run1_constraint(gslope_Run1[%f,-11.0,-2.0],slope_Run1,%f)",w.var("slope_Run1")->getValV(),w.var("slope_Run1")->getError()));
		w.var("gslope_Run1")->setConstant();
		w.factory(Form("Gaussian::slope_Run2_constraint(gslope_Run2[%f,-11.0,-2.0],slope_Run2,%f)",w.var("slope_Run2")->getValV(),w.var("slope_Run2")->getError()));
		w.var("gslope_Run2")->setConstant();

		w.extendSet("globObs","gslope_Run1,gslope_Run2");

		// auto vars1 = w.pdf("Bkg_Run1")->getVariables();
		// /  / auto p1 = (RooRealVar*)vars1->find("slope_Run1");
		// ((RooRealVar*)vars1->find("slope_Run1"))->setConstant(kTRUE);
		// // p1->setConstant(kTRUE);
		// delete vars1;
		//
		// auto vars2 = w.pdf("Bkg_Run2")->getVariables();
		// // auto p2 = (RooRealVar*)vars2->find("slope_Run2");
		// ((RooRealVar*)vars2->find("slope_Run2"))->setConstant(kTRUE);
		// // p2->setConstant(kTRUE);
		// delete vars2;
	}
	else if(bkgType == 1)
	{
		w.factory(Form("Gaussian::c0_Run1_constraint(gc0_Run1[%f,-2.0,2.0],c0_Run1,%f)",w.var("c0_Run1")->getValV(),w.var("c0_Run1")->getError()));
		w.var("gc0_Run1")->setConstant();
		w.factory(Form("Gaussian::c1_Run1_constraint(gc1_Run1[%f,-1.0,1.0],c1_Run1,%f)",w.var("c1_Run1")->getValV(),w.var("c1_Run1")->getError()));
		w.var("gc1_Run1")->setConstant();
		w.factory(Form("Gaussian::c0_Run2_constraint(gc0_Run2[%f,-2.0,2.0],c0_Run2,%f)",w.var("c0_Run2")->getValV(),w.var("c0_Run2")->getError()));
		w.var("gc0_Run2")->setConstant();
		w.factory(Form("Gaussian::c1_Run2_constraint(gc1_Run2[%f,-1.0,1.0],c1_Run2,%f)",w.var("c1_Run2")->getValV(),w.var("c1_Run2")->getError()));
		w.var("gc1_Run2")->setConstant();

		w.extendSet("globObs","gc0_Run1,gc1_Run1,gc0_Run2,gc1_Run2");
		// auto vars1 = w.pdf("Bkg_Run1")->getVariables();
		// // auto p1 = (RooRealVar*)vars1->find("slope_Run1");
		// ((RooRealVar*)vars1->find("c0_Run1"))->setConstant(kTRUE);
		// ((RooRealVar*)vars1->find("c1_Run1"))->setConstant(kTRUE);
		// // p1->setConstant(kTRUE);
		// delete vars1;
		//
		// auto vars2 = w.pdf("Bkg_Run2")->getVariables();
		// // auto p2 = (RooRealVar*)vars2->find("slope_Run2");
		// ((RooRealVar*)vars2->find("c0_Run2"))->setConstant(kTRUE);
		// ((RooRealVar*)vars2->find("c1_Run2"))->setConstant(kTRUE);
		// // p2->setConstant(kTRUE);
		// delete vars2;
	}
	else if(bkgType == 2)
	{
		w.factory(Form("Gaussian::c0_Run1_constraint(gc0_Run1[%f,-2.0,2.0],c0_Run1,%f)",w.var("c0_Run1")->getValV(),w.var("c0_Run1")->getError()));
		w.var("gc0_Run1")->setConstant();
		w.factory(Form("Gaussian::c1_Run1_constraint(gc1_Run1[%f,-1.0,1.0],c1_Run1,%f)",w.var("c1_Run1")->getValV(),w.var("c1_Run1")->getError()));
		w.var("gc1_Run1")->setConstant();
		w.factory(Form("Gaussian::c2_Run1_constraint(gc2_Run1[%f,-1.0,1.0],c2_Run1,%f)",w.var("c2_Run1")->getValV(),w.var("c2_Run1")->getError()));
		w.var("gc2_Run1")->setConstant();
		w.factory(Form("Gaussian::c0_Run2_constraint(gc0_Run2[%f,-2.0,2.0],c0_Run2,%f)",w.var("c0_Run2")->getValV(),w.var("c0_Run2")->getError()));
		w.var("gc0_Run2")->setConstant();
		w.factory(Form("Gaussian::c1_Run2_constraint(gc1_Run2[%f,-1.0,1.0],c1_Run2,%f)",w.var("c1_Run2")->getValV(),w.var("c1_Run2")->getError()));
		w.var("gc1_Run2")->setConstant();
		w.factory(Form("Gaussian::c2_Run2_constraint(gc2_Run2[%f,-1.0,1.0],c2_Run2,%f)",w.var("c2_Run2")->getValV(),w.var("c2_Run2")->getError()));
		w.var("gc2_Run2")->setConstant();

		w.extendSet("globObs","gc0_Run1,gc1_Run1,gc2_Run1,gc0_Run2,gc1_Run2,gc2_Run2");
		// auto vars1 = w.pdf("Bkg_Run1")->getVariables();
		// // auto p1 = (RooRealVar*)vars1->find("slope_Run1");
		// ((RooRealVar*)vars1->find("c0_Run1"))->setConstant(kTRUE);
		// ((RooRealVar*)vars1->find("c1_Run1"))->setConstant(kTRUE);
		// ((RooRealVar*)vars1->find("c2_Run1"))->setConstant(kTRUE);
		//
		// // p1->setConstant(kTRUE);
		// delete vars1;
		//
		// auto vars2 = w.pdf("Bkg_Run2")->getVariables();
		// // auto p2 = (RooRealVar*)vars2->find("slope_Run2");
		// ((RooRealVar*)vars2->find("c0_Run2"))->setConstant(kTRUE);
		// ((RooRealVar*)vars2->find("c1_Run2"))->setConstant(kTRUE);
		// ((RooRealVar*)vars2->find("c2_Run2"))->setConstant(kTRUE);
		//
		// // p2->setConstant(kTRUE);
		// delete vars2;
	}
	else if(bkgType == 3)
	{
		w.factory(Form("Gaussian::c0_Run1_constraint(gc0_Run1[%f,-2.0,2.0],c0_Run1,%f)",w.var("c0_Run1")->getValV(),w.var("c0_Run1")->getError()));
		w.var("gc0_Run1")->setConstant();
		w.factory(Form("Gaussian::c1_Run1_constraint(gc1_Run1[%f,-1.0,1.0],c1_Run1,%f)",w.var("c1_Run1")->getValV(),w.var("c1_Run1")->getError()));
		w.var("gc1_Run1")->setConstant();
		w.factory(Form("Gaussian::c2_Run1_constraint(gc2_Run1[%f,-1.0,1.0],c2_Run1,%f)",w.var("c2_Run1")->getValV(),w.var("c2_Run1")->getError()));
		w.var("gc2_Run1")->setConstant();
		w.factory(Form("Gaussian::c3_Run1_constraint(gc3_Run1[%f,-1.0,1.0],c3_Run1,%f)",w.var("c3_Run1")->getValV(),w.var("c3_Run1")->getError()));
		w.var("gc3_Run1")->setConstant();
		w.factory(Form("Gaussian::c0_Run2_constraint(gc0_Run2[%f,-2.0,2.0],c0_Run2,%f)",w.var("c0_Run2")->getValV(),w.var("c0_Run2")->getError()));
		w.var("gc0_Run2")->setConstant();
		w.factory(Form("Gaussian::c1_Run2_constraint(gc1_Run2[%f,-1.0,1.0],c1_Run2,%f)",w.var("c1_Run2")->getValV(),w.var("c1_Run2")->getError()));
		w.var("gc1_Run2")->setConstant();
		w.factory(Form("Gaussian::c2_Run2_constraint(gc2_Run2[%f,-1.0,1.0],c2_Run2,%f)",w.var("c2_Run2")->getValV(),w.var("c2_Run2")->getError()));
		w.var("gc2_Run2")->setConstant();
		w.factory(Form("Gaussian::c3_Run2_constraint(gc3_Run2[%f,-1.0,1.0],c3_Run2,%f)",w.var("c3_Run2")->getValV(),w.var("c3_Run2")->getError()));
		w.var("gc3_Run2")->setConstant();
		w.extendSet("globObs","gc0_Run1,gc1_Run1,gc2_Run1,gc3_Run1,gc0_Run2,gc1_Run2,gc2_Run2,gc3_Run2");
		// auto vars1 = w.pdf("Bkg_Run1")->getVariables();
		// // auto p1 = (RooRealVar*)vars1->find("slope_Run1");
		// ((RooRealVar*)vars1->find("c0_Run1"))->setConstant(kTRUE);
		// ((RooRealVar*)vars1->find("c1_Run1"))->setConstant(kTRUE);
		// ((RooRealVar*)vars1->find("c2_Run1"))->setConstant(kTRUE);
		// ((RooRealVar*)vars1->find("c3_Run1"))->setConstant(kTRUE);
		//
		// // p1->setConstant(kTRUE);
		// delete vars1;
		//
		// auto vars2 = w.pdf("Bkg_Run2")->getVariables();
		// // auto p2 = (RooRealVar*)vars2->find("slope_Run2");
		// ((RooRealVar*)vars2->find("c0_Run2"))->setConstant(kTRUE);
		// ((RooRealVar*)vars2->find("c1_Run2"))->setConstant(kTRUE);
		// ((RooRealVar*)vars2->find("c2_Run2"))->setConstant(kTRUE);
		// ((RooRealVar*)vars2->find("c3_Run2"))->setConstant(kTRUE);
		// // p2->setConstant(kTRUE);
		// delete vars2;
	}

	if(bkgType == 0)
	{
		w.factory("PROD::model_const1(model1,eff_ratio_constraint1,nXib_constraint1,"
		          "eff_ratio_1405_constraint1,eff_ratio_1520_constraint1,"
		          "eff_ratio_1600_constraint1,slope_Run1_constraint)"); //Multiply model by constraint terms
		w.factory("PROD::model_const2(model2,eff_ratio_constraint2,nXib_constraint2,"
		          "eff_ratio_1405_constraint2,eff_ratio_1520_constraint2,"
		          "eff_ratio_1600_constraint2,slope_Run2_constraint)"); //Multiply model by constraint terms
	}
	else if(bkgType == 1)
	{
		w.factory("PROD::model_const1(model1,eff_ratio_constraint1,nXib_constraint1,"
		          "eff_ratio_1405_constraint1,eff_ratio_1520_constraint1,"
		          "eff_ratio_1600_constraint1,c0_Run1_constraint,c1_Run1_constraint)"); //Multiply model by constraint terms
		w.factory("PROD::model_const2(model2,eff_ratio_constraint2,nXib_constraint2,"
		          "eff_ratio_1405_constraint2,eff_ratio_1520_constraint2,"
		          "eff_ratio_1600_constraint2,c0_Run2_constraint,c1_Run2_constraint)"); //Multiply model by constraint terms
	}
	else if(bkgType == 2)
	{
		w.factory("PROD::model_const1(model1,eff_ratio_constraint1,nXib_constraint1,"
		          "eff_ratio_1405_constraint1,eff_ratio_1520_constraint1,"
		          "eff_ratio_1600_constraint1,c0_Run1_constraint,c1_Run1_constraint,"
		          "c2_Run1_constraint)");                               //Multiply model by constraint terms
		w.factory("PROD::model_const2(model2,eff_ratio_constraint2,nXib_constraint2,"
		          "eff_ratio_1405_constraint2,eff_ratio_1520_constraint2,"
		          "eff_ratio_1600_constraint2,c0_Run2_constraint,c1_Run2_constraint,"
		          "c2_Run2_constraint)");                               //Multiply model by constraint terms
	}
	else if(bkgType == 3)
	{
		w.factory("PROD::model_const1(model1,eff_ratio_constraint1,nXib_constraint1,"
		          "eff_ratio_1405_constraint1,eff_ratio_1520_constraint1,"
		          "eff_ratio_1600_constraint1,c0_Run1_constraint,c1_Run1_constraint,"
		          "c2_Run1_constraint,c3_Run1_constraint)");                               //Multiply model by constraint terms
		w.factory("PROD::model_const2(model2,eff_ratio_constraint2,nXib_constraint2,"
		          "eff_ratio_1405_constraint2,eff_ratio_1520_constraint2,"
		          "eff_ratio_1600_constraint2,c0_Run2_constraint,c1_Run2_constraint,"
		          "c2_Run2_constraint,c3_Run2_constraint)");                               //Multiply model by constraint terms
	}

	RooAbsPdf* model_const1 = w.pdf("model_const1"); // get the model
	RooAbsPdf* model_const2 = w.pdf("model_const2"); // get the model


	//***********************MAKE COMBINED DATASET************************************
	RooCategory sample("sample","sample");
	sample.defineType("run1");
	sample.defineType("run2");

	RooAbsData *combData;
	if(isBinned)
	{
		// Construct combined dataset in (x,sample)
		combData = new RooDataHist("combData","combined data",*(w.var("Lb_DTF_M_JpsiLConstr")),
		                           Index(sample),Import("run1",*myhist[0]),
		                           Import("run2",*myhist[1]));
	}
	else
	{
		// Construct combined dataset in (x,sample)
		combData = new RooDataSet("combData","combined data",*(w.var("Lb_DTF_M_JpsiLConstr")),
		                          Index(sample),Import("run1",*ds_unb[0]),
		                          Import("run2",*ds_unb[1]));
	}

	combData->Print();
	// Construct a simultaneous pdf using category sample as index
	RooSimultaneous simPdf("simPdf","simultaneous pdf",sample);
	// Associate model with the physics state and model_ctl with the control state
	simPdf.addPdf(*model_const1,"run1");
	simPdf.addPdf(*model_const2,"run2");
	w.import(*combData);
	w.import(simPdf);
	// w.import(sample);
	w.defineSet("obs","Lb_DTF_M_JpsiLConstr,sample");
	cout<<"Done importing combData and simPdf"<<endl;
	w.Print();
	//***********************MAKE NLL************************************
	// RooAbsReal* nll = w.pdf("model_const")->createNLL(ds);
	//
	// new TCanvas();
	// RooPlot* frame = w.var("R")->frame();
	//
	// nll->plotOn(frame,ShiftToZero()); //the ShiftToZero option puts the minimum at 0 on the y-axis
	// frame->Draw();
	//
	// RooAbsReal* pll = nll->createProfile(*w.set("poi"));
	// // RooAbsReal* pll = nll->createProfile(*w.var("nLb"));
	// pll->plotOn(frame,LineColor(kRed));
	// frame->Draw();
	// frame->SetAxisRange(0,10,"Y");
	//
	// RooMinuit mini(*nll);
	// // mini.setVerbose(1);
	// mini.setStrategy(2);
	// mini.minos(*w.set("poi")); //could give no arg list and it will calculate minos errors for all the parameters
	// // mini.minos(*w.var("nLb")); //could give no arg list and it will calculate minos errors for all the parameters
	// // mini.improve();
	//
	// RooFitResult* res = mini.save("myResult","My Result");
	// cout<<"STATUS = "<<res->status()<<endl; //should be 0 for success!
	//
	// cout<<"MIN "<<w.var("R")->getVal()<<" ERR LO"<<w.var("R")->getErrorLo()<<" ERR HI"<<w.var("R")->getErrorHi()<<endl;
	//
	// RooFormulaVar p_mu("p_mu","p_{#mu} using asymptotic formulae","TMath::Prob(2*@0,1.)",RooArgList(*pll));
	//
	// TCanvas *c2 = new TCanvas();
	// RooPlot* frame1 = w.var("R")->frame();
	// p_mu.plotOn(frame1,LineColor(kGreen));
	// // frame1->Draw();
	//
	// c2->cd();
	//
	// Double_t limit = p_mu.findRoot(*(w.var("R")),-5.0,0.0,0.05);
	// cout<<"UL ="<<limit<<endl;
	//
	// RooOneSidedProfileLL opll("opll","opll",*pll); //pll is the pointer to your profilell function made earlier
	// opll.plotOn(frame1,LineColor(kMagenta));
	//
	// RooFormulaVar p_mu_q("p_mu_q","p_{#mu} from one-sided pll asymptotic formulae","ROOT::Math::normal_cdf_c(sqrt(2*@0),1.)",RooArgList(opll));
	// p_mu_q.plotOn(frame1,LineColor(kCyan));
	//
	// frame->GetYaxis()->SetRangeUser(-0.1,1.1);
	// frame1->Draw();
	//
	// TLine l; l.SetLineStyle(2);
	// l.DrawLine(-10,0.05,0,0.05);
	//
	// Double_t limit1 = p_mu_q.findRoot(*(w.var("R")),-5.0,0.0,0.05);
	// cout<<"One sided UL ="<<limit1<<endl;
	//
	// yields.push_back(0.0);
	// return yields;

	//************************DO THE FIT***********************

	RooFitResult *res = simPdf.fitTo(*combData,Minos(*w.set("poi")),Extended(), Save(), Hesse(false), Strategy(1), PrintLevel(0));
	//*******************************************************************

	w.factory(Form("Gaussian::nLb_Run1_syst(gnLb_Run1[%f,1000,8000],nLb_Run1,%f)",w.var("nLb_Run1")->getValV(),0.009*(w.var("nLb_Run1")->getValV())));
	w.var("gnLb_Run1")->setConstant();
	w.factory(Form("Gaussian::nLb_Run2_syst(gnLb_Run2[%f,1000,24000],nLb_Run2,%f)",w.var("nLb_Run2")->getValV(),0.008*(w.var("nLb_Run2")->getValV())));
	w.var("gnLb_Run2")->setConstant();

	w.extendSet("globObs","gnLb_Run1,gnLb_Run2");

	w.factory(Form("Gaussian::nLb_Run1_syst(gnLb_Run1[%f,1000,8000],nLb_Run1,%f)",w.var("nLb_Run1")->getValV(),0.009*(w.var("nLb_Run1")->getValV())));
	w.var("gnLb_Run1")->setConstant();
	w.factory(Form("Gaussian::nLb_Run2_syst(gnLb_Run2[%f,1000,24000],nLb_Run2,%f)",w.var("nLb_Run2")->getValV(),0.008*(w.var("nLb_Run2")->getValV())));
	w.var("gnLb_Run2")->setConstant();

	w.extendSet("globObs","gnLb_Run1,gnLb_Run2");

	w.factory("PROD::model_const_new1(model_const1,nLb_Run1_syst)");
	w.factory("PROD::model_const_new2(model_const2,nLb_Run2_syst)");

	RooAbsPdf* model_const_new1 = w.pdf("model_const_new1");
	RooAbsPdf* model_const_new2 = w.pdf("model_const_new2");

	RooSimultaneous simPdf_new("simPdf_new","simultaneous pdf",sample);

	simPdf_new.addPdf(*model_const_new1,"run1");
	simPdf_new.addPdf(*model_const_new2,"run2");

	w.import(simPdf_new);

	cout<<"*********After importing new simPDF**********"<<endl;
	w.Print();
	//*********************PLOTTING STUFF*********************************************
	TCanvas* c_run1 = new TCanvas("Run1","Run1", 1200, 800);

	RooPlot *frame_run1 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),myLow,myHigh,nbins);
	frame_run1->SetTitle("Run1 Fit");
	// frame_run1->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");

	combData->plotOn(frame_run1,Name("data_Run1"),Cut("sample==sample::run1"),DataError(RooAbsData::Poisson));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Name("fit_Run1"));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb_Run1"))),Name("lb_Run1"),LineColor(kMagenta+2));
	// simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb1_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	// simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb2_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Bkg_Run1"))),LineColor(kRed),Name("bkg_Run1"));
	//  if(xibflag!=0)
	if(strncmp(option,"mcFit",5))//if not MC fit
	{
		simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("XIB1"))),LineColor(kGreen),Name("xib_Run1"));
		if(lst1405flag)
			simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("LST1405_Run1"))),LineColor(kGreen+2),LineStyle(kDashed),Name("lst1405_Run1"));
		if(lst1520flag)
			simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("LST1520_Run1"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1520_Run1"));
		if(lst1600flag)
			simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("LST1600_Run1"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1600_Run1"));
		// if(chic1flag)
		//      simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("chic1_Run1"))),LineColor(kMagenta+2),LineStyle(kDashed),Name("chic1_Run1"));
		//simPdf.plotOn(frame_run1,Components(lst1810shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1810"));
		// if(sigmaflag!=0)
		simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("SIG1"))),LineColor(kBlack),Name("sig_Run1"));
		simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("lstLump_Run1"))),LineColor(kRed+2),LineStyle(kDashed),Name("misclst_Run1"));

		frame_run1->GetYaxis()->SetRangeUser(0,60);
	}
	// Double_t chiSquare1 = frame_run1->chiSquare("fit_run1","data_Run1");
	// cout<<"chi square1/dof = "<<chiSquare1<<endl;
	RooArgSet *allpar_run1 = simPdf.getParameters(*(ds[0]));
	RooArgSet *floatpar_run1 = (RooArgSet*)allpar_run1->selectByAttrib("Constant",kFALSE);
	floatpar_run1->Print();
	int floatpars_run1 = (floatpar_run1->selectByAttrib("Constant",kFALSE))->getSize() - 1;//-1 because sample also gets included in this list
	cout<<"run1 float pars = "<<floatpars_run1<<endl;
	Double_t chi2_run1 = frame_run1->chiSquare("fit_Run1","data_Run1",floatpars_run1);
	cout<<"chi square2/dof = "<<chi2_run1<<endl;

	Int_t fit_ndof_run1 = nbins - floatpars_run1;
	cout<<"chi square2 = "<<chi2_run1*fit_ndof_run1<<endl;

	///////////
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.2,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.2);

	pad1->SetGridx();
	pad1->SetGridy();
	pad2->SetGridx();
	pad2->SetGridy();

	pad1->SetBottomMargin(0.04);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	c_run1->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();
	gPad->SetTopMargin(0.06);
	pad1->Update();

	frame_run1->Draw();

	TLatex l_run1;
	l_run1.SetTextSize(0.025);
	l_run1.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run1));

	c_run1->Modified();

	auto legend_run1 = new TLegend(0.1,0.7,0.3,0.9);
	legend_run1->SetTextSize(0.025);
	legend_run1->AddEntry("data_Run1","Data","lp");
	legend_run1->AddEntry("fit_Run1","Total Fit","l");
	legend_run1->AddEntry("lb_Run1","J/#psi #Lambda shape","l");
	legend_run1->AddEntry("bkg_Run1","Comb. Bkg. shape","l");
	if(strncmp(option,"mcFit",5))//if not MC fit
	{
		legend_run1->AddEntry("sig_Run1","J/#psi #Sigma shape","l");
		legend_run1->AddEntry("xib_Run1","J/#psi #Xi shape","l");
		if(lst1405flag)
		{
			legend_run1->AddEntry("lst1405_Run1","J/#psi #Lambda(1405) shape","l");
		}
		if(lst1520flag)
		{
			legend_run1->AddEntry("lst1520_Run1","J/#psi #Lambda(1520) shape","l");
		}
		if(lst1600flag)
		{
			legend_run1->AddEntry("lst1600_Run1","J/#psi #Lambda(1600) shape","l");
		}
		// if(chic1flag)
		// {
		//      legend_run1->AddEntry("chic1_Run1","#chi_{c1} #Lambda shape","l");
		// }
		// if(lst1810flag)
		//      legend_run1->AddEntry("lst1810","J/#psi #Lambda(1810) shape","l");
		legend_run1->AddEntry("misclst_Run1","misc. J/#psi #Lambda* shapes","l");
	}
	legend_run1->Draw("same");

	c_run1->Update();

	// Pull distribution
	RooPlot *frame_run1x2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),myLow,myHigh,nbins);
	// RooPlot *framex2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), myLow,5800,100);
	RooHist* hpull_run1 = frame_run1->pullHist("data_Run1","fit_Run1");
	frame_run1x2->addPlotable(hpull_run1,"P");
	hpull_run1->SetLineColor(kBlack);
	hpull_run1->SetMarkerColor(kBlack);
	frame_run1x2->SetTitle(0);
	frame_run1x2->GetYaxis()->SetTitle("Pull");
	frame_run1x2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame_run1x2->GetYaxis()->SetTitleSize(0.15);
	frame_run1x2->GetYaxis()->SetLabelSize(0.15);
	frame_run1x2->GetXaxis()->SetTitleSize(0.15);
	frame_run1x2->GetXaxis()->SetLabelSize(0.15);
	frame_run1x2->GetYaxis()->CenterTitle();
	frame_run1x2->GetYaxis()->SetTitleOffset(0.25);
	frame_run1x2->GetXaxis()->SetTitleOffset(0.75);
	frame_run1x2->GetYaxis()->SetNdivisions(505);
	frame_run1x2->GetYaxis()->SetRangeUser(-5.0,5.0);
	pad2->cd();
	frame_run1x2->Draw();

	c_run1->cd();
	// pad1->cd();

	cout<<"Run1 Pull Mean Y = "<<hpull_run1->GetMean(2)<<endl;
	cout<<"Run1 Pull RMS Y = "<<hpull_run1->GetRMS(2)<<endl;

	TCanvas* c_run2 = new TCanvas("Run2","Run2", 1200, 800);

	RooPlot *frame_run2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),myLow,myHigh,nbins);

	frame_run2->SetTitle("Run2 Fit");
	frame_run2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	combData->plotOn(frame_run2,Name("data_Run2"),Cut("sample==sample::run2"),DataError(RooAbsData::Poisson));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Name("fit_Run2"));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb_Run2"))),Name("lb_Run2"),LineColor(kMagenta+2));
	// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb1_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb2_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Bkg_Run2"))),LineColor(kRed),Name("bkg_Run2"));
	if(strncmp(option,"mcFit",5))//if not MC fit
	{
		//  if(xibflag!=0)
		simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("XIB2"))),LineColor(kGreen),Name("xib_Run2"));
		if(lst1405flag)
			simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("LST1405_Run2"))),LineColor(kGreen+2),LineStyle(kDashed),Name("lst1405_Run2"));
		if(lst1520flag)
			simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("LST1520_Run2"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1520_Run2"));
		if(lst1600flag)
			simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("LST1600_Run2"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1600_Run2"));
		// if(chic1flag)
		// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("chic1_Run2"))),LineColor(kMagenta+2),LineStyle(kDashed),Name("chic1_Run2"));
		//simPdf.plotOn(frame_run2,Components(lst1810shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1810"));
		// if(sigmaflag!=0)
		simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("SIG2"))),LineColor(kBlack),Name("sig_Run2"));
		simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("lstLump_Run2"))),LineColor(kRed+2),LineStyle(kDashed),Name("misclst_Run2"));
		frame_run2->GetYaxis()->SetRangeUser(0,160);
	}
	// Double_t chiSquare1 = frame_run2->chiSquare("fit_run2","data_run2");
	// cout<<"chi square1/dof = "<<chiSquare1<<endl;
	RooArgSet *floatpar_run2 = simPdf.getParameters(*(ds[1]));
	floatpar_run2->Print();
	int floatpars_run2 = (floatpar_run2->selectByAttrib("Constant",kFALSE))->getSize() -1;//-1 because sample also gets included in this list
	cout<<"run2 float pars = "<<floatpars_run2<<endl;
	Double_t chi2_run2 = frame_run2->chiSquare("fit_Run2","data_Run2",floatpars_run2);
	cout<<"chi square2/dof = "<<chi2_run2<<endl;

	Int_t fit_ndof_run2 = nbins - floatpars_run2;
	cout<<"chi square2 = "<<chi2_run2*fit_ndof_run2<<endl;

	///////////
	TPad *pad3 = new TPad("pad3","pad3",0.0,0.2,1.0,1.0);
	TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,0.2);

	pad3->SetGridx();
	pad3->SetGridy();
	pad4->SetGridx();
	pad4->SetGridy();

	pad3->SetBottomMargin(0.04);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.3);
	pad4->SetBorderMode(0);
	pad3->SetBorderMode(0);
	c_run1->SetBorderMode(0);
	pad4->Draw();
	pad3->Draw();
	pad3->cd();
	gPad->SetTopMargin(0.06);
	pad3->Update();

	frame_run2->Draw();
	// c_run2->cd();

	auto legend_run2 = new TLegend(0.1,0.7,0.3,0.9);
	legend_run2->SetTextSize(0.025);
	legend_run2->AddEntry("data_Run2","Data","lp");
	legend_run2->AddEntry("fit_Run2","Total Fit","l");
	legend_run2->AddEntry("lb_Run2","J/#psi #Lambda shape","l");
	legend_run2->AddEntry("bkg_Run2","Comb. Bkg. shape","l");
	if(strncmp(option,"mcFit",5))//if not MC fit
	{
		legend_run2->AddEntry("sig_Run2","J/#psi #Sigma shape","l");
		legend_run2->AddEntry("xib_Run2","J/#psi #Xi shape","l");
		if(lst1405flag)
		{
			legend_run2->AddEntry("lst1405_Run2","J/#psi #Lambda(1405) shape","l");
		}
		if(lst1520flag)
		{
			legend_run2->AddEntry("lst1520_Run2","J/#psi #Lambda(1520) shape","l");
		}
		if(lst1600flag)
		{
			legend_run2->AddEntry("lst1600_Run2","J/#psi #Lambda(1600) shape","l");
		}
		// if(chic1flag)
		// {
		//      legend_run2->AddEntry("chic1_Run2","#chi_{c1} #Lambda shape","l");
		// }
		// if(lst1810flag)
		//      legend_run2->AddEntry("lst1810","J/#psi #Lambda(1810) shape","l");
		legend_run2->AddEntry("misclst_Run2","misc. J/#psi #Lambda* shapes","l");
	}
	legend_run2->Draw("same");

	TLatex l_run2;
	l_run2.SetTextSize(0.025);
	l_run2.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run2));

	c_run2->Update();

	// Pull distribution
	RooPlot *frame_run2x2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),myLow,myHigh,nbins);
	// RooPlot *framex2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), myLow,5800,100);
	RooHist* hpull_run2 = frame_run2->pullHist("data_Run2","fit_Run2");
	frame_run2x2->addPlotable(hpull_run2,"P");
	hpull_run2->SetLineColor(kBlack);
	hpull_run2->SetMarkerColor(kBlack);
	frame_run2x2->SetTitle(0);
	frame_run2x2->GetYaxis()->SetTitle("Pull");

	frame_run2x2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame_run2x2->GetYaxis()->SetTitleSize(0.15);
	frame_run2x2->GetYaxis()->SetLabelSize(0.15);
	frame_run2x2->GetXaxis()->SetTitleSize(0.15);
	frame_run2x2->GetXaxis()->SetLabelSize(0.15);
	frame_run2x2->GetYaxis()->CenterTitle();
	frame_run2x2->GetYaxis()->SetTitleOffset(0.25);
	frame_run2x2->GetXaxis()->SetTitleOffset(0.75);
	frame_run2x2->GetYaxis()->SetNdivisions(505);
	frame_run2x2->GetYaxis()->SetRangeUser(-5.0,5.0);
	pad4->cd();
	frame_run2x2->Draw();

	c_run2->cd();

	Double_t chi2_global = chi2_run1*fit_ndof_run1 + chi2_run2*fit_ndof_run2;
	Int_t ndof_global = fit_ndof_run1 + fit_ndof_run2;

	Double_t chi2_ndof_global = chi2_global/ndof_global;

	cout<<"****************************"<<endl;
	cout<<"Global Fit chi2/dof = "<<chi2_ndof_global<<endl;
	cout<<"****************************"<<endl;

	if(calcUL)
	{
		//*************ROOSTATS MODEL CONFIG*********************************

		cout<<"Starting Model Config"<<endl;
		RooStats::ModelConfig mc("ModelConfig",&w);
		mc.SetPdf(*(w.pdf("simPdf_new")));
		mc.SetParametersOfInterest(*w.set("poi"));
		mc.SetObservables(*w.set("obs"));
		mc.SetGlobalObservables(*w.set("globObs"));
		mc.SetNuisanceParameters(*w.set("nuisParams"));

		// import model in the workspace
		w.import(mc);
		//*******************************************************************

		//***************ASIMOV DATASET*************************************
		cout<<"Starting Asimov Dataset"<<endl;
		RooDataSet *asimovData = (RooDataSet*)RooStats::AsymptoticCalculator::GenerateAsimovData(*(mc.GetPdf()),*(mc.GetObservables()));
		RooPlot *frame_asim = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), myLow,myHigh,nbins);

		TCanvas *asim = new TCanvas();
		asimovData->plotOn(frame_asim);
		frame_asim->Draw();

		mc.SetSnapshot(*(w.var("R")));

		RooStats::ModelConfig *bkgOnlyModel = mc.Clone();
		bkgOnlyModel->SetName("bkgOnlyModel");
		Double_t oldval = w.var("R")->getVal();

		w.var("R")->setVal(0);
		bkgOnlyModel->SetSnapshot(*(w.var("R")));

		w.var("R")->setVal(oldval);

		w.import(*bkgOnlyModel);
		// w.Write();
		//*******************************************************************

		//**************Fit like Matt****************************************

/*	cout<<"Starting Fit like Matt"<<endl;
        RooAbsReal* nll = w.pdf("simPdf")->createNLL(*combData);
        RooMinimizer min_nll(*nll);
        min_nll.setErrorLevel(0.5);
        min_nll.setStrategy(2);

        Double_t initnll = nll->getVal();

        RooArgSet *pars = w.pdf("simPdf")->getParameters( *combData );
        // RooArgSet *initpars = (RooArgSet*)pars->snapshot();
        // cout<<"YYYY"<<pars<<endl;
        // print "YYYY", initpars

        // TIterator *it = initpars->createIterator();
        // RooRealVar *par = (RooRealVar*)it->Next();
        // cout<<"YYYY All par list"<<endl;
        // while(par)
        // {
        //      cout<<"YYYY"<<par<<par->GetName()<<par->getVal()<<par->getError()<<endl;
        //      par = (RooRealVar*)it->Next();
        // }
        min_nll.minimize("Minuit","Migrad");

        min_nll.minos(*(w.var("R")));

        RooFitResult *res1 = min_nll.save();

        Double_t finalnll = nll->getVal();*/

		Double_t origval = w.var("R")->getVal();
		Double_t errhi   = w.var("R")->getErrorHi();
		Double_t errlo   = w.var("R")->getErrorLo();

		//**************Upper Limit Calculation*****************************

		// cout<<"Starting Upper Limit Calculation"<<endl;

		//**************AsymptoticCalculator********************************
		// RooStats::AsymptoticCalculator hc(*combData, *bkgOnlyModel, mc, true);
		// hc.SetOneSided(true);
		// RooStats::HypoTestInverter calc(hc);
		//
		// calc.SetConfidenceLevel(0.95);
		// // calc.SetFixedScan(100,-7, -1);
		// calc.SetFixedScan(20,0, 1000);
		// calc.UseCLs(true);
		// // calc.SetVerbose(true);
		// HypoTestInverterResult *r = calc.GetInterval();
		//
		// cout<<"origval = "<<origval<<" error myHigh = "
		//     <<errhi<<" error myLow = "<<errlo<<endl;
		// cout<<"UL = "<<r->UpperLimit()<<" +/- "<<r->UpperLimitEstimatedError()<<endl;
		//
		// new TCanvas();
		// RooStats::HypoTestInverterPlot plot("hti_plot", "hti_plot", r);
		// plot.Draw("CLb 2CL");

		//******************************************************************

		//**************FrequentistCalculator*******************************
		// RooStats::FrequentistCalculator hc(combData, *bkgOnlyModel, mc);
		// TestStatSampler *toymcs = hc.GetTestStatSampler();
		//
		// RooStats::ProfileLikelihoodTestStat teststat(*(mc.GetPdf()));
		// teststat.SetOneSided(true);
		//
		// teststat.SetPrintLevel(1);
		// teststat.SetStrategy(0);
		//
		// toymcs->SetTestStatistic(teststat);
		// toymcs->SetGenerateBinned(false);
		// toymcs->SetUseMultiGen(false);
	}

	//Save Canvases
	//
	// if(isBinned)
	// {
	//      c_run1->SaveAs(Form("../plots/data/JpsiLambda/run1/Fit_HypatiaSig_ExpBkg_%d_%d_%dMeVBins.pdf",myLow,myHigh,binwidth));
	//      c_run2->SaveAs(Form("../plots/data/JpsiLambda/run2/Fit_HypatiaSig_ExpBkg_%d_%d_%dMeVBins.pdf",myLow,myHigh,binwidth));
	// }
	// else
	// {
	//      c_run1->SaveAs(Form("../plots/data/JpsiLambda/run1/Fit_HypatiaSig_ExpBkg_%d_%d_unbinned.pdf",myLow,myHigh));
	//      c_run2->SaveAs(Form("../plots/data/JpsiLambda/run2/Fit_HypatiaSig_ExpBkg_%d_%d_unbinned.pdf",myLow,myHigh));
	// }

	// if(isBinned)
	// {
	//      fileName = Form("../rootFiles/dataFiles/JpsiLambda/ModelConfigs/MyModel_HypatiaSig_ExpBkg_%d_%d_%dMeVBins.root",myLow,myHigh,binwidth);
	// }
	// else
	// {
	//      fileName = Form("../rootFiles/dataFiles/JpsiLambda/ModelConfigs/MyModel_HypatiaSig_ExpBkg_%d_%d_unbinned.root",myLow,myHigh);
	// }

	w.writeToFile(fileName,true);
	cout << "workspace written to file " << fileName << endl;

	// gROOT->ProcessLine(Form(".x StandardHypoTestInvDemo.C(\"%s\",\"w\",\"ModelConfig\",\"bkgOnlyModel\",\"combData\",2,2,true,20,100,2000,100,false,0,%d,%d)",fileName,myLow,myHigh));
	// gSystem->RedirectOutput(0);
}
