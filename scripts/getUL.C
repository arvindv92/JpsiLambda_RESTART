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

void getUL(const char *option)
{
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

	Float_t xib_syst             = 0.05;

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

	Int_t Nobs[2]     = {0,0};
	Int_t Ncomb[2]    = {0,0};
	Float_t Nsig[2]   = {0.0,0.0};
	Float_t Nlb[2]    = {0.0,0.0};
	Float_t Nxib[2]   = {0.0,0.0};
	Float_t xibERR[2] = {0.0,0.0};
	Float_t lbERR[2];
	Float_t nsigERR[2] = {0.0,0.0};
	Int_t nentries[2];

	RooAbsReal* xibInt[2];
	Double_t xibINT[2];

	RooAbsReal* lbInt[2];
	Double_t lbINT[2];
	if(!strncmp(option,"best",4)) //Set parameters for best fit
	{
		isoVersion[0] = "v0";//"v1";
		isoVersion[1] = "v0";//"v0";

		isoConf[0] = 2;//1;
		isoConf[1] = 1;//2;

		bdtConf_nonZero[0] = 2;//2;
		bdtConf_nonZero[1] = 2;//1;

		bdtConf_Zero[0] = 2;//1;
		bdtConf_Zero[1] = 1;//1;

		bdtCut_nonZero[0] = 0.475-0.1;//0.375 - 0.1;
		bdtCut_nonZero[1] = 0.555-0.1;//0.535 - 0.1;

		bdtCut_Zero[0] = 0.365;//0.285;
		bdtCut_Zero[1] = 0.455;//0.415;
	}
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
		                   " GetNorm(%d, \"%s\", %d, %d, %d, %f, %f)\'",
		                   run, isoVersion[i], isoConf[i], bdtConf_nonZero[i],
		                   bdtConf_Zero[i], bdtCut_nonZero[i], bdtCut_Zero[i]));

		ifstream infile(Form("../logs/mc/JpsiXi/run%d/xibNorm_log.txt",run));

		infile>>xibnorm_LL[i];         // Get normalization
		infile>>xibnorm_LL_staterr[i]; // Get statistical error on norm. This is the fit error

		infile>>xibnorm_LL_wt[i];         // Get weighted normalization
		infile>>xibnorm_LL_staterr_wt[i]; // Get weighted statistical error on norm. This is the fit error

		xibnorm_LL_systerr[i] = xibnorm_LL[i]*xib_syst;
		xibnorm_LL_err[i]     = sqrt(pow(xibnorm_LL_staterr[i],2) + pow(xibnorm_LL_systerr[i],2)); //Combining stat and syst in quadrature

		xibnorm_LL_systerr_wt[i] = xibnorm_LL_wt[i]*xib_syst;
		xibnorm_LL_err_wt[i]     = sqrt(pow(xibnorm_LL_staterr_wt[i],2) + pow(xibnorm_LL_systerr_wt[i],2)); //Combining stat and syst in quadrature

		xibnorm_LL[i]         = xibnorm_LL[i]*2; //ACCOUNT FOR XIB0
		xibnorm_LL_err[i]     = xibnorm_LL_err[i]*1.414; //ACCOUNT FOR XIB0

		xibnorm_LL_wt[i]         = xibnorm_LL_wt[i]*2; //ACCOUNT FOR XIB0
		xibnorm_LL_err_wt[i]     = xibnorm_LL_err_wt[i]*1.414; //ACCOUNT FOR XIB0

		cout<<"************************************************"<<endl;
		cout<<"The UNWEIGHTED LL Xib normalization for Run "<<run
		    <<" is "<<xibnorm_LL[i]<<" +/- "<<xibnorm_LL_err[i]<<endl;

		cout<<"The WEIGHTED LL Xib normalization for Run "<<run
		    <<" is "<<xibnorm_LL_wt[i]<<" +/- "<<xibnorm_LL_err_wt[i]<<endl;
		cout<<"************************************************"<<endl;
	}
	//************************************************

	RooDataSet* ds_xi[2];
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

		treein_xi_nonZero->SetBranchStatus("*",0);
		treein_xi_nonZero->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		treein_xi_nonZero->SetBranchStatus(Form("BDT%d",bdtConf_nonZero[i]),1);
		treein_xi_nonZero->SetBranchStatus("Lb_BKGCAT",1);
		treein_xi_nonZero->SetBranchStatus("gb_wts",1);

		TFile *tempFile = new TFile("tempFile.root","RECREATE");

		TTree* treein_xi_Zero_cut    = (TTree*)treein_xi_Zero->CopyTree(Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]));
		TTree* treein_xi_nonZero_cut = (TTree*)treein_xi_nonZero->CopyTree(Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]));

		TList *list = new TList;
		list->Add(treein_xi_Zero_cut);
		list->Add(treein_xi_nonZero_cut);

		TTree *combTree = TTree::MergeTrees(list);
		combTree->SetName("combTree");

		ds_xi[i] = new RooDataSet("ds_xi","ds_xi",combTree,RooArgSet(*myVar),0,"gb_wts");
		ds_xi[i]->Print();
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

		XIB_KEYS[i] = new RooKeysPdf(Form("XIB%d",run),Form("XIB%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*(ds_xi[i]),RooKeysPdf::NoMirror);
		// XIB[i] = new RooHistPdf(Form("XIB%d",run),Form("XIB%d",run),*(w.var("Lb_DTF_M_JpsiLConstr")),*ds_xib_smooth,0);

		RooPlot *framexib = (w.var("Lb_DTF_M_JpsiLConstr"))->frame();
		framexib->SetTitle("J/#psi #Xi");
		ds_xi[i]->plotOn(framexib,Name("xibdata"));
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
		xibERR[i] = xibnorm_LL_err_wt[i]*xibINT[i];
	}
	//*********Hypatia signal shape for Lambda_b0************************

	// w.factory("RooHypatia2::Lb_Run1(Lb_DTF_M_JpsiLConstr,lambda_Run1[-2.0,-4.0,0.0],0,0,"
	//           "sigma_Run1[10.,1.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[1.7,1.0,3.0],"
	//           "2 ,a2_Run1[2.0,1.0,3.0], 2)");

	w.factory("RooHypatia2::Lb_Run1(Lb_DTF_M_JpsiLConstr,lambda_Run1[-2.0,-4.0,0.0],0,0,"
	          "sigma_Run1[10.,1.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[1.7,1.0,4.0],"
	          "2 ,a2_Run1[2.0,1.0,3.0], 2)");

	w.factory("RooHypatia2::Lb_Run2(Lb_DTF_M_JpsiLConstr,lambda_Run2[-2.5,-4.0,0.0],0,0,"
	          "sigma_Run2[10.,1.,20.], mean_Run2[5619.6,5619,5621], a1_Run2[1.5,1.0,3.0],"
	          "2 ,a2_Run2[1.5,1.0,3.0], 2)");

	cout<<"Done defining J/psi Lambda Hypatia shapes"<<endl;
	//*******************************************************************

	//*********Continuum backgruond****************

	// if(bkgType == 0)
	// {
	cout<<"*****UING EXPONENTIAL BKG SHAPE*****"<<endl;
	w.factory("Exponential::Bkg_Run1(Lb_DTF_M_JpsiLConstr,tau_Run1[-0.001,-0.01,-0.0000001])");
	w.factory("Exponential::Bkg_Run2(Lb_DTF_M_JpsiLConstr,tau_Run2[-0.001,-0.01,-0.0000001])");
	// }
	// else if(bkgType == 1)
	// {
	//      cout<<"*****UING 2nd ORDER CHEBYCHEV BKG SHAPE*****"<<endl;
	//      w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0], c1_Run1[0.0,-1.0,1.0]})");
	//      w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-2.0,2.0], c1_Run2[0.0,-1.0,1.0]})");
	// }
	// else if(bkgType == 2)
	// {
	//      cout<<"*****UING 3rd ORDER CHEBYCHEV BKG SHAPE*****"<<endl;
	//      w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0], c1_Run1[0.0,-1.0,1.0], c2_Run1[0.0,-1.0,1.0]})");
	//      w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-2.0,2.0], c1_Run2[0.0,-1.0,1.0], c2_Run2[0.0,-1.0,1.0]})");
	// }
	// else if(bkgType == 3)
	// {
	//      cout<<"*****UING 4th ORDER CHEBYCHEV BKG SHAPE*****"<<endl;
	//      w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[0.0,-2.0,2.0], c1_Run1[0.0,-1.0,1.0], c2_Run1[0.0,-1.0,1.0], c3_Run1[0.0,-1.0,1.0]})");
	//      w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[0.0,-2.0,2.0], c1_Run2[0.0,-1.0,1.0], c2_Run2[0.0,-1.0,1.0], c3_Run2[0.0,-1.0,1.0]})");
	// }

	cout<<"Done defining cont. bkg shapes"<<endl;

	// w.factory("Chebychev::Bkg_Run1(Lb_DTF_M_JpsiLConstr, {c0_Run1[-0.5,-2.0,2.0], c1_Run1[0.5,-1.0,1.0], c2_Run1[0.,-1.0,1.0]})");
	// w.factory("Chebychev::Bkg_Run2(Lb_DTF_M_JpsiLConstr, {c0_Run2[-0.5,-2.0,2.0], c1_Run2[-0.5,-1.0,1.0], c2_Run2[0.,-1.0,1.0]})");

	//*******************************************************************

	w.factory("nLb_Run1[5000,1000,7000]");
	w.factory("nLb_Run2[17000,1000,22000]");
	w.factory(Form("nBkg_Run1[2000,1,%d]",nentries[0]));
	w.factory(Form("nBkg_Run2[4000,1,%d]",nentries[1]));

	w.factory("SUM:model1(nLb_Run1*Lb_Run1 + nBkg_Run1*Bkg_Run1)");
	w.factory("SUM:model2(nLb_Run2*Lb_Run2 + nBkg_Run2*Bkg_Run2)");

	RooAbsPdf* model1 = w.pdf("model1"); // get the model
	RooAbsPdf* model2 = w.pdf("model2"); // get the model

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
	// Associate model with the physics state and model_ctl with the control state
	simPdf.addPdf(*model1,"run1");
	simPdf.addPdf(*model2,"run2");
	w.import(*combData);
	w.import(simPdf);

	//************************DO THE FIT***********************
	RooFitResult *res = simPdf.fitTo(*combData,Extended(), Save(), Hesse(false), Strategy(1), PrintLevel(0), Range(5500,5700));
	//*******************************************************************

	//*********************PLOTTING STUFF*********************************************
	TCanvas* c_run1 = new TCanvas("Run1","Run1", 1200, 800);

	RooPlot *frame_run1 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	frame_run1->SetTitle("Run1 Fit");
	// frame_run1->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");

	combData->plotOn(frame_run1,Name("data_Run1"),Cut("sample==sample::run1"),DataError(RooAbsData::Poisson));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Name("fit_Run1"));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb_Run1"))),Name("lb_Run1"),LineColor(kMagenta+2));
	// simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb1_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	// simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Lb2_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,*combData),Components(*(w.pdf("Bkg_Run1"))),LineColor(kRed),Name("bkg_Run1"));

	frame_run1->GetYaxis()->SetRangeUser(0,60);
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
	legend_run1->Draw("same");

	c_run1->Update();

	// Pull distribution
	RooPlot *frame_run1x2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	// RooPlot *framex2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), low,5800,100);
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

	RooPlot *frame_run2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);

	frame_run2->SetTitle("Run2 Fit");
	frame_run2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	combData->plotOn(frame_run2,Name("data_Run2"),Cut("sample==sample::run2"),DataError(RooAbsData::Poisson));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Name("fit_Run2"));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb_Run2"))),Name("lb_Run2"),LineColor(kMagenta+2));
	// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb1_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Lb2_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,*combData),Components(*(w.pdf("Bkg_Run2"))),LineColor(kRed),Name("bkg_Run2"));

	frame_run2->GetYaxis()->SetRangeUser(0,160);
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
	legend_run2->Draw("same");

	TLatex l_run2;
	l_run2.SetTextSize(0.025);
	l_run2.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run2));

	c_run2->Update();

	// Pull distribution
	RooPlot *frame_run2x2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	// RooPlot *framex2 = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")), low,5800,100);
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

	RooAbsPdf *Lb_Run1 = w.pdf("Lb_Run1");
	RooAbsPdf *Lb_Run2 = w.pdf("Lb_Run2");

	lbInt[0] = w.pdf("Lb_Run1")->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
	lbINT[0] = lbInt[0]->getValV();
	Nlb[0]   = lbINT[0]*(w.var("nLb_Run1")->getVal());
	lbInt[1] = w.pdf("Lb_Run2")->createIntegral(*myVar,NormSet(*myVar),Range("signal_window"));
	lbINT[1] = lbInt[1]->getValV();
	Nlb[1]   = lbINT[1]*(w.var("nLb_Run2")->getVal());
	lbERR[0] = (lbInt[0]->getValV())*(w.var("nLb_Run1")->getError());
	lbERR[1] = (lbInt[1]->getValV())*(w.var("nLb_Run2")->getError());

	for(Int_t i=0; i<=1; i++)
	{
		Nsig[i] = Nobs[i] - Nlb[i] - Nxib[i] - Ncomb[i];
		nsigERR[i] = sqrt( Nobs[i] + Ncomb[i] + pow(lbERR[i],2) + pow(xibERR[i],2));

		cout<<"RUN"<<i+1<<" : nsig = "<<Nsig[i]<<" +/- "<<nsigERR[i]<<endl;
	}

}
