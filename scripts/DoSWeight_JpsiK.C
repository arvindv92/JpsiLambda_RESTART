#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TString.h"
#include "TROOT.h"
#include <iostream>

// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel_JpsiK(RooWorkspace*, Int_t, Int_t, Int_t, Int_t);
void AddData_JpsiK(RooWorkspace*, Int_t, TTree*);
void DoSPlot_JpsiK(RooWorkspace*, Int_t, TTree*, TH1D*);
Double_t MakePlots_JpsiK(RooWorkspace*, Int_t);

Double_t DoSWeight_JpsiK(Int_t run = 1)
{
	ROOT::EnableImplicitMT();
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	gSystem->RedirectOutput(Form("logs/data/B2JpsiK/run%d/sPlot.txt",run),"w");

	gROOT->ProcessLine(".x lhcbStyle.C");
	TStopwatch sw;
	sw.Start();

	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"==> Starting DoSWeight_JpsiK: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"********************************"<<endl;

	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;

	const char *rootFolder = "";
	TString inFileName     = "";
	TString outFileName    = "";
	TString trainFileName  = "";

	Int_t entries_init = 0;
	Int_t binwidth     = 2;//trying 2 MeV bins
	Int_t lowRange     = 5200;
	Int_t highRange    = 5360;//5200-5360 MeV window for sWeighting
	Int_t nBins        = (Int_t)(highRange-lowRange)/binwidth;

	rootFolder    = Form("rootFiles/dataFiles/B2JpsiK/run%d", run);
	inFileName    = Form("%s/jpsik.root", rootFolder);
	outFileName   = Form("%s/jpsik_withsw.root", rootFolder);

	fileIn = TFile::Open(inFileName,"READ");
	if (!fileIn)
	{
		cout << "ERROR: could not open data file" << endl;
		exit(1);
	}
	treeIn = (TTree*)fileIn->Get("MyTuple");
	entries_init       = treeIn->GetEntries();
	cout<<"Incoming Entries = "<<entries_init<<endl;

	treeIn->Draw(Form("B_Mass>>hMass(%d,%d,%d)",nBins,lowRange,highRange),"","goff");
	TH1D *hMass = (TH1D*)gDirectory->Get("hMass");

	fileOut = new TFile(outFileName,"RECREATE");
	treeOut = treeIn->CloneTree(0);

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"sWeighted Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	// Create a new workspace to manage the project.
	RooWorkspace *wSpace = new RooWorkspace("myWS");

	// add the signal and background models to the workspace.
	AddModel_JpsiK(wSpace, lowRange, highRange, entries_init, run);

	//add the data
	AddData_JpsiK(wSpace, run, treeIn);

	// do sPlot.
	//This wil make a new dataset with sWeights added for every event.
	DoSPlot_JpsiK(wSpace, run, treeOut, hMass);

	// Make some plots showing the discriminating variable and
	// the control variable after unfolding.
	Double_t chi2 = MakePlots_JpsiK(wSpace, run);

	//Set Aliases so that this file can be used for isolation BDT training
	treeOut->SetAlias("PT","K_PT");
	treeOut->SetAlias("VCHI2DOF","B_ENDVERTEX_CHI2/3");
	treeOut->SetAlias("MINIPCHI2","K_MINIPCHI2");
	treeOut->SetAlias("IPCHI2","K_IPCHI2_ORIVX");

	//inspect the workspace
	wSpace->Print();

	//save workspace in ROOT file
	wSpace->writeToFile(Form("%s/wSpace_sPlot.root",rootFolder));

	fileOut->cd();
	treeOut->Write("",TObject::kOverwrite);
	fileOut->Close();

	// cleanup
	delete wSpace;

	sw.Stop();
	cout<<"*****Done with DoSWeight_JpsiK*****"<<endl;
	sw.Print();
	gSystem->RedirectOutput(0);
	return chi2;
}
void AddModel_JpsiK(RooWorkspace *ws, Int_t lowRange, Int_t highRange,
                    Int_t nEntries, Int_t run)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting AddModel()"<<endl;
	cout<<"******************************************"<<endl;

	TStopwatch sw;
	sw.Start();

	Int_t binwidth     = 2;//trying 2 MeV bins
	Int_t nBins        = (Int_t)(highRange-lowRange)/binwidth;

	// make a RooRealVar for the observables
	RooRealVar B_Mass("B_DTF_M_JpsiConstr", "M_{J/#psi K^{+}}", lowRange, highRange,"MeV");//discriminating variable
	B_Mass.setBins(nBins);

	cout << "make signal model" << endl;

	// SIGNAL MODEL
	RooRealVar mB("mB", "B Mass", 5280.0, 5275.0, 5285.0, "MeV");
	RooRealVar sigmaB1("sigmaB1", "Width of Gaussian1",6.3,1.0,20.0,"MeV");
	RooRealVar sigmaB2("sigmaB2", "Width of Gaussian2",13.0,1.0,20.0,"MeV");
	RooGaussian sig1("sig1","signal pdf1",B_Mass,mB,sigmaB1);
	RooGaussian sig2("sig2","signal pdf2",B_Mass,mB,sigmaB2);
	RooRealVar f1("f1","fraction",0.5,0.01,1.0);

	RooAddPdf mBModel("mBModel","Sig",RooArgList(sig1,sig2),f1);

	// BACKGROUND MODEL

	cout << "make background model" << endl;

	RooRealVar c0("c0","coefficient #0",-0.1,-10.,10.);
	RooRealVar c1("c1","coefficient #1",-0.1,-10.,10.);
	RooChebychev bkg("bkg","background pdf",B_Mass,RooArgList(c0,c1));

	// COMBINED MODEL
	Int_t ul_sig = (run == 1) ? (250000) : (475000);

	RooRealVar sigYield("sigYield","fitted yield for sig",0,ul_sig);
	RooRealVar bkgYield("bkgYield","fitted yield for bkg",0,nEntries);

	cout << "make full model" << endl;
	RooAddPdf model("model","signal+background models",
	                RooArgList(mBModel, bkg),
	                RooArgList(sigYield,bkgYield));

	cout << "Importing model into workspace" << endl;

	ws->import(model);
	ws->Print();

	sw.Stop();
	cout<<"******************************************"<<endl;
	cout<<"AddModel() complete"<<endl; sw.Print();
	cout<<"******************************************"<<endl;
}
void AddData_JpsiK(RooWorkspace *ws, Int_t run, TTree *treeIn)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting AddData()"<<endl;
	cout<<"******************************************"<<endl;

	TStopwatch sw;
	sw.Start();

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("B_DTF_M_JpsiConstr",1);
	RooRealVar B_Mass = *(ws->var("B_DTF_M_JpsiConstr"));

	RooArgSet s(B_Mass);

	cout<<"Importing data into workspace. This could take a while. Sit tight"<<endl;
	RooDataSet data("data","dataset with x",treeIn,s);

	// import data into workspace
	ws->import(data, Rename("data"));
	ws->Print();

	sw.Stop();
	cout<<"******************************************"<<endl;
	cout<<"Finishing AddData()"<<endl; sw.Print();
	cout<<"******************************************"<<endl;
}
void DoSPlot_JpsiK(RooWorkspace *ws, Int_t run, TTree *treeOut, TH1D *hMass)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting DoSPlot()"<<endl;
	cout<<"******************************************"<<endl;

	TStopwatch sw;
	sw.Start();

	Float_t SIGW = 0., BKGW = 0., bMASS = 0.;
	TBranch *sigW = nullptr;
	TBranch *bkgW = nullptr;

	// get what we need out of the workspace to do the fit
	RooAbsPdf *model = ws->pdf("model");

	RooRealVar *sigYield = ws->var("sigYield");
	RooRealVar *bkgYield = ws->var("bkgYield");
	RooRealVar *B_Mass   = ws->var("B_DTF_M_JpsiConstr");
	RooDataSet *data     = (RooDataSet*) ws->data("data");

	RooDataHist *rdh = new RooDataHist("rdh","",*B_Mass,hMass);
	ws->import(*rdh,Rename("data_binned"));

	//fit model to binned data.
	RooFitResult *r = model->fitTo(*rdh,Hesse(kTRUE),Strategy(2),Save(true));//Fit to binned dataset
	cout<<"***** RooFitResult from fit to binned data*****"<<endl;
	r->Print();
	cout<<"***********************************************"<<endl;

	//I don't think I need to do another fit with all non yield params fixed.
	//The SPlot class takes care of fixing non yield params automatically.

	// RooArgSet *myVars;
	// myVars = model->getParameters(*rdh);
	// myVars->Print("v");
	// myVars->setAttribAll("Constant",kTRUE);
	// sigYield->setConstant(kFALSE);
	// bkgYield->setConstant(kFALSE);
	//
	// //Fit to binned dataset again with all params except yields fixed to fitted value
	// RooFitResult *r1 = model->fitTo(*rdh,Extended(),Hesse(kTRUE), Save(true));
	// cout<<"***** RooFitResult from fit to binned data with nonYield params fixed*****"<<endl;
	// r1->Print();
	// cout<<"***********************************************"<<endl;

	cout<<"*****Doing actual SPlot now*****"<<endl;
	//****************************************************************************
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, model,
	                                             RooArgList(*sigYield,*bkgYield));
	//****************************************************************************
	cout<<"***********************************************"<<endl;

	data->Print("v");

	// Check that our weights have the desired properties
	cout << "Check SWeights:" << endl;

	cout << endl <<  "Yield of signal is "
	     << sigYield->getVal() << ".  From sWeights it is "
	     << sData->GetYieldFromSWeight("sigYield") << endl;


	cout << "Yield of background is "
	     << bkgYield->getVal() << ".  From sWeights it is "
	     << sData->GetYieldFromSWeight("bkgYield") << endl;

	for(Int_t i=0; i < 10; i++)
	{
		cout << "Signal Weight   " << sData->GetSWeight(i,"sigYield")
		     << "   Background Weight   " << sData->GetSWeight(i,"bkgYield")
		     << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
		     << endl;
	}
	cout << endl;

	// import this new dataset with sWeights
	cout << "import new dataset with sWeights" << endl;
	ws->import(*data, Rename("dataWithSWeights"));
	cout <<" after splot import data set looks like"<<endl;
	ws->Print();

	//Add SW & BW branches to output tree
	sigW  = treeOut->Branch("SW",&SIGW,"SIGW/F");
	bkgW  = treeOut->Branch("BW",&BKGW,"BKGW/F");

	Int_t dsentries = data->numEntries();
	cout<<"No. of entries in dataset = "<<dsentries<<endl;

	cout<<"************************************************"<<endl;
	cout<<"Writing output tree with unbinned sweighted data"<<endl;
	cout<<"************************************************"<<endl;

	for(Int_t i = 0; i < dsentries; i++)
	{
		bMASS = (data->get(i))->getRealValue("B_DTF_M_JpsiConstr");
		SIGW  = (data->get(i))->getRealValue("sigYield_sw");
		BKGW  = (data->get(i))->getRealValue("bkgYield_sw");

		sigW->Fill();
		bkgW->Fill();
	}

	sw.Stop();
	cout<<"******************************************"<<endl;
	cout<<"Done with DoSPlot()"<<endl; sw.Print();
	cout<<"******************************************"<<endl;
}
Double_t MakePlots_JpsiK(RooWorkspace* ws, Int_t run)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting MakePlots()"<<endl;
	cout<<"******************************************"<<endl;
	TStopwatch sw;
	sw.Start();

	RooAbsPdf *model   = ws->pdf("model");
	RooDataSet *data   = (RooDataSet*) ws->data("data");
	RooRealVar *B_Mass = ws->var("B_DTF_M_JpsiConstr");

	RooAbsPdf *sig1 = ws->pdf("sig1");
	RooAbsPdf *sig2 = ws->pdf("sig2");
	RooAbsPdf *sig  = ws->pdf("sig");
	RooAbsPdf *bkg  = ws->pdf("bkg");

	TCanvas *fitCanvas = new TCanvas("fitCanvas","",1200,800);

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.2,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.2);

	pad1->SetGridx();
	pad1->SetGridy();
	pad2->SetGridx();
	pad2->SetGridy();

	pad1->SetBottomMargin(0.0);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.45);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	fitCanvas->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();
	//	gPad->SetTopMargin(0.06);
	pad1->Update();

	RooPlot* frame = B_Mass->frame();
	frame->GetXaxis()->SetTitle("m_{J/#psi K^{+}}[MeV/#it{c}^{2}]");
	frame->GetYaxis()->SetTitle("Candidates/(4 MeV/#it{c}^{2})");
	data->plotOn(frame,Name("data"),DataError(RooAbsData::SumW2) );
	model->plotOn(frame,Name("fit") );
	model->plotOn(frame,Components(*sig1),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig2),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig),LineStyle(kDashed), LineColor(kRed),Name("sig"));
	model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen),Name("bkg"));

	RooArgSet *allpar = model->getParameters(*data);
	RooArgSet *floatpar = (RooArgSet*)allpar->selectByAttrib("Constant",kFALSE);
	floatpar->Print();
	int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	cout<<"no of float pars = "<<floatpars<<endl;
	Double_t chi2 = frame->chiSquare("fit","data",floatpars);
	cout<<"chi square2/dof = "<<chi2<<endl;

	frame->Draw();

	TLatex l;
	l.SetTextSize(0.04);
	l.DrawLatexNDC(0.7,0.3,Form("#chi^{2}/ndf = %.2f",chi2));

	fitCanvas->Modified();

	auto legend = new TLegend(0.7,0.5,0.9,0.9);
	legend->SetTextSize(0.04);
	legend->AddEntry("data","Data","lp");
	legend->AddEntry("fit","Total Fit","l");
	legend->AddEntry("bkg","Comb. Bkg.","l");
	legend->AddEntry("sig","B^{+} #rightarrow J/#psi K^{+}","l");
	legend->Draw("same");

	fitCanvas->Update();

	// Pull distribution
	RooPlot *framex2 = B_Mass->frame();
	RooHist* hpull = frame->pullHist("data","fit");
	framex2->addPlotable(hpull,"P");
	framex2->GetXaxis()->SetTitle("m_{J/#psi K^{+}}[MeV/#it{c}^{2}]");
	framex2->GetYaxis()->SetTitle("Pull");
	framex2->GetYaxis()->SetTitleSize(0.25);
	framex2->GetYaxis()->SetLabelSize(0.25);
	framex2->GetXaxis()->SetTitleSize(0.25);
	framex2->GetXaxis()->SetLabelSize(0.25);

	hpull->SetLineColor(kBlack);
	hpull->SetMarkerColor(kBlack);
	framex2->SetTitle(0);
	// framex2->GetYaxis()->SetTitleSize(0.15);
	// framex2->GetYaxis()->SetLabelSize(0.15);
	// framex2->GetXaxis()->SetTitleSize(0.15);
	// framex2->GetXaxis()->SetLabelSize(0.15);
	framex2->GetYaxis()->CenterTitle();
	framex2->GetYaxis()->SetTitleOffset(0.25);
	framex2->GetXaxis()->SetTitleOffset(0.75);
	framex2->GetYaxis()->SetNdivisions(505);
	framex2->GetYaxis()->SetRangeUser(-4.0,4.0);
	pad2->cd();
	framex2->Draw();

	fitCanvas->cd();
	// pad1->cd();

	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;

	fitCanvas->SaveAs(Form("plots/ANA/sPlot_JpsiK_Run%d.pdf",run));
	cout<<"******************************************"<<endl;
	cout<<"Finished MakePlots()"<<endl; sw.Print();
	cout<<"******************************************"<<endl;

	return chi2;
}
