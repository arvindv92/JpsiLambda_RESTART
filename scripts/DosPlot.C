/********************************
   Author : Aravindhan V.
   The purpose of this script is to sWeight data in 5200-6000 GeV, and write out weighted data
 *********************************/
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TPaveLabel.h"
#include "TStopwatch.h"
// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*, Int_t, Int_t, Int_t);
void AddData(RooWorkspace*,Int_t,const char*, TString, TTree*);
void DoSPlot(RooWorkspace*,Int_t,const char*, TString, TString, TTree*);
void MakePlots(RooWorkspace*);

void DosPlot(Int_t run = 1, Int_t trackType = 3, Bool_t logFlag = false)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{
	TStopwatch sw;
	sw.Start();

	const char* logFileName = "", *type = "";
	TString inFileName = "", outFileName = "", trainFileName = "";
	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;
	Int_t entries_init = 0, entries_massWindow = 0, lowRange = 5200, highRange = 6000;
	Double_t Lb_DTF_M_JpsiLConstr = 0.0;

	type          = (trackType == 3) ? ("LL") : ("DD");
	logFileName   = (trackType == 3) ? ("sPlot_LL_log.txt") : ("sPlot_DD_log.txt");
	inFileName    = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type);
	outFileName   = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw.root",run,type);
	trainFileName = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_forIsoTraining.root",run,type);

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.
	if(logFlag)
	{
		gROOT->ProcessLine((TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName)).Data());
	}
	cout<<"********************************"<<endl;
	cout<<"==> Starting DosPlot: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"********************************"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Processing Run "<<run<<" "<<type<<endl;
	cout<<"******************************************"<<endl;

	fileIn             = TFile                     ::Open(inFileName,"READ");
	treeIn             = (TTree*)fileIn->Get("MyTuple");

	entries_init       = treeIn->GetEntries();
	entries_massWindow = treeIn->GetEntries(TString::Format("Lb_DTF_M_JpsiLConstr > %d && Lb_DTF_M_JpsiLConstr < %d",lowRange,highRange));

	fileOut = TFile::Open(outFileName,"RECREATE");
	treeOut = (TTree*)treeIn->CloneTree(0);

	//	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	treeIn->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&Lb_DTF_M_JpsiLConstr);
	cout<<"Incoming Entries = "<<entries_init<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"sWeighted Output file = "<<fileOut->GetName()<<endl;
	cout<<"Training Output file = "<<trainFileName<<endl;
	cout<<"******************************************"<<endl;

	cout<<"Making mass cut on "<<entries_init<<" entries"<<endl;
	for (Int_t i = 0; i < entries_init; i++)
	{
		if(i%50000 == 0) cout<<i<<endl;
		treeIn->GetEntry(i);
		if(Lb_DTF_M_JpsiLConstr > lowRange && Lb_DTF_M_JpsiLConstr < highRange)
		{
			treeOut->Fill();
		}
	}

	// Create a new workspace to manage the project.
	RooWorkspace* wSpace = new RooWorkspace("myWS");

	// add the signal and background models to the workspace.
	// Inside this function you will find a discription our model.
	AddModel(wSpace, lowRange, highRange, entries_massWindow);

	// add some toy data to the workspace
	AddData(wSpace, run, type, inFileName, treeIn);

	// inspect the workspace if you wish
	//  wSpace->Print();

	// do sPlot.
	//This wil make a new dataset with sWeights added for every event.
	DoSPlot(wSpace, run, type, outFileName, trainFileName, treeOut);

	fileOut->cd();
	treeOut->Write("",TObject::kOverwrite);
	fileOut->Close();

	// cleanup
	delete wSpace;

	sw.Stop();
	cout << "==> DosPlot is done! Check fit status! sWeights FTW!: "; sw.Print();
	if(logFlag) gROOT->ProcessLine(".>");
}

void AddModel(RooWorkspace* ws = nullptr, Int_t lowRange = 5200, Int_t highRange = 6000, Int_t nEntries = 0)
{
	cout<<"Starting AddModel()"<<endl;
	// Make models for signal and background

	// make a RooRealVar for the observables
	RooRealVar Lb_DTF_M_JpsiLConstr("Lb_DTF_M_JpsiLConstr", "M_{J/#psi#Lambda}", lowRange, highRange,"MeV");//discriminating variable

	Int_t nbins = (Int_t)(highRange-lowRange)/4;
	Lb_DTF_M_JpsiLConstr.setBins(nbins);

	// SIGNAL MODEL
	cout << "Making signal model" << endl;
	RooRealVar mean("mean","Gaussian Mean",5619.42,5619.0,5621.0);
	RooRealVar sigma1("sigma1","Gaussian sigma1",7.0,5.0,20.0);
	RooRealVar sigma2("sigma2","Gaussian sigma2",7.0,3.0,20.0);

	RooGaussian sig1("sig1","Gaussian signal1",Lb_DTF_M_JpsiLConstr,mean,sigma1);
	RooGaussian sig2("sig2","Gaussian signal2",Lb_DTF_M_JpsiLConstr,mean,sigma2);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5,0.01,1.0);

	RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);

	// BACKGROUND MODEL
	cout << "Making background model" << endl;
	// RooRealVar c0("c0","c0",-0.1,-1.0,1.0);
	// RooRealVar c1("c1","c1",0.1,-1.0,1.0);
	//
	// RooChebychev bkg("bkg","bkg",Lb_DTF_M_JpsiLConstr,RooArgList(c0,c1));

	RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
	RooExponential bkg("bkg","Exponential Bkg",Lb_DTF_M_JpsiLConstr,tau);

	// COMBINED MODEL	cout << "Making full model" << endl;

	RooRealVar sigYield("sigYield","fitted yield for sig",5000,0, 15000);
	RooRealVar bkgYield("bkgYield","fitted yield for bkg",(nEntries/2),0, nEntries);

	RooAddPdf model("model","signal+background models", RooArgList(sig, bkg), RooArgList(sigYield,bkgYield));

	cout << "Importing model into workspace" << endl;

	ws->import(model);
	ws->Print();

	cout<<"AddModel() complete"<<endl;
}

void AddData(RooWorkspace* ws = nullptr, Int_t run = 1, const char* type = "LL", TString inFileName = "", TTree *treeIn = nullptr)
{
	cout<<"Starting AddData()"<<endl;

	// get what we need out of the workspace to make toy data
	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	RooRealVar Lb_DTF_M_JpsiLConstr = *(ws->var("Lb_DTF_M_JpsiLConstr"));

	//Import the data
	cout<<"Importing data. This could take a while. Sit tight"<<endl;

	RooArgSet s(Lb_DTF_M_JpsiLConstr);
	RooDataSet data("data","dataset with x",treeIn,s);

	ws->import(data, Rename("data"));
	ws->Print();
	cout<<"Finishing AddData()"<<endl;
}

void DoSPlot(RooWorkspace* ws = nullptr, Int_t run = 1, const char* type = "LL", TString outFileName = "", TString trainFileName = "", TTree *treeOut = nullptr)
{
	TFile *fileOut_training = nullptr, *fileOut_nonZeroTracks = nullptr, *fileOut_ZeroTracks = nullptr;
	TTree *treeOut_training = nullptr, *treeOut_nonZeroTracks = nullptr, *treeOut_ZeroTracks = nullptr;

	Float_t SIGW = 0., BKGW = 0.;
	TBranch *sigW = nullptr,*bkgW = nullptr;

	cout<<"Starting DoSPlot()"<<endl;
	cout << "Calculating sWeights" << endl;

	// get what we need out of the workspace to do the fit
	RooAbsPdf *model                 = ws->pdf("model");
	RooAbsPdf *sig1                  = ws->pdf("sig1");
	RooAbsPdf *sig2                  = ws->pdf("sig2");
	RooAbsPdf *sig                   = ws->pdf("sig");
	RooAbsPdf *bkg                   = ws->pdf("bkg");
	RooRealVar *sigYield             = ws->var("sigYield");
	RooRealVar *bkgYield             = ws->var("bkgYield");
	RooDataSet *data                 = (RooDataSet*) ws->data("data");
	RooRealVar *Lb_DTF_M_JpsiLConstr = ws->var("Lb_DTF_M_JpsiLConstr");

	// fit the model to the data.
	RooFitResult *r = model->fitTo(*data, Extended(), Strategy(2), Save(true)); //unbinned extended ML fit

	cout<<"Starting MakePlots()"<<endl;

	// make our canvas
	TCanvas *fitCanvas = new TCanvas("sPlot","sPlot demo", 1200, 800);

	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
	pad1->SetBottomMargin(0);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.5);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	fitCanvas->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();

	gPad->SetTopMargin(0.06);
	pad1->Update();

	//plot Lb_DTF_M_JpsiLConstr for data with full model and individual componenets overlayed

	RooPlot* frame = Lb_DTF_M_JpsiLConstr->frame();
	data->plotOn(frame,Name("Hist"),DataError(RooAbsData::SumW2) );
	model->plotOn(frame,Name("curvetot") );
	model->plotOn(frame,Components(*sig1),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig2),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig),LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen));
	model->paramOn(frame, Format("NEU"), Layout(0.65,0.85,0.9));
	frame->SetTitle("Fit of model to discriminating variable");
	frame->Draw();

	////Fit chi2

	Double_t chiSquare1 = frame->chiSquare("curvetot","Hist");
	cout<<"chi square1 = "<<chiSquare1<<endl;
	RooArgSet *floatpar = model->getParameters(data);
	floatpar->Print();
	int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	Double_t chi2 = frame->chiSquare("curvetot","Hist",floatpars);
	cout<<"chi square2 = "<<chi2<<endl;

	TPaveLabel *t1 = new TPaveLabel(0.1,0.8,0.3,0.88, Form("#chi^{2}/dof = %f", chi2),"NDC");
	t1->Draw();
	frame->addObject(t1);
	pad1->Update();

	// Pull distribution
	RooPlot *framex2 = Lb_DTF_M_JpsiLConstr->frame();
	RooHist* hpull = frame->pullHist("Hist","curvetot");
	framex2->addPlotable(hpull,"P");
	hpull->SetLineColor(kBlack);
	hpull->SetMarkerColor(kBlack);
	framex2->SetTitle(0);
	framex2->GetYaxis()->SetTitle("Pull");
	framex2->GetYaxis()->SetTitleSize(0.11);
	framex2->GetYaxis()->SetLabelSize(0.15);
	framex2->GetXaxis()->SetTitleSize(0.15);
	framex2->GetXaxis()->SetLabelSize(0.15);
	framex2->GetYaxis()->CenterTitle();
	framex2->GetYaxis()->SetTitleOffset(0.4);
	framex2->GetXaxis()->SetTitleOffset(0.7);
	framex2->GetYaxis()->SetNdivisions(505);
	framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
	pad2->cd();
	framex2->Draw();

	fitCanvas->cd();

	fitCanvas->Update();
	fitCanvas->SaveAs(TString::Format("plots/fit_run%d%s.pdf",run,type));

	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS  Y = "<<hpull->GetRMS(2)<<endl;

	cout<<"***************"<<endl;
	r->Print("v");
	cout<<"***************"<<endl;
	cout<<"covQual = "<<r->covQual()<<endl;
	cout<<"numStatusHistory = "<<r->numStatusHistory()<<endl;
	for(UInt_t i=0; i<r->numStatusHistory(); i++)
	{
		cout<<r->statusCodeHistory(i)<<" "<<r->statusLabelHistory(i)<<endl;
	}

	// The sPlot technique requires that we fix the parameters
	// of the model that are not yields after doing the fit.
	RooRealVar *sigma1 = ws->var("sigma1");
	RooRealVar *sigma2 = ws->var("sigma2");
	RooRealVar *mean   = ws->var("mean");
	// RooRealVar *c0     = ws->var("c0");
	// RooRealVar *c1     = ws->var("c1");
	RooRealVar *tau    = ws->var("tau");

	mean->setConstant(kTRUE);
	sigma1->setConstant(kTRUE);
	sigma2->setConstant(kTRUE);
	// c0->setConstant(kTRUE);
	// c1->setConstant(kTRUE);
	tau->setConstant(kTRUE);

	RooMsgService::instance().setSilentMode(true);

	// Now we use the SPlot class to add SWeights to our data set
	// based on our model and our yield variables
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, model, RooArgList(*sigYield,*bkgYield) );

	// Check that our weights have the desired properties

	cout << "Check SWeights:" << endl;

	cout << endl <<  "Yield of signal is "
	     << sigYield->getVal() << ".  From sWeights it is "
	     << sData->GetYieldFromSWeight("sigYield") << endl;


	cout << "Yield of background is "
	     << bkgYield->getVal() << ".  From sWeights it is "
	     << sData->GetYieldFromSWeight("bkgYield") << endl
	     << endl;

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

	TCanvas *sWeightCanvas = new TCanvas();
	RooRealVar *sigYield_sw = ws->var("sigYield_sw");
	sigYield_sw->setRange(-1.0,1.5);
	RooPlot* frame_sigsw = sigYield_sw->frame();
	frame_sigsw->SetTitle("Signal sweights");
	data->plotOn(frame_sigsw);
	frame_sigsw->Draw();
	sWeightCanvas->Update();
	sWeightCanvas->SaveAs(TString::Format("plots/fit_run%d%s_sWeights.pdf",run,type));

	TCanvas *sWeightVsMass = new TCanvas();
	RooPlot* frame_sigsw_x = Lb_DTF_M_JpsiLConstr->frame();
	frame_sigsw_x->SetTitle("sWeights vs bmass");
	data->plotOnXY(frame_sigsw_x,YVar(*sigYield_sw),MarkerColor(kBlue));
	frame_sigsw_x->Draw();
	sWeightVsMass->Update();
	sWeightVsMass->SaveAs(TString::Format("plots/fit_run%d%s_sWeightsVsMass.pdf",run,type));

	// create weightfed data set
	RooDataSet *dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"sigYield_sw");
	RooDataSet *dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw");

	TCanvas *sWeightMass = new TCanvas();
	sWeightMass->Divide(1,2);

	sWeightMass->cd(1);

	RooPlot* frame2_1 = Lb_DTF_M_JpsiLConstr->frame();
	dataw_sig->plotOn(frame2_1, DataError(RooAbsData::SumW2));
	frame2_1->SetTitle("Lb DTF Mass Distribution (SW)");
	frame2_1->Draw();

	sWeightMass->cd(2);

	RooPlot* frame2_2 = Lb_DTF_M_JpsiLConstr->frame();
	dataw_bkg->plotOn(frame2_2, DataError(RooAbsData::SumW2));
	frame2_2->SetTitle("Lb DTF Mass Distribution (BW)");
	frame2_2->Draw();

	sWeightMass->Update();
	sWeightMass->SaveAs(TString::Format("plots/fit_run%d%s_sWeightedMass.pdf",run,type));

	Int_t dsentries = data->numEntries();
	cout<<"No. of entries in dataset = "<<dsentries<<endl;

	cout<<"************************************************"<<endl;
	cout<<"Writing output tree with unbinned sweighted data"<<endl;
	cout<<"************************************************"<<endl;

	// bMass = treeOut->Branch("B_MASS",&B_MASS,"B_MASS/D");
	sigW  = treeOut->Branch("SW",&SIGW,"SIGW/F");
	bkgW  = treeOut->Branch("BW",&BKGW,"BKGW/F");

	for(Int_t i = 0; i < dsentries; i++)
	{
		// B_MASS = (data->get(i))->getRealValue("Lb_DTF_M_JpsiLConstr");
		SIGW = (data->get(i))->getRealValue("sigYield_sw");
		BKGW = (data->get(i))->getRealValue("bkgYield_sw");

		// bMass->Fill();
		sigW->Fill();
		bkgW->Fill();
		//treeOut->Fill();
	}

	treeOut->SetBranchStatus("*",0);
	treeOut->SetBranchStatus("SW",1);
	treeOut->SetBranchStatus("BW",1);
	treeOut->SetBranchStatus("Added_H_PT",1);
	treeOut->SetBranchStatus("psi_1S_H_MINIPCHI2",1);
	treeOut->SetBranchStatus("psi_1S_H_VERTEXCHI2_NEW",1);
	treeOut->SetBranchStatus("psi_1S_H_IPCHI2_NEW",1);//Maybe also use IP_NEW?
	treeOut->SetBranchStatus("psi_1S_H_IP_NEW",1);
	treeOut->SetBranchStatus("psi_1S_H_FD_NEW",1);
	treeOut->SetBranchStatus("psi_1S_H_FDCHI2_NEW",1);
	treeOut->SetBranchStatus("psi_1S_H_VERTEX_Z_NEW",1);
	treeOut->SetBranchStatus("Added_H_TRACKCHI2",1);
	treeOut->SetBranchStatus("Added_H_GHOST",1);
	treeOut->SetBranchStatus("Lb_H_OPENING",1);
	treeOut->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

	fileOut_training = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_forIsoTraining.root",run,type),"RECREATE");
	cout<<"Copying tree with s-weights to make signal training file for isolation"<<endl;
	treeOut_training = (TTree*)treeOut->CopyTree("");

	treeOut_training->SetAlias("PT","Added_H_PT");
	treeOut_training->SetAlias("MINIPCHI2","psi_1S_H_MINIPCHI2");
	treeOut_training->SetAlias("VCHI2DOF","psi_1S_H_VERTEXCHI2_NEW");
	treeOut_training->SetAlias("IPCHI2","psi_1S_H_IPCHI2_NEW");
	treeOut_training->SetAlias("IP","psi_1S_H_IP_NEW");
	treeOut_training->SetAlias("FD","psi_1S_H_FD_NEW");
	treeOut_training->SetAlias("FDCHI2","psi_1S_H_FDCHI2_NEW");
	treeOut_training->SetAlias("TRACKORIVX_Z","psi_1S_H_VERTEXCHI2_NEW");
	treeOut_training->SetAlias("GHOSTPROB","Added_H_GHOST");
	treeOut_training->SetAlias("TRACKCHI2DOF","Added_H_TRACKCHI2");

	treeOut->SetBranchStatus("*",1);

	cout<<"Done writing output tree with unbinned sWeighted data"<<endl;

	fileOut_training->cd();
	treeOut_training->Write();
	fileOut_training->Close();

	cout<<"Done writing signal training file for isolation"<<endl;

	fileOut_nonZeroTracks = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw_nonZeroTracks.root",run,type),"RECREATE");
	cout<<"Copying tree with s-weights to make nonZeroTracks file"<<endl;
	treeOut_nonZeroTracks = (TTree*)treeOut->CopyTree("Added_n_Particles > 0");

	fileOut_nonZeroTracks->cd();
	treeOut_nonZeroTracks->Write();
	fileOut_nonZeroTracks->Close();

	cout<<"Done writing output tree with unbinned sWeighted data for nonZeroTracks"<<endl;

	fileOut_ZeroTracks = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw_ZeroTracks.root",run,type),"RECREATE");
	cout<<"Copying tree with s-weights to make ZeroTracks file"<<endl;
	treeOut_ZeroTracks = (TTree*)treeOut->CopyTree("Added_n_Particles == 0");

	fileOut_ZeroTracks->cd();
	treeOut_ZeroTracks->Write();
	fileOut_ZeroTracks->Close();

	cout<<"Done writing output tree with unbinned sWeighted data for ZeroTracks"<<endl;

	cout<<"Finishing DoSPlot()"<<endl;
}
