/********************************
   Author : Aravindhan V.
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
// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*);
void AddData(RooWorkspace*,Int_t,const char*);
void DoSPlot(RooWorkspace*,Int_t,const char*);
void MakePlots(RooWorkspace*);

void Lb_sPlot(Int_t run, Int_t trackType)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{

	Bool_t logFlag = false;
	const char* logFileName = "", *type = "";

	if(trackType == 3)
	{
		cout<<"Processing LL"<<endl;
		logFileName = "sPlot_LL_log.txt";
		type = "LL";
	}
	else if(trackType == 5)
	{
		cout<<"Processing DD"<<endl;
		logFileName = "sPlot_DD_log.txt";
		type = "DD";
	}
	if(logFlag)
	{
		gROOT->ProcessLine((TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName)).Data());
	}
	// Create a new workspace to manage the project.
	RooWorkspace* wspace = new RooWorkspace("myWS");

	// add the signal and background models to the workspace.
	// Inside this function you will find a discription our model.
	AddModel(wspace);

	// add some toy data to the workspace
	AddData(wspace, run, type);

	// inspect the workspace if you wish
	//  wspace->Print();

	// do sPlot.
	//This wil make a new dataset with sWeights added for every event.
	DoSPlot(wspace, run, type);

	// Make some plots showing the discriminating variable and
	// the control variable after unfolding.

	//MakePlots(wspace);

	// cleanup
	delete wspace;
	if(logFlag) gROOT->ProcessLine(".>");
}

void AddModel(RooWorkspace* ws)
{
	cout<<"Starting AddModel()"<<endl;
	// Make models for signal and background

	// set range of observable
	Double_t lowRange = 5200.0, highRange = 6000.0;//highRange changed from 6800

	// make a RooRealVar for the observables
	RooRealVar Lb_DTF_M_JpsiLConstr("Lb_DTF_M_JpsiLConstr", "M_{inv}", lowRange, highRange,"MeV");//discriminating variable

	Int_t nbins = (Int_t)(highRange-lowRange)/4;
	Lb_DTF_M_JpsiLConstr.setBins(nbins);

	// SIGNAL MODEL
	cout << "Making signal model" << endl;
	RooRealVar mean("mean","Gaussian Mean",5619.42,5618.0,5621.0);
	RooRealVar sigma1("sigma1","Gaussian sigma1",7.5,1.0,20.0);
	RooRealVar sigma2("sigma2","Gaussian sigma2",10.0,1.0,20.0);

	RooGaussian sig1("sig1","Gaussian signal1",Lb_DTF_M_JpsiLConstr,mean,sigma1);
	RooGaussian sig2("sig2","Gaussian signal2",Lb_DTF_M_JpsiLConstr,mean,sigma2);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5,0.1,0.9);

	RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);

	// BACKGROUND MODEL
	cout << "make background model" << endl;
	RooRealVar c0("c0","c0",-0.1,-1.0,-1.0);
	RooRealVar c1("c1","c1",0.1,-1.0,1.0);

	RooChebychev bkg("bkg","bkg",Lb_DTF_M_JpsiLConstr,RooArgList(c0,c1));

	// COMBINED MODEL
	cout << "Making full model" << endl;

	RooRealVar sigYield("sigYield","fitted yield for sig",5000,0, 15000);
	RooRealVar bkgYield("bkgYield","fitted yield for bkg",100000,0, 1.1e+06);//adjust this based on nentries

	RooAddPdf model("model","signal+background models", RooArgList(sig, bkg), RooArgList(sigYield,bkgYield));

	// interesting for debugging and visualizing the model
	//model.graphVizTree("fullModel.dot");

	cout << "Importing model into workspace" << endl;

	ws->import(model);
	ws->Print();

	cout<<"AddModel() complete"<<endl;
}

void AddData(RooWorkspace* ws, Int_t run, const char* type)
{
	cout<<"Starting AddData()"<<endl;

	TFile *filein(0);
	TTree *treein(0);
	// get what we need out of the workspace to make toy data

	RooRealVar Lb_DTF_M_JpsiLConstr = *(ws->var("Lb_DTF_M_JpsiLConstr"));

	//const char* cuts = "";
	//Import the data
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.
	filein = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type),"READ");

	treein = (TTree*)filein->Get("MyTuple");

	treein->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

	cout<<"Incoming Entries = "<<treein->GetEntries()<<endl;

	RooArgSet s(Lb_DTF_M_JpsiLConstr);

	RooDataSet data("data","dataset with x",treein,s);

	// import data into workspace
	ws->import(data, Rename("data"));
	ws->Print();
	filein->Close();
	cout<<"Finishing AddData()"<<endl;
}

void DoSPlot(RooWorkspace* ws, Int_t run, const char* type)
{
	cout<<"Starting DoSPlot()"<<endl;
	cout << "Calculating sWeights" << endl;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

	TFile *fileout(0), *fileout_training(0);
	TTree *treeout(0), *treeout_training(0);

	// get what we need out of the workspace to do the fit
	RooAbsPdf* model                 = ws->pdf("model");
	RooAbsPdf* sig1                  = ws->pdf("sig1");
	RooAbsPdf* sig2                  = ws->pdf("sig2");
	RooAbsPdf* sig                   = ws->pdf("sig");
	RooAbsPdf* bkg                   = ws->pdf("bkg");
	RooRealVar* sigYield             = ws->var("sigYield");
	RooRealVar* bkgYield             = ws->var("bkgYield");
	RooDataSet* data                 = (RooDataSet*) ws->data("data");
	RooRealVar* Lb_DTF_M_JpsiLConstr = ws->var("Lb_DTF_M_JpsiLConstr");

	// fit the model to the data.
	model->fitTo(*data, Extended(), Strategy(2));

	//Should ideally be plotting here.

	cout<<"Starting MakePlots()"<<endl;

	// make our canvas
	TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 1200, 800);

	///////////
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
	pad1->SetBottomMargin(0);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.5);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	cdata->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();

	gPad->SetTopMargin(0.06);
	pad1->Update();
	////////////

	//plot Lb_DTF_M_JpsiLConstr for data with full model and individual componenets overlayed

	RooPlot* frame = Lb_DTF_M_JpsiLConstr->frame();
	data->plotOn(frame,Name("Hist"),DataError(RooAbsData::SumW2) );
	model->plotOn(frame,Name("curvetot") );
	model->plotOn(frame,Components(*sig1),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig2),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig),LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen));

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

	// Pull distribution
	RooPlot *framex2 = Lb_DTF_M_JpsiLConstr->frame();
	RooHist* hpull = frame->pullHist("Hist","curvetot");
	framex2->addPlotable(hpull,"P");
	hpull->SetLineColor(kBlack);
	hpull->SetMarkerColor(kBlack);
	framex2->SetTitle(0);
	framex2->GetYaxis()->SetTitle("Pull");
	framex2->GetYaxis()->SetTitleSize(0.15);
	framex2->GetYaxis()->SetLabelSize(0.15);
	framex2->GetXaxis()->SetTitleSize(0.2);
	framex2->GetXaxis()->SetLabelSize(0.15);
	framex2->GetYaxis()->CenterTitle();
	framex2->GetYaxis()->SetTitleOffset(0.45);
	framex2->GetXaxis()->SetTitleOffset(1.1);
	framex2->GetYaxis()->SetNdivisions(505);
	framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
	pad2->cd();
	framex2->Draw();

	cdata->cd();

	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS  Y = "<<hpull->GetRMS(2)<<endl;

	// The sPlot technique requires that we fix the parameters
	// of the model that are not yields after doing the fit.
	RooRealVar* sigma1 = ws->var("sigma1");
	RooRealVar* sigma2 = ws->var("sigma2");
	RooRealVar* mean   = ws->var("mean");
	RooRealVar* c0     = ws->var("c0");
	RooRealVar* c1     = ws->var("c1");

	mean->setConstant(kTRUE);
	sigma1->setConstant(kTRUE);
	sigma2->setConstant(kTRUE);
	c0->setConstant(kTRUE);
	c1->setConstant(kTRUE);

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

	new TCanvas();
	RooRealVar* sigYield_sw = ws->var("sigYield_sw");
	sigYield_sw->setRange(-1.0,1.5);
	RooPlot* frame_sigsw = sigYield_sw->frame();
	frame_sigsw->SetTitle("Signal sweights");
	data->plotOn(frame_sigsw);
	frame_sigsw->Draw();

	new TCanvas();
	RooPlot* frame_sigsw_x = Lb_DTF_M_JpsiLConstr->frame();
	frame_sigsw_x->SetTitle("sWeights vs bmass");
	data->plotOnXY(frame_sigsw_x,YVar(*sigYield_sw),MarkerColor(kBlue));
	frame_sigsw_x->Draw();

	// create weightfed data set
	RooDataSet *dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"sigYield_sw");
	RooDataSet *dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw");

	TCanvas *can2 = new TCanvas();
	can2->Divide(1,2);

	can2->cd(1);

	RooPlot* frame2_1 = Lb_DTF_M_JpsiLConstr->frame();
	dataw_sig->plotOn(frame2_1, DataError(RooAbsData::SumW2));
	frame2_1->SetTitle("Lb DTF Mass Distribution (SW)");
	frame2_1->Draw();

	can2->cd(2);

	RooPlot* frame2_2 = Lb_DTF_M_JpsiLConstr->frame();
	dataw_bkg->plotOn(frame2_2, DataError(RooAbsData::SumW2));
	frame2_2->SetTitle("Lb DTF Mass Distribution (BW)");
	frame2_2->Draw();

	// import this new dataset with sWeights
	cout << "import new dataset with sWeights" << endl;
	ws->import(*data, Rename("dataWithSWeights"));
	cout <<" after splot import data set looks like"<<endl;
	ws->Print();

	Int_t dsentries = data->numEntries();
	cout<<"No. of entries in dataset = "<<dsentries<<endl;

	cout<<"Writing output tree with unbinned sweighted data"<<endl;
	fileout = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw.root",run,type),"RECREATE");
	treeout = new TTree("MyTuple","SWeighted Stuff");

	Double_t B_MASS = 0., SIGW = 0., BKGW = 0.;
	TBranch *bMass(0),*sigW(0),*bkgW(0);
	bMass = treeout->Branch("B_MASS",&B_MASS,"B_MASS/D");
	sigW  = treeout->Branch("SIGW",&SIGW,"SIGW/D");
	bkgW  = treeout->Branch("BKGW",&BKGW,"BKGW/D");

	for(Int_t i = 0; i < dsentries; i++)
	{
		B_MASS = (data->get(i))->getRealValue("Lb_DTF_M_JpsiLConstr");
		SIGW = (data->get(i))->getRealValue("sigYield_sw");
		BKGW = (data->get(i))->getRealValue("bkgYield_sw");

		bMass->Fill();
		sigW->Fill();
		bkgW->Fill();
		//	treeout->Fill();
	}

	treeout->SetBranchStatus("*",0);
	treeout->SetBranchStatus("SW",1);
	treeout->SetBranchStatus("BW",1);
	treeout->SetBranchStatus("Added_H_PT",1);
	treeout->SetBranchStatus("psi_1S_H_MINIPCHI2",1);
	treeout->SetBranchStatus("psi_1S_H_VERTEXCHI2_NEW",1);
	treeout->SetBranchStatus("psi_1S_H_IPCHI2_NEW",1);
	treeout->SetBranchStatus("psi_1S_VERTEX_Z_NEW",1);
	treeout->SetBranchStatus("Added_H_TRACKCHI2",1);
	treeout->SetBranchStatus("Added_H_GHOST",1);
	treeout->SetBranchStatus("Lb_H_OPENING",1);
	treeout->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

	fileout_training = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_forIsoTraining.root",run,type),"RECREATE");
	cout<<"Copying tree with s-weights to make signal training file for isolation"<<endl;
	treeout_training = treeout->CopyTree("");

	treeout_training->SetAlias("PT","Added_H_PT");
	treeout_training->SetAlias("MINIPCHI2","Added_H_MINIPCHI2");
	treeout_training->SetAlias("VCHI2DOF","psi_1S_H_VERTEXCHI2_NEW");
	treeout_training->SetAlias("IPCHI2","psi_1S_H_IPCHI2_NEW");
	treeout_training->SetAlias("TRACKORIVX_Z","psi_1S_H_VERTEXCHI2_NEW");
	treeout_training->SetAlias("GHOSTPROB","Added_H_GHOST");
	treeout_training->SetAlias("TRACKCHI2DOF","Added_H_TRACKCHI2");

	treeout->SetBranchStatus("*",1);
	fileout->cd();
	treeout->Write("",TObject::kOverwrite);
	//fileout->Write();
	fileout->Close();
	cout<<"Done writing output tree with unbinned sweighted data"<<endl;

	fileout_training->cd();
	treeout_training->Write();
	fileout_training->Close();
	cout<<"Done writing signal training file for isolation"<<endl;

	cout<<"Finishing DoSPlot()"<<endl;
}

void MakePlots(RooWorkspace* ws)//this function shouldn't be called now. It's function is taken over by DoSPlot.
{

	// Here we make plots of the discriminating variable (Lb_DTF_M_JpsiLConstr) after the fit
	// and of the control variable (isolation) after unfolding with sPlot.
	cout<<"Starting MakePlots()"<<endl;

	// make our canvas
	TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 1200, 800);

	///////////
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
	pad1->SetBottomMargin(0);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.5);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	cdata->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();

	gPad->SetTopMargin(0.06);
	pad1->Update();
	////////////

	// get what we need out of the workspace
	RooAbsPdf* model = ws->pdf("model");
	RooAbsPdf* sig1 = ws->pdf("sig1");
	RooAbsPdf* sig2 = ws->pdf("sig2");
	RooAbsPdf* sig = ws->pdf("sig");
	RooAbsPdf* bkg = ws->pdf("bkg");
	RooRealVar* Lb_DTF_M_JpsiLConstr = ws->var("Lb_DTF_M_JpsiLConstr");

	ws->Print();
	// note, we get the dataset with sWeights
	RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

	// this shouldn't be necessary, need to fix something with workspace
	// do this to set parameters back to their fitted values.
	model->fitTo(*data);

	//plot Lb_DTF_M_JpsiLConstr for data with full model and individual componenets overlayed

	RooPlot* frame = Lb_DTF_M_JpsiLConstr->frame();
	data->plotOn(frame,Name("Hist"),DataError(RooAbsData::SumW2) );
	model->plotOn(frame,Name("curvetot") );
	model->plotOn(frame,Components(*sig1),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig2),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig),LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen));

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

	// Pull distribution
	RooPlot *framex2 = Lb_DTF_M_JpsiLConstr->frame();
	RooHist* hpull = frame->pullHist("Hist","curvetot");
	framex2->addPlotable(hpull,"P");
	hpull->SetLineColor(kBlack);
	hpull->SetMarkerColor(kBlack);
	framex2->SetTitle(0);
	framex2->GetYaxis()->SetTitle("Pull");
	framex2->GetYaxis()->SetTitleSize(0.15);
	framex2->GetYaxis()->SetLabelSize(0.15);
	framex2->GetXaxis()->SetTitleSize(0.2);
	framex2->GetXaxis()->SetLabelSize(0.15);
	framex2->GetYaxis()->CenterTitle();
	framex2->GetYaxis()->SetTitleOffset(0.45);
	framex2->GetXaxis()->SetTitleOffset(1.1);
	framex2->GetYaxis()->SetNdivisions(505);
	framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
	pad2->cd();
	framex2->Draw();

	cdata->cd();

	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;

	new TCanvas();
	RooRealVar* sigYield_sw = ws->var("sigYield_sw");
	sigYield_sw->setRange(-1.0,1.5);
	RooPlot* frame_sigsw = sigYield_sw->frame();
	frame_sigsw->SetTitle("Signal sweights");
	data->plotOn(frame_sigsw);
	frame_sigsw->Draw();

	new TCanvas();
	RooPlot* frame_sigsw_x = Lb_DTF_M_JpsiLConstr->frame();
	frame_sigsw_x->SetTitle("sWeights vs bmass");
	data->plotOnXY(frame_sigsw_x,YVar(*sigYield_sw),MarkerColor(kBlue));
	frame_sigsw_x->Draw();

	// create weightfed data set
	RooDataSet *dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"sigYield_sw");
	RooDataSet *dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw");

	TCanvas *can2 = new TCanvas();
	can2->Divide(1,2);

	can2->cd(1);

	RooPlot* frame2_1= Lb_DTF_M_JpsiLConstr->frame();
	dataw_sig->plotOn(frame2_1, DataError(RooAbsData::SumW2));
	frame2_1->SetTitle("Lb DTF Mass Distribution (SW)");
	frame2_1->Draw();

	can2->cd(2);

	RooPlot* frame2_2= Lb_DTF_M_JpsiLConstr->frame();
	dataw_bkg->plotOn(frame2_2, DataError(RooAbsData::SumW2));
	frame2_2->SetTitle("Lb DTF Mass Distribution (BW)");
	frame2_2->Draw();

} //This shouldn't be called now. It's function is taken over by DoSPlot.
