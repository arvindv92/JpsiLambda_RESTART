/********************************
   Author : Aravindhan V.
   The purpose of this script is to sWeight data in 5200-6000 GeV, and write out weights.
   When running with isoFlag = true, zeroTracks and nonZeroTracks data are sWeighted separately.
   What about running with isoFlag = false? sWeight them separately anyway?
   INPUT: data from CutOutKs.
   OUTPUT: weights for 5200-6000 GeV, isolation training file.
 *********************************/

#include "DoSWeight_Sanity.h"
// use this order for safety on library loading
using namespace RooFit;
using namespace RooStats;
using namespace std;

void AddModel(RooWorkspace *ws = nullptr, Int_t lowRange = 5200,
              Int_t highRange = 6000, Int_t nEntries = 0, Int_t run = 1);

void AddData(RooWorkspace *ws = nullptr, Int_t run = 1,
             TTree *treeIn = nullptr);

Double_t DosPlot(RooWorkspace *ws = nullptr, Int_t run = 1,
                 const char *type = "LL", TTree *treeOut = nullptr,
                 TH1D *hMass = nullptr);

Double_t DoSWeight_Sanity(Int_t run, Int_t trackType, Bool_t logFlag)
/*
   >run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC.
   Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data.

   >trackType = 3 for LL, 5 for DD.

   >zeroFlag = true when running over zeroTracks data coming out of cutOutKs.
   false when running over nonZeroTracks data.
 */
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	const char *type       = "";

	// suffix        = (zeroFlag) ? ("_ZeroTracks") : ("_nonZeroTracks");
	type          = (trackType == 3) ? ("LL")    : ("DD");

	if(logFlag)//Redirect output to log file
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/sPlot_Sanity_%s_noPID_log.txt",
		                             run, type),"w");
	}
	cout<<"********************************"<<endl;
	cout<<"==> Starting DoSWeight: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"********************************"<<endl;

	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;


	const char *rootFolder = "";
	TString inFileName      = "", outFileName = "", trainFileName = "";

	Int_t entries_init = 0, entries_massWindow = 0;
	Int_t lowRange     = 5200, highRange       = 6000;//5200-6000 MeV window for sWeighting
	Int_t nBins = (Int_t)(highRange-lowRange)/4;
	Double_t Lb_Mass   = 0.0;

	rootFolder    = Form("rootFiles/dataFiles/JpsiLambda/run%d", run);

	inFileName    = Form("%s/jpsilambda_sanity_%s_noPID.root", rootFolder, type);
	outFileName   = Form("%s/sWeightSanity/jpsilambda_%s_sanity_withsw_noPID.root", rootFolder, type);
	// trainFileName = Form("%s/jpsilambda_%s_forIsoTraining.root", rootFolder, type);

	fileIn = TFile::Open(inFileName,"READ");
	if (!fileIn)
	{
		cout << "ERROR: could not open data file" << endl;
		exit(1);
	}
	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nBins,lowRange,highRange),"","goff");
	TH1D *hMass = (TH1D*)gDirectory->Get("hMass");

	entries_init       = treeIn->GetEntries();
	entries_massWindow = treeIn->GetEntries(Form("Lb_DTF_M_JpsiLConstr > %d &&"
	                                             " Lb_DTF_M_JpsiLConstr < %d",
	                                             lowRange,highRange));

	fileOut = new TFile(outFileName,"RECREATE");
	treeOut = treeIn->CloneTree(0);

	// if(!zeroFlag)
	// {
	//      fileOut_training = new TFile(trainFileName,"RECREATE");
	//      treeOut_training = treeIn->CloneTree(0);
	// }

	treeIn->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&Lb_Mass);
	cout<<"Incoming Entries = "<<entries_init<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"sWeighted Output file = "<<fileOut->GetName()<<endl;
	// if(!zeroFlag)
	// {
	//      cout<<"Isolation Training Output file = "
	//          <<fileOut_training->GetName()<<endl;
	// }
	cout<<"******************************************"<<endl;

	cout<<"Making mass cut on "<<entries_init<<" entries"<<endl;
	for (Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0) cout<<i<<endl;
		treeIn->GetEntry(i);
		if(Lb_Mass > lowRange && Lb_Mass < highRange)
		{
			treeOut->Fill();
			// if(!zeroFlag)
			// {
			//      // if(Lb_Mass > 5400 && Lb_Mass < 5700)//5400-5700 only for training
			//      //      {
			//      // treeOut_training->Fill();
			//      //	}
			// }
		}
	}

	// Create a new workspace to manage the project.
	RooWorkspace* wSpace = new RooWorkspace("myWS");

	// add the signal and background models to the workspace.
	// Inside this function you will find a discription the model.
	AddModel(wSpace, lowRange, highRange, entries_massWindow, run);

	// add data to the workspace
	AddData(wSpace, run, treeIn);

	// inspect the workspace if you wish
	wSpace->Print();

	// do sPlot.
	//This wil make a new dataset with sWeights added for every event.
	Double_t myChi2 = DosPlot(wSpace, run, type, treeOut, hMass);

	fileOut->cd();
	treeOut->Write("",TObject::kOverwrite);
	fileOut->Close();

	// if(!zeroFlag)
	// {
	//      fileOut_training->cd();
	//      treeOut_training->Write("",TObject::kOverwrite);
	//      fileOut_training->Close();
	// }
	// cleanup
	delete wSpace;

	sw.Stop();
	cout << "==> DoSWeight is done! Check fit status! sWeights FTW!: ";
	sw.Print();
	// if(logFlag) gROOT->ProcessLine(".>"); //end redirect
	if(logFlag) gSystem->RedirectOutput(0); //end redirect

	return myChi2;
}
void AddModel(RooWorkspace *ws, Int_t lowRange, Int_t highRange, Int_t nEntries, Int_t run)
{
	cout<<"Starting AddModel()"<<endl;

	Int_t binSize = 4; /*4 MeV bins.
	                      Bin size only affects visualization and chi square calculation,
	                      since we are doing an unbinned fit.*/

	// Make models for signal and background

	// make a RooRealVar for the observables
	RooRealVar Lb_Mass("Lb_DTF_M_JpsiLConstr", "M_{J/#psi#Lambda}",
	                   lowRange, highRange,"MeV");//discriminating variable

	Int_t nbins = (Int_t)(highRange-lowRange)/binSize;
	Lb_Mass.setBins(nbins);

	// SIGNAL MODEL
	cout << "Making signal model" << endl;
	RooRealVar mean("mean","Gaussian Mean",5619.42,5619.0,5621.0);
	RooRealVar sigma1("sigma1","Gaussian sigma1",7.0,5.0,20.0);
	RooRealVar sigma2("sigma2","Gaussian sigma2",7.0,3.0,20.0);

	RooGaussian sig1("sig1","Gaussian signal1",Lb_Mass,mean,sigma1);
	RooGaussian sig2("sig2","Gaussian signal2",Lb_Mass,mean,sigma2);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5,0.01,1.0);

	RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);

	// BACKGROUND MODEL
	cout << "Making background model" << endl;
	// RooRealVar c0("c0","c0",-0.1,-1.0,1.0);
	// RooRealVar c1("c1","c1",0.1,-1.0,1.0);
	// RooChebychev bkg("bkg","bkg",Lb_Mass,RooArgList(c0,c1));

	RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
	RooExponential bkg("bkg","Exponential Bkg",Lb_Mass,tau);

	// COMBINED MODEL	cout << "Making full model" << endl;
	Int_t init_sigval = 0;
	if(run == 1) init_sigval = 7000;
	else if(run == 2) init_sigval = 27000;

	RooRealVar sigYield("sigYield","fitted yield for sig",init_sigval,0, 35000);
	RooRealVar bkgYield("bkgYield","fitted yield for bkg",(nEntries/2),0, nEntries);

	RooAddPdf model("model","signal+background models",
	                RooArgList(sig, bkg),
	                RooArgList(sigYield,bkgYield));

	cout << "Importing model into workspace" << endl;

	ws->import(model);
	ws->Print();

	cout<<"AddModel() complete"<<endl;
}
void AddData(RooWorkspace *ws, Int_t run, TTree *treeIn)
{
	cout<<"Starting AddData()"<<endl;

	// get what we need out of the workspace to make toy data
	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	RooRealVar Lb_Mass = *(ws->var("Lb_DTF_M_JpsiLConstr"));

	//Import the data
	cout<<"Importing data. This could take a while. Sit tight"<<endl;

	RooArgSet s(Lb_Mass);
	RooDataSet data("data","dataset with x",treeIn,s);

	ws->import(data, Rename("data"));
	ws->Print();
	cout<<"Finishing AddData()"<<endl;
}
Double_t DosPlot(RooWorkspace* ws, Int_t run, const char *type, TTree *treeOut,
                 TH1D *hMass)
{
	// const char *suffix = (zeroFlag) ? ("_ZeroTracks") : ("_nonZeroTracks");
	Float_t SIGW = 0., BKGW = 0., bMASS = 0.;
	TBranch *sigW = nullptr;
	TBranch *bkgW = nullptr;

	cout<<"Starting DoSPlot()"<<endl;
	cout << "Calculating sWeights" << endl;

	// get what we need out of the workspace to do the fit
	RooAbsPdf *model     = ws->pdf("model");
	RooAbsPdf *sig1      = ws->pdf("sig1");
	RooAbsPdf *sig2      = ws->pdf("sig2");
	RooAbsPdf *sig       = ws->pdf("sig");
	RooAbsPdf *bkg       = ws->pdf("bkg");
	RooRealVar *sigYield = ws->var("sigYield");
	RooRealVar *bkgYield = ws->var("bkgYield");
	RooRealVar *Lb_Mass  = ws->var("Lb_DTF_M_JpsiLConstr");
	RooDataSet *data     = (RooDataSet*) ws->data("data");

	RooDataHist *rdh = new RooDataHist("rdh","",*Lb_Mass,hMass);

	// fit the model to the data. unbinned extended ML fit
	//	RooFitResult *r = model->fitTo(*data, Extended(), Strategy(2), Save(true));
	// Modifying this to try something along the lines of what steve does
	RooArgSet *myVars;
	model->fitTo(*rdh,Hesse(kTRUE),Strategy(2));
	myVars = model->getParameters(*rdh);
	myVars->Print("v");
	myVars->setAttribAll("Constant",kTRUE);
	sigYield->setConstant(kFALSE);
	bkgYield->setConstant(kFALSE);
	RooFitResult *r = model->fitTo(*rdh,Extended(),Hesse(kTRUE), Save(true));

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

	//plot Lb_Mass for data with full model and individual componenets overlayed

	RooPlot* frame = Lb_Mass->frame();
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
	RooPlot *framex2 = Lb_Mass->frame();
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
	fitCanvas->SaveAs(Form("plots/fit_sanity_run%d%s_noPID.pdf",run,type));

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
	// RooRealVar *sigma1 = ws->var("sigma1");
	// RooRealVar *sigma2 = ws->var("sigma2");
	// RooRealVar *mean   = ws->var("mean");
	// RooRealVar *tau    = ws->var("tau");
	//
	// mean->setConstant(kTRUE);
	// sigma1->setConstant(kTRUE);
	// sigma2->setConstant(kTRUE);
	// tau->setConstant(kTRUE);
	//
	// RooMsgService::instance().setSilentMode(true);

	// Now we use the SPlot class to add SWeights to our data set
	// based on our model and our yield variables
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, model,
	                                             RooArgList(*sigYield,*bkgYield) );

	data->Print("v");
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
	sWeightCanvas->SaveAs(Form("plots/fit_sanity_run%d%s_sWeights_noPID.pdf",run,type));

	// TCanvas *sWeightVsMass = new TCanvas();
	// RooPlot* frame_sigsw_x = Lb_Mass->frame();
	// frame_sigsw_x->SetTitle("sWeights vs bmass");
	// data->plotOnXY(frame_sigsw_x,YVar(*sigYield_sw),MarkerColor(kBlue));
	// frame_sigsw_x->Draw();
	// sWeightVsMass->Update();
	// sWeightVsMass->SaveAs(Form("plots/fit_sanity_run%d%s_sWeightsVsMass.pdf",run,type));

	// create weightfed data set
	RooDataSet *dataw_sig = new RooDataSet(data->GetName(),data->GetTitle(),
	                                       data,*data->get(),0,"sigYield_sw");
	RooDataSet *dataw_bkg = new RooDataSet(data->GetName(),data->GetTitle(),
	                                       data,*data->get(),0,"bkgYield_sw");

	TCanvas *sWeightMass = new TCanvas();
	sWeightMass->Divide(1,2);

	sWeightMass->cd(1);

	RooPlot* frame2_1 = Lb_Mass->frame();
	dataw_sig->plotOn(frame2_1, DataError(RooAbsData::SumW2));
	frame2_1->SetTitle("Lb DTF Mass Distribution (SW)");
	frame2_1->Draw();

	sWeightMass->cd(2);

	RooPlot* frame2_2 = Lb_Mass->frame();
	dataw_bkg->plotOn(frame2_2, DataError(RooAbsData::SumW2));
	frame2_2->SetTitle("Lb DTF Mass Distribution (BW)");
	frame2_2->Draw();

	sWeightMass->Update();
	sWeightMass->SaveAs(Form("plots/fit_run%d%s_sWeightedMass_noPID.pdf",run,type));

	Int_t dsentries = data->numEntries();
	cout<<"No. of entries in dataset = "<<dsentries<<endl;

	cout<<"************************************************"<<endl;
	cout<<"Writing output tree with unbinned sweighted data"<<endl;
	cout<<"************************************************"<<endl;

	//Add SW & BW branches to output tree
	sigW  = treeOut->Branch("SW",&SIGW,"SIGW/F");
	bkgW  = treeOut->Branch("BW",&BKGW,"BKGW/F");

	//Add SW & BW branches to isolation training tree
	// if(!zeroFlag)
	// {
	//      sigW_training = treeOut_training->Branch("SW",&SIGW,"SIGW/F");
	//      bkgW_training = treeOut_training->Branch("BW",&BKGW,"BKGW/F");
	// }
	for(Int_t i = 0; i < dsentries; i++)
	{
		bMASS = (data->get(i))->getRealValue("Lb_DTF_M_JpsiLConstr");
		SIGW  = (data->get(i))->getRealValue("sigYield_sw");
		BKGW  = (data->get(i))->getRealValue("bkgYield_sw");

		sigW->Fill();
		bkgW->Fill();
		// if(!zeroFlag)
		// {
		//      if (bMASS > 5400 && bMASS < 5700)
		//      {
		//              sigW_training->Fill();
		//              bkgW_training->Fill();
		//      }
		// }
	}

	//Activate only necessary branches in isolation training tree.
	//It only needs to contain training variables + weights
	// if(!zeroFlag)
	// {
	//      treeOut_training->SetBranchStatus("*",0);
	//      treeOut_training->SetBranchStatus("SW",1);
	//      treeOut_training->SetBranchStatus("BW",1);
	//      treeOut_training->SetBranchStatus("Added_H_PT",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_MINIPCHI2",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_VERTEXCHI2_NEW",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_IPCHI2_NEW",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_IP_NEW",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_FD_NEW",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_FDCHI2_NEW",1);
	//      treeOut_training->SetBranchStatus("psi_1S_H_VERTEX_Z_NEW",1);
	//      treeOut_training->SetBranchStatus("Added_H_TRACKCHI2",1);
	//      treeOut_training->SetBranchStatus("Added_H_GHOST",1);
	//      treeOut_training->SetBranchStatus("Lb_H_OPENING",1);
	//      treeOut_training->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	//
	//      //Give aliases to the branches in isolation training tree
	//      treeOut_training->SetAlias("PT","Added_H_PT");
	//      treeOut_training->SetAlias("MINIPCHI2","psi_1S_H_MINIPCHI2");
	//      treeOut_training->SetAlias("VCHI2DOF","psi_1S_H_VERTEXCHI2_NEW");
	//      treeOut_training->SetAlias("IPCHI2","psi_1S_H_IPCHI2_NEW");
	//      treeOut_training->SetAlias("IP","psi_1S_H_IP_NEW");
	//      treeOut_training->SetAlias("FD","psi_1S_H_FD_NEW");
	//      treeOut_training->SetAlias("FDCHI2","psi_1S_H_FDCHI2_NEW");
	//      treeOut_training->SetAlias("TRACKORIVX_Z","psi_1S_H_VERTEXCHI2_NEW");
	//      treeOut_training->SetAlias("GHOSTPROB","Added_H_GHOST");
	//      treeOut_training->SetAlias("TRACKCHI2DOF","Added_H_TRACKCHI2");
	// }
	cout<<"Finishing DoSPlot()"<<endl;

	return chi2;
}
