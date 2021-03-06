/********************************
   Author : Aravindhan V.
   The purpose of this script is to sWeight data in 5200-6000 GeV, and write out weights.
   INPUT: data from Sanity.
   OUTPUT: weights for 5200-6000 GeV.
 *********************************/

#include "DoSWeight_Sanity.h"

void AddModel_Sanity(RooWorkspace *ws = nullptr, Int_t lowRange = 5200,
                     Int_t highRange = 6000, Int_t nEntries = 0, Int_t run = 1);

void AddData_Sanity(RooWorkspace *ws = nullptr, Int_t run = 1,
                    TTree *treeIn = nullptr);

Double_t DosPlot_Sanity(RooWorkspace *ws = nullptr, Int_t run = 1,
                        const char *type = "LL", TTree *treeOut = nullptr,
                        TH1D *hMass = nullptr);

Double_t DoSWeight_Sanity(Int_t run, Bool_t logFlag)
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

	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
	gROOT->ProcessLine(".x scripts/lhcbStyle.C");

	const char *type = "LL";

	if(logFlag)//Redirect output to log file
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/sPlot_Sanity_%s_log.txt",
		                             run, type),"w");
	}
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"==> Starting DoSWeight: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"********************************"<<endl;

	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;

	const char *rootFolder = "";
	TString inFileName     = "";
	TString outFileName    = "";
	TString trainFileName  = "";

	Int_t entries_init       = 0;
	Int_t entries_massWindow = 0;
	Int_t lowRange           = 5200;
	Int_t highRange          = 6000;//5200-6000 MeV window for sWeighting
	Int_t nBins              = (Int_t)(highRange-lowRange)/4;
	Double_t Lb_Mass         = 0.0;

	rootFolder    = Form("rootFiles/dataFiles/JpsiLambda/run%d", run);
	inFileName    = Form("%s/jpsilambda_sanity_%s.root", rootFolder, type);
	outFileName   = Form("%s/sWeightSanity/jpsilambda_%s_sanity_withsw.root",
	                     rootFolder, type);

	fileIn = TFile::Open(inFileName,"READ");
	if (!fileIn)
	{
		cout << "ERROR: could not open data file" << endl;
		exit(1);
	}
	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",
	                  nBins,lowRange,highRange),"","goff");
	TH1D *hMass = (TH1D*)gDirectory->Get("hMass");

	entries_init       = treeIn->GetEntries();
	entries_massWindow = treeIn->GetEntries(Form("Lb_DTF_M_JpsiLConstr > %d &&"
	                                             " Lb_DTF_M_JpsiLConstr < %d",
	                                             lowRange,highRange));

	treeIn->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&Lb_Mass);
	cout<<"Incoming Entries = "<<entries_init<<endl;

	fileOut = new TFile(outFileName,"RECREATE");
	treeOut = treeIn->CloneTree(0);

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"sWeighted Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	cout<<"Making mass cut on "<<entries_init<<" entries"<<endl;
	for (Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0) cout<<i<<endl;
		treeIn->GetEntry(i);
		if(Lb_Mass > lowRange && Lb_Mass < highRange)
		{
			treeOut->Fill();
		}
	}

	// Create a new workspace to manage the project.
	RooWorkspace* wSpace = new RooWorkspace("myWS");

	// add the signal and background models to the workspace.
	// Inside this function you will find a discription the model.
	AddModel_Sanity(wSpace, lowRange, highRange, entries_massWindow, run);

	// add data to the workspace
	AddData_Sanity(wSpace, run, treeIn);

	// do sPlot.
	//This wil make a new dataset with sWeights added for every event.
	Double_t myChi2 = DosPlot_Sanity(wSpace, run, type, treeOut, hMass);

	// inspect the workspace if you wish
	wSpace->Print();

	//save workspace in ROOT file so you don't have to make and load input data each time, dummy
	wSpace->writeToFile(Form("%s/wSpace_sPlot_Sanity.root",rootFolder));

	fileOut->cd();
	treeOut->Write("",TObject::kOverwrite);
	fileOut->Close();

	delete wSpace;

	sw.Stop();
	cout << "==> DoSWeight is done! Check fit status! sWeights FTW!: ";
	sw.Print();
	// if(logFlag) gROOT->ProcessLine(".>"); //end redirect
	if(logFlag) gSystem->RedirectOutput(0); //end redirect

	ofstream myfile;
	myfile.open("chi2_sanity.txt");
	myfile<<myChi2<<"\n";
	myfile.close();
	return myChi2;
}
void AddModel_Sanity(RooWorkspace *ws, Int_t lowRange, Int_t highRange,
                     Int_t nEntries, Int_t run)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting AddModel()"<<endl;
	cout<<"******************************************"<<endl;

	Int_t binSize = 4; //4 MeV bins.

	// Make models for signal and background

	// make a RooRealVar for the observables
	RooRealVar Lb_Mass("Lb_DTF_M_JpsiLConstr", "m_{J/#psi#Lambda}",
	                   lowRange, highRange,"MeV/#it{c}^{2}");//discriminating variable

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
	Int_t init_sigval = 0, ul_sigval = 0;
	if(run == 1)
	{
		init_sigval = 7000;
		ul_sigval   = 10000;
	}
	else if(run == 2)
	{
		init_sigval = 27000;
		ul_sigval   = 29000;
	}

	RooRealVar sigYield("sigYield","fitted yield for sig",init_sigval,0, ul_sigval);
	RooRealVar bkgYield("bkgYield","fitted yield for bkg",0, nEntries);


	RooAddPdf model("model","signal+background models",
	                RooArgList(sig, bkg),
	                RooArgList(sigYield,bkgYield));

	cout << "Importing model into workspace" << endl;

	ws->import(model);
	ws->Print();

	cout<<"******************************************"<<endl;
	cout<<"AddModel() complete"<<endl;
	cout<<"******************************************"<<endl;
}
void AddData_Sanity(RooWorkspace *ws, Int_t run, TTree *treeIn)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting AddData()"<<endl;
	cout<<"******************************************"<<endl;

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	RooRealVar Lb_Mass = *(ws->var("Lb_DTF_M_JpsiLConstr"));

	cout<<"Importing data into workspace. This could take a while. Sit tight"<<endl;

	RooArgSet s(Lb_Mass);
	RooDataSet data("data","dataset with x",treeIn,s);

	ws->import(data, Rename("data"));
	ws->Print();

	cout<<"******************************************"<<endl;
	cout<<"Finishing AddData()"<<endl;
	cout<<"******************************************"<<endl;
}
Double_t DosPlot_Sanity(RooWorkspace* ws, Int_t run, const char *type, TTree *treeOut,
                        TH1D *hMass)
{
	cout<<"******************************************"<<endl;
	cout<<"Starting DoSPlot()"<<endl;
	cout<<"******************************************"<<endl;

	Float_t SIGW = 0., BKGW = 0., bMASS = 0.;
	TBranch *sigW = nullptr;
	TBranch *bkgW = nullptr;

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

	ws->import(*rdh,Rename("data_binned"));
	// fit the model to the data. unbinned extended ML fit
	//	RooFitResult *r = model->fitTo(*data, Extended(), Strategy(2), Save(true));
	// Modifying this to try something along the lines of what steve does
	RooArgSet *myVars;
	model->fitTo(*rdh,Hesse(kTRUE),Strategy(2));//Fit to binned dataset

	myVars = model->getParameters(*rdh);
	myVars->Print("v");
	myVars->setAttribAll("Constant",kTRUE);
	sigYield->setConstant(kFALSE);
	bkgYield->setConstant(kFALSE);

	RooFitResult *r = model->fitTo(*rdh,Extended(),Hesse(kTRUE), Save(true));//Fit to binned dataset again with all params except yields fixed to fitted value

	cout<<"Starting MakePlots()"<<endl;

	// make our canvas
	// TCanvas *fitCanvas = new TCanvas("sPlot","sPlot demo", 1200, 800);
	//
	// TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	// TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
	// pad1->SetBottomMargin(0);
	// pad2->SetTopMargin(0);
	// pad2->SetBottomMargin(0.5);
	// pad2->SetBorderMode(0);
	// pad1->SetBorderMode(0);
	// fitCanvas->SetBorderMode(0);
	// pad2->Draw();
	// pad1->Draw();
	// pad1->cd();
	//
	// gPad->SetTopMargin(0.06);
	// pad1->Update();

	TCanvas *fitCanvas = new TCanvas("fitCanvas","",600,400);

	//plot Lb_Mass for data with full model and individual componenets overlayed

	RooPlot* frame = Lb_Mass->frame();
	frame->GetXaxis()->SetTitle("m_{J/#psi#Lambda}[MeV/#it{c}^{2}]");
	frame->GetYaxis()->SetTitle("Candidates/(4 MeV/#it{c}^{2})");
	data->plotOn(frame,Name("Hist"),DataError(RooAbsData::SumW2) );
	model->plotOn(frame,Name("curvetot") );
	model->plotOn(frame,Components(*sig1),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig2),LineStyle(kDotted),LineColor(kMagenta));
	model->plotOn(frame,Components(*sig),LineStyle(kDashed), LineColor(kRed));
	model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen));
	// model->paramOn(frame, Format("NEU"), Layout(0.65,0.85,0.9));
	// frame->SetTitle("Fit of model to discriminating variable");
	frame->Draw();

	////Fit chi2

	Double_t chiSquare1 = frame->chiSquare("curvetot","Hist");
	cout<<"chi square1 = "<<chiSquare1<<endl;
	RooArgSet *floatpar = model->getParameters(data);
	floatpar->Print();
	int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	Double_t chi2 = frame->chiSquare("curvetot","Hist",floatpars);
	cout<<"chi square2 = "<<chi2<<endl;

	// TPaveLabel *t1 = new TPaveLabel(0.1,0.8,0.3,0.88, Form("#chi^{2}/dof = %f", chi2),"NDC");
	// t1->Draw();
	// frame->addObject(t1);
	// pad1->Update();

	// Pull distribution
	// RooPlot *framex2 = Lb_Mass->frame();
	// RooHist* hpull = frame->pullHist("Hist","curvetot");
	// framex2->addPlotable(hpull,"P");
	// hpull->SetLineColor(kBlack);
	// hpull->SetMarkerColor(kBlack);
	// framex2->SetTitle(0);
	// framex2->GetYaxis()->SetTitle("Pull");
	// framex2->GetYaxis()->SetTitleSize(0.11);
	// framex2->GetYaxis()->SetLabelSize(0.15);
	// framex2->GetXaxis()->SetTitleSize(0.15);
	// framex2->GetXaxis()->SetLabelSize(0.15);
	// framex2->GetYaxis()->CenterTitle();
	// framex2->GetYaxis()->SetTitleOffset(0.4);
	// framex2->GetXaxis()->SetTitleOffset(0.7);
	// framex2->GetYaxis()->SetNdivisions(505);
	// framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
	// pad2->cd();
	// framex2->Draw();
	//
	// fitCanvas->cd();
	//
	// fitCanvas->Update();
	fitCanvas->SaveAs(Form("plots/fit_sPlot_sanity_run%d.pdf",run));

	// cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	// cout<<"Pull RMS  Y = "<<hpull->GetRMS(2)<<endl;

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

	//****************************************************************************
	// Now we use the SPlot class to add SWeights to our data set
	// based on our model and our yield variables
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, model,
	                                             RooArgList(*sigYield,*bkgYield) );

	//****************************************************************************
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

	Int_t dsentries = data->numEntries();
	cout<<"No. of entries in dataset = "<<dsentries<<endl;

	cout<<"************************************************"<<endl;
	cout<<"Writing output tree with unbinned sweighted data"<<endl;
	cout<<"************************************************"<<endl;

	//Add SW & BW branches to output tree
	sigW  = treeOut->Branch("SW",&SIGW,"SIGW/F");
	bkgW  = treeOut->Branch("BW",&BKGW,"BKGW/F");

	for(Int_t i = 0; i < dsentries; i++)
	{
		bMASS = (data->get(i))->getRealValue("Lb_DTF_M_JpsiLConstr");
		SIGW  = (data->get(i))->getRealValue("sigYield_sw");
		BKGW  = (data->get(i))->getRealValue("bkgYield_sw");

		sigW->Fill();
		bkgW->Fill();
	}

	cout<<"******************************************"<<endl;
	cout<<"Finishing DoSPlot()"<<endl;
	cout<<"******************************************"<<endl;
	return chi2;
}
