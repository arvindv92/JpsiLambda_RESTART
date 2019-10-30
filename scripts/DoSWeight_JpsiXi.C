/********************************
   Author : Aravindhan V.
   The purpose of this script is to sWeight Xib -> J/psi Xi data in 5500-6100 GeV, and write out weights.
   INPUT: data from Cuts_JpsiXi.
   OUTPUT: weights for 5500-6100 GeV.
 *********************************/
#include "DoSWeight_JpsiXi.h"

void AddModel(RooWorkspace *ws = nullptr, Int_t lowRange = 5500,
              Int_t highRange = 6100, Int_t nEntries = 0, Int_t run = 1);

void AddData(RooWorkspace *ws = nullptr, Int_t run = 1,
             TTree *treeIn = nullptr);

Double_t DosPlot(RooWorkspace *ws = nullptr, Int_t run = 1,
                 TTree *treeOut = nullptr, TH1D *hMass = nullptr);

Double_t DoSWeight_JpsiXi(Int_t run, Bool_t logFlag)
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	gROOT->ProcessLine(".x scripts/lhcbStyle.C");

	if(logFlag)//Redirect output to log file
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiXi/run%d/sPlot.txt",run),"w");
	}
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"==> Starting DoSWeight_JpsiXi: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"********************************"<<endl;

	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;

	const char *rootFolder = "";
	TString inFileName     = "";
	TString outFileName    = "";

	Int_t entries_init       = 0;
	Int_t entries_massWindow = 0;
	Int_t lowRange           = 5500;
	Int_t highRange          = 6100;//5500-6100 MeV window for sWeighting
	Int_t nBins              = (Int_t)(highRange-lowRange)/6;
	Double_t Xib_Mass         = 0.0;

	rootFolder  = Form("rootFiles/dataFiles/JpsiXi/run%d", run);
	inFileName  = Form("%s/jpsixi_cut_Run%d.root", rootFolder,run);
	outFileName = Form("%s/jpsixi_withsw.root", rootFolder);

	fileIn = TFile::Open(inFileName,"READ");
	if (!fileIn)
	{
		cout << "ERROR: could not open input file" << endl;
		exit(1);
	}
	treeIn = (TTree*)fileIn->Get("MyTuple");
	entries_init       = treeIn->GetEntries();
	entries_massWindow = treeIn->GetEntries(Form("Xib_DTF_M_JpsiXiLConstr > %d &&"
	                                             " Xib_DTF_M_JpsiXiLConstr < %d",
	                                             lowRange,highRange));

	treeIn->SetBranchAddress("Xib_DTF_M_JpsiXiLConstr",&Xib_Mass);
	cout<<"Incoming Entries = "<<entries_init<<endl;

	treeIn->Draw(Form("Xib_DTF_M_JpsiXiLConstr>>hMass(%d,%d,%d)",
	                  nBins,lowRange,highRange),"","goff");
	TH1D *hMass = (TH1D*)gDirectory->Get("hMass");

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
		if(Xib_Mass > lowRange && Xib_Mass < highRange)
		{
			treeOut->Fill();
		}
	}

	// Create a new workspace to manage the project.
	RooWorkspace *wSpace = new RooWorkspace("myWS");

	// add the signal and background models to the workspace.
	// Inside this function you will find a discription the model.
	AddModel(wSpace, lowRange, highRange, entries_massWindow, run);

	// add data to the workspace
	AddData(wSpace, run, treeIn);

	// do sPlot.
	//This wil make a new dataset with sWeights added for every event.
	Double_t myChi2 = DosPlot(wSpace, run, treeOut, hMass);

	// inspect the workspace if you wish
	wSpace->Print();

	//save workspace in ROOT file so you don't have to make and load input data each time, dummy
	wSpace->writeToFile(Form("rootFiles/dataFiles/JpsiXi/run%d/"
	                         "wSpace_sPlot.root",run));

	fileOut->cd();
	treeOut->Write("",TObject::kOverwrite);
	fileOut->Close();

	// cleanup
	delete wSpace;

	sw.Stop();
	cout << "==> DoSWeight_JpsiXi is done! Check fit status! sWeights FTW!: ";
	sw.Print();

	if(logFlag) gSystem->RedirectOutput(0); //end redirect

	return myChi2;
}

void AddModel(RooWorkspace *ws, Int_t lowRange, Int_t highRange, Int_t nEntries, Int_t run)
{
	cout<<"Starting AddModel()"<<endl;

	Int_t binSize = 6; //6 MeV bins.

	// Make models for signal and background

	// make a RooRealVar for the observables
	RooRealVar Xib_Mass("Xib_DTF_M_JpsiXiLConstr", "m_{J/#psi#Xi}",
	                    lowRange, highRange,"MeV/#it{c}^{2}");//discriminating variable

	Int_t nbins = (Int_t)(highRange-lowRange)/binSize;
	Xib_Mass.setBins(nbins);

	// FIT MODEL = Double Gaussian + Exponential
	//SIGNAL MODEL
	cout << "Making signal model" << endl;

	RooRealVar mean("mean","Crystal Ball Mean",5795.2,5790.0,5800.0,"MeV");
	RooRealVar sigma("sigma","Crystal Ball sigma",10.0,0.0,20.0,"MeV");

	RooRealVar alpha1("alpha1","alpha1",1.068,0.0,1.5);
	RooRealVar alpha2("alpha2","alpha2",-1.102,-1.5,0.0);

	RooRealVar CBn("CBn","Crystal Ball n",10.0);//fixed in both sim fit and data fit.

	RooRealVar fs("fs","fs",1.00,0.5,2.0);//scaling factor for MC width
	RooFormulaVar sigma1("sigma1","","sigma*fs",RooArgList(sigma,fs));

	RooCBShape sig1("sig1","Crystal Ball 1",Xib_Mass,mean,sigma1,alpha1,CBn);
	RooCBShape sig2("sig2","Crystal Ball 2",Xib_Mass,mean,sigma1,alpha2,CBn);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5);//fixed in both sim fit and data fit.
	RooAddPdf sig("sig","Total CB signal",RooArgList(sig1,sig2),frac1);

	// BACKGROUND MODEL
	cout << "Making background model" << endl;

	RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
	RooExponential bkg("bkg","Exponential bkg",Xib_Mass,tau);

	// COMBINED MODEL	cout << "Making full model" << endl;

	Int_t init_sigval = 0, ul_sigval = 0;
	if(run == 1)
	{
		init_sigval = 30;
		ul_sigval   = 80;
	}
	else if(run == 2)
	{
		init_sigval = 160;
		ul_sigval   = 400;
	}

	RooRealVar sigYield("sigYield","fitted yield for sig",init_sigval,0, ul_sigval);
	RooRealVar bkgYield("bkgYield","fitted yield for bkg",0, nEntries);

	RooAddPdf model("model","signal+background models",
	                RooArgList(sig, bkg), RooArgList(sigYield, bkgYield));

	cout << "Importing model into workspace" << endl;

	ws->import(model);
	ws->Print();

	cout<<"AddModel() complete"<<endl;
}

void AddData(RooWorkspace *ws, Int_t run, TTree *treeIn)
{
	cout<<"Starting AddData()"<<endl;

	// get what we need out of the workspace
	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Xib_DTF_M_JpsiXiLConstr",1);
	RooRealVar Xib_Mass = *(ws->var("Xib_DTF_M_JpsiXiLConstr"));

	//Import the data into unbinned dataset
	cout<<"Importing data. This could take a while. Sit tight"<<endl;

	RooArgSet s(Xib_Mass);
	RooDataSet data("data","dataset with x",treeIn,s);

	ws->import(data, Rename("data"));
	ws->Print();
	cout<<"Finishing AddData()"<<endl;
}

Double_t DosPlot(RooWorkspace* ws, Int_t run, TTree *treeOut, TH1D *hMass)
{
	Float_t SIGW = 0., BKGW = 0., bMASS = 0.;
	TBranch *sigW = nullptr;//, *sigW_training = nullptr;
	TBranch *bkgW = nullptr;//, *bkgW_training = nullptr;

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
	RooRealVar *Xib_Mass = ws->var("Xib_DTF_M_JpsiXiLConstr");
	RooDataSet *data     = (RooDataSet*) ws->data("data");

	//Make binned dataset
	RooDataHist *rdh = new RooDataHist("rdh","",*Xib_Mass,hMass);

	ws->import(*rdh,Rename("data_binned"));
	// fit the model to the data. unbinned extended ML fit
	// RooFitResult *r = model->fitTo(*data, Extended(), Strategy(2), Save(true));
	// Modifying this to try something along the lines of what steve does
	RooArgSet *myVars;
	model->fitTo(*rdh,Hesse(kTRUE),Strategy(2));//Fit to binned dataset

	myVars = model->getParameters(*rdh);
	myVars->Print("v");
	myVars->setAttribAll("Constant",kTRUE);//set all params except the yields to be constant
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

	//plot Xib_Mass for data with full model and individual componenets overlayed

	RooPlot* frame = Xib_Mass->frame();
	frame->GetXaxis()->SetTitle("m_{J/#psi#Xi}[MeV/#it{c}^{2}]");
	frame->GetYaxis()->SetTitle("Candidates/(6 MeV/#it{c}^{2})");

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
	// RooPlot *framex2 = Xib_Mass->frame();
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

	//****************************************************************************
	// Now we use the SPlot class to add SWeights to our data set
	// based on our model and our yield variables
	RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot", *data, model,
	                                             RooArgList(*sigYield,*bkgYield) );

	//****************************************************************************
	data->Print("v");
	// Check that our weights have the desired properties

	// import this new dataset with sWeights
	cout << "import new dataset with sWeights" << endl;
	ws->import(*data, Rename("dataWithSWeights"));
	cout <<" after splot import data set looks like"<<endl;
	ws->Print();

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
		bMASS = (data->get(i))->getRealValue("Xib_DTF_M_JpsiXiLConstr");
		SIGW  = (data->get(i))->getRealValue("sigYield_sw");
		BKGW  = (data->get(i))->getRealValue("bkgYield_sw");

		sigW->Fill();
		bkgW->Fill();
	}

	cout<<"Finishing DoSPlot_JpsiXi()"<<endl;

	return chi2;
}
