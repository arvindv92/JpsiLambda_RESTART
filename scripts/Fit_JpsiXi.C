#include "TFile.h"
#include "TTree.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TPaveText.h"
#include "RooHist.h"
#include "TAxis.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TROOT.h"
#include "RooMCStudy.h"

using namespace RooFit;
using namespace std;

void Fit_JpsiXi(Int_t run = 1, Bool_t isData = true, Bool_t logFlag = true)
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	gROOT->ProcessLine(".x scripts/lhcbStyle.C");

	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiXi/run%d/Fit_JpsiXi_log.txt",run),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput((Form("logs/mc/JpsiXi/run%d/Fit_JpsiXi_log.txt",run)),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting Fit_JpsiXi: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr;
	TTree *treeIn = nullptr;

	if(isData)
	{
		fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiXi/run%d/jpsixi_cut_Run%d.root",run,run),"READ");
	}
	else
	{
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cut_LL.root",run),"READ");
	}

	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Xib_DTF_M_JpsiXiLConstr",1);

	Int_t nEntries = treeIn->GetEntries();
	cout<<"nEntries = "<<nEntries<<endl;

	Float_t mean_mc = 0., sigma_mc = 0.,  alpha1_mc = 0., alpha2_mc = 0.;

	// ***** I dont like this *****
	if(run == 1)
	{
		mean_mc = 5.79559e+03;
		sigma_mc = 6.10233e+00;
		alpha1_mc = 7.94704e-01;
		alpha2_mc = -8.76836e-01;
	}
	else if(run == 2)
	{
		mean_mc = 5.79515e+03;
		sigma_mc = 7.30515e+00;
		alpha1_mc = 1.03411e+00;
		alpha2_mc = -1.15930e+00;
	}
	//***** *****

	RooRealVar Xib_M("Xib_DTF_M_JpsiXiLConstr","Xib_DTF_M_JpsiXiLConstr",5500.,6100.,"MeV");

	//NOMINAL FIT MODEL - Double Crystal Ball (common mean and sigma) + Exp Background
	RooRealVar mean("mean","Crystal Ball Mean",5795.2,5790.0,5800.0,"MeV");//mean,alpha1,alpha2,sigma are fixed here by fitting to UNREWEIGHTED simulation
	RooRealVar sigma("sigma","Crystal Ball sigma",10.0,0.0,20.0,"MeV");

	RooRealVar alpha1("alpha1","alpha1",1.068,0.0,1.5);
	RooRealVar alpha2("alpha2","alpha2",-1.102,-1.5,0.0);

	RooRealVar CBn("CBn","Crystal Ball n",10.0);//fixed in both sim fit and here

	RooRealVar fs("fs","fs",1.00,0.5,2.0);
	RooFormulaVar sigma1("sigma1","","sigma*fs",RooArgList(sigma,fs));

	RooCBShape sig1("sig1","Crystal Ball 1",Xib_M,mean,sigma1,alpha1,CBn);
	RooCBShape sig2("sig2","Crystal Ball 2",Xib_M,mean,sigma1,alpha2,CBn);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5);//fixed in both sim fit and here
	RooAddPdf sig("sig","Total CB signal",RooArgList(sig1,sig2),frac1);

	RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
	RooExponential bkg("bkg","Exponential bkg",Xib_M,tau);

	RooRealVar nSig("nSig","nSig",5,nEntries);
	RooRealVar nBkg("nBkg","nBkg",0,nEntries);

	RooAddPdf model_nom("model_nom","Nominal Fit Model",RooArgList(sig,bkg),RooArgList(nSig,nBkg));

//      //Alternate Fit Model 1 - DOUBLE GAUSSIAN (common mean) + EXPO
//      RooRealVar mean_gaus("mean_gaus","Gaussian Mean",5795.2,5790.0,5800.0,"MeV");
//      RooRealVar sigma1_gaus("sigma1_gaus","Gaussian sigma1",10.0,1.0,20.0,"MeV");
//      RooRealVar sigma2_gaus("sigma2_gaus","Gaussian sigma2",10.0,1.0,20.0,"MeV");
//
//      RooGaussian sig1_gaus("sig1_gaus","sig1",Xib_M,mean_gaus,sigma1_gaus);
//      RooGaussian sig2_gaus("sig2_gaus","sig2",Xib_M,mean_gaus,sigma2_gaus);
//
//      RooRealVar frac1_gaus("frac1_gaus","Fraction of sig1 in signal",0.5);
//      RooAddPdf sig_gaus("sig_gaus","Total Gaussian signal",RooArgList(sig1_gaus,sig2_gaus),frac1_gaus);
//
//      RooRealVar tau1("tau1","tau1",-0.0007,-0.01,-0.0000001);
//      RooExponential bkg1("bkg1","Exponential bkg1",Xib_M,tau1);
//
//      RooRealVar nSig1("nSig1","nSig1",5,nEntries);
//      RooRealVar nBkg1("nBkg1","nBkg1",0,nEntries);
//
//      RooAddPdf model_alt1("model_alt1","Alternate Fit Model 1",RooArgList(sig_gaus,bkg1),RooArgList(nSig1,nBkg1));
//
// //Alternate Fit Model 2 - DOUBLE GAUSSIAN (common mean) + 1st ORDER CHEBY
//
//      RooRealVar c0("c0","c0",0.,-1.,1.);
//      RooChebychev bkg2("bkg2","Chebychev bkg",Xib_M,RooArgSet(c0));
//
//      RooRealVar nBkg2("nBkg2","nBkg2",0,nEntries);
//      RooAddPdf model_alt2("model_alt2","Alternate Fit Model 2",RooArgList(sig_gaus,bkg2),RooArgList(nSig1,nBkg2));
//
// //Alternate Fit Model 3 - DOUBLE CB (common mean and Sigma) + 1st ORDER CHEBY
//      RooRealVar mean_CB("mean_CB","Crystal Ball Mean",5795.2,5790.0,5800.0,"MeV");
//      RooRealVar alpha1_CB("alpha1_CB","alpha1",1.068,0.2,10.0);
//      RooRealVar alpha2_CB("alpha2_CB","alpha2",-1.102,-10.0,-1.0);
//      RooRealVar sigma_CB("sigma_CB","Gaussian sigma1",10.0,1.0,15.0,"MeV");
//      RooRealVar CBn_CB("CBn","Crystal Ball n",10.0);
//
//      RooCBShape sig1_CB("sig1_CB","Crystal Ball 1",Xib_M,mean_CB,sigma_CB,alpha1_CB,CBn_CB);
//      RooCBShape sig2_CB("sig2_CB","Crystal Ball 2",Xib_M,mean_CB,sigma_CB,alpha2_CB,CBn_CB);
//
//      RooRealVar frac1_CB("frac1_CB","Fraction of sig1 in signal",0.5);
//      RooAddPdf sig_CB("sig_CB","Total CB signal",RooArgList(sig1_CB,sig2_CB),frac1_CB);
//
//      RooRealVar nSig2("nSig2","nSig2",5,nEntries);
//      RooAddPdf model_alt3("model_alt3","Alternate Fit Model 2",RooArgList(sig_CB,bkg2),RooArgList(nSig2,nBkg2));
//      //
//      // RooRealVar c0("c0","c0",0.,-1.,1.);
//      // RooChebychev bkg("bkg","Chebychev bkg",Xib_M,RooArgSet(c0));

	// Import data into unbinned dataset
	RooDataSet ds("ds","ds",Xib_M,Import(*treeIn));
	ds.Print();

	if(!isData)
	{
		fs.setVal(1.0);
		fs.setConstant(kTRUE);
	}
	// else
	// {
	//      mean.setVal(mean_mc); mean.setConstant(kTRUE);
	//      sigma.setVal(sigma_mc); sigma.setConstant(kTRUE);
	//      alpha1.setVal(alpha1_mc); alpha1.setConstant(kTRUE);
	//      alpha2.setVal(alpha2_mc); alpha2.setConstant(kTRUE);
	//
	//      // mean_CB.setVal(mean_mc); mean_CB.setConstant(kTRUE);
	//      // sigma_CB.setVal(sigma_mc); sigma_CB.setConstant(kTRUE);
	//      // alpha1_CB.setVal(alpha1_mc); alpha1_CB.setConstant(kTRUE);
	//      // alpha2_CB.setVal(alpha2_mc); alpha2_CB.setConstant(kTRUE);
	//
	//      if(run == 1)
	//      {
	//              // mean_gaus.setVal(5795.44); mean_gaus.setConstant(kTRUE);
	//              // sigma1_gaus.setVal(13.15); sigma1_gaus.setConstant(kTRUE);
	//              // sigma2_gaus.setVal(4.89); sigma2_gaus.setConstant(kTRUE);
	//      }
	//      else if(run == 2)
	//      {
	//              // mean_gaus.setVal(5795.06); mean_gaus.setConstant(kTRUE);
	//              // sigma1_gaus.setVal(5.58); sigma1_gaus.setConstant(kTRUE);
	//              // sigma2_gaus.setVal(12.29); sigma2_gaus.setConstant(kTRUE);
	//      }
	// }

	//Fit nominal fit model to dataset
	model_nom.fitTo(ds,Extended(),RooFit::Strategy(2));
	//Fit alternate fit model to dataset
	// model_alt1.fitTo(ds,Extended(),RooFit::Strategy(2));
	// model_alt2.fitTo(ds,Extended(),RooFit::Strategy(2));
	// model_alt3.fitTo(ds,Extended(),RooFit::Strategy(2));

	RooAddPdf myModel = model_nom;

	RooPlot *frame = Xib_M.frame();
	ds.plotOn(frame,Name("data"));
	myModel.plotOn(frame,Name("fit"));
	myModel.paramOn(frame);
	myModel.plotOn(frame,Components(sig),LineStyle(kDashed));
	myModel.plotOn(frame,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	myModel.plotOn(frame,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	myModel.plotOn(frame,Components(bkg),LineColor(kRed));

	Double_t chiSquare1 = frame->chiSquare("fit","data");
	cout<<"chi square1 = "<<chiSquare1<<endl;
	RooArgSet *floatpar = myModel.getParameters(ds);
	floatpar->Print();
	int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	cout<<"float pars = "<<floatpars<<endl;
	Double_t chi2 = frame->chiSquare("fit","data",floatpars);
	cout<<"chi square2 = "<<chi2<<endl;

	TCanvas* c = new TCanvas("c","#Xi_{b} #to J/#psi #Xi Fit", 1200, 800);
	// ///////////
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
	pad1->SetBottomMargin(0);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.5);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	c->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();

	gPad->SetTopMargin(0.06);
	pad1->Update();

	frame->Draw();

	//  TPaveText *pt = (TPaveText*)pad1->GetListOfPrimitives()->FindObject("model_paramBox");

	// pt->AddText(Form("#chi^{2}/dof = %f",chi2));
	c->Modified();

	// Pull distribution
	RooPlot *framex2 = Xib_M.frame();
	RooHist* hpull = frame->pullHist("data","fit");
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
	framex2->GetYaxis()->SetTitleOffset(0.25);
	framex2->GetXaxis()->SetTitleOffset(1.1);
	framex2->GetYaxis()->SetNdivisions(505);
	framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
	pad2->cd();
	framex2->Draw();

	c->cd();

	//Check pull mean and RMS. Ideally the pull mean should be 0 and the RMS should be 1
	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;

	if(!isData)
	{
		c->SaveAs(Form("plots/mc/JpsiXi/run%d/NominalFit.pdf",run));
	}
	else
	{
		c->SaveAs(Form("plots/data/JpsiXi/run%d/NominalFit.pdf",run));
	}

	nSig.Print();
	nBkg.Print();

	// RooMCStudy* mcs = new RooMCStudy(model_alt3,Xib_M,FitModel(model_nom),Silence(),Binned(kTRUE),Extended(),FitOptions(Save(kTRUE),PrintEvalErrors(0)));
	// mcs->generateAndFit(5000);//Is it actually fluctuating the generated yield within the errors?
	//
	// // RooFormulaVar par("par","par","(nSig-nSig1)/nSig1",RooArgList(nSig,nSig1));
	// RooPlot* frame_mcs = mcs->plotParam(nSig,Bins(40));
	// new TCanvas();
	// frame_mcs->Draw();

	sw.Stop();
	cout << "==> Fit_JpsiXi is done! Well Fitted!: "; sw.Print();

	if(logFlag) gSystem->RedirectOutput(0);
}
