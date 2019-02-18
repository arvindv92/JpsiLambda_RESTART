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

using namespace RooFit;
using namespace std;
void Fit_JpsiXi(Int_t run = 1, Bool_t isData = true, Bool_t logFlag = true)
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiXi/run%d/Fit_JpsiXi_log.txt",run),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput((Form("logs/mc/JpsiXi/run%d/Fit_JpsiXi_log.txt",run)),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"==> Starting Fit_JpsiXi: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr;
	TTree *treeIn = nullptr;

	Float_t mean_mc = 0., sigma_mc = 0.,  alpha1_mc = 0., alpha2_mc = 0.;

	if(run == 1)
	{
		mean_mc = 5.79559e+03;
//		mean_mc = 5796.70;
		sigma_mc = 6.10233e+00;
		alpha1_mc = 7.94704e-01;
		alpha2_mc = -8.76836e-01;
	}
	else if(run == 2)
	{
		mean_mc = 5.79515e+03;
//		mean_mc = 5796.70;
		sigma_mc = 7.30515e+00;
		alpha1_mc = 1.03411e+00;
		alpha2_mc = -1.15930e+00;
	}
	RooRealVar Xib_DTF_M_JpsiXiLConstr("Xib_DTF_M_JpsiXiLConstr","Xib_DTF_M_JpsiXiLConstr",5500.,6100.,"MeV");

	// RooRealVar mean("mean","Gaussian Mean",5795.0,5790.0,5800.0);
	// RooRealVar sigma1("sigma1","Gaussian sigma1",10.0,1.0,50.0);
	// RooRealVar sigma2("sigma2","Gaussian sigma2",10.0,1.0,50.0);
	// RooRealVar shift("shift","shift",2.0,0.0,5.0);

	// //  RooFormulaVar sigma_shifted("sigma_shifted","sigma + shift", RooArgSet(sigma,shift));

	// RooGaussian sig1("sig1","Gaussian signal1",Xib_DTF_M_JpsiXiLConstr,mean,sigma1);
	// RooGaussian sig2("sig2","Gaussian signal2",Xib_DTF_M_JpsiXiLConstr,mean,sigma2);

	// RooRealVar frac1("frac1","Fraction of sig1 in signal",0.3,0.,0.9);

	// RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);

	RooRealVar mean("mean","Crystal Ball Mean",5795.2,5790.0,5800.0,"MeV");//mean,alpha1,alpha2,sigma are fixed here by fitting to UNREWEIGHTED simulation
	RooRealVar alpha1("alpha1","alpha1",1.068,0.0,1.5);
	RooRealVar alpha2("alpha2","alpha2",-1.102,-1.5,0.0);
	RooRealVar sigma("sigma","Crystal Ball sigma",10.0,0.0,20.0,"MeV");
	RooRealVar CBn("CBn","Crystal Ball n",10.0);//fixed in both sim fit and here

	RooRealVar fs("fs","fs",1.00,0.5,2.0);
	if(!isData)
	{
		fs.setVal(1.0);
		fs.setConstant(kTRUE);
	}
	else
	{
		mean.setVal(mean_mc); mean.setConstant(kTRUE);
		sigma.setVal(sigma_mc); sigma.setConstant(kTRUE);
		alpha1.setVal(alpha1_mc); alpha1.setConstant(kTRUE);
		alpha2.setVal(alpha2_mc); alpha2.setConstant(kTRUE);
	}
	RooFormulaVar sigma1("sigma1","","sigma*fs",RooArgList(sigma,fs));

	RooCBShape sig1("sig1","Crystal Ball 1",Xib_DTF_M_JpsiXiLConstr,mean,sigma1,alpha1,CBn);
	RooCBShape sig2("sig2","Crystal Ball 2",Xib_DTF_M_JpsiXiLConstr,mean,sigma1,alpha2,CBn);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5);//fixed in both sim fit and here
	RooAddPdf sig("sig","Total CB signal",RooArgList(sig1,sig2),frac1);

	// RooRealVar c0("c0","c0",-1.,-50.,50.);
	// RooRealVar c1("c1","c2",0.1,-50.,50.);
	// RooRealVar c2("c2","c2",0.1,-50.,50.);

	RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
	RooExponential bkg("bkg","Exponential bkg",Xib_DTF_M_JpsiXiLConstr,tau);

	//  RooChebychev bkg("bkg","bkg",Xib_DTF_M_JpsiXiLConstr,RooArgList(c0,c1));
	cout<<"poop2"<<endl;
	//  TFile *fileIn = TFile::Open("../../stage1/jpsilambda_LL.root","READ");

	if(isData)
	{
		fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiXi/run%d/jpsixi_cut_Run%d_BachPiLL.root",run,run),"READ");
	}
	else
	{
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cut_LL_BachPiLL.root",run),"READ");
	}

	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Xib_DTF_M_JpsiXiLConstr",1);

	Int_t nentries = treeIn->GetEntries();
	cout<<"nentries = "<<nentries<<endl;

	RooDataSet ds("ds","ds",Xib_DTF_M_JpsiXiLConstr,Import(*treeIn));
	ds.Print();

	cout<<"stage1"<<endl;
	RooRealVar nsig("nsig","nsig",0,nentries);
	RooRealVar nbkg("nbkg","nbkg",0,nentries);

	RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
	cout<<"stage2"<<endl;
	model.fitTo(ds,Extended(),RooFit::Strategy(2));
	cout<<"stage3"<<endl;
	RooPlot *frame = Xib_DTF_M_JpsiXiLConstr.frame();
	ds.plotOn(frame,Name("data"));
	model.plotOn(frame,Name("fit"));
	model.paramOn(frame);
	model.plotOn(frame,Components(sig),LineStyle(kDashed));
	//  model.plotOn(frame,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	//model.plotOn(frame,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	model.plotOn(frame,Components(bkg),LineColor(kRed));



	Double_t chiSquare1 = frame->chiSquare("fit","data");
	cout<<"chi square1 = "<<chiSquare1<<endl;
	RooArgSet *floatpar = model.getParameters(ds);
	floatpar->Print();
	int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	cout<<"float pars = "<<floatpars<<endl;
	Double_t chi2 = frame->chiSquare("fit","data",floatpars);
	cout<<"chi square2 = "<<chi2<<endl;

	TCanvas* c = new TCanvas("Ksfit","Ksfit", 1200, 800);
	cout<<"poop1"<<endl;
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
	cout<<"poop3"<<endl;
	c->Modified();
	// Pull distribution
	RooPlot *framex2 = Xib_DTF_M_JpsiXiLConstr.frame();
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

	if(!isData)
	{
		c->SaveAs(Form("plots/mc/JpsiXi/run%d/fit.pdf",run));
	}
	else
	{
		c->SaveAs(Form("plots/data/JpsiXi/run%d/fit.pdf",run));
	}


	//Check pull mean and RMS. Ideally the pull mean should be 0 and the RMS should be 1
	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;
	nsig.Print();
	nbkg.Print();

	sw.Stop();
	cout << "==> Fit_JpsiXi is done! Well Fitted!: "; sw.Print();

	if(logFlag) gSystem->RedirectOutput(0);
}
