#include "TFile.h"
#include "TTree.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
//#include "RooChebychev.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooHistPdf.h"
#include "TCanvas.h"
#include "RooFitResult.h"
#include "TH1D.h"
#include "RooFFTConvPdf.h"
#include <iostream>
using namespace RooFit;
using namespace std;
void Fit_effs(Int_t mcType = 1, Int_t run = 1, TString stage = "triggered")
{
	Int_t binwidth = 4; //Final fit is a binned fit with this binning
	Int_t low = 5000, high = 6000;

	RooRealVar bMass("Lb_DTF_M_JpsiLConstr","Fit to J/#psi #Lambda inv. mass",low, high,"MeV"); //MASTER VARIABLE

	Int_t nbins = (Int_t)(high-low)/binwidth;
	bMass.setBins(nbins);

	TFile *fileIn = nullptr, *fileIn_gen = nullptr;
	TTree *treeIn = nullptr, *treeIn_gen = nullptr;
	RooFitResult *res = nullptr;

	cout<<"poop1"<<endl;
	if(mcType == 1)
	{
		RooRealVar mean("mean","Crystal Ball Mean",5619.6,5618,5622);
		RooRealVar alpha1("alpha1","alpha1",1.028, 0.8,1.3);
		RooRealVar alpha2("alpha2","alpha2",-1.097,-1.4,-0.7);
		//  RooRealVar sigma("sigma","Crystal Ball sigma",7.069);//,6.4,7.8);
		RooRealVar sigma("sigma","Crystal Ball sigma",0.,20.);//,6.4,7.8);
		RooRealVar CBn("CBn","Crystal Ball n",10.0);//fixed in both sim fit and here

		//RooGaussian mean_constraint("mean_constraint","CB mean constraint",mean,RooFit::RooConst(5620.09),RooFit::RooConst(0.107));
		//	RooGaussian meanpdg_constraint("mean_constraint","CB mean constraint",mean,RooFit::RooConst(5619.6),RooFit::RooConst(0.17));
		// RooGaussian alpha1_constraint("alpha1_constraint","CB alpha1 constraint",mean,RooFit::RooConst(1.028),RooFit::RooConst(0.041));
		// RooGaussian alpha2_constraint("alpha2_constraint","CB alpha2 constraint",mean,RooFit::RooConst(-1.097),RooFit::RooConst(0.046));
		// RooGaussian sigma_constraint("sigma_constraint","CB sigma constraint",mean,RooFit::RooConst(7.069),RooFit::RooConst(0.132));

		//	RooRealVar fs("fs","Scale factor for width",1.0,0.5,2.0);
		//	RooFormulaVar sigma1("sigma1","sigma1","sigma*fs",RooArgList(sigma,fs));

		RooCBShape sig1("sig1","Crystal Ball 1",bMass,mean,sigma,alpha1,CBn);//*****MODIFIED from sigma1 to sigma********
		RooCBShape sig2("sig2","Crystal Ball 2",bMass,mean,sigma,alpha2,CBn);//*****MODIFIED from sigma1 to sigma********

		RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5);//fixed in both sim fit and here
		RooAddPdf sig("sig","Total CB signal",RooArgList(sig1,sig2),frac1);
		//*********Exponential shape for continuum backgruond**************************************************
		RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
		RooExponential bkg("bkg","Exponential bkg",bMass,tau);

		fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_%s",run, stage.Data()),"READ");
		treeIn = (TTree*)fileIn->Get("MyTuple");

		Int_t nentries = treeIn->GetEntries();

		treeIn->SetBranchStatus("*",0);
		treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

		treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>h_lb(%d,%d,%d)",nbins,low,high),"","goff");
		TH1D *h_lb = (TH1D*)gDirectory->Get("h_lb");

		RooDataHist *dataHist = new RooDataHist("dataHist","dataHist",bMass,h_lb);

		dataHist->Print();

		RooRealVar nsig("nsig","nsig",1,nentries);
		RooRealVar nbkg("nbkg","nbkg",1,nentries);

		RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));

		res = model.fitTo(*dataHist,Extended(),RooFit::Save(true), Strategy(2));

		RooPlot *frame = new RooPlot(bMass, low,5800,100);
		//  frame->GetYaxis()->SetRangeUser(0,100);
		dataHist->plotOn(frame,Name("data"));
		model.plotOn(frame,Name("fit"));
		model.plotOn(frame,Components(sig),Name("sig"),LineColor(kMagenta+2));
		model.plotOn(frame,Components(bkg),LineColor(kRed),Name("bkg"));
		model.plotOn(frame,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
		model.plotOn(frame,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));

		new TCanvas();
		frame->Draw();

		Double_t chiSquare1 = frame->chiSquare("fit","data");
		cout<<"chi square1/dof = "<<chiSquare1<<endl;
		RooArgSet *floatpar = model.getParameters(dataHist);
		floatpar->Print();
		int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
		cout<<"float pars = "<<floatpars<<endl;
		Double_t chi2 = frame->chiSquare("fit","data",floatpars);
		cout<<"chi square2/dof = "<<chi2<<endl;
	}

	if(mcType == 2)
	{
		cout<<"poop2"<<endl;
		fileIn_gen = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma.root",run),"READ");
		treeIn_gen = (TTree*)fileIn_gen->Get("MCTuple/MCDecayTree");

		cout<<"poop3"<<endl;
		treeIn_gen->Draw(Form("sqrt(pow(J_psi_1S_TRUEP_E + Lambda0_TRUEP_E,2) - pow(J_psi_1S_TRUEP_X + Lambda0_TRUEP_X,2) - pow(J_psi_1S_TRUEP_Y + Lambda0_TRUEP_Y,2) - pow(J_psi_1S_TRUEP_Z + Lambda0_TRUEP_Z,2))>>h_sigma_gen(%d,%d,%d)",
		                      nbins,low,high),"","goff");

		cout<<"poop4"<<endl;
		TH1D *h_sigma_gen = (TH1D*)gDirectory->Get("h_sigma_gen");

		RooDataHist *genHist = new RooDataHist("genHist","genHist",bMass,h_sigma_gen);
		genHist->Print();

		RooHistPdf myHist("myHist","sigma_gen",bMass,*genHist,0);//Signal shape, hopefully

		RooRealVar mg("mg","mg",0);
		RooRealVar sg("sg","sg",2,0.1,20);
		RooGaussian gauss("gauss","gauss",bMass,mg,sg);

		RooFFTConvPdf sig("sig","hist (X) gauss",bMass,myHist,gauss);

		// RooPlot *frameSigma_gen = bMass.frame();
		// frameSigma_gen->SetTitle("J/#psi #Sigma");
		// genHist->plotOn(frameSigma_gen,Name("Sigmadata"));
		// sigma_gen.plotOn(frameSigma_gen,Name("Sigmafit"),LineColor(kBlue));
		// lst1405shape_smooth.plotOn(framelst1405,Name("lst1405fitsmooth"),LineColor(kRed),LineStyle(kDashed));
		// new TCanvas();
		// framelst1405->Draw();

		RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
		RooExponential bkg("bkg","Exponential bkg",bMass,tau);

		fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_%s.root",run,stage.Data()),"READ");
		treeIn = (TTree*)fileIn->Get("MyTuple");

		Int_t nentries = treeIn->GetEntries();

		treeIn->SetBranchStatus("*",0);
		treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

		treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>h_sigma(%d,%d,%d)",nbins,low,high),"","goff");
		TH1D *h_sigma = (TH1D*)gDirectory->Get("h_sigma");

		RooDataHist *dataHist = new RooDataHist("dataHist","dataHist",bMass,h_sigma);

		dataHist->Print();
		RooRealVar nsig("nsig","nsig",1,nentries);
		RooRealVar nbkg("nbkg","nbkg",1,nentries);

		RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));

		res = model.fitTo(*dataHist,Extended(),RooFit::Save(true), Strategy(2));

		RooPlot *frame = new RooPlot(bMass, low,5800,100);
		//  frame->GetYaxis()->SetRangeUser(0,100);
		dataHist->plotOn(frame,Name("data"));
		model.plotOn(frame,Components(sig),Name("sig"),LineColor(kMagenta+2));
		model.plotOn(frame,Components(bkg),LineColor(kRed),Name("bkg"));
		model.plotOn(frame,Name("fit"));

		new TCanvas();
		frame->Draw();

		Double_t chiSquare1 = frame->chiSquare("fit","data");
		cout<<"chi square1/dof = "<<chiSquare1<<endl;
		RooArgSet *floatpar = model.getParameters(dataHist);
		floatpar->Print();
		int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
		cout<<"float pars = "<<floatpars<<endl;
		Double_t chi2 = frame->chiSquare("fit","data",floatpars);
		cout<<"chi square2/dof = "<<chi2<<endl;
	}
	if(mcType == 3)
	{
		cout<<"poop2"<<endl;
		fileIn_gen = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiXi/run%d/jpsixi.root",run),"READ");
		treeIn_gen = (TTree*)fileIn_gen->Get("MCTuple/MCDecayTree");

		cout<<"poop3"<<endl;
		treeIn_gen->Draw(Form("sqrt(pow(J_psi_1S_TRUEP_E + Lambda0_TRUEP_E,2) - pow(J_psi_1S_TRUEP_X + Lambda0_TRUEP_X,2) - pow(J_psi_1S_TRUEP_Y + Lambda0_TRUEP_Y,2) - pow(J_psi_1S_TRUEP_Z + Lambda0_TRUEP_Z,2))>>h_xi_gen(%d,%d,%d)",
		                      nbins,low,high),"","goff");

		cout<<"poop4"<<endl;
		TH1D *h_xi_gen = (TH1D*)gDirectory->Get("h_xi_gen");

		RooDataHist *genHist = new RooDataHist("genHist","genHist",bMass,h_xi_gen);
		genHist->Print();

		RooHistPdf sig("sig","xi_gen",bMass,*genHist,0);//Signal shape, hopefully

		// RooPlot *frameXi_gen = bMass.frame();
		// frameXi_gen->SetTitle("J/#psi #Xi");
		// genHist->plotOn(frameXi_gen,Name("Xidata"));
		// xi_gen.plotOn(frameXi_gen,Name("Xifit"),LineColor(kBlue));
		// lst1405shape_smooth.plotOn(framelst1405,Name("lst1405fitsmooth"),LineColor(kRed),LineStyle(kDashed));
		// new TCanvas();
		// framelst1405->Draw();

		RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
		RooExponential bkg("bkg","Exponential bkg",bMass,tau);

		fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiXi/run%d/jpsixi_%s.root",run,stage.Data()),"READ");
		treeIn = (TTree*)fileIn->Get("MyTuple");

		Int_t nentries = treeIn->GetEntries();

		treeIn->SetBranchStatus("*",0);
		treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

		treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>h_xi(%d,%d,%d)",nbins,low,high),"","goff");
		TH1D *h_xi = (TH1D*)gDirectory->Get("h_xi");

		RooDataHist *dataHist = new RooDataHist("dataHist","dataHist",bMass,h_xi);

		dataHist->Print();
		RooRealVar nsig("nsig","nsig",1,nentries);
		RooRealVar nbkg("nbkg","nbkg",1,nentries);

		RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));

		res = model.fitTo(*dataHist,Extended(),RooFit::Save(true), Strategy(2));

		RooPlot *frame = new RooPlot(bMass, low,5800,100);
		//  frame->GetYaxis()->SetRangeUser(0,100);
		dataHist->plotOn(frame,Name("data"));
		model.plotOn(frame,Components(sig),Name("sig"),LineColor(kMagenta+2));
		model.plotOn(frame,Components(bkg),LineColor(kRed),Name("bkg"));
		model.plotOn(frame,Name("fit"));

		new TCanvas();
		frame->Draw();

		Double_t chiSquare1 = frame->chiSquare("fit","data");
		cout<<"chi square1/dof = "<<chiSquare1<<endl;
		RooArgSet *floatpar = model.getParameters(dataHist);
		floatpar->Print();
		int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
		cout<<"float pars = "<<floatpars<<endl;
		Double_t chi2 = frame->chiSquare("fit","data",floatpars);
		cout<<"chi square2/dof = "<<chi2<<endl;
	}
}
