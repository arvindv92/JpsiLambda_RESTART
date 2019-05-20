#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooExponential.h"
//#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"


using namespace RooFit;
using namespace std;

#define Open TFile::Open
void Fitscript_dataEffs(Int_t run = 1, TString stage = "Trigger")
{
	gROOT->ProcessLine(".x lhcbStyle.C");
	gSystem->Exec("date");
	gSystem->Load("RooHypatia2_cpp.so"); //Load library for Hypatia shape

	TFile *fileIn = nullptr, *fileIn_Zero = nullptr, *fileIn_sim = nullptr;
	TTree *treeIn = nullptr, *treeIn_Zero = nullptr;
	TH1D  *hMass  = nullptr, *hMass_Zero  = nullptr;

	Int_t nentries = 0;
	if(stage=="Trigger")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(150,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="Sanity")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL_noPID.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(150,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="CutOutKs")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_LL_noPID.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(150,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="finalBDT")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->AddFriend("MyTuple",Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_LL_FinalBDT2_iso1_v0_noPID.root",run));
		if(run == 1)
		{
			treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)","BDT2 > 0.475");
		}
		else if(run == 2)
		{
			treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)","BDT2 > 0.555");
		}
		hMass = (TH1D*)gDirectory->Get("hMass");

		fileIn_Zero = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_LL_ZeroTracks_noPID.root",run));
		treeIn_Zero = (TTree*)fileIn_Zero->Get("MyTuple");
		if(run == 1)
		{
			treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)","BDT2 > 0.365");
		}
		else if(run == 2)
		{
			treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)","BDT2 > 0.495");
		}
		hMass_Zero = (TH1D*)gDirectory->Get("hMass_Zero");

		hMass->Add(hMass_Zero);
	}
	nentries = hMass->GetEntries();

	RooWorkspace w("w");
	w.factory(Form("Lb_DTF_M_JpsiLConstr[%d,%d]",5500,5800));

	RooRealVar *myVar = w.var("Lb_DTF_M_JpsiLConstr");

	RooDataHist *dh = new RooDataHist("dh","",*myVar,hMass);

	// w.factory("RooHypatia2::sig(Lb_DTF_M_JpsiLConstr,lambda[-2.0,-4.0,0.0],0,0,"
	//           "sigma[10.,1.,20.], mean[5619.6,5619,5621], a1[1.7,1.0,3.0],"
	//           "2 ,a2[3.0,1.0,4.0], 2)");

	w.factory("Gaussian::sig1(Lb_DTF_M_JpsiLConstr,mean[5619.6,5619,5621],"
	          "sigma1[10.,1.,20.])");
	w.factory("Gaussian::sig2(Lb_DTF_M_JpsiLConstr,mean,"
	          "sigma2[15.,1.,20.])");
	w.factory("SUM::sig(f1[0.5,0.,1.]*sig1,sig2)");

	cout<<"Done defining J/psi Lambda Hypatia shapes"<<endl;

	cout<<"*****UING EXPONENTIAL BKG SHAPE*****"<<endl;
	w.factory("Exponential::bkg(Lb_DTF_M_JpsiLConstr,tau[-0.0007,-0.01,-0.0000001])");

	if(run == 1)
	{
		w.factory("nsig[6000,5000,10000]");
	}
	if(run == 2)
	{
		w.factory("nsig[20000,10000,25000]");
	}

	w.factory(Form("nbkg[1,%d]",nentries));

	w.factory("SUM::model(nsig*sig , nbkg*bkg)");

	RooAbsPdf* model = w.pdf("model"); // get the model

	RooFitResult *res = model->fitTo(*dh,Extended(), Save(), Hesse(true), Strategy(2));

	TCanvas* c = new TCanvas("c","", 1200, 800);

	RooPlot *frame = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),5500,5800,75);
	frame->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame->GetYaxis()->SetTitle("Candidates/(4 MeV/#it{c}^{2})");

	dh->plotOn(frame,Name("data"),DataError(RooAbsData::Poisson));
	model->plotOn(frame,Name("fit"));
	model->plotOn(frame,Components(*(w.pdf("sig"))),Name("sig"),LineColor(kMagenta+2));
	model->plotOn(frame,Components(*(w.pdf("bkg"))),LineColor(kRed),Name("bkg"));

	RooArgSet *allpar = model->getParameters(*dh);
	RooArgSet *floatpar = (RooArgSet*)allpar->selectByAttrib("Constant",kFALSE);
	floatpar->Print();
	int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	cout<<"run1 float pars = "<<floatpars<<endl;
	Double_t chi2_run1 = frame->chiSquare("fit","data",floatpars);
	cout<<"chi square2/dof = "<<chi2_run1<<endl;

	frame->Draw();
	auto legend = new TLegend(0.65,0.7,0.85,0.9);
	legend->SetTextSize(0.045);
	legend->AddEntry("data","Data","lp");
	legend->AddEntry("fit","Total Fit","l");
	legend->AddEntry("lb","J/#psi #Lambda shape","l");
	legend->AddEntry("bkg","Comb. Bkg. shape","l");
	legend->Draw("same");

	c->Update();
}
