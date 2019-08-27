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
void Fitscript_dataEffs(Int_t run = 1, TString stage = "Trigger", Bool_t isData = true, Bool_t wtMC = false)
{
	gROOT->ProcessLine(".x lhcbStyle.C");
	gSystem->Exec("date");
	gSystem->Load("RooHypatia2_cpp.so"); //Load library for Hypatia shape

	TFile *fileIn = nullptr, *fileIn_Zero = nullptr, *fileIn_sim = nullptr;
	TTree *treeIn = nullptr, *treeIn_Zero = nullptr, *treeIn_sim = nullptr;
	TH1D  *hMass  = nullptr, *hMass_Zero  = nullptr, *hMass_sim  = nullptr;

	Int_t nentries = 0;
	Int_t low = 5400, high = 6000;
	Int_t binwidth = 2;
	Int_t nbins = (Int_t)(high-low)/(binwidth);

	RooWorkspace w("w");
	w.factory(Form("Lb_DTF_M_JpsiLConstr[%d,%d]",low,high));

	RooRealVar *myVar = w.var("Lb_DTF_M_JpsiLConstr");

	const char *prefix   = "";
	const char *wtexp    = "";
	const char *bdtwtexp = "";
	if(isData)
	{
		prefix = "../rootFiles/dataFiles";
	}
	else
	{
		prefix = "../rootFiles/mcFiles/JpsiLambda";
		if(wtMC && (run == 1))
		{
			wtexp  = "(GB_WT*wt_tau)";
			bdtwtexp = "*(GB_WT*wt_tau)";
		}
		else if(wtMC && (run == 2))
		{
			wtexp  = "(GB_WT*wt_tau)";
			bdtwtexp = "*(GB_WT*wt_tau)";
		}
	}
	//*************Fit to simulation to get params****************************

	// fileIn_sim = Open(Form("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_pidgen.root",run));
	// treeIn_sim = (TTree*)fileIn_sim->Get("MyTuple");
	//
	// treeIn_sim->Draw("Lb_DTF_M_JpsiLConstr>>hMass_sim(150,5500,5800)","Lb_BKGCAT==0||Lb_BKGCAT==50","goff");
	// hMass_sim = (TH1D*)gDirectory->Get("hMass_sim");
	//
	// RooDataHist *dh_sim = new RooDataHist("dh_sim","",*myVar,hMass_sim);

	// nentries = hMass_sim->GetEntries();
	//
	// w.factory("RooCBShape::Lb1(Lb_DTF_M_JpsiLConstr,mean_Run1[5619.6,5619,5621],"
	//           "sigma_Run1[10.,0.,20.], alpha1_Run1[1.028, 0.8,1.3], 10.0)" );
	// w.factory("RooCBShape::Lb2(Lb_DTF_M_JpsiLConstr,mean_Run1,sigma_Run1,alpha2_Run1[-1.097,-1.4,-0.7], 10.0)");
	// w.factory("SUM::sig(0.5*Lb1 , 0.5*Lb2)");
	//
	// cout<<"Done defining J/psi Lambda Hypatia shapes"<<endl;
	//
	// cout<<"*****UING EXPONENTIAL BKG SHAPE*****"<<endl;
	// w.factory("Exponential::bkg(Lb_DTF_M_JpsiLConstr,tau[-0.0007,-0.01,-0.0000001])");
	//
	// w.factory(Form("nsig[1,%d]",nentries));
	// w.factory(Form("nbkg[1,%d]",nentries));
	//
	// w.factory("SUM::model(nsig*sig , nbkg*bkg)");
	//
	// RooAbsPdf* model = w.pdf("model"); // get the model
	//
	// RooFitResult *res = model->fitTo(*dh_sim,Extended(), Save(), Hesse(true), Strategy(2));
	//
	// TCanvas* c = new TCanvas("c","", 1200, 800);
	//
	// RooPlot *frame = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),5500,5800,75);
	// frame->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	// frame->GetYaxis()->SetTitle("Candidates/(4 MeV/#it{c}^{2})");
	//
	// dh_sim->plotOn(frame,Name("data"),DataError(RooAbsData::Poisson));
	// model->plotOn(frame,Name("fit"));
	// model->plotOn(frame,Components(*(w.pdf("sig"))),Name("sig"),LineColor(kMagenta+2));
	// model->plotOn(frame,Components(*(w.pdf("bkg"))),LineColor(kRed),Name("bkg"));
	//
	// RooArgSet *allpar = model->getParameters(*dh_sim);
	// RooArgSet *floatpar = (RooArgSet*)allpar->selectByAttrib("Constant",kFALSE);
	// floatpar->Print();
	// int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
	// cout<<"run1 float pars = "<<floatpars<<endl;
	// Double_t chi2_run1 = frame->chiSquare("fit","data",floatpars);
	// cout<<"chi square2/dof = "<<chi2_run1<<endl;
	//
	// frame->Draw();
	// auto legend = new TLegend(0.65,0.7,0.85,0.9);
	// legend->SetTextSize(0.045);
	// legend->AddEntry("data","Data","lp");
	// legend->AddEntry("fit","Total Fit","l");
	// legend->AddEntry("lb","J/#psi #Lambda shape","l");
	// legend->AddEntry("bkg","Comb. Bkg. shape","l");
	// legend->Draw("same");
	//
	// c->Update();
	//*******************************************************************

	if(stage=="Trigger")
	{
		if(run == 1)
		{
			fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_triggered.root",prefix,run));
			treeIn = (TTree*)fileIn->Get("MyTuple");
			treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),wtexp);
			hMass = (TH1D*)gDirectory->Get("hMass");
		}
		else if(run == 2)
		{
			if(isData)
			{
				fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_triggered_2015.root",prefix,run));
				treeIn = (TTree*)fileIn->Get("MyTuple");
				treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),wtexp,"goff");
				hMass = (TH1D*)gDirectory->Get("hMass");

				fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_triggered_2016.root",prefix,run));
				treeIn = (TTree*)fileIn->Get("MyTuple");
				treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass1(%d,%d,%d)",nbins,low,high),wtexp,"goff");
				TH1D *hMass1 = (TH1D*)gDirectory->Get("hMass1");
				hMass->Add(hMass1);

				fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_triggered_2017.root",prefix,run));
				treeIn = (TTree*)fileIn->Get("MyTuple");
				treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass2(%d,%d,%d)",nbins,low,high),wtexp,"goff");
				TH1D *hMass2 = (TH1D*)gDirectory->Get("hMass2");
				hMass->Add(hMass2);

				fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_triggered_2018.root",prefix,run));
				treeIn = (TTree*)fileIn->Get("MyTuple");
				treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass3(%d,%d,%d)",nbins,low,high),wtexp,"goff");
				TH1D *hMass3 = (TH1D*)gDirectory->Get("hMass3");
				hMass->Add(hMass3);
			}
			else
			{
				fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_triggered.root",prefix,run));
				treeIn = (TTree*)fileIn->Get("MyTuple");
				treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),wtexp);
				hMass = (TH1D*)gDirectory->Get("hMass");
			}
			hMass->Draw();
		}
	}
	else if(stage=="Sanity")
	{
		fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_sanity_LL.root",prefix,run));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),wtexp);
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="CutOutKs")
	{
		fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_cutoutks_LL.root",prefix,run));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),wtexp);
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="finalBDT")
	{
		fileIn = Open(Form("%s/JpsiLambda/run%d/jpsilambda_cutoutks_LL_nonZeroTracks.root",prefix,run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		if(run == 1)
		{
			treeIn->AddFriend("MyTuple",Form("%s/JpsiLambda/run%d/jpsilambda_LL_FinalBDT1_iso2_v1.root",prefix,run));
			treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),Form("(BDT1 > 0.405)%s",bdtwtexp));
		}
		else if(run == 2)
		{
			treeIn->AddFriend("MyTuple",Form("%s/JpsiLambda/run%d/jpsilambda_LL_FinalBDT1_iso2_v1.root",prefix,run));
			treeIn->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass(%d,%d,%d)",nbins,low,high),Form("(BDT1 > 0.475)%s",bdtwtexp));
		}

		hMass = (TH1D*)gDirectory->Get("hMass");

		fileIn_Zero = Open(Form("%s/JpsiLambda/run%d/jpsilambda_cutoutks_LL_ZeroTracks.root",prefix,run));
		treeIn_Zero = (TTree*)fileIn_Zero->Get("MyTuple");
		treeIn_Zero->AddFriend("MyTuple",Form("%s/JpsiLambda/run%d/jpsilambda_zeroTracksLL_FinalBDT1.root",prefix,run));

		if(run == 1)
		{
			treeIn_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass_Zero(%d,%d,%d)",nbins,low,high),Form("(BDT1 > 0.335)%s",bdtwtexp));
		}
		else if(run == 2)
		{
			treeIn_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hMass_Zero(%d,%d,%d)",nbins,low,high),Form("(BDT1 > 0.405)%s",bdtwtexp));
		}
		hMass_Zero = (TH1D*)gDirectory->Get("hMass_Zero");

		hMass->Add(hMass_Zero);
	}
	nentries = (Int_t)(hMass->Integral()) + 1;
	cout<<"nentries = "<<nentries<<endl;

	RooDataHist *dh = new RooDataHist("dh","",*myVar,hMass);

	if(run == 1)
	{
		w.factory("RooCBShape::Lb1(Lb_DTF_M_JpsiLConstr,mean[5619.6,5619,5621],"
		          "sigma[7.28,5.0,9.0], alpha1[1.08453], 10.0)" );
		w.factory("RooCBShape::Lb2(Lb_DTF_M_JpsiLConstr,mean,sigma,alpha2[-1.11446], 10.0)");
	}
	if(run == 2)
	{
		w.factory("RooCBShape::Lb1(Lb_DTF_M_JpsiLConstr,mean[5619.6,5619,5621],"
		          "sigma[7.34,5.0,9.0], alpha1[0.98616], 10.0)" );
		w.factory("RooCBShape::Lb2(Lb_DTF_M_JpsiLConstr,mean,sigma,alpha2[-1.036], 10.0)");
	}

	w.factory("SUM::sig(0.5*Lb1 , 0.5*Lb2)");

	cout<<"Done defining J/psi Lambda Hypatia shapes"<<endl;

	cout<<"*****UING EXPONENTIAL BKG SHAPE*****"<<endl;
	w.factory("Exponential::bkg(Lb_DTF_M_JpsiLConstr,tau[-0.0007,-0.01,-0.0000001])");

	w.factory(Form("nsig[1,%d]",nentries));
	w.factory(Form("nbkg[1,%d]",nentries));

	w.factory("SUM::model(nsig*sig , nbkg*bkg)");

	RooAbsPdf* model = w.pdf("model"); // get the model

	// w.factory("RooHypatia2::sig(Lb_DTF_M_JpsiLConstr,lambda[-2.0,-4.0,0.0],0,0,"
	//           "sigma[10.,1.,20.], mean[5619.6,5619,5621], a1[1.7,1.0,3.0],"
	//           "2 ,a2[3.0,1.0,4.0], 2)");
	if(isData || !wtMC)
	{
		RooFitResult *res = model->fitTo(*dh,Extended(), Save(), Hesse(true), Strategy(2));
	}
	else
	{
		RooFitResult *res = model->fitTo(*dh,Extended(), Save(), Hesse(true), Strategy(2), SumW2Error(kTRUE));
	}
	TCanvas* c = new TCanvas("c","", 1200, 800);

	RooPlot *frame = new RooPlot(*(w.var("Lb_DTF_M_JpsiLConstr")),low,high,nbins);
	frame->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame->GetYaxis()->SetTitle("Candidates/(2 MeV/#it{c}^{2})");

	if(isData)
	{
		dh->plotOn(frame,Name("data"),DataError(RooAbsData::Poisson));
	}
	else
	{
		dh->plotOn(frame,Name("data"));
	}
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
	legend->AddEntry("sig","J/#psi #Lambda shape","l");
	legend->AddEntry("bkg","Comb. Bkg. shape","l");
	legend->Draw("same");

	c->Update();
}
