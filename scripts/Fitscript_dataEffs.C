#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TH1D.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"

using namespace std;
using namespace RooFit;

#define Open TFile::Open
void Fitscript_dataEffs(Int_t run = 1, TString stage = "Trigger")
{
	gROOT->ProcessLine(".x lhcbStyle.C");
	gSystem->Exec("date");
	gSystem->Load("RooHypatia2_cpp.so"); //Load library for Hypatia shape

	TFile *fileIn = nullptr, *fileIn_Zero = nullptr;
	TTree *treeIn = nullptr, *treeIn_Zero = nullptr;
	TH1D  *hMass  = nullptr, *hMass_Zero  = nullptr;

	if(stage=="Trigger")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="Sanity")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL_noPID.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="CutOutKs")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_LL_noPID.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");
	}
	else if(stage=="finalBDT")
	{
		fileIn = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root",run));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass(75,5500,5800)");
		hMass = (TH1D*)gDirectory->Get("hMass");

		fileIn_Zero = Open(Form("../rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_LL_ZeroTracks_noPID.root",run));
		treeIn_Zero = (TTree*)fileIn_Zero->Get("MyTuple");
		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>hMass_Zero(75,5500,5800)");
		hMass_Zero = (TH1D*)gDirectory->Get("hMass_Zero");

		hMass->Add(hMass_Zero);
	}

	RooWorkspace w("w");
	w.factory(Form("Lb_DTF_M_JpsiLConstr[%d,%d]",5500,5800));

	RooRealVar *myVar = w.var("Lb_DTF_M_JpsiLConstr");

	RooDataHist *dh = new RooDataHist("dh","dh",*myVar,hMass);


}
