#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"

void addGraphics(TH1F *h, TString Xtitle, TString Ytitle, int iCol){
	h->SetXTitle(Xtitle);
	//h->SetFillColor(30);
	//int bw = h->GetBinWidth(1);
	h->SetYTitle(Ytitle);
	h->SetStats(kFALSE);
	h->SetMinimum(0.1);
	h->SetMaximum(1.3*h->GetMaximum());
	// h->SetTitleSize(0.1);
	// h->SetLineColor(iCol);
	// h->SetMarkerColor(iCol);
	// h->SetMarkerSize(0.7);
	// h->SetMarkerStyle(20);
	// h->GetXaxis()->SetTitleOffset(1.0);
	// h->GetYaxis()->SetTitleOffset(1.15);
	// h->GetXaxis()->SetTitleSize(0.055);
	// h->GetYaxis()->SetTitleSize(0.055);
	// h->GetXaxis()->SetLabelSize(0.045);
	// h->GetYaxis()->SetLabelSize(0.045);
	// h->SetNdivisions(404,"X");
	// h->SetNdivisions(505,"Y");
	// h->SetLineWidth(2);
	h->SetTitle("");

}

void MakePlots()
{
	gROOT->ProcessLine(".x lhcbStyle.C");
	TString m_jpsiL = "m_{J/#psi#Lambda}[MeV/#it{c}^{2}]";
	TString m_pipi   = "m_{#pi#pi}[MeV/#it{c}^{2}]";
	TString bin_4   = "Candidates/(4 MeV/#it{c}^{2})";

	TLatex *myLatex = new TLatex();
	myLatex->SetTextFont(42);
	myLatex->SetTextColor(1);
	myLatex->SetTextAlign(12);
	myLatex->SetNDC(kTRUE);

	myLatex->SetTextSize(0.065);

	{
		//jpsi lambda mass in simulation for Lb-> JpsiLambda and Lb-> Jpsi Sigma
		//Run 2 used here.
		TFile *fileIn = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run2/jpsisigma.root");
		TTree *treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");

		TFile *fileIn1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/jpsilambda.root");
		TTree *treeIn1 = (TTree*)fileIn1->Get("Lb2JpsiLTree/MyTuple");

		treeIn->Draw("Lb_DTF_M_JpsiLConstr>>h0(250,5000,6000)","","goff");
		treeIn1->Draw("Lb_DTF_M_JpsiLConstr>>h1(100,5400,5800)","","goff");

		TH1F *m_sigma = (TH1F*)gDirectory->Get("h0");
		TH1F *m_lambda = (TH1F*)gDirectory->Get("h1");

		addGraphics(m_sigma,m_jpsiL,bin_4,1);
		addGraphics(m_lambda,m_jpsiL,bin_4,1);

		TCanvas *c = new TCanvas("c","",1200,400);
		c->Divide(2,1);

		c->cd(1);
		m_sigma->Draw();
		myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
		c->cd(2);
		m_lambda->Draw();
		myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");

		c->SaveAs("../plots/ANA/signal_shapes.pdf");
	}
	{
		//Ks -> Lambda misID veto plot
		TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_LL_noPID.root");
		TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");

		TFile *fileIn1_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_LL_noPID.root");
		TTree *treeIn1_cut = (TTree*)fileIn1_cut->Get("MyTuple");

		TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_LL_noPID.root");
		TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");

		TFile *fileIn2_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_LL_noPID.root");
		TTree *treeIn2_cut = (TTree*)fileIn2_cut->Get("MyTuple");

		treeIn1->Draw("L_WMpipi>>h2(100,300,700)","","goff");
		treeIn1_cut->Draw("L_WMpipi>>h3(100,300,700)","","goff");

		treeIn2->Draw("L_WMpipi>>h4(100,300,700)","","goff");
		treeIn2_cut->Draw("L_WMpipi>>h5(100,300,700)","","goff");

		TH1F *h2 = (TH1F*)gDirectory->Get("h2");
		TH1F *h3 = (TH1F*)gDirectory->Get("h3");
		TH1F *h4 = (TH1F*)gDirectory->Get("h4");
		TH1F *h5 = (TH1F*)gDirectory->Get("h5");

		addGraphics(h2,m_pipi,bin_4,1);
		addGraphics(h3,m_pipi,bin_4,2);
		addGraphics(h4,m_pipi,bin_4,1);
		addGraphics(h5,m_pipi,bin_4,2);

		TCanvas *c = new TCanvas("c","",1200,400);
		c->Divide(2,1);

		c->cd(1);
		h2->Draw();
		h3->Draw("same");
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");

		c->cd(2);
		h4->Draw();
		h5->Draw("same");
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");

		c->SaveAs("../plots/ANA/lambda_veto.pdf");
	}
}
