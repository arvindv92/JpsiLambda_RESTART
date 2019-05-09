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
	h->SetTitleSize(0.1);
	h->SetLineColor(iCol);
	h->SetMarkerColor(iCol);
	h->SetMarkerSize(0.7);
	h->SetMarkerStyle(20);
	h->GetXaxis()->SetTitleOffset(1.0);
	h->GetYaxis()->SetTitleOffset(1.15);
	h->GetXaxis()->SetTitleSize(0.055);
	h->GetYaxis()->SetTitleSize(0.055);
	h->GetXaxis()->SetLabelSize(0.045);
	h->GetYaxis()->SetLabelSize(0.045);
	h->SetNdivisions(404,"X");
	h->SetNdivisions(505,"Y");
	h->SetLineWidth(2);
	h->SetTitle("");

}

void MakePlots()
{
	gROOT->ProcessLine(".x lhcbStyle.C");
	TString m_jpsiL = "m_{#jpsi#PLambda}[#rm{MeV}/c^{2}]";
	TString bin_4   = "Candidates/(4 #rm{MeV}/c^{2})";

	{
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

		TLatex *myLatex = new TLatex();
		myLatex->SetTextFont(42);
		myLatex->SetTextColor(1);
		myLatex->SetTextAlign(12);
		myLatex->SetNDC(kTRUE);

		myLatex->SetTextSize(0.065);


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
}
