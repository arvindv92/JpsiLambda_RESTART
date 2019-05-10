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
	h->SetLineColor(iCol);
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

		TCanvas *c1 = new TCanvas("c1","",1200,400);
		c1->Divide(2,1);

		c1->cd(1);
		h2->Draw();
		h3->Draw("same");
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");

		c1->cd(2);
		h4->Draw();
		h5->Draw("same");
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");

		c1->SaveAs("../plots/ANA/lambda_veto.pdf");
	}
	{
		//mass distributions after cutoutks
		TFile *fileIn1_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_LL_noPID.root");
		TTree *treeIn1_cut = (TTree*)fileIn1_cut->Get("MyTuple");

		TFile *fileIn2_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_LL_noPID.root");
		TTree *treeIn2_cut = (TTree*)fileIn2_cut->Get("MyTuple");

		treeIn1_cut->Draw("Lb_DTF_M_JpsiLConstr>>h6(250,5000,6000)","","goff");
		treeIn2_cut->Draw("Lb_DTF_M_JpsiLConstr>>h7(250,5000,6000)","","goff");

		TH1F *h6 = (TH1F*)gDirectory->Get("h6");
		TH1F *h7 = (TH1F*)gDirectory->Get("h7");

		addGraphics(h6,m_jpsiL,bin_4,1);
		addGraphics(h7,m_jpsiL,bin_4,1);

		TCanvas *c2 = new TCanvas("c2","",1200,400);
		c2->Divide(2,1);

		c2->cd(1);
		h6->Draw();
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");
		c2->cd(2);
		h7->Draw();
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");

		c2->SaveAs("../plots/ANA/mass_cutoutks.pdf");
	}
	{
		// isolation BDT output before final BDT
		TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LL_iso1_v0_noPID.root");
		TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");

		TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_LL_iso1_v0_noPID.root");
		TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");

		treeIn1->Draw("BDTkMin_v0>>h8(100,-1.0,1.0)","","goff");
		treeIn2->Draw("BDTkMin_v0>>h9(100,-1.0,1.0)","","goff");

		TH1F *h8 = (TH1F*)gDirectory->Get("h8");
		TH1F *h9 = (TH1F*)gDirectory->Get("h9");

		addGraphics(h8,"isolation BDT",bin_4,1);
		addGraphics(h9,"isolation BDT",bin_4,1);

		TCanvas *c3 = new TCanvas("c3","",1200,400);
		c3->Divide(2,1);

		c3->cd(1);
		h8->Draw();
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");
		c3->cd(2);
		h9->Draw();
		myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");

		c3->SaveAs("../plots/ANA/isolation.pdf");
	}
	// {
	//      TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/"
	//                                   "/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root");
	//      TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");
	//      TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/"
	//                                   "/jpsilambda_cutoutks_LL_nonZeroTracks_noPID.root");
	//      TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");
	//
	//      TFile *fileIn1_zero = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/"
	//                                        "jpsilambda_cutoutks_LL_ZeroTracks_noPID.root");
	//      TTree *treeIn1_zero = (TTree*)fileIn1_zero->Get("MyTuple");
	//      TFile *fileIn2_zero = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/"
	//                                        "jpsilambda_cutoutks_LL_ZeroTracks_noPID.root");
	//      TTree *treeIn2_zero = (TTree*)fileIn2_zero->Get("MyTuple");
	//      treeIn1->AddFriend("MyTuple","../rooFiles/dataFiles/JpsiLambda/run1/"
	//                         "jpsilambda_LL_FinalBDT2_iso1_v0_noPID.root");
	//      treeIn2->AddFriend("MyTuple","../rooFiles/dataFiles/JpsiLambda/run2/"
	//                         "jpsilambda_LL_FinalBDT2_iso1_v0_noPID.root");
	//      treeIn1_zero->AddFriend("MyTuple","../rooFiles/dataFiles/JpsiLambda/run1/"
	//                              "jpsilambda_zeroTracksLL_FinalBDT2_noPID.root");
	//      treeIn2_zero->AddFriend("MyTuple","../rooFiles/dataFiles/JpsiLambda/run2/"
	//                              "jpsilambda_zeroTracksLL_FinalBDT2_noPID.root");
	// }
}
