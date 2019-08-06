#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TList.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
using namespace RooFit;
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
	TString bin_half   = "Candidates/(0.5 MeV/#it{c}^{2})";
	TString bin_1   = "Candidates/(1 MeV/#it{c}^{2})";

	TLatex *myLatex = new TLatex();
	myLatex->SetTextFont(42);
	myLatex->SetTextColor(1);
	myLatex->SetTextAlign(12);
	myLatex->SetNDC(kTRUE);

	myLatex->SetTextSize(0.065);

	TPaveText* lhcbName = new TPaveText(gStyle->GetPadLeftMargin() + 0.05,
	                                    0.87 - gStyle->GetPadTopMargin(),
	                                    gStyle->GetPadLeftMargin() + 0.20,
	                                    0.95 - gStyle->GetPadTopMargin(),
	                                    "BRNDC");
	lhcbName->AddText("LHCb");
	lhcbName->SetFillColor(0);
	lhcbName->SetTextAlign(12);
	lhcbName->SetBorderSize(0);

	TText *lhcbLabel = new TText();
	lhcbLabel->SetTextFont(132);
	lhcbLabel->SetTextColor(1);
	lhcbLabel->SetTextSize(0.06);
	lhcbLabel->SetTextAlign(12);

	TLatex *lhcbLatex = new TLatex();
	lhcbLatex->SetTextFont(132);
	lhcbLatex->SetTextColor(1);
	lhcbLatex->SetTextSize(0.06);
	lhcbLatex->SetTextAlign(12);

	// {
	//      //jpsi lambda mass in simulation for Lb-> JpsiLambda and Lb-> Jpsi Sigma
	//      //Run 2 used here.
	//      TFile *fileIn = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run2/jpsisigma.root");
	//      TTree *treeIn = (TTree*)fileIn->Get("Lb2JpsiLTree/MyTuple");
	//
	//      TFile *fileIn1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/jpsilambda.root");
	//      TTree *treeIn1 = (TTree*)fileIn1->Get("Lb2JpsiLTree/MyTuple");
	//
	//      treeIn->Draw("Lb_DTF_M_JpsiLConstr>>h0(150,5200,5800)","","goff");
	//      treeIn1->Draw("Lb_DTF_M_JpsiLConstr>>h1(100,5400,5800)","","goff");
	//
	//      TH1F *m_sigma = (TH1F*)gDirectory->Get("h0");
	//      TH1F *m_lambda = (TH1F*)gDirectory->Get("h1");
	//
	//      addGraphics(m_sigma,m_jpsiL,bin_4,1);
	//      addGraphics(m_lambda,m_jpsiL,bin_4,1);
	//
	//      TCanvas *c = new TCanvas("c","",1200,400);
	//      c->Divide(2,1);
	//
	//      c->cd(1);
	//      m_sigma->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//      // lhcbName->Draw();
	//      c->cd(2);
	//      m_lambda->Draw();
	//      // lhcbName->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      c->SaveAs("../plots/ANA/signal_shapes.pdf");
	// }
	// {
	//      //Ks -> Lambda misID veto plot
	//      TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_LL.root");
	//      TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");
	//
	//      TFile *fileIn1_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_LL.root");
	//      TTree *treeIn1_cut = (TTree*)fileIn1_cut->Get("MyTuple");
	//
	//      TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root");
	//      TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");
	//
	//      TFile *fileIn2_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_LL.root");
	//      TTree *treeIn2_cut = (TTree*)fileIn2_cut->Get("MyTuple");
	//
	//      treeIn1->Draw("L_WMpipi>>h2(100,300,700)","","goff");
	//      treeIn1_cut->Draw("L_WMpipi>>h3(100,300,700)","","goff");
	//
	//      treeIn2->Draw("L_WMpipi>>h4(100,300,700)","","goff");
	//      treeIn2_cut->Draw("L_WMpipi>>h5(100,300,700)","","goff");
	//
	//      TH1F *h2 = (TH1F*)gDirectory->Get("h2");
	//      TH1F *h3 = (TH1F*)gDirectory->Get("h3");
	//      TH1F *h4 = (TH1F*)gDirectory->Get("h4");
	//      TH1F *h5 = (TH1F*)gDirectory->Get("h5");
	//
	//      addGraphics(h2,m_pipi,bin_4,1);
	//      addGraphics(h3,m_pipi,bin_4,2);
	//      addGraphics(h4,m_pipi,bin_4,1);
	//      addGraphics(h5,m_pipi,bin_4,2);
	//
	//      TCanvas *c1 = new TCanvas("c1","",600,400);
	//      c1->SetTopMargin(0.07);
	//      h2->Draw();
	//      h3->Draw("same");
	//      lhcbName->Draw();
	//
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");
	//      TCanvas *d1 = new TCanvas("d1","",600,400);
	//      d1->SetTopMargin(0.07);
	//      h4->Draw();
	//      h5->Draw("same");
	//      lhcbName->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");
	//
	//      c1->SaveAs("../plots/ANA/lambda_veto_run1.pdf");
	//      d1->SaveAs("../plots/ANA/lambda_veto_run2.pdf");
	//
	// }
	// {
	//      //mass distributions after cutoutks
	//      TFile *fileIn1_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_LL.root");
	//      TTree *treeIn1_cut = (TTree*)fileIn1_cut->Get("MyTuple");
	//
	//      TFile *fileIn2_cut = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_LL.root");
	//      TTree *treeIn2_cut = (TTree*)fileIn2_cut->Get("MyTuple");
	//
	//      treeIn1_cut->Draw("Lb_DTF_M_JpsiLConstr>>h6(150,5200,5800)","","goff");
	//      treeIn2_cut->Draw("Lb_DTF_M_JpsiLConstr>>h7(150,5200,5800)","","goff");
	//
	//      TH1F *h6 = (TH1F*)gDirectory->Get("h6");
	//      TH1F *h7 = (TH1F*)gDirectory->Get("h7");
	//
	//      addGraphics(h6,m_jpsiL,bin_4,1);
	//      addGraphics(h7,m_jpsiL,bin_4,1);
	//
	//      TCanvas *c2 = new TCanvas("c2","",600,400);
	//
	//      h6->Draw();
	//      lhcbName->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");
	//
	//      TCanvas *d2 = new TCanvas("d2","",600,400);
	//
	//      h7->Draw();
	//      h7->GetYaxis()->SetTitleOffset(1.0);
	//      lhcbName->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");
	//
	//      c2->SaveAs("../plots/ANA/mass_cutoutks_run1.pdf");
	//      d2->SaveAs("../plots/ANA/mass_cutoutks_run2.pdf");
	//
	// }
	// {
	//      // isolation BDT output before final BDT
	//      TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LL_iso2_v0.root");
	//      TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");
	//
	//      TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_LL_iso1_v0.root");
	//      TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");
	//
	//      treeIn1->Draw("BDTkMin_v0>>h8(100,-1.0,1.0)","","goff");
	//      treeIn2->Draw("BDTkMin_v0>>h9(100,-1.0,1.0)","","goff");
	//
	//      TH1F *h8 = (TH1F*)gDirectory->Get("h8");
	//      TH1F *h9 = (TH1F*)gDirectory->Get("h9");
	//
	//      addGraphics(h8,"isolation BDT",bin_4,1);
	//      addGraphics(h9,"isolation BDT",bin_4,1);
	//
	//      TCanvas *c3 = new TCanvas("c3","",600,400);
	//
	//      h8->Draw();
	//      lhcbName->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Run 1");
	//
	//      TCanvas *d3 = new TCanvas("d3","",600,400);
	//      h9->Draw();
	//      lhcbName->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Run 2");
	//
	//      c3->SaveAs("../plots/ANA/isolation_run1.pdf");
	//      d3->SaveAs("../plots/ANA/isolation_run2.pdf");
	// }
	// {
	//      //mass distributions after finalBDT cut
	//      TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambdaoutks_LL_nonZeroTracks.root");
	//      TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");
	//
	//      treeIn1->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LL_FinalBDT2_iso2_v0.root");
	//
	//      TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambdaoutks_LL_nonZeroTracks.root");
	//      TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");
	//
	//      treeIn2->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_LL_FinalBDT2_iso1_v0.root");
	//
	//      TFile *fileIn1_zero = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambdaoutks_LL_ZeroTracks.root");
	//      TTree *treeIn1_zero = (TTree*)fileIn1_zero->Get("MyTuple");
	//
	//      treeIn1_zero->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_zeroTracksLL_FinalBDT2.root");
	//
	//      TFile *fileIn2_zero = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambdaoutks_LL_ZeroTracks.root");
	//      TTree *treeIn2_zero = (TTree*)fileIn2_zero->Get("MyTuple");
	//
	//      treeIn2_zero->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_zeroTracksLL_FinalBDT2.root");
	//
	//      treeIn1->Draw("Lb_DTF_M_JpsiLConstr>>g1(150,5200,5800)","BDT2 > 0.475","goff");
	//      treeIn2->Draw("Lb_DTF_M_JpsiLConstr>>g2(150,5200,5800)","BDT2 > 0.555","goff");
	//
	//      treeIn1_zero->Draw("Lb_DTF_M_JpsiLConstr>>g3(150,5200,5800)","BDT2 > 0.365","goff");
	//      treeIn2_zero->Draw("Lb_DTF_M_JpsiLConstr>>g4(150,5200,5800)","BDT2 > 0.495","goff");
	//
	//      treeIn1->Draw("Lb_DTF_M_JpsiLConstr>>f1(750,4000,7000)","BDT2 > 0.475","goff");
	//      treeIn2->Draw("Lb_DTF_M_JpsiLConstr>>f2(750,4000,7000)","BDT2 > 0.555","goff");
	//
	//      treeIn1_zero->Draw("Lb_DTF_M_JpsiLConstr>>f3(750,4000,7000)","BDT2 > 0.365","goff");
	//      treeIn2_zero->Draw("Lb_DTF_M_JpsiLConstr>>f4(750,4000,7000)","BDT2 > 0.495","goff");
	//
	//      TH1F *g1 = (TH1F*)gDirectory->Get("g1");
	//      TH1F *g2 = (TH1F*)gDirectory->Get("g2");
	//      TH1F *g3 = (TH1F*)gDirectory->Get("g3");
	//      TH1F *g4 = (TH1F*)gDirectory->Get("g4");
	//
	//      TH1F *f1 = (TH1F*)gDirectory->Get("f1");
	//      TH1F *f2 = (TH1F*)gDirectory->Get("f2");
	//      TH1F *f3 = (TH1F*)gDirectory->Get("f3");
	//      TH1F *f4 = (TH1F*)gDirectory->Get("f4");
	//
	//      TH1F *tot1 = (TH1F*)g1->Clone();
	//      tot1->Add(g3);//Add histograms for nonZeroTracks and ZeroTracks
	//
	//      TH1F *tot2 = (TH1F*)g2->Clone();
	//      tot2->Add(g4);//Add histograms for nonZeroTracks and ZeroTracks
	//
	//      TH1F *tot1_wide = (TH1F*)f1->Clone();
	//      tot1_wide->Add(f3);//Add histograms for nonZeroTracks and ZeroTracks
	//
	//      TH1F *tot2_wide = (TH1F*)f2->Clone();
	//      tot2_wide->Add(f4);//Add histograms for nonZeroTracks and ZeroTracks
	//
	//      TH1F *tot1_clone = (TH1F*)tot1->Clone();
	//      TH1F *tot2_clone = (TH1F*)tot2->Clone();
	//
	//      addGraphics(g1,m_jpsiL,bin_4,1);
	//      addGraphics(g2,m_jpsiL,bin_4,1);
	//      addGraphics(g3,m_jpsiL,bin_4,1);
	//      addGraphics(g4,m_jpsiL,bin_4,1);
	//
	//      addGraphics(f1,m_jpsiL,bin_4,1);
	//      addGraphics(f2,m_jpsiL,bin_4,1);
	//      addGraphics(f3,m_jpsiL,bin_4,1);
	//      addGraphics(f4,m_jpsiL,bin_4,1);
	//
	//      addGraphics(tot1,m_jpsiL,bin_4,1);
	//      addGraphics(tot2,m_jpsiL,bin_4,1);
	//      addGraphics(tot1_wide,m_jpsiL,bin_4,1);
	//      addGraphics(tot2_wide,m_jpsiL,bin_4,1);
	//
	//      addGraphics(tot1_clone,m_jpsiL,bin_4,1);
	//      addGraphics(tot2_clone,m_jpsiL,bin_4,1);
	//
	//      tot1_clone->SetMaximum(40);
	//      tot2_clone->SetMaximum(100);
	//
	//      TCanvas *c4 = new TCanvas("c4","",600,400);
	//      g1->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *c5 = new TCanvas("c5","",600,400);
	//      g2->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *c6 = new TCanvas("c6","",600,400);
	//      g3->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *c7 = new TCanvas("c7","",600,400);
	//      g4->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *d4 = new TCanvas("d4","",600,400);
	//      f1->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *d5 = new TCanvas("d5","",600,400);
	//      f2->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *d6 = new TCanvas("d6","",600,400);
	//      f3->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *d7 = new TCanvas("d7","",600,400);
	//      f4->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *e4 = new TCanvas("e4","",600,400);
	//      tot1->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *e5 = new TCanvas("e5","",600,400);
	//      tot2->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *e6 = new TCanvas("e6","",600,400);
	//      tot1_wide->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *e7 = new TCanvas("e7","",600,400);
	//      tot2_wide->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *e8 = new TCanvas("e8","",600,400);
	//      tot1_clone->Draw();
	//      lhcbName->Draw();
	//
	//      TCanvas *e9 = new TCanvas("e9","",600,400);
	//      tot2_clone->Draw();
	//      lhcbName->Draw();
	//
	//      c4->SaveAs("../plots/ANA/mass_run1_finalBDT_nonZeroTracks.pdf");
	//      c5->SaveAs("../plots/ANA/mass_run2_finalBDT_nonZeroTracks.pdf");
	//      c6->SaveAs("../plots/ANA/mass_run1_finalBDT_ZeroTracks.pdf");
	//      c7->SaveAs("../plots/ANA/mass_run2_finalBDT_ZeroTracks.pdf");
	//
	//      d4->SaveAs("../plots/ANA/mass_run1_finalBDT_nonZeroTracks_wide.pdf");
	//      d5->SaveAs("../plots/ANA/mass_run2_finalBDT_nonZeroTracks_wide.pdf");
	//      d6->SaveAs("../plots/ANA/mass_run1_finalBDT_ZeroTracks_wide.pdf");
	//      d7->SaveAs("../plots/ANA/mass_run2_finalBDT_ZeroTracks_wide.pdf");
	//
	//      e4->SaveAs("../plots/ANA/mass_run1_finalBDT_total.pdf");
	//      e5->SaveAs("../plots/ANA/mass_run2_finalBDT_total.pdf");
	//      e6->SaveAs("../plots/ANA/mass_run1_finalBDT_total_wide.pdf");
	//      e7->SaveAs("../plots/ANA/mass_run2_finalBDT_total_wide.pdf");
	//
	//      e8->SaveAs("../plots/ANA/mass_run1_finalBDT_total_zoomY.pdf");
	//      e9->SaveAs("../plots/ANA/mass_run2_finalBDT_total_zoomY.pdf");
	// }
	// {
	//      //J/psi Xi mass reco'd as J/psi Lambda, after BDT selection
	//      TFile *file1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiXi/run1/jpsixi_cutoutks_LL_nonZeroTracks.root");
	//      TTree *tree1 = (TTree*)file1->Get("MyTuple");
	//      tree1->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiXi/run1/jpsixi_LL_FinalBDT2_iso2_v0.root");
	//
	//      TFile *file2 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiXi/run1/jpsixi_cutoutks_LL_ZeroTracks.root");
	//      TTree *tree2 = (TTree*)file2->Get("MyTuple");
	//      tree2->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiXi/run1/jpsixi_zeroTracksLL_FinalBDT2.root");
	//
	//      TFile *file3 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiXi/run2/jpsixi_cutoutks_LL_nonZeroTracks.root");
	//      TTree *tree3 = (TTree*)file3->Get("MyTuple");
	//      tree3->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiXi/run2/jpsixi_LL_FinalBDT2_iso1_v0.root");
	//
	//      TFile *file4 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiXi/run2/jpsixi_cutoutks_LL_ZeroTracks.root");
	//      TTree *tree4 = (TTree*)file4->Get("MyTuple");
	//      tree4->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiXi/run2/jpsixi_zeroTracksLL_FinalBDT2.root");
	//
	//      TTree *tree1_cut = (TTree*)tree1->CopyTree("BDT2 > 0.475");
	//      TTree *tree2_cut = (TTree*)tree2->CopyTree("BDT2 > 0.365");
	//
	//      TTree *tree3_cut = (TTree*)tree3->CopyTree("BDT2 > 0.555");
	//      TTree *tree4_cut = (TTree*)tree4->CopyTree("BDT2 > 0.495");
	//
	//      TList *list1 = new TList;
	//      list1->Add(tree1_cut);
	//      list1->Add(tree2_cut);
	//
	//      TFile *tempFile = new TFile("tempFile.root","RECREATE");
	//
	//      TTree *combTree1 = TTree::MergeTrees(list1);
	//      combTree1->SetName("combTree1");
	//
	//      TList *list2 = new TList;
	//      list2->Add(tree3_cut);
	//      list2->Add(tree4_cut);
	//
	//      TTree *combTree2 = TTree::MergeTrees(list2);
	//      combTree2->SetName("combTree2");
	//
	//      RooRealVar *gbWtVar  = new RooRealVar("GB_WT","gb Weight Var",-100.,100.);
	//
	//      RooRealVar *Lb_Mass = new RooRealVar("Lb_DTF_M_JpsiLConstr","",5200,5800);
	//
	//      RooDataSet *ds1_nowt = new RooDataSet("ds1_nowt","ds1_nowt",combTree1,RooArgSet(*Lb_Mass,*gbWtVar));
	//      RooDataSet *ds2_nowt = new RooDataSet("ds2_nowt","ds2_nowt",combTree2,RooArgSet(*Lb_Mass,*gbWtVar));
	//
	//      RooDataSet *ds1 = new RooDataSet("ds1","ds1",RooArgSet(*Lb_Mass,*gbWtVar),Import(*ds1_nowt),WeightVar(*gbWtVar));
	//      RooDataSet *ds2 = new RooDataSet("ds2","ds2",RooArgSet(*Lb_Mass,*gbWtVar),Import(*ds2_nowt),WeightVar(*gbWtVar));
	//
	//      RooKeysPdf *xibFit1 = new RooKeysPdf("xibFit1","xibFit1",*Lb_Mass,*ds1,RooKeysPdf::NoMirror);
	//      RooKeysPdf *xibFit2 = new RooKeysPdf("xibFit2","xibFit2",*Lb_Mass,*ds2,RooKeysPdf::NoMirror);
	//
	//      RooPlot *frame1 = Lb_Mass->frame();
	//      ds1->plotOn(frame1);
	//      xibFit1->plotOn(frame1);
	//      frame1->GetXaxis()->SetTitle(m_jpsiL);
	//      frame1->GetYaxis()->SetTitle(bin_4);
	//
	//      RooPlot *frame2 = Lb_Mass->frame();
	//      ds2->plotOn(frame2);
	//      xibFit2->plotOn(frame2);
	//      frame2->GetXaxis()->SetTitle(m_jpsiL);
	//      frame2->GetYaxis()->SetTitle(bin_4);
	//
	//      // tree1->Draw("Lb_DTF_M_JpsiLConstr>>xi1(150,5200,5800)","GB_WT*(BDT2 > 0.475)","goff");
	//      // tree2->Draw("Lb_DTF_M_JpsiLConstr>>xi2(150,5200,5800)","GB_WT*(BDT2 > 0.365)","goff");
	//      // tree3->Draw("Lb_DTF_M_JpsiLConstr>>xi3(150,5200,5800)","GB_WT*(BDT2 > 0.555)","goff");
	//      // tree4->Draw("Lb_DTF_M_JpsiLConstr>>xi4(150,5200,5800)","GB_WT*(BDT2 > 0.495)","goff");
	//      //
	//      // TH1F *xi1 = (TH1F*)gDirectory->Get("xi1");
	//      // TH1F *xi2 = (TH1F*)gDirectory->Get("xi2");
	//      // TH1F *xi3 = (TH1F*)gDirectory->Get("xi3");
	//      // TH1F *xi4 = (TH1F*)gDirectory->Get("xi4");
	//      //
	//      // xi1->Add(xi2);
	//      // xi3->Add(xi4);
	//      //
	//      // addGraphics(xi1,m_jpsiL,bin_4,1);
	//      // addGraphics(xi3,m_jpsiL,bin_4,1);
	//      //
	//      // TCanvas *can_xi1 = new TCanvas();
	//      // xi1->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//      //
	//      // TCanvas *can_xi2 = new TCanvas();
	//      // xi3->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      TCanvas *can_xi1 = new TCanvas();
	//      frame1->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      TCanvas *can_xi2 = new TCanvas();
	//      frame2->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      can_xi1->SaveAs("../plots/ANA/jpsixi_jpsilambda_run1.pdf");
	//      can_xi2->SaveAs("../plots/ANA/jpsixi_jpsilambda_run2.pdf");
	// }
	//
	// {
	//      //J/psi Ks mass reco'd as J/psi Lambda, after BDT selection
	//      TFile *file1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_cutoutks_LL_nonZeroTracks.root");
	//      TTree *tree1 = (TTree*)file1->Get("MyTuple");
	//      tree1->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_LL_FinalBDT2_iso2_v0.root");
	//
	//      TFile *file2 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_cutoutks_LL_ZeroTracks.root");
	//      TTree *tree2 = (TTree*)file2->Get("MyTuple");
	//      tree2->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_zeroTracksLL_FinalBDT2.root");
	//
	//      TFile *file3 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_cutoutks_LL_nonZeroTracks.root");
	//      TTree *tree3 = (TTree*)file3->Get("MyTuple");
	//      tree3->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_LL_FinalBDT2_iso1_v0.root");
	//
	//      TFile *file4 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_cutoutks_LL_ZeroTracks.root");
	//      TTree *tree4 = (TTree*)file4->Get("MyTuple");
	//      tree4->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_zeroTracksLL_FinalBDT2.root");
	//
	//      // TTree *tree1_cut = (TTree*)tree1->CopyTree("BDT2 > 0.475");
	//      // TTree *tree2_cut = (TTree*)tree2->CopyTree("BDT2 > 0.365");
	//      //
	//      // TTree *tree3_cut = (TTree*)tree3->CopyTree("BDT2 > 0.555");
	//      // TTree *tree4_cut = (TTree*)tree4->CopyTree("BDT2 > 0.495");
	//      //
	//      // TList *list1 = new TList;
	//      // list1->Add(tree1_cut);
	//      // list1->Add(tree2_cut);
	//      //
	//      // TTree *combTree1 = TTree::MergeTrees(list1);
	//      // combTree1->SetName("combTree1");
	//      //
	//      // TList *list2 = new TList;
	//      // list2->Add(tree3_cut);
	//      // list2->Add(tree4_cut);
	//      //
	//      // TTree *combTree2 = TTree::MergeTrees(list2);
	//      // combTree2->SetName("combTree2");
	//
	//      // RooRealVar *Lb_Mass = new RooRealVar("Lb_DTF_M_JpsiLConstr","",5200,5800);
	//      //
	//      // RooDataSet *ds1 = new RooDataSet("ds1","ds1",combTree1,RooArgSet(*Lb_Mass));
	//      // RooDataSet *ds2 = new RooDataSet("ds2","ds2",combTree2,RooArgSet(*Lb_Mass));
	//      //
	//      // RooKeysPdf *xibFit1 = new RooKeysPdf("xibFit1","xibFit1",*Lb_Mass,*ds1,RooKeysPdf::NoMirror);
	//      // RooKeysPdf *xibFit2 = new RooKeysPdf("xibFit2","xibFit2",*Lb_Mass,*ds2,RooKeysPdf::NoMirror);
	//      //
	//      // RooPlot *frame1 = Lb_Mass->frame();
	//      // ds1->plotOn(frame1);
	//      // xibFit1->plotOn(frame1);
	//      // frame1->GetXaxis()->SetTitle(m_jpsiL);
	//      // frame1->GetYaxis()->SetTitle(bin_4);
	//      //
	//      // RooPlot *frame2 = Lb_Mass->frame();
	//      // ds2->plotOn(frame2);
	//      // xibFit2->plotOn(frame2);
	//      // frame2->GetXaxis()->SetTitle(m_jpsiL);
	//      // frame2->GetYaxis()->SetTitle(bin_4);
	//
	//      tree1->Draw("Lb_DTF_M_JpsiLConstr>>ks1(150,5200,5800)","(BDT2 > 0.475)","goff");
	//      tree2->Draw("Lb_DTF_M_JpsiLConstr>>ks2(150,5200,5800)","(BDT2 > 0.365)","goff");
	//      tree3->Draw("Lb_DTF_M_JpsiLConstr>>ks3(150,5200,5800)","(BDT2 > 0.555)","goff");
	//      tree4->Draw("Lb_DTF_M_JpsiLConstr>>ks4(150,5200,5800)","(BDT2 > 0.495)","goff");
	//
	//      TH1F *ks1 = (TH1F*)gDirectory->Get("ks1");
	//      TH1F *ks2 = (TH1F*)gDirectory->Get("ks2");
	//      TH1F *ks3 = (TH1F*)gDirectory->Get("ks3");
	//      TH1F *ks4 = (TH1F*)gDirectory->Get("ks4");
	//
	//      ks1->Add(ks2);
	//      ks3->Add(ks4);
	//
	//      addGraphics(ks1,m_jpsiL,bin_4,1);
	//      addGraphics(ks3,m_jpsiL,bin_4,1);
	//
	//      TCanvas *can_ks1 = new TCanvas();
	//      ks1->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      TCanvas *can_ks2 = new TCanvas();
	//      ks3->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      // TCanvas *can_xi1 = new TCanvas();
	//      // frame1->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//      //
	//      // TCanvas *can_xi2 = new TCanvas();
	//      // frame2->Draw();
	//      // myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      can_ks1->SaveAs("../plots/ANA/jpsiks_jpsilambda_run1.pdf");
	//      can_ks2->SaveAs("../plots/ANA/jpsiks_jpsilambda_run2.pdf");
	// }
	//
	// {
	//      //J/psi Ks mass reco'd as J/psi Lambda, WMpipi after BDT selection
	//      TFile *file1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_cutoutks_LL_nonZeroTracks.root");
	//      TTree *tree1 = (TTree*)file1->Get("MyTuple");
	//      tree1->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_LL_FinalBDT2_iso2_v0.root");
	//
	//      TFile *file2 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_cutoutks_LL_ZeroTracks.root");
	//      TTree *tree2 = (TTree*)file2->Get("MyTuple");
	//      tree2->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_zeroTracksLL_FinalBDT2.root");
	//
	//      TFile *file3 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_cutoutks_LL_nonZeroTracks.root");
	//      TTree *tree3 = (TTree*)file3->Get("MyTuple");
	//      tree3->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_LL_FinalBDT2_iso1_v0.root");
	//
	//      TFile *file4 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_cutoutks_LL_ZeroTracks.root");
	//      TTree *tree4 = (TTree*)file4->Get("MyTuple");
	//      tree4->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_zeroTracksLL_FinalBDT2.root");
	//
	//      TFile *file5 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run1/jpsiks_cutoutks_LL.root");
	//      TTree *tree5 = (TTree*)file5->Get("MyTuple");
	//
	//      TFile *file6 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiKs/run2/jpsiks_cutoutks_LL.root");
	//      TTree *tree6 = (TTree*)file6->Get("MyTuple");
	//
	//      TTree *tree1_cut = (TTree*)tree1->CopyTree("BDT2 > 0.475 && Lb_BKGCAT==30");
	//      TTree *tree2_cut = (TTree*)tree2->CopyTree("BDT2 > 0.365 && Lb_BKGCAT==30");
	//
	//      TTree *tree3_cut = (TTree*)tree3->CopyTree("BDT2 > 0.555 && Lb_BKGCAT==30");
	//      TTree *tree4_cut = (TTree*)tree4->CopyTree("BDT2 > 0.495 && Lb_BKGCAT==30");
	//
	//      TList *list1_cut = new TList;
	//      list1_cut->Add(tree1_cut);
	//      list1_cut->Add(tree2_cut);
	//
	//      TList *list2_cut = new TList;
	//      list2_cut->Add(tree3_cut);
	//      list2_cut->Add(tree4_cut);
	//
	//      // TList *list1 = new TList;
	//      // list1->Add(tree1);
	//      // list1->Add(tree2);
	//      //
	//      // TList *list2 = new TList;
	//      // list2->Add(tree3);
	//      // list2->Add(tree4);
	//
	//      TFile *tempFile = new TFile("tempFile.root","RECREATE");
	//
	//      TTree *combTree1_cut = TTree::MergeTrees(list1_cut);
	//      combTree1_cut->SetName("combTree1_cut");
	//
	//      TTree *combTree2_cut = TTree::MergeTrees(list2_cut);
	//      combTree2_cut->SetName("combTree2_cut");
	//
	//      // TTree *combTree1 = TTree::MergeTrees(list1);
	//      // combTree1->SetName("combTree1");
	//      //
	//      // TTree *combTree2 = TTree::MergeTrees(list2);
	//      // combTree2->SetName("combTree2");
	//
	//      Int_t nentries_run1_cut = combTree1_cut->GetEntries();
	//      Int_t nentries_run2_cut = combTree2_cut->GetEntries();
	//
	//      Int_t nentries_run1 = tree5->GetEntries();
	//      Int_t nentries_run2 = tree6->GetEntries();
	//
	//      RooRealVar *Lb_Mass = new RooRealVar("Lb_DTF_L_WMpipi_JpsiConstr","",450,550);
	//      RooRealVar *Lb_BKGCAT = new RooRealVar("Lb_BKGCAT","",0,200);
	//
	//      RooDataSet *ds_run1_cut = new RooDataSet("ds_run1_cut","ds_run1_cut",combTree1_cut,RooArgSet(*Lb_Mass));
	//      RooDataSet *ds_run2_cut = new RooDataSet("ds_run2_cut","ds_run2_cut",combTree2_cut,RooArgSet(*Lb_Mass));
	//
	//      RooDataSet *ds_run1 = new RooDataSet("ds_run1","ds_run1",tree5,RooArgSet(*Lb_Mass,*Lb_BKGCAT),"Lb_BKGCAT==30");
	//      RooDataSet *ds_run2 = new RooDataSet("ds_run2","ds_run2",tree6,RooArgSet(*Lb_Mass,*Lb_BKGCAT),"Lb_BKGCAT==30");
	//
	//      RooRealVar mean_gaus_run1("mean_gaus_run1","Gaussian Mean",497.0,490.0,500.0,"MeV");
	//      RooRealVar sigma1_gaus_run1("sigma1_gaus_run1","Gaussian sigma1",10.0,1.0,20.0,"MeV");
	//      RooRealVar sigma2_gaus_run1("sigma2_gaus_run1","Gaussian sigma2",10.0,1.0,20.0,"MeV");
	//
	//      RooGaussian sig1_gaus_run1("sig1_gaus_run1","sig",*Lb_Mass,mean_gaus_run1,sigma1_gaus_run1);
	//      RooGaussian sig2_gaus_run1("sig2_gaus_run1","sig",*Lb_Mass,mean_gaus_run1,sigma2_gaus_run1);
	//
	//      RooRealVar frac1_gaus_run1("frac1_gaus_run1","Fraction of sig1 in signal",0.5,0.0,1.0);
	//      RooAddPdf sig_gaus_run1("sig_gaus_run1","Total Gaussian signal",RooArgList(sig1_gaus_run1,sig2_gaus_run1),frac1_gaus_run1);
	//
	//      RooRealVar c0_run1("c0_run1","c0_run1",0.,-5.,5.);
	//      RooRealVar c1_run1("c1_run1","c1_run1",0.,-1.,1.);
	//
	//      RooChebychev bkg_run1("bkg_run1","Chebychev bkg",*Lb_Mass,RooArgSet(c0_run1,c1_run1));
	//
	//      // RooRealVar tau_run1("tau_run1","tau_run1",-0.0007,-0.01,-0.0000001);
	//      // RooExponential bkg_run1("bkg_run1","Exponential bkg_run1",*Lb_Mass,tau_run1);
	//
	//      RooRealVar nsig_run1_cut("nsig_run1_cut","nsig_run1_cut",0,nentries_run1_cut);
	//      RooRealVar nbkg_run1_cut("nbkg_run1_cut","nbkg_run1_cut",0,nentries_run1_cut);
	//
	//      RooRealVar nsig_run1("nsig_run1","nsig_run1",0,nentries_run1);
	//      RooRealVar nbkg_run1("nbkg_run1","nbkg_run1",0,nentries_run1);
	//
	//      RooAddPdf model_run1_cut("model_run1_cut","Nominal Fit Model_run1",
	//                               RooArgList(sig_gaus_run1,bkg_run1),RooArgList(nsig_run1_cut,nbkg_run1_cut));
	//
	//      RooRealVar mean_gaus_run2("mean_gaus_run2","Gaussian Mean",497.0,490.0,500.0,"MeV");
	//      RooRealVar sigma1_gaus_run2("sigma1_gaus_run2","Gaussian sigma1",10.0,1.0,20.0,"MeV");
	//      RooRealVar sigma2_gaus_run2("sigma2_gaus_run2","Gaussian sigma2",10.0,1.0,20.0,"MeV");
	//
	//      RooGaussian sig1_gaus_run2("sig1_gaus_run2","sig",*Lb_Mass,mean_gaus_run2,sigma1_gaus_run2);
	//      RooGaussian sig2_gaus_run2("sig1_gaus_run2","sig",*Lb_Mass,mean_gaus_run2,sigma2_gaus_run2);
	//
	//      RooRealVar frac1_gaus_run2("frac1_gaus_run2","Fraction of sig1 in signal",0.5,0.0,1.0);
	//      RooAddPdf sig_gaus_run2("sig_gaus_run2","Total Gaussian signal",RooArgList(sig1_gaus_run2,sig2_gaus_run2),frac1_gaus_run2);
	//
	//      RooRealVar c0_run2("c0_run2","c0_run2",0.,-5.,5.);
	//      RooRealVar c1_run2("c1_run2","c1_run2",0.,-1.,1.);
	//
	//      RooChebychev bkg_run2("bkg_run2","Chebychev bkg",*Lb_Mass,RooArgSet(c0_run2,c1_run2));
	//
	//      // RooRealVar tau_run2("tau_run2","tau_run2",-0.0007,-0.01,-0.0000001);
	//      // RooExponential bkg_run2("bkg_run2","Exponential bkg_run2",*Lb_Mass,tau_run2);
	//
	//      RooRealVar nsig_run2_cut("nsig_run2_cut","nsig_run2_cut",0,nentries_run2_cut);
	//      RooRealVar nbkg_run2_cut("nbkg_run2_cut","nbkg_run2_cut",0,nentries_run2_cut);
	//
	//      RooRealVar nsig_run2("nsig_run2","nsig_run2",0,nentries_run2);
	//      RooRealVar nbkg_run2("nbkg_run2","nbkg_run2",0,nentries_run2);
	//
	//      RooAddPdf model_run2_cut("model_run2_cut","Nominal Fit Model_run2",
	//                               RooArgList(sig_gaus_run2,bkg_run2),RooArgList(nsig_run2_cut,nbkg_run2_cut));
	//      RooAddPdf model_run1("model_run1","Nominal Fit Model_run1",
	//                           RooArgList(sig_gaus_run1,bkg_run1),RooArgList(nsig_run1,nbkg_run1));
	//      RooAddPdf model_run2("model_run2","Nominal Fit Model_run2",
	//                           RooArgList(sig_gaus_run2,bkg_run2),RooArgList(nsig_run2,nbkg_run2));
	//
	//      model_run1_cut.fitTo(*ds_run1_cut,Extended());
	//      model_run2_cut.fitTo(*ds_run2_cut,Extended());
	//
	//      model_run1.fitTo(*ds_run1,Extended());
	//      model_run2.fitTo(*ds_run2,Extended());
	//
	//
	//      RooPlot *frame_cut_run1 = Lb_Mass->frame();
	//      frame_cut_run1->GetXaxis()->SetTitle(m_jpsiL);
	//      frame_cut_run1->GetYaxis()->SetTitle(bin_4);
	//      ds_run1_cut->plotOn(frame_cut_run1,Name("data_run1"));
	//      model_run1_cut.plotOn(frame_cut_run1,Name("fit"));
	//      model_run1_cut.plotOn(frame_cut_run1,Components(sig_gaus_run1),LineStyle(kDashed));
	//      //  myModel.plotOn(frame_cut_run1,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	//      //myModel.plotOn(frame_cut_run1,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	//      model_run1_cut.plotOn(frame_cut_run1,Components(bkg_run1),LineColor(kRed));
	//
	//      RooPlot *frame_cut_run2 = Lb_Mass->frame();
	//      frame_cut_run2->GetXaxis()->SetTitle(m_jpsiL);
	//      frame_cut_run2->GetYaxis()->SetTitle(bin_4);
	//      ds_run2_cut->plotOn(frame_cut_run2,Name("data_run2"));
	//      model_run2_cut.plotOn(frame_cut_run2,Name("fit"));
	//      model_run2_cut.plotOn(frame_cut_run2,Components(sig_gaus_run2),LineStyle(kDashed));
	//      //  myModel.plotOn(frame_cut_run2,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	//      //myModel.plotOn(frame_cut_run2,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	//      model_run2_cut.plotOn(frame_cut_run2,Components(bkg_run2),LineColor(kRed));
	//
	//      RooPlot *frame_run1 = Lb_Mass->frame();
	//      frame_run1->GetXaxis()->SetTitle(m_jpsiL);
	//      frame_run1->GetYaxis()->SetTitle(bin_4);
	//      ds_run1->plotOn(frame_run1,Name("data_run1"));
	//      model_run1.plotOn(frame_run1,Name("fit"));
	//      model_run1.plotOn(frame_run1,Components(sig_gaus_run1),LineStyle(kDashed));
	//      //  myModel.plotOn(frame_run1,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	//      //myModel.plotOn(frame_run1,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	//      model_run1.plotOn(frame_run1,Components(bkg_run1),LineColor(kRed));
	//
	//      RooPlot *frame_run2 = Lb_Mass->frame();
	//      frame_run2->GetXaxis()->SetTitle(m_jpsiL);
	//      frame_run2->GetYaxis()->SetTitle(bin_4);
	//      ds_run2->plotOn(frame_run2,Name("data_run2"));
	//      model_run2.plotOn(frame_run2,Name("fit"));
	//      model_run2.plotOn(frame_run2,Components(sig_gaus_run2),LineStyle(kDashed));
	//      //  myModel.plotOn(frame_run2,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	//      //myModel.plotOn(frame_run2,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	//      model_run2.plotOn(frame_run2,Components(bkg_run2),LineColor(kRed));
	//
	//      // RooKeysPdf *xibFit1 = new RooKeysPdf("xibFit1","xibFit1",*Lb_Mass,*ds1,RooKeysPdf::NoMirror);
	//      // RooKeysPdf *xibFit2 = new RooKeysPdf("xibFit2","xibFit2",*Lb_Mass,*ds2,RooKeysPdf::NoMirror);
	//      //
	//      // RooPlot *frame1 = Lb_Mass->frame();
	//      // ds1->plotOn(frame1);
	//      // xibFit1->plotOn(frame1);
	//      // frame1->GetXaxis()->SetTitle(m_jpsiL);
	//      // frame1->GetYaxis()->SetTitle(bin_4);
	//      //
	//      // RooPlot *frame2 = Lb_Mass->frame();
	//      // ds2->plotOn(frame2);
	//      // xibFit2->plotOn(frame2);
	//      // frame2->GetXaxis()->SetTitle(m_jpsiL);
	//      // frame2->GetYaxis()->SetTitle(bin_4);
	//
	//      TCanvas *ks_cut_run1 = new TCanvas();
	//      frame_cut_run1->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      TCanvas *ks_cut_run2 = new TCanvas();
	//      frame_cut_run2->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      ks_cut_run1->SaveAs("../plots/ANA/jpsiks_wmpipi_bdtcut_run1.pdf");
	//      ks_cut_run2->SaveAs("../plots/ANA/jpsiks_wmpipi_bdtcut_run2.pdf");
	//
	//      TCanvas *ks_run1 = new TCanvas();
	//      frame_run1->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      TCanvas *ks_run2 = new TCanvas();
	//      frame_run2->Draw();
	//      myLatex->DrawLatex(0.18,0.85,"LHCb Simulation");
	//
	//      ks_run1->SaveAs("../plots/ANA/jpsiks_wmpipi_run1.pdf");
	//      ks_run2->SaveAs("../plots/ANA/jpsiks_wmpipi_run2.pdf");
	// }
	{
		// J\psi and Lambda mass distributions after all selections. Data and MC overlaid.
		TFile *fileIn1 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_LL_nonZeroTracks.root");
		TTree *treeIn1 = (TTree*)fileIn1->Get("MyTuple");

		treeIn1->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LL_FinalBDT2_iso2_v0.root");

		TFile *fileIn2 = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_LL_nonZeroTracks.root");
		TTree *treeIn2 = (TTree*)fileIn2->Get("MyTuple");

		treeIn2->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_LL_FinalBDT2_iso1_v0.root");

		TFile *fileIn1_zero = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_LL_ZeroTracks.root");
		TTree *treeIn1_zero = (TTree*)fileIn1_zero->Get("MyTuple");

		treeIn1_zero->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_zeroTracksLL_FinalBDT2.root");

		TFile *fileIn2_zero = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_LL_ZeroTracks.root");
		TTree *treeIn2_zero = (TTree*)fileIn2_zero->Get("MyTuple");

		treeIn2_zero->AddFriend("MyTuple","../rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_zeroTracksLL_FinalBDT2.root");

		treeIn1->Draw("Jpsi_M>>jpsi1(110,3040,3150)","BDT2 > 0.475","goff");
		treeIn2->Draw("Jpsi_M>>jpsi2(110,3040,3150)","BDT2 > 0.555","goff");

		treeIn1_zero->Draw("Jpsi_M>>jpsi3(110,3040,3150)","BDT2 > 0.365","goff");
		treeIn2_zero->Draw("Jpsi_M>>jpsi4(110,3040,3150)","BDT2 > 0.495","goff");

		treeIn1->Draw("L_M>>lambda1(50,1105,1130)","BDT2 > 0.475","goff");
		treeIn2->Draw("L_M>>lambda2(50,1105,1130)","BDT2 > 0.555","goff");

		treeIn1_zero->Draw("L_M>>lambda3(50,1105,1130)","BDT2 > 0.365","goff");
		treeIn2_zero->Draw("L_M>>lambda4(50,1105,1130)","BDT2 > 0.495","goff");

		TH1F *jpsi1 = (TH1F*)gDirectory->Get("jpsi1");
		TH1F *jpsi2 = (TH1F*)gDirectory->Get("jpsi2");
		TH1F *jpsi3 = (TH1F*)gDirectory->Get("jpsi3");
		TH1F *jpsi4 = (TH1F*)gDirectory->Get("jpsi4");

		TH1F *lambda1 = (TH1F*)gDirectory->Get("lambda1");
		TH1F *lambda2 = (TH1F*)gDirectory->Get("lambda2");
		TH1F *lambda3 = (TH1F*)gDirectory->Get("lambda3");
		TH1F *lambda4 = (TH1F*)gDirectory->Get("lambda4");

		jpsi1->Add(jpsi3);
		jpsi2->Add(jpsi4);

		lambda1->Add(lambda3);
		lambda2->Add(lambda4);

		addGraphics(jpsi1,"m_{#mu#mu}[MeV/#it{c}^{2}]",bin_1,1);
		addGraphics(jpsi2,"m_{#mu#mu}[MeV/#it{c}^{2}]",bin_1,1);

		addGraphics(lambda1,"m_{p#pi}[MeV/#it{c}^{2}]",bin_half,1);
		addGraphics(lambda2,"m_{p#pi}[MeV/#it{c}^{2}]",bin_half,1);

		TFile *mcfileIn1 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_cutoutks_LL_nonZeroTracks.root");
		TTree *mcTreeIn1 = (TTree*)mcfileIn1->Get("MyTuple");

		mcTreeIn1->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_LL_FinalBDT2_iso2_v0.root");

		TFile *mcfileIn2 = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/jpsilambda_cutoutks_LL_nonZeroTracks.root");
		TTree *mcTreeIn2 = (TTree*)mcfileIn2->Get("MyTuple");

		mcTreeIn2->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/jpsilambda_LL_FinalBDT2_iso1_v0.root");

		TFile *mcfileIn1_zero = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_cutoutks_LL_ZeroTracks.root");
		TTree *mcTreeIn1_zero = (TTree*)mcfileIn1_zero->Get("MyTuple");

		mcTreeIn1_zero->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_zeroTracksLL_FinalBDT2.root");

		TFile *mcfileIn2_zero = TFile::Open("../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/jpsilambda_cutoutks_LL_ZeroTracks.root");
		TTree *mcTreeIn2_zero = (TTree*)mcfileIn2_zero->Get("MyTuple");

		mcTreeIn2_zero->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiLambda/run2/jpsilambda_zeroTracksLL_FinalBDT2.root");

		mcTreeIn1->Draw("Jpsi_M>>mcjpsi1(110,3040,3150)","GB_WT*wt_tau*(BDT2 > 0.475)","goff");
		mcTreeIn2->Draw("Jpsi_M>>mcjpsi2(110,3040,3150)","GB_WT*wt_tau*(BDT2 > 0.555)","goff");

		mcTreeIn1_zero->Draw("Jpsi_M>>mcjpsi3(110,3040,3150)","GB_WT*wt_tau*(BDT2 > 0.365)","goff");
		mcTreeIn2_zero->Draw("Jpsi_M>>mcjpsi4(110,3040,3150)","GB_WT*wt_tau*(BDT2 > 0.495)","goff");

		mcTreeIn1->Draw("L_M>>mclambda1(50,1105,1130)","GB_WT*wt_tau*(BDT2 > 0.475)","goff");
		mcTreeIn2->Draw("L_M>>mclambda2(50,1105,1130)","GB_WT*wt_tau*(BDT2 > 0.555)","goff");

		mcTreeIn1_zero->Draw("L_M>>mclambda3(50,1105,1130)","GB_WT*wt_tau*(BDT2 > 0.365)","goff");
		mcTreeIn2_zero->Draw("L_M>>mclambda4(50,1105,1130)","GB_WT*wt_tau*(BDT2 > 0.495)","goff");

		TH1F *mcjpsi1 = (TH1F*)gDirectory->Get("mcjpsi1");
		TH1F *mcjpsi2 = (TH1F*)gDirectory->Get("mcjpsi2");
		TH1F *mcjpsi3 = (TH1F*)gDirectory->Get("mcjpsi3");
		TH1F *mcjpsi4 = (TH1F*)gDirectory->Get("mcjpsi4");

		TH1F *mclambda1 = (TH1F*)gDirectory->Get("mclambda1");
		TH1F *mclambda2 = (TH1F*)gDirectory->Get("mclambda2");
		TH1F *mclambda3 = (TH1F*)gDirectory->Get("mclambda3");
		TH1F *mclambda4 = (TH1F*)gDirectory->Get("mclambda4");

		mcjpsi1->Add(mcjpsi3);
		mcjpsi2->Add(mcjpsi4);

		mclambda1->Add(mclambda3);
		mclambda2->Add(mclambda4);

		addGraphics(mcjpsi1,"m_{#mu#mu}[MeV/#it{c}^{2}]",bin_1,2);
		addGraphics(mcjpsi2,"m_{#mu#mu}[MeV/#it{c}^{2}]",bin_1,2);

		addGraphics(mclambda1,"m_{p#pi}[MeV/#it{c}^{2}]",bin_half,2);
		addGraphics(mclambda2,"m_{p#pi}[MeV/#it{c}^{2}]",bin_half,2);

		TCanvas *jpsi_can1 = new TCanvas();
		jpsi1->DrawNormalized();
		mcjpsi1->DrawNormalized("same");
		lhcbName->Draw();

		TCanvas *jpsi_can2 = new TCanvas();
		jpsi2->DrawNormalized();
		mcjpsi2->DrawNormalized("same");
		lhcbName->Draw();

		TCanvas *lambda_can1 = new TCanvas();
		lambda1->DrawNormalized();
		mclambda1->DrawNormalized("same");
		lhcbName->Draw();

		TCanvas *lambda_can2 = new TCanvas();
		lambda2->DrawNormalized();
		mclambda2->DrawNormalized("same");
		lhcbName->Draw();

		jpsi_can1->SaveAs("../plots/ANA/mass_jpsi_run1.pdf");
		jpsi_can2->SaveAs("../plots/ANA/mass_jpsi_run2.pdf");

		lambda_can1->SaveAs("../plots/ANA/mass_lambda_run1.pdf");
		lambda_can2->SaveAs("../plots/ANA/mass_lambda_run2.pdf");

	}
}
