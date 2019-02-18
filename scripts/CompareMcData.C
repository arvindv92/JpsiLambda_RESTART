#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "iostream"

using namespace std;
void CompareMcData(Int_t run = 1, Int_t mcType = 0, Int_t trackType = 3, const char *varName = "", Float_t low = 0.0, Float_t high = 0.0, Int_t nBins = 50, Bool_t tmFlag = false)
{
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

  gStyle->SetOptStat(0);

  TFile *fileIn_data = nullptr, *fileIn_mc = nullptr;
  TTree *treeIn_data = nullptr, *treeIn_mc = nullptr;
  TString tmCut = "";

  if(tmFlag)
    tmCut = "(Lb_BKGCAT==0||Lb_BKGCAT==50)";
  Double_t myChi2 = 0.0;

  const char *type = (trackType == 3) ? ("LL") : ("DD");

  const char *folder = "", *part = "";
  switch(mcType)
    {
    case 0:
      {
	folder = "";
	part = "";
	break;
      }
    case 1:
      {
	folder = "JpsiLambda";
	part = "jpsilambda";
	break;
      }
    case 2:
      {
	folder = "JpsiSigma";
	part = "jpsisigma";
	break;
      }
    case 3:
      {
	folder = "JpsiXi";
	part = "jpsixi";
	break;
      }
    case 4:
      {
	folder = "Bu_JpsiX";
	part = "bu_jpsix";
	break;
      }
    case 5:
      {
	folder = "Bd_JpsiX";
	part = "bd_jpsix";
	break;
      }
    }

  fileIn_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/sWeightSanity/jpsilambda_%s_sanity_withsw.root",run,type),"READ");
  treeIn_data = (TTree*)fileIn_data->Get("MyTuple");

  fileIn_mc   = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_cutoutks_%s.root",folder,run,part,type),"READ");
  treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

  TCanvas *c1 = new TCanvas();
  // c1->Divide(2,1);
  // c1->cd(1);

  treeIn_data->Draw(Form("%s>>dataHist(%d,%f,%f)",varName,nBins,low,high),"SW","goff");
  treeIn_mc->Draw(Form("%s>>mcHist(%d,%f,%f)",varName,nBins,low,high),tmCut,"goff");

  TH1D *dataHist = (TH1D*)gDirectory->Get("dataHist");
  TH1D *mcHist = (TH1D*)gDirectory->Get("mcHist");

  dataHist->Scale(1.0/dataHist->Integral());
  mcHist->Scale(1.0/mcHist->Integral());
  
  for(Int_t i=1; i<=dataHist->GetNbinsX(); i++) {
    if(dataHist->GetBinContent(i) < 0)
      dataHist->SetBinContent(i,0);
  }

  myChi2 = mcHist->Chi2Test(dataHist,"UW CHI2/NDF");
  cout<<"myChi2 = "<<myChi2<<endl;
  Float_t maxY = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()) , mcHist->GetBinContent(mcHist->GetMaximumBin()));

  mcHist->SetLineColor(kRed);
  mcHist->SetLineWidth(2.0);
  dataHist->SetLineColor(kBlack);
  dataHist->SetLineWidth(2.0);
  dataHist->SetMarkerStyle(20);
  dataHist->SetMarkerSize(0.8);

  mcHist->GetYaxis()->SetRangeUser(0,maxY*1.2);
  mcHist->Draw("HIST");
  dataHist->Draw("E0same");

  TLatex chi2;
  chi2.SetTextSize(0.05);
  chi2.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2));

  // c1->cd(2);

  // treeIn_data->Draw("L_TAU>>l_tau_data(100,0,0.5)","SW","goff");
  // treeIn_mc->Draw("L_TAU>>l_tau_mc(100,0,0.5)",tmCut,"goff");

  // TH1D *l_tau_data = (TH1D*)gDirectory->Get("l_tau_data");
  // TH1D *l_tau_mc = (TH1D*)gDirectory->Get("l_tau_mc");

  // l_tau_data->Scale(1.0/l_tau_data->Integral());
  // l_tau_mc->Scale(1.0/l_tau_mc->Integral());

  // maxY = TMath::Max(l_tau_data->GetBinContent(l_tau_data->GetMaximumBin()) , l_tau_mc->GetBinContent(l_tau_mc->GetMaximumBin()));

  // l_tau_mc->SetLineColor(kRed);
  // l_tau_mc->SetLineWidth(2.0);
  // l_tau_data->SetLineColor(kBlack);
  // l_tau_data->SetLineWidth(2.0);
  // l_tau_data->SetMarkerStyle(20);
  // l_tau_data->SetMarkerSize(0.8);

  // l_tau_mc->GetYaxis()->SetRangeUser(0,maxY*1.2);
  // l_tau_mc->Draw("HIST");
  // l_tau_data->Draw("E0same");

  // TCanvas *c2 = new TCanvas();
  
  // treeIn_data->Draw("p_PIDp>>p_PIDp_data(100,-10,150)","SW","goff");
  // treeIn_mc->Draw("p_PIDp>>p_PIDp_mc(100,-10,150)",tmCut,"goff");

  // TH1D *p_PIDp_data = (TH1D*)gDirectory->Get("p_PIDp_data");
  // TH1D *p_PIDp_mc = (TH1D*)gDirectory->Get("p_PIDp_mc");

  // p_PIDp_data->Scale(1.0/p_PIDp_data->Integral());
  // p_PIDp_mc->Scale(1.0/p_PIDp_mc->Integral());

  // p_PIDp_mc->SetLineColor(kRed);
  // p_PIDp_mc->SetLineWidth(2.0);
  // p_PIDp_data->SetLineColor(kBlack);
  // p_PIDp_data->SetLineWidth(2.0);
  // p_PIDp_data->SetMarkerStyle(20);
  // p_PIDp_data->SetMarkerSize(0.8);
  
  // maxY = TMath::Max(p_PIDp_data->GetBinContent(p_PIDp_data->GetMaximumBin()) , p_PIDp_mc->GetBinContent(p_PIDp_mc->GetMaximumBin()));
  // // cout<<"maxY = "<<maxY<<endl;
  // //  gPad->DrawFrame(0.0,0.0,150.0,maxY*1.5);

  // p_PIDp_mc->GetYaxis()->SetRangeUser(0,maxY*1.2);
  // p_PIDp_mc->Draw("HIST");
  // p_PIDp_data->Draw("E0same");

  // TCanvas *c3 = new TCanvas();
  // c3->Divide(3,1);
  // c3->cd(1);

  // treeIn_data->Draw("Lb_ETA>>Lb_ETA_data(20,2,6)","SW","goff");
  // treeIn_mc->Draw("Lb_ETA>>Lb_ETA_mc(20,2,6)",tmCut,"goff");

  // TH1D *Lb_ETA_data = (TH1D*)gDirectory->Get("Lb_ETA_data");
  // TH1D *Lb_ETA_mc = (TH1D*)gDirectory->Get("Lb_ETA_mc");

  // Lb_ETA_data->Scale(1.0/Lb_ETA_data->Integral());
  // Lb_ETA_mc->Scale(1.0/Lb_ETA_mc->Integral());

  // Lb_ETA_mc->SetLineColor(kRed);
  // Lb_ETA_mc->SetLineWidth(2.0);
  // Lb_ETA_data->SetLineColor(kBlack);
  // Lb_ETA_data->SetLineWidth(2.0);
  // Lb_ETA_data->SetMarkerStyle(20);
  // Lb_ETA_data->SetMarkerSize(0.8);
  
  // maxY = TMath::Max(Lb_ETA_data->GetBinContent(Lb_ETA_data->GetMaximumBin()) , Lb_ETA_mc->GetBinContent(Lb_ETA_mc->GetMaximumBin()));

  // Lb_ETA_mc->GetYaxis()->SetRangeUser(0,maxY*1.2);
  // Lb_ETA_mc->Draw("HIST");
  // Lb_ETA_data->Draw("E0same");

  // c3->cd(2);

  // treeIn_data->Draw("Lb_PT>>Lb_PT_data(50,0,20000)","SW","goff");
  // treeIn_mc->Draw("Lb_PT>>Lb_PT_mc(50,0,20000)",tmCut,"goff");

  // TH1D *Lb_PT_data = (TH1D*)gDirectory->Get("Lb_PT_data");
  // TH1D *Lb_PT_mc = (TH1D*)gDirectory->Get("Lb_PT_mc");

  // Lb_PT_data->Scale(1.0/Lb_PT_data->Integral());
  // Lb_PT_mc->Scale(1.0/Lb_PT_mc->Integral());

  // Lb_PT_mc->SetLineColor(kRed);
  // Lb_PT_mc->SetLineWidth(2.0);
  // Lb_PT_data->SetLineColor(kBlack);
  // Lb_PT_data->SetLineWidth(2.0);
  // Lb_PT_data->SetMarkerStyle(20);
  // Lb_PT_data->SetMarkerSize(0.8);
  
  // maxY = TMath::Max(Lb_PT_data->GetBinContent(Lb_PT_data->GetMaximumBin()) , Lb_PT_mc->GetBinContent(Lb_PT_mc->GetMaximumBin()));
 
  // Lb_PT_mc->GetYaxis()->SetRangeUser(0,maxY*1.2);
  // Lb_PT_mc->Draw("HIST");
  // Lb_PT_data->Draw("E0same");

  // c3->cd(3);

  // treeIn_data->Draw("Lb_P>>Lb_P_data(50,0,500000)","SW","goff");
  // treeIn_mc->Draw("Lb_P>>Lb_P_mc(50,0,500000)",tmCut,"goff");

  // TH1D *Lb_P_data = (TH1D*)gDirectory->Get("Lb_P_data");
  // TH1D *Lb_P_mc = (TH1D*)gDirectory->Get("Lb_P_mc");

  // Lb_P_data->Scale(1.0/Lb_P_data->Integral());
  // Lb_P_mc->Scale(1.0/Lb_P_mc->Integral());

  // Lb_P_mc->SetLineColor(kRed);
  // Lb_P_mc->SetLineWidth(2.0);
  // Lb_P_data->SetLineColor(kBlack);
  // Lb_P_data->SetLineWidth(2.0);
  // Lb_P_data->SetMarkerStyle(20);
  // Lb_P_data->SetMarkerSize(0.8);
  
  // maxY = TMath::Max(Lb_P_data->GetBinContent(Lb_P_data->GetMaximumBin()) , Lb_P_mc->GetBinContent(Lb_P_mc->GetMaximumBin()));
 
  // Lb_P_mc->GetYaxis()->SetRangeUser(0,maxY*1.2);
  // Lb_P_mc->Draw("HIST");
  // Lb_P_data->Draw("E0same");

  // TCanvas *c4 = new TCanvas();
  // // c3->Divide(3,1);
  // // c3->cd(1);

  // treeIn_data->Draw("(Lb_ConsLb_chi2/Lb_ConsLb_nDOF)>>Lb_ConsLb_chi2_data(100,0,50)","SW","goff");
  // treeIn_mc->Draw("(Lb_ConsLb_chi2/Lb_ConsLb_nDOF)>>Lb_ConsLb_chi2_mc(100,0,50)",tmCut,"goff");

  // TH1D *Lb_ConsLb_chi2_data = (TH1D*)gDirectory->Get("Lb_ConsLb_chi2_data");
  // TH1D *Lb_ConsLb_chi2_mc = (TH1D*)gDirectory->Get("Lb_ConsLb_chi2_mc");

  // Lb_ConsLb_chi2_data->Scale(1.0/Lb_ConsLb_chi2_data->Integral());
  // Lb_ConsLb_chi2_mc->Scale(1.0/Lb_ConsLb_chi2_mc->Integral());

  // Lb_ConsLb_chi2_mc->SetLineColor(kRed);
  // Lb_ConsLb_chi2_mc->SetLineWidth(2.0);
  // Lb_ConsLb_chi2_data->SetLineColor(kBlack);
  // Lb_ConsLb_chi2_data->SetLineWidth(2.0);
  // Lb_ConsLb_chi2_data->SetMarkerStyle(20);
  // Lb_ConsLb_chi2_data->SetMarkerSize(0.8);
  
  // maxY = TMath::Max(Lb_ConsLb_chi2_data->GetBinContent(Lb_ConsLb_chi2_data->GetMaximumBin()) , Lb_ConsLb_chi2_mc->GetBinContent(Lb_ConsLb_chi2_mc->GetMaximumBin()));

  // Lb_ConsLb_chi2_mc->GetYaxis()->SetRangeUser(0,maxY*1.2);
  // Lb_ConsLb_chi2_mc->Draw("HIST");
  // Lb_ConsLb_chi2_data->Draw("E0same");

}
