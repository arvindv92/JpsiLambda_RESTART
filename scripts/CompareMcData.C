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
//Compare sWeighted data weight generated MC

void routine(Int_t run = 1, Int_t mcType = 0,const char *varName = "",Float_t low = 0.0, Float_t high = 0.0, Int_t nBins = 50);
void CompareMcData(Int_t run = 1, Int_t mcType = 0)
{
	const char* varNameArray[26] = {"Lb_PT","Lb_P","Lb_ETA","Jpsi_PT","Jpsi_P","Jpsi_ETA","L_PT","L_P","L_ETA","p_PT","p_P","p_ETA","pi_PT","pi_P","pi_ETA","(Jpsi_P/Lb_P)","(p_P/L_P)","p_ProbNNp","pi_ProbNNpi","p_PIDp","pi_PIDK","(p_P-pi_P)/(p_P+pi_P)","L_TAU","L_FD_ORIVX","L_ENDVERTEX_Z","ntracks"};
	Float_t lowArray[26]         = {0.0,0.0,2.0,0.0,0.0,1.5,0.0,0.0,1.5,0.0,0.0,1.5,0.0,0.0,1.5,0.3,0.7,0.0,0.0,-50,-100,0.5,0.0,0.0,0,0};
	Float_t highArray[26]        = {20000,300000,6.0,20000,300000,6.5,15000,200000,5.5,10000,100000,5.5,3000,50000,5.5,1.0,1.0,1.0,1.0,150,40,1.0,0.3,1000,1000,1000 };
	Float_t nBinArray[26]        = {50,50,20,50,50,20,50,50,20,50,50,20,50,50,20,15,15,50,50,50,50,20,40,25,25,50 };

	for(Int_t i=17; i<20; i++)
	{
		cout<<"********"<<endl;
		cout<<"VARNAME = "<<varNameArray[i]<<endl;
		routine(run,mcType,varNameArray[i],lowArray[i],highArray[i],nBinArray[i]);
	}

}
void routine(Int_t run, Int_t mcType,const char *varName,Float_t low, Float_t high, Int_t nBins)
{
  cout<<"Entered routine"<<endl;
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	gStyle->SetOptStat(0);

	TFile *fileIn_data = nullptr, *fileIn_mc = nullptr;
	TTree *treeIn_data = nullptr, *treeIn_mc = nullptr, *treeIn_mcgen = nullptr;
	TString tmCut = "";

	// if(tmFlag)
	//      tmCut = "(Lb_BKGCAT==0||Lb_BKGCAT==50)";
	Double_t myChi2 = 0.0, myChi2_rw = 0.0, myChi2_corr_rw = 0.0;
	Double_t genChi2 = 0.0, genChi2_rw = 0.0, genChi2_corr_rw = 0.0;

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
	cout<<"poop1"<<endl;

	fileIn_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/sWeightSanity/jpsilambda_LL_sanity_withsw_noPID.root",run),"READ");
	treeIn_data = (TTree*)fileIn_data->Get("MyTuple");

	fileIn_mc   = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_pidgen.root",folder,run,part),"READ");
	// fileIn_mc   = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s.root",folder,run,part),"READ");

	// treeIn_mcgen = (TTree*)fileIn_mc->Get("MCTuple/MCDecayTree");
	treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");
	cout<<"poop2"<<endl;
	// treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");
	// treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/gbWeights_pid_rec.root",folder,run));
	treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/gbWeights_rec.root",folder,run));
	cout<<"poop3"<<endl;
	treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_rec.root",folder,run));
	cout<<"poop4"<<endl;
	// treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_pidgen.root",folder,run,part));
	cout<<"poop5"<<endl;
	treeIn_data->Draw(Form("%s>>dataHist(%d,%f,%f)",varName,nBins,low,high),"SW","goff");
	cout<<"poop6"<<endl;
	treeIn_mc->Draw(Form("%s>>mcHist(%d,%f,%f)",varName,nBins,low,high),"(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");
	cout<<"poop7"<<endl;
	treeIn_mc->Draw(Form("%s_corr>>mcHist_corr(%d,%f,%f)",varName,nBins,low,high),"(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");
	cout<<"poop8"<<endl;
	treeIn_mc->Draw(Form("%s_corr>>mcHist_corr_rw(%d,%f,%f)",varName,nBins,low,high),"gb_wts*wt_tau*(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");
	cout<<"poop9"<<endl;
	// treeIn_mcgen->Draw(Form("%s>>genHist_rw(%d,%f,%f)",mcvarName,nBins,low,high),"","goff");
	// treeIn_mcgen->Draw(Form("%s>>genHist(%d,%f,%f)",mcvarName,nBins,low,high),"","goff");

	TH1D *dataHist = (TH1D*)gDirectory->Get("dataHist");
	TH1D *mcHist = (TH1D*)gDirectory->Get("mcHist");
	TH1D *mcHist_corr = (TH1D*)gDirectory->Get("mcHist_corr");
	TH1D *mcHist_corr_rw = (TH1D*)gDirectory->Get("mcHist_corr_rw");
	cout<<"poop10"<<endl;
	// TH1D *genHist = (TH1D*)gDirectory->Get("genHist");
	// TH1D *genHist_rw = (TH1D*)gDirectory->Get("genHist_rw");

	dataHist->Scale(1.0/dataHist->Integral());
	mcHist->Scale(1.0/mcHist->Integral());
	mcHist_corr->Scale(1.0/mcHist_corr->Integral());
	cout<<"poop10.5"<<endl;
	mcHist_corr_rw->Scale(1.0/mcHist_corr_rw->Integral());
	cout<<"poop11"<<endl;
	// genHist->Scale(1.0/genHist->Integral());
	// genHist_rw->Scale(1.0/genHist_rw->Integral());

	for(Int_t i=1; i<=dataHist->GetNbinsX(); i++)
	{
		if(dataHist->GetBinContent(i) < 0)
			dataHist->SetBinContent(i,0);
	}

	myChi2    = mcHist->Chi2Test(dataHist,"WW CHI2/NDF");
	myChi2_rw = mcHist_corr->Chi2Test(dataHist,"WW CHI2/NDF");
	myChi2_corr_rw = mcHist_corr_rw->Chi2Test(dataHist,"WW CHI2/NDF");
	// genChi2    = genHist->Chi2Test(dataHist,"WW CHI2/NDF");
	// genChi2_rw = genHist_rw->Chi2Test(dataHist,"WW CHI2/NDF");
	cout<<"poop12"<<endl;
	cout<<"myChi2 = "<<myChi2<<endl;
	cout<<"myChi2_rw = "<<myChi2_rw<<endl;
	cout<<"myChi2_corr_rw = "<<myChi2_corr_rw<<endl;

	// cout<<"gen myChi2 = "<<genChi2<<endl;
	// cout<<"gen myChi2_rw = "<<genChi2_rw<<endl;

	Float_t maxY = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist->GetBinContent(mcHist->GetMaximumBin()));
	Float_t maxY1 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist_corr->GetBinContent(mcHist_corr->GetMaximumBin()));
	Float_t maxY2 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist_corr_rw->GetBinContent(mcHist_corr_rw->GetMaximumBin()));
	// Float_t maxY2 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), genHist->GetBinContent(genHist->GetMaximumBin()));
	// Float_t maxY3 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), genHist_rw->GetBinContent(genHist_rw->GetMaximumBin()));

	mcHist->SetLineColor(kRed);
	mcHist->SetLineWidth(2.0);

	mcHist_corr->SetLineColor(kRed);
	mcHist_corr->SetLineWidth(2.0);

	mcHist_corr_rw->SetLineColor(kRed);
	mcHist_corr_rw->SetLineWidth(2.0);

	// genHist->SetLineColor(kRed);
	// genHist->SetLineWidth(2.0);

	// genHist_rw->SetLineColor(kRed);
	// genHist_rw->SetLineWidth(2.0);

	dataHist->SetLineColor(kBlack);
	dataHist->SetLineWidth(2.0);
	dataHist->SetMarkerStyle(20);
	dataHist->SetMarkerSize(0.8);

	mcHist->GetYaxis()->SetRangeUser(0,maxY*1.2);
	mcHist_corr->GetYaxis()->SetRangeUser(0,maxY1*1.2);
	mcHist_corr_rw->GetYaxis()->SetRangeUser(0,maxY2*1.2);
	// genHist->GetYaxis()->SetRangeUser(0,maxY2*1.2);
	// genHist_rw->GetYaxis()->SetRangeUser(0,maxY3*1.2);

	TCanvas *c1 = new TCanvas();
	c1->Divide(3,1);
	c1->cd(1);

	mcHist->Draw("HIST");
	dataHist->Draw("E0same");

	TLatex chi2;
	chi2.SetTextSize(0.05);
	chi2.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2));

	c1->cd(2);

	mcHist_corr->Draw("HIST");
	dataHist->Draw("E0same");

	TLatex chi2_rw;
	chi2_rw.SetTextSize(0.05);
	chi2_rw.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2_rw));

	c1->cd(3);

	mcHist_corr_rw->Draw("HIST");
	dataHist->Draw("E0same");

	TLatex chi2_corr_rw;
	chi2_corr_rw.SetTextSize(0.05);
	chi2_corr_rw.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2_corr_rw));

}
