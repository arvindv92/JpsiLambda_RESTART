#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TStyle.h"
#include "TString.h"
#include "TLatex.h"
#include "TROOT.h"
#include "iostream"

using namespace std;
//Compare sWeighted data weight generated MC

TCanvas* routine(Int_t run = 1, Int_t mcType = 0,const char *varName = "",Float_t low = 0.0,
                 Float_t high = 0.0, Int_t nBins = 50, TString unit = "", TString option = "kinematic");

void addGraphics(TH1F *h, TString Xtitle, TString Ytitle, int iCol){
	h->SetXTitle(Xtitle);
	//h->SetFillColor(30);
	//int bw = h->GetBinWidth(1);
	h->SetYTitle(Ytitle);
	h->SetStats(kFALSE);
	h->SetMinimum(0.1);
	h->SetMaximum(1.2*h->GetMaximum());
	h->SetMinimum(0);
	// h->SetTitleSize(0.1);
	h->SetLineColor(iCol);
	h->SetMarkerColor(iCol);
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

void CompareMcData(Int_t run = 1, Int_t mcType = 0, TString option = "kinematic")
//option = "kinematic" or "pid" or "finalBDT"
{
	gROOT->ProcessLine(".x lhcbStyle.C");
	if(option == "kinematic")
	{
		const char* varNameArray[15] = {"Lb_PT","Lb_P","Lb_ETA","Jpsi_PT","Jpsi_P",
			                        "Jpsi_ETA","L_PT","L_P","L_ETA","p_PT","p_P",
			                        "p_ETA","pi_PT","pi_P","pi_ETA"};
		Float_t lowArray[15]         = {0.0,0.0,2.0,0.0,0.0,1.5,0.0,0.0,1.5,0.0,0.0,1.5,0.0,0.0,1.5};
		Float_t highArray[15]        = {20000,300000,6.0,20000,300000,6.5,15000,200000,5.5,10000,100000,5.5,3000,50000,5.5};
		Float_t nBinArray[15]        = {50,50,20,50,50,20,50,50,20,50,50,20,50,50,20};
		TString units[15]            = {"MeV/#it{c}^{2}","MeV/#it{c}^{2}","","MeV/#it{c}^{2}",
			                        "MeV/#it{c}^{2}","","MeV/#it{c}^{2}","MeV/#it{c}^{2}","",
			                        "MeV/#it{c}^{2}","MeV/#it{c}^{2}","","MeV/#it{c}^{2}","MeV/#it{c}^{2}",""};

		TCanvas *canArray[15];

		for(Int_t i=0; i<15; i++)
		{
			cout<<"********"<<endl;
			cout<<"VARNAME = "<<varNameArray[i]<<endl;
			canArray[i] = routine(run,mcType,varNameArray[i],lowArray[i],highArray[i],nBinArray[i],units[i],option);
		}
		gROOT->SetBatch(kFALSE);
		TCanvas *c1 = new TCanvas("c1","",1200,1200);
		c1->Divide(2,3);
		TCanvas *c2 = new TCanvas("c2","",1200,1200);
		c2->Divide(2,3);
		TCanvas *c3 = new TCanvas("c3","",1200,800);
		c3->Divide(2,2);

		for(Int_t i=0; i<6; i++)
		{
			c1->cd(i+1);
			canArray[i]->DrawClonePad();
		}
		for(Int_t i=6; i<12; i++)
		{
			c2->cd(i-5);
			canArray[i]->DrawClonePad();
		}
		for(Int_t i=12; i<15; i++)
		{
			c3->cd(i-11);
			canArray[i]->DrawClonePad();
		}

		gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
		c1->SaveAs(Form("plots/ANA/MC_Data_kinematic_run%d_page1.pdf",run));
		c2->SaveAs(Form("plots/ANA/MC_Data_kinematic_run%d_page2.pdf",run));
		c3->SaveAs(Form("plots/ANA/MC_Data_kinematic_run%d_page3.pdf",run));

	}
	else if(option == "pid")
	{
		const char* varNameArray[2] = {"p_ProbNNp","p_PIDp"};
		Float_t lowArray[2]         = {0.5,-20.0};
		Float_t highArray[2]        = {1.0,100.0};
		Float_t nBinArray[2]        = {25,50};
		TString units[2]            = {"",""};

		TCanvas *canArray[2];
		for(Int_t i=0; i<2; i++)
		{
			cout<<"********"<<endl;
			cout<<"VARNAME = "<<varNameArray[i]<<endl;
			canArray[i] = routine(run,mcType,varNameArray[i],lowArray[i],highArray[i],nBinArray[i],units[i],option);
		}
		gROOT->SetBatch(kFALSE);
		TCanvas *c1 = new TCanvas("c1","",1200,400);
		c1->Divide(2,1);

		for(Int_t i=0; i<2; i++)
		{
			c1->cd(i+1);
			canArray[i]->DrawClonePad();
		}

		gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
		c1->SaveAs(Form("plots/ANA/MC_Data_pid_run%d_page1.pdf",run));
	}
	else if(option == "finalBDT")
	{
		const char* varNameArray[20] = {"log10((Lb_ConsLb_chi2))","log10(Lb_MINIPCHI2)","log10(acos(Lb_DIRA_OWNPV))","log10(Lb_FD_OWNPV)",
			                        "log10(Jpsi_MINIPCHI2)","log10(Jpsi_M)","log10(L_FDCHI2_ORIVX)","log10(acos(L_DIRA_ORIVX))",
			                        "log10(L_FD_ORIVX)","log10(acos(L_DIRA_OWNPV))","L_dm","log10(L_MINIPCHI2)",
			                        "log10(p_MINIPCHI2)","p_ProbNNghost","log10(p_PT)","p_ProbNNp","pi_ProbNNghost",
			                        "log10(pi_MINIPCHI2)","log10(pi_PT)","BDTkMin_v0"};
		Float_t lowArray[20]         = {-1,-4,-5,0,-2,3.48,-2,-6,-1,-5,0,-2,0,0,2,0,0,1,1.0,-1};
		Float_t highArray[20]        = {3,4,0,2.5,5,3.5,6,-1,4,0,15,6,5,0.3,4,1,0.3,4.5,3.5,1};
		Float_t nBinArray[20]        = {50,50,50,50,50,50,50,50,50,50,50,50,50,30,50,50,30,50,50,50};
		TString units[20]            = {"","","","mm","","MeV/#it{c}^{2}","","","mm","","MeV/#it{c}^{2}","","","","MeV/#it{c}^{2}","","","","MeV/#it{c}^{2}",""};

		TCanvas *canArray[20];
		for(Int_t i=0; i<20; i++)
		{
			cout<<"********"<<endl;
			cout<<"VARNAME = "<<varNameArray[i]<<endl;
			canArray[i] = routine(run,mcType,varNameArray[i],lowArray[i],highArray[i],nBinArray[i],units[i],option);
		}
		gROOT->SetBatch(kFALSE);

		TCanvas *c1 = new TCanvas("c1","",1200,1200);
		c1->Divide(2,3);
		TCanvas *c2 = new TCanvas("c2","",1200,1200);
		c2->Divide(2,3);
		TCanvas *c3 = new TCanvas("c3","",1200,1200);
		c3->Divide(2,3);
		TCanvas *c4 = new TCanvas("c4","",1200,400);
		c4->Divide(2,1);

		for(Int_t i=0; i<6; i++)
		{
			c1->cd(i+1);
			canArray[i]->DrawClonePad();
		}
		for(Int_t i=6; i<12; i++)
		{
			c2->cd(i-5);
			canArray[i]->DrawClonePad();
		}
		for(Int_t i=12; i<18; i++)
		{
			c3->cd(i-11);
			canArray[i]->DrawClonePad();
		}
		for(Int_t i=18; i<20; i++)
		{
			c4->cd(i-17);
			canArray[i]->DrawClonePad();
		}

		gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
		c1->SaveAs(Form("plots/ANA/MC_Data_finalBDT_run%d_page1.pdf",run));
		c2->SaveAs(Form("plots/ANA/MC_Data_finalBDT_run%d_page2.pdf",run));
		c3->SaveAs(Form("plots/ANA/MC_Data_finalBDT_run%d_page3.pdf",run));
		c4->SaveAs(Form("plots/ANA/MC_Data_finalBDT_run%d_page4.pdf",run));
	}
}
TCanvas* routine(Int_t run, Int_t mcType,const char *varName,Float_t low, Float_t high, Int_t nBins, TString unit, TString option)
{
	gROOT->SetBatch(kTRUE);
	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
	// gStyle->SetOptStat(0);

	TFile *fileIn_data = nullptr, *fileIn_mc = nullptr;
	TTree *treeIn_data = nullptr, *treeIn_mc = nullptr, *treeIn_mcgen = nullptr;
	TString tmCut = "";

	// if(tmFlag)
	//      tmCut = "(Lb_BKGCAT==0||Lb_BKGCAT==50)";
	Double_t myChi2 = 0.0, myChi2_rw = 0.0;
	Double_t myChi2_corr_rw = 0.0, myChi2_corr = 0.0;
	Double_t myChi2_uncorr = 0.0;
	Double_t genChi2 = 0.0, genChi2_rw = 0.0, genChi2_corr_rw = 0.0;

	const char *folder = "", *part = "";
	const char *gbWt = "";
	if(run == 1)
	{
		gbWt = "GB_WT";
	}
	else if(run == 2)
	{
		gbWt = "GB_WT";
	}
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
	default:
	{
		cout<<"$$$ MC Type doesn't match any of the allowed cases. Exiting! $$$"<<endl;
		exit(1);
	}
	}

	if(option == "kinematic" or option == "pid")
	{
		fileIn_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/sWeightSanity/jpsilambda_LL_sanity_withsw.root",run),"READ");
		treeIn_data = (TTree*)fileIn_data->Get("MyTuple");
		fileIn_mc   = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_sanity_LL.root",folder,run,part),"READ");
		treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");
	}
	else if(option == "finalBDT")
	{
		fileIn_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_LL_withsw_nonZeroTracks.root",run),"READ");
		treeIn_data = (TTree*)fileIn_data->Get("MyTuple");

		fileIn_mc   = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_cutoutks_LL_nonZeroTracks.root",folder,run,part),"READ");
		treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

		if(run == 1)
		{
			treeIn_data->AddFriend("MyTuple","rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LLsig_iso2_v0.root");
			treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run1/%s_LL_iso2_v0.root",folder,part));
		}
		else if(run == 2)
		{
			treeIn_data->AddFriend("MyTuple","rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_LLsig_iso1_v0.root");
			treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run2/%s_LL_iso1_v0.root",folder,part));
		}
	}

	// treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/gbWeights_rec.root",folder,run));
	// treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/tauWeights_rec.root",folder,run));

	treeIn_data->Draw(Form("%s>>dataHist(%d,%f,%f)",varName,nBins,low,high),"SW"," goff");

	if(!strncmp(varName,"p_PIDp",6) || !strncmp(varName,"p_ProbNNp",9))
	{
		treeIn_mc->Draw(Form("%s>>mcHist_uncorr(%d,%f,%f)",varName,nBins,low,high),"(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");
		treeIn_mc->Draw(Form("%s_corr>>mcHist(%d,%f,%f)",varName,nBins,low,high),"(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");
		treeIn_mc->Draw(Form("%s_corr>>mcHist_rw(%d,%f,%f)",varName,nBins,low,high),Form("%s*wt_tau*(Lb_BKGCAT==0||Lb_BKGCAT==50)",gbWt),"goff");
	}
	else
	{
		treeIn_mc->Draw(Form("%s>>mcHist(%d,%f,%f)",varName,nBins,low,high),"(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");
		treeIn_mc->Draw(Form("%s>>mcHist_rw(%d,%f,%f)",varName,nBins,low,high),Form("%s*wt_tau*(Lb_BKGCAT==0||Lb_BKGCAT==50)",gbWt),"goff");
	}

	TH1F *dataHist    = (TH1F*)gDirectory->Get("dataHist");
	TH1F *mcHist      = (TH1F*)gDirectory->Get("mcHist");
	TH1F *mcHist_rw   = (TH1F*)gDirectory->Get("mcHist_rw");
	TH1F *mcHist_uncorr = nullptr;

	if(!strncmp(varName,"p_PIDp",6) || !strncmp(varName,"p_ProbNNp",9))
	{
		mcHist_uncorr = (TH1F*)gDirectory->Get("mcHist_uncorr");
		mcHist_uncorr->Scale(1.0/mcHist_uncorr->Integral());
		myChi2_uncorr = mcHist_uncorr->Chi2Test(dataHist,"WW CHI2/NDF");
		cout<<"myChi2_uncorr = "<<myChi2_uncorr<<endl;
	}

	dataHist->Scale(1.0/dataHist->Integral());
	mcHist->Scale(1.0/mcHist->Integral());
	mcHist_rw->Scale(1.0/mcHist_rw->Integral());

	for(Int_t i=1; i<=dataHist->GetNbinsX(); i++)
	{
		if(dataHist->GetBinContent(i) < 0)
			dataHist->SetBinContent(i,0);
	}

	myChi2    = mcHist->Chi2Test(dataHist,"WW CHI2/NDF");
	myChi2_rw = mcHist_rw->Chi2Test(dataHist,"WW CHI2/NDF");
	// myChi2_corr = mcHist_corr->Chi2Test(dataHist,"WW CHI2/NDF");
	// myChi2_corr_rw = mcHist_corr_rw->Chi2Test(dataHist,"WW CHI2/NDF");
	// genChi2    = genHist->Chi2Test(dataHist,"WW CHI2/NDF");
	// genChi2_rw = genHist_rw->Chi2Test(dataHist,"WW CHI2/NDF");

	cout<<"myChi2    = "<<myChi2<<endl;
	cout<<"myChi2_rw = "<<myChi2_rw<<endl;
	// cout<<"myChi2_corr = "<<myChi2_corr<<endl;
	// cout<<"myChi2_corr_rw = "<<myChi2_corr_rw<<endl;

	// cout<<"gen myChi2 = "<<genChi2<<endl;
	// cout<<"gen myChi2_rw = "<<genChi2_rw<<endl;

	// Float_t maxY = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist->GetBinContent(mcHist->GetMaximumBin()));
	// Float_t maxY1 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist_rw->GetBinContent(mcHist_rw->GetMaximumBin()));
	// // Float_t maxY1 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist_corr->GetBinContent(mcHist_corr->GetMaximumBin()));
	// // Float_t maxY2 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), mcHist_corr_rw->GetBinContent(mcHist_corr_rw->GetMaximumBin()));
	// // Float_t maxY2 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), genHist->GetBinContent(genHist->GetMaximumBin()));
	// // Float_t maxY3 = TMath::Max(dataHist->GetBinContent(dataHist->GetMaximumBin()), genHist_rw->GetBinContent(genHist_rw->GetMaximumBin()));
	//
	// mcHist->SetLineColor(kRed);
	// mcHist->SetLineWidth(2.0);
	//
	// mcHist_rw->SetLineColor(kRed);
	// mcHist_rw->SetLineWidth(2.0);
	// // mcHist_corr->SetLineColor(kRed);
	// // mcHist_corr->SetLineWidth(2.0);
	//
	// mcHist_corr_rw->SetLineColor(kRed);
	// mcHist_corr_rw->SetLineWidth(2.0);
	//
	// // genHist->SetLineColor(kRed);
	// // genHist->SetLineWidth(2.0);
	//
	// // genHist_rw->SetLineColor(kRed);
	// // genHist_rw->SetLineWidth(2.0);
	//
	// dataHist->SetLineColor(kBlack);
	// dataHist->SetLineWidth(2.0);
	// dataHist->SetMarkerStyle(20);
	// dataHist->SetMarkerSize(0.8);

	TString Xtit;
	if(unit!="")
		Xtit = TString(varName)+"["+unit+"]";
	else
		Xtit = TString(varName);

	if(!strncmp(varName,"BDTkMin_v0",10))
		Xtit = "isolation BDT";

	addGraphics(dataHist,Xtit,"Candidates(normalized)",1);
	addGraphics(mcHist,Xtit,"Candidates(normalized)",4);
	addGraphics(mcHist_rw,Xtit,"Candidates(normalized)",2);
	if(!strncmp(varName,"p_PIDp",6) || !strncmp(varName,"p_ProbNNp",9))
	{
		addGraphics(mcHist_uncorr,Xtit,"Candidates(normalized)",8);
	}

	TCanvas *c1 = new TCanvas(varName,"",600,400);

	if(!strncmp(varName,"p_PIDp",6) || !strncmp(varName,"p_ProbNNp",9))
	{
		Double_t max1 = mcHist_uncorr->GetMaximum();
		Double_t max2 = dataHist->GetMaximum();
		Double_t max3 = mcHist->GetMaximum();
		Double_t max4 = mcHist_rw->GetMaximum();

		Double_t max5 = TMath::Max(max1,max2);
		Double_t max6 = TMath::Max(max5,max3);
		Double_t max7 = TMath::Max(max6,max4);

		mcHist_uncorr->SetMaximum(max7*1.05);

		mcHist_uncorr->Draw("HIST");
		dataHist->Draw("E0same");
		mcHist->Draw("HISTsame");
		mcHist_rw->Draw("HISTsame");
	}
	else
	{
		Double_t max2 = dataHist->GetMaximum();
		Double_t max3 = mcHist->GetMaximum();
		Double_t max4 = mcHist_rw->GetMaximum();

		Double_t max5 = TMath::Max(max2,max3);
		Double_t max6 = TMath::Max(max5,max4);

		dataHist->SetMaximum(max6*1.05);

		dataHist->Draw("E0");
		mcHist->Draw("HISTsame");
		mcHist_rw->Draw("HISTsame");
	}

	// if(!strncmp(varName,"p_PIDp",6) || !strncmp(varName,"p_ProbNNp",9))
	// {
	//      TLatex chi2_uncorr;
	//      chi2_uncorr.SetTextSize(0.06);
	//      chi2_uncorr.DrawLatexNDC(.6,.85,Form("Orig. #chi^{2}/ndf = %.1f",myChi2_uncorr));
	//
	//      TLatex chi2;
	//      chi2.SetTextSize(0.06);
	//      chi2.DrawLatexNDC(.6,.75,Form("PIDCorr. #chi^{2}/ndf = %.1f",myChi2));
	//
	//      TLatex chi2_rw;
	//      chi2_rw.SetTextSize(0.06);
	//      chi2_rw.DrawLatexNDC(.6,.65,Form("Kin RW. #chi^{2}/ndf = %.1f",myChi2_rw));
	// }
	// else
	// {
	//      TLatex chi2;
	//      chi2.SetTextSize(0.06);
	//      chi2.DrawLatexNDC(.65,.85,Form("Orig. #chi^{2}/ndf = %.1f",myChi2));
	//
	//      TLatex chi2_rw;
	//      chi2_rw.SetTextSize(0.06);
	//      chi2_rw.DrawLatexNDC(.65,.75,Form("RW #chi^{2}/ndf = %.1f",myChi2_rw));
	// }
	return c1;

	gROOT->SetBatch(kFALSE);

	// c1->cd(2);

	// mcHist_rw->Draw("HIST");
	// dataHist->Draw("E0same");

	// TLatex chi2_rw;
	// chi2_rw.SetTextSize(0.05);
	// chi2_rw.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2_rw));

	// mcHist_corr->Draw("HIST");
	// dataHist->Draw("E0same");
	//
	// TLatex chi2_corr;
	// chi2_corr.SetTextSize(0.05);
	// chi2_corr.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2_corr));
	//
	// c1->cd(3);
	//
	// mcHist_corr_rw->Draw("HIST");
	// dataHist->Draw("E0same");
	//
	// TLatex chi2_corr_rw;
	// chi2_corr_rw.SetTextSize(0.05);
	// chi2_corr_rw.DrawLatexNDC(.6,.85,Form("#chi^{2}/ndf = %.3f",myChi2_corr_rw));

}
