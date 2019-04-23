#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
using namespace std;

//This code does the calculation of weights, not the application
//The weights are calculated by comparing sWeighted data (post Sanity) and raw generator JpsiLambda MC (before any selections)
//2D reweighting in Lb_PT and Lb_ETA
void pteta(Int_t run = 1)
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	gSystem->Exec(Form("cp rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda.root"
	                   " rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt.root",run,run));

	TFile *filein_data,*filein_mc,*filein_wts,*filein_rw;
	TTree *treein_data,*treein_mc,*treein_rw_rec,*treein_rw_gen;

	Int_t nentries_mc_rec, nentries_mc_gen;
	Int_t nbins_lbeta = 20;
	Int_t nbins_lbpt = 50;
	Int_t etabin, ptbin;

	Float_t weights[50][20];
	Float_t binwidth_lbeta, binwidth_lbpt;
	Float_t ptetawt_rec, ptetawt_gen;

	Double_t myeta_rec, mypt_rec, mypz_rec;
	Double_t myeta_gen, mypt_gen, mypz_gen;
	Int_t lbetalow, lbetahigh, lbptlow, lbpthigh;

	TH2D *pteta_data, *pteta_mc, *wts;

	lbetalow = 2;
	lbetahigh = 6;
	lbptlow = 0;
	lbpthigh = 20000;

	binwidth_lbeta = (Float_t)(lbetahigh-lbetalow)/nbins_lbeta;
	binwidth_lbpt = (Float_t)(lbpthigh-lbptlow)/nbins_lbpt;

	filein_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/sWeightSanity/jpsilambda_LL_sanity_withsw_new.root",run),"READ");
	filein_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda.root",run),"READ");

	treein_data = (TTree*)filein_data->Get("MyTuple");
	treein_mc = (TTree*)filein_mc->Get("Lb2JpsiLTree/MyTuple");

	// nentries_mc = treein_mc->GetEntries();
	// nentries_data = treein_data->GetEntries();

	treein_data->Draw(Form("Lb_ETA:Lb_PT>>pteta_data(%d,%d,%d,%d,%d,%d)",nbins_lbpt,lbptlow,lbpthigh,nbins_lbeta,lbetalow,lbetahigh),"SW","goff");
	// treein_mc->Draw(Form("(-log(tan(0.5*atan(Lambda_b0_TRUEPT/Lambda_b0_TRUEP_Z)))):Lambda_b0_TRUEPT>>pteta_mc(%d,%d,%d,%d,%d,%d)",nbins_lbpt,lbptlow,lbpthigh,nbins_lbeta,lbetalow,lbetahigh),"","goff");
	treein_mc->Draw(Form("(-log(tan(0.5*atan(Lb_PT/Lb_PZ)))):Lb_PT>>pteta_mc(%d,%d,%d,%d,%d,%d)",nbins_lbpt,lbptlow,lbpthigh,nbins_lbeta,lbetalow,lbetahigh),"(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");

	// treein_mc->Draw(Form("Lb_ETA>>lbeta_mc(%d,%d,%d)",nbins_lbeta,lbetalow,lbetahigh),"","goff");
	// treein_mc->Draw(Form("Lb_PT>>lbpt_mc(%d,%d,%d)",nbins_lbpt,lbptlow,lbpthigh),"","goff");

	pteta_data = (TH2D*)gDirectory->Get("pteta_data");
	pteta_mc = (TH2D*)gDirectory->Get("pteta_mc");

	// pteta_data->Smooth(1);//NOTE SMOOTHING HERE
	// pteta_mc->Smooth(1); //NOTE SMOOTHING HERE

	cout<<"Data integral = "<<pteta_data->Integral()<<endl;
	cout<<"MC integral = "<<pteta_mc->Integral()<<endl;

	pteta_data->Scale(1.0/pteta_data->Integral());
	pteta_mc->Scale(1.0/pteta_mc->Integral());

	wts = (TH2D*)pteta_data->Clone("wts");
	wts->Divide(pteta_mc);

	// for(Int_t i=0; i<nbins_lbpt; i++)
	// {
	//      for(Int_t j=0; j<nbins_lbeta; j++)
	//      {
	//              binc_data=pteta_data->GetBinContent(i+1,j+1);
	//              binc_mc=pteta_mc->GetBinContent(i+1,j+1);
	//              // if(binc_data!=0 && binc_mc!=0)
	//              //      cout<<"data = "<<binc_data<<" mc = "<<binc_mc<<endl;
	//
	//              if(binc_mc!=0 && binc_data > 0)
	//                      weights[i][j] = (Float_t)binc_data/binc_mc;
	//
	//              else if(binc_mc!=0 && binc_data <= 0)
	//                      weights[i][j] = 0;
	//
	//              else if(binc_mc==0)
	//                      weights[i][j] = 1;
	//
	//              if(weights[i][j] > 5)
	//                      cout<<"REDFLAG data = "<<binc_data<<" mc = "<<binc_mc<<" i = "<<i<<" j = "<<j<<endl;
	//
	//              wts->SetBinContent(i+1,j+1,weights[i][j]);
	//      }
	// }
	cout<<"poop2"<<endl;
	// filein_wts = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/ptetawt.root",run),
	//                          "RECREATE");
	//
	// wts->Write();
	// filein_wts->Close();

	filein_rw = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt.root",run),
	                        "UPDATE");
	treein_rw_rec = (TTree*)filein_rw->Get("Lb2JpsiLTree/MyTuple");
	treein_rw_gen = (TTree*)filein_rw->Get("MCTuple/MCDecayTree");

	nentries_mc_rec = treein_rw_rec->GetEntries();
	nentries_mc_gen = treein_rw_gen->GetEntries();

	treein_rw_rec->SetBranchAddress("Lb_PZ",&mypz_rec);
	treein_rw_rec->SetBranchAddress("Lb_PT",&mypt_rec);

	treein_rw_gen->SetBranchAddress("Lambda_b0_TRUEP_Z",&mypz_gen);
	treein_rw_gen->SetBranchAddress("Lambda_b0_TRUEPT",&mypt_gen);

	TBranch *wtbranch_rec = treein_rw_rec->Branch("ptetawt", &ptetawt_rec, "ptetawt/F");
	cout<<"Looping over reconstructed tree"<<endl;
	for(Int_t i=0; i<nentries_mc_rec; i++) {
		if(i%10000 == 0)
			cout<<i<<endl;

		treein_rw_rec->GetEntry(i);
		if(mypz_rec > 0 && mypt_rec > 0)
		{
			myeta_rec = -log(tan(0.5*atan(mypt_rec/mypz_rec)));

			Int_t globBin = wts->FindBin(mypt_rec,myeta_rec);
			ptetawt_rec   = wts->GetBinContent(globBin);

			// etabin = TMath::CeilNint((Double_t)(myeta_rec-lbetalow)/binwidth_lbeta);
			// ptbin = TMath::CeilNint((Double_t)(mypt_rec-lbptlow)/binwidth_lbpt);
			//
			// if(myeta_rec > lbetahigh)
			//         etabin = nbins_lbeta;
			// if(mypt_rec > lbpthigh)
			//         ptbin = nbins_lbpt;
			// if(myeta_rec < lbetalow)
			//         etabin = 1;
			// if(mypt_rec < lbptlow)
			//         ptbin = 1;
			//
			// ptetawt_rec = weights[ptbin-1][etabin-1];
		}
		else
		{
			ptetawt_rec = 1.0;
		}
		wtbranch_rec->Fill();
	}
	treein_rw_rec->Write("",TObject::kOverwrite);

	TBranch *wtbranch_gen = treein_rw_gen->Branch("ptetawt", &ptetawt_gen, "ptetawt/F");
	cout<<"Looping over generated tree"<<endl;
	for(Int_t i=0; i<nentries_mc_gen; i++) {
		if(i%100000 == 0)
			cout<<i<<endl;

		treein_rw_gen->GetEntry(i);

		if(mypz_gen > 0 && mypt_gen > 0)
		{
			myeta_gen = -log(tan(0.5*atan(mypt_gen/mypz_gen)));

			Int_t globBin = wts->FindBin(mypt_gen,myeta_gen);
			ptetawt_gen   = wts->GetBinContent(globBin);

			// etabin = TMath::CeilNint((Double_t)(myeta_gen-lbetalow)/binwidth_lbeta);
			// ptbin = TMath::CeilNint((Double_t)(mypt_gen-lbptlow)/binwidth_lbpt);
			//
			// if(myeta_gen > lbetahigh)
			//         etabin = nbins_lbeta;
			// if(mypt_gen > lbpthigh)
			//         ptbin = nbins_lbpt;
			// if(myeta_gen < lbetalow)
			//         etabin = 1;
			// if(mypt_gen < lbptlow)
			//         ptbin = 1;
			//
			// ptetawt_gen = weights[ptbin-1][etabin-1];
		}
		else
		{
			ptetawt_gen = 1.0;
		}
		wtbranch_gen->Fill();
	}
	treein_rw_gen->Write("",TObject::kOverwrite);

	filein_rw->Close();
	cout<<"DONE"<<endl;
	//
	// gDirectory->Add(pteta_data);
	// gDirectory->Add(pteta_mc);
}
