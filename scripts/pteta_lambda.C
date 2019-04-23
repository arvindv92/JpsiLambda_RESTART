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
//2D reweighting in L_PT and L_ETA
void pteta_lambda(Int_t run = 1)
{
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

  gSystem->Exec(Form("cp rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt.root"
		     " rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt_lambda.root",run,run));

  TFile *filein_data,*filein_mc,*filein_wts,*filein_rw;
  TTree *treein_data,*treein_mc,*treein_rw_rec,*treein_rw_gen;

  Int_t nentries_mc_rec, nentries_mc_gen;
  Int_t nbins_leta = 20;
  Int_t nbins_lbpt = 50;
  Int_t etabin, ptbin;

  Float_t weights[50][20];
  Float_t binwidth_leta, binwidth_lbpt;
  Float_t ptetawt_rec, ptetawt_gen;

  Double_t myeta_rec, mypt_rec, mypz_rec;
  Double_t myeta_gen, mypt_gen, mypz_gen;
  Float_t letalow, letahigh, lbptlow, lbpthigh;

  TH2D *pteta_data, *pteta_mc, *wts;

  letalow = 1.5;
  letahigh = 5.5;
  lbptlow = 0;
  lbpthigh = 15000;

  binwidth_leta = (Float_t)(letahigh-letalow)/nbins_leta;
  binwidth_lbpt = (Float_t)(lbpthigh-lbptlow)/nbins_lbpt;

  filein_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/sWeightSanity/jpsilambda_LL_sanity_withsw_new.root",run),"READ");
  filein_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt.root",run),"READ");

  treein_data = (TTree*)filein_data->Get("MyTuple");
  treein_mc = (TTree*)filein_mc->Get("MyTuple");

  // nentries_mc = treein_mc->GetEntries();
  // nentries_data = treein_data->GetEntries();

  treein_data->Draw(Form("L_ETA:L_PT>>pteta_data(%d,%f,%f,%d,%f,%f)",nbins_lbpt,lbptlow,lbpthigh,nbins_leta,letalow,letahigh),"SW","goff");
  // treein_mc->Draw(Form("(-log(tan(0.5*atan(Lambda_b0_TRUEPT/Lambda_b0_TRUEP_Z)))):Lambda_b0_TRUEPT>>pteta_mc(%d,%d,%d,%d,%d,%d)",nbins_lbpt,lbptlow,lbpthigh,nbins_leta,letalow,letahigh),"","goff");
  treein_mc->Draw(Form("L_ETA:L_PT>>pteta_mc(%d,%f,%f,%d,%f,%f)",nbins_lbpt,lbptlow,lbpthigh,nbins_leta,letalow,letahigh),"ptetawt*(L_BKGCAT==0||L_BKGCAT==50)","goff");

  // treein_mc->Draw(Form("L_ETA>>leta_mc(%d,%d,%d)",nbins_leta,letalow,letahigh),"","goff");
  // treein_mc->Draw(Form("L_PT>>lbpt_mc(%d,%d,%d)",nbins_lbpt,lbptlow,lbpthigh),"","goff");

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
  //      for(Int_t j=0; j<nbins_leta; j++)
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
  filein_wts = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/ptetawt_lambda.root",run),
                           "RECREATE");
  
  wts->Write();
  filein_wts->Close();

  filein_rw = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt_lambda.root",run),
			  "UPDATE");
  treein_rw_rec = (TTree*)filein_rw->Get("MyTuple");
  treein_rw_gen = (TTree*)filein_rw->Get("MCDecayTree");

  nentries_mc_rec = treein_rw_rec->GetEntries();
  nentries_mc_gen = treein_rw_gen->GetEntries();

  treein_rw_rec->SetBranchAddress("L_PZ",&mypz_rec);
  treein_rw_rec->SetBranchAddress("L_PT",&mypt_rec);

  treein_rw_gen->SetBranchAddress("Lambda_b0_TRUEP_Z",&mypz_gen);
  treein_rw_gen->SetBranchAddress("Lambda_b0_TRUEPT",&mypt_gen);

  TBranch *wtbranch_rec = treein_rw_rec->Branch("ptetawt_L", &ptetawt_rec, "ptetawt_L/F");
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

	// etabin = TMath::CeilNint((Double_t)(myeta_rec-letalow)/binwidth_leta);
	// ptbin = TMath::CeilNint((Double_t)(mypt_rec-lbptlow)/binwidth_lbpt);
	//
	// if(myeta_rec > letahigh)
	//         etabin = nbins_leta;
	// if(mypt_rec > lbpthigh)
	//         ptbin = nbins_lbpt;
	// if(myeta_rec < letalow)
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

  TBranch *wtbranch_gen = treein_rw_gen->Branch("ptetawt_L", &ptetawt_gen, "ptetawt_L/F");
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

	// etabin = TMath::CeilNint((Double_t)(myeta_gen-letalow)/binwidth_leta);
	// ptbin = TMath::CeilNint((Double_t)(mypt_gen-lbptlow)/binwidth_lbpt);
	//
	// if(myeta_gen > letahigh)
	//         etabin = nbins_leta;
	// if(mypt_gen > lbpthigh)
	//         ptbin = nbins_lbpt;
	// if(myeta_gen < letalow)
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
