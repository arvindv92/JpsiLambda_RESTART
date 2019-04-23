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
//Reweighting in Jpsi_P/Lb_P. Done ON TOP of pteta reweighting.
//pjpsi in this code is actually Jpsi_P/Lb_P
void pjpsi(Int_t run = 1)
{
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

  gSystem->Exec(Form("cp rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt.root"
		     " rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_pjpsiwt.root",run,run));

  TFile *filein_data,*filein_mc,*filein_wts,*filein_rw;
  TTree *treein_data,*treein_mc,*treein_rw_rec,*treein_rw_gen;

  Int_t nentries_mc_rec, nentries_mc_gen;
  Int_t nbins_pjpsi = 12;

  Float_t binwidth_pjpsi;
  Float_t pjpsiwt_rec, pjpsiwt_gen;

  Double_t myjpsipt_gen, myjpsipz_gen;
  Double_t mylbpt_gen, mylbpz_gen;
  Double_t myjpsip_rec, myjpsip_gen;
  Double_t mylbp_rec, mylbp_gen;
  Double_t mypjpsi_rec, mypjpsi_gen;

  Int_t pjpsilow, pjpsihigh;

  TH1D *pjpsi_data, *pjpsi_mc, *wts;

  bool interpolate = false;

  pjpsilow = 0;
  pjpsihigh = 1;

  binwidth_pjpsi = (Float_t)(pjpsihigh-pjpsilow)/nbins_pjpsi;

  filein_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/sWeightSanity/jpsilambda_LL_sanity_withsw_new.root",run),"READ");
  filein_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_ptetawt.root",run),"READ");

  treein_data = (TTree*)filein_data->Get("MyTuple");
  treein_mc = (TTree*)filein_mc->Get("MyTuple");

  // nentries_mc = treein_mc->GetEntries();
  // nentries_data = treein_data->GetEntries();

  treein_data->Draw(Form("(Jpsi_P/Lb_P)>>pjpsi_data(%d,%d,%d)",nbins_pjpsi,pjpsilow,pjpsihigh),"SW","goff");
  treein_mc->Draw(Form("(Jpsi_P/Lb_P)>>pjpsi_mc(%d,%d,%d)",nbins_pjpsi,pjpsilow,pjpsihigh),"ptetawt*(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");

  pjpsi_data = (TH1D*)gDirectory->Get("pjpsi_data");
  pjpsi_mc = (TH1D*)gDirectory->Get("pjpsi_mc");

  // pjpsi_data->Smooth(1);//NOTE SMOOTHING HAPPENING HERE
  // pjpsi_mc->Smooth(1);//NOTE SMOOTHING HAPPENING HERE

  wts = (TH1D*)pjpsi_data->Clone("wts");
  wts->Divide(pjpsi_mc);

  cout<<"Data integral = "<<pjpsi_data->Integral()<<endl;
  cout<<"MC integral = "<<pjpsi_mc->Integral()<<endl;

  pjpsi_data->Scale(1.0/pjpsi_data->Integral());
  pjpsi_mc->Scale(1.0/pjpsi_mc->Integral());

  // for(Int_t j=0; j<nbins_pjpsi; j++)
  // {
  //      binc_data=pjpsi_data->GetBinContent(j+1);
  //      binc_mc=pjpsi_mc->GetBinContent(j+1);
  //
  //      if(binc_mc!=0 && binc_data > 0)
  //              weights[j] = (Float_t)binc_data/binc_mc;
  //
  //      else if(binc_mc!=0 && binc_data <= 0)
  //              weights[j] = 0;
  //
  //      else if(binc_mc==0)
  //              weights[j] = 1;
  //
  //      if(weights[j] > 5)
  //              cout<<"REDFLAG data = "<<binc_data<<" mc = "<<binc_mc<<" j = "<<j<<endl;
  //
  //      wts->SetBinContent(j+1,weights[j]);
  // }
  cout<<"poop2"<<endl;
  filein_wts = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/pjpsiwt.root",run),
			   "RECREATE");

  wts->Write();
  filein_wts->Close();

  filein_rw = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_pjpsiwt.root",run),
			  "UPDATE");
  treein_rw_rec = (TTree*)filein_rw->Get("MyTuple");
  treein_rw_gen = (TTree*)filein_rw->Get("MCDecayTree");

  nentries_mc_rec = treein_rw_rec->GetEntries();
  nentries_mc_gen = treein_rw_gen->GetEntries();

  treein_rw_rec->SetBranchAddress("Lb_P",&mylbp_rec);
  treein_rw_rec->SetBranchAddress("Jpsi_P",&myjpsip_rec);

  treein_rw_gen->SetBranchAddress("Lambda_b0_TRUEP_Z",&mylbpz_gen);
  treein_rw_gen->SetBranchAddress("Lambda_b0_TRUEPT",&mylbpt_gen);
  treein_rw_gen->SetBranchAddress("J_psi_1S_TRUEP_Z",&myjpsipz_gen);
  treein_rw_gen->SetBranchAddress("J_psi_1S_TRUEPT",&myjpsipt_gen);

  TBranch *wtbranch_rec = treein_rw_rec->Branch("pjpsiwt", &pjpsiwt_rec, "pjpsiwt/F");
  cout<<"Looping over reconstructed tree"<<endl;
  Int_t bin_neighbour;
  Double_t wt_center;
  for(Int_t i=0; i<nentries_mc_rec; i++) {
    if(i%10000 == 0)
      cout<<i<<endl;

    treein_rw_rec->GetEntry(i);

    if(mylbp_rec > 0 && myjpsip_rec > 0)
      {
	mypjpsi_rec = myjpsip_rec/mylbp_rec;
	Int_t globBin  = wts->FindBin(mypjpsi_rec);

	if(interpolate)
	  {
	    
	    // pjpsiwt_rec   = wts->GetBinContent(globBin);
	    wt_center   = wts->GetBinContent(globBin);

	    Double_t x_center = wts->GetBinCenter(globBin);

	    if(mypjpsi_rec < x_center)
	      {
		if(globBin!=1)
		  {
		    bin_neighbour = globBin-1;
		  }
		else
		  {
		    pjpsiwt_rec = wt_center;
		    continue;
		  }
	      }
	    else
	      {
		if(globBin != nbins_pjpsi)
		  {
		    bin_neighbour = globBin+1;
		  }
		else
		  {
		    pjpsiwt_rec = wt_center;
		    continue;
		  }
	      }

	    Double_t x_neighbour = wts->GetBinCenter(bin_neighbour);

	    Double_t wt_neighbour = wts->GetBinContent(bin_neighbour);

	    Double_t slope = (wt_center-wt_neighbour)/(x_center-x_neighbour);
	    Double_t c = wt_center - slope*x_center;

	    pjpsiwt_rec = slope*(mypjpsi_rec-x_neighbour) + c;
	  }
	else
	  {
	    pjpsiwt_rec   = wts->GetBinContent(globBin);
	  }
	// pjpsibin = TMath::CeilNint((Double_t)(mypjpsi_rec-pjpsilow)/binwidth_pjpsi);

	// if(mypjpsi_rec > pjpsihigh)
	//      pjpsibin = nbins_pjpsi;
	// if(mypjpsi_rec < pjpsilow)
	//      pjpsibin = 1;
	//
	// pjpsiwt_rec = weights[pjpsibin-1];
      }
    else
      {
	pjpsiwt_rec = 1.0;
      }
    wtbranch_rec->Fill();
  }
  treein_rw_rec->Write("",TObject::kOverwrite);

  TBranch *wtbranch_gen = treein_rw_gen->Branch("pjpsiwt", &pjpsiwt_gen, "pjpsiwt/F");
  cout<<"Looping over generated tree"<<endl;
  for(Int_t i=0; i<nentries_mc_gen; i++) {
    if(i%100000 == 0)
      cout<<i<<endl;

    treein_rw_gen->GetEntry(i);

    myjpsip_gen = sqrt(pow(myjpsipz_gen,2)+pow(myjpsipt_gen,2));
    mylbp_gen = sqrt(pow(mylbpz_gen,2)+pow(mylbpt_gen,2));

    if(mylbp_gen > 0 && myjpsip_gen > 0)
      {
	mypjpsi_gen = myjpsip_gen/mylbp_gen;
	Int_t globBin = wts->FindBin(mypjpsi_gen);
	pjpsiwt_gen   = wts->GetBinContent(globBin);

	// pjpsibin = TMath::CeilNint((Double_t)(mypjpsi_gen-pjpsilow)/binwidth_pjpsi);
	//
	// if(mypjpsi_gen > pjpsihigh)
	//      pjpsibin = nbins_pjpsi;
	// if(mypjpsi_gen < pjpsilow)
	//      pjpsibin = 1;
	//
	// pjpsiwt_gen = weights[pjpsibin-1];
      }
    else
      {
	pjpsiwt_gen = 1.0;
      }
    wtbranch_gen->Fill();
  }
  treein_rw_gen->Write("",TObject::kOverwrite);

  filein_rw->Close();
  cout<<"DONE"<<endl;

  // gDirectory->Add(pjpsi_data);
  // gDirectory->Add(pjpsi_mc);
}
