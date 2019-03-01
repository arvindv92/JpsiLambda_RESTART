#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
using namespace std;
void ppimomasym()
{
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

  //  if(flag==1)
  gSystem->Exec("cp rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_pidgen_ptetawt.root rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_pidgen_ppiwt.root");

  TFile *filein_data,*filein_mc,*filein_wts;
  TTree *treein_data,*treein_mc,*treein_wts;
  
  Int_t nentries_mc, nentries_data;
  Int_t nbins_ppimomasym = 20;
  Int_t ppimomasymbin;

  Float_t weights[20];
  Float_t binwidth_ppimomasym;
  Float_t ppimomasymwt;

  Double_t myppimomasym, myproton_p, mypion_p;
  Float_t ppimomasymlow, ppimomasymhigh;

  TH1D *ppimomasym_data, *ppimomasym_mc;

  ppimomasymlow = 0.5;
  ppimomasymhigh = 0.9;

  binwidth_ppimomasym = (Float_t)(ppimomasymhigh-ppimomasymlow)/nbins_ppimomasym;

  filein_data = TFile::Open("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LL_withsw.root","READ");
  filein_mc = TFile::Open("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_pidgen_ptetawt.root","READ");

  treein_data = (TTree*)filein_data->Get("MyTuple");
  treein_mc = (TTree*)filein_mc->Get("MyTuple");

  nentries_mc = treein_mc->GetEntries();
  nentries_data = treein_data->GetEntries();

  treein_data->Draw(Form("((p_P-pi_P)/(p_P+pi_P))>>ppimomasym_data(%d,%f,%f)",nbins_ppimomasym,ppimomasymlow,ppimomasymhigh),"SW","goff");
  treein_mc->Draw(Form("((p_P-pi_P)/(p_P+pi_P))>>ppimomasym_mc(%d,%f,%f)",nbins_ppimomasym,ppimomasymlow,ppimomasymhigh),"ptetawt*(Lb_BKGCAT==0||Lb_BKGCAT==50)","goff");

  ppimomasym_data = (TH1D*)gDirectory->Get("ppimomasym_data");
  ppimomasym_mc = (TH1D*)gDirectory->Get("ppimomasym_mc");

  cout<<"Data integral = "<<ppimomasym_data->Integral()<<endl;
  cout<<"MC integral = "<<ppimomasym_mc->Integral()<<endl;

  ppimomasym_data->Scale(1.0/ppimomasym_data->Integral());
  ppimomasym_mc->Scale(1.0/ppimomasym_mc->Integral());

  Float_t binc_data, binc_mc;
  cout<<"poop1"<<endl;
  for(Int_t i=1;i<=nbins_ppimomasym;i++) {

      binc_data = ppimomasym_data->GetBinContent(i);
      binc_mc = ppimomasym_mc->GetBinContent(i);

      if(binc_mc != 0 && binc_data > 0)
	weights[i-1] = (Float_t)binc_data/binc_mc;
   
      else if(binc_mc != 0 && binc_data <= 0) 
	weights[i-1] = 0;

      else if(binc_mc == 0) 
	weights[i-1] = 1;

      if(weights[i-1] > 10)
	cout<<"REDFLAG data = "<<binc_data<<" mc = "<<binc_mc<<" i = "<<i<<endl;
  }
  cout<<"poop2"<<endl;
  
  filein_wts = TFile::Open("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run1/jpsilambda_pidgen_ppiwt.root","UPDATE");
  
  treein_wts = (TTree*)filein_wts->Get("MyTuple");

  treein_wts->SetBranchAddress("p_P",&myproton_p);
  treein_wts->SetBranchAddress("pi_P",&mypion_p);

  TBranch *wtbranch = treein_wts->Branch("ppimomasymwt", &ppimomasymwt, "ppimomasymwt/F");
  cout<<"poop3"<<endl;
  for(Int_t i=0; i<nentries_mc; i++) {
    if(i%1000 == 0)
      cout<<i<<endl;

    treein_wts->GetEntry(i);
    
    myppimomasym = ((myproton_p-mypion_p)/(myproton_p+mypion_p));

    ppimomasymbin = TMath::CeilNint((Double_t)(myppimomasym-ppimomasymlow)/binwidth_ppimomasym);

    if(myppimomasym > ppimomasymhigh) 
      ppimomasymbin = nbins_ppimomasym;
   
    if(myppimomasym < ppimomasymlow) 
      ppimomasymbin = 1;

    ppimomasymwt = weights[ppimomasymbin-1];

    wtbranch->Fill();
  }
  treein_wts->Write("",TObject::kOverwrite);
  cout<<"DONE"<<endl;
}
