#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <iostream>

using namespace std;

void docut();

void addswbranch()
{
  TFile *swfile, *filein;
  TTree *swtree, *treein;
  TBranch *swbranch, *bwbranch;
  Int_t nentries;

  Double_t SW, BW;

  TLorentzVector B, K;

  //first cutdown the file to 5200 - 5360 MeV in B mass
  docut();
  
  swfile = TFile::Open("../jpsik_withsw.root","READ");
  swtree = (TTree*)swfile->Get("MyTuple");

  swtree->SetBranchAddress("SIGW",&SW);
  swtree->SetBranchAddress("BKGW",&BW);

  filein = TFile::Open("jpsik_forTMVAisoTraining_temp.root","update");
  treein = (TTree*)filein->Get("MyTuple");
  
  swbranch = treein->Branch("SW",&SW,"SW/D");
  bwbranch = treein->Branch("BW",&BW,"BW/D");
  nentries = treein->GetEntries();
 
  cout<<"Looping over "<<nentries<<" jpsik entries"<<endl;

  for(Int_t i=0; i<nentries; i++){
    if(i%10000 == 0) cout<<i<<endl;
    treein->GetEntry(i);
    swtree->GetEntry(i);

    swbranch->Fill();
    bwbranch->Fill();
  }

  treein->Write("",TObject::kOverwrite);
  filein->Close();

}

void docut()
{
  TFile *filein_og, *fileout_cut;
  TTree *treein_og, *treeout_cut;

  filein_og = TFile::Open("../../jpsik.root","read");
  treein_og = (TTree*)filein_og->Get("MyTuple");

  fileout_cut = new TFile("jpsik_forTMVAisoTraining_temp.root","recreate");

  cout<<"Making mass cut on jpsik file"<<endl;

  treeout_cut = (TTree*)treein_og->CopyTree("B_DTF_M_JpsiConstr > 5200 && B_DTF_M_JpsiConstr < 5360");

  cout<<"Done making mass cut on jpsik file"<<endl;

  fileout_cut->Write();
  fileout_cut->Close();
}
 
