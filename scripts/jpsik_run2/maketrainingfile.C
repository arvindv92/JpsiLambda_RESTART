#include <TFile.h>
#include <TTree.h>
#include <iostream>

using namespace std;
void maketrainingfile()
{
  TFile *filein, *fileout;
  TTree *treein, *treeout;
  
  Int_t nentries;
  
  Double_t SW, BW, B_ENDVERTEX_CHI2, B_FD_OWNPV, B_FDCHI2_OWNPV, K_PT, K_MINIPCHI2, K_IPCHI2_ORIVX, K_IP_ORIVX, B_VCHI2DOF, K_TRACK_GhostProb, K_TRACK_CHI2NDOF;
  ULong64_t evtNum = 0;
    
  cout<<"I am making the files needed to train the isolation BDT for background"<<endl;
  cout<<"I am taking in the masscut file for jpsik that have the signal weights branch added to them, and then making training file with only the necessary branches with the right names"<<endl;

  filein = TFile::Open("jpsik_forTMVAisoTraining_temp.root","read");
  treein = (TTree*)filein->Get("MyTuple");

  nentries = treein->GetEntries();
  
  treein->SetBranchAddress("SW",&SW);
  treein->SetBranchAddress("BW",&BW);
  treein->SetBranchAddress("B_ENDVERTEX_CHI2",&B_ENDVERTEX_CHI2);
  treein->SetBranchAddress("B_FD_OWNPV",&B_FD_OWNPV);
  treein->SetBranchAddress("B_FDCHI2_OWNPV",&B_FDCHI2_OWNPV);
  treein->SetBranchAddress("K_MINIPCHI2",&K_MINIPCHI2);
  treein->SetBranchAddress("K_IP_ORIVX",&K_IP_ORIVX);
  treein->SetBranchAddress("K_IPCHI2_ORIVX",&K_IPCHI2_ORIVX);
  treein->SetBranchAddress("K_PT",&K_PT);
  treein->SetBranchAddress("K_TRACK_GhostProb",&K_TRACK_GhostProb);
  treein->SetBranchAddress("K_TRACK_CHI2NDOF",&K_TRACK_CHI2NDOF);
  treein->SetBranchAddress("eventNumber",&evtNum);

  fileout = new TFile("jpsik_forTMVAisoTraining.root","recreate");
  treeout = new TTree("MyTuple","Background training tree for isolation BDT");

  treeout->Branch("SW",&SW,"SW/D");
  treeout->Branch("BW",&BW,"BW/D");
  treeout->Branch("VCHI2DOF",&B_VCHI2DOF,"VCHI2DOF/D");
  treeout->Branch("FD",&B_FD_OWNPV,"FD/D");
  treeout->Branch("FDCHI2",&B_FDCHI2_OWNPV,"FDCHI2/D");
  treeout->Branch("PT",&K_PT,"PT/D");
  treeout->Branch("MINIPCHI2",&K_MINIPCHI2,"MINIPCHI2/D");
  treeout->Branch("IP",&K_IP_ORIVX,"IP/D");
  treeout->Branch("IPCHI2",&K_IPCHI2_ORIVX,"IPCHI2/D");
  treeout->Branch("GHOSTPROB",&K_TRACK_GhostProb,"GHOSTPROB/D");
  treeout->Branch("TRACKCHI2DOF",&K_TRACK_CHI2NDOF,"TRACKCHI2DOF/D");
  treeout->Branch("eventNumber",&evtNum);

  for(Int_t i=0; i<nentries; i++) {

    if(i%100000==0) cout<<i<<endl;
    treein->GetEntry(i);
    B_VCHI2DOF = B_ENDVERTEX_CHI2/3.0;
    
    treeout->Fill();
  }

  fileout->Write();
  fileout->Close();
}
  
