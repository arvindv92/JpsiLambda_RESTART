#define cuts_cxx
#include "cuts.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TCut.h>
#include <TSystem.h>
#include <iostream>

using namespace std;
void cuts::Loop()
{
//   In a ROOT session, you can do:
//      root> .L cuts.C
//      root> cuts t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
  if (fChain == 0) return;
      
  gROOT->ProcessLine(".> trigger_log.txt");
   
  TCut bcuts = "B_DTF_VCHI2NDOF_JpsiConstr > 0 && B_DTF_VCHI2NDOF_JpsiConstr < 20 && B_MINIPCHI2 < 25 && (Jpsi_Hlt1DiMuonHighMassDecision_TOS==1||Jpsi_Hlt1TrackMuonDecision_TOS==1||Jpsi_Hlt1TrackAllL0Decision_TOS==1)&&(Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS==1) &&B_DIRA_OWNPV > 0.999 && B_TAU > 0 && B_DTF_M_JpsiConstr > 5000 && B_DTF_M_JpsiConstr < 5600";

  TFile *fileout = new TFile("jpsik.root","RECREATE");
  TTree *treeout = (TTree*)fChain->CloneTree(0);

  Long64_t entries_init, entries_final;

  entries_init = fChain->GetEntries();

  cout<<"Incoming Entries = "<<entries_init<<endl;
   
  cout<<"I am making the following trigger cuts. Sit tight"<<endl;
  bcuts.Print();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<entries_init;jentry++) {
    if(jentry%100000 == 0)
      cout<<jentry<<endl;
  
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    if(B_DTF_VCHI2NDOF_JpsiConstr > 0 && B_DTF_VCHI2NDOF_JpsiConstr < 20) {
      if(B_MINIPCHI2 < 25) {
	if((Jpsi_Hlt1DiMuonHighMassDecision_TOS==1||Jpsi_Hlt1TrackMuonDecision_TOS==1||Jpsi_Hlt1TrackAllL0Decision_TOS==1)&&(Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS==1)) {
	    if(B_DIRA_OWNPV > 0.999) {
	      if(B_TAU > 0) {
		if(B_DTF_M_JpsiConstr > 5000 && B_DTF_M_JpsiConstr < 5600) {
		  treeout->Fill();
		}
	      }
	    }
	  }
      }
    }
  }
  entries_final = treeout->GetEntries();
  cout<<"Outgoing entries = "<<entries_final<<endl;

  gROOT->ProcessLine(".>");
  gSystem->Exec("cat cuts_log.txt");

  fileout->Write();
  fileout->Close();

}
