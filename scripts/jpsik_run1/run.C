#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TCut.h>
#include <iostream>

using namespace std;
void run()
{
   TChain *h1 = new TChain("B2JpsiKTree/MyTuple");

    for(Int_t i=0; i<=272;i++)
    {
      h1->Add(TString::Format("../2011_MagDown/1158_%d/jpsik.root",i),-1);
    }

  cout<<"DONE ATTACHING ROOT FILES"<<endl;

  TCut dtfchi2ndofcut = "(B_DTF_VCHI2NDOF_JpsiConstr > 0 && B_DTF_VCHI2NDOF_JpsiConstr < 20)";
  TCut triggercut = "(Jpsi_Hlt1DiMuonHighMassDecision_TOS==1||Jpsi_Hlt1TrackMuonDecision_TOS==1||Jpsi_Hlt1TrackAllL0Decision_TOS==1)&&(Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS==1)";
  TCut bminipchi2cut = "(B_MINIPCHI2 < 25)";
  TCut lifetimecut = "(B_TAU > 0)";
  TCut diracut = "(B_DIRA_OWNPV > 0.999)";
  TCut pidcut = "(K_PIDK != 0 && K_PIDp != 0)";
  TCut masscut = "(B_DTF_M_JpsiConstr > 5000 && B_DTF_M_JpsiConstr < 5600)";
  
  TCut totalcut = lifetimecut && triggercut && dtfchi2ndofcut && diracut && masscut && bminipchi2cut && pidcut;

  TFile *outfile = new TFile("jpsik.root","RECREATE");
  cout<<"Making outfile"<<endl;
  TTree *outtree = (TTree*)h1->CopyTree(totalcut);

  outfile->Write();
  outfile->Close();
}
