#define checkDup_cxx
#include "checkDup.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void checkDup::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L checkDup.C
//      Root > checkDup t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
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

   Long64_t nentries = fChain->GetEntriesFast();

   TFile *fDup = new TFile("DupXib","RECREATE");
   TH1F* hdup = new TH1F("hdup","Duplicates",10,0.0,10.0);
   
   int ndup = 0;
   ULong64_t lastEv = 0;
   UInt_t    lastRun = 0;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(runNumber==lastRun && eventNumber==lastEv) ndup++;
      //if(runNumber==lastRun) ndup++;

      if(!(runNumber==lastRun && eventNumber==lastEv)){
        //if(!(runNumber==lastRun)){        
        lastEv = eventNumber;
        lastRun = runNumber;
        hdup->Fill(ndup);
        ndup = 1;
      }
      if(jentry < 100) cout << runNumber << " " << eventNumber << endl;
      

      // if (Cut(ientry) < 0) continue;
   }

   hdup->Draw();
}
