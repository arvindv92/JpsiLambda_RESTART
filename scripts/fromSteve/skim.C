#include "CodeInput.h"
#include <iostream>
#include <fstream>
using std::ifstream;

using namespace std;
using namespace RooFit;

void skim(){

  TTree* treeXib[6];
  TTree* treeXibNew[6];

  // Get data files
  year = "Run1";
  getDataFiles();
  treeXib[4] = treeDD;
  year = "Run2";
  getDataFiles();
  treeXib[5] = treeDD;
  
  // Get cuts
  getCuts();  

  // Create trees with cuts
  TCut longCut = Xib_Cut && XiL_Cut;
  TCut downCut = Xib_Cut && XiD_Cut;

  // Run 1 Lambda_DD pi_D
  TString fileName = "/data1/sblusk/jpsilambdaDD_Run1.root";
  TFile *fout = new TFile(fileName,"RECREATE");
  treeXibNew[4] = treeXib[4]->CopyTree(LambdabDD_Cut);
  treeXibNew[4]->Write();
  fout->Write();
  fout->Close();

  // Run  2 Lambda_DD pi_D
  TString fileName = "/data1/sblusk/jpsilambdaDD_Run2.root";
  TFile *fout = new TFile(fileName,"RECREATE");
  treeXibNew[5] = treeXib[5]->CopyTree(LambdabDD_Cut);
  treeXibNew[5]->Write();
  fout->Write();
  fout->Close();  



  return;

}
