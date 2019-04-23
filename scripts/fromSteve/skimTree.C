#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include "CodeInput.h"

void skimXib_Run1(){

  TString year = "Run1";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiXi/data/";
  TString yd[4] = {"2011_MagDown", "2011_MagUp", "2012_MagDown", "2012_MagUp"}
  TString iJob[4] = {"1059", "1060", "1061", "1062"};
  //int nJob[4] = {10, 10, 10, 10};
  int nJob[4] = {267, 226, 336, 393};
  TString fullpath, file, fullfile, filename;
  bool ex;
  

  for(int iy=0; iy<4; iy++){
    for(int j=0; j<=nJob[iy]; j++){
      fullpath = dir+yd[iy]+"/"+iJob[iy];
      file = Form("_%d/jpsixi.root",j);
      filename = fullpath+file;
      ex = exists(filename);
      cout << fullfile << " " << ex << endl;
      if(ex){
        fullfile = fullpath+file+"/Xib2JpsiXiTree/MyTuple";      
        ch->Add(fullfile);                  
      }
    }
  }

  getLooseCuts();
  
  TTree *tree = (TTree*)ch;
  TFile *fout = new TFile("JpsiXi_LooseCutData_"+year+".root","RECREATE");  
  TTree* newtree = tree->CopyTree(Xib_Cut);
  newtree->SetName("MyTuple");
  newtree->Write();
  fout->Write();
  
  return;
  


}

void skimXib_Run2(){

  TString year = "Run2";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiXi/data/";
  TString yd[2] = {"2016_MagDown", "2016_MagUp"};
  TString iJob[2] = {"1053", "1054"};
  int nJob[2] = {617, 615};
  TString fullpath, file, fullfile, filename;
  bool ex;
  

  for(int iy=0; iy<2; iy++){
    for(int j=0; j<=nJob[iy]; j++){
      fullpath = dir+yd[iy]+"/"+iJob[iy];
      file = Form("_%d/jpsixi.root",j);
      filename = fullpath+file;
      ex = exists(filename);
      cout << fullfile << " " << ex << endl;
      if(ex){
        fullfile = fullpath+file+"/Xib2JpsiXiTree/MyTuple";      
        ch->Add(fullfile);                  
      }
    }
  }

  getLooseCuts();
  
  TTree *tree = (TTree*)ch;
  TFile *fout = new TFile("JpsiXi_LooseCutData_"+year+".root","RECREATE");  
  TTree* newtree = tree->CopyTree(Xib_Cut);
  newtree->SetName("MyTuple");
  newtree->Write();
  fout->Write();
  
  return;
  


}

void skimLb_Run2(){

  TString year = "Run2";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiLambda/massdump/data/";
  TString yd[4] = {"2016_MagDown", "2016_MagDown", "2016_MagUp", "2016_MagUp"};
  //TString iJob[2] = {"1044", "1045"};
  //int nJob[2] = {627, 614};
  TString iJob[4] = {"1044", "1149", "1045", "1150"};
  int nJob[4] = {637, 153, 614, 169};
  TString fullpath, file, fullfile, filename;
  bool ex;
  

  for(int iy=0; iy<4; iy++){
    for(int j=0; j<=nJob[iy]; j++){
      fullpath = dir+yd[iy]+"/"+iJob[iy];
      file = Form("_%d/jpsilambda.root",j);
      filename = fullpath+file;
      ex = exists(filename);
      cout << fullfile << " " << ex << endl;
      if(ex){
        fullfile = fullpath+file+"/Lb2JpsiLTree/MyTuple";      
        ch->Add(fullfile);                  
      }
    }
  }

  getLooseCuts();
  
  TTree *tree = (TTree*)ch;
  TFile *fout = new TFile("/data1/sblusk/XibFrag/data/JpsiL0_LooseCutData_"+year+"_2.root","RECREATE");  
  TTree* newtree = tree->CopyTree(Lambdab_Cut);
  newtree->SetName("MyTuple");
  newtree->Write();
  fout->Write();
  
  return;
  


}

inline bool exists(TString name) {
  struct stat buffer;   
  return (stat (name, &buffer) == 0); 

}
