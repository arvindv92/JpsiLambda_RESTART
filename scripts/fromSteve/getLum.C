#include <sys/stat.h>
#include <unistd.h>
#include <string>
#include "CodeInput.h"

void getLumRun1_Xib(){

  TString year = "Run1";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiXi/data/";
  TString yd[4] = {"2011_MagDown", "2011_MagUp", "2012_MagDown", "2012_MagUp"}
  TString iJob[4] = {"1059", "1060", "1061", "1062"};
  //int nJob[4] = {10, 10, 10, 10};
  int nJob[4] = {267, 226, 336, 393};
  TString fullpath, file, fullfile, filename;
  bool ex;

  TFile *f;
  TTree* tr;

  double IntegratedLuminosity;
  double lum;

  for(int iy=0; iy<4; iy++){
    for(int j=0; j<=nJob[iy]; j++){
      fullpath = dir+yd[iy]+"/"+iJob[iy];
      file = Form("_%d/jpsixi.root",j);
      filename = fullpath+file;
      ex = exists(filename);
      //cout << filename << " " << ex << endl;
      if(ex){
        f = new TFile(filename);
        tr = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
        tr->SetBranchAddress("IntegratedLuminosity", &IntegratedLuminosity);
        tr->GetEntry(0);
        lum = lum + IntegratedLuminosity;
      }
    }
  }

  cout << lum << endl;
  
  return;

}

void getLumRun1_Lb(){

  TString year = "Run1";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiLambda/massdump/data/";
  TString yd[4] = {"2011_MagDown", "2011_MagUp", "2012_MagDown", "2012_MagUp"}
  //TString iJob[4] = {"1059", "1060", "1061", "1062"};
  int nJob[4][2] = {
    {253, 20},
    {213, 27},
    {342, 16},
    {376, 45}
  };
  
  int iJob[4][2] = {
    {887, 981},
    {888, 982},
    {889, 952},
    {890, 983}
  };

  TString fullpath, file, fullfile, filename;
  bool ex;

  TFile *f;
  TTree* tr;

  double IntegratedLuminosity;
  double lum;
  TString js;

  for(int iy=0; iy<4; iy++){
    for(int k=0; k<2; k++){
      js = Form("%d",iJob[iy][k]);
      fullpath = dir+yd[iy]+"/"+js;//iJob[iy][k];
      cout << fullpath << endl;
      for(int j=0; j<=nJob[iy][k]; j++){
        file = Form("_%d/jpsilambda.root",j);
        filename = fullpath+file;
        ex = exists(filename);
        //cout << filename << " " << ex << endl;
        if(ex){
          f = new TFile(filename);
          tr = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
          tr->SetBranchAddress("IntegratedLuminosity", &IntegratedLuminosity);
          tr->GetEntry(0);
          lum = lum + IntegratedLuminosity;
        }
      }
    }
    
  }

  cout << lum << endl;
  
  return;

}

void getLumRun2_Lb(){

  TString year = "Run2";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiLambda/massdump/data/";
  TString yd[2] = {"2016_MagDown", "2016_MagUp"};
  int nJob[2][2] = {
    {637, 153},
    {614, 169}
  };
  
  int iJob[2][2] = {
    {1044, 1149},
    {1045, 1150}
  };

  TString fullpath, file, fullfile, filename;
  bool ex;

  TFile *f;
  TTree* tr;

  double IntegratedLuminosity;
  double lum;
  TString js;

  for(int iy=0; iy<2; iy++){
    for(int k=0; k<2; k++){
      js = Form("%d",iJob[iy][k]);
      fullpath = dir+yd[iy]+"/"+js;//iJob[iy][k];
      cout << fullpath << endl;
      for(int j=0; j<=nJob[iy][k]; j++){
        file = Form("_%d/jpsilambda.root",j);
        filename = fullpath+file;
        ex = exists(filename);
        //cout << filename << " " << ex << endl;
        if(ex){
          f = new TFile(filename);
          tr = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
          tr->SetBranchAddress("IntegratedLuminosity", &IntegratedLuminosity);
          tr->GetEntry(0);
          lum = lum + IntegratedLuminosity;
        }
      }
    }
    
  }

  cout << lum << endl;
  
  return;

}

void getLumRun2_Xib(){

  TString year = "Run2";

  TChain* ch = new TChain("","MyTuple");
  TString dir = "/data1/avenkate/JpsiXi/data/";
  TString yd[2] = {"2016_MagDown", "2016_MagUp"};
  int nJob[2][2] = {
    {617, 0},
    {615, 0}
  };
  
  int iJob[2][2] = {
    {1053, 0},
    {1054, 0}
  };

  TString fullpath, file, fullfile, filename;
  bool ex;

  TFile *f;
  TTree* tr;

  double IntegratedLuminosity;
  double lum;
  TString js;

  for(int iy=0; iy<2; iy++){
    for(int k=0; k<2; k++){
      js = Form("%d",iJob[iy][k]);
      fullpath = dir+yd[iy]+"/"+js;//iJob[iy][k];
      cout << fullpath << endl;
      for(int j=0; j<=nJob[iy][k]; j++){
        file = Form("_%d/jpsixi.root",j);
        filename = fullpath+file;
        ex = exists(filename);
        //cout << filename << " " << ex << endl;
        if(ex){
          f = new TFile(filename);
          tr = (TTree*)f->Get("GetIntegratedLuminosity/LumiTuple");
          tr->SetBranchAddress("IntegratedLuminosity", &IntegratedLuminosity);
          tr->GetEntry(0);
          lum = lum + IntegratedLuminosity;
        }
      }
    }
    
  }

  cout << lum << endl;
  
  return;

}

inline bool exists(TString name) {
  struct stat buffer;   
  return (stat (name, &buffer) == 0); 

}
