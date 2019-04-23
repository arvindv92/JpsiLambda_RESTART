#include "CodeInput.h"

void test(){

  TFile *f = new TFile("/data1/avenkate/JpsiXi/data/total/jpsixi.root");
  TTree *tree = (TTree*)f->Get("Xib2JpsiXiTree/MyTuple");

  getCuts();
  


  TCut Xib_Base = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_IPCHI2_OWNPV<12 && Xib_TAU*1000>0.2";
  if(applyEtaCutToHb) Xib_Base = Xib_Base && "(-log(tan(0.5*asin(Xib_PT/Xib_P)))>2 && -log(tan(0.5*asin(Xib_PT/Xib_P)))<6)";
  Xib_CutLLL = Xib_Base && Xib_Trig && Jpsi_Cut && LambdaLL_Cut && Xi_Cut && XiL_Cut ;
  Xib_CutDDL = Xib_Base && Xib_Trig && Jpsi_Cut && LambdaDD_Cut && Xi_Cut && XiL_Cut ;


  TH1F *h1 = new TH1F("h1","#Xi_{b}^{-}#rightarrow J/#psi#Xi^{-}_{LLL}",62,5640,5950);
  TH1F *h2 = new TH1F("h2","#Xi_{b}^{-}#rightarrow J/#psi#Xi^{-}_{DDL}",62,5640,5950);
  TCanvas *c = new TCanvas("c","",800,400);
  c->Divide(2,1);
  c->cd(1);
  tree->Draw("Xib_DTF_M_JpsiXiLConstr>>h1",Xib_CutLLL);
  c->cd(2);
  tree->Draw("Xib_DTF_M_JpsiXiLConstr>>h2",Xib_CutDDL);
  
    
  
}
