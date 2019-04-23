#include "CodeInput.h"
//#include "CodeInput.July17.mod.h"

using namespace std;
using namespace RooFit;

TString tagCond = "_def";

void compLbLL(){
  lambdaType = "LL";
  parentPart = "Lb";  
  makePDF = false;
  mcOpt = 5;  //def = 4
  year = "Run2";

  // Get Cuts
  getCuts();
  simCut = MCLambdabLL_Cut;

  // Get needed files
  getLbLLFiles();

  // Get the histogram weights
  getWeightHistos();  

  // Get the MC tree for plotting MC results
  if(mcOpt>=2){
    getNewTree(chainl);
  }else{
    newtree = (TTree*)chainl;
  }
  
  // print out cuts
  cout << simCut.Print() << endl;  

  // Book histograms
  bookLbBaseHistos();

  // Fill histograms
  fillLbHistos(TreeWithWeight,newtree);

  // make Plots
  makeLbBasePlots();
  makeLbPlotsLL();

  //makeWeightFile()
}

void compLbDD(){
  lambdaType = "DD";
  parentPart = "Lb";  
  makePDF = true;
  mcOpt = 5;  //def = 5

  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);

  // Get Cuts
  getCuts();
  simCut = MCLambdabDD_Cut;

  // Get needed files
  getLbDDFiles();

  // Get the histogram weights
  getWeightHistos();  

  // Get the MC tree for plotting MC results
  if(mcOpt>=2){
    getNewTree(chainl);
  }else{
    newtree = (TTree*)chainl;
  }
  
  // print out cuts
  cout << simCut.Print() << endl;  

  // Book histograms
  bookLbBaseHistos();

  // Fill histograms
  fillLbHistos(TreeWithWeight,newtree);

  // make Plots
  makeLbBasePlots();

  if(MakeWeightFile && mcOpt<=3) makeWeightFile();
  
}


void compXib(){
  lambdaType = "DD";
  parentPart = "Xib";
  year = "Run1";
  
  makePDF = true;
  mcOpt = 5;  //def = 5

  // Get Cuts
  getCuts();
  simCut = MCXib_Cut&&XiL_Cut;
  //return;
  
  
  // Get needed files
  getXibFiles();

  // Get the histogram weights
  getWeightHistos();  

  // Get the MC tree for plotting MC results
  if(mcOpt>=2){
    getNewTreeXib(chainl);
  }else{
    newtree = (TTree*)chainl;
  }
  
  // print out cuts
  cout << simCut.Print() << endl;  

  // Book histograms
  //bookLbBaseHistos();
  bookXibBaseHistos();

  // Fill histograms
  fillLbHistos(TreeWithWeight,newtree);
  //fillLbHistos();

  //return;
  
  // make Plots
  makeLbBasePlots();
  //fsigdata->Close("R");
  makeXibPlots();

  
  
  return;
  
}


void compLbXib(){
  lambdaType = "DD";
  parentPart = "Xib";
  
  makePDF = false;
  mcOpt = 4;  //def = 4

  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);

  // Get Cuts
  getCuts();

  // Get needed files
  TChain* xibChain = getXibFiles();
  TChain* lbChain = getLbDDFiles();

  cout << xibChain->GetEntries() << " " << lbChain->GetEntries() << endl;
  
  //return;
  

  // Get the histogram weights
  getWeightHistos();  

  // Get the MC tree for plotting MC results
  if(mcOpt>=2){
    simCut = MCXib_Cut;
    TTree* treeXib = getNewTreeXib(xibChain);
    MCXib_Cut = simCut;
    simCut = MCLambdabDD_Cut;
    TTree* treeLb = getNewTree(lbChain);
    MCLambdabDD_Cut = simCut;
  }else{
    newtree = (TTree*)chainl;
  }

  
  TCanvas *c1 = new TCanvas("c1","p",800,800);
  TH1F *h1 = new TH1F("h1","pion mom",100,0,20);
  TH1F *h2 = new TH1F("h2","pion mom",100,0,20);  
  treeXib->Draw("0.001*pi_P>>h1",MCXib_Cut);
  treeLb->Draw("0.001*pi_P>>h2",MCLambdabDD_Cut);
  h2->SetLineColor(2);
  h2->Scale(h1->Integral()/h2->Integral());
  h1->Draw();
  h2->Draw("same");
  



  return;
  

  // print out cuts
  cout << simCut.Print() << endl;  

  // Book histograms
  bookXibBaseHistos();

  // Fill histograms
  fillLbHistos(TreeWithWeight,newtree);

  // make Plots
  makeLbBasePlots();
  makeXibPlots();

}



