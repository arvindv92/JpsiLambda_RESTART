#include <TFile.h>
#include <TTree.h>
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <RooHist.h>
#include <TAxis.h>

using namespace RooFit;
using namespace std;
void fitscript()
{
  RooRealVar B_DTF_M_JpsiConstr("B_DTF_M_JpsiConstr","B_DTF_M_JpsiConstr",5200.,5360.);
  // RooRealVar mean("mean","Gaussian Mean",5.27908e+03);
  // RooRealVar sigma1("sigma1","Gaussian sigma1",7.51971e+00);
  // RooRealVar sigma2("sigma2","Gaussian sigma2",1.36722e+01);

  RooRealVar mean("mean","Gaussian Mean",5280.,5200.0,5400.0);
  RooRealVar sigma1("sigma1","Gaussian sigma1",6.3,0.0,50.0);
  RooRealVar sigma2("sigma2","Gaussian sigma2",13.0,0.0,50.0);

  RooGaussian sig1("sig1","Gaussian signal1",B_DTF_M_JpsiConstr,mean,sigma1);
  RooGaussian sig2("sig2","Gaussian signal2",B_DTF_M_JpsiConstr,mean,sigma2);
  
  RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5,0.0,1.0);
  
  RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);

  //  RooRealVar c0("c0","c0",-1.45520e-01);
  //  RooRealVar c1("c1","c2",-2.27969e-02);

  RooRealVar c0("c0","c0",-0.1,-10.,10.);
  RooRealVar c1("c1","c2",-0.1,-10.,10.);

  RooChebychev bkg("bkg","bkg",B_DTF_M_JpsiConstr,RooArgList(c0,c1));
  
  //  RooRealVar tau("tau","tau",-0.0007,-0.1,-0.0000001);
  //  RooExponential bkg("bkg","Exponential bkg",B_DTF_M_JpsiConstr,tau);

  TFile *filein = TFile::Open("./jpsik.root","READ");
  TTree *treein = (TTree*)filein->Get("MyTuple");

  treein->SetBranchStatus("*",0);
  treein->SetBranchStatus("B_DTF_M_JpsiConstr",1);

  Int_t nentries = treein->GetEntries();
  cout<<"nentries = "<<nentries<<endl;

  RooDataSet ds("ds","ds",B_DTF_M_JpsiConstr,Import(*treein));
  ds.Print();

  cout<<"stage1"<<endl;
  RooRealVar nsig("nsig","nsig",0,nentries);
  RooRealVar nbkg("nbkg","nbkg",0,nentries);

  // RooRealVar nsig("nsig","nsig",2.13176e+05);
  // RooRealVar nbkg("nbkg","nbkg",2.97442e+05);

  RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
  cout<<"stage2"<<endl;
  model.fitTo(ds,Extended());
  cout<<"stage3"<<endl;
  RooPlot *frame = B_DTF_M_JpsiConstr.frame();
  ds.plotOn(frame,Name("data"));
  model.plotOn(frame,Name("fit"));
  model.plotOn(frame,Components(sig),LineStyle(kDashed));
  model.plotOn(frame,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
  model.plotOn(frame,Components(sig2),LineStyle(kDashed),LineColor(kMagenta));
  model.plotOn(frame,Components(bkg),LineColor(kRed));
  
  Double_t chiSquare1 = frame->chiSquare("fit","data");
  cout<<"chi square1 = "<<chiSquare1<<endl;
  RooArgSet *floatpar = model.getParameters(ds);
  floatpar->Print();
  int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
  cout<<"float pars = "<<floatpars<<endl;
  Double_t chi2 = frame->chiSquare("fit","data",floatpars);
  cout<<"chi square2 = "<<chi2<<endl;

  TCanvas* c = new TCanvas("Ksfit","Ksfit", 1200, 800);
  cout<<"poop1"<<endl;
  // ///////////
  TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
  TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.5);
  pad2->SetBorderMode(0);
  pad1->SetBorderMode(0);
  c->SetBorderMode(0);
  pad2->Draw();
  pad1->Draw();
  pad1->cd();

  gPad->SetTopMargin(0.06);
  pad1->Update();
  cout<<"poop2"<<endl;

  frame->Draw();

  //  TPaveText *pt = (TPaveText*)pad1->GetListOfPrimitives()->FindObject("model_paramBox");

  // pt->AddText(Form("#chi^{2}/dof = %f",chi2));
  cout<<"poop3"<<endl;
  c->Modified();
  // Pull distribution
  RooPlot *framex2 = B_DTF_M_JpsiConstr.frame();
  RooHist* hpull = frame->pullHist("data","fit");
  framex2->addPlotable(hpull,"P");
  hpull->SetLineColor(kBlack);
  hpull->SetMarkerColor(kBlack);
  framex2->SetTitle(0);
  framex2->GetYaxis()->SetTitle("Pull");
  framex2->GetYaxis()->SetTitleSize(0.15);
  framex2->GetYaxis()->SetLabelSize(0.15);
  framex2->GetXaxis()->SetTitleSize(0.2);
  framex2->GetXaxis()->SetLabelSize(0.15);
  framex2->GetYaxis()->CenterTitle();
  framex2->GetYaxis()->SetTitleOffset(0.25);
  framex2->GetXaxis()->SetTitleOffset(1.1);
  framex2->GetYaxis()->SetNdivisions(505);
  framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
  pad2->cd();
  framex2->Draw();

  c->cd();

  //Check pull mean and RMS. Ideally the pull mean should be 0 and the RMS should be 1
  cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
  cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;
  nsig.Print();
  nbkg.Print();
}
  
