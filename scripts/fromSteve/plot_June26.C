#include "cuts_June26.h"

using namespace std;
using namespace RooFit;

bool makePDF = false;
float mass,err_mass;

// Systematic shift of relative tracking efficiency
//  0 = no shift  (nominal)
// 1 = +1 sigtma
// 1 = -1 sigma
float shift = 0;

// Do efficiency calculation and compute f_Xib / f_Lb
bool doEfficiency_And_Rate = true;

int do_splot = 0;
  // do_splot = 0 ==> No sPlot
  // do_splot = 1 ==> Do Lambda_b --> J/psi Lambda0 (DD)
  // do_splot = 2 ==> Do Lambda_b --> J/psi Lambda0 (LL)
  // do_splot = 11 ==> Do Xi_b --> J/psi Xi-


void addGraphics(TH1F *h, TString Xtitle, TString Ytitle, int iCol){
  h->SetXTitle(Xtitle);
  //h->SetFillColor(30);
  int bw = h->GetBinWidth(1);
  h->SetYTitle(Ytitle);
  h->SetStats(kFALSE);
  h->SetMinimum(0.1);
  h->SetMaximum(1.3*h->GetMaximum());
  h->SetTitleSize(0.1);
  h->SetLineColor(iCol);
  h->SetMarkerColor(iCol);
  h->SetMarkerSize(0.7);
  h->SetMarkerStyle(20);
  h->GetXaxis()->SetTitleOffset(1.0);  
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetTitleSize(0.055);  
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetLabelSize(0.045);  
  h->GetYaxis()->SetLabelSize(0.045);  
  h->SetNdivisions(404,"X");
  h->SetNdivisions(505,"Y");
  h->SetLineWidth(2);
  h->SetTitle("");
  
}

void setFrameAtt(RooPlot *frame, TString xtit, TString ytit)
{
  frame->SetMaximum(1.1*frame->GetMaximum());
  //frame->SetMaximum(35);
  frame->SetMinimum(0.1);    
  frame->SetLabelSize(0.045,"X");
  frame->SetLabelSize(0.045,"Y");
  frame->SetTitleSize(0.050,"X");
  frame->SetTitleSize(0.050,"Y");
  frame->SetTitleOffset(0.95,"X");
  frame->SetTitleOffset(1.45,"Y");
  frame->SetXTitle(xtit);
  frame->SetYTitle(ytit);
  frame->SetNdivisions(505,"X");
  frame->SetNdivisions(505,"Y");
  //frame->GetXaxis()->SetTitleFont(myconfigs->m_font);
  //frame->GetYaxis()->SetTitleFont(myconfigs->m_font);
}


TCanvas* makeRooPlot(RooAbsPdf* pdf, RooAbsData* rdh, RooRealVar *mvar, int mode, RooRealVar* nS,
                     TString tag="def",double mlow=2420, double mhi=2520, int nbin=50){
  
  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.045);
  double xt = 0.18;
  int bw = (mhi-mlow)/nbin;


  TString comps = "sig*";
  TString titleY = Form("Candidates / [%d MeV/#it{c}^{2}]",bw);
  TString legTitle, titleX;

  if(mode==0){
    titleX = "#it{J/#psi#Lambda} mass [MeV/#it{c}^{2}]";
    legTitle = "#it{#Lambda_{b}^{0}}#rightarrow #it{J/#psi#Lambda}";
  }else if(mode==1){
     titleX = "#it{J/#psi#Xi^{#font[122]{-}}} mass [MeV/#it{c}^{2}]";
     legTitle = "#it{#Xi_{b}^{#font[122]{-}}}#rightarrow #it{J/#psi#Xi^{#font[122]{-}}}";
  }  


  TCanvas *c = new TCanvas("c_"+tag, "Mass Plot "+tag,0,0,600,600);
  c->SetRightMargin(0.05);
  c->SetLeftMargin(0.15);
  c->SetBottomMargin(0.11);
  RooPlot* plot2 = mvar->frame(mlow,mhi,nbin);
  plot2->SetTitle("");
  rdh->plotOn(plot2,MarkerStyle(20),MarkerSize(0.75),Name("Data"));
  pdf->plotOn(plot2,LineColor(kBlue),LineWidth(2),Name("Full")); 
  pdf->plotOn(plot2,Components(comps),LineColor(kRed),LineWidth(2),Name("Signal"));
  setFrameAtt(plot2,titleX, titleY);
  plot2->SetMinimum(0.0001);
  plot2->SetMaximum(1.2*plot2->GetMaximum());
  plot2->Draw();
  //hh->SetLineColor(3);
  //hh->Draw("same");

  float ns = nS->getVal();
  float ens = nS->getError();
  TString y = Form("N_{sig} =  %4.0f #pm %2.0f",ns,ens);
  TString y1 = Form("M =  %7.2f #pm %3.2f MeV/#it{c}^{2}",mass,err_mass);
  myLatex->SetTextSize(0.065);
  myLatex->DrawLatex(xt,0.85,"LHCb");
  myLatex->SetTextSize(0.04);
  myLatex->DrawLatex(xt,0.77,y);
  //myLatex->DrawLatex(xt,0.73,y1);
  //myLatex->DrawLatex(xt,0.8,lt);

  //c->Print("./MassPlot_"+tag+".png");

  TLegend* legend1 = new TLegend(0.64,0.65,0.88,0.88);
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->SetTextSize(0.05);
  legend1->AddEntry("Data","Data","LEP");
  legend1->AddEntry("Full","Full PDF","L");
  legend1->AddEntry("Signal",legTitle,"L");
  legend1->Draw();

  return c;
}


void plot_June26(){

  int mcOpt = 1;

  nbinLb = (mHiLb - mLowLb) / bwLb;
  nbinXib = (mHiXib - mLowXib) / bwXib;
  
  TFile *fw = new TFile("fsig_Lb.root");
  TH2F *whist = (TH2F*)fw->Get("hratio");
  TFile *ftr = new TFile("ratio2012S20.root");
  TH2D *trackEff = (TH2D*)ftr->Get("Ratio");


  TFile *fL_ll = new TFile("/data1/avenkate/JpsiLambda/massdump/data/total/jpsilambda_LL.root");
  TFile *fL_dd = new TFile("/data1/avenkate/JpsiLambda/massdump/data/total/jpsilambda_DD.root");
  TFile *fXi = new TFile("/data1/avenkate/JpsiXi/data/total/jpsixi.root");
  
  TTree* treeLL = (TTree*)fL_ll->Get("MyTuple");
  TTree* treeDD = (TTree*)fL_dd->Get("MyTuple");
  TTree* treeXi = (TTree*)fXi->Get("Xib2JpsiXiTree/MyTuple");
  
  getCuts();
  

  // Data
  
  TH1F *h1 = new TH1F("h1","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
  TH1F *h1LL = new TH1F("h1LL","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
  TH1F *h2 = new TH1F("h2","#Xi_{b}^{0} mass",nbinXib,mLowXib,mHiXib);
  TH1F *h2L = new TH1F("h2L","#Xi_{b}^{0} mass, Long",nbinXib,mLowXib,mHiXib);
  TH1F *h2D = new TH1F("h2D","#Xi_{b}^{0} mass, Down",nbinXib,mLowXib,mHiXib);
  h1->SetMinimum(0);
  h2->SetMinimum(0);
  h2L->SetMinimum(0);
  h2D->SetMinimum(0);


  TCanvas *c = new TCanvas("c","",800,400);
  c->Divide(2,1);
  c->cd(1);
  treeDD->Draw("Lb_DTF_M_JpsiLConstr>>h1",Lambdab_Cut);
  c->cd(2);
  treeLL->Draw("Lb_DTF_M_JpsiLConstr>>h1LL",LambdabLL_Cut);
  
  //return;

  
  // Fit Lambda_b signal

  RooRealVar* mvar1 = new RooRealVar("Lb_DTF_M_JpsiLConstr","",mLowLb,mHiLb);
  RooRealVar* mg1 = new RooRealVar("mg1","",5620,5610,5630);
  RooRealVar* sg1A = new RooRealVar("sg1A","",5,2,12);
  RooRealVar* sg1scale = new RooRealVar("sg1scale","",1.5,1.0,2.5);
  RooFormulaVar* sg1B = new RooFormulaVar("sg1B","","@0*@1",RooArgList(*sg1A,*sg1scale));
  //RooRealVar* sg1B = new RooRealVar("sg1B","",8,2,15);
  RooGaussian* gauss1A = new RooGaussian("gauss1A","Signal",*mvar1,*mg1,*sg1A);
  RooGaussian* gauss1B = new RooGaussian("gauss1B","Signal",*mvar1,*mg1,*sg1B);
  RooRealVar* fr1 = new RooRealVar("fr1","fraction",0.8,0.2,1.0);
  RooAddPdf* sig1 = new RooAddPdf("sig1","Signal",RooArgList(*gauss1A,*gauss1B),RooArgList(*fr1));
  RooRealVar* b1  = new RooRealVar("b1", "b1", -0.0003,-0.05,0.05);
  RooExponential* comb1 = new RooExponential("comb1","f_{back}",*mvar1,*b1);
  RooRealVar* nsig1 = new RooRealVar("nsig1","",12000,0,30000000);
  RooRealVar* ncomb1 = new RooRealVar("ncomb1","",6000,0,30000000);
  RooAbsPdf* pdf1 = new RooAddPdf("pdf1","",RooArgList(*sig1,*comb1),RooArgList(*nsig1,*ncomb1));
  RooDataHist *rdh1 = new RooDataHist("rdh1","", *mvar1, h1);
  RooDataHist *rdh1LL = new RooDataHist("rdh1LL","", *mvar1, h1LL);
  
  TString modes[2];
  modes[0] = "#Lambda_{b}^{0}#rightarrow#J/#psi#Lambda";

  pdf1->fitTo(*rdh1,Hesse(kTRUE),Strategy(1));
  mass = mg1->getVal();
  err_mass = mg1->getError();
  TCanvas* cc1 = makeRooPlot(pdf1,rdh1,mvar1,0,nsig1,modes[0],mLowLb,mHiLb,nbinLb);
  cc1->Draw();
  cc1->Update();
  double fitmass_Lb_DD = mg1->getVal();
  double err_fitmass_Lb_DD = mg1->getError();
  double nsigLb = nsig1->getVal();
  double err_nsigLb = nsig1->getError();
  
  if(makePDF){
    cc1->Print("./Plots/LbMass_DD.C");
    cc1->Print("./Plots/LbMass_DD.eps");
    cc1->Print("./Plots/LbMass_DD.png");
    cc1->Print("./Plots/LbMass_DD.pdf");
  }

  pdf1->fitTo(*rdh1LL,Hesse(kTRUE),Strategy(1));
  mass = mg1->getVal();
  err_mass = mg1->getError();
  TCanvas* cc1L = makeRooPlot(pdf1,rdh1LL,mvar1,0,nsig1,modes[0]+"LL",mLowLb,mHiLb,nbinLb);
  cc1L->Draw();
  cc1L->Update();
  double fitmass_Lb_LL = mg1->getVal();
  double err_fitmass_Lb_LL = mg1->getError();
  if(makePDF){
    cc1L->Print("./Plots/LbMass_LL.png");
    cc1L->Print("./Plots/LbMass_LL.pdf");
  }
  

  TCanvas *c1 = new TCanvas("c1","",600,600);
  c1->Divide(2,2);
  c1->cd(1);
  treeXi->Draw("Xib_DTF_M_JpsiXiLConstr>>h2",Xib_Cut);
  c1->cd(2);
  treeXi->Draw("Xib_DTF_M_JpsiXiLConstr>>h2L",Xib_Cut&&XiL_Cut);  
  c1->cd(3);
  treeXi->Draw("Xib_DTF_M_JpsiXiLConstr>>h2D",Xib_Cut&&XiD_Cut);
  c1->Update();

  //return;
  


  // Now Xib
  RooRealVar* mvar2 = new RooRealVar("Xib_DTF_M_JpsiXiLConstr","",mLowXib,mHiXib);
  RooRealVar* mg2 = new RooRealVar("mg2","",5797,5780,5810);
  RooRealVar* sg2A = new RooRealVar("sg2A","",5,2,12);
  RooRealVar* sg2scale = new RooRealVar("sg2scale","",sg1scale->getVal());
  RooFormulaVar* sg2B = new RooFormulaVar("sg2B","","@0*@1",RooArgList(*sg2A,*sg2scale));
  //RooRealVar* sg2B = new RooRealVar("sg2B","",8,2,15);
  RooGaussian* sig2 = new RooGaussian("sig2","Signal",*mvar2,*mg2,*sg2A);
  //RooGaussian* gauss2A = new RooGaussian("gauss2A","Signal",*mvar2,*mg2,*sg2A);
  RooGaussian* gauss2B = new RooGaussian("gauss2B","Signal",*mvar2,*mg2,*sg2B);
  RooRealVar* fr2 = new RooRealVar("fr2","fraction",fr1->getVal());//8,0.2,1.0);
  //RooAddPdf* sig2 = new RooAddPdf("sig2","Signal",RooArgList(*gauss2A,*gauss2B),RooArgList(*fr2));
  RooRealVar* b2  = new RooRealVar("b2", "b2", -0.0003,-2,2);
  RooExponential* comb2 = new RooExponential("comb2","f_{back}",*mvar2,*b2);
  //RooChebychev* comb2 = new RooChebychev("comb2","f_{back}",*mvar2,*b2);
  RooRealVar* nsig2 = new RooRealVar("nsig2","",12000,0,30000000);
  RooRealVar* ncomb2 = new RooRealVar("ncomb2","",6000,0,30000000);
  RooAbsPdf* pdf2 = new RooAddPdf("pdf2","",RooArgList(*sig2,*comb2),RooArgList(*nsig2,*ncomb2));

  RooDataHist *rdh2 = new RooDataHist("rdh2","", *mvar2, h2);
  RooDataHist *rdh2L = new RooDataHist("rdh2L","", *mvar2, h2L);
  RooDataHist *rdh2D = new RooDataHist("rdh2D","", *mvar2, h2D);
  
  modes[1] = "#Xi_{b}^{0}#rightarrow#J/#psi#Xi^{#font[122]{-}}";

  double yieldXib, yieldXibL, yieldXibD;
  double err_yieldXib, err_yieldXibL, err_yieldXibD;

  pdf2->fitTo(*rdh2,Hesse(kTRUE),Strategy(1));
  mass = mg2->getVal();
  err_mass = mg2->getError();
  TCanvas* cc2 = makeRooPlot(pdf2,rdh2,mvar2,1,nsig2,modes[1],mLowXib,mHiXib,nbinXib);
  cc2->Draw();
  cc2->Update();
  yieldXib = nsig2->getVal();
  err_yieldXib = nsig2->getError();
  double fitmass_Xib = mg2->getVal();
  double err_fitmass_Xib = mg2->getError();


  pdf2->fitTo(*rdh2D,Hesse(kTRUE),Strategy(1));
  mass = mg2->getVal();
  err_mass = mg2->getError();
  TCanvas* cc2D = makeRooPlot(pdf2,rdh2D,mvar2,1,nsig2,modes[1]+"D",mLowXib,mHiXib,nbinXib);
  cc2D->Draw();
  cc2D->Update();
  yieldXibD = nsig2->getVal();
  err_yieldXibD = nsig2->getError();
  double fitmassD_Xib = mg2->getVal();
  double err_fitmassD_Xib = mg2->getError();

  pdf2->fitTo(*rdh2L,Hesse(kTRUE),Strategy(1));
  mass = mg2->getVal();
  err_mass = mg2->getError();
  TCanvas* cc2L = makeRooPlot(pdf2,rdh2L,mvar2,1,nsig2,modes[1]+"L",mLowXib,mHiXib,nbinXib);
  cc2L->Draw();
  cc2L->Update();
  yieldXibL = nsig2->getVal();
  err_yieldXibL = nsig2->getError();
  double fitmassL_Xib = mg2->getVal();
  double err_fitmassL_Xib = mg2->getError();  
  
  if(makePDF){
    cc2->Print("./Plots/XibMass_LambdaDD_XiLorD.png");
    cc2L->Print("./Plots/XibMass_LambdaDD_XiL.png");
    cc2D->Print("./Plots/XibMass_LambdaDD_XiD.png");

    cc2->Print("./Plots/XibMass_LambdaDD_XiLorD.pdf");
    cc2L->Print("./Plots/XibMass_LambdaDD_XiL.pdf");
    cc2D->Print("./Plots/XibMass_LambdaDD_XiD.pdf");

    cc2->Print("./Plots/XibMass_LambdaDD_XiLorD.eps");
    cc2L->Print("./Plots/XibMass_LambdaDD_XiL.eps");
    cc2D->Print("./Plots/XibMass_LambdaDD_XiD.eps");

    cc2->Print("./Plots/XibMass_LambdaDD_XiLorD.C");
    cc2L->Print("./Plots/XibMass_LambdaDD_XiL.C");
    cc2D->Print("./Plots/XibMass_LambdaDD_XiD.C");

  }


  //return;

  if(doEfficiency_And_Rate)
  {
    


  // Lambda_b MC
  TChain *chainlt = new TChain("","MCDecayTree");
  TChain *chainl = new TChain("","MyTuple");
  // MC Truth, all events
  chainlt->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_0/jpsilambda.root/MCTuple/MCDecayTree");
  chainlt->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_0/jpsilambda.root/MCTuple/MCDecayTree");
  chainlt->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_0/jpsilambda.root/MCTuple/MCDecayTree");
  chainlt->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_0/jpsilambda.root/MCTuple/MCDecayTree");
  // reco'd
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  TCanvas *c2 = new TCanvas("c2","",600,600);

  TH1F *h1mcr = new TH1F("h1mcr","#Lambda_{b}^{0} pT",100,0,50);
  TH1F *h1mcg = new TH1F("h1mcg","#Lambda_{b}^{0} pT",100,0,50);
  double LbNum, LbDen;
  double XibNumL, XibNumD, XibNum, XibDen;
  ULong64_t eventNumber;

  if(mcOpt==2){
    int ib;
    double etaLb, w, wl;
    double Lb_P, Lb_PT, Lb_TAU, Lb_TRUETAU;
    TFile *fout = new TFile("/data1/sblusk/junk.root","RECREATE");
    TTree* newtree = chainl->CopyTree(MCLambdab_Cut);
    newtree->SetBranchAddress("Lb_P", &Lb_P);
    newtree->SetBranchAddress("Lb_PT", &Lb_PT);
    newtree->SetBranchAddress("Lb_TAU", &Lb_TAU);
    newtree->SetBranchAddress("Lb_TRUETAU", &Lb_TRUETAU);
    newtree->SetBranchAddress("eventNumber", &eventNumber);
    Long64_t nentries = newtree->GetEntriesFast();
    TProfile *ha = new TProfile("ha","weight vs p",50,0,300,0,5);
    cout << "Number of reco'd events in MC = " << nentries << endl;
    Long64_t ientry ;
    //nentries = 100;
    // Loop over reconstructed MC tuple
    for (Long64_t i=0;i<nentries; i++) {
      ientry = newtree->LoadTree(i);
      if (ientry < 0) break;
      newtree->GetEntry(i);
      if(i%10000 == 0) cout << "At entry = " << i << endl;
      etaLb = -log(tan(0.5*asin(Lb_PT/Lb_P))); 
      int ib = whist->FindBin(0.001*Lb_PT,etaLb);
      w = whist->GetBinContent(ib);
      if(w<0 || w>5) w = 1.0;
      if(Lb_TRUETAU>0) {
        wl = exp(-1000*Lb_TRUETAU/1.470)/exp(-1000*Lb_TRUETAU/1.425);
      }else{
        wl = exp(-1000*Lb_TAU/1.470)/exp(-1000*Lb_TAU/1.425);
      }
      h1mcr->Fill(0.001*Lb_PT,w*wl);
      if(eventNumber==247539 ||eventNumber==247623 ||eventNumber==206593 ){
        cout << eventNumber << " " << Lb_TRUETAU << " " << etaLb << " " << Lb_PT << " " << w << " " << wl << endl;  
      }      
    }    
    //return;
    
    // Loop over generated MC tuple    
    nentries = chainlt->GetEntriesFast();
    cout << "Number of gen'd events in MC = " << nentries << endl;
    double Lambda_b0_TRUEP_Z, Lambda_b0_TRUEPT, Lambda_b0_TRUETAU;
    chainlt->SetBranchAddress("Lambda_b0_TRUEP_Z", &Lambda_b0_TRUEP_Z);
    chainlt->SetBranchAddress("Lambda_b0_TRUEPT", &Lambda_b0_TRUEPT);
    chainlt->SetBranchAddress("Lambda_b0_TRUETAU", &Lambda_b0_TRUETAU);
    chainlt->SetBranchAddress("eventNumber", &eventNumber);
    
    for (Long64_t i=0;i<nentries; i++) {
      ientry = chainlt->LoadTree(i);
      if (ientry < 0) break;
      chainlt->GetEntry(i);
      if(i%10000 == 0) cout << "At entry = " << i << endl;
      etaLb = -log(tan(0.5*atan(Lambda_b0_TRUEPT/Lambda_b0_TRUEP_Z))); 
      int ib = whist->FindBin(0.001*Lambda_b0_TRUEPT,etaLb);
      w = whist->GetBinContent(ib);
      if(w<0 || w>5) w = 1.0;
      wl = exp(-1000*Lambda_b0_TRUETAU/1.470)/exp(-1000*Lambda_b0_TRUETAU/1.425);
      if(eventNumber==247539 ||eventNumber==247623 ||eventNumber==206593 ){
        cout << eventNumber << " "  << Lb_TRUETAU << " " << etaLb << " " << Lambda_b0_TRUEPT << " " << w << " " << wl << endl;  
      }      

      h1mcg->Fill(0.001*Lambda_b0_TRUEPT,w*wl);
    }
    LbNum = h1mcr->Integral();
    LbDen = h1mcg->Integral();
  }else{
    chainl->Draw("0.001*Lb_PT>>h1mcr",MCLambdab_Cut);
    LbNum = h1mcr->Integral();
    LbDen = chainlt->GetEntries();    
  }
  
  //cout << chainlt->GetEntries() << endl;
  double effLambdab = LbNum / LbDen;
  double err_effLambdab = sqrt(effLambdab*(1-effLambdab) / LbDen);
  
  TCanvas * c = new TCanvas("c","",800,400);
  c->Divide(2,1);
  c->cd(1);
  h1mcr->Draw();
  c->cd(2);
  h1mcg->Draw();
  

  cout << "=========================================================="<<endl;
  cout << " Reco Lb = " << h1mcr->GetEntries() << endl;
  cout << " Reco Lb (weighted) = " << LbNum << endl;
  cout << "Lambda_b efficiency = " << effLambdab << " +- " << err_effLambdab << endl;
  cout << "=========================================================="<<endl;
  
  // return;
  

  // Xi_b MC
  TChain *chainxt = new TChain("","MCDecayTree");
  TChain *chainx = new TChain("","MyTuple");
  // MC Truth, all events
  chainxt->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagDown/1001_0/jpsixi.root/MCTuple/MCDecayTree");
  chainxt->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagUp/1002_0/jpsixi.root/MCTuple/MCDecayTree");
  chainxt->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagDown/1003_0/jpsixi.root/MCTuple/MCDecayTree");
  chainxt->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagUp/1004_0/jpsixi.root/MCTuple/MCDecayTree");
  // reco'd
  chainx->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagDown/1001_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainx->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagUp/1002_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainx->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagDown/1003_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainx->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagUp/1004_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");

  TH1F *bk1 = new TH1F("bk1","Xib BKGCAT",30,0,150);
  TH1F *m1 = new TH1F("m1","Mass",50,5700,5900);
  TH1F *m2 = new TH1F("m2","Mass",50,5700,5900);
  TH1F *m3 = new TH1F("m3","Mass",50,5700,5900);
  TH1F *m4 = new TH1F("m4","Mass",50,5700,5900);

  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);

  TCanvas *c4 = new TCanvas("c4","",800,800);
  c4->Divide(2,2);
  c4->cd(1);
  chainx->Draw("Xib_BKGCAT>>bk1",Xib_Cut);
  addGraphics(bk1,"#Xi_{b} BKGCAT","",1);
  c4->cd(2);
  chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m1",Xib_Cut&&"Xib_BKGCAT<20");
  addGraphics(m1,"#Xi_{b} mass","",1);
  m1->Draw();
  myLatex->DrawLatex(0.18,0.85,"BKGCAT == 0");
  c4->cd(3);
  chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m2",Xib_Cut&&"Xib_BKGCAT==60");
  //chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m3",Xib_Cut&&"(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5795)<25)");
  addGraphics(m2,"#Xi_{b} mass","",1);
  //addGraphics(m3,"#Xi_{b} mass","",1);
  //m3->SetFillColor(4);
  m2->Draw();
  myLatex->DrawLatex(0.18,0.85,"BKGCAT == 60 (Ghost)");
  //m3->Draw("same");
  c4->cd(4);
  chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m4",Xib_Cut&&"!((Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5795)<25)||(Xib_BKGCAT<30))");
  addGraphics(m4,"#Xi_{b} mass","",1);
  m4->Draw();
  myLatex->DrawLatex(0.18,0.85,"BKGCAT == 50 (Low mass BG)");
  c4->Update();

  //return;
  


  TH1F *h2mcg = new TH1F("h2mcg","#Xi_{b}^{-} pT",100,0,50);
  TH1F *h2mcr = new TH1F("h2mcr","#Xi_{b}^{-} pT",100,0,50);
  TH1F *h2mcrL = new TH1F("h2mcrL","#Xi_{b}^{-} pT",100,0,50);
  TH1F *h2mcrD = new TH1F("h2mcrD","#Xi_{b}^{-} pT",100,0,50);

  



  TCanvas *c3 = new TCanvas("c3","",600,600);
  c3->Divide(2,2);
  c3->cd(1);

  if(mcOpt==2){
    int ib, BachPi_TRACK_Type;
    double etaXib, w, wt, ewt, etaBach, wl;
    double Xib_PT, Xib_P, Xib_TAU, Xib_TRUETAU, BachPi_P, BachPi_PT;
    //TCut ccc = MCXib_Cut&&XiL_Cut;
    TTree* newtree = chainx->CopyTree(MCXib_Cut);
    newtree->SetBranchAddress("BachPi_TRACK_Type", &BachPi_TRACK_Type);
    newtree->SetBranchAddress("BachPi_P", &BachPi_P);
    newtree->SetBranchAddress("BachPi_PT", &BachPi_PT);
    newtree->SetBranchAddress("Xib_P", &Xib_P);
    newtree->SetBranchAddress("Xib_PT", &Xib_PT);
    newtree->SetBranchAddress("Xib_TAU", &Xib_TAU);
    newtree->SetBranchAddress("Xib_TRUETAU", &Xib_TRUETAU);
    Long64_t nentries = newtree->GetEntriesFast();
    Long64_t ientry ;
    //nentries = 100;
    // Loop over reconstructed MC tuple
    for (Long64_t i=0;i<nentries; i++) {
      ientry = newtree->LoadTree(i);
      if (ientry < 0) break;
      newtree->GetEntry(i);
      if(i%10000 == 0) cout << "At entry = " << i << endl;
      etaXib = -log(tan(0.5*asin(Xib_PT/Xib_P))); 
      int ib = whist->FindBin(0.001*Xib_PT,etaXib);
      w = whist->GetBinContent(ib);
      if(w<0 || w>5) w = 1.0;
      if(useEtaWeight) w = w*(exp(-0.5*(etaXib-XibWtPar[0])**2/XibWtPar[1]**2));
      if(useBachPiWeight) w = w*(BachPiPar[0]+BachPiPar[1]*exp(-0.001*BachPi_P/BachPiPar[2]));
      // Tracking correction
      wt = 1.0;
      if(BachPi_TRACK_Type==3){
        wt = 1.00 + shift * 0.1;
        etaBach = -log(tan(0.5*asin(BachPi_PT/BachPi_P))); 
        //cout << etaBach << " " << BachPi_P << endl;
        if(etaBach>1.9 && etaBach<4.9 && BachPi_P>5000 && BachPi_P<200000){
          ib = trackEff->FindBin(0.001*BachPi_P,etaBach);
          wt = trackEff->GetBinContent(ib);
          ewt = trackEff->GetBinError(ib);
          wt = wt + shift*ewt;
        }
      }else{
        // Downstream
        wt = 1.0;//0.95 + shift * 0.05;
      }
      //if(i<100){ cout << etaBach << " " << BachPi_P << " " << wt << " " << ewt << endl; }

      if(wt<0 || wt>5) wt = 1.0;
      // systematic shift +- 0.04 ps
      if(Xib_TRUETAU>0){
        wl = exp(-1000*Xib_TRUETAU/(1.566+XibLifetimeWeight*0.04))/exp(-1000*Xib_TRUETAU/1.566);
      }else{
        wl = exp(-1000*Xib_TAU/(1.566+XibLifetimeWeight*0.04))/exp(-1000*Xib_TAU/1.566);        
      }
      
      if(BachPi_TRACK_Type==3) h2mcrL->Fill(0.001*Xib_PT,w*wt*wl);
      if(BachPi_TRACK_Type==5) h2mcrD->Fill(0.001*Xib_PT,w*wl);
      if(BachPi_TRACK_Type==3) h2mcr->Fill(0.001*Xib_PT,w*wt*wl);
      if(BachPi_TRACK_Type==5) h2mcr->Fill(0.001*Xib_PT,w*wl);
    }    



    // Loop over generated MC tuple    
    nentries = chainxt->GetEntriesFast();
    cout << "Number of gen'd events in MC = " << nentries << endl;
    double Xi_bminus_TRUEP_Z, Xi_bminus_TRUEPT, Xi_bminus_TRUETAU, piminus0_TRUEPT, piminus0_TRUEP_Z;
    chainxt->SetBranchAddress("Xi_bminus_TRUEP_Z", &Xi_bminus_TRUEP_Z);
    chainxt->SetBranchAddress("Xi_bminus_TRUEPT", &Xi_bminus_TRUEPT);
    chainxt->SetBranchAddress("Xi_bminus_TRUETAU", &Xi_bminus_TRUETAU);
    chainxt->SetBranchAddress("piminus0_TRUEPT", &piminus0_TRUEPT);
    chainxt->SetBranchAddress("piminus0_TRUEP_Z", &piminus0_TRUEP_Z);
    
    for (Long64_t i=0;i<nentries; i++) {
      ientry = chainxt->LoadTree(i);
      if (ientry < 0) break;
      chainxt->GetEntry(i);
      if(i%10000 == 0) cout << "At entry = " << i << endl;
      etaXib = -log(tan(0.5*atan(Xi_bminus_TRUEPT/Xi_bminus_TRUEP_Z))); 
      int ib = whist->FindBin(0.001*Xi_bminus_TRUEPT,etaXib);
      w = whist->GetBinContent(ib);
      if(w<0 || w>5) w = 1.0;
      if(useEtaWeight) w = w*(exp(-0.5*(etaXib-XibWtPar[0])**2/XibWtPar[1]**2));
      if(useBachPiWeight) w = w*(BachPiPar[0]+BachPiPar[1]*exp(-0.001*sqrt(piminus0_TRUEPT**2+piminus0_TRUEP_Z**2)/BachPiPar[2]));
      wl = exp(-1000*Xi_bminus_TRUETAU/(1.566+XibLifetimeWeight*0.04))/exp(-1000*Xi_bminus_TRUETAU/1.566);
      h2mcg->Fill(0.001*Xi_bminus_TRUEPT,w*wl);
    }    
  }else{
    chainx->Draw("0.001*Xib_PT>>h2mcr",MCXib_Cut);  
    c3->cd(2);
    chainx->Draw("0.001*Xib_PT>>h2mcrL",MCXib_Cut&&XiL_Cut);  
    c3->cd(3);
    chainx->Draw("0.001*Xib_PT>>h2mcrD",MCXib_Cut&&XiD_Cut);    
    chainxt->Draw("0.001*Xi_bminus_TRUEPT>>h2mcg");
  }
  XibNum = h2mcr->Integral();
  XibNumL = h2mcrL->Integral();
  XibNumD = h2mcrD->Integral();
  XibDen = h2mcg->Integral();

  c3->cd(1);
  h2mcg->Draw();
  c3->cd(2);
  h2mcr->Draw();
  c3->cd(3);
  h2mcrL->Draw();
  c3->cd(4);
  h2mcrD->Draw();  
  c3->Update();
  
  cout << chainxt->GetEntries() << endl;
  double effXib = XibNum/XibDen;//h2mc->GetEntries() / chainxt->GetEntries() ;
  double effXibL = XibNumL/XibDen;//h2mcL->GetEntries() / chainxt->GetEntries() ;
  double effXibD = XibNumD/XibDen;//h2mcD->GetEntries() / chainxt->GetEntries() ;
  double err_effXib = sqrt(effXib*(1-effXib)/ XibDen) ;
  double err_effXibL = sqrt(effXibL*(1-effXibL)/ XibDen) ;
  double err_effXibD = sqrt(effXibD*(1-effXibD)/ XibDen) ;


  double corr = (2.0/3.0)*(tauLb/tauXib)*(genEffLb/genEffXib);
  double err_corr = corr*sqrt((err_tauLb/tauLb)**2 +(err_tauXib/tauXib)**2 + (err_genEffLb/genEffLb)**2 + (err_genEffXib/genEffXib)**2);

  double lumRatio = lumJpsiXi / lumJpsiL0;
  
  double yieldRatio = yieldXib / nsigLb;
  double effRatio = effXib / effLambdab;
  double result = yieldRatio / effRatio;
  double err_yieldRatio = yieldRatio*sqrt((err_yieldXib/yieldXib)**2 + (err_nsigLb/nsigLb)**2);
  double err_result = err_yieldRatio / effRatio;

  result = result / lumRatio;
  err_result = err_result / lumRatio;

  double frag =  result * corr;
  double err_frag = frag*sqrt((err_result/result)**2 + (err_corr/corr)**2);
  double err_effRatio = effRatio*sqrt( (err_effXib/effXib)**2 + (err_effLambdab/effLambdab)**2);

  cout << "Long + Downstream Pion from Xib" << endl;
  cout << "-------------------------------" << endl;
  cout << "Yield ratio =     " << yieldRatio << " +- " << err_yieldRatio << endl;
  cout << "Xib efficiency = " << effXib << "+-" << err_effXib << endl;
  cout << "Eff   ratio =     " << effRatio << " +- " << err_effRatio << endl;
  cout << "Corr yield ratio: " << result << " +- " << err_result << endl;
  cout << "f(Xib)/f(Lb)    : " << frag << " +- " << err_frag << "+-" << 0.3*frag << endl;

  double yieldRatioL = yieldXibL / nsigLb;
  double effRatioL = effXibL / effLambdab;
  double resultL = yieldRatioL / effRatioL;
  double err_yieldRatioL = yieldRatioL*sqrt((err_yieldXibL/yieldXibL)**2 + (err_nsigLb/nsigLb)**2);
  double err_resultL = err_yieldRatioL / effRatioL;
  double err_effRatioL = effRatioL*sqrt( (err_effXibL/effXibL)**2 + (err_effLambdab/effLambdab)**2);

  resultL = resultL / lumRatio;
  err_resultL = err_resultL / lumRatio;
  
  cout << endl;
  
  double fragL =  resultL * corr;
  double err_fragL = fragL*sqrt((err_resultL/resultL)**2 + (err_corr/corr)**2);

  cout << "Long Pion from Xib" << endl;
  cout << "------------------" << endl;
  cout << "Yield ratio = " << yieldRatioL << " +- " << err_yieldRatioL << endl;
  cout << "Xib efficiency = " << effXibL << "+-" << err_effXibL << endl;
  cout << "Eff   ratio = " << effRatioL << " +- " << err_effRatioL << endl;
  cout << "Corr yield ratio: " << resultL << " +- " << err_resultL << endl;
  cout << "f(Xib)/f(Lb)    : " << fragL << " +- " << err_fragL << "+-" << 0.3*fragL << endl;
  cout << endl;

  

  double yieldRatioD = yieldXibD / nsigLb;
  double effRatioD = effXibD / effLambdab;
  double resultD = yieldRatioD / effRatioD;
  double err_yieldRatioD = yieldRatioD*sqrt((err_yieldXibD/yieldXibD)**2 + (err_nsigLb/nsigLb)**2);
  double err_resultD = err_yieldRatioD / effRatioD;

  resultD = resultD / lumRatio;
  err_resultD = err_resultD / lumRatio;

  double fragD =  resultD * corr;
  double err_fragD = fragD*sqrt((err_resultD/resultD)**2 + (err_corr/corr)**2);
  double err_effRatioD = effRatioD*sqrt( (err_effXibD/effXibD)**2 + (err_effLambdab/effLambdab)**2);

  cout << "Downstream Pion from Xib" << endl;
  cout << "------------------------" << endl;
  cout << "Yield ratio = " << yieldRatioD << " +- " << err_yieldRatioD << endl;
  cout << "Xib efficiency = " << effXibD << "+-" << err_effXibD << endl;
  cout << "Eff   ratio = " << effRatioD << " +- " << err_effRatioD << endl;
  cout << "Corr yield ratio: " << resultD << " +- " << err_resultD << endl;
  cout << "f(Xib)/f(Lb)    : " << fragD << " +- " << err_fragD << "+-" << 0.3*fragD << endl;  

  cout << endl;
  cout << "Fitted masses" << endl;
  cout << "   Lambda_b mass(DD)   : " << fitmass_Lb_DD << " +- " << err_fitmass_Lb_DD << endl;
  cout << "   Lambda_b mass(LL)   : " << fitmass_Lb_LL << " +- " << err_fitmass_Lb_LL << endl;
  cout << "   Xi_b mass(L+D)  : " << fitmass_Xib << " +- " << err_fitmass_Xib << endl;
  cout << "   Xi_b mass(L)    : " << fitmassL_Xib << " +- " << err_fitmassL_Xib << endl;
  cout << "   Xi_b mass(D)    : " << fitmassD_Xib << " +- " << err_fitmassD_Xib << endl;
  
  }
  

  
  
  if(do_splot == 0) return;

  TString sLb_Trig = Lb_Trig;
  TString tag = "DD";
  if(do_splot == 2) tag = "LL";

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if(sLb_Trig.Contains("L0Global")>0 && sLb_Trig.Contains("L0DiMuon")>0) {
    tag = tag + "_L0JPsiTOSorGlobalTIS";
    cout << "Trigger is (JpsiL0DiMuon TOS || Lb L0 Global TIS)" << endl;
  }else if(sLb_Trig.Contains("L0DiMuon")>0 && sLb_Trig.Contains("L0Global")==0  ){
    tag = tag + "_L0JPsiTOS";
    cout << "Trigger is (JpsiL0DiMuon TOS)" << endl;
  }else if(sLb_Trig.Contains("L0DiMuon")==0 && sLb_Trig.Contains("L0Global")>0 ){
    tag = tag + "_L0GlobalTIS";
    cout << "Trigger is (Lb L0 Global TIS)" << endl;
  }else{
    cout << "Shoot, something wrong here!!" << endl;
  }    
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  
  if(do_splot == 1 || do_splot==2){
    

    TString fileName = "Lb2JpsiL0_"+tag+".root";
    TFile *fout = new TFile(fileName,"RECREATE");
    
    if(do_splot==1) TTree* treeNew1 = treeDD->CopyTree(Lambdab_Cut);
    if(do_splot==2) TTree* treeNew1 = treeLL->CopyTree(LambdabLL_Cut);
    
    //--------------------------------------------
    // Variables that are needed in the sWeighting
    //--------------------------------------------
    RooRealVar* L_ENDVERTEX_Z = new RooRealVar("L_ENDVERTEX_Z","Lambda Z",-1.0e10,1.0e10);
    RooRealVar* L_TAU = new RooRealVar("L_TAU","Lambda decay time",-1.0e10,1.0e10);
    RooRealVar* L_P = new RooRealVar("L_P","Lambda mom",-1.0e10,1.0e10);
    RooRealVar* L_PT = new RooRealVar("L_PT","Lambda mom",-1.0e10,1.0e10);
    RooRealVar* Lb_p = new RooRealVar("Lb_P","Lambda mom",-1.0e10,1.0e10);
    RooRealVar* Lb_pT = new RooRealVar("Lb_PT","Lambda mom",-1.0e10,1.0e10);
    RooRealVar* Lb_OWNPV_Z = new RooRealVar("Lb_OWNPV_Z","Lambdab origin Z",-1.0e10,1.0e10);
    RooRealVar* Jpsi_pT = new RooRealVar("Jpsi_PT","Jpsi pT",-1.0e10,1.0e10);
    RooRealVar* Jpsi_p = new RooRealVar("Jpsi_P","Jpsi mom",-1.0e10,1.0e10);
    RooRealVar* L_FD_OWNPV = new RooRealVar("L_FD_OWNPV","Lambda FD",-1.0e10,1.0e10);
    RooRealVar* p_IP_OWNPV = new RooRealVar("p_IP_OWNPV","proton IP",-1.0e10,1.0e10);
    RooRealVar* pi_IP_OWNPV = new RooRealVar("pi_IP_OWNPV","pion IP",-1.0e10,1.0e10);
    RooRealVar* p_IPCHI2_OWNPV = new RooRealVar("p_IPCHI2_OWNPV","proton IP",-1.0e10,1.0e10);
    RooRealVar* pi_IPCHI2_OWNPV = new RooRealVar("pi_IPCHI2_OWNPV","pion IP",-1.0e10,1.0e10);
    RooRealVar* pi_ORIVX_X = new RooRealVar("pi_ORIVX_X","pion X prod",-1.0e10,1.0e10);
    RooRealVar* pi_ORIVX_Y = new RooRealVar("pi_ORIVX_Y","pion Y prod",-1.0e10,1.0e10);
    RooRealVar* pi_ORIVX_Z = new RooRealVar("pi_ORIVX_Z","pion Z prod",-1.0e10,1.0e10);
    RooRealVar* pi_PX = new RooRealVar("pi_PX","pion PX",-1.0e10,1.0e10);
    RooRealVar* pi_PY = new RooRealVar("pi_PY","pion PY",-1.0e10,1.0e10);
    RooRealVar* pi_PZ = new RooRealVar("pi_PZ","pion PZ",-1.0e10,1.0e10);
    RooRealVar* pi_ID = new RooRealVar("pi_ID","pion PZ",-1.0e10,1.0e10);
    RooRealVar* p_PX = new RooRealVar("p_PX","proton PX",-1.0e10,1.0e10);
    RooRealVar* p_PY = new RooRealVar("p_PY","proton PY",-1.0e10,1.0e10);
    RooRealVar* p_PZ = new RooRealVar("p_PZ","proton PZ",-1.0e10,1.0e10);
    RooRealVar* runNumber = new RooRealVar("runNumber","run #",-1.0e10,1.0e10);
    
    RooArgSet args;
    args.add(*mvar1);
    args.add(*L_TAU);
    args.add(*L_P);
    args.add(*L_PT);
    args.add(*L_ENDVERTEX_Z);
    args.add(*Lb_p);
    args.add(*Lb_pT);
    args.add(*Lb_OWNPV_Z);
    args.add(*Jpsi_pT);
    args.add(*Jpsi_p);
    args.add(*L_FD_OWNPV);
    args.add(*p_IP_OWNPV);
    args.add(*pi_IP_OWNPV);
    args.add(*p_IPCHI2_OWNPV);
    args.add(*pi_IPCHI2_OWNPV);
    args.add(*pi_ORIVX_X);
    args.add(*pi_ORIVX_Y);
    args.add(*pi_ORIVX_Z);
    args.add(*pi_PX);
    args.add(*pi_PY);
    args.add(*pi_PZ);
    args.add(*pi_ID);
    args.add(*p_PX);
    args.add(*p_PY);
    args.add(*p_PZ);
    args.add(*runNumber);

    RooDataSet* data_red = new RooDataSet("data_red","dataset with reduced vars",treeNew1,args);
    
    TString wsig = "nsig1_sw";
    
    RooArgSet *myVars1;
    if(do_splot==1) {
      pdf1->fitTo(*rdh1,Hesse(kTRUE),Strategy(2));
      myVars1 = pdf1->getParameters(*rdh1) ;
      myVars1->Print("v") ;
      myVars1->setAttribAll("Constant",kTRUE) ;
      nsig1->setConstant(kFALSE);
      ncomb1->setConstant(kFALSE);
      pdf1->fitTo(*rdh1,Extended(),Hesse(kTRUE));
    }else{
      pdf1->fitTo(*rdh1LL,Hesse(kTRUE),Strategy(2));
      myVars1 = pdf1->getParameters(*rdh1LL) ;
      myVars1->Print("v") ;
      myVars1->setAttribAll("Constant",kTRUE) ;
      nsig1->setConstant(kFALSE);
      ncomb1->setConstant(kFALSE);
      pdf1->fitTo(*rdh1LL,Extended(),Hesse(kTRUE));
    }
    
      //--------------
    // Get sWeights
    //--------------
    RooStats::SPlot* sDataA1 = new RooStats::SPlot("sDataA1","An SPlot", *data_red, pdf1 , RooArgList(*nsig1,*ncomb1) );
    
    data_red->Print("v") ; 
    gStyle->SetOptStat(011000);
    //--------------------------
    // create weighted data set 
    //--------------------------
    RooDataSet * dataw_signal = new RooDataSet(data_red->GetName(),data_red->GetTitle(),data_red,*data_red->get(),0,wsig) ;
    dataw_signal->Print("v") ; 
    //-------------------------------------------
    // Project out weighted data onto histograms
    //-------------------------------------------
    
    TTree *newtree2 = dataw_signal->tree();
    newtree2->SetName("TreeWithWeight");
    newtree2->Write();
    
    
    fout->Write();

  }else{
    
    TString fileName = "Xib2JpsiXi.root";
    TFile *fout = new TFile(fileName,"RECREATE");
    TCut longCut = Xib_Cut && XiL_Cut;    

    TTree* treeNew1 = treeXi->CopyTree(longCut);
    
    //--------------------------------------------
    // Variables that are needed in the sWeighting
    //--------------------------------------------
    RooRealVar* L_TAU = new RooRealVar("L_TAU","Lambda decay time",-1.0e10,1.0e10);
    RooRealVar* L_P = new RooRealVar("L_P","Lambda mom",-1.0e10,1.0e10);
    RooRealVar* L_PT = new RooRealVar("L_PT","Lambda mom",-1.0e10,1.0e10);
    RooRealVar* Jpsi_pT = new RooRealVar("Jpsi_PT","Jpsi pT",-1.0e10,1.0e10);
    RooRealVar* Jpsi_p = new RooRealVar("Jpsi_P","Jpsi mom",-1.0e10,1.0e10);
    RooRealVar* L_FD_OWNPV = new RooRealVar("L_FD_OWNPV","Lambda FD",-1.0e10,1.0e10);
    RooRealVar* L_FD_ORIVX = new RooRealVar("L_FD_ORIVX","Lambda FD",-1.0e10,1.0e10);

    RooRealVar* Xi_TAU = new RooRealVar("Xi_TAU","Xi decay time",-1.0e10,1.0e10);
    RooRealVar* Xi_P = new RooRealVar("Xi_P","Xi mom",-1.0e10,1.0e10);
    RooRealVar* Xi_PT = new RooRealVar("Xi_PT","Xi mom",-1.0e10,1.0e10);
    RooRealVar* Xi_ENDVERTEX_CHI2 = new RooRealVar("Xi_ENDVERTEX_CHI2","Xi chisq_vtx",-1.0e10,1.0e10);
    RooRealVar* L_ENDVERTEX_Z = new RooRealVar("L_ENDVERTEX_Z","Lambda Z",-1.0e10,1.0e10);
    RooRealVar* Xi_DIRA_ORIVX = new RooRealVar("Xi_DIRA_ORIVX","Xi DIRA",-1.0e10,1.0e10);
    //RooRealVar* Xi_FDCHI2_ORIVX = new RooRealVar("Xi_DIRA_ORIVX","Xi DIRA",-1.0e10,1.0e10);

    
    RooRealVar* Xib_p = new RooRealVar("Xib_P","Xib mom",-1.0e10,1.0e10);
    RooRealVar* Xib_pT = new RooRealVar("Xib_PT","Xib pT",-1.0e10,1.0e10);
    RooRealVar* Xi_FD_OWNPV = new RooRealVar("Xi_FD_OWNPV","Xi FD",-1.0e10,1.0e10);
    RooRealVar* Xib_OWNPV_Z = new RooRealVar("Xib_OWNPV_Z","Xib origin Z",-1.0e10,1.0e10);


    RooRealVar* BachPi_pT = new RooRealVar("BachPi_PT","BachPi PT",-1.0e10,1.0e10);
    RooRealVar* BachPi_p = new RooRealVar("BachPi_P","BachPi P",-1.0e10,1.0e10);
    RooRealVar* BachPi_IPCHI2_OWNPV = new RooRealVar("BachPi_IPCHI2_OWNPV","BachPi chisqip",-1.0e10,1.0e10);


    RooArgSet args;
    args.add(*mvar2);
    args.add(*L_TAU);
    args.add(*L_P);
    args.add(*L_PT);
    args.add(*Jpsi_pT);
    args.add(*Jpsi_p);
    args.add(*L_FD_OWNPV);
    args.add(*L_FD_ORIVX);
    args.add(*Xi_TAU);
    args.add(*Xi_P);
    args.add(*Xi_PT);
    args.add(*Xib_p);
    args.add(*Xib_pT);
    args.add(*Xib_OWNPV_Z);
    args.add(*Xi_FD_OWNPV);
    args.add(*L_ENDVERTEX_Z);
    args.add(*Xi_ENDVERTEX_CHI2);
    args.add(*Xi_DIRA_ORIVX);
    args.add(*BachPi_pT);
    args.add(*BachPi_p);
    args.add(*BachPi_IPCHI2_OWNPV);

    RooDataSet* data_red = new RooDataSet("data_red","dataset with reduced vars",treeNew1,args);
    
    TString wsig = "nsig2_sw";
    
    RooArgSet *myVars1;
    pdf2->fitTo(*rdh2L,Hesse(kTRUE),Strategy(2));
    myVars1 = pdf2->getParameters(*rdh2L) ;
    myVars1->Print("v") ;
    myVars1->setAttribAll("Constant",kTRUE) ;
    nsig2->setConstant(kFALSE);
    ncomb2->setConstant(kFALSE);
    pdf2->fitTo(*rdh2L,Extended(),Hesse(kTRUE),Strategy(2));
    //--------------
    // Get sWeights
    //--------------
    cout << "Going to call SPlot now"<< endl;    
    RooStats::SPlot* sDataA1 = new RooStats::SPlot("sDataA1","An SPlot", *data_red, pdf2, RooArgList(*nsig2,*ncomb2) );
    cout << "Done calling SPlot now"<< endl;    
    
    data_red->Print("v") ; 
    gStyle->SetOptStat(011000);
    //--------------------------
    // create weighted data set 
    //--------------------------
    RooDataSet * dataw_signal = new RooDataSet(data_red->GetName(),data_red->GetTitle(),data_red,*data_red->get(),0,wsig) ;
    dataw_signal->Print("v") ; 
    //-------------------------------------------
    // Project out weighted data onto histograms
    //-------------------------------------------
    
    TTree *newtree2 = dataw_signal->tree();
    newtree2->SetName("TreeWithWeight");
    newtree2->Write();
    
    
    fout->Write();


  }
  
  


  
  return;
  






  
}
