#include <sys/stat.h>
#include <unistd.h>
#include <string>

using namespace std;
using namespace RooFit;

// add Xib- weight, w(eta)
bool applyEtaWeight = false;   
bool MakeWeightFile = false;
bool useNewStrip  = true;
int preScale = 100;

TString year = "Run1";
//TString year = "Run2";
TString modes[2];

bool mcFit = false;
bool applyHlt = true;

// Defaults
/*
bool applyJpsiLWeight = true;
bool applyLbKinWeight = true; // run with mcopt = 1, then set to true
bool applyLKinWeight = true; // run with mcopt = 2, then set to true
bool applyPrPiKinWeight = true;
bool applyJpsiXiWeight = true;
bool applyXiWeight = true;
*/

bool applyJpsiLWeight = false;
bool applyLbKinWeight = false; // run with mcopt = 1, then set to true
bool applyLKinWeight = false; // run with mcopt = 2, then set to true
bool applyPrPiKinWeight = false;
bool applyJpsiXiWeight = false;
bool applyXiWeight = false;

float LbMassWin = 8.0;

bool applyXiPiTrackingAcceptanceCut = false;

float XibLifetimeWeight = 0;

double tauLb = 1.470;
double err_tauLb = 0.010;

double tauXib = 1.571;
double err_tauXib = 0.041;

double genEffLb = 0.181;
double err_genEffLb = 0.003;
double genEffXib = 0.171;
double err_genEffXib = 0.001;

double lumJpsiL0 = 2962.0;
double lumJpsiXi = 2811.0;
double lumJpsiL0_2016 = 1551.63;
double lumJpsiXi_2016 = 1622.41;

TCut Jpsi_Cut;
TCut Lambda_Cut;
TCut LambdaLL_Cut;
TCut LambdaDD_Cut;
TCut Xi_Cut;
TCut XiL_Cut;
TCut XiD_Cut;

TCut MCXib_Cut;
TCut MCXib_Cut_NoWin;
TCut Xib_Cut;
TCut Lb_Cut;
TCut Lambdab_Cut;
TCut MCLambdab_Cut;
TCut LambdabLL_Cut;
TCut MCLambdabLL_Cut;
TCut MCLambdabLL_Cut_NoWin;
TCut LambdabDD_Cut;
TCut MCLambdabDD_Cut;
TCut MCLambdabDD_Cut_NoWin;
TCut XiWeight;


double mLowLb=5500;
double mHiLb = 5750;
double mLowXib= 5640;
double mHiXib = 5950;
//double mLowXib= 5700;
//double mHiXib = 5900;
float bwLb = 2.0;
float bwXib = 5.0;
int nbinLb, nbinXib;

TCut Lb_Trig;
TCut Xib_Trig;

TCut XibWt;
TCut BachPiWt;

TCut JpsiLWeight;
TCut JpsiXiWeight;


//float XibWtPar[2] = {6.33, 2.78};   // exponential correction
float XibWtPar[2] = {1.0, 0.3};   // linear correction
float BachPiPar[3] = {0.746, 1.30, 2.12};
float JpsiWeightPar[3];// = {0.75, 0.38, 0.0};

float JpsiXiWeightPar[3] = {0.47, 0.834, 0.0};  // was 0.6
//float JpsiXiWeightPar[3] = {0.0, 4.5, 0.0};
float XiWeightPar[3] = {1.59, -0.35, 0.029};

int nHistsLb;
int histLbf;
int nHistsL;
int histLf;
int nHistsp;
int histpf;
int nHistspi;
int histpif;
int nHistsR;
int histRf;
int histJpsif;
int nHistsJpsi;
int nHistsEx;
int histExf;
int nHistsSys;
int histSys;
int nHistsBachPi;
int histBachPif;
int nHistsXi;
int histXif;
int nHistsPV;
int histPV;
int nHistsXiVtx;
int histXiVtxf;


int nRatio = 0;
int nHists;
int nHists2D;
TString histNames[200];
TString histExp[200];
TString histExp2[200];
TString hTit[200];
int NB[200];
double LO[200];
double HI[200];

TH1F *hda[200];
TH1F *hmc[200];
TH1F *hr[20];

TH2F *whistLb ;
TH2F *whistL;
TH2F *whistPi ;
TH1F *whistPrPi;
TProfile *prPi;
TH2D *trackEff;


TTree* TreeWithWeight ;
TTree* newtree;
TTree *treeLL;
TTree *treeDD;
TTree *treeXi;

TString signalWeight;
TCut simCut;

TLegend* legend1;


TChain* chainl;
TChain* chainlt;
TChain* chainx;
TChain* chainxt;

int mcOpt;
bool makePDF;

TString parentString;
TString parentPart;
TString lambdaType;
TString tag;
TString tagCondLb;
TString tagCondXib;
TString tagCond;

double bottomMarginOffset = 0.13;

ULong64_t eventNumber;
double Jpsi_TRUEPT, Jpsi_TRUEP_Z, Jpsi_P, L_P, L_PT;
double L_TRUEPT, L_TRUEP_Z, pi_TRUEP_Z, pi_TRUEPT, p_TRUEP_Z, p_TRUEPT;
double J_psi_1S_TRUEPT,J_psi_1S_TRUEP_Z,  Lambda_b0_TRUETAU, Lambda0_TRUEP_Z, Lambda0_TRUEPT;
double piminus_TRUEP_Z, piminus_TRUEPT, pplus_TRUEPT, pplus_TRUEP_Z;
double pi_PT, pi_PZ, p_PT, p_PZ, ptPi, pPi, pPr, asym;      
double Lb_TRUEPT, Lb_TRUEP_Z;
double Lb_P, Lb_PT, Lb_TAU, Lb_TRUETAU;
double Lambda_b0_TRUEP_Z, Lambda_b0_TRUEPT; 
double Ximinus_TRUEPT;
double Xi_bminus_TRUEP_Z, Xi_bminus_TRUEPT, Xi_bminus_TRUETAU, piminus0_TRUEPT, piminus0_TRUEP_Z;
int BachPi_TRACK_Type;
double Xib_PT, Xib_P, Xib_TAU, Xib_TRUETAU, BachPi_P, BachPi_PT, Xib_TRUEPT,Xib_TRUEP_Z, Xi_PT, Xi_TRUEPT;
float mass, err_mass;

TH1F *h1;// = new TH1F("h1","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
TH1F *h1LL;// = new TH1F("h1LL","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
TH1F *h2;// = new TH1F("h2","#Xi_{b}^{0} mass",nbinXib,mLowXib,mHiXib);
TH1F *h2L;// = new TH1F("h2L","#Xi_{b}^{0} mass, Long",nbinXib,mLowXib,mHiXib);
TH1F *h2D;// = new TH1F("h2D","#Xi_{b}^{0} mass, Down",nbinXib,mLowXib,mHiXib);

RooRealVar* mvar1;
RooRealVar* mg1;
RooRealVar* sg1A;
RooRealVar* sg1scale;
RooFormulaVar* sg1B;
RooAbsPdf* gauss1A;
RooAbsPdf* gauss1B;
RooRealVar* fr1;
RooAddPdf* sig1;
RooRealVar* b1 ;
RooExponential* comb1;
RooRealVar* nsig1;
RooRealVar* ncomb1;
RooAbsPdf* pdf1;
RooDataHist *rdh1;
RooDataHist *rdh1LL;


RooRealVar* mvar2;
RooRealVar* mg2;
RooRealVar* sg2A;
RooRealVar* sg2scale;
RooFormulaVar* sg2B;
RooAbsPdf* sig2;
RooAbsPdf* gauss2A;
RooAbsPdf* gauss2B;
RooRealVar* fr2;
RooRealVar* alpha1A;
RooRealVar* alpha1B;
RooRealVar* alpha2A;
RooRealVar* alpha2B;
RooRealVar* nCB;
RooRealVar* b2;
RooExponential* comb2;
RooRealVar* nsig2;
RooRealVar* ncomb2;
RooAbsPdf* pdf2;

RooDataHist *rdh2;
RooDataHist *rdh2L;
RooDataHist *rdh2D;


double fitmass_Xib;
double err_fitmass_Xib;
double fitmassD_Xib;
double err_fitmassD_Xib;
double fitmassL_Xib;
double err_fitmassL_Xib;

double fitmass_Lb_LL;
double err_fitmass_Lb_LL;
double fitmass_Lb_DD;
double err_fitmass_Lb_DD;
double nsigLb;
double err_nsigLb;

double yieldXib, yieldXibL, yieldXibD;
double err_yieldXib, err_yieldXibL, err_yieldXibD;

///////////////////////////////////////////////////////////////////////////////

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

  RooAbsPdf* sig;// = sig1;
  if(mode==0){
    titleX = "#it{J/#psi#Lambda} mass [MeV/#it{c}^{2}]";
    legTitle = "#it{#Lambda_{b}^{0}}#rightarrow #it{J/#psi#Lambda}";
    sig = sig1;
  }else if(mode==1){
     titleX = "#it{J/#psi#Xi^{#font[122]{-}}} mass [MeV/#it{c}^{2}]";
     legTitle = "#it{#Xi_{b}^{#font[122]{-}}}#rightarrow #it{J/#psi#Xi^{#font[122]{-}}}";
     sig = sig2;
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
  plot2->SetMinimum(0.1);
  plot2->SetMaximum(1.2*plot2->GetMaximum());
  plot2->Draw();
  if(mcFit){
    c->SetLogy();
    sig->paramOn(plot2,RooFit::Layout(0.6,0.95,0.9),Format("NE",AutoPrecision(1)));
    //plot2->getAttText()->SetTextSize(0.4) ;
  }
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

  if(!mcFit){
    TLegend* legend1 = new TLegend(0.64,0.65,0.88,0.88);
    legend1->SetFillStyle(0);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->SetTextSize(0.05);
    legend1->AddEntry("Data","Data","LEP");
    legend1->AddEntry("Full","Full PDF","L");
    legend1->AddEntry("Signal",legTitle,"L");
    legend1->Draw();
  }
  
  return c;
}


void fillLbMass(){

  h1 = new TH1F("h1","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
  h1LL = new TH1F("h1LL","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
  h1->SetMinimum(0);
  h1LL->SetMinimum(0);


  TCanvas *cLam = new TCanvas("cLam","",800,400);
  cLam->Divide(2,1);
  cLam->cd(1);
  treeDD->Draw("Lb_DTF_M_JpsiLConstr>>h1",LambdabDD_Cut);
  cLam->cd(2);
  treeLL->Draw("Lb_DTF_M_JpsiLConstr>>h1LL",LambdabLL_Cut);


}

void fillXibMass(){

  h2 = new TH1F("h2","#Xi_{b}^{0} mass",nbinXib,mLowXib,mHiXib);
  h2L = new TH1F("h2L","#Xi_{b}^{0} mass, Long",nbinXib,mLowXib,mHiXib);
  h2D = new TH1F("h2D","#Xi_{b}^{0} mass, Down",nbinXib,mLowXib,mHiXib);
  h2->SetMinimum(0);
  h2L->SetMinimum(0);
  h2D->SetMinimum(0);

  TCanvas *cXib = new TCanvas("cXib","",600,600);
  cXib->Divide(2,2);
  cXib->cd(1);
  treeXi->Draw("Xib_DTF_M_JpsiXiLConstr>>h2",Xib_Cut);
  cXib->cd(2);
  treeXi->Draw("Xib_DTF_M_JpsiXiLConstr>>h2L",Xib_Cut&&XiL_Cut);  
  cXib->cd(3);
  treeXi->Draw("Xib_DTF_M_JpsiXiLConstr>>h2D",Xib_Cut&&XiD_Cut);
  cXib->Update();


}

void doXibFitAndPlot(){
  /*
  pdf2->fitTo(*rdh2,Hesse(kTRUE),Strategy(1));
  mass = mg2->getVal();
  err_mass = mg2->getError();
  TCanvas* cc2 = makeRooPlot(pdf2,rdh2,mvar2,1,nsig2,modes[1],mLowXib,mHiXib,nbinXib);
  cc2->Draw();
  cc2->Update();
  yieldXib = nsig2->getVal();
  err_yieldXib = nsig2->getError();
  fitmass_Xib = mg2->getVal();
  err_fitmass_Xib = mg2->getError();
  */

  pdf2->fitTo(*rdh2D,Hesse(kTRUE),Strategy(1));
  mass = mg2->getVal();
  err_mass = mg2->getError();
  TCanvas* cc2D = makeRooPlot(pdf2,rdh2D,mvar2,1,nsig2,modes[1]+"D",mLowXib,mHiXib,nbinXib);
  cc2D->Draw();
  cc2D->Update();
  yieldXibD = nsig2->getVal();
  err_yieldXibD = nsig2->getError();
  fitmassD_Xib = mg2->getVal();
  err_fitmassD_Xib = mg2->getError();

  pdf2->fitTo(*rdh2L,Hesse(kTRUE),Strategy(1));
  mass = mg2->getVal();
  err_mass = mg2->getError();
  TCanvas* cc2L = makeRooPlot(pdf2,rdh2L,mvar2,1,nsig2,modes[1]+"L",mLowXib,mHiXib,nbinXib);
  cc2L->Draw();
  cc2L->Update();
  yieldXibL = nsig2->getVal();
  err_yieldXibL = nsig2->getError();
  fitmassL_Xib = mg2->getVal();
  err_fitmassL_Xib = mg2->getError();  
  

  if(makePDF){
    TString suff = "XiLorD";
    //if(mcFit) suff = suff + "_MC";  
    //cc2->Print("./Plots/XibMass_LambdaDD_"+suff+".C");
    //cc2->Print("./Plots/XibMass_LambdaDD_"+suff+".pdf");
    //cc2->Print("./Plots/XibMass_LambdaDD_"+suff+".png");
    //cc2->Print("./Plots/XibMass_LambdaDD_"+suff+".eps");
    suff = "XiL";
    if(mcFit) suff = suff + "_MC";  
    cc2L->Print("./Plots/XibMass_LambdaDD_"+suff+".C");
    cc2L->Print("./Plots/XibMass_LambdaDD_"+suff+".pdf");
    cc2L->Print("./Plots/XibMass_LambdaDD_"+suff+".png");
    cc2L->Print("./Plots/XibMass_LambdaDD_"+suff+".eps");
    suff = "XiD";
    if(mcFit) suff = suff + "_MC";  
    cc2D->Print("./Plots/XibMass_LambdaDD_"+suff+".C");
    cc2D->Print("./Plots/XibMass_LambdaDD_"+suff+".pdf");
    cc2D->Print("./Plots/XibMass_LambdaDD_"+suff+".png");
    cc2D->Print("./Plots/XibMass_LambdaDD_"+suff+".eps");
  }



}


 
void fillLbMassMC(){

  h1 = new TH1F("h1","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
  h1LL = new TH1F("h1LL","#Lambda_{b}^{0} mass",nbinLb,mLowLb,mHiLb);
  h1->SetMinimum(0);
  h1LL->SetMinimum(0);


  TCanvas *cLam = new TCanvas("cLam","",800,400);
  cLam->Divide(2,1);
  cLam->cd(1);
  chainl->Draw("Lb_DTF_M_JpsiLConstr>>h1",MCLambdabDD_Cut_NoWin);
  cLam->cd(2);
  chainl->Draw("Lb_DTF_M_JpsiLConstr>>h1LL",MCLambdabLL_Cut_NoWin);


}

void fillXibMassMC(){

  h2 = new TH1F("h2","#Xi_{b}^{0} mass",nbinXib,mLowXib,mHiXib);
  h2L = new TH1F("h2L","#Xi_{b}^{0} mass, Long",nbinXib,mLowXib,mHiXib);
  h2D = new TH1F("h2D","#Xi_{b}^{0} mass, Down",nbinXib,mLowXib,mHiXib);
  h2->SetMinimum(0);
  h2L->SetMinimum(0);
  h2D->SetMinimum(0);

  TCanvas *cXib = new TCanvas("cXib","",600,600);
  cXib->Divide(2,2);
  cXib->cd(1);
  chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>h2",Xib_Cut);
  cXib->cd(2);
  chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>h2L",Xib_Cut&&XiL_Cut);  
  cXib->cd(3);
  chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>h2D",Xib_Cut&&XiD_Cut);
  cXib->Update();


}

void getLbMassFitFunctions(){

  mvar1 = new RooRealVar("Lb_DTF_M_JpsiLConstr","",mLowLb,mHiLb);
  mg1 = new RooRealVar("m_{1}","",5620,5610,5630);
  sg1A = new RooRealVar("#sigma_{1A}","",7.56,2,12);
  sg1scale = new RooRealVar("#sigma_{1B}/#sigma_{1A}","",2.32,0.5,3.5);
  sg1B = new RooFormulaVar("sg1B","","@0*@1",RooArgList(*sg1A,*sg1scale));
  alpha1A = new RooRealVar("#alpha_{1A}","",1.5,0.5,3.0);
  alpha1B = new RooRealVar("#alpha_{1B}","",-2.0,-5.0,-0.2);
  nCB = new RooRealVar("nCB","",10);
  gauss1A = new RooCBShape("gauss1A","Signal",*mvar1,*mg1,*sg1A,*alpha1A,*nCB);
  gauss1B = new RooCBShape("gauss1B","Signal",*mvar1,*mg1,*sg1B,*alpha1B,*nCB);
  //gauss1A = new RooGaussian("gauss1A","Signal",*mvar1,*mg1,*sg1A);
  //gauss1B = new RooGaussian("gauss1B","Signal",*mvar1,*mg1,*sg1B);
  fr1 = new RooRealVar("f_{1A}","fraction",0.5);//,0.3,0.7);
  sig1 = new RooAddPdf("sig1","Signal",RooArgList(*gauss1A,*gauss1B),RooArgList(*fr1));
  b1  = new RooRealVar("b_{1}", "b1", -0.0003,-0.05,0.05);
  comb1 = new RooExponential("comb1","f_{back}",*mvar1,*b1);
  nsig1 = new RooRealVar("n_{sig1}","",12000,0,30000000);
  ncomb1 = new RooRealVar("n_{comb1}","",6000,0,30000000);
  pdf1 = new RooAddPdf("pdf1","",RooArgList(*sig1,*comb1),RooArgList(*nsig1,*ncomb1));
  rdh1 = new RooDataHist("rdh1","", *mvar1, h1);
  rdh1LL = new RooDataHist("rdh1LL","", *mvar1, h1LL);
}

void getXibMassFitFunctions(){

  mvar2 = new RooRealVar("Xib_DTF_M_JpsiXiLConstr","",mLowXib,mHiXib);
  mg2 = new RooRealVar("m_{2}","",5797,5780,5810);
  sg2A = new RooRealVar("#sigma_{2A}","",sg1A->getVal(),3,12);
  sg2scale = new RooRealVar("#sigma_{2B}/#sigma_{2A}","",sg1scale->getVal(),1.0,3.5);
  sg2B = new RooFormulaVar("sg2B","","@0*@1",RooArgList(*sg2A,*sg2scale));
  alpha2A = new RooRealVar("#alpha_{2A}","",1.5,0.5,3.0);
  alpha2B = new RooRealVar("#alpha_{2B}","",-2.0,-5.0,-0.2);
  gauss2A = new RooCBShape("gauss2A","Signal",*mvar2,*mg2,*sg2A,*alpha2A,*nCB);
  gauss2B = new RooCBShape("gauss2B","Signal",*mvar2,*mg2,*sg2B,*alpha2B,*nCB);
  //gauss2A = new RooGaussian("gauss2A","Signal",*mvar2,*mg2,*sg2A);
  //gauss2B = new RooGaussian("gauss2B","Signal",*mvar2,*mg2,*sg2B);
  fr2 = new RooRealVar("f_{2A}","fraction",0.5);//,0.3,0.7);//8,0.2,1.0);
  sig2 = new RooAddPdf("sig2","Signal",RooArgList(*gauss2A,*gauss2B),RooArgList(*fr2));
  //sig2 = new RooGaussian("sig2","Signal",*mvar2,*mg2,*sg2A);
  b2  = new RooRealVar("b2", "b2", -0.0003,-2,2);
  comb2 = new RooExponential("comb2","f_{back}",*mvar2,*b2);
  nsig2 = new RooRealVar("n_{sig2}","",12000,0,30000000);
  ncomb2 = new RooRealVar("n_{comb2}","",6000,0,30000000);
  pdf2 = new RooAddPdf("pdf2","",RooArgList(*sig2,*comb2),RooArgList(*nsig2,*ncomb2));

  rdh2 = new RooDataHist("rdh2","", *mvar2, h2);
  rdh2L = new RooDataHist("rdh2L","", *mvar2, h2L);
  rdh2D = new RooDataHist("rdh2D","", *mvar2, h2D);


}


inline bool exists(TString name) {
  struct stat buffer;   
  return (stat (name, &buffer) == 0); 

}

long GetFileSize(TString filename)
{
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    return rc == 0 ? stat_buf.st_size : -1;
}


double get2DKinWt(TH2F* h, double pt, double p, int iopt = 0){
  // iopt = 0 --> provide pt, pz
  // iopt = 1 --> provide pt, pmag 


  double eta;
  double w = 1.0;
  if(p<=0 || pt<=0) return 1.0;
  if(iopt==0) eta = -log(tan(0.5*atan(pt/p)));
  if(iopt==1) eta = -log(tan(0.5*asin(pt/p)));
  int ib = h->FindBin(0.001*pt,eta);
  w = h->GetBinContent(ib);
  if(w<=0 || w>5) w = 1.0;

  return w;


}

double get1DKinWt(TH1F* h, double x){

  
  double w = 1.0;
  int ib = h->FindBin(x);
  w = h->GetBinContent(ib);
  if(w<=0 || w>5) w = 1.0;

  return w;
}

void makePlots(TCanvas *ccc, TString filename){

  if(parentPart=="Lb") tagCond = tagCondLb;
  if(parentPart=="Xib") tagCond = tagCondXib;

//if( !applyLbKinWeight ) tagCond = "_unweighted";
  
  ccc->Print("./Plots/"+filename+tagCond+".pdf");
  ccc->Print("./Plots/"+filename+tagCond+".png");
  ccc->Print("./Plots/"+filename+tagCond+".eps");
  ccc->Print("./Plots/"+filename+tagCond+".C");
  
}

void setConditions(){

  if(mcOpt==1){
    applyJpsiLWeight = false;
    applyLbKinWeight = false; 
    applyLKinWeight = false; 
    applyPrPiKinWeight = false;
    applyJpsiXiWeight = false;
    applyXiWeight = false;
  }else if(mcOpt==2){
    applyJpsiLWeight = false;
    applyLbKinWeight = true; 
    applyLKinWeight = false; 
    applyPrPiKinWeight = false;
    applyJpsiXiWeight = false;
    applyXiWeight = false;
  }else if(mcOpt==3){
    applyJpsiLWeight = true;
    applyLbKinWeight = true; 
    applyLKinWeight = false; 
    applyPrPiKinWeight = false;
    applyJpsiXiWeight = true;
    applyXiWeight = false;
  }else if(mcOpt>=4 ){
    applyJpsiLWeight = true;
    applyLbKinWeight = true; 
    applyLKinWeight = false; 
    applyPrPiKinWeight = true;
    applyJpsiXiWeight = true;
    applyXiWeight = false;
  }
  
  parentString = "#Lambda_{b}^{0}";
  if(parentPart=="Xib")  parentString = "#Xi_{b}^{#font[122]{-}}";

}


void getTagCond(){

  int ieta = 0;
  int iBach = 0;
  int iJpsiL = 0;
  int iLbKin = 0;
  int iLKin = 0;
  int iPrPiKin = 0;
  int iJpsiXi = 0;
  int iXi = 0;
  if(applyEtaWeight) ieta = 1;
  if(applyJpsiLWeight) iJpsiL = 1;
  if(applyLbKinWeight) iLbKin = 1;
  if(applyLKinWeight) iLKin = 1;
  if(applyPrPiKinWeight) iPrPiKin = 1;
  if(applyJpsiXiWeight) iJpsiXi = 1;
  if(applyXiWeight) iXi = 1;
  
  int iopt = mcOpt;

  TString tagCondLb = Form("_Lb%d_Jpsi%d_Lam%d_PrPi%d_mcOpt%d",iLbKin,iJpsiL,iLKin,iPrPiKin,iopt);
  TString tagCondXib = Form("_Lb%d_Jpsi%d_JpsiXi%d_Lam%d_PrPi%d_XiWt%d_XiEta%d_mcOpt%d",iLbKin,iJpsiL,iJpsiXi,iLKin,iPrPiKin,iXi,ieta,iopt);
  tagCondLb = tagCondLb+"_"+year;
  tagCondXib = tagCondXib+"_"+year;
  

  cout << "tagCond_Lb = " << tagCondLb << endl;
  cout << "tagCond_Xib = " << tagCondXib << endl;
  

}


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

void getCuts(){

  if(year=="Run1"){
    JpsiWeightPar[0] = 0.75;// = {0.75, 0.38, 0.0};
    JpsiWeightPar[1] = 0.38;// = {0.75, 0.38, 0.0};
    JpsiWeightPar[2] = 0.0;// = {0.75, 0.38, 0.0};
  }else  if(year=="Run2"){
    JpsiWeightPar[0] = 0.92;// = {0.75, 0.38, 0.0};
    JpsiWeightPar[1] = 0.12;// = {0.75, 0.38, 0.0};
    JpsiWeightPar[2] = 0.0;// = {0.75, 0.38, 0.0};
  }

  
  
  //TCut Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Lb_L0Global_TIS>0)";
  //Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Xib_L0Global_TIS>0)";
  Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Jpsi_L0MuonDecision_TOS>0)";
  Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Jpsi_L0MuonDecision_TOS>0)";
  //Lb_Trig = "( Lb_L0Global_TIS>0>0 && Jpsi_L0DiMuonDecision_TOS==0 )";
  //Xib_Trig = "( Xib_L0Global_TIS>0 && Jpsi_L0DiMuonDecision_TOS==0 )";
  if(applyHlt){
    //Lb_Trig = Lb_Trig && "Jpsi_Hlt1Phys_TOS>0 && Jpsi_Hlt2Phys_TOS>0";  // def
    //Xib_Trig = Xib_Trig && "Jpsi_Hlt1Phys_TOS>0 && Jpsi_Hlt2Phys_TOS>0";  //def
    Lb_Trig = Lb_Trig && "(Jpsi_Hlt1TrackMuonDecision_TOS>0 || Jpsi_Hlt1DiMuonHighMassDecision_TOS>0)";
    Xib_Trig = Xib_Trig && "(Jpsi_Hlt1TrackMuonDecision_TOS>0 || Jpsi_Hlt1DiMuonHighMassDecision_TOS>0)";
    Lb_Trig = Lb_Trig && "Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS>0";
    Xib_Trig = Xib_Trig && "Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS>0";
    //Lb_Trig = Lb_Trig && "(Jpsi_Hlt1SingleMuonNoIPDecision_TOS>0 || Jpsi_Hlt1TrackMuonDecision_TOS>0 || Jpsi_Hlt1DiMuonHighMassDecision_TOS>0 )";
    // Xib_Trig = Xib_Trig && "(Jpsi_Hlt1SingleMuonNoIPDecision_TOS>0 || Jpsi_Hlt1TrackMuonDecision_TOS>0 || Jpsi_Hlt1DiMuonHighMassDecision_TOS>0 )";
  }
  

  
  // default mass window on Xi = 11  MeV

  Jpsi_Cut = "abs(Jpsi_MM-3098)<40";// && Jpsi_L0DiMuonDecision_TOS>0";
  Jpsi_Cut = Jpsi_Cut && "muplus_IPCHI2_OWNPV>4&&muminus_IPCHI2_OWNPV>4&&muplus_PT>550&&muminus_PT>550&&Jpsi_ENDVERTEX_CHI2<16 && Jpsi_FDCHI2_OWNPV>12";
  Lambda_Cut  = "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100 && L_FD_OWNPV/L_ENDVERTEX_ZERR>5 && L_TAU>0";
  TString lamWin = Form("abs(L_M-1115.7)<%f",LbMassWin);
  TCut cutLamMass = lamWin;
  Lambda_Cut = Lambda_Cut && cutLamMass;
  LambdaLL_Cut  = Lambda_Cut && "p_TRACK_Type==3 && pi_TRACK_Type==3";
  LambdaDD_Cut  = Lambda_Cut && "p_TRACK_Type==5 && pi_TRACK_Type==5";

  //Lambda_Cut = Lambda_Cut && "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100&&L_FD_OWNPV/L_ENDVERTEX_ZERR>5"; 
  //LambdaLL_Cut = LambdaLL_Cut && "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100&&L_FD_OWNPV/L_ENDVERTEX_ZERR>5"; 

  Xi_Cut = "abs(Xi_M-L_M+1115.7-1321.9)<10 && Xi_DIRA_ORIVX>0.995 && Xi_ENDVERTEX_CHI2<25 && Xi_TAU>0"; // Use Long Track for pion from Xi- decay
  //XiL_Cut = "BachPi_TRACK_Type==3 && BachPi_P>5000"; // Use Long Track for pion from Xi- decay
  //XiD_Cut = "BachPi_TRACK_Type==5 && BachPi_P>5000"; // Use DownStream Track for pion from Xi- decay
  XiL_Cut = "BachPi_TRACK_Type==3&&BachPi_IPCHI2_OWNPV>16"; // Use Long Track for pion from Xi- decay
  XiD_Cut = "BachPi_TRACK_Type==5"; // Use DownStream Track for pion from Xi- decay

  // Apply p and eta cut, for comparison to region where tracking corrections are valid
  if(applyXiPiTrackingAcceptanceCut) XiL_Cut = XiL_Cut && "BachPi_P>5000 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))>1.9 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))<4.9";
  

  //Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_MINIPCHI2<12 && Xib_TAU*1000>0.2";
  Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_IPCHI2_OWNPV<12 && Xib_TAU*1000>0.2";
  Xib_Cut = Xib_Cut && Xib_Trig && Jpsi_Cut && LambdaDD_Cut && Xi_Cut && (XiL_Cut||XiD_Cut) ;
  float ml = mLowXib;
  float mh = mHiXib;
  TString mCut = Form("Xib_DTF_M_JpsiXiLConstr>%f && Xib_DTF_M_JpsiXiLConstr<%f",ml,mh);
  Xib_Cut = Xib_Cut && mCut;
  MCXib_Cut_NoWin = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<300)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<0))";
  //MCXib_Cut_NoWin = Xib_Cut && "(Xib_BKGCAT<20)";
  MCXib_Cut = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25))";

  Lb_Cut = "(Lb_ENDVERTEX_CHI2/Lb_ENDVERTEX_NDOF < 10) && L_TAU>0 && Lb_MINIPCHI2 < 12 && Lb_TAU*1000>0.2";
  LambdabDD_Cut = Jpsi_Cut && LambdaDD_Cut && Lb_Cut && Lb_Trig;
  LambdabLL_Cut = Jpsi_Cut && LambdaLL_Cut && Lb_Cut && Lb_Trig;

  ml = mLowLb;
  mh = mHiLb;
  mCut = Form("Lb_DTF_M_JpsiLConstr>%f&& Lb_DTF_M_JpsiLConstr<%f",ml,mh);
  LambdabDD_Cut = LambdabDD_Cut && mCut;
  LambdabLL_Cut = LambdabLL_Cut && mCut;

  MCLambdabDD_Cut_NoWin = LambdabDD_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<300)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<0))";
  MCLambdabLL_Cut_NoWin = LambdabLL_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<300)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<0))";

  MCLambdabDD_Cut = LambdabDD_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";
  MCLambdabLL_Cut = LambdabLL_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";

  //TString ss = Form("(1.0*exp(-0.5*((-log(tan(0.5*asin(Xib_PT/Xib_P))))-%4.2f)**2/%4.2f**2))",XibWtPar[0], XibWtPar[1]); // exponential
  TString ss = Form("%4.2f+%4.2f*(-log(tan(0.5*asin(Xib_PT/Xib_P))))",XibWtPar[0], XibWtPar[1]);
  XibWt = ss;//"(1.0*exp(-0.5*((-log(tan(0.5*asin(Xib_PT/Xib_P))))-6.33)**2/2.78**2))";
  ss = Form("(%6.3f+%6.3f*exp(-0.001*BachPi_P/%6.3f))",BachPiPar[0],BachPiPar[1],BachPiPar[2]);
  BachPiWt = ss;
  ss = Form("(%5.3f+%5.3f*(Jpsi_P/Lb_P)+%5.3f*(Jpsi_P/Lb_P)**2)",JpsiWeightPar[0],JpsiWeightPar[1],JpsiWeightPar[2]);
  JpsiLWeight = ss;
  ss = Form("(%5.3f+%5.3f*(Jpsi_P/Xib_P)+%5.3f*(Jpsi_P/Xib_P)**2)",JpsiXiWeightPar[0],JpsiXiWeightPar[1],JpsiXiWeightPar[2]);
  JpsiXiWeight = ss;
  ss = Form("(%6.3f+%6.3f*(Xi_PT/1000.0)+%6.3f*(Xi_PT/1000.0)**2)",XiWeightPar[0],XiWeightPar[1],XiWeightPar[2]);
  XiWeight = ss;//"(2.6-1.09*(Xi_PT/1000.0)+0.14*(Xi_PT/1000.)**2)";

  setConditions();
  getTagCond();
  

}

void getLooseCuts(){

  Jpsi_Cut = "abs(Jpsi_MM-3098)<80";// && Jpsi_L0DiMuonDecision_TOS>0";
  Jpsi_Cut = Jpsi_Cut && "muplus_IPCHI2_OWNPV>4&&muminus_IPCHI2_OWNPV>4&&muplus_PT>550&&muminus_PT>550&&Jpsi_ENDVERTEX_CHI2<16 && Jpsi_FDCHI2_OWNPV>12";
  Lambda_Cut  = "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100 && L_FD_OWNPV/L_ENDVERTEX_ZERR>5 && abs(L_M-1115.7)<20 && L_TAU>0";
  LambdaLL_Cut  = Lambda_Cut && "p_TRACK_Type==3 && pi_TRACK_Type==3";
  LambdaDD_Cut  = Lambda_Cut && "p_TRACK_Type==5 && pi_TRACK_Type==5";


  Xi_Cut = "abs(Xi_M-L_M+1115.7-1321.9)<30 && Xi_DIRA_ORIVX>0.995 && Xi_ENDVERTEX_CHI2<25 && Xi_TAU>0"; // Use Long Track for pion from Xi- decay
  XiL_Cut = "BachPi_TRACK_Type==3"; // Use Long Track for pion from Xi- decay
  XiD_Cut = "BachPi_TRACK_Type==5"; // Use DownStream Track for pion from Xi- decay

  // Apply p and eta cut, for comparison to region where tracking corrections are valid
  if(applyXiPiTrackingAcceptanceCut) XiL_Cut = XiL_Cut && "BachPi_P>5000 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))>1.9 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))<4.9";
  

  Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_IPCHI2_OWNPV<16 && Xib_TAU*1000>0.2";
  Xib_Cut = Xib_Cut && Xib_Trig && Jpsi_Cut && LambdaDD_Cut && Xi_Cut && (XiL_Cut||XiD_Cut) ;
  float ml = mLowXib;
  float mh = mHiXib;
  TString mCut = Form("Xib_DTF_M_JpsiXiLConstr>%f && Xib_DTF_M_JpsiXiLConstr<%f",ml-50,mh+50);
  Xib_Cut = Xib_Cut && mCut;
  MCXib_Cut = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25))";

  Lb_Cut = "(Lb_ENDVERTEX_CHI2/Lb_ENDVERTEX_NDOF < 10) && L_TAU>0 && Lb_MINIPCHI2 < 12 && Lb_TAU*1000>0.2";
  LambdabDD_Cut = Jpsi_Cut && LambdaDD_Cut && Lb_Cut && Lb_Trig;
  LambdabLL_Cut = Jpsi_Cut && LambdaLL_Cut && Lb_Cut && Lb_Trig;
  Lambdab_Cut = Jpsi_Cut && Lb_Cut && Lb_Trig && Lambda_Cut;

  ml = mLowLb;
  mh = mHiLb;
  mCut = Form("Lb_DTF_M_JpsiLConstr>%f&& Lb_DTF_M_JpsiLConstr<%f",ml-50,mh+50);
  LambdabDD_Cut = LambdabDD_Cut && mCut;
  LambdabLL_Cut = LambdabLL_Cut && mCut;
  Lambdab_Cut = Lambdab_Cut && mCut;

  MCLambdabDD_Cut = LambdabDD_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";
  MCLambdabLL_Cut = LambdabLL_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";

  

}


void bookLbBaseHistos(){

  signalWeight = "nsig1_sw";
  
 TString HB = "#Lambda";

 int i1 = 0;
 int i2 = 0;
 
 double maxDecayTime = 1.0;
 double maxZ = 2600;
 if(lambdaType=="LL"){
   maxDecayTime = 0.1;
   maxZ = 800;
 }
 

 histLf = 0;
 // Lambda0
 histNames[i1] = "taulam"; histExp[i1] = "L_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=50; 
 LO[i1]=0.0; HI[i1]=maxDecayTime; i1++;
 histNames[i1] = "plam"; histExp[i1] = "0.001*L_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
 LO[i1]=0.0; HI[i1]=120.0; i1++;
 histNames[i1] = "ptlam"; histExp[i1] = "0.001*L_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=50; 
 LO[i1]=0.0; HI[i1]=10.0; i1++;
 histNames[i1] = "fdlam"; histExp[i1] = "L_FD_OWNPV"; hTit[i1] = HB+" flight distance [mm]"; NB[i1]=26; 
 LO[i1]=0.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "vzlam"; histExp[i1] = "L_ENDVERTEX_Z"; hTit[i1] = HB+" z decay position [mm]"; NB[i1]=26; 
 LO[i1]=-100.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "etalam"; histExp[i1] = "-log(tan(0.5*asin(L_PT/L_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=14; 
 LO[i1]=1.5; HI[i1]=5.0; i1++;
 nHistsL = i1;

 if(parentPart=="Lb"){ 
   HB = "#Lambda_{b}^{0}";
   histLbf = i1;
   histNames[i1] = "plamb"; histExp[i1] = "0.001*Lb_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
   LO[i1]=0.0; HI[i1]=300.0; i1++;
   histNames[i1] = "ptlamb"; histExp[i1] = "0.001*Lb_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=20.0; i1++;
   histNames[i1] = "etalamb"; histExp[i1] = "-log(tan(0.5*asin(Lb_PT/Lb_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=20; 
   LO[i1]=1.5; HI[i1]=6.5; i1++;
 }else{
   signalWeight = "nsig2_sw";
   HB = "#Xi_{b}^{#font[122]{-}}";
   histLbf = i1;
   histNames[i1] = "pxib"; histExp[i1] = "0.001*Xib_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
   LO[i1]=0.0; HI[i1]=300.0; i1++;
   histNames[i1] = "ptxib"; histExp[i1] = "0.001*Xib_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=20.0; i1++;
   histNames[i1] = "etaxib"; histExp[i1] = "-log(tan(0.5*asin(Xib_PT/Xib_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=20; 
   LO[i1]=1.5; HI[i1]=6.5; i1++;
 }
 

 nHistsLb = i1 - histLbf;

 HB = "proton (from #Lambda)";
 histpf = i1;
 histNames[i1] = "pp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2+p_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=120.0; i1++;
 histNames[i1] = "ptp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 nHistsp = i1 - histpf;
 HB = "pion (from #Lambda)";
 histpif = i1;
 histNames[i1] = "ppi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2+pi_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=30; 
 LO[i1]=0.0; HI[i1]=30.0; i1++;
 histNames[i1] = "ptpi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=2.0; i1++;
 histNames[i1] = "ppiasym"; histExp[i1] = "(sqrt(p_PX**2+p_PY**2+p_PZ**2)-sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))/(sqrt(p_PX**2+p_PY**2+p_PZ**2)+sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))"; hTit[i1] = "p-#pi^{-} momentum asymmetry "; NB[i1]=20; LO[i1]=0.5; HI[i1]=1.0; i1++;
 nHists = i1;
 nHistspi = i1 - histpif ;

 histRf = i1;
 histNames[i1] = "jpsioverLb"; histExp[i1] = "Jpsi_P/Lb_P"; hTit[i1] = "p_{J/#psi} / p_{#Lambda_{b}}"; NB[i1]=12; 
 LO[i1]=0.25; HI[i1]=1.0; i1++;
 histNames[i1] = "LoverLb"; histExp[i1] = "L_P/Lb_P"; hTit[i1] = "p_{#Lambda} / p_{#Lambda_{b}}"; NB[i1]=25; 
 LO[i1]=0.0; HI[i1]=0.75; i1++;

 nHistsR = i1 - histRf ;
 
 histJpsif = i1;
 histNames[i1] = "jpsipt"; histExp[i1] = "0.001*Jpsi_PT"; hTit[i1] = "J/#psi p_{T} [GeV/#it{c}]"; NB[i1]=36; 
 LO[i1]=0.0; HI[i1]=18.0; i1++;
 histNames[i1] = "jpsieta"; histExp[i1] = "-log(tan(0.5*asin(Jpsi_PT/Jpsi_P)))"; hTit[i1] = "J/#psi #eta"; NB[i1]=18; 
 LO[i1]=1.5; HI[i1]=6.0; i1++;
 histNames[i1] = "muppt"; histExp[i1] = "0.001*muplus_PT"; hTit[i1] = "#mu^{+} p_{T} [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=10.0; i1++;
 histNames[i1] = "mumpt"; histExp[i1] = "0.001*muminus_PT"; hTit[i1] = "#mu^{#font[122]{-}} p_{T} [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=10.0; i1++;
 nHistsJpsi = i1 - histJpsif ;


 histPV = i1;
 HB = "#Lambda_{b}^{0}";
 histNames[i1] = "pvz"; histExp[i1] = "Lb_OWNPV_Z"; hTit[i1] = HB+" PV"; NB[i1]=60; 
 LO[i1]=-300.0; HI[i1]=300.0; i1++;
 nHistsPV = i1 - histPV;

 

 if(lambdaType=="LL"){
   HB = "proton (from #Lambda)";
   histExf = i1;
   histNames[i1] = "pIP"; histExp[i1] = "p_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=4.0; i1++;
   histNames[i1] = "pIPCHI2"; histExp[i1] = "log(p_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "pDOCAz"; histExp[i1] = "abs(p_PY*(pi_ORIVX_X-0.459)-p_PX*(pi_ORIVX_Y+0.014))/sqrt(p_PX**2+p_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=30; LO[i1]=0.0; HI[i1]=3.0; i1++;

   HB = "pion (from #Lambda)";
   histNames[i1] = "piIP"; histExp[i1] = "pi_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=10.0; i1++;
   histNames[i1] = "piIPCHI2"; histExp[i1] = "log(pi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "piDOCAz"; histExp[i1] = "abs(pi_PY*(pi_ORIVX_X-0.459)-pi_PX*(pi_ORIVX_Y+0.014))/sqrt(pi_PX**2+pi_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=32; LO[i1]=0.0; HI[i1]=8.0; i1++;
   nHistsEx = i1 - histExf;
 } 


 // Do these last
 HB = "proton (from #Lambda)";
 histSys = i1;
 histNames[i1] = "pIPCHI2_sys"; histExp[i1] = "p_IPCHI2_OWNPV"; hTit[i1] = HB+" #chi^{2}_{IP}"; NB[i1]=100; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 HB = "pion (from #Lambda)";
 histNames[i1] = "piIPCHI2_sys"; histExp[i1] = "pi_IPCHI2_OWNPV"; hTit[i1] = HB+" #chi^{2}_{IP}"; NB[i1]=100; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 nHistsSys = i1 - histSys;
 

 nHists = i1;

 cout << "# Histograms = " << nHists << endl;
 
 
 // Book the histograms
 for(int i=0; i<nHists; i++){
   hda[i] = new TH1F(histNames[i]+"_da",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i] = new TH1F(histNames[i]+"_mc",hTit[i],NB[i],LO[i],HI[i]);
 }


}

void bookXibBaseHistos(){

  signalWeight = "nsig2_sw";
  
 TString HB = "#Lambda";

 int i1 = 0;
 int i2 = 0;
 
 double maxDecayTime = 0.6;
 double maxZ = 2600;
 if(lambdaType=="LL"){
   maxDecayTime = 0.1;
   maxZ = 800;
 }
 

 histLf = 0;
 // Lambda0
 histNames[i1] = "taulam"; histExp[i1] = "L_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=12; 
 LO[i1]=0.0; HI[i1]=maxDecayTime; i1++;
 histNames[i1] = "plam"; histExp[i1] = "0.001*L_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 histNames[i1] = "ptlam"; histExp[i1] = "0.001*L_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=12; 
 LO[i1]=0.0; HI[i1]=6.0; i1++;
 histNames[i1] = "fdlam"; histExp[i1] = "L_FD_OWNPV"; hTit[i1] = HB+" flight distance [mm]"; NB[i1]=13; 
 LO[i1]=0.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "vzlam"; histExp[i1] = "L_ENDVERTEX_Z"; hTit[i1] = HB+" z decay position [mm]"; NB[i1]=13; 
 LO[i1]=-100.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "etalam"; histExp[i1] = "-log(tan(0.5*asin(L_PT/L_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=14; 
 LO[i1]=1.5; HI[i1]=5.0; i1++;
 nHistsL = i1;

 HB = "#Xi_{b}^{#font[122]{-}}";
 histLbf = i1;
 histNames[i1] = "pxib"; histExp[i1] = "0.001*Xib_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=300.0; i1++;
 histNames[i1] = "ptxib"; histExp[i1] = "0.001*Xib_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=20.0; i1++;
 histNames[i1] = "etaxib"; histExp[i1] = "-log(tan(0.5*asin(Xib_PT/Xib_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=10; 
 LO[i1]=1.5; HI[i1]=6.5; i1++;
 

 nHistsLb = i1 - histLbf;

 HB = "proton (from #Lambda)";
 histpf = i1;
 histNames[i1] = "pp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2+p_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=80.0; i1++;
 histNames[i1] = "ptp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=15; 
 LO[i1]=0.0; HI[i1]=6.0; i1++;
 nHistsp = i1 - histpf;
 HB = "pion (from #Lambda)";
 histpif = i1;
 histNames[i1] = "ppi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2+pi_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=20.0; i1++;
 histNames[i1] = "ptpi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=15; 
 LO[i1]=0.0; HI[i1]=1.5; i1++;
 histNames[i1] = "ppiasym"; histExp[i1] = "(sqrt(p_PX**2+p_PY**2+p_PZ**2)-sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))/(sqrt(p_PX**2+p_PY**2+p_PZ**2)+sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))"; hTit[i1] = "p-#pi^{-} momentum asymmetry "; NB[i1]=10; LO[i1]=0.5; HI[i1]=1.0; i1++;
 nHistspi = i1 - histpif ;

 histRf = i1;
 histNames[i1] = "jpsioverLb"; histExp[i1] = "Jpsi_P/Xib_P"; hTit[i1] = "p_{J/#psi} / p_{#Xi_{b}}"; NB[i1]=6; 
 LO[i1]=0.25; HI[i1]=1.0; i1++;
 histNames[i1] = "LoverLb"; histExp[i1] = "L_P/Xib_P"; hTit[i1] = "p_{#Lambda} / p_{#Xi_{b}}"; NB[i1]=6; 
 LO[i1]=0.0; HI[i1]=0.75; i1++;

 nHistsR = i1 - histRf ;
 
 histJpsif = i1;
 histNames[i1] = "jpsipt"; histExp[i1] = "0.001*Jpsi_PT"; hTit[i1] = "J/#psi p_{T} [GeV/#it{c}]"; NB[i1]=18; 
 LO[i1]=0.0; HI[i1]=18.0; i1++;
 nHistsJpsi = i1 - histJpsif ;

 if(lambdaType=="LL"){
   HB = "proton (from #Lambda)";
   histExf = i1;
   histNames[i1] = "pIP"; histExp[i1] = "p_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=4.0; i1++;
   histNames[i1] = "pIPCHI2"; histExp[i1] = "log(p_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "pDOCAz"; histExp[i1] = "abs(p_PY*(pi_ORIVX_X-0.459)-p_PX*(pi_ORIVX_Y+0.014))/sqrt(p_PX**2+p_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=30; LO[i1]=0.0; HI[i1]=3.0; i1++;

   HB = "pion (from #Lambda)";
   histNames[i1] = "piIP"; histExp[i1] = "pi_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=10.0; i1++;
   histNames[i1] = "piIPCHI2"; histExp[i1] = "log(pi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "piDOCAz"; histExp[i1] = "abs(pi_PY*(pi_ORIVX_X-0.459)-pi_PX*(pi_ORIVX_Y+0.014))/sqrt(pi_PX**2+pi_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=32; LO[i1]=0.0; HI[i1]=8.0; i1++;
   nHistsEx = i1 - histExf;
 }

   
 HB = "#pi^{#font[122]{-}} (from #Xi^{#font[122]{-}})";
 histBachPif = i1; 
 histNames[i1] = "ppib"; histExp[i1] = "0.001*BachPi_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=18; 
 LO[i1]=0.0; HI[i1]=18.0; i1++;
 histNames[i1] = "ptpib"; histExp[i1] = "0.001*BachPi_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=15; 
 LO[i1]=0.0; HI[i1]=1.5; i1++;
 histNames[i1] = "pibIPCHI2"; histExp[i1] = "log(BachPi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=20; 
 LO[i1]=1.0; HI[i1]=10.0; i1++;
 nHistsBachPi = i1 - histBachPif ;


 HB = "#Xi^{#font[122]{-}}";
 histXif = i1; 
 histNames[i1] = "pXi"; histExp[i1] = "0.001*Xi_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=10; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 histNames[i1] = "ptXi"; histExp[i1] = "0.001*Xi_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=16; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 histNames[i1] = "etaXi"; histExp[i1] = "-log(tan(0.5*asin(Xi_PT/Xi_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=14; 
 LO[i1]=1.5; HI[i1]=5.0; i1++;
 nHistsXi = i1 - histXif ;

 histXiVtxf = i1; 
 histNames[i1] = "XiTau"; histExp[i1] = "Xi_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=12; 
 LO[i1]=0.0; HI[i1]=0.12; i1++;
 histNames[i1] = "XiChisq"; histExp[i1] = "Xi_ENDVERTEX_CHI2"; hTit[i1] = HB+" #chi^{2}_{vtx}"; NB[i1]=25; 
 LO[i1]=0.0; HI[i1]=25.0; i1++;
 histNames[i1] = "XiDira"; histExp[i1] = "1000*acos(Xi_DIRA_ORIVX)"; hTit[i1] = HB+" #theta_{DIRA, ORIVX}"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 nHistsXiVtx = i1 - histXiVtxf ;

 nHists = i1;

 cout << "Number of histograms " << nHists << endl;
 


 nHists = i1;
 histSys = nHists;
 
 
 // Book the histograms
 for(int i=0; i<nHists; i++){
   hda[i] = new TH1F(histNames[i]+"_da",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i] = new TH1F(histNames[i]+"_mc",hTit[i],NB[i],LO[i],HI[i]);
 }

}

void fillLbHistos(TTree* tree_data, TTree *tree_sim){
  
  TString hist;
  double vint;
  //simCut = MCLambdabDD_Cut;
  TCanvas *ctest = new TCanvas("ctest","",600,600);
  
  for(int i=0; i<nHists; i++){
    hist = histNames[i]+"_da";
    cout << "Filling histogram: " << i << " " << hTit[i] << " " << histExp[i] << " 0 " << " " << signalWeight << " " << tree_data << endl;
    tree_data->Draw(histExp[i]+">>"+hist,signalWeight);
    if(i==1){
      hda[i]->Draw();
      ctest->Update();      
    }
    cout << "Integral: " << hda[1]->Integral() << endl;
    
    addGraphics(hda[i],hTit[i],"",1);
    hist = histNames[i]+"_mc";
    cout << "Filling histogram: " << i << " " << hTit[i] << " 1 " << endl;
    tree_sim->Draw(histExp[i]+">>"+hist,simCut);    
    addGraphics(hmc[i],hTit[i],"",2);
    hmc[i]->Sumw2();
    if(i<histSys) hmc[i]->Scale(hda[i]->Integral()/hmc[i]->Integral());
    cout << "Filling histogram: " << i << " " << hTit[i] << " 2 " << endl;
  }

  legend1 = new TLegend(0.70,0.70,0.88,0.88);
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->SetTextSize(0.055);
  legend1->AddEntry(hda[0],"Data","LEP");
  legend1->AddEntry(hmc[0],"Sim","L");
  legend1->Draw();

  
}

TChain* getLbDDFiles(){

  bool ex;
  long siz;
  chainl = new TChain("","MyTuple");
  // reco'd
  if(year=="Run1"){
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  }else if(year=="Run2"){
    if(year=="Run2"){
      TString file, ffile;
      // Lb --> Jpsi L0
      int nFile = 20;
      // Get reconstructed MC -- MagDown
      for(int i=0; i<nFile; i++){
        file = Form("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2016_MagDown/1068_%d/jpsilambda.root",i);
        ffile = file + "/Lb2JpsiLTree/MyTuple";
        ex = exists(file);
        if(!ex) continue;
        siz = GetFileSize(file);
        cout << "file: " << file << ", size = " << siz << endl;
        chainl->Add(ffile);
      }
      // -- MagUp
      for(int i=0; i<nFile; i++){
        file = Form("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2016_MagUp/1069_%d/jpsilambda.root",i);
        ffile = file + "/Lb2JpsiLTree/MyTuple";
        ex = exists(file);
        if(!ex) continue;
        siz = GetFileSize(file);
        cout << "file: " << file << ", size = " << siz << endl;
        chainl->Add(ffile);
      }
    }
  }


  TString sLb_Trig = Lb_Trig.GetTitle();
  tag = "DD";

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
  

  TString filename  = "Lb2JpsiL0_"+year+"_"+tag+".root";
  TFile *f = new TFile(filename);
  TreeWithWeight = (TTree*)f->Get("TreeWithWeight");

  cout << filename << " " << ",  TreeWithWeight = " << TreeWithWeight  << endl;
  
  return chainl;
}

TChain* getLbLLFiles(){
  bool ex;
  long siz;

  chainl = new TChain("","MyTuple");
  // reco'd
  if(year=="Run1"){
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
    chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  }else if(year=="Run2"){
    if(year=="Run2"){
      TString file, ffile;
      // Lb --> Jpsi L0
      int nFile = 20;
      // Get reconstructed MC -- MagDown
      for(int i=0; i<nFile; i++){
        file = Form("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2016_MagDown/1068_%d/jpsilambda.root",i);
        ffile = file + "/Lb2JpsiLTree/MyTuple";
        ex = exists(file);
        if(!ex) continue;
        siz = GetFileSize(file);
        cout << "file: " << file << ", size = " << siz << endl;
        chainl->Add(ffile);
      }
      // -- MagUp
      for(int i=0; i<nFile; i++){
        file = Form("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2016_MagUp/1069_%d/jpsilambda.root",i);
        ffile = file + "/Lb2JpsiLTree/MyTuple";
        ex = exists(file);
        if(!ex) continue;
        siz = GetFileSize(file);
        cout << "file: " << file << ", size = " << siz << endl;
        chainl->Add(ffile);
      }
    }
  }
    


  TString sLb_Trig = Lb_Trig.GetTitle();
  tag = "LL";

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
  

  TString filename  = "Lb2JpsiL0_"+year+"_"+tag+".root";
  TFile *f = new TFile(filename);
  TreeWithWeight = (TTree*)f->Get("TreeWithWeight");

  cout << filename << " " << ",  TreeWithWeight = " << TreeWithWeight  << endl;

  
  return chainl;

}

TChain* getXibFiles(){

  chainl = new TChain("","MyTuple");
  // reco'd
  TString mcdir = "/data1/avenkate/JpsiXi/mc/Pythia8/old";
  mcdir = "/data1/avenkate/JpsiXi/mc/Pythia8/";
  TString ids[4] = {"1085","1086","1087","1088"};
  
  // reco'd
  chainl->Add(mcdir+"2011_MagDown/"+ids[0]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add(mcdir+"2011_MagUp/"+ids[1]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add(mcdir+"2012_MagDown/"+ids[2]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add(mcdir+"2012_MagUp/"+ids[3]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  //chainl->Add(mcdir+"2011_MagDown/"+ids[0]+"_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  //chainl->Add(mcdir+"2011_MagUp/"+ids[1]+"_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  //chainl->Add(mcdir+"2012_MagDown/"+ids[2]+"_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  //chainl->Add(mcdir+"2012_MagUp/"+ids[3]+"_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  tag = "DD";
  TString sLb_Trig = Xib_Trig.GetTitle();
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


  TString fileName = "Xib2JpsiXi_"+year+"_"+tag;
  fileName = fileName + "_NewStrip";;
  if(applyXiPiTrackingAcceptanceCut) fileName = fileName+"_withTrackingAcc";
  fileName = fileName +".root";

  /*
   TString fileName = "Xib2JpsiXi_L0JPsiTOS";;
   if(useNewStrip) fileName = "Xib2JpsiXi_L0JPsiTOS_NewStrip";;
   if(applyXiPiTrackingAcceptanceCut) fileName = fileName+"_withTrackingAcc";
   fileName = fileName +".root";
  */

   TFile *f = new TFile(fileName);
   TreeWithWeight = (TTree*)f->Get("TreeWithWeight");
   cout << "================================="<<endl;
   cout << TreeWithWeight << endl;
   cout << "================================="<<endl;
  
  return chainl;
  
}


void getWeightHistos(){
  TFile *fw = new TFile("fsig_Lb_new.save.root");
  //TH2F *whistLb = (TH2F*)fw->Get("hratio");
  whistLb = (TH2F*)fw->Get("hratioLb")->Clone("whistLb");
  TFile *fw2 = new TFile("fsig_Lb_new2.save.root");
  whistL = (TH2F*)fw2->Get("hratioL")->Clone("whistL");
  //if(applyPrPiKinWeight){
  TFile *fw3 = new TFile("fsig_Lb_new3.save.root");
  whistPi = (TH2F*)fw3->Get("hratioPi")->Clone("whistPi");
  whistPrPi = (TH1F*)fw3->Get("hratioPrPi")->Clone("whistPrPi");
  if(year=="Run1"){
    TFile *ftr = new TFile("ratio2012S20.root");
    trackEff = (TH2D*)ftr->Get("Ratio");
  }else if(year=="Run2"){
    TFile *ftr = new TFile("Ratio_Long_P-ETA_2016_25ns.root");
    trackEff = (TH2D*)ftr->Get("Ratio");
  }
  
  

  //}
  /*
    TCanvas *cw = new TCanvas("cw","weights",600,600);
    cw->Divide(2,2);
    cw->cd(1);
    whistLb->Draw();
    cw->cd(2);
    whistL->Draw();
    cw->cd(3);
    whistPi->Draw();
    cw->cd(4);
    whistPrPi->Draw();
  */


    //fw->Close();
    //fw2->Close();
    //fw3->Close();


}

TTree* getNewTree(TChain *chain){

  cout << "=============================================="<<endl;  
  cout << whistPi << " " << whistLb << " " << whistL << endl;
  cout << "=============================================="<<endl;

  int ib;
  double eta, w, wPi, etaPi, pi_PX, pi_PY, pi_PZ, p_PX, p_PY, p_PZ, ptPi;
  double w1, w2, w3, pi_PT, p_PT;
  double L_P, L_PT, Lb_P, Lb_PT, L_TAU;
  chain->SetBranchAddress("Lb_P", &Lb_P);
  chain->SetBranchAddress("Lb_PT", &Lb_PT);
  chain->SetBranchAddress("L_P", &L_P);
  chain->SetBranchAddress("L_PT", &L_PT);
  chain->SetBranchAddress("pi_PT", &pi_PT);
  chain->SetBranchAddress("pi_PZ", &pi_PZ);
  chain->SetBranchAddress("p_PT", &p_PT);
  chain->SetBranchAddress("p_PZ", &p_PZ);
  Long64_t nentries = chain->GetEntriesFast();
  cout << nentries << endl;
  Long64_t ientry ;
  //nentries = 100;
  TFile *ff = new TFile("/data1/sblusk/junk2.root","RECREATE");
  prPi = new TProfile("prPi","",20,0,2.0,0.0,4.0);
  double swt, swtPi, swtLb, swtL;
  newtree = (TTree*) chain->CloneTree(0);
  newtree->Branch("swtLb", &swtLb);
  newtree->Branch("swtL", &swtL);
  newtree->Branch("swtPi", &swtPi);
  for (Long64_t ii=0;ii<nentries; ii++) {
    ientry = chain->LoadTree(ii);
    if (ientry < 0) break;
    chain->GetEntry(ii);
    if(ii%10000 == 0) cout << "At entry = " << ii << endl;
    swtPi = 1.0;
    //    if(applyPrPiKinWeight){
      if(mcOpt==4 ){
        swtPi = get2DKinWt(whistPi,pi_PT,pi_PZ,0);
        prPi->Fill(0.001*pi_PT,swtPi);
        if(ii<200){
          eta = -log(tan(0.5*atan(pi_PT/pi_PZ)));
          cout << ii << " " << 0.001*pi_PT << " " << eta << " " << swtPi << endl;
        }
      }else if(mcOpt==5){
        pPr = sqrt(p_PT*p_PT+p_PZ*p_PZ);
        pPi = sqrt(pi_PT*pi_PT+pi_PZ*pi_PZ);
        asym = (pPr-pPi) / (pPr+pPi);      
        swtPi = get1DKinWt(whistPrPi,asym);
      }
      // }
    
    swtLb = get2DKinWt(whistLb,Lb_PT,Lb_P,1);
    swtL = get2DKinWt(whistL,L_PT,L_P,1);

    newtree->Fill(); 
  }
  newtree->AutoSave();
  TCut swLb = "(swtLb)";
  TCut swL = "(swtL)";
  TCut swPi = "(swtPi)";
  if(mcOpt>=2 && applyLbKinWeight) simCut = simCut*swLb;
  if(mcOpt>=3 && applyLKinWeight) simCut = simCut*swL;
  if(mcOpt>=4 && applyPrPiKinWeight) simCut = simCut*swPi;  
  if(applyJpsiLWeight) simCut = simCut*JpsiLWeight;
  
  return newtree;
  
}

TTree* getNewTreeXib(TChain* chain){

  int ib;
  double eta, w, wPi, etaPi, pi_PT, pi_PZ, p_PT, p_PZ, ptPi, pPi, pPr, asym;
  double L_P, L_PT, Xib_P, Xib_PT, Xi_P, Xi_PT, L_TAU, BachPi_P;
  chain->SetBranchAddress("Xib_P", &Xib_P);
  //chain->SetBranchAddress("BachPi_P", &BachPi_P);
  chain->SetBranchAddress("Xi_PT", &Xi_PT);
  chain->SetBranchAddress("Xib_PT", &Xib_PT);
  chain->SetBranchAddress("L_P", &L_P);
  chain->SetBranchAddress("L_PT", &L_PT);
  //chain->SetBranchAddress("Xi_P", &Xi_P);
  //chain->SetBranchAddress("Xi_PT", &Xi_PT);
  chain->SetBranchAddress("pi_PT", &pi_PT);
  chain->SetBranchAddress("pi_PZ", &pi_PZ);
  chain->SetBranchAddress("p_PT", &p_PT);
  chain->SetBranchAddress("p_PZ", &p_PZ);
  Long64_t nentries = chain->GetEntriesFast();
  cout << nentries << endl;
  Long64_t ientry ;
  TFile *ff = new TFile("/data1/sblusk/junk2.root","RECREATE");
  double swt, swtPi, swtLb, swtL;
  newtree = (TTree*) chain->CloneTree(0);
  newtree->Branch("swtLb", &swtLb);
  newtree->Branch("swtL", &swtL);
  newtree->Branch("swtPi", &swtPi);
  for (Long64_t i=0;i<nentries; i++) {
    ientry = chain->LoadTree(i);
    if (ientry < 0) break;
    chain->GetEntry(i);
    if(i%10000 == 0) cout << "At entry = " << i << endl;
    swtPi = 1.0;

    if(applyPrPiKinWeight){
      if(mcOpt==4){
        swtPi = get2DKinWt(whistPi,pi_PT,pi_PZ,0);
      }else if(mcOpt==5){
        pPr = sqrt(p_PT*p_PT+p_PZ*p_PZ);
        pPi = sqrt(pi_PT*pi_PT+pi_PZ*pi_PZ);
        asym = (pPr-pPi) / (pPr+pPi);
        swtPi = get1DKinWt(whistPrPi,asym);
      }
    }
    
    swtLb = get2DKinWt(whistLb,Xib_PT,Xib_P,1);
    swtL = get2DKinWt(whistL,L_PT,L_P,1);

    newtree->Fill(); 
    //if(i<100) cout << ib << " " << Xib_PT << " " << etaXib << " " << w << endl;
  }
  newtree->AutoSave();
  TCut swLb = "(swtLb)";
  TCut swL = "(swtL)";
  TCut swPi = "(swtPi)";
  simCut = (simCut && XiL_Cut);
  // Weights from Lb samples
  if(mcOpt>=2 && applyLbKinWeight) simCut = simCut*swLb;
  if(mcOpt>=3 && applyLKinWeight) simCut = simCut*swL;
  if(mcOpt>=4 && applyPrPiKinWeight) simCut = simCut*swPi;
  // Extra weight for Xib- decay
  if(applyEtaWeight) simCut = simCut*XibWt;
  if(applyJpsiXiWeight) simCut = simCut*JpsiXiWeight ;
  if(applyXiWeight) simCut = simCut*XiWeight ;

  return newtree;


}


void makeWeightFile(){

  if(mcOpt <= 3){    
    TFile *fout;
    if(mcOpt==1) fout = new TFile("fsig_Lb_new.root","RECREATE");
    if(mcOpt==2) fout = new TFile("fsig_Lb_new2.root","RECREATE");
    if(mcOpt==3) fout = new TFile("fsig_Lb_new3.root","RECREATE");
    TCanvas *c2 = new TCanvas("c2","",400,400);

    TH2F *h201 = new TH2F("h201","eta vs pT",25,0,25,10,1.5,6.5);
    TH2F *h301 = new TH2F("h301","eta vs pT",25,0,25,10,1.5,6.5);
    TreeWithWeight->Draw("-log(tan(0.5*asin(Lb_PT/Lb_P))):0.001*Lb_PT>>h201","nsig1_sw");
    newtree->Draw("-log(tan(0.5*asin(Lb_PT/Lb_P))):0.001*Lb_PT>>h301",simCut);
    TH2F *hratioLb = (TH2F*)h201->Clone("hratioLb");
    h301->Scale(h201->Integral()/h301->Integral());    
    hratioLb->Divide(h301);    
    h201->Write();
    h301->Write();
    hratioLb->Write();

    TH2F *h202 = new TH2F("h202","eta vs pT",25,0,10,7,1.5,5.0);
    TH2F *h302 = new TH2F("h302","eta vs pT",25,0,10,7,1.5,5.0);
    TreeWithWeight->Draw("-log(tan(0.5*asin(L_PT/L_P))):0.001*L_PT>>h202","nsig1_sw");
    newtree->Draw("-log(tan(0.5*asin(L_PT/L_P))):0.001*L_PT>>h302",simCut);
    TH2F *hratioL = (TH2F*)h202->Clone("hratioL");
    h302->Scale(h202->Integral()/h302->Integral());    
    hratioL->Divide(h302);    
    h202->Write();
    h302->Write();
    hratioL->Write();

    TH2F *h203 = new TH2F("h203","eta vs pT",20,0,2,16,1.5,5.5);
    TH2F *h303 = new TH2F("h303","eta vs pT",20,0,2,16,1.5,5.5);
    TreeWithWeight->Draw("-log(tan(0.5*atan(sqrt(pi_PX**2+pi_PY**2)/pi_PZ))):0.001*sqrt(pi_PX**2+pi_PY**2)>>h203","nsig1_sw");
    newtree->Draw("-log(tan(0.5*atan(sqrt(pi_PX**2+pi_PY**2)/pi_PZ))):0.001*sqrt(pi_PX**2+pi_PY**2)>>h303",simCut);
    TH2F *hratioPi = (TH2F*)h203->Clone("hratioPi");
    h303->Scale(h203->Integral()/h303->Integral());    
    hratioPi->Divide(h303);    
    h203->Write();
    h303->Write();
    hratioPi->Write();

    TH1F *h204 = new TH1F("h204","Momentum asymmetry",100,0.5,1);
    TH1F *h304 = new TH1F("h304","Momentum asymmetry",100,0.5,1);
    TreeWithWeight->Draw("(sqrt(p_PX**2+p_PY**2+p_PZ**2)-sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))/(sqrt(p_PX**2+p_PY**2+p_PZ**2)+sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))>>h204","nsig1_sw");
    newtree->Draw("(sqrt(p_PX**2+p_PY**2+p_PZ**2)-sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))/(sqrt(p_PX**2+p_PY**2+p_PZ**2)+sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))>>h304",simCut);
    TH1F *hratioPrPi = (TH1F*)h204->Clone("hratioPrPi");
    h304->Scale(h204->Integral()/h304->Integral());    
    hratioPrPi->Divide(h304);    
    h204->Write();
    h304->Write();
    hratioPrPi->Write();

    fout->Write();
    fout->Close();
  }


}

void makeLbBasePlots(){
  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);


  cout << "============================="<<endl;
  cout << "I am here - 0" << endl;
  cout << "============================="<<endl;

  // Plot results
  TCanvas *c0 = new TCanvas("c0","Lambda0",1200,800);
  cout << "I made it past Tcanvas c0" << endl;
  
  c0->Divide(3,2);
  c0->cd(1);
  int j,k;
  for(int i=0;i<nHistsL;i++){
    c0->cd(i+1);
    c0->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histLf + i;
    cout << i << " " << j << endl;    
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }


  cout << "============================="<<endl;
  cout << "I am here - 1" << endl;
  cout << "============================="<<endl;

  TCanvas *ca = new TCanvas("ca","b-hadron",1200,400);
  ca->Divide(3,1);
  for(int i=0;i<nHistsLb;i++){
    ca->cd(i+1);
    ca->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histLbf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }

  
  cout << "============================="<<endl;
  cout << "I am here - 2" << endl;
  cout << "============================="<<endl;


  TCanvas *cb = new TCanvas("cb","Proton/pion",600,800);
  cb->Divide(2,3);
  for(int i=0;i<nHistsp+nHistspi;i++){
    cb->cd(i+1);
    cb->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histpf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  // Form asymmetry
    if(i==nHistsp+nHistspi-1){      
      cb->cd(i+2);
      cb->cd(i+2)->SetBottomMargin(bottomMarginOffset);
      TH1F *hasym = (TH1F*)hda[j]->Clone("hasym");
      hasym->Divide(hmc[j]);
      hasym->GetYaxis()->SetRangeUser(0,2);
      hasym->GetYaxis()->SetTitle("Momentum asymmetry");
      hasym->Draw("e");   
    }                               
  }


  cout << "============================="<<endl;
  cout << "I am here - 3" << endl;
  cout << "============================="<<endl;

  double v1[2];
  TF1 *func = new TF1("func","[0]+[1]*x",0.2,1);
  TCanvas *cc = new TCanvas("cc","Ratio (Jpsi, L) / b-hadron",800,800);
  cc->Divide(2,2);
  for(int i=0;i<nHistsR;i++){
    cc->cd(i+1);
    cc->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histRf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==2) legend1->Draw();

    k = nRatio+i;
    cc->cd(i+3)->SetBottomMargin(bottomMarginOffset);;
    hr[k] = (TH1F*)hda[j]->Clone(Form("hr_%d",k));
    hr[k]->Divide(hmc[j]);
    hr[k]->SetMaximum(2);
    if(i==0) { 
      func->SetParameters(0.2,1.4);
      func->SetRange(0.24,1.0);
    }else if(i==1) {
      func->SetRange(0.00,0.66);
      func->SetParameters(1.3,-0.8);
    }  
    hr[k]->Fit("func","RM");
    v1[0] = func->GetParameter(0);
    v1[1] = func->GetParameter(1);
    if(parentPart=="Lb"){
      if(i==0) myLatex->DrawLatex(0.18,0.8,Form("%4.2f+%5.3f (p_{J/#psi}/p_{#Lambda_{b}})",v1[0],v1[1]));
      if(i==1) myLatex->DrawLatex(0.18,0.8,Form("%4.2f+%5.3f (p_{#Lambda}/p_{#Lambda_{b}})",v1[0],v1[1]));
    }else{
      if(i==0) myLatex->DrawLatex(0.18,0.8,Form("%4.2f+%5.3f (p_{J/#psi}/p_{#Xi_{b}^{#font[122]{-}}})",v1[0],v1[1]));
      if(i==1) myLatex->DrawLatex(0.18,0.8,Form("%4.2f+%5.3f (p_{#Lambda}/p_{#Xi_{b}^{#font[122]{-}}})",v1[0],v1[1]));
    }
    
  }
  nRatio = nRatio + 2;


  cout << "============================="<<endl;
  cout << "I am here - 4" << endl;
  cout << "============================="<<endl;

  TCanvas *cd = new TCanvas("cd","Jpsi",800,800);
  //cd->Divide(2,2);
  for(int i=0;i<nHistsJpsi;i++){
    //cd->cd(i+1);
    cd->SetBottomMargin(bottomMarginOffset);
    j = histJpsif + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }  


  TCanvas *cz = new TCanvas("cz","PV",800,400);
  cz->Divide(2,1);
  TF1 *gfunc = new TF1("gfunc","[0]*exp(-0.5*(x-[1])**2/[2]**2)",-200,200);
  gfunc->SetParNames("Con","Mean","Sigma");
  for(int i=0;i<nHistsPV;i++){
    cz->cd(i+1);
    cz->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histPV + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
    cz->cd(i+2)->SetBottomMargin(bottomMarginOffset);;
    k = nRatio + i;
    hr[k] = (TH1F*)hda[j]->Clone(Form("hr_%d",k));
    hr[k]->SetStats(kTRUE);
    hr[k]->Divide(hmc[j]);
    hr[k]->SetMaximum(2);
    hr[k]->SetMinimum(0);
    hr[k]->GetYaxis()->SetTitle("ratio");
    hr[k]->Fit("gaus","RM");
  }  
  nRatio++;

  cout << "============================="<<endl;
  cout << "I am here - 5" << endl;
  cout << "============================="<<endl;

  TString name = "Lambdab_from_"+parentPart+tag;
  if(parentPart=="Xib") name = "Xib_from_"+parentPart+tag;
  if(makePDF){
    makePlots(c0,"Lambda_from_"+parentPart+tag);
    makePlots(ca,name);
    makePlots(cb,"ProtonPionMomenta_from_"+parentPart+tag);
    makePlots(cc,"MomentumRatio_from_"+parentPart+tag);
    makePlots(cd,"Jpsi_from_"+parentPart+tag);
    makePlots(cz,"PV_from_"+parentPart+tag);
  }


  cout << "============================="<<endl;
  cout << "I am here - 6" << endl;
  cout << "============================="<<endl;

}

void makeLbPlotsLL(){
  int j,k;

  TCanvas *ce = new TCanvas("ce","",800,600);
  ce->Divide(3,2);
  for(int i=0;i<nHistsEx;i++){
    ce->cd(i+1);
    ce->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histExf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }
  
  if(makePDF){
    makePlots(ce,"ExtraLambda_from_"+parentPart+tag);
  } 
}

void makeXibPlots(){

  cout << histBachPif << " " << nHistsBachPi << endl;
  int j,k;
  cout << "I am here = 000 " << endl;
  TCanvas *cf = new TCanvas("cf","BachPi",900,300);
  cout << "I am here = 00 " << endl;
  cf->Divide(3,1);
  cout << "I am here = 0 " << endl;
  for(int i=0;i<nHistsBachPi;i++){
    cf->cd(i+1);
    cout << "I am here = 1 " << endl;
    cf->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    cout << "I am here = 2 " << endl;
    j = histBachPif + i;
    //cout << " j =  " << j << " " << hTit[j] << endl;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    cout << "I am here = 3 " << endl;
    if(i==0) legend1->Draw();
  }

  TH1F* hrXi[3];
  TCanvas *cg = new TCanvas("cg","Xi",900,600);

  TF1 *func2 = new TF1("func2","[0]+[1]*x+[2]*x*x",1.0,7.5);
  func2->SetParameters(1.6,-0.35,0.03);

  cg->Divide(3,2);
  for(int i=0;i<nHistsXi;i++){
    cg->cd(i+1);
    cg->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histXif + i;
    //cout << " j =  " << j << " " << hTit[j] << endl;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
    cg->cd(i+1+3);
    cg->cd(i+1+3)->SetBottomMargin(bottomMarginOffset);    
    hrXi[i] = (TH1F*)hda[j]->Clone(Form("hrXi_%d",i));
    hrXi[i]->Divide(hmc[j]);
    hrXi[i]->GetYaxis()->SetRangeUser(0,4);
    if(i==1){
      hrXi[i]->Fit("func2","RM");
    }
    
    hrXi[i]->Draw();    
    
    
  }


  double maxim1, maxim2, v1, v2, dv;
  TCanvas *cj = new TCanvas("cj","Xi Vtx",1200,400);
  cj->Divide(3,1);
  for(int i=0;i<nHistsXiVtx;i++){
    cj->cd(i+1);
    cj->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histXiVtxf + i;
    //cout << " j =  " << j << " " << hTit[j] << endl;
    maxim1 = hmc[j]->GetMaximum();
    maxim2 = hda[j]->GetMaximum();
    if(maxim1 > maxim2){
      hmc[j]->GetYaxis()->SetRangeUser(-5,1.2*maxim1);
      hmc[j]->Draw("hist");  
      hda[j]->Draw("e,same");
    }else{
      hda[j]->GetYaxis()->SetRangeUser(-5,1.4*maxim1);
      hda[j]->Draw("e");
      hmc[j]->Draw("hist,same");  
    }    
    v1 = hmc[j]->Integral(13,25)/hmc[j]->Integral();
    v2 = hda[j]->Integral(13,25)/hda[j]->Integral();
    dv = v2 - v1;
    cout << "Difference between data and mc = " << v1 << " " << v2 << " " << dv << endl;
    
    if(i==0) legend1->Draw();
  }

  TH1F *hrp[4];
  TCanvas *ch = new TCanvas("ch","Proton/pion ratios",800,800);
  ch->Divide(2,2);
  for(int i=0;i<nHistsp+nHistspi;i++){
    ch->cd(i+1);
    ch->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histpf + i;

    hrp[i] = (TH1F*)hda[j]->Clone(Form("hrp_%d",i));
    hrp[i]->Divide(hmc[j]);
    hrp[i]->GetYaxis()->SetRangeUser(0,4);
    hrp[i]->Draw();    

  }


  
  if(makePDF){
    makePlots(cf,"BachPi_from_"+parentPart+tag);
    makePlots(cg,"Xi_from_"+parentPart+tag);
    makePlots(cj,"XiVtx_from_"+parentPart+tag);
    makePlots(ch,"ProtonPionRatios_from_"+parentPart+tag);
  } 
}


void getSyst(){

  double frd[10];
  double frm[10];
  int j;
  
  for(int i=0;i<nHistsSys;i++){
    j = histSys + i;
    frd[i] = hda[j]->Integral(9,25)/hda[j]->GetEntries();
    frm[i] = hmc[j]->Integral(9,25)/hmc[j]->GetEntries();
    cout << hTit[j] << endl;
    cout << " --> Data fraction: " << frd[i] << endl;
    cout << " --> MC fraction: " << frm[i] << endl;
  }
  
}




void getLbMCFilesForPlot(){

  TString file, ffile;

  chainl = new TChain("","MyTuple");
  chainlt = new TChain("","MCDecayTree");
  chainx = new TChain("","MyTuple");
  chainxt = new TChain("","MCDecayTree");

  TString mcdir = "/data1/avenkate/JpsiXi/mc/Pythia8/old";
  mcdir = "/data1/avenkate/JpsiXi/mc/Pythia8/";
  TString ids[4] = {"1085","1086","1087","1088"};

  TString Run1Dir[4] = {"2011_MagDown","2011_MagUp","2012_MagDown","2012_MagUp" };
  TString Run2Dir[2] = {"2016_MagDown","2016_MagUp"};
  
  
  bool ex;
  if(year=="Run1"){
    int idLb = 920;
    int numFiles = 1;
    if(mcFit) numFiles = 4;
    TString mcdirLb = "/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/";
    for(int j=0;j<numFiles;j++){
      for(int k=0; k<4;k++){
        file = Form("%d_%d",idLb+k,j);
        file = mcdirLb+Run1Dir[k]+"/"+file+"/jpsilambda.root";
        ex = exists(file);
        if(!ex) continue;
        ffile = file + "/MCTuple/MCDecayTree";
        chainlt->Add(ffile);
        ffile = file + "/Lb2JpsiLTree/MyTuple";
        chainl->Add(ffile);
      }
    }   

    idLb = 1085;
    mcdirLb = "/data1/avenkate/JpsiXi/mc/Pythia8/";
    for(int j=0;j<numFiles;j++){
      for(int k=0; k<4;k++){
        file = Form("%d_%d",idLb+k,j);
        file = mcdirLb+Run1Dir[k]+"/"+file+"/jpsixi.root";
        ex = exists(file);
        if(!ex) continue;
        ffile = file + "/MCTuple/MCDecayTree";
        chainxt->Add(ffile);
        ffile = file + "/Xib2JpsiXiTree/MyTuple";
        chainx->Add(ffile);
      }
    }   

    /*

    // MC Truth, all events
    chainxt->Add(mcdir+"2011_MagDown/"+ids[0]+"_0/jpsixi.root/MCTuple/MCDecayTree");
    chainxt->Add(mcdir+"2011_MagUp/"+ids[1]+"_0/jpsixi.root/MCTuple/MCDecayTree");
    chainxt->Add(mcdir+"2012_MagDown/"+ids[2]+"_0/jpsixi.root/MCTuple/MCDecayTree");
    chainxt->Add(mcdir+"2012_MagUp/"+ids[3]+"_0/jpsixi.root/MCTuple/MCDecayTree");
    // reco'd
    chainx->Add(mcdir+"2011_MagDown/"+ids[0]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
    chainx->Add(mcdir+"2011_MagUp/"+ids[1]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
    chainx->Add(mcdir+"2012_MagDown/"+ids[2]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
    chainx->Add(mcdir+"2012_MagUp/"+ids[3]+"_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
    */

  }

  if(year=="Run2"){

    // Lb --> Jpsi L0
    int nFile = 30;
    // Get reconstructed MC -- MagDown
    for(int i=0; i<nFile; i++){
      file = Form("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2016_MagDown/1068_%d/jpsilambda.root",i);
      ffile = file + "/Lb2JpsiLTree/MyTuple";
      ex = exists(file);
      if(!ex) continue;
      chainl->Add(ffile);
      ffile = file + "/MCTuple/MCDecayTree";
      chainlt->Add(ffile);
    }
    // -- MagUp
    for(int i=0; i<nFile; i++){
      file = Form("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2016_MagUp/1069_%d/jpsilambda.root",i);
      ffile = file + "/Lb2JpsiLTree/MyTuple";
      ex = exists(file);
      if(!ex) continue;
      chainl->Add(ffile);
      ffile = file + "/MCTuple/MCDecayTree";
      chainlt->Add(ffile);
    }

    // Xib --> Jpsi Xi
    nFile = 30;
    // Get reconstructed MC -- MagDown
    for(int i=0; i<nFile; i++){
      file = Form("/data1/avenkate/JpsiXi/mc/Pythia8/2016_MagDown/1091_%d/jpsixi.root",i);
      ffile = file + "/Xib2JpsiXiTree/MyTuple";
      ex = exists(file);
      if(!ex) continue;
      chainx->Add(ffile);
      ffile = file + "/MCTuple/MCDecayTree";
      chainxt->Add(ffile);
    }
    // -- MagUp
    for(int i=0; i<nFile; i++){
      file = Form("/data1/avenkate/JpsiXi/mc/Pythia8/2016_MagUp/1092_%d/jpsixi.root",i);
      ffile = file + "/Xib2JpsiXiTree/MyTuple";
      ex = exists(file);
      if(!ex) continue;
      chainx->Add(ffile);
      ffile = file + "/MCTuple/MCDecayTree";
      chainxt->Add(ffile);
    }
  }
  
    
  return;
  

  }

void setLbBranchAddresses(){
  newtree->SetBranchAddress("Jpsi_P", &Jpsi_P);
      newtree->SetBranchAddress("Jpsi_TRUEPT", &Jpsi_TRUEPT);
      newtree->SetBranchAddress("Jpsi_TRUEP_Z", &Jpsi_TRUEP_Z);
      newtree->SetBranchAddress("Lb_TAU", &Lb_TAU);
      newtree->SetBranchAddress("Lb_TRUETAU", &Lb_TRUETAU);
      newtree->SetBranchAddress("Lb_P", &Lb_P);
      newtree->SetBranchAddress("Lb_PT", &Lb_PT);
      newtree->SetBranchAddress("L_P", &L_P);
      newtree->SetBranchAddress("L_PT", &L_PT);
      newtree->SetBranchAddress("pi_PT", &pi_PT);
      newtree->SetBranchAddress("pi_PZ", &pi_PZ);
      newtree->SetBranchAddress("p_PT", &p_PT);
      newtree->SetBranchAddress("p_PZ", &p_PZ);

      newtree->SetBranchAddress("Jpsi_TRUEP_Z", &Jpsi_TRUEP_Z);
      newtree->SetBranchAddress("Jpsi_TRUEPT", &Jpsi_TRUEPT);
      newtree->SetBranchAddress("Lb_TRUEP_Z", &Lb_TRUEP_Z);
      newtree->SetBranchAddress("Lb_TRUEPT", &Lb_TRUEPT);
      newtree->SetBranchAddress("L_TRUEP_Z", &L_TRUEP_Z);
      newtree->SetBranchAddress("L_TRUEPT", &L_TRUEPT);
      newtree->SetBranchAddress("pi_TRUEPT", &pi_TRUEPT);
      newtree->SetBranchAddress("pi_TRUEP_Z", &pi_TRUEP_Z);
      newtree->SetBranchAddress("p_TRUEPT", &p_TRUEPT);
      newtree->SetBranchAddress("p_TRUEP_Z", &p_TRUEP_Z);
      //newtree->SetBranchAddress("eventNumber", &eventNumber);
      

      chainlt->SetBranchAddress("J_psi_1S_TRUEP_Z", &J_psi_1S_TRUEP_Z);
      chainlt->SetBranchAddress("J_psi_1S_TRUEPT", &J_psi_1S_TRUEPT);
      chainlt->SetBranchAddress("Lambda_b0_TRUEP_Z", &Lambda_b0_TRUEP_Z);
      chainlt->SetBranchAddress("Lambda_b0_TRUEPT", &Lambda_b0_TRUEPT);
      chainlt->SetBranchAddress("Lambda0_TRUEP_Z", &Lambda0_TRUEP_Z);
      chainlt->SetBranchAddress("Lambda0_TRUEPT", &Lambda0_TRUEPT);
      chainlt->SetBranchAddress("Lambda_b0_TRUETAU", &Lambda_b0_TRUETAU);
      chainlt->SetBranchAddress("piminus_TRUEPT", &piminus_TRUEPT);
      chainlt->SetBranchAddress("piminus_TRUEP_Z", &piminus_TRUEP_Z);
      chainlt->SetBranchAddress("pplus_TRUEPT", &pplus_TRUEPT);
      chainlt->SetBranchAddress("pplus_TRUEP_Z", &pplus_TRUEP_Z);
      //chainlt->SetBranchAddress("eventNumber", &eventNumber);

      return;
      


}

void setXibBranchAddresses(){

      newtree->SetBranchAddress("BachPi_TRACK_Type", &BachPi_TRACK_Type);
      newtree->SetBranchAddress("BachPi_P", &BachPi_P);
      newtree->SetBranchAddress("BachPi_PT", &BachPi_PT);
      newtree->SetBranchAddress("Jpsi_P", &Jpsi_P);
      newtree->SetBranchAddress("Xib_P", &Xib_P);
      newtree->SetBranchAddress("Xib_PT", &Xib_PT);
      newtree->SetBranchAddress("Xib_TAU", &Xib_TAU);
      newtree->SetBranchAddress("Xib_TRUETAU", &Xib_TRUETAU);
      newtree->SetBranchAddress("L_P", &L_P);
      newtree->SetBranchAddress("L_PT", &L_PT);
      newtree->SetBranchAddress("pi_PT", &pi_PT);
      newtree->SetBranchAddress("pi_PZ", &pi_PZ);
      newtree->SetBranchAddress("p_PT", &p_PT);
      newtree->SetBranchAddress("p_PZ", &p_PZ);
      newtree->SetBranchAddress("Xi_PT", &Xi_PT);

      newtree->SetBranchAddress("Jpsi_TRUEP_Z", &Jpsi_TRUEP_Z);
      newtree->SetBranchAddress("Jpsi_TRUEPT", &Jpsi_TRUEPT);
      newtree->SetBranchAddress("Xib_TRUEP_Z", &Xib_TRUEP_Z);
      newtree->SetBranchAddress("Xib_TRUEPT", &Xib_TRUEPT);
      newtree->SetBranchAddress("Xi_TRUEPT", &Xi_TRUEPT);
      newtree->SetBranchAddress("L_TRUEP_Z", &L_TRUEP_Z);
      newtree->SetBranchAddress("L_TRUEPT", &L_TRUEPT);
      newtree->SetBranchAddress("pi_TRUEPT", &pi_TRUEPT);
      newtree->SetBranchAddress("pi_TRUEP_Z", &pi_TRUEP_Z);
      newtree->SetBranchAddress("p_TRUEPT", &p_TRUEPT);
      newtree->SetBranchAddress("p_TRUEP_Z", &p_TRUEP_Z);
      //newtree->SetBranchAddress("eventNumber", &eventNumber);
      

      chainlt->SetBranchAddress("J_psi_1S_TRUEP_Z", &J_psi_1S_TRUEP_Z);
      chainlt->SetBranchAddress("J_psi_1S_TRUEPT", &J_psi_1S_TRUEPT);
      chainlt->SetBranchAddress("Lambda_b0_TRUEP_Z", &Lambda_b0_TRUEP_Z);
      chainlt->SetBranchAddress("Lambda_b0_TRUEPT", &Lambda_b0_TRUEPT);
      chainlt->SetBranchAddress("Lambda0_TRUEP_Z", &Lambda0_TRUEP_Z);
      chainlt->SetBranchAddress("Lambda0_TRUEPT", &Lambda0_TRUEPT);
      chainlt->SetBranchAddress("Lambda_b0_TRUETAU", &Lambda_b0_TRUETAU);
      chainlt->SetBranchAddress("piminus_TRUEPT", &piminus_TRUEPT);
      chainlt->SetBranchAddress("piminus_TRUEP_Z", &piminus_TRUEP_Z);
      chainlt->SetBranchAddress("pplus_TRUEPT", &pplus_TRUEPT);
      chainlt->SetBranchAddress("pplus_TRUEP_Z", &pplus_TRUEP_Z);
      //chainlt->SetBranchAddress("eventNumber", &eventNumber);


      chainxt->SetBranchAddress("J_psi_1S_TRUEP_Z", &J_psi_1S_TRUEP_Z);
      chainxt->SetBranchAddress("J_psi_1S_TRUEPT", &J_psi_1S_TRUEPT);
      chainxt->SetBranchAddress("Xi_bminus_TRUEP_Z", &Xi_bminus_TRUEP_Z);
      chainxt->SetBranchAddress("Xi_bminus_TRUEPT", &Xi_bminus_TRUEPT);
      chainxt->SetBranchAddress("Xi_bminus_TRUETAU", &Xi_bminus_TRUETAU);
      chainxt->SetBranchAddress("Ximinus_TRUEPT", &Ximinus_TRUEPT);
      //chainxt->SetBranchAddress("piminus0_TRUEPT", &piminus0_TRUEPT);
      //chainxt->SetBranchAddress("piminus0_TRUEP_Z", &piminus0_TRUEP_Z);
      chainxt->SetBranchAddress("piminus_TRUEPT", &piminus_TRUEPT);
      chainxt->SetBranchAddress("piminus_TRUEP_Z", &piminus_TRUEP_Z);
      chainxt->SetBranchAddress("pplus_TRUEPT", &pplus_TRUEPT);
      chainxt->SetBranchAddress("pplus_TRUEP_Z", &pplus_TRUEP_Z);
      chainxt->SetBranchAddress("Lambda0_TRUEP_Z", &Lambda0_TRUEP_Z);
      chainxt->SetBranchAddress("Lambda0_TRUEPT", &Lambda0_TRUEPT);      
      //chainxt->SetBranchAddress("eventNumber", &eventNumber);

      return;
      


}


void getDataFiles(){

  if(year=="Run1"){
    TFile *fL_ll = new TFile("/data1/avenkate/JpsiLambda/massdump/data/total/jpsilambda_LL.root");
    TFile *fL_dd = new TFile("/data1/avenkate/JpsiLambda/massdump/data/total/jpsilambda_DD.root");
    treeLL = (TTree*)fL_ll->Get("MyTuple");
    treeDD = (TTree*)fL_dd->Get("MyTuple");
    
    TFile *fXi;
    fXi = new TFile("/data1/sblusk/XibFrag/data/JpsiXi_LooseCutData_Run1.root");
    treeXi = (TTree*)fXi->Get("MyTuple");
  }
  
  if(year=="Run2"){
    TFile *fL = new TFile("/data1/sblusk/XibFrag/data/JpsiL0_LooseCutData_Run2.root");
    treeLL = (TTree*)fL->Get("MyTuple");
    treeDD = (TTree*)fL->Get("MyTuple");
    TFile *fXi;
    fXi = new TFile("/data1/sblusk/XibFrag/data/JpsiXi_LooseCutData_Run2.root");
    treeXi = (TTree*)fXi->Get("MyTuple");
  }  


}


void doLbFitAndPlot(){


  pdf1->fitTo(*rdh1,Hesse(kTRUE),Strategy(1));
  mass = mg1->getVal();
  err_mass = mg1->getError();
  TCanvas* cc1 = makeRooPlot(pdf1,rdh1,mvar1,0,nsig1,modes[0],mLowLb,mHiLb,nbinLb);
  cc1->Draw();
  cc1->Update();
  fitmass_Lb_DD = mg1->getVal();
  err_fitmass_Lb_DD = mg1->getError();
  nsigLb = nsig1->getVal();
  err_nsigLb = nsig1->getError();
  
  TString suff = "DD";
  if(mcFit) suff = suff + "_MC";  

  if(makePDF){
    cc1->Print("./Plots/LbMass_"+suff+".C");
    cc1->Print("./Plots/LbMass_"+suff+".pdf");
    cc1->Print("./Plots/LbMass_"+suff+".eps");
    cc1->Print("./Plots/LbMass_"+suff+".png");
  }

  pdf1->fitTo(*rdh1LL,Hesse(kTRUE),Strategy(1));
  mass = mg1->getVal();
  err_mass = mg1->getError();
  TCanvas* cc1L = makeRooPlot(pdf1,rdh1LL,mvar1,0,nsig1,modes[0]+"LL",mLowLb,mHiLb,nbinLb);
  cc1L->Draw();
  cc1L->Update();
  fitmass_Lb_LL = mg1->getVal();
  err_fitmass_Lb_LL = mg1->getError();

  TString suff = "LL";
  if(mcFit) suff = suff + "_MC";
  if(makePDF){
    cc1L->Print("./Plots/LbMass_"+suff+".C");
    cc1L->Print("./Plots/LbMass_"+suff+".pdf");
    cc1L->Print("./Plots/LbMass_"+suff+".eps");
    cc1L->Print("./Plots/LbMass_"+suff+".png");
  }

}

