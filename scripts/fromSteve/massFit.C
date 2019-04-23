#include "CodeInput.h"
#include <iostream>
#include <fstream>
using std::ifstream;

using namespace std;
using namespace RooFit;

float mass,err_mass;

bool defFit = true;
bool useCommonMassXib = true;
bool useCommonMassLb = true;

static const int nSample = 6;
static const int nSampleFit = 6;
TTree *treeXib[nSample];
TTree* treeXibNew[nSample];
RooDataSet* dataXib[nSample];
RooRealVar* mvar;
RooRealVar* mvarLb;
RooRealVar* dm[nSample];
RooRealVar* mg[nSample];
RooRealVar* sg_mc[nSample];
RooFormulaVar* sgA[nSample];
RooFormulaVar* sgB[nSample];
RooRealVar* sgScale1[nSample];
RooRealVar* sgScale2[nSample];
RooRealVar* sgScalenSample[nSample];
RooRealVar* alphaA[nSample];
RooRealVar* alphaB[nSample];
RooAbsPdf* gaussA[nSample];
RooAbsPdf* gaussB[nSample];
RooRealVar* fr[nSample];
RooRealVar* nCB;
RooRealVar* bp[nSample] ;
RooAddPdf* signal[nSample];
RooAbsPdf* comb[nSample];
RooRealVar* nsig[nSample];
RooRealVar* ncomb[nSample];
RooAbsPdf* pdf[nSample];
RooAbsReal* nlls[nSample];
RooFormulaVar* sgB[nSample];
RooRealVar* sgL[nSample];
RooRealVar* sgH[nSample];

RooRealVar* mLb[2];
RooFormulaVar* mXib[4];



// Parameters
// 0 = Run 1, Xib--> Jpsi L0(DD) pi(L)
// 1 = Run 2, Xib--> Jpsi L0(DD) pi(L)
// 2 = Run 1, Xib--> Jpsi L0(DD) pi(D)
// 3 = Run 2, Xib--> Jpsi L0(DD) pi(D)
// 4 = Run 1, Lb-->  Jpsi L0(DD) 
// 5 = Run 2, Lb-->  Jpsi L0(DD) 

double sigL[nSample] = { 20.0, 22.0, 24, 20, 20.5, 22.6 };
double sigH[nSample] = { 17.2, 21.0, 21, 18, 20.5, 22.6 };
double frbf[nSample] = {0.85, 0.84, 0.90, 0.83, 0.897,0.892};
double fracCB[nSample] = {0.5, 0.5, 0.5, 0.5, 0.35, 0.40};

  
  

double widmc[nSample] = { 7.1, 7.45, 8.0, 8.4, 7.5, 7.5 };
double alpha_a[nSample] = {0.96, 0.91, 1.14, 1.08, 1.02, 0.99};
double alpha_b[nSample] = {-1.155, -1.05, -1.25, -1.154, -1.31, -1.216};  
double mean[nSample] = {5797, 5797, 5797, 5797, 5619, 5619};
double nsigEst[nSample] = {500, 500, 500, 500, 13000, 13000};
double nbakEst[nSample] = {500, 500, 500, 500, 10000, 10000};
  
int iFit[6] = {1, 1, 1, 1, 1, 1};

int meanOpt = 1; 

void makePlots2(TCanvas *ccc, TString filename){

  
  ccc->Print("./Plots/"+filename+tagCond+".pdf");
  ccc->Print("./Plots/"+filename+tagCond+".png");
  ccc->Print("./Plots/"+filename+tagCond+".eps");
  ccc->Print("./Plots/"+filename+tagCond+".C");
  
}

void getMeans(){

  int k;
  if(meanOpt==1){
    // Lb mass, Run 1, Run 2
    for(int i=0; i<2;i++){
      mLb[i] = new RooRealVar(Form("mLb_%d",i),"",5619, 5610, 5630);
    }
    // Xib mass, Run 1, Run 2
    for(int i=0; i<4;i++){
      dm[i] = new RooRealVar(Form("dm_%d",i),"",178,160,190);
      k = i%2;
      mXib[i] = new RooFormulaVar(Form("mXib_%d",i),"","@0+@1",RooArgList(*dm[i],*mLb[k]));
      cout << "Xib initial mass: " << dm[i]->getVal() << " " << mLb[k]->getVal() << " " << mXib[i]->getVal() << endl;
    }
  }

  return;
}


void getXibMassFitStuff(int i, RooRealVar *mmvar, RooAbsReal* mean){
  
  cout << "i = " << i << "initial mass: " << mean->getVal() << endl;
  int k = i;

  sg_mc[i] = new RooRealVar(Form("sg_mc_%d",i),"",widmc[i]);
  sgScale1[i] = new RooRealVar(Form("sgScale1_%d",i),"",1.0,0.5,2.0);//,0.5,1.5);
  sgA[i] = new RooFormulaVar(Form("sgA_%d",i),"","@0*@1",RooArgList(*sg_mc[i],*sgScale1[i]));
  sgScale2[i] = new RooRealVar(Form("sgScale2_%d",i),"",1.0);//,0.5,4.0);
  sgB[i] = new RooFormulaVar(Form("sgB_%d",i),"","@0*@1",RooArgList(*sgA[i],*sgScale2[i]));
  alphaA[i] = new RooRealVar(Form("alphaA_%d",i),"",alpha_a[i],0.9,3.0);
  alphaB[i] = new RooRealVar(Form("alphaB_%d",i),"",alpha_b[i],-5.0,-0.2);
  fr[i] = new RooRealVar(Form("fr_%d",i),"fraction",fracCB[i],0,1);

  if(defFit){
    //gaussA[i] = new RooCBShape(Form("gaussA_%d",i),"Signal",*mmvar,*mg[k],*sgA[i],*alphaA[i],*nCB);
    //gaussB[i] = new RooCBShape(Form("gaussB_%d",i),"Signal",*mmvar,*mg[k],*sgB[i],*alphaB[i],*nCB);
    gaussA[i] = new RooCBShape(Form("gaussA_%d",i),"Signal",*mmvar,*mean,*sgA[i],*alphaA[i],*nCB);
    gaussB[i] = new RooCBShape(Form("gaussB_%d",i),"Signal",*mmvar,*mean,*sgB[i],*alphaB[i],*nCB);
  }else{
    gaussA[i] = new RooGaussian(Form("gaussA_%d",i),"Signal",*mmvar,*mean,*sgA[i]);
    fr[i]->setVal(frbf[i]);
    sgL[i] = new RooRealVar(Form("sgL_%d",i),"",sigL[i]);
    sgH[i] = new RooRealVar(Form("sgH_%d",i),"",sigH[i]);
    gaussB[i] = new RooBifurGauss(Form("gaussB_%d",i),"Signal",*mmvar,*mean,*sgL[i],*sgH[i]);
  }
  
  
     //gauss1A = new RooGaussian("gauss1A","Signal",*mvar1,*mg1,*sg1A);
    //gauss1B = new RooGaussian("gauss1B","Signal",*mvar1,*mg1,*sg1B);

  signal[i] = new RooAddPdf(Form("signal_%d",i),"Signal",RooArgList(*gaussA[i],*gaussB[i]),RooArgList(*fr[i]));
  bp[i]  = new RooRealVar(Form("bp_%d",i), "b1", -0.0003,-2,2);

  if(!altBackFunction){
    comb[i] = new RooExponential(Form("comb_%d",i),"Back function",*mmvar,*bp[i]);
  }else{
    comb[i] = new RooChebychev(Form("comb_%d",i),"Back function",*mmvar,*bp[i]);
  }

  nsig[i] = new RooRealVar(Form("nsig_%d",i),"",nsigEst[i],0,30000000);
  ncomb[i] = new RooRealVar(Form("ncomb_%d",i),"",nbakEst[i],0,30000000);
  pdf[i] = new RooAddPdf(Form("pdf_%d",i),"",RooArgList(*signal[i],*comb[i]),RooArgList(*nsig[i],*ncomb[i]));

  fr[i]->setConstant();
  if(defFit){
    alphaA[i]->setConstant();
    alphaB[i]->setConstant();
    sgScale2[i]->setConstant();
  }
  

}



void massFit(){
  makePDF = false;
  altBackFunction = false;
  bool makePDF2 = true;
  forPaper = true;
  

  modes[0] = "#Lambda_{b}^{0}#rightarrow#J/#psi#Lambda";
  modes[1] = "#Xi_{b}^{0}#rightarrow#J/#psi#Xi^{#font[122]{-}}";

  // Apply systematic changes, as needed.

  nbinLb = (mHiLb - mLowLb) / bwLb;
  nbinXib = (mHiXib - mLowXib) / bwXib;
  

  // Get data files
  year = "Run1";
  getDataFiles();
  treeXib[0] = treeXi;
  treeXib[4] = treeDD;
  year = "Run2";
  getDataFiles();
  treeXib[1] = treeXi;
  treeXib[5] = treeDD;
  
  // Get cuts
  getCuts();  

  // Get Means for mass fits
  getMeans();
  //return;

  TString fileName = "/data1/sblusk/junk.root";
  TFile *fout = new TFile(fileName,"RECREATE"); 

  // CReate trees with cuts
  TCut longCut = Xib_Cut && XiL_Cut;
  TCut downCut = Xib_Cut && XiD_Cut;
  // Run 1 Lambda_DD pi_L
  treeXibNew[0] = treeXib[0]->CopyTree(longCut);
  // Run  2 Lambda_DD pi_L
  treeXibNew[1] = treeXib[1]->CopyTree(longCut);
  // Run 1 Lambda_DD pi_D
  treeXibNew[2] = treeXib[0]->CopyTree(downCut);
  // Run  2Lambda_DD pi_D
  treeXibNew[3] = treeXib[1]->CopyTree(downCut);
  // Run 1 Lambda_DD pi_D
  treeXibNew[4] = treeXib[4]->CopyTree(LambdabDD_Cut);
  // Run  2Lambda_DD pi_D
  treeXibNew[5] = treeXib[5]->CopyTree(LambdabDD_Cut);

  

  TCanvas *c = new TCanvas("c","",800,800);
  c->Divide(3,2);
  c->cd(1);
  treeXibNew[0]->Draw("Xib_DTF_M_JpsiXiLConstr");
  c->cd(2);
  treeXibNew[1]->Draw("Xib_DTF_M_JpsiXiLConstr");
  c->cd(3);
  treeXibNew[2]->Draw("Xib_DTF_M_JpsiXiLConstr");
  c->cd(4);
  treeXibNew[3]->Draw("Xib_DTF_M_JpsiXiLConstr");
  c->cd(5);
  treeXibNew[4]->Draw("Lb_DTF_M_JpsiLConstr");
  c->cd(6);
  treeXibNew[5]->Draw("Lb_DTF_M_JpsiLConstr");
  
  //return;

  nCB = new RooRealVar("nCB","",10);
  mvar = new RooRealVar("Xib_DTF_M_JpsiXiLConstr","",mLowXib,mHiXib);
  mvarLb = new RooRealVar("Lb_DTF_M_JpsiLConstr","",mLowLb,mHiLb);
  RooArgSet args;
  RooArgSet argsLb;
  args.add(*mvar);
  argsLb.add(*mvarLb);

  
  int k;
  RooArgSet myArgs;
  for(int j=0; j<nSampleFit; j++){
    if(iFit[j]<0.5) continue;
    if(j<4) {
      dataXib[j] = new RooDataSet(Form("dataXib_%d",j),"dataset with reduced vars",treeXibNew[j],args);
      k = j;
      if(useCommonMassXib) k = 0;
      getXibMassFitStuff(j, mvar, mXib[k]);
    }else{
      dataXib[j] = new RooDataSet(Form("dataXib_%d",j),"dataset with reduced vars",treeXibNew[j],argsLb);
      k = j - 4;
      //if(useCommonMassLb) k = 0;
      getXibMassFitStuff(j, mvarLb, mLb[k]);
    }    
    cout << " Sample: " << j << ", number of entries = " << dataXib[j]->sumEntries() << endl;
    nlls[j] = pdf[j]->createNLL(*dataXib[j],Extended());
    myArgs.add(*nlls[j]);
  }
  
  RooAbsReal* nll = new RooAddition("nll","total nll",myArgs);
  RooMinuit* minuit = new RooMinuit(*nll);
  RooFitResult* fit_start = minuit->save("start","start");
  minuit->setErrorLevel(0.5);
  minuit->setStrategy(2);
  minuit->migrad();
  minuit->migrad();
  minuit->hesse();

  TString tags[6];
  tags[0] = "XibDDL_Run1";  
  tags[1] = "XibDDL_Run2";  
  tags[2] = "XibDDD_Run1";  
  tags[3] = "XibDDD_Run2";
  tags[4] = "LbDD_Run1";
  tags[5] = "LbDD_Run2";
  TString file;

  TCanvas* cc[nSample];
  for(int j=0; j<nSampleFit; j++){
    TString tag = Form("year_%d",j);
    if(j<4) cc[j] = makeRooPlot(pdf[j],dataXib[j],mvar,1,nsig[j],tag,mLowXib,mHiXib,nbinXib);
    if(j>=4) cc[j] = makeRooPlot(pdf[j],dataXib[j],mvarLb,0,nsig[j],tag,mLowLb,mHiLb,nbinLb);
    if(makePDF2){
      file = "MassFit_"+tags[j];
      if(useCommonMassLb) file = file + "_constr";
      makePlots2(cc[j], file) ;
    }
  }
  

  return;

}
