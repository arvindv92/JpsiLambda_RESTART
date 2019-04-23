#include "CodeInput.h"
#include <iostream>
#include <fstream>
using std::ifstream;

using namespace std;
using namespace RooFit;

float mass,err_mass;


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



void sig(){
  plotPull = false;
  makePDFPull = false;  

  mcFit = true;
  makePDF = false;
  year = "Run1";
  altFitFunction = false;

  massExpLb = "Lb_M-Jpsi_M-L_M+3096.90+1115.683";
  massExpXib = "Xib_M-Jpsi_M-Xi_M+3096.90+1321.71";
  

  mLowLb=5400;
  mHiLb = 5850;
  mLowXib= 5540;
  mHiXib = 6050;

  mcOpt = 5;
  setConditions();
  //cout << applyPrPiKinWeight << endl;
  //return;
  modes[0] = "#Lambda_{b}^{0}#rightarrow#J/#psi#Lambda";
  modes[1] = "#Xi_{b}^{0}#rightarrow#J/#psi#Xi^{#font[122]{-}}";

  bwLb = 2.0;
  bwXib = 2.0;
  

  nbinLb = (mHiLb - mLowLb) / bwLb;
  nbinXib = (mHiXib - mLowXib) / bwXib;

  getLbMCFilesForPlot();
  


  // Get cuts
  getCuts();
  
  // Fill Lb mass histograms
  fillLbMassMC();

  // Get mass fit functions
  getLbMassFitFunctions();
  

  // Fit and Plot result
  doLbFitAndPlot();

  // Get Xib files
  getXibFiles();

  // Fill Xib mass histograms
  fillXibMassMC();

  // Get Xib mass fit functions
  getXibMassFitFunctions();

  // Fit and Plot result
  doXibFitAndPlot();  

  
  return;

    

}
