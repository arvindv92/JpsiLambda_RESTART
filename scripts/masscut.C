#include <TFile.h>
#include <TTree.h>
#include <TString.h>
void masscut(Int_t run = 1, Int_t trackType = 3, Int_t low = 5200, Int_t high = 6000)
/*
run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
trackType = 3 for LL, 5 for DD.
high and low refer to the mass cut you want to make
*/
{
  Bool_t logFlag = 0;
  const char* logFileName;
  TFile *filein(0), *fileout(0);
  TTree *treein(0), *treeout(0);

  if(trackType == 3)
  {
    cout<<"Processing LL"<<endl;
    logFileName = "masscut_LL_log.txt";
  }
  else if(trackType == 5)
  {
    cout<<"Processing DD"<<endl;
    logFileName = "masscut_DD_log.txt";
  }
  if(logFlag == 1)
  {
    gROOT->ProcessLine((TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName)).Data());
  }

  filein = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type),"READ");
  fileout = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_masscut_%s.root",run,type),"RECREATE");
  treein = (TTree*)filein->Get("MyTuple");
}
