#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include <iostream>

using namespace std;

void MASTER_JpsiXi_Run2()
{
  TStopwatch sw;
  sw.Start();

  gSystem->Exec("date");

  Bool_t testing         = false; // when true, analysis will only run over a subset of data
  Bool_t loose           = true;  // when true, analysis will run over data/MC from "loose" stripping line. Only LL
  Int_t run              = 2;     // Run 1 or Run 2?
  Bool_t isData          = false;  // Data or MC?
  Bool_t logFlag         = true;  // set to false only while testing.
  Int_t  year            = 2015;

  /*
  //  Data- Run 2
  cout<<"Processing Run 2"<<endl;
  cout<<"Processing 2015"<<endl;
  year = 2015;
  //  Trigger Cut
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
  gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
                          run, year, isData, testing, loose, logFlag));
  gROOT->Reset();

  //Sanity Cuts
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
  gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
  gROOT->Reset();


  cout<<"Processing 2016"<<endl;
  year = 2016;
  //  Trigger Cut
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
  gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
                          run, year, isData, testing, loose, logFlag));
  gROOT->Reset();

  //Sanity Cuts
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
  gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
  gROOT->Reset();


  cout<<"Processing 2017"<<endl;
  year = 2017;
  //  Trigger Cut
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
  gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
                          run, year, isData, testing, loose, logFlag));
  gROOT->Reset();

  //Sanity Cuts
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
  gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
  gROOT->Reset();


  cout<<"Processing 2018"<<endl;
  year = 2018;
  //  Trigger Cut
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
  gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
                          run, year, isData, testing, loose, logFlag));
  gROOT->Reset();

  //Sanity Cuts
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
  gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
  gROOT->Reset();
*/

  //  Simulation- Run 2
  cout<<"Processing Run 2 simulation"<<endl;
  //  Trigger Cut
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Trigger_JpsiXi.C+");
  gROOT->ProcessLine(Form("Trigger_JpsiXi(%d, %d, %d, %d, %d, %d)",
                          run, year, isData, testing, loose, logFlag));
  gROOT->Reset();

  //Sanity Cuts
  gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
  gROOT->ProcessLine(".L Cuts_JpsiXi.C+");
  gROOT->ProcessLine(Form("Cuts_JpsiXi(%d, %d, %d, %d)", run, year, isData, logFlag));
  gROOT->Reset();

  
  sw.Stop();
  cout << "==> MASTER_JpsiXi is done! Huzzah!: "; sw.Print();
}
