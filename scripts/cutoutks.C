#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TCut.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>

using namespace std;

void cutoutks(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0, Int_t trackType = 3)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{
  TFile *filein(0), *fileout(0);
  TTree *treein(0), *treeout(0);

  TCut L_dmcut, pidcut,masswindow, totalcut;

  Int_t entries_init = 0, entries_final = 0, entries_gen = 0;
  Double_t L_dm = 0.0, p_PIDp = 0.0, Lb_DTF_L_WMpipi_JpsiConstr = 0.0;
  Bool_t logFlag = 0, genFlag = 0;

  fstream genfile;
  masswindow = "Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520";

  const char* type = (trackType == 3)?("LL"):("DD");

  if(trackType == 3) cout<<"Processing LL"<<endl;
  if(trackType == 5) cout<<"Processing DD"<<endl;

  switch(isData)
  {
  case 0: // MC
      switch(mcType)
      {
      case 1: //JpsiLambda
          switch(run)
          {
          case 1:
              if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/cutoutks_log.txt");
              if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/gen_log.txt"))
              {
                  genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/gen_log.txt");
                  genFlag = 1;
              }
              filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_sanity_%s.root",type),"READ");
              fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_cutoutks_%s.root",type),"RECREATE");
              break;

          case 2:
              if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/cutoutks_log.txt");
              if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/gen_log.txt"))
              {
                  genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/gen_log.txt");
                  genFlag = 1;
              }
              filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambd_sanity_%s.root",type),"READ");
              fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_cutoutks_%s.root",type),"RECREATE");
              break;
          }
          break;

      case 2: //JpsiSigma
          switch(run)
          {
          case 1:
              if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/cutoutks_log.txt");
              if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/gen_log.txt"))
              {
                  genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/gen_log.txt");
                  genFlag = 1;
              }
              filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_%s.root",type),"READ");
              fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_cutoutks_%s.root",type),"RECREATE");
              break;

          case 2:
              if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/cutoutks_log.txt");
              if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/gen_log.txt"))
              {
                  genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/gen_log.txt");
                  genFlag = 1;
              }
              filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_%s.root",type),"READ");
              fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_cutoutks_%s.root",type),"RECREATE");
              break;
          }
          break;

      case 3: //JpsiXi
          switch(run)
          {
          case 1:
              if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/cutoutks_log.txt");
              if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/gen_log.txt"))
              {
                  genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/gen_log.txt");
                  genFlag = 1;
              }
              filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_%s.root",type),"READ");
              fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi_cutoutks_%s.root",type),"RECREATE");
              break;

          case 2:
              if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/cutoutks_log.txt");
              if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/gen_log.txt"))
              {
                  genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/gen_log.txt");
                  genFlag = 1;
              }
              filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_%s.root",type),"READ");
              fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi_cutoutks_%s.root",type),"RECREATE");
              break;
          }
          break;
      }
      treein = (TTree*)filein->Get("MyTuple");
      break;

  case 1: // Data
      switch(run)
      {
      case 1:
          if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run1/cutoutks_log.txt");
          filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_%s.root",type),"READ");
          fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_cutoutks_%s.root",type),"RECREATE");
          break;

      case 2:
          if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run2/cutoutks_log.txt");
          filein = TFile::Open(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_%s.root",type),"READ");
          fileout = new TFile(TString::Format("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_cutoutks_%s.root",type),"RECREATE");
          break;
      }
      treein = (TTree*)filein->Get("MyTuple");
      break;
  }
  cout<<"******************************************"<<endl;
  cout<<"Input file = "<<filein->GetName()<<endl;
  cout<<"Output file = "<<fileout->GetName()<<endl;
  cout<<"******************************************"<<endl;

  entries_init = treein->GetEntries();
  cout<<"Incoming entries = "<<entries_init<<endl;

  treein->SetBranchAddress("L_dm",&L_dm);
  treein->SetBranchAddress("p_PIDp",&p_PIDp);
  treein->SetBranchAddress("Lb_DTF_L_WMpipi_JpsiConstr",&Lb_DTF_L_WMpipi_JpsiConstr);
  treeout = (TTree*)treein->CloneTree(0);

  if(trackType == 3)
  {
      L_dmcut = "L_dm < 7.5";
      pidcut = "p_PIDp > 10";
      cout<<"Making LL file"<<endl;
  }
  else if(trackType == 5)
  {
      L_dmcut = "L_dm < 10";
      pidcut = "p_PIDp > 15";
      cout<<"Making DD file"<<endl;
  }

  totalcut = (L_dmcut && pidcut && masswindow) || (!masswindow);

  cout<<"I am making the following cutoutks cuts. Sit tight"<<endl;
  totalcut.Print();

  for (Int_t i = 0; i < entries_init; i++)
  {
      if(i % 100000 == 0) cout<<i<<endl;
      treein->GetEntry(i);
      if(trackType == 3)
      {
          if((L_dm < 7.5 && p_PIDp > 10 && Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520) || !(Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520))
          {
              treeout->Fill();
          }
      }
      else if(trackType == 5)
      {
          if((L_dm < 10 && p_PIDp > 15 && Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520) || !(Lb_DTF_L_WMpipi_JpsiConstr > 480 && Lb_DTF_L_WMpipi_JpsiConstr < 520))
          {
              treeout->Fill();
          }
      }
  }

  entries_final = treeout->GetEntries();
  cout<<"Outgoing Entries = "<<entries_final<<endl;

  if(isData==0)
  {
      if(entries_init != 0)
      {
          if(trackType == 3)
          {
              eff_excl_LL = (Float_t)entries_final_LL*100/entries_init;
              eff_excl_LL_err = sqrt( eff_excl_LL*(100.0-eff_excl_LL)/entries_init);
          }
          else if(trackType == 5)
          {
              eff_excl_DD = (Float_t)entries_final_DD*100/entries_init;
              eff_excl_DD_err = sqrt( eff_excl_DD*(100.0-eff_excl_DD)/entries_init);
          }

      }
      cout<<"******************************************"<<endl;
      cout<<"LL Sanity cuts made with exclusive efficiency = "<<eff_excl_LL<<"% +/- " <<eff_excl_LL_err<<" %"<<endl;
      cout<<"DD Sanity cuts made with exclusive efficiency = "<<eff_excl_DD<<"% +/- " <<eff_excl_DD_err<<" %"<<endl;
      cout<<"******************************************"<<endl;

      if(genFlag)
      {
          genfile>>entries_gen;//NEEDS TO BE TESTED.
          cout<<"Original generated number = "<<entries_gen<<endl;
          if(entries_gen != 0)
          {
              eff_incl_LL = (Float_t)entries_final_LL*100/entries_gen;
              eff_incl_LL_err = sqrt(eff_incl_LL*(100.0-eff_incl_LL)/entries_gen);

              eff_incl_DD = (Float_t)entries_final_DD*100/entries_gen;
              eff_incl_DD_err = sqrt(eff_incl_DD*(100.0-eff_incl_DD)/entries_gen);
          }
          cout<<"******************************************"<<endl;
          cout<<"LL Sanity cuts made with inclusive efficiency = "<<eff_incl_LL<<"% +/- " <<eff_incl_LL_err<<" %"<<endl;
          cout<<"DD Sanity cuts made with inclusive efficiency = "<<eff_incl_DD<<"% +/- " <<eff_incl_DD_err<<" %"<<endl;
          cout<<"******************************************"<<endl;
      }
  }

  gROOT->ProcessLine(".>");
  if(flag==1)
    gSystem->Exec("cat cutoutks_LL_log.txt");
  else if(flag==2)
    gSystem->Exec("cat cutoutks_DD_log.txt");

  fileout->Write();
  fileout->Close();
  filein->Close();
}
