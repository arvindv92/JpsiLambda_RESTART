/********************************
Author : Aravindhan V.
*********************************/
#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TH1D.h>
#include <TROOT.h>
#include <TSystem.h>
#include "collateFiles.h"

using namespace std;
void trigger(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0)
/*
run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
isData = 1 for data, 0 for MC
mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
*/
{
    gROOT->ProcessLine(".L collateFiles.C++");
    Bool_t logFlag = 0; //This should be 0 only while testing.
    Bool_t collateFlag = 1; //If you don't want to re-collate MC, set this to zero. For example, if only the trigger condition changes.
    TFile *fileout(0);
    TTree *treeout(0), *mytree(0);
    char *triggerCut(0);

    Int_t entries_init = 0, entries_final = 0, entries_gen = 0;
    Bool_t hlt1DiMuonHighMass = 0,hlt1TrackMuon = 0,hlt1TrackAllL0 = 0,hlt2DiMuonDetached = 0;
    Float_t eff_excl = 0.0, eff_excl_err = 0.0,eff_incl = 0.0,eff_incl_err = 0.0;

    if(isData == 1)
    {
        TH1D *lumihist, *lumierrhist;
        Double_t lumi, lumierr;

        TChain *h1 = new TChain("Lb2JpsiLTree/MyTuple");
        TChain *h2 = new TChain("GetIntegratedLuminosity/LumiTuple");
        if(run == 1)
        {
            if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run1/trigger_log.txt");

            collateFiles(run,isData,mcType,&h1,&h2);
            fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_triggered.root","RECREATE");
        }
        else if(run == 2)
        {
            if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run2/trigger_log.txt");

            collateFiles(run,isData,mcType,&h1,&h2);
            fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_triggered.root","RECREATE");
        }

        entries_init = h1->GetEntries();

        h2->Draw("IntegratedLuminosity>>lumihist","","goff");
        h2->Draw("IntegratedLuminosityErr>>lumierrhist","","goff");

        lumihist = (TH1D*)gDirectory->Get("lumihist");
        lumierrhist = (TH1D*)gDirectory->Get("lumierrhist");

        cout<<lumihist->GetMean()<<" "<<lumihist->GetEntries()<<endl;
        lumi = (lumihist->GetMean())*(lumihist->GetEntries());

        lumierr = (lumierrhist->GetMean())*(lumierrhist->GetEntries());

        cout<<"Processing "<<lumi<<" +/- "<<lumierr<< "luminosity"<<endl;

        mytree = (TTree*)h1;
    }

    if(isData == 0)
    {
        TFile *filein(0);
        TTree *treein(0), *treein_gen(0);
        ofstream genfile;//The number of generated MC entries in every MC case will be written out to this file, so that it can be accessed for calculating exclusive efficiencies later

        if(mcType == 1)//Jpsi Lambda
        {
            if(run == 1)
            {
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/trigger_log.txt");
                genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/gen_log.txt");
                cout<<"PROCESSING MC for Run 1 Jpsi Lambda"<<endl;
                if(collateFlag) collateFiles(run,isData,mcType);
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda.root");
                treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");
                treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
                fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_triggered.root","RECREATE");
            }
            if(run == 2)
            {
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/trigger_log.txt");
                genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/gen_log.txt");
                cout<<"PROCESSING MC for Run 2 Jpsi Lambda"<<endl;
                if(collateFlag) collateFiles(run,isData,mcType);
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda.root");
                treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");
                treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
                fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_triggered.root","RECREATE");
            }
        }
        if(mcType == 2)//Jpsi Sigma
        {
            if(run == 1)
            {
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/trigger_log.txt");
                genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/gen_log.txt");
                cout<<"PROCESSING MC for Run 1 JpsiSigma"<<endl;
                if(collateFlag) collateFiles(run,isData,mcType);
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma.root");
                treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");
                treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
                fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_triggered.root","RECREATE");
            }
            if(run == 2)
            {
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/trigger_log.txt");
                genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/gen_log.txt");
                cout<<"PROCESSING MC for Run 2 JpsiSigma"<<endl;
                if(collateFlag) collateFiles(run,isData,mcType);
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma.root");
                treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");
                treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
                fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_triggered.root","RECREATE");
            }
        }
        if(mcType == 3)//Jpsi Xi
        {
            if(run == 1)
            {
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/trigger_log.txt");
                genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/gen_log.txt");
                cout<<"PROCESSING MC for Run 1 JpsiXi"<<endl;
                if(collateFlag) collateFiles(run,isData,mcType);
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi.root");
                treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");
                treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
                fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi_triggered.root","RECREATE");
            }
            if(run == 2)
            {
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/trigger_log.txt");
                genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/gen_log.txt");
                cout<<"PROCESSING MC for Run 2 JpsiXi"<<endl;
                if(collateFlag) collateFiles(run,isData,mcType);
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi.root");
                treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");
                treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
                fileout = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi_triggered.root","RECREATE");
            }
        }
        entries_gen = treein_gen->GetEntries();
        genfile<<entries_gen<<endl;
        genfile.close();
        entries_init = treein->GetEntries();
        mytree = treein;
    }

    cout<<mytree->GetEntries()<<endl;
    treeout = mytree->CloneTree(0);

    cout<<"Incoming entries = "<<entries_init<<endl;

    mytree->SetBranchAddress("Lb_Hlt1DiMuonHighMassDecision_TOS",&hlt1DiMuonHighMass);
    mytree->SetBranchAddress("Lb_Hlt1TrackMuonDecision_TOS",&hlt1TrackMuon);
    mytree->SetBranchAddress("Lb_Hlt1TrackAllL0Decision_TOS",&hlt1TrackAllL0);
    mytree->SetBranchAddress("Lb_Hlt2DiMuonDetachedJPsiDecision_TOS",&hlt2DiMuonDetached);

    triggerCut = (char*)"(Lb_Hlt1DiMuonHighMassDecision_TOS==1||Lb_Hlt1TrackMuonDecision_TOS==1||Lb_Hlt1TrackAllL0Decision_TOS==1)&&(Lb_Hlt2DiMuonDetachedJPsiDecision_TOS==1)";

    cout<<"I am making the following trigger cuts. Sit tight"<<endl;
    cout<<triggerCut<<endl;

    for(Int_t i = 0; i < entries_init; i++)
    {
        //if(i==100000) break;//TESTING ONLY
        if(i%100000 == 0)
        {
            cout<<i<<endl;
        }
        mytree->GetEntry(i);
        if(hlt2DiMuonDetached)
        {
            if(hlt1DiMuonHighMass||hlt1TrackMuon||hlt1TrackAllL0)
            {
                treeout->Fill();
            }
        }
    }

    entries_final = treeout->GetEntries();
    cout<<"Outgoing entries = "<<entries_final<<endl;

    if(isData == 0)
    {
        if(entries_init != 0)
        {
            eff_excl = (Float_t)entries_final*100/entries_init;
            eff_excl_err = sqrt( eff_excl*(100.0-eff_excl)/entries_init);
        }
        cout<<"******************************************"<<endl;
        cout<<"Trigger cut made with exclusive efficiency = "<<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
        cout<<"******************************************"<<endl;

        if(entries_gen != 0)
        {
            eff_incl = (Float_t)entries_final*100/entries_gen;
            eff_incl_err = sqrt( eff_incl*(100.0-eff_incl)/entries_gen);
        }
        cout<<"******************************************"<<endl;
        cout<<"Trigger cut made with inclusive efficiency = "<<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
        cout<<"******************************************"<<endl;
    }

    fileout->cd();
    treeout->Write();
    fileout->Close();

    if(logFlag) gROOT->ProcessLine(".>");
    //if(logFlag) gSystem->Exec("cat trigger_log.txt");
}
