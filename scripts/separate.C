#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TCut.h>
#include <TROOT.h>
#include <TSystem.h>

using namespace std;
//DEPRECATED. THE FUNCTION OF THIS SCRIPT HAS BEEN TAKEN OVER Y SANITY.C
void separate(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0)
/*
After some thought, it doesn't make sense to have a separate script that does the separation after applying sanity cuts. It seems inefficient to loop over the data again.
Instead, when looping over the data in sanity.C, I should also check the track type in the same loop and write out to fileout_LL/DD accordingly.
*/

/*
run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
isData = 1 for data, 0 for MC
mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
*/
{
    Bool_t logFlag = 0, genFlag = 0;

//    if(logFlag) gROOT->ProcessLine(".> separate_log.txt");

    TFile *filein(0), *fileout_LL(0), *fileout_DD(0);
    TTree *treein(0), *treeout_LL(0), *treeout_DD(0);

    TCut LLcut, DDcut;

    Int_t entries_init = 0, entries_final_LL = 0, entries_final_DD = 0, entries_gen = 0, p_TRACK_Type = 0;

    fstream genfile;

    switch(isData)
    {
        case 0: // MC
        switch(mcType)
        {
            case 1: //JpsiLambda
            switch(run)
            {
                case 1:
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/separate_log.txt");
                if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/gen_log.txt"))
                {
                    genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/gen_log.txt");
                    genFlag = 1;
                }
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_sanity.root");
                fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_LL.root","RECREATE");
                fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_DD.root","RECREATE");
                break;

                case 2:
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/separate_log.txt");
                if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/gen_log.txt"))
                {
                    genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/gen_log.txt");
                    genFlag = 1;
                }
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity.root");
                fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_LL.root","RECREATE");
                fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_DD.root","RECREATE");
                break;
            }
            break;

            case 2: //JpsiSigma
            switch(run)
            {
                case 1:
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/separate_log.txt");
                if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/gen_log.txt"))
                {
                    genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/gen_log.txt");
                    genFlag = 1;
                }
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity.root");
                fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_LL.root","RECREATE");
                fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_DD.root","RECREATE");
                break;

                case 2:
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/separate_log.txt");
                if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/gen_log.txt"))
                {
                    genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/gen_log.txt");
                    genFlag = 1;
                }
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity.root");
                fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_LL.root","RECREATE");
                fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_DD.root","RECREATE");
                break;
            }
            break;

            case 3: //JpsiXi
            switch(run)
            {
                case 1:
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/separate_log.txt");
                if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/gen_log.txt"))
                {
                    genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/gen_log.txt");
                    genFlag = 1;
                }
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity.root");
                fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi_LL.root","RECREATE");
                fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi_DD.root","RECREATE");
                break;

                case 2:
                if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/separate_log.txt");
                if(!gSystem->AccessPathName("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/gen_log.txt"))
                {
                    genfile.open("/data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/gen_log.txt");
                    genFlag = 1;
                }
                filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity.root");
                fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi_LL.root","RECREATE");
                fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi_DD.root","RECREATE");
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
            if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run1/separate_log.txt");
            filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity.root");
            fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_LL.root","RECREATE");
            fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_DD.root","RECREATE");
            break;

            case 2:
            if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run2/separate_log.txt");
            filein = TFile::Open("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity.root");
            fileout_LL = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_LL.root","RECREATE");
            fileout_DD = new TFile("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_DD.root","RECREATE");
            break;
        }
        treein = (TTree*)filein->Get("MyTuple");
        break;
    }
    cout<<"******************************************"<<endl;
    cout<<"Input file = "<<filein->GetName()<<endl;
    cout<<"LL Output file = "<<fileout_LL->GetName()<<endl;
    cout<<"DD Output file = "<<fileout_DD->GetName()<<endl;
    cout<<"******************************************"<<endl;

    entries_init = treein->GetEntries();
    cout<<"Incoming entries = "<<entries_init<<endl;

    treein->SetBranchAddress("p_TRACK_Type",&p_TRACK_Type);

    treeout_LL = (TTree*)treein->CloneTree(0);
    treeout_DD = (TTree*)treein->CloneTree(0);

    LLcut = "(p_TRACK_Type == 3 && pi_TRACK_Type == 3)";
    DDcut = "(p_TRACK_Type == 5 && pi_TRACK_Type == 5)";

    cout<<"Separating. Sit tight."<<endl;

    for (Int_t i = 0; i < entries_init; i++)
    {
        if (p_TRACK_Type == 3) //This ensures pi_TRACK_Type = 3.
        {
            treeout_LL->Fill();
        }
        else
        {
            treeout_DD->Fill();
        }
    }

    entries_final_LL = treeout_LL->GetEntries();
    cout<<"Outgoing LL entries = "<<entries_final_LL<<endl;

    fileout_LL->Write();
    fileout_LL->Close();

    fileout_DD = new TFile("jpsilambda_DD.root","RECREATE");

    cout<<"Making DD file. Sit tight"<<endl;

    entries_final_DD = treeout_DD->GetEntries();

    cout<<"Outgoing DD entries = "<<entries_final_DD<<endl;

    filein->Close();
    fileout_DD->Write();
    fileout_DD->Close();

    if(logFlag) gROOT->ProcessLine(".>");
    if(logFlag) gSystem->Exec("cat separate_log.txt");
}
