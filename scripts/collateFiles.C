#include "collateFiles.h"

using namespace std;

void collateFiles(Int_t run, Int_t isData, Int_t mcType, TChain** h1, TChain** h2)
/*
run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
isData = 1 for data, 0 for MC
mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
*/
{
    cout<<"Staring to collate"<<endl;
//    gROOT->ProcessLine(".> collate_log.txt");
    if(isData == 1)
    {
        cout<<"PROCESSING DATA"<<endl;
        gSystem->cd("/data1/avenkate/JpsiLambda/massdump/data");
        if(run==1)
        {
            cout<<"Adding Run 1 ROOT files to TChain. Sit tight"<<endl;

            gSystem->Exec("ls 2011_Mag*/*/jpsilambda.root 2012_Mag*/*/jpsilambda.root > run1Files.txt");
            TFileCollection fc_run1("run1","","run1Files.txt");
            TCollection *run1_list = (TCollection*)fc_run1.GetList();
            (*h1)->AddFileInfoList(run1_list);
            (*h2)->AddFileInfoList(run1_list);
        }
        else if(run == 2)
        {
            cout<<"Adding Run 2 ROOT files to TChain. Sit tight"<<endl;

            gSystem->Exec("ls 2015_Mag*/*/jpsilambda.root 2016_Mag*/*/jpsilambda.root 2017_Mag*/*/jpsilambda.root 2018_Mag*/*/jpsilambda.root > run2Files.txt");
            TFileCollection fc_run2("run2","","run2Files.txt");
            TCollection *run2_list = (TCollection*)fc_run2.GetList();
            (*h1)->AddFileInfoList(run2_list);
            (*h2)->AddFileInfoList(run2_list);
        }
        cout<<"DONE ATTACHING ROOT FILES"<<endl;
    }

    if(isData == 0)
    {
        gSystem->cd("/data1/avenkate/JpsiLambda/massdump/mc");
        if(mcType == 1)//Jpsi Lambda MC
        {
            if(run == 1)
            {
                gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda.root JpsiLambda/Pythia8/2011*/*/*.root JpsiLambda/Pythia8/2012*/*/*.root");
            }
            if(run == 2)
            {
                gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda.root JpsiLambda/Pythia8/2015*/*/*.root JpsiLambda/Pythia8/2016*/*/*.root");
            }
        }
        if(mcType == 2)//Jpsi Sigma MC
        {
            if(run == 1)
            {
                gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma.root  JpsiSigma/Pythia8/2012*/*/*.root");
            }
            if(run == 2)
            {
                gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma.root JpsiSigma/Pythia8/2015*/*/*.root JpsiSigma/Pythia8/2016*/*/*.root");
            }
        }
        if(mcType == 3)//Jpsi Xi MC
        {
            if(run == 1)
            {
                gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi.root JpsiXi/Pythia8/2011*/*/*.root JpsiXi/Pythia8/2012*/*/*.root");
            }
            if(run == 2)
            {
                gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi.root JpsiXi/Pythia8/2015*/*/*.root JpsiXi/Pythia8/2016*/*/*.root");
            }
        }
    }
    // gROOT->ProcessLine(".>");
    // gSystem->Exec("cat collate_log.txt");
}
