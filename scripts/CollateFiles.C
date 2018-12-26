/********************************
   Author : Aravindhan V.
 *********************************/
#include "CollateFiles.h"

using namespace std;

void CollateFiles(Int_t run, Int_t year, Bool_t isData, Int_t mcType, TChain** h1, TChain** h2, Bool_t testing, Bool_t loose, Bool_t logFlag)//defaults provided in header
/*DOC
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   testing = true to run only over a subset of data
   loose = true to run over data from loose stripping line. Only LL for loose stripping line
 */
{
	if(isData)
	{
		// if(loose)
		// {
		//      gSystem->cd("/data1/avenkate/JpsiLambda/massdump/data/LooseLL");
		// }
		// else
		// {
		//      gSystem->cd("/data1/avenkate/JpsiLambda/massdump/data");
		// }
		cout<<"WD = "<<gSystem->pwd()<<endl;

		if(logFlag) gROOT->ProcessLine(Form(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run%d/collate_%d_log.txt",run,year));
		cout<<"Adding Run "<<run<<" Data ROOT files for "<<year<<" to TChain. Sit tight"<<endl;

		if(loose)
		{
			gSystem->Exec(Form("ls /data1/avenkate/JpsiLambda/massdump/data/LooseLL/%d_Mag*/*/jpsilambda.root > run%dFiles_%d.txt",year,run,year));
		}
		else
		{
			gSystem->Exec(Form("ls /data1/avenkate/JpsiLambda/massdump/data/%d_Mag*/*/jpsilambda.root > run%dFiles_%d.txt",year,run,year));
		}
		TFileCollection fc_run(Form("run%d",run),"",Form("run%dFiles_%d.txt",run,year));
		TCollection *run_list = (TCollection*)fc_run.GetList();
		if(testing)
		{
			(*h1)->AddFileInfoList(run_list,100);//add only files from 100 subjobs for testing purposes
			(*h2)->AddFileInfoList(run_list,100);
		}
		else
		{
			(*h1)->AddFileInfoList(run_list);
			(*h2)->AddFileInfoList(run_list);
		}

		// if(run == 1)
		// {
		//      if(logFlag) gROOT->ProcessLine(Form(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run1/collate_%d_log.txt",year));
		//      cout<<"Adding Run 1 Data ROOT files for "<<year<<" to TChain. Sit tight"<<endl;
		//
		//      if(loose)
		//      {
		//              gSystem->Exec(Form("ls /data1/avenkate/JpsiLambda/massdump/data/LooseLL/%d_Mag*/*/jpsilambda.root > run1Files_%d.txt",year,year));
		//      }
		//      else
		//      {
		//              gSystem->Exec(Form("ls /data1/avenkate/JpsiLambda/massdump/data/%d_Mag*/*/jpsilambda.root > run1Files_%d.txt",year,year));
		//      }
		//      TFileCollection fc_run1("run1","",Form("run1Files_%d.txt",year));
		//      TCollection *run1_list = (TCollection*)fc_run1.GetList();
		//      if(testing)
		//      {
		//              (*h1)->AddFileInfoList(run1_list,100);//add only files from 100 subjobs for testing purposes
		//              (*h2)->AddFileInfoList(run1_list,100);
		//      }
		//      else
		//      {
		//              (*h1)->AddFileInfoList(run1_list);
		//              (*h2)->AddFileInfoList(run1_list);
		//      }
		// }
		// else if(run == 2)
		// {
		//      if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/data/JpsiLambda/run2/collate_log.txt");
		//      cout<<"Adding Run 2 Data ROOT files to TChain. Sit tight"<<endl;
		//
		//      if(loose)
		//      {
		//              gSystem->Exec("ls /data1/avenkate/JpsiLambda/massdump/data/LooseLL/2015_Mag*/*/jpsilambda.root /data1/avenkate/JpsiLambda/massdump/data/LooseLL/2016_Mag*/*/jpsilambda.root /data1/avenkate/JpsiLambda/massdump/data/LooseLL/2017_Mag*/*/jpsilambda.root /data1/avenkate/JpsiLambda/massdump/data/LooseLL/2018_Mag*/*/jpsilambda.root > run2Files.txt");
		//      }
		//      else
		//      {
		//              gSystem->Exec("ls /data1/avenkate/JpsiLambda/massdump/data/2015_Mag*/*/jpsilambda.root /data1/avenkate/JpsiLambda/massdump/data/2016_Mag*/*/jpsilambda.root	/data1/avenkate/JpsiLambda/massdump/data/2017_Mag*/*/jpsilambda.root /data1/avenkate/JpsiLambda/massdump/data/2018_Mag*/*/jpsilambda.root > run2Files.txt");
		//      }
		//      TFileCollection fc_run2("run2","","run2Files.txt");
		//      TCollection *run2_list = (TCollection*)fc_run2.GetList();
		//      if(testing)
		//      {
		//              (*h1)->AddFileInfoList(run2_list,100);//add only files from 100 subjobs for testing purposes
		//              (*h2)->AddFileInfoList(run2_list,100);
		//      }
		//      else
		//      {
		//              (*h1)->AddFileInfoList(run2_list);
		//              (*h2)->AddFileInfoList(run2_list);
		//      }
		// }
		cout<<"DONE ATTACHING ROOT FILES"<<endl;
	}//end Data loop

	else //MC
	{
		cout<<"Collating MC for ";
		if(loose)
		{
			gSystem->cd("/data1/avenkate/JpsiLambda/massdump/mc/LooseLL");
		}
		else
		{
			gSystem->cd("/data1/avenkate/JpsiLambda/massdump/mc");
		}
		cout<<"WD = "<<gSystem->pwd()<<endl;

		if(mcType == 1)//Jpsi Lambda MC
		{
			cout<<"JpsiLambda Run "<<run<<endl;
			if(run == 1)
			{
				if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run1/collate_log.txt");
				gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run1/jpsilambda.root JpsiLambda/Pythia8/2011*/*/*.root JpsiLambda/Pythia8/2012*/*/*.root");
			}
			else if(run == 2)
			{
				if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiLambda/run2/collate_log.txt");
				gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/run2/jpsilambda.root JpsiLambda/Pythia8/2015*/*/*.root JpsiLambda/Pythia8/2016*/*/*.root");
			}
		}
		else if(mcType == 2)//Jpsi Sigma MC
		{
			cout<<"JpsiSigma Run "<<run<<endl;
			if(run == 1)
			{
				if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run1/collate_log.txt");
				gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run1/jpsisigma.root  JpsiSigma/Pythia8/2012*/*/*.root");
			}
			else if(run == 2)
			{
				if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiSigma/run2/collate_log.txt");
				gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiSigma/run2/jpsisigma.root JpsiSigma/Pythia8/2015*/*/*.root JpsiSigma/Pythia8/2016*/*/*.root");
			}
		}
		else if(mcType == 3)//Jpsi Xi MC
		{
			cout<<"JpsiXi Run "<<run<<endl;
			if(run == 1)
			{
				if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run1/collate_log.txt");
				gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run1/jpsixi.root JpsiXi/Pythia8/2011*/*/*.root JpsiXi/Pythia8/2012*/*/*.root");
			}
			else if(run == 2)
			{
				if(logFlag) gROOT->ProcessLine(".> /data1/avenkate/JpsiLambda_RESTART/logs/mc/JpsiXi/run2/collate_log.txt");
				gSystem->Exec("hadd -f /data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiXi/run2/jpsixi.root JpsiXi/Pythia8/2015*/*/*.root JpsiXi/Pythia8/2016*/*/*.root");
			}
		}
	}//end MC loop
	cout<<"DONE COLLATING"<<endl;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	cout<<"WD = "<<gSystem->pwd()<<endl;

	if(logFlag) gROOT->ProcessLine(".>");
	// gSystem->Exec("cat collate_log.txt");
}
