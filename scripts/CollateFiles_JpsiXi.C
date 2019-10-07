/********************************
   Author : Aravindhan V.
 *********************************/
#include "CollateFiles_JpsiXi.h"

using namespace std;

void CollateFiles_JpsiXi(Int_t run, Int_t year, Bool_t isData, TChain** h1,
                         TChain** h2, Bool_t testing, Bool_t logFlag)
/*DOC
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   testing = true to run only over a subset of data
 */
{
	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");

	//Set up logging
	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiXi/run%d/Collate_%d_log.txt",
		                             run,year),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiXi/run%d/Collate_log.txt",
		                             run),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting JpsiXi CollateFiles: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	if(isData)
	{
		cout<<"Adding Run "<<run<<" Data ROOT files for "<<year
		    <<" to TChain. Sit tight"<<endl;

		const char *dir = "/data1/avenkate/JpsiXi/data";

		gSystem->Exec(Form("ls %s/LooseLL/%d_Mag*/*/jpsixi.root > fileList/run%dFiles_Xi_%d.txt",
		                   dir,year,run,year));

		TFileCollection fc_run(Form("run%d",run),"",Form("fileList/run%dFiles_Xi_%d.txt",
		                                                 run,year));
		TCollection *run_list = (TCollection*)fc_run.GetList();
		if(testing)
		{
			//add only files from 50 subjobs for testing purposes
			(*h1)->AddFileInfoList(run_list,50);
			(*h2)->AddFileInfoList(run_list,50);
		}
		else
		{
			(*h1)->AddFileInfoList(run_list);
			(*h2)->AddFileInfoList(run_list);
		}
		cout<<"DONE ATTACHING ROOT FILES"<<endl;
	}//end Data block

	else //MC
	{
		cout<<"Collating MC for ";

		gSystem->cd("/data1/avenkate/JpsiXi/mc/LooseLL/Pythia8");

		cout<<"WD = "<<gSystem->pwd()<<endl;

		cout<<"JpsiXi Run "<<run<<endl;

		//No PIDGEN files here because no PID cuts are made
		if(run == 1)
		{
			gSystem->Exec("hadd -f rootFiles/mcFiles/JpsiXi/run1/jpsixi.root 2011*/*/*.root 2012*/*/*.root");
		}
		else if(run == 2)
		{
			gSystem->Exec("hadd -f rootFiles/mcFiles/JpsiXi/run2/jpsixi.root 2015*/*/*.root 2016*/*/*.root");
		}
	}//end MC block

	cout<<"DONE COLLATING"<<endl;

	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
	cout<<"WD = "<<gSystem->pwd()<<endl;

	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiXi/run%d/Trigger_%d_log.txt",
		                             run,year),"a");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiXi/run%d/Trigger_log.txt",
		                             run),"a");
	}
}
