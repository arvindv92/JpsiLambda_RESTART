/********************************
   Author : Aravindhan V.
 *********************************/
#include "CollateFiles.h"

using namespace std;

void CollateFiles(Int_t run, Int_t year, Bool_t isData,
                  Int_t mcType, TChain** h1, TChain** h2, Bool_t testing,
                  Bool_t loose, Bool_t logFlag)//defaults provided in header
/*DOC
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   testing = true to run only over a subset of data
   loose = true to run over data from loose stripping line. Only LL for loose stripping line
 */
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	cout<<"******************************************"<<endl;
	cout<<"==> Starting CollateFiles: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	const char *folder = "", *part = "";
	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Bu_JpsiX";
		part = "bu_jpsix";
		break;
	}
	case 5:
	{
		folder = "Bd_JpsiX";
		part = "bd_jpsix";
		break;
	}
	case 6:
	{
		folder = "Lst1405";
		part = "lst1405";
		break;
	}
	case 7:
	{
		folder = "Lst1520";
		part = "lst1520";
		break;
	}
	case 8:
	{
		folder = "Lst1600";
		part = "lst1600";
		break;
	}
	case 9:
	{
		folder = "chiC1";
		part = "chic1";
		break;
	}
	case 10:
	{
		folder = "JpsiKs";
		part = "jpsiks";
		break;
	}
	case 11:
	{
		folder = "Xib0";
		part = "xib0";
	}
	}

	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/collate_%d_log.txt",run,year),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/collate_log.txt",folder,run),"w");
	}

	if(isData)
	{
		cout<<"Adding Run "<<run<<" Data ROOT files for "<<year<<" to TChain. Sit tight"<<endl;

		const char *dir = "/data1/avenkate/JpsiLambda/massdump/data";
		if(loose)
		{
			gSystem->Exec(Form("ls %s/LooseLL/%d_Mag*/*/jpsilambda.root > run%dFiles_%d.txt",
			                   dir,year,run,year));
		}
		else
		{
			gSystem->Exec(Form("ls %s/%d_Mag*/*/jpsilambda.root > run%dFiles_%d.txt",
			                   dir,year,run,year));
		}
		TFileCollection fc_run(Form("run%d",run),"",Form("run%dFiles_%d.txt",run,year));
		TCollection *run_list = (TCollection*)fc_run.GetList();
		if(testing)
		{
			(*h1)->AddFileInfoList(run_list,50);//add only files from 50 subjobs for testing purposes
			(*h2)->AddFileInfoList(run_list,50);
		}
		else
		{
			(*h1)->AddFileInfoList(run_list);
			(*h2)->AddFileInfoList(run_list);
		}
		cout<<"DONE ATTACHING ROOT FILES"<<endl;
	}//end Data loop

	else //MC
	{
		cout<<"Collating MC for "<<folder<< "run "<<run<<endl;
		if(loose)
		{
			gSystem->cd("/data1/avenkate/JpsiLambda/massdump/mc/LooseLL");
		}
		else
		{
			gSystem->cd("/data1/avenkate/JpsiLambda/massdump/mc");
		}
		cout<<"WD = "<<gSystem->pwd()<<endl;
		const char *path = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda";
		if(mcType!=2)
		{
			if(run == 1)
			{
				// gSystem->Exec(Form("hadd -f %s/JpsiLambda/run1/jpsilambda.root JpsiLambda/2011*/*/*.root JpsiLambda/2012*/*/*.root",path));
				gSystem->Exec(Form("hadd -f %s/%s/run1/%s_pidgen.root %s/2011*_pidgen.root %s/2012*_pidgen.root",path,folder,part,folder,folder));
			}
			else if(run == 2)
			{
				// gSystem->Exec(Form("hadd -f %s/JpsiLambda/run2/jpsilambda.root JpsiLambda/2015*/*/*.root JpsiLambda/2016*/*/*.root",path));
				gSystem->Exec(Form("hadd -f %s/%s/run2/%s_pidgen.root %s/2015*_pidgen.root %s/2016*_pidgen.root",path,folder,part,folder,folder));
			}
		}
		else
		{
			if(run == 1)
			{
				// gSystem->Exec(Form("hadd -f %s/JpsiLambda/run1/jpsilambda.root JpsiLambda/2011*/*/*.root JpsiLambda/2012*/*/*.root",path));
				gSystem->Exec(Form("hadd -f %s/%s/run1/%s_pidgen.root %s/2012*_pidgen.root",path,folder,part,folder));
			}
			else if(run == 2)
			{
				// gSystem->Exec(Form("hadd -f %s/JpsiLambda/run2/jpsilambda.root JpsiLambda/2015*/*/*.root JpsiLambda/2016*/*/*.root",path));
				gSystem->Exec(Form("hadd -f %s/%s/run2/%s_pidgen.root %s/2015*_pidgen.root %s/2016*_pidgen.root",path,folder,part,folder,folder));
			}
		}
		// if(mcType == 1)//Jpsi Lambda MC
		// {
		//      if(run == 1)
		//      {
		//              // gSystem->Exec(Form("hadd -f %s/JpsiLambda/run1/jpsilambda.root JpsiLambda/2011*/*/*.root JpsiLambda/2012*/*/*.root",path));
		//              gSystem->Exec(Form("hadd -f %s/JpsiLambda/run1/jpsilambda_pidgen.root JpsiLambda/2011*_pidgen.root JpsiLambda/2012*_pidgen.root",path));
		//      }
		//      else if(run == 2)
		//      {
		//              // gSystem->Exec(Form("hadd -f %s/JpsiLambda/run2/jpsilambda.root JpsiLambda/2015*/*/*.root JpsiLambda/2016*/*/*.root",path));
		//              gSystem->Exec(Form("hadd -f %s/JpsiLambda/run2/jpsilambda_pidgen.root JpsiLambda/2015*_pidgen.root JpsiLambda/2016*_pidgen.root",path));
		//      }
		// }
		// else if(mcType == 2)//Jpsi Sigma MC
		// {
		//      if(run == 1)
		//      {
		//              // gSystem->Exec(Form("hadd -f %s/JpsiSigma/run1/jpsisigma.root JpsiSigma/2012*/*/*.root",path));
		//              gSystem->Exec(Form("hadd -f %s/JpsiSigma/run1/jpsisigma_pidgen.root JpsiSigma/2012*_pidgen.root",path));
		//      }
		//      else if(run == 2)
		//      {
		//              // gSystem->Exec(Form("hadd -f %s/JpsiSigma/run2/jpsisigma.root JpsiSigma/2015*/*/*.root JpsiSigma/2016*/*/*.root",path));
		//              gSystem->Exec(Form("hadd -f %s/JpsiSigma/run2/jpsisigma_pidgen.root JpsiSigma/2015*_pidgen.root JpsiSigma/2016*_pidgen.root",path));
		//      }
		// }
		// else if(mcType == 3)//Jpsi Xi MC
		// {
		//      if(run == 1)
		//      {
		//              // gSystem->Exec(Form("hadd -f %s/JpsiXi/run1/jpsixi.root JpsiXi/2011*/*/*.root JpsiXi/2012*/*/*.root",path));
		//              gSystem->Exec(Form("hadd -f %s/JpsiXi/run1/jpsixi_pidgen.root JpsiXi/2011*_pidgen.root JpsiXi/2012*_pidgen.root",path));
		//      }
		//      else if(run == 2)
		//      {
		//              // gSystem->Exec(Form("hadd -f %s/JpsiXi/run2/jpsixi.root JpsiXi/2015*/*/*.root JpsiXi/2016*/*/*.root",path));
		//              gSystem->Exec(Form("hadd -f %s/JpsiXi/run2/jpsixi_pidgen.root JpsiXi/2015*_pidgen.root JpsiXi/2016*_pidgen.root",path));
		//      }
		// }
		// else if(mcType == 6||mcType==7||mcType==8||mcType==9||mcType==11) //Lambda* and chiC1 MC
		// {
		//      if(run == 1)
		//      {
		//              gSystem->Exec(Form("hadd -f %s/%s/run1/%s.root %s/2011*/*/*.root %s/2012*/*/*.root",path,folder,part,folder,folder));
		//      }
		//      else if(run == 2)
		//      {
		//              gSystem->Exec(Form("hadd -f %s/%s/run2/%s.root %s/2015*/*/*.root %s/2016*/*/*.root",path,folder,part,folder,folder));
		//      }
		// }
	}//end MC loop
	cout<<"DONE COLLATING"<<endl;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	cout<<"WD = "<<gSystem->pwd()<<endl;

	if(logFlag && isData)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/Trigger_%d_log.txt",
		                             run,year),"a");
	}
	else if(logFlag && !isData)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Trigger_log.txt",
		                             folder,run),"a");
	}
}
