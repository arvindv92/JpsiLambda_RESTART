/********************************
   Author : Aravindhan V.
 *********************************/
#include "CollateFiles.h"

using namespace std;

void CollateFiles(Int_t run, Int_t year, Bool_t isData,
                  Int_t mcType, TChain** h1, TChain** h2, Bool_t testing,
                  Bool_t loose, Bool_t logFlag)//defaults provided in header
/*DOC
   run = 1 or 2
   isData = true for data, false for MC

   mcType = 0 when processing data. non zero for processing MC.
   mcType = 1 for Lb -> J/psi Lambda MC
   mcType = 2 for Lb -> J/psi Sigma MC        (reco'd JpsiLambda)
   mcType = 3 for Xib -> Jpsi Xi MC           (reco'd JpsiLambda)
   mcType = 4 for Lb -> J/psi Lambda*(1405)MC (reco'd JpsiLambda)
   mcType = 5 for Lb -> J/psi Lambda*(1520)MC (reco'd JpsiLambda)
   mcType = 6 for Lb -> J/psi Lambda*(1600)MC (reco'd JpsiLambda)
   mcType = 7 forB0 -> Jpsi Ks MC             (reco'd JpsiLambda)
   mcType = 8 for Xib0 -> J/psi Lambda MC
   mcType = 9 for Xib- -> J/psi Xi- MC          (reco'd J/psi Xi-)

   testing = true to run only over a subset of data
   loose   = true to run over data from loose stripping line. Only LL for loose stripping line
   logFlag = true if you want the output to be piped to a log file.
 */
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
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
		folder = "Lst1405";
		part = "lst1405";
		break;
	}
	case 5:
	{
		folder = "Lst1520";
		part = "lst1520";
		break;
	}
	case 6:
	{
		folder = "Lst1600";
		part = "lst1600";
		break;
	}
	case 7:
	{
		folder = "JpsiKs";
		part = "jpsiks";
		break;
	}
	case 8:
	{
		folder = "Xib0";
		part = "xib0";
	}
	default:
	{
		cout<<"$$$ MC Type doesn't match any of the allowed cases. Exiting! $$$"<<endl;
		exit(1);
	}
	}

	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/Collate_%d_log.txt",run,year),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Collate_log.txt",folder,run),"w");
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

		if(run == 1)
		{
			if(mcType!=2)
			{
				gSystem->Exec(Form("hadd -f %s/%s/run1/%s.root %s/2011*/*/*.root %s/2012*/*/*.root",path,folder,part,folder,folder));
				gSystem->Exec(Form("hadd -f %s/%s/run1/%s_pidgen.root %s/2011*_pidgen.root %s/2012*_pidgen.root",path,folder,part,folder,folder));
			}
			else //Lb -> J/psi Sigma Run 1 MC doesn't exist for 2011.
			{
				gSystem->Exec(Form("hadd -f %s/%s/run1/%s.root %s/2012*/*/*.root",path,folder,part,folder));
				gSystem->Exec(Form("hadd -f %s/%s/run1/%s_pidgen.root %s/2012*_pidgen.root",path,folder,part,folder));
			}
		}
		else if(run == 2)
		{
			gSystem->Exec(Form("hadd -f %s/%s/run1/%s.root %s/2015*/*/*.root %s/2016*/*/*.root",path,folder,part,folder,folder));
			gSystem->Exec(Form("hadd -f %s/%s/run2/%s_pidgen.root %s/2015*_pidgen.root %s/2016*_pidgen.root",path,folder,part,folder,folder));
		}
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
