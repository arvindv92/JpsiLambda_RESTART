/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply loose sanity cuts to triggered data/MC
   Split by year is only done for Run 2 data, to allow parallel processing of different years.
   For MC, all years are processed together.
 *********************************/
#include "Sanity.h"

void Sanity(Int_t run, Int_t year, Bool_t isData, Int_t mcType, Bool_t logFlag)
/*
   run    = 1 or 2.
   isData = 1 for data, 0 for MC

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

   logFlag = true if you want the output to be piped to a log file.
 */
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	const char *folder = "", *part = "";
	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part   = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part   = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part   = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part   = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Lst1405";
		part   = "lst1405";
		break;
	}
	case 5:
	{
		folder = "Lst1520";
		part   = "lst1520";
		break;
	}
	case 6:
	{
		folder = "Lst1600";
		part   = "lst1600";
		break;
	}
	case 7:
	{
		folder = "JpsiKs";
		part   = "jpsiks";
		break;
	}
	case 8:
	{
		folder = "Xib0";
		part   = "xib0";
		break;
	}
	default:
	{
		cout<<"$$$ MC Type doesn't match any of the allowed cases. Exiting! $$$"<<endl;
		exit(1);
	}
	}

	//Set up logging
	if(isData && logFlag && run == 1)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/Sanity.txt",
		                             run),"w");
	}
	else if(isData && logFlag && run == 2)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/Sanity_%d.txt",
		                             run,year),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Sanity.txt",
		                             folder,run),"w");
	}
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting Sanity: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	Int_t entries_init = 0, entries_final_LL = 0;
	Int_t entries_gen  = 0, Lb_BKGCAT        = 0;

	Float_t Lb_ConsLb_chi2 = 0.0, Lb_ConsLb_nDOF  = 0.0;
	Float_t eff_excl_LL    = 0.0, eff_excl_LL_err = 0.0;
	Float_t eff_incl_LL    = 0.0, eff_incl_LL_err = 0.0;

	Double_t Lb_TAU  = 0.0, Lb_ETA  = 0.0, Lb_PT  = 0.0;
	Double_t L_TAU   = 0.0, L_M = 0.0;

	Bool_t genFlag = false;

	TCut dtfCut   = "", lifetimeCut = "", pidCut = "";
	TCut accCut   = "", tmCut       = "";
	TCut LMassCut = "";

	TFile *fileIn = nullptr, *fileOut_LL = nullptr;
	TTree *treeIn = nullptr, *treeOut_LL = nullptr;

	fstream genFile;//contains no of generated candidates in MC.

	// Set up input, output
	if(!isData) // MC
	{
		if(!gSystem->AccessPathName((Form("logs/mc/JpsiLambda/%s/run%d/gen_log.txt",folder,run))))
		{
			genFile.open((Form("logs/mc/JpsiLambda/%s/run%d/gen_log.txt",folder,run)));
			genFlag = true;
		}
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_triggered.root",folder,run,part));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		fileOut_LL = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_sanity_LL.root",folder,run,part),"RECREATE");
		treeOut_LL = (TTree*)treeIn->CloneTree(0);

		treeIn->SetBranchAddress("Lb_BKGCAT",&Lb_BKGCAT);
	} // end MC block
	else //Data
	{
		if(run == 1)
		{
			fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered.root",run));
			treeIn = (TTree*)fileIn->Get("MyTuple");

			fileOut_LL = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL.root",run),"RECREATE");
			treeOut_LL = (TTree*)treeIn->CloneTree(0);
		}
		else if(run == 2)
		{
			fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered_%d.root",run,year));
			treeIn = (TTree*)fileIn->Get("MyTuple");

			fileOut_LL = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL_%d.root",run,year),"RECREATE");
			treeOut_LL = (TTree*)treeIn->CloneTree(0);
		}
	}//end Data block
	 //end setup of input, output

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"LL Output file = "<<fileOut_LL->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("Lb_TAU",&Lb_TAU);
	treeIn->SetBranchAddress("Lb_PT",&Lb_PT);
	treeIn->SetBranchAddress("L_TAU",&L_TAU);
	treeIn->SetBranchAddress("L_M",&L_M);
	treeIn->SetBranchAddress("Lb_ConsLb_chi2",&Lb_ConsLb_chi2);
	treeIn->SetBranchAddress("Lb_ConsLb_nDOF",&Lb_ConsLb_nDOF);
	treeIn->SetBranchAddress("Lb_ETA",&Lb_ETA);

	lifetimeCut = "(Lb_TAU > 0)&&(L_TAU > 0)";
	dtfCut = "(Lb_ConsLb_chi2/Lb_ConsLb_nDOF > 0 && Lb_ConsLb_chi2/Lb_ConsLb_nDOF < 50)";
	accCut = "(Lb_ETA > 2 && Lb_ETA < 6 && Lb_PT < 20000)";
	LMassCut = "(L_M > 1104 && L_M < 1129)";

	cout<<"I am making the following sanity cuts, and then separating in LL and DD files. Sit tight"<<endl;

	lifetimeCut.Print();
	dtfCut.Print();
	accCut.Print();
	LMassCut.Print();

	for(Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn->GetEntry(i);
		if(L_M > 1104 && L_M < 1129)//new cut added
		{
			if(Lb_TAU > 0 && L_TAU > 0) //require lifetimes to be positive, universal
			{
				if((Lb_ConsLb_chi2/Lb_ConsLb_nDOF) > 0 && (Lb_ConsLb_chi2/Lb_ConsLb_nDOF) < 50)// DTF chi2 cut, universal
				{
					if(Lb_ETA > 2 && Lb_ETA < 6 && Lb_PT < 20000) //acceptance cut, universal, pt < 20000 cut needs to be added here.
					{
						treeOut_LL->Fill();
					}
				}
			}
		}
	}

	entries_final_LL = treeOut_LL->GetEntries();
	cout<<"Outgoing LL entries = "<<entries_final_LL<<endl;

	if(isData==0)//Efficiency calculation for MC.
	{
		entries_init     = treeIn->GetEntries("(Lb_BKGCAT==0||Lb_BKGCAT==50)");
		entries_final_LL = treeOut_LL->GetEntries("(Lb_BKGCAT==0||Lb_BKGCAT==50)");
		if(entries_init != 0)//Exclusive efficiency calculation
		{
			eff_excl_LL = (Float_t)entries_final_LL*100/entries_init;
			eff_excl_LL_err = sqrt(eff_excl_LL*(100.0-eff_excl_LL)/entries_init);

			cout<<"******************************************"<<endl;
			cout<<"LL Sanity cuts made with exclusive efficiency = "<<
			        eff_excl_LL<<"% +/- " <<eff_excl_LL_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}

		if(genFlag)
		{
			genFile>>entries_gen;//Inclusive efficiency calculation.
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl_LL     = (Float_t)entries_final_LL*100/entries_gen;
				eff_incl_LL_err = sqrt(eff_incl_LL*(100.0-eff_incl_LL)/entries_gen);

				cout<<"******************************************"<<endl;
				cout<<"LL Sanity cuts made with inclusive efficiency = "<<
				        eff_incl_LL<<"% +/- " <<eff_incl_LL_err<<" %"<<endl;
				cout<<"******************************************"<<endl;
			}
		}
	}

	fileOut_LL->Write();
	fileOut_LL->Close();

	sw.Stop();
	cout << "==> Sanity is done! Mazel Tov!: "; sw.Print();

	if(logFlag) gSystem->RedirectOutput(0);
}
