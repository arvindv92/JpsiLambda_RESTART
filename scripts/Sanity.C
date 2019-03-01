/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply loose sanity cuts to triggered data/MC, and then separate out into LL and DD files.
   Split by year is only done for data, to allow parallel processing of different years.
   For MC, all years are processed together.
 *********************************/
#include "Sanity.h"

void Sanity(Int_t run, Int_t year, Bool_t isData, Int_t mcType, Bool_t logFlag)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
         mcType = 0 when running over data.
   When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   mcType = 4 for Bu_JpsiX, 5 for Bd_JpsiX
 */
//***MODIFIED TM CUTS REMOVED. I HATE TM*******
/*
   NB: BKGCAT applies to MC, but not data. Read  Steve's analysis and think about exactly what TM condition is to be applied
   Lb_DTF_CTAU_L -> Lb_TAU
   DTFchi2 cut is now 0-50
   PID cuts are applied to reject junk events with PID's of 0 or -1000 for Run1 and reject -1000 for Run2
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
	}

	//Set up logging
	if(isData && logFlag)
	{
		//gROOT->ProcessLine(Form(".> logs/data/JpsiLambda/run%d/sanity_%d_log.txt",run,year));
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/Sanity_%d.txt",
		                             run,year),"w");
	}
	else if(!isData && logFlag)
	{
		//gROOT->ProcessLine(Form(".> logs/mc/JpsiLambda/%s/run%d/sanity_log.txt",folder,run));
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Sanity.txt",
		                             folder,run),"w");
	}
	cout<<"******************************************"<<endl;
	cout<<"==> Starting Sanity: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	Int_t entries_init   = 0, entries_final_LL = 0, entries_final_DD = 0;
	Int_t entries_gen    = 0, Lb_BKGCAT        = 0, p_TRACK_Type     = 0;
	Int_t p_MOTHER_ID    = 0, pi_MOTHER_ID     = 0, L_MOTHER_ID      = 0;
	Int_t L_GD_MOTHER_ID = 0, Jpsi_MOTHER_ID   = 0;

	Float_t eff_excl_LL    = 0.0, eff_excl_LL_err = 0.0;
	Float_t eff_incl_LL    = 0.0, eff_incl_LL_err = 0.0;
	Float_t eff_excl_DD    = 0.0, eff_excl_DD_err = 0.0;
	Float_t eff_incl_DD    = 0.0,eff_incl_DD_err  = 0.0;
	Float_t Lb_ConsLb_chi2 = 0.0, Lb_ConsLb_nDOF  = 0.0;

	Double_t Lb_TAU  = 0.0, Lb_ETA  = 0.0, Lb_PT  = 0.0;
	Double_t pi_PIDp = 0.0, pi_PIDK = 0.0, p_PIDp = 0.0;
	Double_t p_PIDK  = 0.0, L_TAU   = 0.0;
	Bool_t genFlag   = false;

	TCut dtfCut = "", lifetimeCut = "", pidCut = "";
	TCut accCut = "", tmCut       = "";

	TFile *fileIn = nullptr, *fileOut_LL = nullptr, *fileOut_DD = nullptr;
	TTree *treeIn = nullptr, *treeOut_LL = nullptr, *treeOut_DD = nullptr;

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

		fileOut_DD = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_sanity_DD.root",folder,run,part),"RECREATE");
		treeOut_DD = (TTree*)treeIn->CloneTree(0);

		treeIn->SetBranchAddress("Lb_BKGCAT",&Lb_BKGCAT);
		treeIn->SetBranchAddress("p_MC_MOTHER_ID",&p_MOTHER_ID);
		treeIn->SetBranchAddress("pi_MC_MOTHER_ID",&pi_MOTHER_ID);
		treeIn->SetBranchAddress("L_MC_MOTHER_ID",&L_MOTHER_ID);
		treeIn->SetBranchAddress("L_MC_GD_MOTHER_ID",&L_GD_MOTHER_ID);
		treeIn->SetBranchAddress("Jpsi_MC_MOTHER_ID",&Jpsi_MOTHER_ID);

	} // end MC block
	else //Data
	{
		fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_triggered_%d.root",run,year));
		treeIn = (TTree*)fileIn->Get("MyTuple");

		fileOut_LL = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_LL_%d.root",run,year),"RECREATE");
		treeOut_LL = (TTree*)treeIn->CloneTree(0);

		fileOut_DD = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_DD_%d.root",run,year),"RECREATE");
		treeOut_DD = (TTree*)treeIn->CloneTree(0);
	}//end Data block
	 //end setup of input, output

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"LL Output file = "<<fileOut_LL->GetName()<<endl;
	cout<<"DD Output file = "<<fileOut_DD->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("Lb_TAU",&Lb_TAU);
	treeIn->SetBranchAddress("Lb_PT",&Lb_PT);
	treeIn->SetBranchAddress("L_TAU",&L_TAU);
	treeIn->SetBranchAddress("Lb_ConsLb_chi2",&Lb_ConsLb_chi2);
	treeIn->SetBranchAddress("Lb_ConsLb_nDOF",&Lb_ConsLb_nDOF);
	treeIn->SetBranchAddress("Lb_ETA",&Lb_ETA);
	treeIn->SetBranchAddress("p_TRACK_Type",&p_TRACK_Type);

	if(isData)
	{
		treeIn->SetBranchAddress("pi_PIDp",&pi_PIDp);
		treeIn->SetBranchAddress("pi_PIDK",&pi_PIDK);
		treeIn->SetBranchAddress("p_PIDp",&p_PIDp);
		treeIn->SetBranchAddress("p_PIDK",&p_PIDK);
	}
	else
	{
		treeIn->SetBranchAddress("pi_PIDp_corr",&pi_PIDp);
		treeIn->SetBranchAddress("pi_PIDK_corr",&pi_PIDK);
		treeIn->SetBranchAddress("p_PIDp_corr",&p_PIDp);
		treeIn->SetBranchAddress("p_PIDK_corr",&p_PIDK);
	}

	lifetimeCut = "(Lb_TAU > 0)&&(L_TAU > 0)";
	dtfCut = "(Lb_ConsLb_chi2/Lb_ConsLb_nDOF > 0 && Lb_ConsLb_chi2/Lb_ConsLb_nDOF < 50)";
	accCut = "(Lb_ETA > 2 && Lb_ETA < 6 && Lb_PT < 20000)";

	if(run == 1)
	{
		//This also ensures that pi_PIDK != -1000 && p_PIDK != -1000
		pidCut = "(pi_PIDp != 0 && pi_PIDp != -1000 && pi_PIDK != 0 && p_PIDp !=0 && p_PIDp!= -1000 && p_PIDK !=0)";
	}
	else if(run == 2)
	{
		pidCut = "(pi_PIDp != -1000 && p_PIDp!= -1000)";
	}
	// if(!isData)
	// {
	//      if(mcType == 1)
	//      {
	//              // tmCut  = "(abs(p_MOTHER_ID) == 3122 && abs(pi_MOTHER_ID) == 3122 && abs(Jpsi_MOTHER_ID) == 5122 && abs(L_MOTHER_ID) == 5122)";
	//              tmCut = "(abs(Lb_TRUEID == 5122) || (Lb_BKGCAT == 60 && abs(Lb_DTF_M_JpsiLConstr - 5620) < 35))";
	//      }
	//      if(mcType == 2)
	//      {
	//              tmCut  = "(abs(p_MOTHER_ID) == 3122 && abs(pi_MOTHER_ID) == 3122 && abs(Jpsi_MOTHER_ID) == 5122 && abs(L_MOTHER_ID) == 3212 && abs(L_GD_MOTHER_ID) == 5122)";
	//      }
	//      if(mcType == 2)
	//      {
	//              tmCut  = "(abs(p_MOTHER_ID) == 3122 && abs(pi_MOTHER_ID) == 3122 && abs(Jpsi_MOTHER_ID) == 5132 && abs(L_MOTHER_ID) == 3312 && abs(L_GD_MOTHER_ID) == 5122)";
	//      }
	// }

	cout<<"I am making the following sanity cuts, and then separating in LL and DD files. Sit tight"<<endl;

	lifetimeCut.Print();
	dtfCut.Print();
	pidCut.Print();
	accCut.Print();
	//if(!isData) tmCut.Print();

	for(Int_t i = 0; i < entries_init; i++)
	{
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn->GetEntry(i);
		if(Lb_TAU > 0 && L_TAU > 0) //require lifetimes to be positive, universal
		{
			if(Lb_ETA > 2 && Lb_ETA < 6 && Lb_PT < 20000) //acceptance cut, universal, pt < 20000 cut needs to be added here.
			{
				if((Lb_ConsLb_chi2/Lb_ConsLb_nDOF) > 0 && (Lb_ConsLb_chi2/Lb_ConsLb_nDOF) < 50)// DTF chi2 cut, universal
				{
					if((run == 1 && pi_PIDp != 0 && pi_PIDp != -1000 && pi_PIDK != 0 && p_PIDp !=0 && p_PIDp!= -1000 && p_PIDK !=0)|| (run == 2 && p_PIDp != -1000 && pi_PIDp != -1000))
					{
						// if( (isData == 1) || (isData == 0 && abs(p_MOTHER_ID) == 3122 && abs(pi_MOTHER_ID) == 3122 && ((mcType == 1 && abs(Jpsi_MOTHER_ID) == 5122 && abs(L_MOTHER_ID) == 5122) || (mcType == 2 && abs(Jpsi_MOTHER_ID) == 5122 && abs(L_MOTHER_ID) == 3212 && abs(L_GD_MOTHER_ID) == 5122) || (mcType == 3 && abs(Jpsi_MOTHER_ID) == 5132 && abs(L_MOTHER_ID) == 3312 && abs(L_GD_MOTHER_ID) == 5122))))
						// {
						if(p_TRACK_Type == 3)
						{
							treeOut_LL->Fill();
						}
						else
						{
							treeOut_DD->Fill();
						}
						// }
					}
				}
			}
		}
	}

	entries_final_LL = treeOut_LL->GetEntries();
	entries_final_DD = treeOut_DD->GetEntries();
	cout<<"Outgoing LL entries = "<<entries_final_LL<<endl;
	cout<<"Outgoing DD entries = "<<entries_final_DD<<endl;

	if(isData==0)//Efficiency calculation for MC.
	{
		if(entries_init != 0)//Exclusive efficiency calculation
		{
			eff_excl_LL = (Float_t)entries_final_LL*100/entries_init;
			eff_excl_LL_err = sqrt(eff_excl_LL*(100.0-eff_excl_LL)/entries_init);

			eff_excl_DD = (Float_t)entries_final_DD*100/entries_init;
			eff_excl_DD_err = sqrt(eff_excl_DD*(100.0-eff_excl_DD)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<"LL Sanity cuts made with exclusive efficiency = "<<
		        eff_excl_LL<<"% +/- " <<eff_excl_LL_err<<" %"<<endl;
		cout<<"DD Sanity cuts made with exclusive efficiency = "<<
		        eff_excl_DD<<"% +/- " <<eff_excl_DD_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genFile>>entries_gen;//Inclusive efficiency calculation.
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl_LL = (Float_t)entries_final_LL*100/entries_gen;
				eff_incl_LL_err = sqrt(eff_incl_LL*(100.0-eff_incl_LL)/entries_gen);

				eff_incl_DD = (Float_t)entries_final_DD*100/entries_gen;
				eff_incl_DD_err = sqrt(eff_incl_DD*(100.0-eff_incl_DD)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<"LL Sanity cuts made with inclusive efficiency = "<<
			        eff_incl_LL<<"% +/- " <<eff_incl_LL_err<<" %"<<endl;
			cout<<"DD Sanity cuts made with inclusive efficiency = "<<
			        eff_incl_DD<<"% +/- " <<eff_incl_DD_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
	}

	fileOut_LL->Write();
	fileOut_LL->Close();

	fileOut_DD->Write();
	fileOut_DD->Close();

	sw.Stop();
	cout << "==> Sanity is done! Mazel Tov!: "; sw.Print();

	//if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);
	//if(logFlag) gSystem->Exec("cat sanity_log.txt");
}
