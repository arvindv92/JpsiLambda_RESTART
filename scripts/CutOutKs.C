/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply cuts to remove the Ks0 -> Lambda0 mis-ID reflection.
   Input: Files coming from sanity.
   Output: LL and DD files, separated for zeroTracks and nonZeroTracks.
   zeroTracks contains those events for which TupleToolAllTracks did not add any
   charged tracks to the J/psi. nonZeroTracks contains everything else.
   LL and DD are processed separately.
 *********************************/

#include "CutOutKs.h"

void CutOutKs(Int_t run, Int_t year, Bool_t isData, Int_t mcType, Int_t trackType, Bool_t logFlag)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data.
   When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   mcType = 4 for Bu_JpsiX, 5 for Bd_JpsiX
   trackType = 3 for LL, 5 for DD.
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
		break;
	}
	}
	//Set up logging
	if(isData && logFlag && run == 1)
	{
		//		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/CutOutKs_%d.txt",run,year),"w");
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/CutOutKs.txt",run),"w");
	}
	else if(isData && logFlag && run == 2)
	{
		//		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/CutOutKs_%d.txt",run,year),"w");
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/CutOutKs_%d.txt",run,year),"w");
	}
	else if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/CutOutKs.txt",folder,run),"w");
	}
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting CutOutKs: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileOut_nonZero = nullptr, *fileOut_Zero = nullptr;
	TFile *fileIn          = nullptr;

	TTree *treeOut_nonZero = nullptr, *treeOut_Zero = nullptr;
	TTree *treeIn          = nullptr;

	const char *L_dmCut       = "", *pidCut           = "", *type        = "";
	const char *logFolder     = "", *rootFolder       = "";
	const char *fileName_Zero = "", *fileName_nonZero = "";

	Int_t entries_init       = 0, entries_final = 0, entries_final_nonZero = 0;
	Int_t entries_final_Zero = 0, entries_gen   = 0, nTracks               = 0;

	Float_t eff_excl = 0.0, eff_excl_err = 0.0;
	Float_t eff_incl = 0.0, eff_incl_err = 0.0;

	Double_t L_dm  = 0.0, p_PIDp = 0.0, WMpipi = 0.0;
	Bool_t genFlag = false;
	fstream genFile;

	type = (trackType == 3) ? ("LL") : ("DD");

	//Setup input and output.
	if(!isData)  // MC
	{
		logFolder        = Form("logs/mc/JpsiLambda/%s/run%d",folder,run);
		rootFolder       = Form("rootFiles/mcFiles/JpsiLambda/%s/run%d",folder,run);
		fileName_nonZero = Form("%s_cutoutks_%s_nonZeroTracks.root",part,type);
		fileName_Zero    = Form("%s_cutoutks_%s_ZeroTracks.root",part,type);

		if(!gSystem->AccessPathName((Form("%s/gen_log.txt",logFolder))))
		{
			genFile.open((Form("%s/gen_log.txt",logFolder)));
			genFlag = true;
		}
		else cout<<"*****GenFile not accessible!!*****"<<endl;

		fileIn  = TFile::Open(Form("%s/%s_sanity_%s.root",rootFolder,part,type),"READ");
		treeIn  = (TTree*)fileIn->Get("MyTuple");

		fileOut_nonZero = new TFile(Form("%s/%s",rootFolder,fileName_nonZero),"RECREATE");
		treeOut_nonZero = (TTree*)treeIn->CloneTree(0);
		fileOut_Zero    = new TFile(Form("%s/%s",rootFolder,fileName_Zero),"RECREATE");
		treeOut_Zero    = (TTree*)treeIn->CloneTree(0);
	}
	else   // Data
	{
		rootFolder       = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);

		if(run == 1)
		{
			fileIn  = TFile::Open(Form("%s/jpsilambda_sanity_%s.root",rootFolder,type),"READ");
			treeIn  = (TTree*)fileIn->Get("MyTuple");

			fileName_nonZero = Form("jpsilambda_cutoutks_%s_nonZeroTracks.root",type);
			fileName_Zero    = Form("jpsilambda_cutoutks_%s_ZeroTracks.root",type);
		}
		else if(run == 2)
		{
			fileIn  = TFile::Open(Form("%s/jpsilambda_sanity_%s_%d.root",rootFolder,type,year),"READ");
			treeIn  = (TTree*)fileIn->Get("MyTuple");

			fileName_nonZero = Form("jpsilambda_cutoutks_%s_%d_nonZeroTracks.root",type,year);
			fileName_Zero    = Form("jpsilambda_cutoutks_%s_%d_ZeroTracks.root",type,year);
		}

		fileOut_nonZero = new TFile(Form("%s/%s",rootFolder,fileName_nonZero),"RECREATE");
		treeOut_nonZero = (TTree*)treeIn->CloneTree(0);
		fileOut_Zero    = new TFile(Form("%s/%s",rootFolder,fileName_Zero),"RECREATE");
		treeOut_Zero    = (TTree*)treeIn->CloneTree(0);
	}
	//end setup of input, output
	cout<<"******************************************"<<endl;
	cout<<"==> Starting CutOutKs: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	if(run == 1)
	{
		cout<<"******************************************"<<endl;
		cout<<"Processing Run "<<run<<" "<<type
		    <<((isData) ? (" Data") : (" MC type "))<<mcType<<endl;
		cout<<"******************************************"<<endl;
	}
	else if(run == 2)
	{
		cout<<"******************************************"<<endl;
		cout<<"Processing Run "<<run<<" YEAR "<<year<<" "<<type
		    <<((isData) ? (" Data") : (" MC type "))<<mcType<<endl;
		cout<<"******************************************"<<endl;
	}

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"nonZeroTracks Output file = "<<fileOut_nonZero->GetName()<<endl;
	cout<<"ZeroTracks Output file = "<<fileOut_Zero->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("L_dm",&L_dm);
	treeIn->SetBranchAddress("Lb_DTF_L_WMpipi_JpsiConstr",&WMpipi);
	treeIn->SetBranchAddress("Added_n_Particles",&nTracks);

	if(isData)
	{
		treeIn->SetBranchAddress("p_PIDp",&p_PIDp);
	}
	else
	{
		treeIn->SetBranchAddress("p_PIDp_corr",&p_PIDp);
	}

	L_dmCut = "L_dm < 7.5";        //These cuts are not optimized. RC might pain.
	pidCut  = "p_PIDp > 10";
	cout<<"Making LL file"<<endl;

	cout<<"I am making the following cutoutks cuts only on events within "
	    <<"480-520 MeV in Lb_DTF_L_WMpipi_JpsiConstr. "
	    <<"Events outside this window aren't cut on. Sit tight"<<endl;
	cout<<L_dmCut<<endl;
	cout<<pidCut<<endl;

	for (Int_t i = 0; i < entries_init; i++)
	{
		if(i % 100000 == 0) cout<<i<<endl;
		treeIn->GetEntry(i);
		if(WMpipi < 480 || WMpipi > 520)
		{
			if(nTracks > 0) treeOut_nonZero->Fill();
			else if(nTracks == 0) treeOut_Zero->Fill();
			continue;
		}
		else
		{
			if (WMpipi > 480 && WMpipi < 520)
			{
				if (p_PIDp > 10)
				{
					if(L_dm < 7.5)
					{
						if(nTracks > 0) treeOut_nonZero->Fill();
						else if(nTracks == 0) treeOut_Zero->Fill();
					}
				}
			}
		}
	}
	entries_final_nonZero = treeOut_nonZero->GetEntries();
	entries_final_Zero    = treeOut_Zero->GetEntries();
	entries_final         = entries_final_nonZero + entries_final_Zero;

	cout<<"Outgoing Entries with nonZeroTracks = "<<entries_final_nonZero<<endl;
	cout<<"Outgoing Entries with ZeroTracks = "<<entries_final_Zero<<endl;

	if(!isData)//Calculate inclusive and exclusive efficiencies for MC
	{
		entries_init = treeIn->GetEntries("(Lb_BKGCAT==0||Lb_BKGCAT==50)");
		entries_final_nonZero = treeOut_nonZero->GetEntries("(Lb_BKGCAT==0||Lb_BKGCAT==50)");
		entries_final_Zero    = treeOut_Zero->GetEntries("(Lb_BKGCAT==0||Lb_BKGCAT==50)");
		entries_final         = entries_final_nonZero + entries_final_Zero;

		if(entries_init != 0)
		{
			eff_excl     = (Float_t)entries_final*100/entries_init;
			eff_excl_err = sqrt(eff_excl*(100.0-eff_excl)/entries_init);

			cout<<"******************************************"<<endl;
			cout<<type<<" CutOutKs cuts made with exclusive efficiency = "
			    <<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;
		}
		if(genFlag)
		{
			genFile>>entries_gen;
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl     = (Float_t)entries_final*100/entries_gen;
				eff_incl_err = sqrt(eff_incl*(100.0-eff_incl)/entries_gen);

				cout<<"******************************************"<<endl;
				cout<<type<<" CutOutKs cuts made with inclusive efficiency = "
				    <<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
				cout<<"******************************************"<<endl;
			}
			genFile.close();
		}
	}

	if(!isData)
	{
		treeOut_nonZero->SetAlias("PT","Added_H_PT");
		treeOut_nonZero->SetAlias("MINIPCHI2","psi_1S_H_MINIPCHI2");
		treeOut_nonZero->SetAlias("VCHI2DOF","psi_1S_H_VERTEXCHI2_NEW");
		treeOut_nonZero->SetAlias("IPCHI2","psi_1S_H_IPCHI2_NEW");
		treeOut_nonZero->SetAlias("IP","psi_1S_H_IP_NEW");
		treeOut_nonZero->SetAlias("FD","psi_1S_H_FD_NEW");
		treeOut_nonZero->SetAlias("FDCHI2","psi_1S_H_FDCHI2_NEW");
		treeOut_nonZero->SetAlias("TRACKORIVX_Z","psi_1S_H_VERTEXCHI2_NEW");
		treeOut_nonZero->SetAlias("GHOSTPROB","Added_H_GHOST");
		treeOut_nonZero->SetAlias("TRACKCHI2DOF","Added_H_TRACKCHI2");
	}
	fileOut_nonZero->Write();
	fileOut_nonZero->Close();

	fileOut_Zero->Write();
	fileOut_Zero->Close();

	fileIn->Close();

	sw.Stop();
	cout << "==> CutOutKs is done! Death to Ks0!: "; sw.Print();

	// if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);
}
