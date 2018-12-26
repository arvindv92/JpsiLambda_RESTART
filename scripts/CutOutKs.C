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
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{
	cout<<"***********Starting CutOutKs***********"<<endl;

	TStopwatch sw;
	sw.Start();

	TFile *fileIn = nullptr;
	TFile *fileOut_nonZero = nullptr, *fileOut_Zero = nullptr;
	TTree *treeIn = nullptr;
	TTree *treeOut_nonZero = nullptr, *treeOut_Zero = nullptr;

	const char *L_dmCut = "", *pidCut = "", *type = "";
	const char *logFolder = "", *rootFolder = "", *logFileName = "";
	const char *fileName_Zero = "", *fileName_nonZero = "";
	Int_t entries_init = 0, entries_final = 0, entries_final_nonZero = 0;
	Int_t entries_final_Zero = 0, entries_gen = 0, nTracks = 0;
	Float_t eff_excl = 0.0, eff_excl_err = 0.0;
	Float_t eff_incl = 0.0, eff_incl_err = 0.0;

	Double_t L_dm = 0.0, p_PIDp = 0.0, WMpipi = 0.0;
	Bool_t genFlag = false;
	fstream genFile;

	type        = (trackType == 3) ? ("LL") : ("DD");
	logFileName = Form("cutoutks_%d_%s_log.txt",year,type);

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	cout<<"WD = "<<gSystem->pwd()<<endl;

	//Setup logging, input and output files.
	if(!isData)  // MC
	{
		switch(mcType)
		{
		case 1:                                 //JpsiLambda
			logFolder        = Form("logs/mc/JpsiLambda/run%d",run);
			rootFolder       = Form("rootFiles/mcFiles/JpsiLambda/run%d",run);
			fileName_nonZero = Form("jpsilambda_cutoutks_%s_nonZeroTracks.root",type);
			fileName_Zero    = Form("jpsilambda_cutoutks_%s_ZeroTracks.root",type);

			if(logFlag)
			{
				gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
			}
			if(!gSystem->AccessPathName((Form("%s/gen_log.txt",logFolder))))
			{
				genFile.open((Form("%s/gen_log.txt",logFolder)));
				genFlag = true;
			}
			else cout<<"*****GenFile not accessible!!*****"<<endl;

			fileIn  = TFile::Open(Form("%s/jpsilambda_sanity_%s.root",rootFolder,type),"READ");
			treeIn  = (TTree*)fileIn->Get("MyTuple");
			// fileOut = new TFile(Form("%s/jpsilambda_cutoutks_%s.root",rootFolder,type),"RECREATE");
			// treeOut = (TTree*)treeIn->CloneTree(0);

			fileOut_nonZero = new TFile(Form("%s/%s",rootFolder,fileName_nonZero),"RECREATE");
			treeOut_nonZero = (TTree*)treeIn->CloneTree(0);
			fileOut_Zero    = new TFile(Form("%s/%s",rootFolder,fileName_Zero),"RECREATE");
			treeOut_Zero    = (TTree*)treeIn->CloneTree(0);
			break;

		case 2:                                 //JpsiSigma
			logFolder        = Form("logs/mc/JpsiSigma/run%d",run);
			rootFolder       = Form("rootFiles/mcFiles/JpsiSigma/run%d",run);
			fileName_nonZero = Form("jpsisigma_cutoutks_%s_nonZeroTracks.root",type);
			fileName_Zero    = Form("jpsisigma_cutoutks_%s_ZeroTracks.root",type);
			if(logFlag)
			{
				gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
			}
			if(!gSystem->AccessPathName((Form("%s/gen_log.txt",logFolder))))
			{
				genFile.open((Form("%s/gen_log.txt",logFolder)));
				genFlag = 1;
			}
			else cout<<"*****GenFile not accessible!!*****"<<endl;

			fileIn  = TFile::Open(Form("%s/jpsisigma_sanity_%s.root",rootFolder,type),"READ");
			treeIn  = (TTree*)fileIn->Get("MyTuple");
			// fileOut = new TFile(Form("%s/jpsisigma_cutoutks_%s.root",rootFolder,type),"RECREATE");
			// treeOut = (TTree*)treeIn->CloneTree(0);

			fileOut_nonZero = new TFile(Form("%s/%s",rootFolder,fileName_nonZero),"RECREATE");
			treeOut_nonZero = (TTree*)treeIn->CloneTree(0);
			fileOut_Zero    = new TFile(Form("%s/%s",rootFolder,fileName_Zero),"RECREATE");
			treeOut_Zero    = (TTree*)treeIn->CloneTree(0);
			break;

		case 3:                                 //JpsiXi
			logFolder        = Form("logs/mc/JpsiXi/run%d",run);
			rootFolder       = Form("rootFiles/mcFiles/JpsiXi/run%d",run);
			fileName_nonZero = Form("jpsixi_cutoutks_%s_nonZeroTracks.root",type);
			fileName_Zero    = Form("jpsixi_cutoutks_%s_ZeroTracks.root",type);
			if(logFlag)
			{
				gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
			}
			if(!gSystem->AccessPathName((Form("%s/gen_log.txt",logFolder))))
			{
				genFile.open((Form("%s/gen_log.txt",logFolder)));
				genFlag = true;
			}
			else cout<<"*****GenFile not accessible!!*****"<<endl;

			fileIn  = TFile::Open(Form("%s/jpsixi_sanity_%s.root",rootFolder,type),"READ");
			treeIn  = (TTree*)fileIn->Get("MyTuple");
			// fileOut = new TFile(Form("%s/jpsixi_cutoutks_%s.root",rootFolder,type),"RECREATE");
			// treeOut = (TTree*)treeIn->CloneTree(0);

			fileOut_nonZero = new TFile(Form("%s/%s",rootFolder,fileName_nonZero),"RECREATE");
			treeOut_nonZero = (TTree*)treeIn->CloneTree(0);
			fileOut_Zero    = new TFile(Form("%s/%s",rootFolder,fileName_Zero),"RECREATE");
			treeOut_Zero    = (TTree*)treeIn->CloneTree(0);
			break;
		}
	}
	else   // Data
	{
		logFolder        = Form("logs/data/JpsiLambda/run%d",run);
		rootFolder       = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);
		fileName_nonZero = Form("jpsilambda_cutoutks_%s_%d_nonZeroTracks.root",type,year);
		fileName_Zero    = Form("jpsilambda_cutoutks_%s_%d_ZeroTracks.root",type,year);
		if(logFlag)
		{
			gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
		}
		fileIn  = TFile::Open(Form("%s/jpsilambda_sanity_%s_%d.root",rootFolder,type,year),"READ");
		treeIn  = (TTree*)fileIn->Get("MyTuple");
		// fileOut = new TFile(Form("%s/jpsilambda_cutoutks_%s_%d.root",rootFolder,type,year),"RECREATE");
		// treeOut = (TTree*)treeIn->CloneTree(0);

		fileOut_nonZero = new TFile(Form("%s/%s",rootFolder,fileName_nonZero),"RECREATE");
		treeOut_nonZero = (TTree*)treeIn->CloneTree(0);
		fileOut_Zero    = new TFile(Form("%s/%s",rootFolder,fileName_Zero),"RECREATE");
		treeOut_Zero    = (TTree*)treeIn->CloneTree(0);
	}
	//end setup of input, output, logging
	cout<<"******************************************"<<endl;
	cout<<"==> Starting CutOutKs: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Processing Run "<<run<<" YEAR "<<year<<" "<<type
	    <<((isData) ? (" Data") : (" MC type "))<<mcType<<endl;
	cout<<"******************************************"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	// cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"nonZeroTracks Output file = "<<fileOut_nonZero->GetName()<<endl;
	cout<<"ZeroTracks Output file = "<<fileOut_Zero->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	treeIn->SetBranchAddress("L_dm",&L_dm);
	treeIn->SetBranchAddress("p_PIDp",&p_PIDp);
	treeIn->SetBranchAddress("Lb_DTF_L_WMpipi_JpsiConstr",&WMpipi);
	treeIn->SetBranchAddress("Added_n_Particles",&nTracks);

	if(trackType == 3)
	{
		L_dmCut = "L_dm < 7.5";//These cuts are not optimized. RC might pain.
		pidCut  = "p_PIDp > 10";
		cout<<"Making LL file"<<endl;
	}
	else if(trackType == 5)
	{
		L_dmCut = "L_dm < 10";
		pidCut  = "p_PIDp > 15";
		cout<<"Making DD file"<<endl;
	}

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
			// treeOut->Fill();
			continue;
		}
		else
		{
			if(trackType == 3)
			{
				if (WMpipi > 480 && WMpipi < 520)
				{
					if (p_PIDp > 10)
					{
						if(L_dm < 7.5)
						{
							if(nTracks > 0) treeOut_nonZero->Fill();
							else if(nTracks == 0) treeOut_Zero->Fill();
							// treeOut->Fill();
						}
					}
				}
			}
			if(trackType == 5)
			{
				if (WMpipi > 480 && WMpipi < 520)
				{
					if (p_PIDp > 15)
					{
						if(L_dm < 10)
						{
							if(nTracks > 0) treeOut_nonZero->Fill();
							if(nTracks == 0) treeOut_Zero->Fill();
							// treeOut->Fill();
						}
					}
				}
			}
		}
	}

	// entries_final = treeOut->GetEntries();
	entries_final_nonZero = treeOut_nonZero->GetEntries();
	entries_final_Zero    = treeOut_Zero->GetEntries();
	entries_final         = entries_final_nonZero + entries_final_Zero;
	// cout<<"Outgoing Entries = "<<entries_final<<endl;

	cout<<"Outgoing Entries with nonZeroTracks = "<<entries_final_nonZero<<endl;
	cout<<"Outgoing Entries with ZeroTracks = "<<entries_final_Zero<<endl;

	if(!isData)//Calculate inclusive and exclusive efficiencies for MC
	{
		if(entries_init != 0)
		{
			eff_excl     = (Float_t)entries_final*100/entries_init;
			eff_excl_err = sqrt(eff_excl*(100.0-eff_excl)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<type<<" CutOutKs cuts made with exclusive efficiency = "
		    <<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genFile>>entries_gen;//NEEDS TO BE TESTED.
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl     = (Float_t)entries_final*100/entries_gen;
				eff_incl_err = sqrt(eff_incl*(100.0-eff_incl)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<type<<" CutOutKs cuts made with inclusive efficiency = "
			    <<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;

			genFile.close();
		}
	}

	// fileOut->Write();
	// fileOut->Close();

	fileOut_nonZero->Write();
	fileOut_nonZero->Close();

	fileOut_Zero->Write();
	fileOut_Zero->Close();

	fileIn->Close();

	sw.Stop();
	cout << "==> CutOutKs is done! Death to Ks0!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
}
