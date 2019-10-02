#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>

using namespace std;
void CutFinalBDT(Int_t run, Bool_t isData, Int_t mcType,
                 const char* isoVersion, Int_t isoConf, Int_t bdtConf_nonZero,
                 Int_t bdtConf_Zero, Float_t bdtCut, Float_t bdtCut_ZeroTracks, Bool_t isoFlag,
                 Bool_t logFlag, const char *FOM, const char *Part)
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");

	const char *folder = "", *part = "";
	const char *type = (trackType == 3) ? ("LL") : ("DD");

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
	}
	const char* logFileName = "";

	// if(isoFlag) logFileName = Form("CutFinalBDT%d_%s_iso%d_%s.txt",bdtConf,type,isoConf,isoVersion);
	// else logFileName = Form("CutFinalBDT%d_%s_noIso.txt",bdtConf,type);
	//
	// //****Set up logging*****
	// if(isData && logFlag)
	// {
	//      gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/%s",run,logFileName),"a");
	// }
	// else if(!isData && logFlag)
	// {
	//      gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/%s",folder,run,logFileName),"a");
	// }
	cout<<"******************************************"<<endl;
	cout<<"==> Starting CutFinalBDT: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr, *fileIn_zeroTracks = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeIn_zeroTracks = nullptr, *treeOut = nullptr;
	TTree *myTree = nullptr;
	Bool_t genFlag = false;

	Float_t eff_excl = 0.0, eff_excl_err = 0.0;
	Float_t eff_incl = 0.0, eff_incl_err = 0.0;

	Int_t nEntries_nonZero = 0, nEntries_zeroTracks = 0, entries_init = 0, entries_gen;
	Int_t entries_final = 0;
	Double_t myBDT = 0., myBDT_zeroTracks = 0.;
	fstream genFile;

	//*****Set up input, output*****
	if(!isData) // MC
	{
		if(!gSystem->AccessPathName((Form("logs/mc/JpsiLambda/%s/run%d/gen_log.txt",
		                                  folder,run))))
		{
			genFile.open((Form("logs/mc/JpsiLambda/%s/run%d/gen_log.txt",
			                   folder,run)));
			genFlag = true;
		}
		if(isoFlag)
		{
			fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_cutoutks_%s_nonZeroTracks.root",
			                          folder,run,part,type));
			treeIn = (TTree*)fileIn->Get("MyTuple");

			// treeIn->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_%s_iso%d_%s.root",
			//                                  folder,run,part,type,isoConf,isoVersion));

			fileIn_zeroTracks = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_cutoutks_%s_ZeroTracks.root",
			                                     folder,run,part,type));
			treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");

			treeIn_zeroTracks->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_zeroTracks%s_FinalBDT%d.root",
			                                            folder,run,part,type,bdtConf_Zero));

			// fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/BDTcut/%s_%s/%s_%s_BDT%dcut_iso%d_%s.root",
			//                          folder,run,FOM,Part,part,type,bdtConf,isoConf,isoVersion),"RECREATE");
			fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/BDTcut/best.root",
			                         folder,run),"RECREATE");
			treeOut = (TTree*)treeIn->CloneTree(0);
			treeIn->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_%s_FinalBDT%d_iso%d_%s.root",
			                                 folder,run,part,type,bdtConf_nonZero,isoConf,isoVersion));

			treeIn->SetBranchAddress(Form("BDT%d",bdtConf_nonZero),&myBDT);
			treeIn_zeroTracks->SetBranchAddress(Form("BDT%d",bdtConf_Zero),&myBDT_zeroTracks);
		}
		// else
		// {
		//      fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_cutoutks_%s.root",
		//                                folder,run,part,type));
		//      treeIn = (TTree*)fileIn->Get("MyTuple");
		//
		//      fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/BDTcut/%s_%s/%s_%s_BDT%dcut_noIso.root",
		//                               folder,run,FOM,Part,part,type,bdtConf),"RECREATE");
		//
		//      treeOut = (TTree*)treeIn->CloneTree(0);
		//      treeIn->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_%s_FinalBDT%d_noIso.root",
		//                                       folder,run,part,type,bdtConf));
		//
		//      treeIn->SetBranchAddress(Form("BDT%d",bdtConf),&myBDT);
		// }

	}
	else   // Data
	{
		if(isoFlag)
		{
			fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_nonZeroTracks.root",
			                          run,type));
			treeIn = (TTree*)fileIn->Get("MyTuple");

			fileIn_zeroTracks = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_ZeroTracks.root",
			                                     run,type));
			treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");
			treeIn_zeroTracks->AddFriend("MyTuple",Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_zeroTracks%s_FinalBDT%d.root",
			                                            run,type,bdtConf_Zero));

			// fileOut = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/BDTcut/%s_%s/jpsilambda_%s_BDT%dcut_iso%d_%s.root",
			//                          run,FOM,Part,type,bdtConf,isoConf,isoVersion),"RECREATE");
			fileOut = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/BDTcut/best.root",
			                         run),"RECREATE");

			treeOut = (TTree*)treeIn->CloneTree(0);

			treeIn->AddFriend("MyTuple",Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_FinalBDT%d_iso%d_%s.root",
			                                 run,type,bdtConf_nonZero,isoConf,isoVersion));

			treeIn->SetBranchAddress(Form("BDT%d",bdtConf_nonZero),&myBDT);
			treeIn_zeroTracks->SetBranchAddress(Form("BDT%d",bdtConf_Zero),&myBDT_zeroTracks);
		}
		// else
		// {
		//      fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",
		//                                run,type));
		//      treeIn = (TTree*)fileIn->Get("MyTuple");
		//
		//      fileOut = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/BDTcut/%s_%s/jpsilambda_%s_BDT%dcut_noIso.root",
		//                               run,FOM,Part,type,bdtConf),"RECREATE");
		//
		//      treeOut = (TTree*)treeIn->CloneTree(0);
		//
		//      treeIn->AddFriend("MyTuple",Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_FinalBDT%d_noIso.root",
		//                                       run,type,bdtConf));
		//      treeIn->SetBranchAddress(Form("BDT%d",bdtConf),&myBDT);
		// }
	}
	nEntries_nonZero = treeIn->GetEntries();
	if(isoFlag) nEntries_zeroTracks = treeIn_zeroTracks->GetEntries();

	entries_init = nEntries_nonZero + nEntries_zeroTracks;

	for (Int_t i = 0; i < nEntries_nonZero; i++)
	{
		//	if(i==200000) break;
		if(i%100000 == 0)
		{
			cout<<i<<endl;
		}
		treeIn->GetEntry(i);

		if(myBDT > bdtCut)
		{
			treeOut->Fill();
		}
	}
	if(isoFlag)
	{
		treeOut->CopyAddresses(treeIn_zeroTracks);
		for (Int_t i = 0; i < nEntries_zeroTracks; i++)
		{
			if(i%100000 == 0)
			{
				cout<<i<<endl;
			}
			treeIn_zeroTracks->GetEntry(i);

			if(myBDT_zeroTracks > bdtCut_ZeroTracks)
			{
				treeOut->Fill();
			}
		}
	}

	entries_final = treeOut->GetEntries();

	cout<<"Outgoing Entries = "<<entries_final<<endl;

	if(!isData && mcType != 4 && mcType != 5)//Calculate inclusive and exclusive efficiencies for MC
	{
		if(entries_init != 0)
		{
			eff_excl     = (Float_t)entries_final*100/entries_init;
			eff_excl_err = sqrt(eff_excl*(100.0-eff_excl)/entries_init);
		}
		cout<<"******************************************"<<endl;
		cout<<type<<" CutFinalBDT cuts made with exclusive efficiency = "
		    <<eff_excl<<"% +/- " <<eff_excl_err<<" %"<<endl;
		cout<<"******************************************"<<endl;

		if(genFlag)
		{
			genFile>>entries_gen;
			cout<<"Original generated number = "<<entries_gen<<endl;
			if(entries_gen != 0)
			{
				eff_incl     = (Float_t)entries_final*100/entries_gen;
				eff_incl_err = sqrt(eff_incl*(100.0-eff_incl)/entries_gen);
			}
			cout<<"******************************************"<<endl;
			cout<<type<<" CutFinalBDT cuts made with inclusive efficiency = "
			    <<eff_incl<<"% +/- " <<eff_incl_err<<" %"<<endl;
			cout<<"******************************************"<<endl;

			genFile.close();
		}
	}

	fileOut->Write();
	fileOut->Close();

	fileIn->Close();
	if(isoFlag) fileIn_zeroTracks->Close();

	sw.Stop();
	cout << "==> CutFinalBDT is done! Death to Background!: "; sw.Print();

	// if(logFlag) gSystem->RedirectOutput(0);
}
