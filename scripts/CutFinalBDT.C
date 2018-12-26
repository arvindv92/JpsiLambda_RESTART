#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <fstream>

using namespace std;
//HAS ERRORS. READ AGAIN
void CutFinalBDT(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0, Int_t trackType = 3, const char* isoVersion = "v1", Int_t isoConf = 1, Int_t bdtConf = 1, Float_t bdtCut = -1.0, Bool_t isoFlag = true, Bool_t logFlag = false)
{
	TFile *fileIn = nullptr, *fileIn_zeroTracks = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeIn_zeroTracks = nullptr, *treeOut = nullptr;

	const char *type = "";
	TString logFileName = "";
	Bool_t genFlag = false;
	fstream genFile;

	if(trackType == 3)
	{
		cout<<"Processing LL"<<endl;
		type = "LL";
	}
	else if(trackType == 5)
	{
		cout<<"Processing DD"<<endl;
		type = "DD";
	}
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	if(isoFlag) logFileName = Form("CutFinalBDT%d_%s_iso%d_%s.txt",bdtConf,type,isoConf,isoVersion);
	else logFileName = Form("CutFinalBDT%d_%s_noIso.txt",bdtConf,type);

	if(!isData) // MC
	{
		switch(mcType)
		{
		case 1:         //JpsiLambda
			if(logFlag)
			{
				gROOT->ProcessLine(Form(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName.Data()));
			}
			if(!gSystem->AccessPathName((Form("logs/mc/JpsiLambda/run%d/gen_log.txt",run))))
			{
				genFile.open((Form("logs/mc/JpsiLambda/run%d/gen_log.txt",run)));
				genFlag = true;
			}
			if(isoFlag)
			{
				fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_%s_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
				treeIn = (TTree*)fileIn->Get("MyTuple");
				fileIn_zeroTracks = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_zeroTracks%s_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
				treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");
				fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/run%d/BDTcut/jpsilambda_%s_BDT%dcut_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"RECREATE");
			}
			else
			{
				fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_%s_FinalBDT%d_noIso.root",run,type,bdtConf),"READ");
				treeIn = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("rootFiles/mcFiles/JpsiLambda/run%d/BDTcut/jpsilambda_%s_BDT%dcut_noIso.root",run,type,bdtConf),"RECREATE");
			}
			treeOut = (TTree*)treeIn->CloneTree(0);
			break;

		case 2:         //JpsiSigma
			if(logFlag)
			{
				gROOT->ProcessLine(Form(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName.Data()));
			}
			if(!gSystem->AccessPathName((Form("logs/mc/JpsiSigma/run%d/gen_log.txt",run))))
			{
				genFile.open((Form("logs/mc/JpsiSigma/run%d/gen_log.txt",run)));
				genFlag = true;
			}
			if(isoFlag)
			{
				fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_%s_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
				treeIn = (TTree*)fileIn->Get("MyTuple");
				fileIn_zeroTracks = TFile::Open(Form("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_zeroTracks%s_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
				treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");
				fileOut = new TFile(Form("rootFiles/mcFiles/JpsiSigma/run%d/BDTcut/jpsisigma_%s_BDT%dcut_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"RECREATE");
			}
			else
			{
				fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_%s_FinalBDT%d_noIso.root",run,type,bdtConf),"READ");
				treeIn = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("rootFiles/mcFiles/JpsiSigma/run%d/BDTcut/jpsisigma_%s_BDT%dcut_noIso.root",run,type,bdtConf),"RECREATE");
			}
			treeOut = (TTree*)treeIn->CloneTree(0);
			break;

		case 3:         //JpsiXi
			if(logFlag)
			{
				gROOT->ProcessLine(Form(".> logs/mc/JpsiXi/run%d/%s",run,logFileName.Data()));
			}
			if(!gSystem->AccessPathName((Form("logs/mc/JpsiXi/run%d/gen_log.txt",run))))
			{
				genFile.open((Form("logs/mc/JpsiXi/run%d/gen_log.txt",run)));
				genFlag = true;
			}
			if(isoFlag)
			{
				fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_%s_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
				treeIn = (TTree*)fileIn->Get("MyTuple");
				fileIn_zeroTracks = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_zeroTracks%s_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
				treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");
				fileOut = new TFile(Form("rootFiles/mcFiles/JpsiXi/run%d/BDTcut/jpsixi_%s_BDT%dcut_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"RECREATE");
			}
			else
			{
				fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_%s_FinalBDT%d_noIso.root",run,type,bdtConf),"READ");
				treeIn = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("rootFiles/mcFiles/JpsiXi/run%d/BDTcut/jpsixi_%s_BDT%dcut_noIso.root",run,type,bdtConf),"RECREATE");
			}
			treeOut = (TTree*)treeIn->CloneTree(0);
			break;
		}
	}
	else   // Data
	{
		if(logFlag)
		{
			gROOT->ProcessLine(Form(".> logs/data/JpsiLambda/run%d/%s",run,logFileName.Data()));
		}
		fileIn = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_sanity_%s_%d.root",run,type,year),"READ");
		treeIn = (TTree*)fileIn->Get("MyTuple");
		fileOut = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_%d.root",run,type,year),"RECREATE");
		treeOut = (TTree*)treeIn->CloneTree(0);
		fileOut_nonZeroTracks = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_%d_nonZeroTracks.root",run,type,year),"RECREATE");
		treeOut_nonZeroTracks = (TTree*)treeIn->CloneTree(0);
		fileOut_ZeroTracks = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_%d_ZeroTracks.root",run,type,year),"RECREATE");
		treeOut_ZeroTracks = (TTree*)treeIn->CloneTree(0);
	}
	cout<<"******************************************"<<endl;
	cout<<"Starting optimizeFinalBDT "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;
}
