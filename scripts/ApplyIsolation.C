/********************************
   Author : Aravindhan V.
 *********************************/
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

//#include "TMVAGui.C"

//#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
//#endif

using namespace TMVA;
using namespace std;

void ApplyIsolation(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0, Int_t trackType = 3, Int_t flag = 1, const char* version = "v1", Int_t isoConf = 1)
/* run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = true for data, false for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
   flag = 1 when applying on all data, flag = 2 when applying only on signal training sample
   version = "v0","v1", "v2" or "v3"
 */
{
	TStopwatch sw;
	sw.Start();

	// #ifdef __CINT__
	// gROOT->ProcessLine( ".O0" );         // turn off optimization in CINT
	// #endif

	const Int_t Max = 200;
	Int_t entries_init = 0, nTracks = 0;
	Float_t PT[Max], IPCHI2[Max], VCHI2DOF[Max], MINIPCHI2[Max], FD[Max], FDCHI2[Max];
	Float_t ipChi2 = 0., vChi2Dof = 0., log_minIpChi2 = 0., log_PT = 0., log_FD = 0., log_fdChi2 = 0.;
	Double_t BDTk[Max], BDTkMin = 1.1;
	Bool_t logFlag = false, genFlag = false; //logFlag should be 0 only while testing.
	const char *type = "", *logFileName = "";
	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;
	fstream genFile;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	cout<<"WD = "<<gSystem->pwd()<<endl;

	logFileName = (trackType == 3) ? ("isoApplication_LL_log.txt") : ("isoApplication_DD_log.txt");

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

	//set up input, output, logging
	if(!isData)  // MC
	{
		switch(mcType)
		{
		case 1:  //JpsiLambda
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_nonZeroTracks.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_%s_withiso_%s.root",run,type,version),"RECREATE");
			break;

		case 2: //JpsiSigma
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_cutoutks_%s_nonZeroTracks.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_%s_withiso_%s.root",run,type,version),"RECREATE");
			break;

		case 3: //JpsiXi
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiXi/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cutoutks_%s_nonZeroTracks.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_%s_withiso_%s.root",run,type,version),"RECREATE");
			break;
		}
	}//end MC block
	else  // Data
	{
		if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%s.txt",run,logFileName));
		if(flag == 1)
		{
			fileIn      = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s_nonZeroTracks.root",run,type));
			fileOut     = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withiso_%s.root",run,type,version),"RECREATE");
		}
		else if(flag == 2)
		{
			fileIn  = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw_nonZeroTracks.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_withiso_%s.root",run,type,version),"RECREATE");
		}
	}
	//end Data block
	//end set up of input, output, logging

	cout<<"******************************************"<<endl;
	cout<< "==> Start of ApplyIsolation"<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	treeIn = (TTree*)fileIn->Get("MyTuple");
	treeOut = (TTree*)treeIn->CloneTree(0);

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	TMVA::Tools::Instance();// This loads the library

	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used

	if (!strncmp(version,"v0",2) || !strncmp(version,"v1",2) || !strncmp(version,"v2",2) || !strncmp(version,"v3",2) )
	{
		reader->AddVariable("IPCHI2", &ipChi2);
		reader->AddVariable("VCHI2DOF", &vChi2Dof);
		reader->AddVariable("log_minIpChi2 := log10(MINIPCHI2)", &log_minIpChi2);
	}
	if(!strncmp(version,"v1",2))
	{
		reader->AddVariable("log_PT := log10(PT)", &log_PT);
	}
	if(!strncmp(version,"v2",2))
	{
		reader->AddVariable("log_FD := log10(FD)", &log_FD);
	}
	if(!strncmp(version,"v3",2))
	{
		reader->AddVariable("log_fdChi2 := log10(FDCHI2)", &log_fdChi2);
	}
	if(!strncmp(version,"v4",2))
	{
		reader->AddVariable("log_PT := log10(PT)", &log_PT);
		reader->AddVariable("log_FD := log10(FD)", &log_FD);
		reader->AddVariable("log_fdChi2 := log10(FDCHI2)", &log_fdChi2);
	}

	TString dir        = "dataset/weights/";
	TString prefix     = TString::Format("TMVAClassification-isok%s_dataRun%d_%s",type,run,version);
	TString methodName = "BDT method";

	TString weightfile = dir + prefix + TString("_") + TString::Format("BDTconf%d",isoConf) + TString(".weights.xml");// change  configuration to conf1 or conf5 depending on which is better
	reader->BookMVA(methodName, weightfile);

	// treeIn->SetBranchStatus("*",0);
	// treeIn->SetBranchStatus("psi_1S_H_MINIPCHI2",1);
	// treeIn->SetBranchStatus("psi_1S_H_IPCHI2_NEW",1);
	// treeIn->SetBranchStatus("psi_1S_H_VERTEXCHI2_NEW",1);
	// treeIn->SetBranchStatus("Added_H_PT",1);
	// treeIn->SetBranchStatus("Added_n_Particles",1);
	// treeIn->SetBranchStatus("psi_1S_H_FD_NEW",1);
	// treeIn->SetBranchStatus("psi_1S_H_FDCHI2_NEW",1);
	// treeIn->SetBranchStatus("Added_H_TRACKCHI2",1);
	// treeIn->SetBranchStatus("Added_H_GHOST",1);

	treeIn->SetBranchAddress( "Added_n_Particles", &nTracks);
	treeIn->SetBranchAddress( "psi_1S_H_IPCHI2_NEW", IPCHI2);
	treeIn->SetBranchAddress( "psi_1S_H_VERTEXCHI2_NEW", VCHI2DOF);
	treeIn->SetBranchAddress( "psi_1S_H_MINIPCHI2", MINIPCHI2);

	if(!strncmp(version,"v1",2) || !strncmp(version,"v4",2))
	{
		treeIn->SetBranchAddress( "Added_H_PT", PT);
	}
	if(!strncmp(version,"v2",2) || !strncmp(version,"v4",2))
	{
		treeIn->SetBranchAddress( "psi_1S_H_FD_NEW", FD);
	}
	if(!strncmp(version,"v3",2) || !strncmp(version,"v4",2))
	{
		treeIn->SetBranchAddress( "psi_1S_H_FDCHI2_NEW", FDCHI2);
	}

	//treeOut->Branch("ntracks",&ntracks,"ntracks/I");
	treeOut->Branch(TString::Format("BDTk_%s",version),BDTk,TString::Format("BDTk_%s[Added_n_Particles]/D",version));
	treeOut->Branch(TString::Format("BDTkMin_%s",version),&BDTkMin,TString::Format("BDTkMin_%s/D",version));

	cout << "--- Processing: " << entries_init << " events" <<endl;


	Int_t tripFlag = 1;
	for (Long64_t ievt = 0; ievt < entries_init; ievt++)
	{
		if (ievt%50000 == 0) cout << "--- ... Processing event: " << ievt << endl;

		treeIn->GetEntry(ievt);
		BDTkMin = 1.1;

		//if(nTracks == 0) BDTkMin = 1.0; //some events don't have any tracks added to them by TTAT

		for(Int_t j = 0; j < nTracks; j++) //Loop over all the added tracks in a given event
		{
			tripFlag = 1;
			ipChi2 = IPCHI2[j];
			vChi2Dof = VCHI2DOF[j];

			(MINIPCHI2[j] > 0) ? (log_minIpChi2 = log10(MINIPCHI2[j])) : (tripFlag = -1);

			if(!strncmp(version,"v1",2) || !strncmp(version,"v4",2))
			{
				(PT[j] > 0) ? (log_PT = log10(PT[j])) : (tripFlag = -1);
			}
			else if(!strncmp(version,"v2",2) || !strncmp(version,"v4",2))
			{
				(FD[j] > 0) ? (log_FD = log10(FD[j])) : (tripFlag = -1);
			}
			else if(!strncmp(version,"v3",2) || !strncmp(version,"v4",2))
			{
				(FDCHI2[j] > 0) ? (log_fdChi2 = log10(FDCHI2[j])) : (tripFlag = -1);
			}
			BDTk[j] = (tripFlag == 1) ? (reader->EvaluateMVA(methodName)) : (1.1);

			if(BDTk[j] < BDTkMin)
			{
				BDTkMin = BDTk[j];
			}
		}
		treeOut->Fill();
	}

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	cout<< "--- Created root file: "<<fileOut->GetName()<<" containing the MVA output histograms" <<endl;

	delete reader;

	sw.Stop();// Get elapsed time
	cout<< "==> End of ApplyIsolation! Isolated forever!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
}
