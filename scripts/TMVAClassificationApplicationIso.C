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

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
/****MODIFIED****/
void TMVAClassificationApplicationIso(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0, Int_t trackType = 3, Int_t flag = 1, const char* version = "v1")
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
   flag = 1 when applying on all data, flag = 2 when applying only on signal training sample
   version = "v1" or "v2"
 */
{
	#ifdef __CINT__
	gROOT->ProcessLine( ".O0" );         // turn off optimization in CINT
	#endif

	Bool_t logFlag = false, genFlag = false; //logFlag should be 0 only while testing.
	fstream genFile;
	const char *type = "", *logFileName = "";
	TFile *fileIn(0), *fileOut(0);
	TTree *treeIn(0), *treeOut(0);
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

	logFileName = (trackType == 3) ? ("isolationTraining_LL_log.txt") : ("isolationTraining_DD_log.txt");

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

	switch(isData)
	{
	case 0:                                                                 // MC
		switch(mcType)
		{
		case 1:                                                                 //JpsiLambda
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_%s_withiso.root",run,type),"RECREATE");
			/*  switch(run)
			   {
			   case 1:
			   if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run1/sanity_log.txt");
			   if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run1/gen_log.txt"))
			   {
			   genFile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
			   genFlag = 1;
			   }
			   fileIn = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type));
			   fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_%s_withiso.root",run,type),"RECREATE");
			   break;

			   case 2:
			   if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run2/sanity_log.txt");
			   if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run2/gen_log.txt"))
			   {
			   genFile.open("logs/mc/JpsiLambda/run2/gen_log.txt");
			   genFlag = 1;
			   }
			   fileIn = TFile::Open("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_triggered.root");
			   fileOut_LL = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root","RECREATE");
			   fileOut_DD = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity_DD.root","RECREATE");
			   break;
			   }
			 */
			break;

		case 2: //JpsiSigma
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_cutoutks_%s.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_%s_withiso.root",run,type),"RECREATE");
			/*	switch(run)
			   {
			   case 1:
			   if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run1/sanity_log.txt");
			   if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run1/gen_log.txt"))
			   {
			   genFile.open("logs/mc/JpsiSigma/run1/gen_log.txt");
			   genFlag = 1;
			   }
			   fileIn = TFile::Open("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_triggered.root");
			   fileOut_LL = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_LL.root","RECREATE");
			   fileOut_DD = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_DD.root","RECREATE");
			   break;

			   case 2:
			   if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run2/sanity_log.txt");
			   if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run2/gen_log.txt"))
			   {
			   genFile.open("logs/mc/JpsiSigma/run2/gen_log.txt");
			   genFlag = 1;
			   }
			   fileIn = TFile::Open("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_triggered.root");
			   fileOut_LL = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_LL.root","RECREATE");
			   fileOut_DD = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_DD.root","RECREATE");
			   break;
			   }
			 */
			break;

		case 3: //JpsiXi
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiXi/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_cutoutks_%s.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_%s_withiso.root",run,type),"RECREATE");
			/*	switch(run)
			   {
			   case 1:
			   if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run1/sanity_log.txt");
			   if(!gSystem->AccessPathName("logs/mc/JpsiXi/run1/gen_log.txt"))
			   {
			   genFile.open("logs/mc/JpsiXi/run1/gen_log.txt");
			   genFlag = 1;
			   }
			   fileIn = TFile::Open("rootFiles/mcFiles/JpsiXi/run1/jpsixi_triggered.root");
			   fileOut_LL = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_LL.root","RECREATE");
			   fileOut_DD = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_DD.root","RECREATE");
			   break;

			   case 2:
			   if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run2/sanity_log.txt");
			   if(!gSystem->AccessPathName("logs/mc/JpsiXi/run2/gen_log.txt"))
			   {
			   genFile.open("logs/mc/JpsiXi/run2/gen_log.txt");
			   genFlag = 1;
			   }
			   fileIn = TFile::Open("rootFiles/mcFiles/JpsiXi/run2/jpsixi_triggered.root");
			   fileOut_LL = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_LL.root","RECREATE");
			   fileOut_DD = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_DD.root","RECREATE");
			   break;
			   }
			 */
			break;
		}
		treeIn = (TTree*)fileIn->Get("MyTuple");
		treeOut = (TTree*)treeIn->CloneTree(0);
		break;

	case 1: // Data
		if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%s.txt",run,logFileName));
		if(flag == 1)
		{
			fileIn      = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type));
			fileOut     = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withiso.root",run,type),"RECREATE");
		}
		else if(flag == 2)
		{
			fileIn  = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_withiso.root",run,type),"RECREATE");
		}
		/*	switch(run)
		   {
		   case 1:
		   if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run1/sanity_log.txt");
		   fileIn = TFile::Open("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_triggered.root");
		   fileOut_LL = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_LL.root","RECREATE");
		   fileOut_DD = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_DD.root","RECREATE");
		   break;

		   case 2:
		   if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run2/sanity_log.txt");
		   fileIn = TFile::Open("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_triggered.root");
		   fileOut_LL = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root","RECREATE");
		   fileOut_DD = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_DD.root","RECREATE");
		   break;
		   }
		 */
		treeIn      = (TTree*)fileIn->Get("MyTuple");
		treeOut     = (TTree*)treeIn->CloneTree(0);
		break;
	}
	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

// This loads the library
	TMVA::Tools::Instance();

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;

	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

// Create a set of variables and declare them to the reader
// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t var[6];

	reader->AddVariable( "PT", &var[0] );
	reader->AddVariable( "IPCHI2", &var[1] );
	reader->AddVariable( "VCHI2DOF", &var[2] );
	reader->AddVariable( "MINIPCHI2", &var[3] );

	if(!strncmp(version,"v2",2))
	{
		reader->AddVariable( "GHOSTPROB", &var[4] );
		reader->AddVariable( "TRACKCHI2DOF", &var[5] );
	}

	TString dir        = "weights/";
	TString prefix     = TString::Format("TMVAClassification-isok%s_data_%s",type,version);
	TString methodName = "BDT method";

	TString weightfile = dir + prefix + TString("_") + TString("BDTconf1") + TString(".weights.xml");// change  configuration to conf1 or conf5 depending on which is better
	reader->BookMVA(methodName, weightfile);

// Prepare input tree (this must be replaced by your data source)
// in this example, there is a toy tree with signal and one with background events
// we'll later on use only the "signal" events for the test in this example.
//
/*
   TFile *target(0);
   TString targetname;

   if(!strncmp(version,"v1",2)) {
   switch(flag) {
   case 1:
   targetname = "applicationfiles/TMVApp-isok_LL_BDTconf1_temp.root";
   break;
   case 2:
   targetname = "applicationfiles/TMVApp-isok_DD_BDTconf1_temp.root";
   break;
   case 3:
   targetname = "applicationfiles/TMVApp-isok_LLsig_BDTconf1.root";
   break;
   case 4:
   targetname = "applicationfiles/TMVApp-isok_DDsig_BDTconf1.root";
   break;
   }
   }

   else if(version == "v2") {
   switch(flag) {
   case 1:
   targetname = "applicationfiles/TMVApp-isok_LL_BDTconf1_v2.root";
   break;
   case 2:
   targetname = "applicationfiles/TMVApp-isok_DD_BDTconf1_v2.root";
   break;
   case 3:
   targetname = "applicationfiles/TMVApp-isok_LLsig_BDTconf1_v2.root";
   break;
   case 4:
   targetname = "applicationfiles/TMVApp-isok_DDsig_BDTconf1_v2.root";
   break;
   }
   }
   target  = new TFile( targetname,"RECREATE" );
   TTree *treeOut = new TTree("MyTuple","BDT stuff");
 */
// --- Event loop

// Prepare the event tree
// - here the variable names have to corresponds to your tree
// - you can use the same variables as above which is slightly faster,
//   but of course you can use different ones and copy the values inside the event loop
//
	std::cout << "--- Select signal sample" << std::endl;

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("psi_1S_H_MINIPCHI2",1);
	treeIn->SetBranchStatus("psi_1S_H_IPCHI2_NEW",1);
	treeIn->SetBranchStatus("psi_1S_H_VERTEXCHI2_NEW",1);
	treeIn->SetBranchStatus("Added_H_PT",1);
	treeIn->SetBranchStatus("Added_n_Particles",1);
	treeIn->SetBranchStatus("Added_H_TRACKCHI2",1);
	treeIn->SetBranchStatus("Added_H_GHOST",1);

	const Int_t Max = 200;
	Int_t ntracks;
	Float_t PT[Max], IPCHI2[Max], VCHI2DOF[Max], MINIPCHI2[Max], myTRACKCHI2[Max], GHOST[Max], Spec;
	Double_t BDT[Max];

	treeIn->SetBranchAddress( "Added_H_PT", PT );
	treeIn->SetBranchAddress( "psi_1S_H_IPCHI2_NEW", IPCHI2 );
	treeIn->SetBranchAddress( "psi_1S_H_VERTEXCHI2_NEW", VCHI2DOF );
	treeIn->SetBranchAddress( "psi_1S_H_MINIPCHI2", MINIPCHI2 );
	treeIn->SetBranchAddress( "Added_n_Particles",&ntracks);
	if(flag == 2)
	{
		treeIn->SetBranchAddress( "SW", &Spec );
	}
	if(!strncmp(version,"v2",2))
	{
		treeIn->SetBranchAddress( "Added_H_TRACKCHI2", myTRACKCHI2 );
		treeIn->SetBranchAddress( "Added_H_GHOST", GHOST );
	}

	treeOut->Branch("ntracks",&ntracks,"ntracks/I");
	treeOut->Branch("BDT",BDT,"BDT[ntracks]/D");
	treeOut->Branch("PT",PT,"PT[ntracks]/F");
	treeOut->Branch("IPCHI2",IPCHI2,"IPCHI2[ntracks]/F");
	treeOut->Branch("VCHI2DOF",VCHI2DOF,"VCHI2DOF[ntracks]/F");
	treeOut->Branch("MINIPCHI2",MINIPCHI2,"MINIPCHI2[ntracks]/F");

	if(!strncmp(version,"v2",2)) {
		treeOut->Branch("TRACKCHI2",myTRACKCHI2,"TRACKCHI2[ntracks]/F");
		treeOut->Branch("GHOST",GHOST,"GHOST[ntracks]/F");
	}

	std::cout << "--- Processing: " << treeIn->GetEntries() << " events" << std::endl;
	TStopwatch sw;
	sw.Start();
	for (Long64_t ievt=0; ievt<treeIn->GetEntries(); ievt++) {

		if (ievt%50000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

		treeIn->GetEntry(ievt);
		for(Int_t j=0; j<ntracks; j++)
		{
			var[0] = PT[j];
			var[1] = IPCHI2[j];
			var[2] = VCHI2DOF[j];
			var[3] = MINIPCHI2[j];
			if(!strncmp(version,"v2",2)) {
				var[4] = myTRACKCHI2[j];
				var[5] = GHOST[j];
			}
			// --- Return the MVA outputs and fill into histograms
			BDT[j] = reader->EvaluateMVA(methodName);
		}
		treeOut->Fill();
	}
	sw.Stop();// Get elapsed time
	std::cout << "--- End of event loop: ";
	sw.Print();
	//treeOut->Print();

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	std::cout << "--- Created root file: "<<fileOut->GetName()<<" containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;

	if(logFlag) gROOT->ProcessLine(".>");
}
