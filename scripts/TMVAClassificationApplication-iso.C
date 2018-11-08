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
void TMVAClassificationApplication_iso(Int_t run = 1, Int_t isData = 1, Int_t mcType = 0, Int_t trackType = 3 /*, Int_t flag = 1*/, TString version = "v1")
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
   flag = 1 when applying on all data, flag = 2 when applying only on signal training sample
   version = "v1" or "v2"
 */
{
	Bool_t logFlag = false, genFlag = false; //logFlag should be 0 only while testing.
	fstream genfile;
	TFile *filein(0), *filein_sig(0), *fileout(0), *fileout_sig(0);
	switch(isData)
	{
	case 0: // MC
		switch(mcType)
		{
		case 1: //JpsiLambda
			switch(run)
			{
			case 1:
				if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run1/sanity_log.txt");
				if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run1/gen_log.txt"))
				{
					genfile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
					genFlag = 1;
				}
				filein = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type));
				fileout = new TFile("rootFiles/mcFiles/JpsiLambda/run1/jpsilambda_LL_withiso.root","RECREATE");
				break;

			case 2:
				if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiLambda/run2/sanity_log.txt");
				if(!gSystem->AccessPathName("logs/mc/JpsiLambda/run2/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiLambda/run2/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_triggered.root");
				fileout_LL = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root","RECREATE");
				fileout_DD = new TFile("rootFiles/mcFiles/JpsiLambda/run2/jpsilambda_sanity_DD.root","RECREATE");
				break;
			}
			break;

		case 2: //JpsiSigma
			switch(run)
			{
			case 1:
				if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run1/sanity_log.txt");
				if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run1/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiSigma/run1/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_triggered.root");
				fileout_LL = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_LL.root","RECREATE");
				fileout_DD = new TFile("rootFiles/mcFiles/JpsiSigma/run1/jpsisigma_sanity_DD.root","RECREATE");
				break;

			case 2:
				if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiSigma/run2/sanity_log.txt");
				if(!gSystem->AccessPathName("logs/mc/JpsiSigma/run2/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiSigma/run2/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_triggered.root");
				fileout_LL = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_LL.root","RECREATE");
				fileout_DD = new TFile("rootFiles/mcFiles/JpsiSigma/run2/jpsisigma_sanity_DD.root","RECREATE");
				break;
			}
			break;

		case 3: //JpsiXi
			switch(run)
			{
			case 1:
				if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run1/sanity_log.txt");
				if(!gSystem->AccessPathName("logs/mc/JpsiXi/run1/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiXi/run1/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open("rootFiles/mcFiles/JpsiXi/run1/jpsixi_triggered.root");
				fileout_LL = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_LL.root","RECREATE");
				fileout_DD = new TFile("rootFiles/mcFiles/JpsiXi/run1/jpsixi_sanity_DD.root","RECREATE");
				break;

			case 2:
				if(logFlag) gROOT->ProcessLine(".> logs/mc/JpsiXi/run2/sanity_log.txt");
				if(!gSystem->AccessPathName("logs/mc/JpsiXi/run2/gen_log.txt"))
				{
					genfile.open("logs/mc/JpsiXi/run2/gen_log.txt");
					genFlag = 1;
				}
				filein = TFile::Open("rootFiles/mcFiles/JpsiXi/run2/jpsixi_triggered.root");
				fileout_LL = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_LL.root","RECREATE");
				fileout_DD = new TFile("rootFiles/mcFiles/JpsiXi/run2/jpsixi_sanity_DD.root","RECREATE");
				break;
			}
			break;
		}
		treein = (TTree*)filein->Get("MyTuple");
		treein->SetBranchAddress("Lb_BKGCAT",&Lb_BKGCAT);
		break;

	case 1: // Data
		switch(run)
		{
		case 1:
			if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run1/sanity_log.txt");
			filein = TFile::Open("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_triggered.root");
			fileout_LL = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_LL.root","RECREATE");
			fileout_DD = new TFile("rootFiles/dataFiles/JpsiLambda/run1/jpsilambda_sanity_DD.root","RECREATE");
			break;

		case 2:
			if(logFlag) gROOT->ProcessLine(".> logs/data/JpsiLambda/run2/sanity_log.txt");
			filein = TFile::Open("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_triggered.root");
			fileout_LL = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_LL.root","RECREATE");
			fileout_DD = new TFile("rootFiles/dataFiles/JpsiLambda/run2/jpsilambda_sanity_DD.root","RECREATE");
			break;
		}
		treein = (TTree*)filein->Get("MyTuple");
		break;
	}

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

  #ifdef __CINT__
	gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
  #endif

	if((flag == 1) && (version=="v1"))
		gROOT->ProcessLine(".> isoapplication_LL_v1_log.txt");
	else if((flag == 2) && (version=="v1"))
		gROOT->ProcessLine(".> isoapplication_DD_v1_log.txt");
	else if((flag == 1) && (version=="v2"))
		gROOT->ProcessLine(".> isoapplication_LL_v2_log.txt");
	else if((flag == 2) && (version=="v2"))
		gROOT->ProcessLine(".> isoapplication_DD_v2_log.txt");

	//---------------------------------------------------------------

	// This loads the library
	TMVA::Tools::Instance();

	// Default MVA methods to be trained + tested

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;

	// --------------------------------------------------------------------------------------------------

	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	Float_t var[6];

	reader->AddVariable( "PT", &var[0] );
	reader->AddVariable( "IPCHI2", &var[1] );
	reader->AddVariable( "VCHI2DOF", &var[2] );
	reader->AddVariable( "MINIPCHI2", &var[3] );

	if(version == "v2") {
		reader->AddVariable( "GHOSTPROB", &var[4] );
		reader->AddVariable( "TRACKCHI2DOF", &var[5] );
	}
	//   reader->AddVariable( "COSTHETA := cos(Lb_H_OPENING)", &var[4] );

	// Spectator variables declared in the training have to be added to the reader, too
	Float_t spec1;

	// --- Book the MVA methods

	TString dir    = "weights/";
	TString prefix;

	if(flag == 1 || flag == 3) {
		if(version == "v1")
			prefix = "TMVAClassification-isokLL_data";
		else if(version == "v2")
			prefix = "TMVAClassification-isokLL_data_v2";
	}
	else if(flag == 2 || flag == 4) {
		if(version == "v1")
			prefix = "TMVAClassification-isokDD_data";
		else if(version =="v2")
			prefix = "TMVAClassification-isokDD_data_v2";
	}

	TString methodName = TString("BDT method");
	TString weightfile = dir + prefix + TString("_") + TString("BDTconf1") + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile );

	UInt_t nbin = 100;

	// Prepare input tree (this must be replaced by your data source)
	// in this example, there is a toy tree with signal and one with background events
	// we'll later on use only the "signal" events for the test in this example.
	//
	TFile *input(0);
	TString fname;

	switch(flag) {
	case 1:
		//     fname = "./jpsilambda_LL_forTMVAisoApp.root";
		fname = "../cutoutks/jpsilambda_LL.root";
		break;
	case 2:
		//     fname = "./jpsilambda_DD_forTMVAisoApp.root";
		fname = "../cutoutks/jpsilambda_DD.root";
		break;
	case 3:
		fname = "sweightstuff/jpsilambda_LL_forTMVAisoTraining.root";
		break;
	case 4:
		fname = "sweightstuff/jpsilambda_DD_forTMVAisoTraining.root";
		break;
	}

	input = TFile::Open( fname, "READ");

	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}
	std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

	TFile *target(0);
	TString targetname;

	if(version == "v1") {
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
	TTree *treeout = new TTree("MyTuple","BDT stuff");

	// --- Event loop

	// Prepare the event tree
	// - here the variable names have to corresponds to your tree
	// - you can use the same variables as above which is slightly faster,
	//   but of course you can use different ones and copy the values inside the event loop
	//
	std::cout << "--- Select signal sample" << std::endl;
	TTree* theTree = (TTree*)input->Get("MyTuple");

	theTree->SetBranchStatus("*",0);
	theTree->SetBranchStatus("psi_1S_H_MINIPCHI2",1);
	theTree->SetBranchStatus("psi_1S_H_IPCHI2_NEW",1);
	theTree->SetBranchStatus("psi_1S_H_VERTEXCHI2_NEW",1);
	theTree->SetBranchStatus("Added_H_PT",1);
	theTree->SetBranchStatus("Added_n_Particles",1);
	theTree->SetBranchStatus("Added_H_TRACKCHI2",1);
	theTree->SetBranchStatus("Added_H_GHOST",1);
	//   theTree->SetBranchStatus("Lb_H_OPENING",1);

	const Int_t Max = 200;

	Int_t ntracks;
	Float_t PT[Max], IPCHI2[Max], VCHI2DOF[Max], MINIPCHI2[Max], ORIVX_Z[Max], TRACKCHI2[Max], GHOST[Max], Spec;
	Double_t BDT[Max];

	// Double_t PT, IPCHI2, VCHI2DOF, MINIPCHI2, Spec;
	// Double_t BDT;

	theTree->SetBranchAddress( "Added_H_PT", PT );
	theTree->SetBranchAddress( "psi_1S_H_IPCHI2_NEW", IPCHI2 );
	theTree->SetBranchAddress( "psi_1S_H_VERTEXCHI2_NEW", VCHI2DOF );
	theTree->SetBranchAddress( "psi_1S_H_MINIPCHI2", MINIPCHI2 );
	//   theTree->SetBranchAddress( "SW", &Spec );
	//   theTree->SetBranchAddress( "Added_n_Tracks",&ntracks);
	theTree->SetBranchAddress( "Added_n_Particles",&ntracks);
	if(version == "v2") {
		theTree->SetBranchAddress( "Added_H_TRACKCHI2", TRACKCHI2 );
		theTree->SetBranchAddress( "Added_H_GHOST", GHOST );
	}
	// theTree->SetBranchAddress( "PT", &PT );
	// theTree->SetBranchAddress( "IPCHI2", &IPCHI2 );
	// theTree->SetBranchAddress( "VCHI2DOF", &VCHI2DOF );
	// theTree->SetBranchAddress( "MINIPCHI2", &MINIPCHI2 );
	// theTree->SetBranchAddress( "SW", &Spec );

	treeout->Branch("ntracks",&ntracks,"ntracks/I");
	treeout->Branch("BDT",BDT,"BDT[ntracks]/D");
	treeout->Branch("PT",PT,"PT[ntracks]/F");
	treeout->Branch("IPCHI2",IPCHI2,"IPCHI2[ntracks]/F");
	treeout->Branch("VCHI2DOF",VCHI2DOF,"VCHI2DOF[ntracks]/F");
	treeout->Branch("MINIPCHI2",MINIPCHI2,"MINIPCHI2[ntracks]/F");

	if(version == "v2") {
		treeout->Branch("TRACKCHI2",TRACKCHI2,"TRACKCHI2[ntracks]/F");
		treeout->Branch("GHOST",GHOST,"GHOST[ntracks]/F");
	}
	//   treeout->Branch("Lb_H_OPENING",Lb_H_OPENING,"Lb_H_OPENING[ntracks]/F");
	//   treeout->Branch("SW",&Spec,"SW/F");

	// treeout->Branch("BDT",&BDT,"BDT/D");
	// treeout->Branch("PT",&PT,"PT/D");
	// treeout->Branch("IPCHI2",&IPCHI2,"IPCHI2/D");
	// treeout->Branch("VCHI2DOF",&VCHI2DOF,"VCHI2DOF/D");
	// treeout->Branch("MINIPCHI2",&MINIPCHI2,"MINIPCHI2/D");
	// treeout->Branch( "SW", &Spec,"SW/D" );

	std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

	std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
	TStopwatch sw;
	sw.Start();
	for (Long64_t ievt=0; ievt<theTree->GetEntries(); ievt++) {

		if (ievt%50000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

		theTree->GetEntry(ievt);
		for(Int_t j=0; j<ntracks; j++)
		{
			var[0] = PT[j];
			var[1] = IPCHI2[j];
			var[2] = VCHI2DOF[j];
			var[3] = MINIPCHI2[j];
			if(version == "v2") {
				var[4] = TRACKCHI2[j];
				var[5] = GHOST[j];
			}
			//	 var[4] = cos(Lb_H_OPENING[j]);

			// var[0] = PT;
			// var[1] = IPCHI2;
			// var[2] = VCHI2DOF;
			// var[3] = MINIPCHI2;

			// --- Return the MVA outputs and fill into histograms

			BDT[j] = reader->EvaluateMVA("BDT method");

			//     BDT = reader->EvaluateMVA( "BDT method"          );

		}
		treeout->Fill();
	}
	treeout->Write();
	treeout->Print();
	// Get elapsed time
	sw.Stop();
	std::cout << "--- End of event loop: ";
	sw.Print();

	// --- Write histograms

	target->Close();

	std::cout << "--- Created root file: "<<target->GetName()<<" containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;

	gROOT->ProcessLine(".>");
	if((flag == 1||flag==3) && (version=="v1"))
		gSystem->Exec("cat isoapplication_LL_v1_log.txt");
	else if((flag == 2||flag==4) && (version=="v1"))
		gSystem->Exec("cat isoapplication_DD_v1_log.txt");
	else if((flag == 1||flag==3) && (version=="v2"))
		gSystem->Exec("cat isoapplication_LL_v2_log.txt");
	else if((flag == 2||flag==4) && (version=="v2"))
		gSystem->Exec("cat isoapplication_DD_v2_log.txt");

}
