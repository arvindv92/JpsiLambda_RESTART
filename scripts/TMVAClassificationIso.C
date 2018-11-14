/********************************
   Author : Aravindhan V.
 *********************************/
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
#endif

void TMVAClassificationIso(Int_t run = 1,Int_t trackType = 3, TString version = "v1")
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   version = "v0","v1","v2" or "v3"
 */
{
	Bool_t logFlag = false;
	TString outfileName = "", fname_sig = "", fname_bkg = "";
	TCut mycuts = "", mycutb = "";
	const char *logFileName = "", *type = "";
	TFile *outputFile = NULL, *input_sig = NULL, *input_bkg = NULL;
	TTree *sigTree = NULL,*bkgTree = NULL;
	TMVA::Factory *factory = NULL;
	TMVA::DataLoader *dataloader = NULL;
	TMVA::Tools::Instance(); // This loads the library

	logFileName = (trackType == 3) ? (TString::Format("isolationTraining_LL_%s_log.txt",version.Data())) : (TString::Format("isolationTraining_DD_%s_log.txt",version.Data()));

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.

	if(logFlag)
	{
		gROOT->ProcessLine((TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName)).Data());
	}
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

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassification" << std::endl;

	outfileName = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/TMVA-isok%s_data_%s.root",run,type,version.Data());
	outputFile  = TFile::Open(outfileName, "RECREATE");
	factory     = new TMVA::Factory( TString::Format("TMVAClassification-isok%s_data_%s",type,version.Data()), outputFile,
	                                 "!V:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );

	dataloader = new TMVA::DataLoader("dataset");

	if(version == "v0" || version == "v1" || version == "v2" || version == "v3")
	{
		dataloader->AddVariable("IPCHI2",'F');
		dataloader->AddVariable("VCHI2DOF",'F');
		dataloader->AddVariable("MINIPCHI2",'F');
	}
	if(version == "v1")
	{
		dataloader->AddVariable("PT",'F');
	}
	if(version == "v2")
	{
		dataloader->AddVariable("FD",'F');
	}
	if(version == "v3")
	{
		dataloader->AddVariable("FDCHI2",'F');
	}

	fname_sig = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_forIsoTraining.root",run,type);
	input_sig = TFile::Open(fname_sig);

	fname_bkg = TString::Format("/data1/avenkate/B+toJpsiK+/RealData_AllKaons/total_run%d/splot/maketuples/jpsik_forTMVAisoTraining.root",run);//Background 1 Training Sample
	input_bkg = TFile::Open( fname_bkg );

	// TString fname2 = "./jpsikpi.root";//Background 2 Training Sample
	// input2 = TFile::Open( fname2 );

	std::cout << "--- TMVAClassification       : Using background input file: " << input_bkg->GetName() << std::endl;
	std::cout << "--- TMVAClassification       : Using signal input file: " << input_sig->GetName() << std::endl;

	// --- Register the training and test trees

	sigTree = (TTree*)input_sig->Get("MyTuple");
	bkgTree = (TTree*)input_bkg->Get("MyTuple");

	// global event weights per tree (see below for setting event-wise weights)
	Double_t signalWeight     = 1.0;
	Double_t backgroundWeight = 1.0;

	// You can add an arbitrary number of signal or background trees
	dataloader->AddSignalTree    (sigTree,signalWeight);
	dataloader->AddBackgroundTree(bkgTree,backgroundWeight);

	dataloader->SetSignalWeightExpression("SW");
	dataloader->SetBackgroundWeightExpression("SW");
	//dataloader->SetBackgroundWeightExpression("BW");

	// Apply additional cuts on the signal and background samples (can be different)
	mycuts = "PT > 0 && MINIPCHI2 > 0 && VCHI2DOF > 0 && PT < 10000 && MINIPCHI2 < 50";//CHANGED FROM MINIPCHI2 < 10000
	mycutb = "PT > 0 && MINIPCHI2 > 0 && VCHI2DOF > 0 && PT < 10000 && MINIPCHI2 < 50";//CHANGED FROM MINIPCHI2 < 10000
	if(version == "v2")
	{
		mycuts = "PT > 0 && MINIPCHI2 > 0 && VCHI2DOF > 0 && PT < 10000 && MINIPCHI2 < 100 && FD < 10000";
		mycutb = "PT > 0 && MINIPCHI2 > 0 && VCHI2DOF > 0 && PT < 10000 && MINIPCHI2 < 100 && FD < 10000";
	}
	else if(version == "v3")
	{
		mycuts = "PT > 0 && MINIPCHI2 > 0 && VCHI2DOF > 0 && PT < 10000 && MINIPCHI2 < 100 && FDCHI2 < 1000";
		mycutb = "PT > 0 && MINIPCHI2 > 0 && VCHI2DOF > 0 && PT < 10000 && MINIPCHI2 < 100 && FDCHI2 < 1000";
	}

	dataloader->PrepareTrainingAndTestTree(mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTconf1",
	                    "!H:!V:NTrees=850:MinNodeSize=1.25%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDTconf5",
	                    "!H:!V:NTrees=850:MinNodeSize=0.75%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	/*	factory->BookMethod( TMVA::Types::kBDT, "BDTconf2",
	   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.7:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf3",
	   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.7:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf4",
	   "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
	   "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTG",
	   "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
	 */
	// ---- STILL EXPERIMENTAL and only implemented for BDT's !
	// factory->OptimizeAllMethods("SigEffAt001","Scan");
	// factory->OptimizeAllMethods("ROCIntegral","FitGA");

	// --------------------------------------------------------------------------------------------------

	// ---- Now you can tell the factory to train, test, and evaluate the MVAs

	factory->TrainAllMethods();

	factory->TestAllMethods();

	factory->EvaluateAllMethods();

	// --------------------------------------------------------------

	// Save the output
	outputFile->Close();

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;

	delete factory;

	//      Launch the GUI for the root macros
	//   if (!gROOT->IsBatch()) TMVAGui( outfileName );

	if(logFlag) gROOT->ProcessLine(".>");
}
