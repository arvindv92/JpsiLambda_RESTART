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
#include "TStopwatch.h"
#include "TROOT.h"
//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVAGui.h"
#endif

using namespace std;

void TrainFinalBDT(Int_t run = 1,Int_t trackType = 3, TString isoVersion = "v1", Int_t isoConf = 1, Bool_t isoFlag = true, Bool_t logFlag = false)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   isoVersion = "v0","v1","v2" or "v3"
   isoConf = 1 or 5
   isoFlag = true if you want to use isolation in the final BDT.
 */
{
	TStopwatch sw;
	sw.Start();

	TString outfileName = "";
	const char *logFileName = "", *type = "";
	TFile *input = nullptr, *outputFile = nullptr;//Signal and Background training files
	TTree *treeIn = nullptr;
	TCut signalCut = "", bkgCut = "";
	TMVA::DataLoader *dataloader = nullptr;
	TMVA::Factory *factory = nullptr;
	TMVA::Tools::Instance();// This loads the library

	logFileName = (trackType == 3) ? (TString::Format("finalBDTTraining_LL_iso%d_%s_log.txt",isoConf,isoVersion.Data())) : (TString::Format("finalBDTTraining_DD_iso%d_%s_log.txt",isoConf,isoVersion.Data()));

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	if(logFlag) //set up logging
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

	cout<<"*****************************"<<endl;
	cout << "==> Starting TrainFinalBDT: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"*****************************"<<endl;


	gROOT->ProcessLine("(TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 10;");

	if(isoFlag)
	{
		outfileName = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/TMVAtraining/iso/TMVALite-JpsiLambda_%s_data_iso%d_%s.root",run,type,isoConf,isoVersion.Data());
		outputFile = TFile::Open( outfileName, "RECREATE" );
		factory = new TMVA::Factory( TString::Format("TMVAClassificationLite-JpsiLambda%s_dataRun%d_iso%d_%s",type,run,isoConf,isoVersion.Data()), outputFile,
		                             "!V:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );
	}
	else
	{
		outfileName = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/TMVAtraining/noIso/TMVALite-JpsiLambda%s_data_noIso.root",run,type);
		outputFile = TFile::Open( outfileName, "RECREATE");
		factory = new TMVA::Factory( TString::Format("TMVAClassificationLite-JpsiLambda%s_dataRun%d_noIso",type,run), outputFile,
		                             "!V:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );
	}

	dataloader = new TMVA::DataLoader("dataset");

	dataloader->AddVariable( "log_dtfchi2 := log10(Lb_ConsLb_chi2)", 'F' );
	dataloader->AddVariable( "log_lbminipchi2 := log10(Lb_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "logacos_lbdira := log10(acos(Lb_DIRA_OWNPV))", 'F' );
	dataloader->AddVariable( "log_lbfd_ownpv := log10(Lb_FD_OWNPV)", 'F' );
//	dataloader->AddVariable( "log_ltau := log10(L_TAU)", 'F' );
	//dataloader->AddVariable( "Lb_DTF_CTAUS_L", 'F' );

	dataloader->AddVariable( "log_jpsiminipchi2 := log10(Jpsi_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "log_jpsimass := log10(Jpsi_M)", 'F' );
//	dataloader->AddVariable( "Jpsi_CosTheta", 'F' );
//	dataloader->AddVariable( "Jpsi_PT", 'F' );

	dataloader->AddVariable( "log_lfdchi2 := log10(L_FDCHI2_ORIVX)", 'F' );
	dataloader->AddVariable( "logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))", 'F' );
	dataloader->AddVariable( "log_lfd_orivx := log10(L_FD_ORIVX)", 'F' );
	dataloader->AddVariable( "logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))", 'F' );
	dataloader->AddVariable( "L_dm", 'F' );
	dataloader->AddVariable( "log_lminipchi2 := log10(L_MINIPCHI2)", 'F' );
//	dataloader->AddVariable( "L_PT", 'F' );
	dataloader->AddVariable( "L_ENDVERTEX_CHI2", 'F' );
//	dataloader->AddVariable( "L_CosTheta", 'F' );

	dataloader->AddVariable( "p_PIDp", 'F' );
	dataloader->AddVariable( "log_pminipchi2 := log10(p_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "p_TRACK_GhostProb", 'F' );
	dataloader->AddVariable( "log_p_PT := log10(p_PT)", 'F' );
	dataloader->AddVariable( "p_ProbNNp", 'F' );

//	dataloader->AddVariable( "pi_TRACK_GhostProb", 'F' );
//	dataloader->AddVariable( "pi_PIDK", 'F');
	dataloader->AddVariable( "log_piminipchi2 := log10(pi_MINIPCHI2)", 'F');
	dataloader->AddVariable( "log_pi_PT := log10(pi_PT)", 'F' );
	dataloader->AddVariable( "pi_ProbNNpi", 'F' );

	if(isoFlag) dataloader->AddVariable( TString::Format("BDTkMin_%s",isoVersion.Data()), 'F' );

	if(isoFlag)
	{
		input = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_withiso%d_%s.root",run,type,isoConf,isoVersion.Data()));
	}
	else
	{
		input = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw.root",run,type));
	}

	cout << "--- TMVAClassification       : Using input file: " << input->GetName() << endl;
	//cout << "--- TMVAClassification       : Using signal input file: " << input->GetName() << endl;

	treeIn = (TTree*)input->Get("MyTuple");

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	treeIn->SetBranchStatus("Lb_ConsLb_chi2",1);
	treeIn->SetBranchStatus("Lb_MINIPCHI2",1);
	treeIn->SetBranchStatus("Lb_DIRA_OWNPV",1);
	treeIn->SetBranchStatus("Lb_FD_OWNPV",1);
//	treeIn->SetBranchStatus("Lb_DTF_CTAUS_L",1);
	treeIn->SetBranchStatus("L_TAU",1);

	treeIn->SetBranchStatus("Jpsi_MINIPCHI2",1);
	treeIn->SetBranchStatus("Jpsi_M",1);
	treeIn->SetBranchStatus("Jpsi_CosTheta",1);
	treeIn->SetBranchStatus("Jpsi_PT",1);

	treeIn->SetBranchStatus("L_FDCHI2_ORIVX",1);
	treeIn->SetBranchStatus("L_DIRA_ORIVX",1);
	treeIn->SetBranchStatus("L_FD_ORIVX",1);
	treeIn->SetBranchStatus("L_DIRA_OWNPV",1);
	treeIn->SetBranchStatus("L_dm",1);
	treeIn->SetBranchStatus("L_MINIPCHI2",1);
	treeIn->SetBranchStatus("L_PT",1);
	treeIn->SetBranchStatus("L_ENDVERTEX_CHI2",1);
	treeIn->SetBranchStatus("L_CosTheta",1);

	treeIn->SetBranchStatus("muplus_MINIPCHI2",1);
	treeIn->SetBranchStatus("muminus_MINIPCHI2",1);

	treeIn->SetBranchStatus("p_PIDp",1);
	treeIn->SetBranchStatus("p_MINIPCHI2",1);
	treeIn->SetBranchStatus("p_TRACK_GhostProb",1);
	treeIn->SetBranchStatus("p_PT",1);
	treeIn->SetBranchStatus("p_ProbNNp",1);

	treeIn->SetBranchStatus("pi_PIDK",1);
	treeIn->SetBranchStatus("pi_MINIPCHI2",1);
	treeIn->SetBranchStatus("pi_TRACK_GhostProb",1);
	treeIn->SetBranchStatus("pi_PT",1);
	treeIn->SetBranchStatus("pi_ProbNNpi",1);

	if(isoFlag) treeIn->SetBranchStatus(TString::Format("BDTkMin_%s",isoVersion.Data()),1);

	treeIn->SetBranchStatus("SW",1);
	treeIn->SetBranchStatus("BW",1);

	signalCut = "";
	bkgCut    = "";

	dataloader->SetInputTrees(treeIn, signalCut, bkgCut);

	// global event weights per tree (see below for setting event-wise weights)
	// Double_t signalWeight     = 1.0;
	// Double_t backgroundWeight = 1.0;

	// You can add an arbitrary number of signal or background trees
	// dataloader->AddSignalTree    ( signal,     signalWeight     );
	// dataloader->AddBackgroundTree( background, backgroundWeight );

	dataloader->SetSignalWeightExpression("SW");
	dataloader->SetBackgroundWeightExpression("BW");

	// Apply additional cuts on the signal and background samples (can be different)
	TCut mycuts = "Lb_ConsLb_chi2 > 0 && Lb_MINIPCHI2 > 0 && Lb_FD_OWNPV > 0 && L_TAU > 0 && Jpsi_MINIPCHI2 > 0 && Jpsi_M > 0 && L_FDCHI2_ORIVX > 0 && L_FD_ORIVX > 0 && L_MINIPCHI2 > 0 && p_MINIPCHI2 > 0 && p_PT > 0 && pi_MINIPCHI2 > 0 && pi_PT > 0 && Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000 && Lb_DTF_M_JpsiLConstr < 5700";//Maybe add Lbmass < 5700 here? Maybe dtfchi2 cut?
	TCut mycutb = "Lb_ConsLb_chi2 > 0 && Lb_MINIPCHI2 > 0 && Lb_FD_OWNPV > 0 && L_TAU > 0 && Jpsi_MINIPCHI2 > 0 && Jpsi_M > 0 && L_FDCHI2_ORIVX > 0 && L_FD_ORIVX > 0 && L_MINIPCHI2 > 0 && p_MINIPCHI2 > 0 && p_PT > 0 && pi_MINIPCHI2 > 0 && pi_PT > 0 && Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000 && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 7000";
	cout<<"3"<<endl;
	dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

	cout<<"3.5"<<endl;
	factory->BookMethod(dataloader,TMVA::Types::kBDT, "BDTconf1",
	                    "!H:!V:NTrees=300:MinNodeSize=1.0%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
	cout<<"4"<<endl;
	factory->BookMethod(dataloader,TMVA::Types::kBDT, "BDTconf2",
	                    "!H:!V:NTrees=300:MinNodeSize=0.5%:MaxDepth=5:BoostType=AdaBoost:AdaBoostBeta=0.2:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	/*factory->BookMethod( TMVA::Types::kBDT, "BDTconf6",
	   "!H:!V:NTrees=1000:MinNodeSize=0.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf2",
	                     "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.7:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf3",
	                     "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.7:SeparationType=GiniIndex:nCuts=20" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf4",
	                     "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
	   factory->BookMethod( TMVA::Types::kBDT, "BDTGconf1",
	                     "!H:!V:NTrees=200:MinNodeSize=0.625%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
	   factory->BookMethod( TMVA::Types::kBDT, "BDTGconf2",
	                     "!H:!V:NTrees=400:MinNodeSize=0.625%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
	   factory->BookMethod( TMVA::Types::kBDT, "BDTGconf3",
	                     "!H:!V:NTrees=800:MinNodeSize=0.625%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3" );
	   Support Vector Machine
	   factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );*/
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

	cout << "==> Wrote root file: " << outputFile->GetName() << endl;
	cout << "==> TMVAClassification is done!" << endl;

	delete factory;

	sw.Stop();
	cout << "==> TrainFinalBDT is done! Home streeeeetch!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
	// Launch the GUI for the root macros
	//   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
