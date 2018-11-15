// @(#)root/tmva $Id$
/**********************************************************************************
* Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
* Package   : TMVA                                                               *
* Root Macro: TMVAClassification                                                 *
*                                                                                *
* This macro provides examples for the training and testing of the               *
* TMVA classifiers.                                                              *
*                                                                                *
* As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
* and linearly correlated input variables.                                       *
*                                                                                *
* The methods to be used can be switched on and off by means of booleans, or     *
* via the prompt command, for example:                                           *
*                                                                                *
*    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
*                                                                                *
* (note that the backslashes are mandatory)                                      *
* If no method given, a default set of classifiers is used.                      *
*                                                                                *
* The output file "TMVA.root" can be analysed with the use of dedicated          *
* macros (simply say: root -l <macro.C>), which can be conveniently              *
* invoked through a GUI that will appear at the end of the run of this macro.    *
* Launch the GUI via the command:                                                *
*                                                                                *
*    root -l ./TMVAGui.C                                                         *
*                                                                                *
**********************************************************************************/

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

void TMVAClassification_JpsiLambda(Int_t run = 1,Int_t trackType = 3, TString version = "v1")
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   version = "v1","v2" or "v3"
 */
{
	Bool_t logFlag = false;
	const char *logFileName = "", *type = "";
	TString outfileName = "";
	TMVA::DataLoader *dataloader = NULL;
	TMVA::Factory *factory = NULL;
	logFileName = (trackType == 3) ? (TString::Format("finalBDTTraining_LL_%s_log.txt",version.Data())) : (TString::Format("finalBDTTraining_DD_%s_log.txt",version.Data()));

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
	//---------------------------------------------------------------
	// This loads the library
	TMVA::Tools::Instance();

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassification" << std::endl;

	gROOT->ProcessLine("(TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 10;");
	// --------------------------------------------------------------------------------------------------

	// --- Here the preparation phase begins

	// Create a ROOT output file where TMVA will store ntuples, histograms, etc.


	outfileName = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/TMVA-JpsiLambda%s_data_%s.root",run,type,version.Data());

	TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

	factory = new TMVA::Factory( TString::Format("TMVAClassification-JpsiLambda%s_data_%s",type,version.Data()), outputFile,
	                             "!V:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );

	dataloader = new TMVA::DataLoader("dataset");

	dataloader->AddVariable( "log_dtfchi2 := log10(Lb_ConsLb_chi2)", 'F' );
	dataloader->AddVariable( "log_lbminipchi2 := log10(Lb_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "logacos_lbdira := log10(acos(Lb_DIRA_OWNPV))", 'F' );
	dataloader->AddVariable( "log_lbfd_ownpv := log10(Lb_FD_OWNPV)", 'F' );
	dataloader->AddVariable( "log_dtfctaul := log10(Lb_DTF_CTAU_L)", 'F' );
	dataloader->AddVariable( "Lb_DTF_CTAUS_L", 'F' );

	dataloader->AddVariable( "log_jpsiminipchi2 := log10(Jpsi_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "log_jpsimass := log10(Jpsi_M)", 'F' );
	dataloader->AddVariable( "Jpsi_CosTheta", 'F' );
	dataloader->AddVariable( "Jpsi_PT", 'F' );

	dataloader->AddVariable( "log_lfdchi2 := log10(L_FDCHI2_ORIVX)", 'F' );
	dataloader->AddVariable( "logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))", 'F' );
	dataloader->AddVariable( "log_lfd_orivx := log10(L_FD_ORIVX)", 'F' );
	dataloader->AddVariable( "logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))", 'F' );
	dataloader->AddVariable( "L_dm", 'F' );
	dataloader->AddVariable( "log_lminipchi2 := log10(L_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "L_PT", 'F' );
	dataloader->AddVariable( "L_ENDVERTEX_CHI2", 'F' );
	dataloader->AddVariable( "L_CosTheta", 'F' );

	dataloader->AddVariable( "p_PIDp", 'F' );
	dataloader->AddVariable( "log_pminipchi2 := log10(p_MINIPCHI2)", 'F' );
	dataloader->AddVariable( "p_TRACK_GhostProb", 'F' );
	dataloader->AddVariable( "log_p_PT := log10(p_PT)", 'F' );
	dataloader->AddVariable( "p_ProbNNp", 'F' );

	dataloader->AddVariable( "pi_TRACK_GhostProb", 'F' );
	dataloader->AddVariable( "pi_PIDK", 'F');
	dataloader->AddVariable( "log_piminipchi2 := log10(pi_MINIPCHI2)", 'F');
	dataloader->AddVariable( "log_pi_PT := log10(pi_PT)", 'F' );
	dataloader->AddVariable( "pi_ProbNNpi", 'F' );

	dataloader->AddVariable( TString::Format("BDTkMin_%s",version.Data()), 'F' );

	TFile *input(0);//Signal and Background training files
	TTree *treeIn(0);
	TString fname, fname1;
	TCut signalCut, bkgCut;

	input = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withiso_%s.root",run,type,version.Data()));

	std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
	//std::cout << "--- TMVAClassification       : Using signal input file: " << input->GetName() << std::endl;

	treeIn = (TTree*)input->Get("MyTuple");

	treeIn->SetBranchStatus("*",0);
	treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	treeIn->SetBranchStatus("Lb_ConsLb_chi2",1);
	treeIn->SetBranchStatus("Lb_MINIPCHI2",1);
	treeIn->SetBranchStatus("Lb_DIRA_OWNPV",1);
	treeIn->SetBranchStatus("Lb_FD_OWNPV",1);
	treeIn->SetBranchStatus("Lb_DTF_CTAUS_L",1);
	treeIn->SetBranchStatus("Lb_DTF_CTAU_L",1);

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

	treeIn->SetBranchStatus(TString::Format("BDTkMin_%s",version.Data()),1);

	signalCut = "Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000";//Maybe add Lbmass < 5700 here?
	bkgCut    = "Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000 && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 7000";

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
	TCut mycuts = "";
	TCut mycutb = "";

	dataloader->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

	factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTconf1",
	                     "!H:!V:NTrees=850:MinNodeSize=1.25%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );
	factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTconf5",
	                     "!H:!V:NTrees=850:MinNodeSize=0.625%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

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

	std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
	std::cout << "==> TMVAClassification is done!" << std::endl;

	delete factory;

	if(logFlag) gROOT->ProcessLine(".>");
	// Launch the GUI for the root macros
	//   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
