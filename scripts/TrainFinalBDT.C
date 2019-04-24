/********************************
   Author : Aravindhan V.
 *********************************/
#include "TrainFinalBDT.h"
void TrainFinalBDT(Int_t run, Int_t trackType, const char* isoVersion,
                   Int_t isoConf, Bool_t isoFlag, Bool_t logFlag)
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

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	const char *type = "";

	type = (trackType == 3) ? ("LL") : ("DD");

	if(isoFlag && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainFinalBDT_%s_iso%d_%s.txt",
		                             run,type,isoConf,isoVersion),"w");
	}
	else if(!isoFlag && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainFinalBDT_%s_noIso.txt",
		                             run,type),"w");
	}

	cout<<"*****************************"<<endl;
	cout << "==> Starting TrainFinalBDT: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"*****************************"<<endl;

	TFile *input        = nullptr, *input_iso  = nullptr, *outputFile = nullptr;
	TTree *treeIn       = nullptr, *treeIn_iso = nullptr;
	TString outfileName = "";
	TCut signalCut      = "", bkgCut           = "";
	TCut baseCut        = "", myCutS           = "", myCutB           = "";

	const char *rootFolder = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);

	TMVA::DataLoader *dataLoader = nullptr;
	TMVA::Factory *factory       = nullptr;
	TMVA::Tools::Instance();// This loads the library

	gROOT->ProcessLine("(TMVA::gConfig().GetVariablePlotting())."
	                   "fMaxNumOfAllowedVariablesForScatterPlots = 10;");

	if(isoFlag)
	{
		outfileName = Form("%s/TMVAtraining/iso/"
		                   "TMVA-JpsiLambda_%s_data_iso%d_%s.root",
		                   rootFolder,type,isoConf,isoVersion);
		outputFile  = TFile::Open(outfileName, "RECREATE");
		factory     = new TMVA::Factory(Form("TMVAClassification-"
		                                     "JpsiLambda%s_dataRun%d_iso%d_%s",
		                                     type,run,isoConf,isoVersion),outputFile,
		                                "!V:!Silent:Color:!DrawProgressBar:"
		                                "AnalysisType=Classification");
	}
	else
	{
		outfileName = Form("%s/TMVAtraining/noIso/"
		                   "TMVA-JpsiLambda%s_data_noIso.root",
		                   rootFolder,type);
		outputFile  = TFile::Open( outfileName, "RECREATE");
		factory = new TMVA::Factory( Form("TMVAClassification-"
		                                  "JpsiLambda%s_dataRun%d_noIso",
		                                  type,run), outputFile,
		                             "!V:!Silent:Color:!DrawProgressBar:"
		                             "AnalysisType=Classification" );
	}

	dataLoader = new TMVA::DataLoader("dataset");

	dataLoader->AddVariable("log_dtfchi2         := log10(Lb_ConsLb_chi2)",'F');
	dataLoader->AddVariable("log_lbminipchi2     := log10(Lb_MINIPCHI2)",'F');
	dataLoader->AddVariable("logacos_lbdira      := log10(acos(Lb_DIRA_OWNPV))",'F');
	dataLoader->AddVariable("log_lbfd_ownpv      := log10(Lb_FD_OWNPV)",'F');
	//	dataLoader->AddVariable("log_ltau        := log10(L_TAU)",'F');
	//dataLoader->AddVariable("Lb_DTF_CTAUS_L",'F');

	dataLoader->AddVariable("log_jpsiminipchi2   := log10(Jpsi_MINIPCHI2)",'F');
	dataLoader->AddVariable("log_jpsimass        := log10(Jpsi_M)",'F');
	//	dataLoader->AddVariable("Jpsi_CosTheta",'F');
	//	dataLoader->AddVariable("Jpsi_PT",'F');

	dataLoader->AddVariable("log_lfdchi2         := log10(L_FDCHI2_ORIVX)",'F');
	dataLoader->AddVariable("logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))",'F');
	dataLoader->AddVariable("log_lfd_orivx       := log10(L_FD_ORIVX)",'F');
	dataLoader->AddVariable("logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))",'F');
	dataLoader->AddVariable("L_dm",'F');
	dataLoader->AddVariable("log_lminipchi2      := log10(L_MINIPCHI2)",'F');
	//	dataLoader->AddVariable("L_PT",'F');
	// dataLoader->AddVariable("L_ENDVERTEX_CHI2",'F');
	//	dataLoader->AddVariable("L_CosTheta",'F');

	//	dataLoader->AddVariable("p_PIDp",'F');
	dataLoader->AddVariable("log_pminipchi2      := log10(p_MINIPCHI2)",'F');
	dataLoader->AddVariable("p_ProbNNghost",'F');
	dataLoader->AddVariable("log_p_PT            := log10(p_PT)",'F');
	dataLoader->AddVariable("p_ProbNNp",'F');

	dataLoader->AddVariable("pi_ProbNNghost",'F');
	//	dataLoader->AddVariable("pi_PIDK",'F');
	dataLoader->AddVariable("log_piminipchi2     := log10(pi_MINIPCHI2)",'F');
	dataLoader->AddVariable("log_pi_PT           := log10(pi_PT)",'F');
	dataLoader->AddVariable("pi_ProbNNpi",'F');

	if(isoFlag) dataLoader->AddVariable(Form("BDTkMin_%s",isoVersion),'F');

	if(isoFlag)
	{
		input      = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks.root",
		                              rootFolder,type),"READ");
		treeIn     = (TTree*)input->Get("MyTuple");
		input_iso  = TFile::Open(Form("%s/jpsilambda_%ssig_iso%d_%s.root",
		                              rootFolder,type,isoConf,isoVersion));
		treeIn_iso = (TTree*)input_iso->Get("MyTuple");
	}
	else
	{
		gSystem->Exec(Form("hadd -f %s/jpsilambda_%s_withsw.root "
		                   "%s/jpsilambda_%s_withsw_nonZeroTracks.root"
		                   " %s/jpsilambda_%s_withsw_ZeroTracks.root",
		                   rootFolder,type,rootFolder,type,rootFolder,type));
		input  = TFile::Open(Form("%s/jpsilambda_%s_withsw.root",
		                          rootFolder,type));
		treeIn = (TTree*)input->Get("MyTuple");
	}

	cout << "--- TMVAClassification : Using input file: "
	     << input->GetName() << endl;
	if(isoFlag)
	{
		cout << "--- TMVAClassification : Using isolation file: "
		     << input_iso->GetName() << endl;
	}

	// treeIn = (TTree*)input->Get("MyTuple");
	if(isoFlag) treeIn->AddFriend(treeIn_iso);

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
	//	treeIn->SetBranchStatus("L_ENDVERTEX_CHI2",1);
	treeIn->SetBranchStatus("L_CosTheta",1);

	treeIn->SetBranchStatus("muplus_MINIPCHI2",1);
	treeIn->SetBranchStatus("muminus_MINIPCHI2",1);

	//	treeIn->SetBranchStatus("p_PIDp",1);
	treeIn->SetBranchStatus("p_MINIPCHI2",1);
	treeIn->SetBranchStatus("p_ProbNNghost",1);
	treeIn->SetBranchStatus("p_PT",1);
	treeIn->SetBranchStatus("p_ProbNNp",1);

	//	treeIn->SetBranchStatus("pi_PIDK",1);
	treeIn->SetBranchStatus("pi_MINIPCHI2",1);
	treeIn->SetBranchStatus("pi_ProbNNghost",1);
	treeIn->SetBranchStatus("pi_PT",1);
	treeIn->SetBranchStatus("pi_ProbNNpi",1);

	treeIn->SetBranchStatus("SW",1);
	treeIn->SetBranchStatus("BW",1);

	signalCut = "";
	bkgCut    = "";

	dataLoader->SetInputTrees(treeIn, signalCut, bkgCut);

	dataLoader->SetSignalWeightExpression("SW");
	//	dataLoader->SetBackgroundWeightExpression("BW");

	// Apply additional cuts on the signal and background samples (can be different)
	baseCut = "Lb_ConsLb_chi2 > 0 && Lb_MINIPCHI2 > 0 && Lb_FD_OWNPV > 0 &&"
	          "L_TAU > 0 && Jpsi_MINIPCHI2 > 0 && Jpsi_M > 0 &&"
	          "L_FDCHI2_ORIVX > 0 && L_FD_ORIVX > 0 && L_MINIPCHI2 > 0 &&"
	          "p_MINIPCHI2 > 0 && p_PT > 0 && pi_MINIPCHI2 > 0 &&"
	          "pi_PT > 0";

	myCutS = baseCut;//removed M < 5700 cut as it was not helping
	myCutB = baseCut && "Lb_DTF_M_JpsiLConstr > 5700";

	Int_t nEntries_S = treeIn->GetEntries(myCutS);
	Int_t nEntries_B = treeIn->GetEntries(myCutB);

	Int_t nTrain_S = (Int_t)nEntries_S*0.8;// 80/20 split
	Int_t nTest_S  = nEntries_S - nTrain_S;

	Int_t nTrain_B = (Int_t)nEntries_B*0.8;// 80/20 split
	Int_t nTest_B  = nEntries_B - nTrain_B;

	if(isoFlag)
	{
		myCutS = myCutS && Form("BDTkMin_%s < 1.0",isoVersion);
		myCutB = myCutB && Form("BDTkMin_%s < 1.0",isoVersion);
	}
	dataLoader->PrepareTrainingAndTestTree( myCutS, myCutB,
	                                        Form("nTrain_Signal=%d:nTest_Signal=%d:"
	                                             "nTrain_Background=%d:nTest_Background=%d"
	                                             "SplitMode=Random:NormMode=NumEvents:!V",
	                                             nTrain_S,nTest_S,nTrain_B,nTest_B));

	factory->BookMethod(dataLoader,TMVA::Types::kBDT, "BDTconf1",
	                    "!H:!V:NTrees=500:MinNodeSize=1.0%:MaxDepth=4:"
	                    "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
	                    "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
	                    "nCuts=100" );

	factory->BookMethod(dataLoader,TMVA::Types::kBDT, "BDTconf2",
	                    "!H:!V:NTrees=500:MinNodeSize=0.5%:MaxDepth=4:"
	                    "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
	                    "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
	                    "nCuts=100" );

	/*factory->BookMethod( TMVA::Types::kBDT, "BDTconf6",
	   "!H:!V:NTrees=1000:MinNodeSize=0.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.3:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf2",
	   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.7:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf3",
	   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.7:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf4",
	   "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1" );
	   factory->BookMethod( TMVA::Types::kBDT, "BDTGconf1",
	   "!H:!V:NTrees=200:MinNodeSize=0.625%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=-1:MaxDepth=3" );
	   factory->BookMethod( TMVA::Types::kBDT, "BDTGconf2",
	   "!H:!V:NTrees=400:MinNodeSize=0.625%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=-1:MaxDepth=3" );
	   factory->BookMethod( TMVA::Types::kBDT, "BDTGconf3",
	   "!H:!V:NTrees=800:MinNodeSize=0.625%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=-1:MaxDepth=3" );
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

	// if(!isoFlag)
	// {
	//      gSystem->Exec(Form("rm -f %s/jpsilambda_%s_withsw.root",rootFolder,type));
	// }

	sw.Stop();
	cout << "==> TrainFinalBDT is done! Home streeeeetch!: "; sw.Print();

	// if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);
	// Launch the GUI for the root macros
	//   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}