/********************************
   Author : Aravindhan V.
   The purpose of this script is to train a BDT to isolate background coming
   from b -> Jpsi + charged tracks
   Signal Training: sWeighted Lb -> J/psi L data coming from DoSWeight
   Background Training: sWeighted B+ -> J/psi K+ data.
   Output: Weights file
   TODO: Also train against B0 -> J/psi K+ pi-
 *********************************/
 #include "TrainIsolation.h"

void TrainIsolation(Int_t run, Int_t trackType,
                    const char* isoVersion, Bool_t logFlag)
/*
   >run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC.
   Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data

   >trackType = 3 for LL, 5 for DD.

   >isoVersion = "v0","v1","v2","v3". Different versions correspond
   to different training variables.

   v0: IPCHI2, VCHI2DOF, log_MINIPCHI2
   v1: IPCHI2, VCHI2DOF, log_MINIPCHI2, log_PT
   v2: IPCHI2, VCHI2DOF, log_MINIPCHI2, log_FD, log_FDCHI2
   v3: IPCHI2, VCHI2DOF, log_MINIPCHI2, log_FD, log_FDCHI2, log_PT
 */
{
	cout<<"***********Starting TrainIsolation***********"<<endl;

	TStopwatch sw;
	sw.Start();

	TFile *input_sig = nullptr, *input_bkg = nullptr, *outputFile = nullptr;
	TTree *sigTree   = nullptr, *bkgTree   = nullptr;

	TString outFileName = "", fname_sig    = "", fname_bkg = "";
	TCut myCut          = "";
	const char *type    = "", *logFileName = "";

	TMVA::Factory *factory       = nullptr;
	TMVA::DataLoader *dataloader = nullptr;
	TMVA::Tools::Instance(); // This loads the library

	type        = (trackType == 3) ? ("LL") : ("DD");
	logFileName = Form("isolationTraining_%s_%s_log.txt",
	                   type,isoVersion);
	cout<<"logFile = "<<logFileName<<endl;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	if(logFlag) //set up logging
	{
		gROOT->ProcessLine(Form(".> logs/data/JpsiLambda/run%d/%s",
		                        run,logFileName));
	}
	gSystem->Exec("date");
	cout<<endl;
	cout<<"******************************************"<<endl;
	cout<<"Processing Run "<<run<<" "<<type<<endl;
	cout<<"******************************************"<<endl;

	cout<<"*****************************"<<endl;
	cout<<"==> Starting TrainIsolation: "<<isoVersion<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"*****************************"<<endl;

	outFileName = Form("rootFiles/dataFiles/JpsiLambda/run%d/"
	                   "TMVAtraining/iso/TMVA300-isok%s_data_%s.root",
	                   run,type,isoVersion);
	outputFile  = TFile::Open(outFileName, "RECREATE");

	factory     = new TMVA::Factory(Form("TMVAClassification300-isok%s_dataRun%d_%s",
	                                     type,run,isoVersion), outputFile,
	                                "!V:!Silent:Color:!DrawProgressBar:"
	                                "AnalysisType=Classification" );

	dataloader  = new TMVA::DataLoader("dataset");

	dataloader->AddVariable("IPCHI2",'F');
	dataloader->AddVariable("VCHI2DOF",'F');
	dataloader->AddVariable("log_MINIPCHI2 := log10(MINIPCHI2)",'F');

	myCut = "MINIPCHI2 > 0 && VCHI2DOF > 0";

	if(strncmp(isoVersion,"v1",2)==0)
	{
		dataloader->AddVariable("log_PT := log10(PT)",'F');

		myCut = myCut && "PT > 0 ";
	}
	else if(strncmp(isoVersion,"v2",2)==0)
	{
		dataloader->AddVariable("log_FD     := log10(FD)",'F');
		dataloader->AddVariable("log_FDCHI2 := log10(FDCHI2)",'F');

		myCut = myCut && "FD > 0 && FDCHI2 > 0";
	}
	else if(strncmp(isoVersion,"v3",2)==0)
	{
		dataloader->AddVariable("log_PT     := log10(PT)",'F');
		dataloader->AddVariable("log_FD     := log10(FD)",'F');
		dataloader->AddVariable("log_FDCHI2 := log10(FDCHI2)",'F');

		myCut = myCut && "PT > 0 && FDCHI2 > 0 && FD > 0";
	}

	fname_sig = Form("rootFiles/dataFiles/JpsiLambda/run%d/"
	                 "jpsilambda_%s_forIsoTraining.root",run,type);
	input_sig = TFile::Open(fname_sig);
	if (!input_sig)
	{
		cout << "ERROR: could not open signal training file" << endl;
		exit(1);
	}

	fname_bkg = Form("/data1/avenkate/B+toJpsiK+/RealData_AllKaons/"
	                 "total_run%d/splot/maketuples/jpsik_forTMVAisoTraining.root",
	                 run);
	input_bkg = TFile::Open(fname_bkg);
	if (!input_bkg)
	{
		cout << "ERROR: could not open background training file" << endl;
		exit(1);
	}

	// TString fname2 = "./jpsikpi.root";//Background 2 Training Sample
	// input2 = TFile::Open( fname2 );

	cout << "--- TMVAClassification : Using background input file: "
	     << input_bkg->GetName() << endl;
	cout << "--- TMVAClassification : Using signal input file: "
	     << input_sig->GetName() << endl;

	sigTree = (TTree*)input_sig->Get("MyTuple");
	bkgTree = (TTree*)input_bkg->Get("MyTuple");

	// global event weights per tree (see below for setting event-wise weights)
	dataloader->AddSignalTree(sigTree,1.0);
	dataloader->AddBackgroundTree(bkgTree,1.0);

	dataloader->SetSignalWeightExpression("SW");
	dataloader->SetBackgroundWeightExpression("SW");

	dataloader->PrepareTrainingAndTestTree(myCut, myCut,
	                                       "nTrain_Signal=0:nTrain_Background=0:"
	                                       "SplitMode=Random:NormMode=NumEvents:"
	                                       "!V");

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "isoConf1_300",
	                    "!H:!V:NTrees=300:MinNodeSize=1.25%:MaxDepth=3:"
	                    "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
	                    "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
	                    "nCuts=-1" );

	factory->BookMethod(dataloader, TMVA::Types::kBDT, "isoConf2_300",
	                    "!H:!V:NTrees=300:MinNodeSize=0.75%:MaxDepth=3:"
	                    "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
	                    "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
	                    "nCuts=-1" );

	/*	factory->BookMethod( TMVA::Types::kBDT, "BDTconf2",
	   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.7:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf3",
	   "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.7:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTconf4",
	   "!H:!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
	   "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=-1" );

	   factory->BookMethod( TMVA::Types::kBDT, "BDTG",
	   "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=-1:MaxDepth=2" );
	 */
	// ---- STILL EXPERIMENTAL and only implemented for BDT's !
	// factory->OptimizeAllMethods("SigEffAt001","Scan");
	// factory->OptimizeAllMethods("ROCIntegral","FitGA");

	// -----------------------------------------------------------------------

	// ---- train, test, and evaluate the MVAs

	factory->TrainAllMethods();

	factory->TestAllMethods();

	factory->EvaluateAllMethods();

	// Save the output
	outputFile->Close();

	cout << "==> Wrote root file: " << outputFile->GetName() << endl;
	cout << "==> TMVAClassification is done!" << endl;

	delete factory;

	sw.Stop();
	cout << "==> End of TrainIsolation! Cheers!: "; sw.Print();

	if(logFlag) gROOT->ProcessLine(".>");
}
