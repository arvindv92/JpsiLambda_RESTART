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

   >isoVersion = "v0","v1". Different versions correspond
   to different training variables.

   v0: IPCHI2, VCHI2DOF, log_MINIPCHI2
   v1: IPCHI2, VCHI2DOF, log_MINIPCHI2, log_PT
 */
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	const char *type = "";
	type        = (trackType == 3) ? ("LL") : ("DD");

	if(logFlag) //set up logging
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainIsolation_%s_%s.txt",
		                             run,type,isoVersion),"w");
	}
	cout<<"*****************************"<<endl;
	cout<<"==> Starting TrainIsolation: "<<isoVersion<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"*****************************"<<endl;

	TFile *input_sig = nullptr, *input_bkg = nullptr, *outputFile = nullptr;
	TTree *sigTree   = nullptr, *bkgTree   = nullptr;

	TString outFileName = "", fname_sig    = "", fname_bkg = "";
	TCut myCut          = "";

	TMVA::DataLoader *dataLoader = nullptr;
	TMVA::Factory *factory       = nullptr;
	TMVA::Tools::Instance(); // This loads the library

	outFileName = Form("rootFiles/dataFiles/JpsiLambda/run%d/"
	                   "TMVAtraining/iso/TMVA300-isok%s_data_%s.root",
	                   run,type,isoVersion);
	outputFile  = TFile::Open(outFileName, "RECREATE");

	factory     = new TMVA::Factory(Form("TMVAClassification300-isok%s_dataRun%d_%s",
	                                     type,run,isoVersion), outputFile,
	                                "!V:!Silent:Color:!DrawProgressBar:"
	                                "AnalysisType=Classification" );

	dataLoader  = new TMVA::DataLoader("dataset");

	dataLoader->AddVariable("IPCHI2",'F');
	dataLoader->AddVariable("VCHI2DOF",'F');
	dataLoader->AddVariable("log_MINIPCHI2 := log10(MINIPCHI2)",'F');

	myCut = "MINIPCHI2 > 0 && VCHI2DOF > 0";

	if(strncmp(isoVersion,"v1",2)==0)
	{
		dataLoader->AddVariable("log_PT := log10(PT)",'F');

		myCut = myCut && "PT > 0 ";
	}
	// else if(strncmp(isoVersion,"v2",2)==0)
	// {
	//      dataLoader->AddVariable("log_FD     := log10(FD)",'F');
	//      dataLoader->AddVariable("log_FDCHI2 := log10(FDCHI2)",'F');
	//
	//      myCut = myCut && "FD > 0 && FDCHI2 > 0";
	// }
	// else if(strncmp(isoVersion,"v3",2)==0)
	// {
	//      dataLoader->AddVariable("log_PT     := log10(PT)",'F');
	//      dataLoader->AddVariable("log_FD     := log10(FD)",'F');
	//      dataLoader->AddVariable("log_FDCHI2 := log10(FDCHI2)",'F');
	//
	//      myCut = myCut && "PT > 0 && FDCHI2 > 0 && FD > 0";
	// }

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
	dataLoader->AddSignalTree(sigTree,1.0);
	dataLoader->AddBackgroundTree(bkgTree,1.0);

	dataLoader->SetSignalWeightExpression("SW");
	dataLoader->SetBackgroundWeightExpression("SW");

	Int_t nEntries_S = sigTree->GetEntries(myCut);
	Int_t nEntries_B = bkgTree->GetEntries(myCut);

	Int_t nTrain_S = (Int_t)nEntries_S*0.8;// 80/20 split
	Int_t nTest_S  = nEntries_S - nTrain_S;

	Int_t nTrain_B = (Int_t)nEntries_B*0.8;// 80/20 split
	Int_t nTest_B  = nEntries_B - nTrain_B;

	dataLoader->PrepareTrainingAndTestTree( myCut, myCut,
	                                        Form("nTrain_Signal=%d:nTest_Signal=%d:"
	                                             "nTrain_Background=%d:nTest_Background=%d"
	                                             "SplitMode=Random:NormMode=NumEvents:!V",
	                                             nTrain_S,nTest_S,nTrain_B,nTest_B));

	factory->BookMethod(dataLoader, TMVA::Types::kBDT, "isoConf1_300",
	                    "!H:!V:NTrees=300:MinNodeSize=1.25%:MaxDepth=3:"
	                    "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
	                    "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
	                    "nCuts=200" );

	factory->BookMethod(dataLoader, TMVA::Types::kBDT, "isoConf2_300",
	                    "!H:!V:NTrees=300:MinNodeSize=0.75%:MaxDepth=3:"
	                    "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
	                    "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
	                    "nCuts=200" );

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

	// if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);

}
