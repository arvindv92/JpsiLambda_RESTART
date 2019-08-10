/********************************
   Author : Aravindhan V.
 *********************************/
#include "TrainFinalBDT.h"
void TrainFinalBDT(Int_t run, Int_t trackType, const char* isoVersion,
                   Int_t isoConf, Bool_t isoFlag, Int_t bdtConf, Bool_t logFlag,
                   Bool_t simFlag)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   isoVersion = "v0","v1","v2" or "v3"
   isoConf = 1 or 5
   isoFlag = true if you want to use isolation in the final BDT.
   simFlag = false to train BDT on sWeighted data (default) for signal. true to train on MC for signal instead. both options use data for background.
 */
//NB MODIFIED: Removing pi_ProbNNpi from variables because MC doesn't model it well
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	const char *type = "";

	type = (trackType == 3) ? ("LL") : ("DD");

	if(!simFlag) //training on data
	{
		if(isoFlag && logFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainFinalBDT%d_%s_iso%d_%s.txt",
			                             run,bdtConf,type,isoConf,isoVersion),"w");
		}
		else if(!isoFlag && logFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainFinalBDT%d_%s_noIso.txt",
			                             run,bdtConf,type),"w");
		}
	}
	else //training on MC
	{
		if(isoFlag && logFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainFinalBDT%d_MC_%s_iso%d_%s.txt",
			                             run,bdtConf,type,isoConf,isoVersion),"w");
		}
		else if(!isoFlag && logFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/TrainFinalBDT%d_MC_%s_noIso.txt",
			                             run,bdtConf,type),"w");
		}
	}
	cout<<"*****************************"<<endl;
	cout<<"*****************************"<<endl;
	cout<<"*****************************"<<endl;
	cout<<"*****************************"<<endl;
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

	const char *rootFolder   = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);
	const char *MCrootFolder = Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d",run);

	TMVA::DataLoader *dataLoader = nullptr;
	// TMVA::Factory *factory       = nullptr;
	TMVA::Tools::Instance();// This loads the library

	gROOT->ProcessLine("(TMVA::gConfig().GetVariablePlotting())."
	                   "fMaxNumOfAllowedVariablesForScatterPlots = 10;");

	// if(!simFlag)
	// {
	//      if(isoFlag)
	//      {
	//              outfileName = Form("%s/TMVAtraining/iso/"
	//                                 "data_iso%d_%s_BDT%d.root",
	//                                 rootFolder,isoConf,isoVersion,bdtConf);
	//              outputFile  = TFile::Open(outfileName, "RECREATE");
	//              factory     = new TMVA::Factory(Form("dataRun%d_iso%d_%s_BDT%d",
	//                                                   run,isoConf,isoVersion,bdtConf),outputFile,
	//                                              "!V:!Silent:Color:!DrawProgressBar:"
	//                                              "AnalysisType=Classification");
	//      }
	//      else
	//      {
	//              outfileName = Form("%s/TMVAtraining/noIso/"
	//                                 "data_noIso_BDT%d.root",
	//                                 rootFolder,bdtConf);
	//              outputFile  = TFile::Open( outfileName, "RECREATE");
	//              factory = new TMVA::Factory( Form("dataRun%d_noIso_BDT%d",
	//                                                run,bdtConf), outputFile,
	//                                           "!V:!Silent:Color:!DrawProgressBar:"
	//                                           "AnalysisType=Classification" );
	//      }
	// }
	// else
	// {
	//      if(isoFlag)
	//      {
	//              outfileName = Form("%s/TMVAtraining/iso/"
	//                                 "MC_iso%d_%s_BDT%d.root",
	//                                 rootFolder,isoConf,isoVersion,bdtConf);
	//              outputFile  = TFile::Open(outfileName, "RECREATE");
	//              factory     = new TMVA::Factory(Form("MCRun%d_iso%d_%s_BDT%d",
	//                                                   run,isoConf,isoVersion,bdtConf),outputFile,
	//                                              "!V:!Silent:Color:!DrawProgressBar:"
	//                                              "AnalysisType=Classification");
	//      }
	//      else
	//      {
	//              outfileName = Form("%s/TMVAtraining/noIso/"
	//                                 "MC_noIso_BDT%d.root",
	//                                 rootFolder,bdtConf);
	//              outputFile  = TFile::Open( outfileName, "RECREATE");
	//              factory = new TMVA::Factory( Form("MCRun%d_noIso_BDT%d",
	//                                                run,bdtConf), outputFile,
	//                                           "!V:!Silent:Color:!DrawProgressBar:"
	//                                           "AnalysisType=Classification" );
	//      }
	// }
	dataLoader = new TMVA::DataLoader("dataset");

	dataLoader->AddSpectator("eventNumber := eventNumber % 4096", 'I');
	dataLoader->AddSpectator("runNumber", 'I');
	dataLoader->AddSpectator("ntracks", 'I');

	dataLoader->AddVariable("log_dtfchi2     := log10(Lb_ConsLb_chi2)",'F');
	dataLoader->AddVariable("log_lbminipchi2 := log10(Lb_MINIPCHI2)",'F');
	dataLoader->AddVariable("logacos_lbdira  := log10(acos(Lb_DIRA_OWNPV))",'F');
	dataLoader->AddVariable("log_lbfd_ownpv  := log10(Lb_FD_OWNPV)",'F');
	//	dataLoader->AddVariable("log_ltau        := log10(L_TAU)",'F');
	//dataLoader->AddVariable("Lb_DTF_CTAUS_L",'F');

	dataLoader->AddVariable("log_jpsiminipchi2 := log10(Jpsi_MINIPCHI2)",'F');
	dataLoader->AddVariable("log_jpsimass      := log10(Jpsi_M)",'F');
	//	dataLoader->AddVariable("Jpsi_CosTheta",'F');
	//	dataLoader->AddVariable("Jpsi_PT",'F');

	dataLoader->AddVariable("log_lfdchi2         := log10(L_FDCHI2_ORIVX)",'F');
	dataLoader->AddVariable("logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))",'F');
	dataLoader->AddVariable("log_lfd_orivx       := log10(L_FD_ORIVX)",'F');
	dataLoader->AddVariable("logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))",'F');
	dataLoader->AddVariable("L_dm",'F');
	dataLoader->AddVariable("log_lminipchi2 := log10(L_MINIPCHI2)",'F');
	dataLoader->AddVariable("L_ENDVERTEX_CHI2",'F');
	//	dataLoader->AddVariable("L_PT",'F');

	//	dataLoader->AddVariable("L_CosTheta",'F');

	dataLoader->AddVariable("log_pminipchi2      := log10(p_MINIPCHI2)",'F');
	dataLoader->AddVariable("p_ProbNNghost",'F');
	dataLoader->AddVariable("log_p_PT            := log10(p_PT)",'F');
	dataLoader->AddVariable("p_ProbNNp",'F');

	dataLoader->AddVariable("pi_ProbNNghost",'F');
	dataLoader->AddVariable("log_piminipchi2     := log10(pi_MINIPCHI2)",'F');
	dataLoader->AddVariable("log_pi_PT           := log10(pi_PT)",'F');

	if(isoFlag) dataLoader->AddVariable(Form("BDTkMin_%s",isoVersion),'F');

	TFile *input_sig = nullptr, *input_bkg = nullptr;
	TTree *treeIn_sig = nullptr, *treeIn_bkg = nullptr;
	TFile *input_iso_sig = nullptr, *input_iso_bkg = nullptr;
	TTree *treeIn_iso_sig = nullptr, *treeIn_iso_bkg = nullptr;

	if(!simFlag)//training on data
	{
		if(isoFlag)
		{
			input      = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks.root",
			                              rootFolder,type),"READ");
			treeIn     = (TTree*)input->Get("MyTuple");
			input_iso  = TFile::Open(Form("%s/jpsilambda_%ssig_iso%d_%s.root",
			                              rootFolder,type,isoConf,isoVersion));//This is the isolation applied to sWeighted data
			treeIn_iso = (TTree*)input_iso->Get("MyTuple");
		}
		else
		{
			// gSystem->Exec(Form("hadd -f %s/jpsilambda_%s_withsw.root "
			//                    "%s/jpsilambda_%s_withsw_nonZeroTracks.root"
			//                    " %s/jpsilambda_%s_withsw_ZeroTracks.root",
			//                    rootFolder,type,rootFolder,type,rootFolder,type));
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
	}
	else //training on simulation
	{
		if(isoFlag)
		{
			input_sig      = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks.root",
			                                  MCrootFolder,type),"READ");
			treeIn_sig     = (TTree*)input_sig->Get("MyTuple");
			input_iso_sig  = TFile::Open(Form("%s/jpsilambda_%s_iso%d_%s.root",
			                                  MCrootFolder,type,isoConf,isoVersion));//This is the isolation applied to sWeighted data
			treeIn_iso_sig = (TTree*)input_iso_sig->Get("MyTuple");

			input_bkg      = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks.root",
			                                  rootFolder,type),"READ");
			treeIn_bkg     = (TTree*)input_bkg->Get("MyTuple");
			input_iso_bkg  = TFile::Open(Form("%s/jpsilambda_%ssig_iso%d_%s.root",
			                                  rootFolder,type,isoConf,isoVersion));//This is the isolation applied to sWeighted data
			treeIn_iso_bkg = (TTree*)input_iso_bkg->Get("MyTuple");
		}
		else
		{
			input_sig      = TFile::Open(Form("%s/jpsilambda_cutoutks_%s.root",
			                                  MCrootFolder,type),"READ");
			treeIn_sig     = (TTree*)input_sig->Get("MyTuple");

			input_bkg      = TFile::Open(Form("%s/jpsilambda_%s_withsw.root",
			                                  rootFolder,type),"READ");
			treeIn_bkg     = (TTree*)input_bkg->Get("MyTuple");
		}
	}
	// treeIn = (TTree*)input->Get("MyTuple");
	if(!simFlag)// train on data
	{
		if(isoFlag) treeIn->AddFriend(treeIn_iso);

		treeIn->SetBranchStatus("*",0);
		treeIn->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		treeIn->SetBranchStatus("Lb_ConsLb_chi2",1);
		treeIn->SetBranchStatus("Lb_MINIPCHI2",1);
		treeIn->SetBranchStatus("Lb_DIRA_OWNPV",1);
		treeIn->SetBranchStatus("Lb_FD_OWNPV",1);
		//	treeIn->SetBranchStatus("Lb_DTF_CTAUS_L",1);
		// treeIn->SetBranchStatus("L_TAU",1);

		treeIn->SetBranchStatus("Jpsi_MINIPCHI2",1);
		treeIn->SetBranchStatus("Jpsi_M",1);
		// treeIn->SetBranchStatus("Jpsi_CosTheta",1);
		// treeIn->SetBranchStatus("Jpsi_PT",1);

		treeIn->SetBranchStatus("L_FDCHI2_ORIVX",1);
		treeIn->SetBranchStatus("L_DIRA_ORIVX",1);
		treeIn->SetBranchStatus("L_FD_ORIVX",1);
		treeIn->SetBranchStatus("L_DIRA_OWNPV",1);
		treeIn->SetBranchStatus("L_dm",1);
		treeIn->SetBranchStatus("L_MINIPCHI2",1);
		// treeIn->SetBranchStatus("L_PT",1);
		treeIn->SetBranchStatus("L_ENDVERTEX_CHI2",1);
		// treeIn->SetBranchStatus("L_CosTheta",1);

		//	treeIn->SetBranchStatus("p_PIDp",1);
		treeIn->SetBranchStatus("p_MINIPCHI2",1);
		treeIn->SetBranchStatus("p_ProbNNghost",1);
		treeIn->SetBranchStatus("p_PT",1);
		treeIn->SetBranchStatus("p_ProbNNp",1);

		//	treeIn->SetBranchStatus("pi_PIDK",1);
		treeIn->SetBranchStatus("pi_MINIPCHI2",1);
		treeIn->SetBranchStatus("pi_ProbNNghost",1);
		treeIn->SetBranchStatus("pi_PT",1);
		// treeIn->SetBranchStatus("pi_ProbNNpi",1);

		treeIn->SetBranchStatus("SW",1);
		treeIn->SetBranchStatus("BW",1);

		treeIn->SetBranchStatus("eventNumber",1);
		treeIn->SetBranchStatus("runNumber",1);
		treeIn->SetBranchStatus("ntracks",1);
	}
	else// train on MC
	{
		if(isoFlag)
		{
			treeIn_sig->AddFriend(treeIn_iso_sig);
			treeIn_bkg->AddFriend(treeIn_iso_bkg);
		}
		treeIn_bkg->SetBranchStatus("*",0);
		treeIn_bkg->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		treeIn_bkg->SetBranchStatus("Lb_ConsLb_chi2",1);
		treeIn_bkg->SetBranchStatus("Lb_MINIPCHI2",1);
		treeIn_bkg->SetBranchStatus("Lb_DIRA_OWNPV",1);
		treeIn_bkg->SetBranchStatus("Lb_FD_OWNPV",1);
		//	treeIn_bkg->SetBranchStatus("Lb_DTF_CTAUS_L",1);
		// treeIn_bkg->SetBranchStatus("L_TAU",1);

		treeIn_bkg->SetBranchStatus("Jpsi_MINIPCHI2",1);
		treeIn_bkg->SetBranchStatus("Jpsi_M",1);
		// treeIn_bkg->SetBranchStatus("Jpsi_CosTheta",1);
		// treeIn_bkg->SetBranchStatus("Jpsi_PT",1);

		treeIn_bkg->SetBranchStatus("L_FDCHI2_ORIVX",1);
		treeIn_bkg->SetBranchStatus("L_DIRA_ORIVX",1);
		treeIn_bkg->SetBranchStatus("L_FD_ORIVX",1);
		treeIn_bkg->SetBranchStatus("L_DIRA_OWNPV",1);
		treeIn_bkg->SetBranchStatus("L_dm",1);
		treeIn_bkg->SetBranchStatus("L_MINIPCHI2",1);
		// treeIn_bkg->SetBranchStatus("L_PT",1);
		treeIn_bkg->SetBranchStatus("L_ENDVERTEX_CHI2",1);
		// treeIn_bkg->SetBranchStatus("L_CosTheta",1);

		//	treeIn_bkg->SetBranchStatus("p_PIDp",1);
		treeIn_bkg->SetBranchStatus("p_MINIPCHI2",1);
		treeIn_bkg->SetBranchStatus("p_ProbNNghost",1);
		treeIn_bkg->SetBranchStatus("p_PT",1);
		treeIn_bkg->SetBranchStatus("p_ProbNNp",1);

		//	treeIn_bkg->SetBranchStatus("pi_PIDK",1);
		treeIn_bkg->SetBranchStatus("pi_MINIPCHI2",1);
		treeIn_bkg->SetBranchStatus("pi_ProbNNghost",1);
		treeIn_bkg->SetBranchStatus("pi_PT",1);
		// treeIn_bkg->SetBranchStatus("pi_ProbNNpi",1);

		treeIn_bkg->SetBranchStatus("eventNumber",1);
		treeIn_bkg->SetBranchStatus("runNumber",1);
		treeIn_bkg->SetBranchStatus("ntracks",1);

		treeIn_sig->SetBranchStatus("*",0);
		treeIn_sig->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		treeIn_sig->SetBranchStatus("Lb_ConsLb_chi2",1);
		treeIn_sig->SetBranchStatus("Lb_MINIPCHI2",1);
		treeIn_sig->SetBranchStatus("Lb_DIRA_OWNPV",1);
		treeIn_sig->SetBranchStatus("Lb_FD_OWNPV",1);
		//	treeIn_sig->SetBranchStatus("Lb_DTF_CTAUS_L",1);
		// treeIn_sig->SetBranchStatus("L_TAU",1);

		treeIn_sig->SetBranchStatus("Jpsi_MINIPCHI2",1);
		treeIn_sig->SetBranchStatus("Jpsi_M",1);
		// treeIn_sig->SetBranchStatus("Jpsi_CosTheta",1);
		// treeIn_sig->SetBranchStatus("Jpsi_PT",1);

		treeIn_sig->SetBranchStatus("L_FDCHI2_ORIVX",1);
		treeIn_sig->SetBranchStatus("L_DIRA_ORIVX",1);
		treeIn_sig->SetBranchStatus("L_FD_ORIVX",1);
		treeIn_sig->SetBranchStatus("L_DIRA_OWNPV",1);
		treeIn_sig->SetBranchStatus("L_dm",1);
		treeIn_sig->SetBranchStatus("L_MINIPCHI2",1);
		// treeIn_sig->SetBranchStatus("L_PT",1);
		treeIn_sig->SetBranchStatus("L_ENDVERTEX_CHI2",1);
		// treeIn_sig->SetBranchStatus("L_CosTheta",1);

		//	treeIn_sig->SetBranchStatus("p_PIDp",1);
		treeIn_sig->SetBranchStatus("p_MINIPCHI2",1);
		treeIn_sig->SetBranchStatus("p_ProbNNghost",1);
		treeIn_sig->SetBranchStatus("p_PT",1);
		treeIn_sig->SetBranchStatus("p_ProbNNp",1);

		//	treeIn_sig->SetBranchStatus("pi_PIDK",1);
		treeIn_sig->SetBranchStatus("pi_MINIPCHI2",1);
		treeIn_sig->SetBranchStatus("pi_ProbNNghost",1);
		treeIn_sig->SetBranchStatus("pi_PT",1);
		// treeIn_sig->SetBranchStatus("pi_ProbNNpi",1);

		treeIn_sig->SetBranchStatus("GB_WT",1);
		treeIn_sig->SetBranchStatus("Lb_BKGCAT",1);

		treeIn_sig->SetBranchStatus("eventNumber",1);
		treeIn_sig->SetBranchStatus("runNumber",1);
		treeIn_sig->SetBranchStatus("ntracks",1);
	}
	signalCut = "";
	bkgCut    = "";

	if(!simFlag)
	{
		dataLoader->SetInputTrees(treeIn, signalCut, bkgCut);
		dataLoader->SetSignalWeightExpression("SW");
	}
	else
	{
		dataLoader->SetSignalTree(treeIn_sig);
		dataLoader->SetBackgroundTree(treeIn_bkg);
		dataLoader->SetSignalWeightExpression("GB_WT");
	}

	//	dataLoader->SetBackgroundWeightExpression("BW");
	Int_t nEntries_S = 0, nEntries_B = 0;

	// Apply additional cuts on the signal and background samples (can be different)
	baseCut = "Lb_ConsLb_chi2 > 0 && Lb_MINIPCHI2 > 0 && Lb_FD_OWNPV > 0 &&"
	          "Jpsi_MINIPCHI2 > 0 && Jpsi_M > 0 &&"
	          "L_FDCHI2_ORIVX > 0 && L_FD_ORIVX > 0 && L_MINIPCHI2 > 0 &&"
	          "p_MINIPCHI2 > 0 && p_PT > 0 && pi_MINIPCHI2 > 0 &&"
	          "pi_PT > 0";

	myCutS = baseCut;//removed M < 5700 cut as it was not helping
	myCutB = baseCut && "Lb_DTF_M_JpsiLConstr > 5700 && "
	         "!(Lb_DTF_M_JpsiLConstr>5770 && Lb_DTF_M_JpsiLConstr<5810)";//cutting out portion where Xib0->J/psi Lambda hides

	if(isoFlag)
	{
		myCutS = myCutS && Form("BDTkMin_%s < 1.0",isoVersion);
		myCutB = myCutB && Form("BDTkMin_%s < 1.0",isoVersion);
	}

	if(simFlag)
	{
		myCutS = myCutS && "(Lb_BKGCAT==0||Lb_BKGCAT==50)";
		nEntries_S = treeIn_sig->GetEntries(myCutS);
		nEntries_B = treeIn_bkg->GetEntries(myCutB);
	}
	else
	{
		nEntries_S = treeIn->GetEntries(myCutS);
		nEntries_B = treeIn->GetEntries(myCutB);
	}

	Int_t nTrain_S = (Int_t)nEntries_S*0.8;// 80/20 split
	Int_t nTest_S  = nEntries_S - nTrain_S;

	Int_t nTrain_B = (Int_t)nEntries_B*0.8;// 80/20 split
	Int_t nTest_B  = nEntries_B - nTrain_B;

	// dataLoader->PrepareTrainingAndTestTree( myCutS, myCutB,
	//                                         Form("nTrain_Signal=%d:nTest_Signal=%d:"
	//                                              "nTrain_Background=%d:nTest_Background=%d"
	//                                              "SplitMode=Random:NormMode=NumEvents:!V",
	//                                              nTrain_S,nTest_S,nTrain_B,nTest_B));

	dataLoader->PrepareTrainingAndTestTree( myCutS, myCutB,
	                                        "nTest_Signal=1"
	                                        ":nTest_Background=1"
	                                        ":SplitMode=Random"
	                                        ":NormMode=NumEvents"
	                                        ":!V");

	// TString splitExpr = " ( 19*int([eventNumber]) + 29*int([runNumber]) + 37*int([ntracks]) )% int([NumFolds]) ";
	TString splitExpr = "int(fabs([eventNumber]))%int([NumFolds])";

	TString cvOptions = Form("!V"
	                         ":!Silent"
	                         ":AnalysisType=Classification"
	                         ":SplitType=Deterministic"
	                         ":NumFolds=5"
	                         ":SplitExpr=%s",
	                         splitExpr.Data());
	TMVA::CrossValidation *cv = nullptr;

	if(!simFlag)
	{
		if(isoFlag)
		{
			outfileName = Form("%s/TMVAtraining/iso/"
			                   "CVdata_iso%d_%s_BDT%d.root",
			                   rootFolder,isoConf,isoVersion,bdtConf);
			outputFile  = TFile::Open(outfileName, "RECREATE");
			cv          = new TMVA::CrossValidation{Form("CVdataRun%d_iso%d_%s_BDT%d",
				                                     run,isoConf,isoVersion,bdtConf),dataLoader, outputFile, cvOptions};
		}
		else
		{
			outfileName = Form("%s/TMVAtraining/noIso/"
			                   "CVdata_noIso_BDT%d.root",
			                   rootFolder,bdtConf);
			outputFile  = TFile::Open( outfileName, "RECREATE");
			cv          = new TMVA::CrossValidation{Form("CVdataRun%d_noIso_BDT%d",
				                                     run,bdtConf),dataLoader, outputFile, cvOptions};
		}
	}
	else
	{
		if(isoFlag)
		{
			outfileName = Form("%s/TMVAtraining/iso/"
			                   "CVMC_iso%d_%s_BDT%d.root",
			                   rootFolder,isoConf,isoVersion,bdtConf);
			outputFile  = TFile::Open(outfileName, "RECREATE");
			cv     = new TMVA::CrossValidation{Form("CVMCRun%d_iso%d_%s_BDT%d",
				                                run,isoConf,isoVersion,bdtConf),dataLoader, outputFile, cvOptions};
		}
		else
		{
			outfileName = Form("%s/TMVAtraining/noIso/"
			                   "CVMC_noIso_BDT%d.root",
			                   rootFolder,bdtConf);
			outputFile  = TFile::Open( outfileName, "RECREATE");
			cv = new TMVA::CrossValidation{ Form("CVMCRun%d_noIso_BDT%d",
				                             run,bdtConf), dataLoader, outputFile, cvOptions};
		}
	}

	if(bdtConf == 1)
	{
		cv->BookMethod(TMVA::Types::kBDT, "BDTconf1",
		               "!H:!V:NTrees=500:MinNodeSize=1.0%:MaxDepth=4:"
		               "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
		               "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
		               "nCuts=100" );
		// factory->BookMethod(dataLoader,TMVA::Types::kBDT, "BDTconf1",
		//                "!H:!V:NTrees=500:MinNodeSize=1.0%:MaxDepth=4:"
		//                "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
		//                "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
		//                "nCuts=100" );

		// factory->BookMethod(dataLoader,TMVA::Types::kBDT, "BDTconf1_REAL",
		//                     "!H:!V:NTrees=500:MinNodeSize=1.0%:MaxDepth=4:"
		//                     "BoostType=RealAdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
		//                     "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
		//                     "nCuts=100" );
	}
	else if(bdtConf == 2)

	{
		cv->BookMethod(TMVA::Types::kBDT, "BDTconf2",
		               "!H:!V:NTrees=500:MinNodeSize=0.5%:MaxDepth=4:"
		               "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
		               "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
		               "nCuts=100" );
		// factory->BookMethod(dataLoader,TMVA::Types::kBDT, "BDTconf2",
		//                "!H:!V:NTrees=500:MinNodeSize=0.5%:MaxDepth=4:"
		//                "BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
		//                "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
		//                "nCuts=100" );


		// factory->BookMethod(dataLoader,TMVA::Types::kBDT, "BDTconf2_REAL",
		//                     "!H:!V:NTrees=500:MinNodeSize=0.5%:MaxDepth=4:"
		//                     "BoostType=RealAdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:"
		//                     "BaggedSampleFraction=0.5:SeparationType=GiniIndex:"
		//                     "nCuts=100" );
	}

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

	// cv->TrainAllMethods();
	//
	// cv->TestAllMethods();
	//
	// cv->EvaluateAllMethods();
	cv->Evaluate();
	// --------------------------------------------------------------

	// Save the output
	outputFile->Close();

	cout << "==> Wrote root file: " << outputFile->GetName() << endl;
	cout << "==> TMVAClassification is done!" << endl;

	delete cv;

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
