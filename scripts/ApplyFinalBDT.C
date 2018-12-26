/********************************
   Author : Aravindhan V.
 *********************************/

#include "ApplyFinalBDT.h"

void ApplyFinalBDT(Int_t run, Bool_t isData, Int_t mcType, Int_t trackType,
                   const char* isoVersion, Int_t isoConf, Int_t bdtConf,
                   Int_t flag, Bool_t isoFlag, Bool_t zeroFlag, Bool_t logFlag,
                   Bool_t newFlag)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
   isoVersion = "v1","v2", "v3" or "v4". Refers to isolation isoVersion.
   isoConf = 1 or 2. Refers to config of isolation BDT.
   bdtConf = 1 or 2. refers to finalBDT training configs
   flag = 1 when applying on all data, flag = 2 when applying only on signal training sample
   isoFlag = true if you want to use isolation in the final BDT.
   zeroFlag = true to apply final BDT on zeroTrack data, false to apply on nonZeroTracks data.
 */
{
	cout<<"***********Starting ApplyFinalBDT***********"<<endl;

	TStopwatch sw;
	sw.Start();

	TFile *fileIn = nullptr, *fileIn_iso = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeIn_iso = nullptr, *treeOut = nullptr;

	Int_t entries_init    = 0;
	const char *type      = "", *mynew = "", *logFileName = "";
	const char *logFolder = "", *rootFolder = "";

	Float_t log_dtfChi2       = 0., log_lbMinIpChi2 = 0., logAcos_lbDira_ownPV = 0.;
	Float_t log_lbFd_ownPV    = 0., log_lTau        = 0.;
	Float_t log_jpsiMinIpChi2 = 0., log_jpsiMass    = 0., jpsi_cosTheta        = 0.;
	Float_t jpsi_pt           = 0., pJpsi_lb        = 0., ptJpsi_lb            = 0.;

	Float_t log_lFdChi2_orivX   = 0., logAcos_lDira_orivX = 0., log_lFd_orivX  = 0.;
	Float_t logAcos_lDira_ownPV = 0., l_dm                = 0., log_lMinIpChi2 = 0., l_pt = 0.;
	Float_t l_endVertex_Chi2    = 0., l_cosTheta          = 0.;
	Float_t pLambda_lb          = 0., ptLambda_lb         = 0.;
	Float_t p_PIDp              = 0., log_pMinIpChi2      = 0., p_ghostProb    = 0.;

	Float_t log_p_pt  = 0., p_probNNp       = 0.;
	Float_t pi_PIDK   = 0., log_piMinIpChi2 = 0., pi_ghostProb = 0.;
	Float_t log_pi_pt = 0., pi_probNNpi     = 0.;
	Float_t BDTK      = 0.;
	Float_t chi2array[200] {0.};

	Double_t lbMinIpChi2   = 0., lbDira_ownPV = 0., lbFd_ownPV   = 0.;
	Double_t lTau          = 0., lbPT         = 0., lbP          = 0.;
	Double_t jpsiMinIpChi2 = 0., jpsiMass     = 0., jpsiCosTheta = 0.;
	Double_t jpsiPT        = 0., jpsiP        = 0.;
	Double_t lFdChi2_orivX = 0., lDira_orivX  = 0., lFd_orivX    = 0.;

	Double_t lDira_ownPV = 0., lDm            = 0., lMinIpChi2  = 0., lPT = 0.;
	Double_t lP          = 0., lEndVertexChi2 = 0., lCosTheta   = 0.;
	Double_t pPIDp       = 0., pMinIpChi2     = 0., pGhostProb  = 0.;
	Double_t pPT         = 0., pProbNNp       = 0.;
	Double_t piPIDK      = 0., piMinIpChi2    = 0., piGhostProb = 0.;

	Double_t piPT   = 0., piProbNNpi = 0.;
	Double_t myBDTK = 0., BDT        = 0.;

	TMVA::Tools::Instance();// This loads the library
	TMVA::Reader *reader = nullptr;

	if(newFlag) mynew = "_new";
	type        = (trackType == 3) ? ("LL") : ("DD");
	logFileName = Form("finalBDT%dApplication_%s_iso%d_%s%s.txt",
	                   bdtConf,type,isoConf,isoVersion,mynew);
	cout<<"logFile = "<<logFileName<<endl;

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	cout<<"WD = "<<gSystem->pwd()<<endl;

	//set up logging, input, output
	if(!isData) // MC
	{
		switch(mcType)
		{
		case 1: //JpsiLambda
			logFolder  = Form("logs/mc/JpsiLambda/run%d",run);
			rootFolder = Form("rootFiles/mcFiles/JpsiLambda/run%d",run);
			if(logFlag) gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn      = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks.root",
					                               rootFolder,type));
					treeIn      = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso  = TFile::Open(Form("%s/jpsilambda_%s_iso%d_%s.root",
					                               rootFolder,type,isoConf,isoVersion));
					treeIn_iso  = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut     = new TFile(Form("%s/jpsilambda_%s_FinalBDT%d_iso%d_%s%s.root",
					                             rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                        "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_ZeroTracks.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsilambda_zeroTracks%s_FinalBDT%d_iso%d_%s%s.root",
					                         rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s.root",
				                           rootFolder,type));
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsilambda_%s_FinalBDT%d_noIso%s.root",
				                         rootFolder,type,bdtConf,mynew),"RECREATE");
			}
			break;
		case 2: //JpsiSigma
			logFolder  = Form("logs/mc/JpsiSigma/run%d",run);
			rootFolder = Form("rootFiles/mcFiles/JpsiSigma/run%d",run);
			if(logFlag) gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn     = TFile::Open(Form("%s/jpsisigma_cutoutks_%s_nonZeroTracks.root",
					                              rootFolder,type));
					treeIn     = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso = TFile::Open(Form("%s/jpsisigma_%s_iso%d_%s.root",
					                              rootFolder,type,isoConf,isoVersion));
					treeIn_iso = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut    = new TFile(Form("%s/jpsisigma_%s_FinalBDT%d_iso%d_%s%s.root",
					                            rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                       "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsisigma_cutoutks_%s_ZeroTracks.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsisigma_zeroTracks%s_FinalBDT%d_iso%d_%s%s.root",
					                         rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsisigma_cutoutks_%s.root",
				                           rootFolder,type));
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsisigma_%s_FinalBDT%d_noIso%s.root",
				                         rootFolder,type,bdtConf,mynew),"RECREATE");
			}
			break;
		case 3: //JpsiXi
			logFolder  = Form("logs/mc/JpsiXi/run%d",run);
			rootFolder = Form("rootFiles/mcFiles/JpsiXi/run%d",run);
			if(logFlag) gROOT->ProcessLine(Form(".> %s/%s",logFolder,logFileName));
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn      = TFile::Open(Form("%s/jpsixi_cutoutks_%s_nonZeroTracks.root",
					                               rootFolder,type));
					treeIn      = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso  = TFile::Open(Form("%s/jpsixi_%s_iso%d_%s.root",
					                               rootFolder,type,isoConf,isoVersion));
					treeIn_iso  = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut     = new TFile(Form("%s/jpsixi_%s_FinalBDT%d_iso%d_%s%s.root",
					                             rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                        "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsixi_cutoutks_%s_ZeroTracks.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsixi_zeroTracks%s_FinalBDT%d_iso%d_%s%s.root",
					                         rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsixi_cutoutks_%s.root",
				                           rootFolder,type));
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsixi_%s_FinalBDT%d_noIso%s.root",
				                         rootFolder,type,bdtConf,mynew),"RECREATE");
			}
			break;
		}
	}//end MC block
	else // Data
	{
		logFolder  = Form("logs/data/JpsiLambda/run%d",run);
		rootFolder = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);
		if(logFlag) gSystem->RedirectOutput(Form("%s/%s",logFolder,logFileName),"a");
		if(flag == 1)
		{
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn      = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks.root",
					                               rootFolder,type));
					treeIn      = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso  = TFile::Open(Form("%s/jpsilambda_%s_iso%d_%s.root",
					                               rootFolder,type,isoConf,isoVersion));
					treeIn_iso  = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut     = new TFile(Form("%s/jpsilambda_%s_FinalBDT%d_iso%d_%s%s.root",
					                             rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                        "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_ZeroTracks.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsilambda_zeroTracks%s_FinalBDT%d_iso%d_%s%s.root",
					                         rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s.root",
				                           rootFolder,type));
				//has to be hadded file doesnt exist
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsilambda_%s_FinalBDT%d_noIso%s.root",
				                         rootFolder,type,bdtConf,mynew),
				                    "RECREATE");
			}
		}
		else if(flag == 2)
		{
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn     = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks.root",
					                              rootFolder,type));
					treeIn     = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso = TFile::Open(Form("%s/jpsilambda_%ssig_iso%d_%s.root",
					                              rootFolder,type,isoConf,isoVersion));
					treeIn_iso = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut    = new TFile(Form("%s/jpsilambda_%ssig_FinalBDT%d_iso%d_%s%s.root",
					                            rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                       "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsilambda_%s_withsw_ZeroTracks.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsilambda_zeroTracks%ssig_FinalBDT%d_iso%d_%s%s.root",
					                         rootFolder,type,bdtConf,isoConf,isoVersion,mynew),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsilambda_%s_withsw.root",
				                           rootFolder,type));
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsilambda_%ssig_FinalBDT%d_noIso%s.root",
				                         rootFolder,type,bdtConf,mynew),"RECREATE");
			}
		}
	}//end Data block
	//end setup of input, output, logging
	cout<<"******************************************"<<endl;
	cout<<"Starting ApplyFinalBDT "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;
	if (!fileIn)
	{
		cout << "ERROR: could not open data file" << endl;
		exit(1);
	}
	cout<<"******************************************"<<endl;
	cout<<"Processing Run "<<run<<" "<<type
	    <<((isData) ? (" Data") : (" MC type "))
	    <<mcType<<((flag == 2) ? ("sWeighted") : (""))
	    <<((zeroFlag) ? ("zeroTracks") : (""))<<endl;
	cout<<"******************************************"<<endl;

	treeOut = new TTree("MyTuple","");

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();

	cout<<"Incoming Entries = "<<entries_init;

	reader = new TMVA::Reader("!Color:!Silent");

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type
	//to those given in the weight file(s) used

	reader->AddVariable("log_dtfchi2         := log10(Lb_ConsLb_chi2)", &log_dtfChi2);
	reader->AddVariable("log_lbminipchi2     := log10(Lb_MINIPCHI2)", &log_lbMinIpChi2);
	reader->AddVariable("logacos_lbdira      := log10(acos(Lb_DIRA_OWNPV))", &logAcos_lbDira_ownPV);
	reader->AddVariable("log_lbfd_ownpv      := log10(Lb_FD_OWNPV)", &log_lbFd_ownPV);
	//	reader->AddVariable("log_ltau        := log10(L_TAU)", &log_lTau);
	//reader->AddVariable("Lb_DTF_CTAUS_L", &l_dtfCtauS);

	reader->AddVariable("log_jpsiminipchi2   := log10(Jpsi_MINIPCHI2)", &log_jpsiMinIpChi2);
	reader->AddVariable("log_jpsimass        := log10(Jpsi_M)", &log_jpsiMass);
	//	reader->AddVariable("Jpsi_CosTheta", &jpsi_cosTheta);
	//	reader->AddVariable("Jpsi_PT", &jpsi_pt);

	reader->AddVariable("log_lfdchi2         := log10(L_FDCHI2_ORIVX)", &log_lFdChi2_orivX);
	reader->AddVariable("logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))", &logAcos_lDira_orivX);
	reader->AddVariable("log_lfd_orivx       := log10(L_FD_ORIVX)", &log_lFd_orivX);
	reader->AddVariable("logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))", &logAcos_lDira_ownPV);
	reader->AddVariable("L_dm", &l_dm);
	reader->AddVariable("log_lminipchi2      := log10(L_MINIPCHI2)", &log_lMinIpChi2);
	//	reader->AddVariable("L_PT", &l_pt);
	reader->AddVariable("L_ENDVERTEX_CHI2", &l_endVertex_Chi2);
	//	reader->AddVariable("L_CosTheta", &l_cosTheta);

	reader->AddVariable("p_PIDp", &p_PIDp);
	reader->AddVariable("log_pminipchi2      := log10(p_MINIPCHI2)", &log_pMinIpChi2);
	reader->AddVariable("p_TRACK_GhostProb", &p_ghostProb);
	reader->AddVariable("log_p_PT            := log10(p_PT)", &log_p_pt);
	reader->AddVariable("p_ProbNNp", &p_probNNp);

	//	reader->AddVariable("pi_TRACK_GhostProb", &pi_ghostProb);
	//	reader->AddVariable("pi_PIDK", &pi_PIDK);
	reader->AddVariable("log_piminipchi2     := log10(pi_MINIPCHI2)", &log_piMinIpChi2);
	reader->AddVariable("log_pi_PT           := log10(pi_PT)", &log_pi_pt);
	reader->AddVariable("pi_ProbNNpi", &pi_probNNpi);

	if(newFlag)
	{
		reader->AddVariable("pJpsi_lb    := (Jpsi_P/Lb_P)",&pJpsi_lb);
		reader->AddVariable("pLambda_lb  := (L_P/Lb_P)",&pLambda_lb);
		reader->AddVariable("ptJpsi_lb   := (Jpsi_PT/Lb_PT)",&ptJpsi_lb);
		reader->AddVariable("ptLambda_lb := (L_PT/Lb_PT)",&ptLambda_lb);
	}

	if(isoFlag && !zeroFlag)
	{
		reader->AddVariable(Form("BDTkMin_%s",isoVersion), &BDTK);
	}

	//***********Load weights file**********************************
	TString dir    = "dataset/weights/";
	TString prefix;
	if(isoFlag && !zeroFlag)
	{
		prefix = Form("TMVAClassificationLite-JpsiLambda%s_dataRun%d_iso%d_%s%s",
		              type,run,isoConf,isoVersion,mynew);
	}
	else if(!isoFlag || (isoFlag && zeroFlag))
	{
		prefix = Form("TMVAClassificationLite-JpsiLambda%s_dataRun%d_noIso%s",
		              type,run,mynew);
	}
	TString methodName = TString("BDT method");
	TString weightfile = dir + prefix + TString("_") +
	                     Form("BDTconf%d",bdtConf) + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile);
	//*************************************************************

	cout << "--- TMVAClassificationApp    : Using input file: "
	     << fileIn->GetName() << endl;

	treeIn->SetBranchAddress("Lb_ConsLb_chi2", chi2array);
	treeIn->SetBranchAddress("Lb_MINIPCHI2", &lbMinIpChi2);
	treeIn->SetBranchAddress("Lb_DIRA_OWNPV", &lbDira_ownPV );
	treeIn->SetBranchAddress("Lb_FD_OWNPV", &lbFd_ownPV);
	treeIn->SetBranchAddress("L_TAU", &lTau);
	//treeIn->SetBranchAddress("Lb_DTF_CTAUS_L", &l_dtfCtauS);
	treeIn->SetBranchAddress("Lb_P", &lbP);
	treeIn->SetBranchAddress("Lb_PT", &lbPT);

	treeIn->SetBranchAddress("Jpsi_MINIPCHI2", &jpsiMinIpChi2);
	treeIn->SetBranchAddress("Jpsi_M", &jpsiMass);
	treeIn->SetBranchAddress("Jpsi_CosTheta", &jpsiCosTheta);
	treeIn->SetBranchAddress("Jpsi_PT", &jpsiPT);
	treeIn->SetBranchAddress("Jpsi_P", &jpsiP);

	treeIn->SetBranchAddress("L_FDCHI2_ORIVX", &lFdChi2_orivX);
	treeIn->SetBranchAddress("L_DIRA_ORIVX", &lDira_orivX);
	treeIn->SetBranchAddress("L_FD_ORIVX", &lFd_orivX);
	treeIn->SetBranchAddress("L_DIRA_OWNPV", &lDira_ownPV);
	treeIn->SetBranchAddress("L_dm", &lDm);
	treeIn->SetBranchAddress("L_MINIPCHI2", &lMinIpChi2);
	treeIn->SetBranchAddress("L_PT", &lPT);
	treeIn->SetBranchAddress("L_P", &lP);
	treeIn->SetBranchAddress("L_ENDVERTEX_CHI2", &lEndVertexChi2);
	treeIn->SetBranchAddress("L_CosTheta", &lCosTheta);

	treeIn->SetBranchAddress("p_PIDp", &pPIDp);
	treeIn->SetBranchAddress("p_MINIPCHI2", &pMinIpChi2);
	treeIn->SetBranchAddress("p_TRACK_GhostProb", &pGhostProb);
	treeIn->SetBranchAddress("p_PT", &pPT);
	treeIn->SetBranchAddress("p_ProbNNp", &pProbNNp);

	treeIn->SetBranchAddress("pi_TRACK_GhostProb", &piGhostProb);
	treeIn->SetBranchAddress("pi_PIDK", &piPIDK);
	treeIn->SetBranchAddress("pi_MINIPCHI2", &piMinIpChi2);
	treeIn->SetBranchAddress("pi_PT", &piPT);
	treeIn->SetBranchAddress("pi_ProbNNpi", &piProbNNpi);

	if(isoFlag && !zeroFlag)
	{
		treeIn_iso->SetBranchAddress(Form("BDTkMin_%s",isoVersion), &myBDTK);
	}

	treeOut->Branch( Form("BDT%d",bdtConf), &BDT, Form("BDT%d/D",bdtConf));

	cout << "--- Processing: " << entries_init
	     << " events" << endl;

	for (Long64_t ievt = 0; ievt < entries_init; ievt++) {

		if (ievt%50000 == 0) cout << "--- ... Processing event: " << ievt << endl;

		treeIn->GetEntry(ievt);
		if(isoFlag && !zeroFlag) treeIn_iso->GetEntry(ievt);

		log_dtfChi2          = log10(chi2array[0]);
		log_lbMinIpChi2      = log10(lbMinIpChi2);
		logAcos_lbDira_ownPV = log10(acos(lbDira_ownPV));
		log_lbFd_ownPV       = log10(lbFd_ownPV);
		//	log_lTau         = log10(lTau);

		log_jpsiMinIpChi2    = log10(jpsiMinIpChi2);
		log_jpsiMass         = log10(jpsiMass);
		//	jpsi_cosTheta    = jpsiCosTheta;
		//	jpsi_pt          = jpsiPT;

		log_lFdChi2_orivX    = log10(lFdChi2_orivX);
		logAcos_lDira_orivX  = log10(acos(lDira_orivX));
		log_lFd_orivX        = log10(lFd_orivX);
		logAcos_lDira_ownPV  = log10(acos(lDira_ownPV));
		l_dm                 = lDm;
		log_lMinIpChi2       = log10(lMinIpChi2);
		//	l_pt             = lPT;
		l_endVertex_Chi2     = lEndVertexChi2;
		//	l_cosTheta       = lCosTheta;

		pPIDp                = p_PIDp;
		log_pMinIpChi2       = log10(pMinIpChi2);
		p_ghostProb          = pGhostProb;
		log_p_pt             = log10(pPT);
		p_probNNp            = pProbNNp;

		// pi_ghostProb      = piGhostProb;
		//	pi_PIDK          = piPIDK;
		log_piMinIpChi2      = log10(piMinIpChi2);
		log_pi_pt            = log10(piPT);
		pi_probNNpi          = piProbNNpi;

		if(isoFlag && !zeroFlag) BDTK = myBDTK;

		if(newFlag) {
			pJpsi_lb    = jpsiP/lbP;
			ptJpsi_lb   = jpsiPT/lbPT;
			pLambda_lb  = lP/lbP;
			ptLambda_lb = lPT/lbPT;
		}
		BDT = reader->EvaluateMVA("BDT method");

		treeOut->Fill();
	}
	// Get elapsed time
	sw.Stop();
	cout << "--- End of event loop: ";

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	cout << "--- Created root file: "<<fileOut->GetName()
	     <<" containing the MVA output histograms" << endl;

	delete reader;

	cout << "==> ApplyFinalBDT is done! Go optimize and fit!: "; sw.Print();
	//if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);
}
