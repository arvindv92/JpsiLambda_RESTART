/********************************
   Author : Aravindhan V.
 *********************************/

#include "ApplyFinalBDT.h"

void ApplyFinalBDT(Int_t run, Bool_t isData, Int_t mcType, Int_t trackType,
                   const char* isoVersion, Int_t isoConf, Int_t bdtConf,
                   Int_t flag, Bool_t isoFlag, Bool_t zeroFlag, Bool_t logFlag,
                   Bool_t simFlag)
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
//NB MODIFIED: Removing pi_ProbNNpi from variables because MC doesn't model it well
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	const char *folder = "", *part = "", *type = "";

	type        = (trackType == 3) ? ("LL") : ("DD");
	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Bu_JpsiX";
		part = "bu_jpsix";
		break;
	}
	case 5:
	{
		folder = "Bd_JpsiX";
		part = "bd_jpsix";
		break;
	}
	case 6:
	{
		folder = "Lst1405";
		part = "lst1405";
		break;
	}
	case 7:
	{
		folder = "Lst1520";
		part = "lst1520";
		break;
	}
	case 8:
	{
		folder = "Lst1600";
		part = "lst1600";
		break;
	}
	case 9:
	{
		folder = "chiC1";
		part = "chic1";
		break;
	}
	case 10:
	{
		folder = "JpsiKs";
		part = "jpsiks";
		break;
	}
	case 11:
	{
		folder = "Xib0";
		part = "xib0";
	}
	}
	//Set up logging
	if(isData && isoFlag && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/ApplyFinalBDT%d_%s_iso%d_%s_noPID.txt",
		                             run,bdtConf,type,isoConf,isoVersion),"a");
	}
	else if(!isData && isoFlag && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/ApplyFinalBDT%d_%s_iso%d_%s_noPID.txt",
		                             folder,run,bdtConf,type,isoConf,isoVersion),"w");
	}
	else if(isData && !isoFlag && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/ApplyFinalBDT%d_%s_noIso_noPID.txt",
		                             run,bdtConf,type),"a");
	}
	else if(!isData && !isoFlag && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/ApplyFinalBDT%d_%s_noIso_noPID.txt",
		                             folder,run,bdtConf,type),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"Starting ApplyFinalBDT "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr, *fileIn_iso = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeIn_iso = nullptr, *treeOut = nullptr;

	Int_t entries_init     = 0;
	const char *rootFolder = "";

	Float_t log_dtfChi2       = 0., log_lbMinIpChi2 = 0., logAcos_lbDira_ownPV = 0.;
	Float_t log_lbFd_ownPV    = 0.;
	Float_t log_jpsiMinIpChi2 = 0., log_jpsiMass    = 0.;

	Float_t log_lFdChi2_orivX   = 0., logAcos_lDira_orivX = 0., log_lFd_orivX  = 0.;
	Float_t logAcos_lDira_ownPV = 0., l_dm                = 0., log_lMinIpChi2 = 0.;
	Float_t l_EndVtx_Chi2 = 0.;
	Float_t log_pMinIpChi2      = 0., p_ghostProbNN    = 0.;

	Float_t log_p_pt  = 0., p_ProbNNp       = 0.;
	Float_t log_piMinIpChi2 = 0., pi_ghostProbNN = 0.;
	Float_t log_pi_pt = 0.;
	Float_t BDTK      = 0.;
	Float_t chi2array[200] {0.};

	Double_t lbMinIpChi2   = 0., lbDira_ownPV = 0., lbFd_ownPV = 0.;
	Double_t lEndVtxChi2    = 0.;
	Double_t jpsiMinIpChi2 = 0., jpsiMass     = 0.;
	Double_t lFdChi2_orivX = 0., lDira_orivX  = 0., lFd_orivX  = 0.;

	Double_t lDira_ownPV = 0., lDm           = 0., lMinIpChi2 = 0.;
	Double_t pMinIpChi2  = 0., pGhostProbNN  = 0.;
	Double_t pPT         = 0., pProbNNp      = 0.;
	Double_t piMinIpChi2 = 0., piGhostProbNN = 0.;

	Double_t piPT   = 0.;
	Double_t myBDTK = 0., BDT        = 0.;

	TMVA::Tools::Instance();// This loads the library
	TMVA::Reader *reader = nullptr;

	//set up input, output
	if(!isData) // Apply on MC
	{
		rootFolder = Form("rootFiles/mcFiles/JpsiLambda/%s/run%d",folder,run);
		if(isoFlag)
		{
			if(!zeroFlag)
			{
				fileIn      = TFile::Open(Form("%s/%s_cutoutks_%s_nonZeroTracks_noPID.root",
				                               rootFolder,part,type));
				treeIn      = (TTree*)fileIn->Get("MyTuple");
				fileIn_iso  = TFile::Open(Form("%s/%s_%s_iso%d_%s_noPID.root",
				                               rootFolder,part,type,isoConf,isoVersion));
				treeIn_iso  = (TTree*)fileIn_iso->Get("MyTuple");
				fileOut     = new TFile(Form("%s/%s_%s_FinalBDT%d_iso%d_%s_noPID.root",
				                             rootFolder,part,type,bdtConf,isoConf,isoVersion),
				                        "RECREATE");
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/%s_cutoutks_%s_ZeroTracks_noPID.root",
				                           rootFolder,part,type));
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/%s_zeroTracks%s_FinalBDT%d_noPID.root",
				                         rootFolder,part,type,bdtConf),
				                    "RECREATE");
			}
		}
		else
		{
			fileIn  = TFile::Open(Form("%s/%s_cutoutks_%s_noPID.root",
			                           rootFolder,part,type));
			treeIn  = (TTree*)fileIn->Get("MyTuple");
			fileOut = new TFile(Form("%s/%s_%s_FinalBDT%d_noIso_noPID.root",
			                         rootFolder,part,type,bdtConf),"RECREATE");
		}
	}//end MC block
	else // Apply on Data
	{
		rootFolder = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);
		if(flag == 1)
		{
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn      = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks_noPID.root",
					                               rootFolder,type));
					treeIn      = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso  = TFile::Open(Form("%s/jpsilambda_%s_iso%d_%s_noPID.root",
					                               rootFolder,type,isoConf,isoVersion));
					treeIn_iso  = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut     = new TFile(Form("%s/jpsilambda_%s_FinalBDT%d_iso%d_%s_noPID.root",
					                             rootFolder,type,bdtConf,isoConf,isoVersion),
					                        "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_ZeroTracks_noPID.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsilambda_zeroTracks%s_FinalBDT%d_noPID.root",
					                         rootFolder,type,bdtConf),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_noPID.root",
				                           rootFolder,type));
				//has to be hadded file doesnt exist
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsilambda_%s_FinalBDT%d_noIso_noPID.root",
				                         rootFolder,type,bdtConf),
				                    "RECREATE");
			}
		}
		else if(flag == 2)
		{
			if(isoFlag)
			{
				if(!zeroFlag)
				{
					fileIn     = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks_noPID.root",
					                              rootFolder,type));
					treeIn     = (TTree*)fileIn->Get("MyTuple");
					fileIn_iso = TFile::Open(Form("%s/jpsilambda_%ssig_iso%d_%s_noPID.root",
					                              rootFolder,type,isoConf,isoVersion));
					treeIn_iso = (TTree*)fileIn_iso->Get("MyTuple");
					fileOut    = new TFile(Form("%s/jpsilambda_%ssig_FinalBDT%d_iso%d_%s_noPID.root",
					                            rootFolder,type,bdtConf,isoConf,isoVersion),
					                       "RECREATE");
				}
				else
				{
					fileIn  = TFile::Open(Form("%s/jpsilambda_%s_withsw_ZeroTracks_noPID.root",
					                           rootFolder,type));
					treeIn  = (TTree*)fileIn->Get("MyTuple");
					fileOut = new TFile(Form("%s/jpsilambda_zeroTracks%ssig_FinalBDT%d_noPID.root",
					                         rootFolder,type,bdtConf),
					                    "RECREATE");
				}
			}
			else
			{
				fileIn  = TFile::Open(Form("%s/jpsilambda_%s_withsw_noPID.root",
				                           rootFolder,type));
				treeIn  = (TTree*)fileIn->Get("MyTuple");
				fileOut = new TFile(Form("%s/jpsilambda_%ssig_FinalBDT%d_noIso_noPID.root",
				                         rootFolder,type,bdtConf),"RECREATE");
			}
		}
	}//end Data block
	//end setup of input, output
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

	reader->AddVariable("log_jpsiminipchi2   := log10(Jpsi_MINIPCHI2)", &log_jpsiMinIpChi2);
	reader->AddVariable("log_jpsimass        := log10(Jpsi_M)", &log_jpsiMass);

	reader->AddVariable("log_lfdchi2         := log10(L_FDCHI2_ORIVX)", &log_lFdChi2_orivX);
	reader->AddVariable("logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))", &logAcos_lDira_orivX);
	reader->AddVariable("log_lfd_orivx       := log10(L_FD_ORIVX)", &log_lFd_orivX);
	reader->AddVariable("logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))", &logAcos_lDira_ownPV);
	reader->AddVariable("L_dm", &l_dm);
	reader->AddVariable("log_lminipchi2      := log10(L_MINIPCHI2)", &log_lMinIpChi2);
	reader->AddVariable("L_ENDVERTEX_CHI2",&l_EndVtx_Chi2);

	reader->AddVariable("log_pminipchi2      := log10(p_MINIPCHI2)", &log_pMinIpChi2);
	reader->AddVariable("p_ProbNNghost", &p_ghostProbNN);
	reader->AddVariable("log_p_PT            := log10(p_PT)", &log_p_pt);
	reader->AddVariable("p_ProbNNp", &p_ProbNNp);

	reader->AddVariable("pi_ProbNNghost", &pi_ghostProbNN);
	reader->AddVariable("log_piminipchi2     := log10(pi_MINIPCHI2)", &log_piMinIpChi2);
	reader->AddVariable("log_pi_PT           := log10(pi_PT)", &log_pi_pt);

	if(isoFlag && !zeroFlag)
	{
		reader->AddVariable(Form("BDTkMin_%s",isoVersion), &BDTK);
	}

	//***********Load weights file**********************************
	TString dir    = "dataset/weights/";
	TString prefix;
	if(!simFlag)//train on data
	{
		if(isoFlag && !zeroFlag)
		{
			prefix = Form("dataRun%d_iso%d_%s_BDT%d_noPID",
			              run,isoConf,isoVersion,bdtConf);
		}
		else if(!isoFlag || (isoFlag && zeroFlag))
		{
			prefix = Form("dataRun%d_noIso_BDT%d_noPID",
			              run,bdtConf);
		}
	}
	else //train on MC
	{
		if(isoFlag && !zeroFlag)
		{
			prefix = Form("MCRun%d_iso%d_%s_BDT%d_noPID",
			              run,isoConf,isoVersion,bdtConf);
		}
		else if(!isoFlag || (isoFlag && zeroFlag))
		{
			prefix = Form("MCRun%d_noIso_BDT%d_noPID",
			              run,bdtConf);
		}
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
	// treeIn->SetBranchAddress("L_TAU", &lTau);
	//treeIn->SetBranchAddress("Lb_DTF_CTAUS_L", &l_dtfCtauS);
	// treeIn->SetBranchAddress("Lb_P", &lbP);
	// treeIn->SetBranchAddress("Lb_PT", &lbPT);

	treeIn->SetBranchAddress("Jpsi_MINIPCHI2", &jpsiMinIpChi2);
	treeIn->SetBranchAddress("Jpsi_M", &jpsiMass);
	// treeIn->SetBranchAddress("Jpsi_CosTheta", &jpsiCosTheta);
	// treeIn->SetBranchAddress("Jpsi_PT", &jpsiPT);
	// treeIn->SetBranchAddress("Jpsi_P", &jpsiP);

	treeIn->SetBranchAddress("L_FDCHI2_ORIVX", &lFdChi2_orivX);
	treeIn->SetBranchAddress("L_DIRA_ORIVX", &lDira_orivX);
	treeIn->SetBranchAddress("L_FD_ORIVX", &lFd_orivX);
	treeIn->SetBranchAddress("L_DIRA_OWNPV", &lDira_ownPV);
	treeIn->SetBranchAddress("L_dm", &lDm);
	treeIn->SetBranchAddress("L_MINIPCHI2", &lMinIpChi2);
	treeIn->SetBranchAddress("L_ENDVERTEX_CHI2", &l_EndVtx_Chi2);
	// treeIn->SetBranchAddress("L_CosTheta", &lCosTheta);

	treeIn->SetBranchAddress("p_MINIPCHI2", &pMinIpChi2);
	treeIn->SetBranchAddress("p_ProbNNghost", &pGhostProbNN);
	treeIn->SetBranchAddress("p_PT", &pPT);

	if(isData)
	{
		treeIn->SetBranchAddress("p_ProbNNp", &pProbNNp);
	}
	else
	{
		treeIn->SetBranchAddress("p_ProbNNp_corr", &pProbNNp);
	}

	treeIn->SetBranchAddress("pi_ProbNNghost", &piGhostProbNN);
	treeIn->SetBranchAddress("pi_MINIPCHI2", &piMinIpChi2);
	treeIn->SetBranchAddress("pi_PT", &piPT);

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
		l_EndVtx_Chi2        = lEndVtxChi2;
		//	l_cosTheta       = lCosTheta;

		log_pMinIpChi2       = log10(pMinIpChi2);
		p_ghostProbNN        = pGhostProbNN;
		log_p_pt             = log10(pPT);
		p_ProbNNp            = pProbNNp;

		pi_ghostProbNN       = piGhostProbNN;
		log_piMinIpChi2      = log10(piMinIpChi2);
		log_pi_pt            = log10(piPT);

		if(isoFlag && !zeroFlag) BDTK = myBDTK;

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
	// if(logFlag) gROOT->ProcessLine(".>");
	if(logFlag) gSystem->RedirectOutput(0);
}
