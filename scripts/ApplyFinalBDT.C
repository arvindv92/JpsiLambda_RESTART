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

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
using namespace std;

void ApplyFinalBDT(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0, Int_t trackType = 3, const char* version = "v1", Int_t bdtConf = 1, Int_t flag = 1)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
   version = "v1","v2" or "v3"
   flag = 1 when applying on all data, flag = 2 when applying only on signal training sample
 */
{
	TStopwatch sw;
	sw.Start();

	Bool_t logFlag = false, genFlag = false;
	const char *logFileName = "", *type = "";
	TString outfileName = "";
	TFile *fileIn = NULL, *fileOut = NULL;
	TTree *treeIn = NULL, *treeOut = NULL;

	Int_t entries_init = 0;
	fstream genFile;

	logFileName = (trackType == 3) ? (TString::Format("finalBDTApplication_LL_%s_log.txt",version)) : (TString::Format("finalBDTApplication_DD_%s_log.txt",version));

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");//This could be problematic when putting all scripts together in a master script.
	cout<<"WD = "<<gSystem->pwd()<<endl;

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


	if(!isData) // MC
	{
		switch(mcType)
		{
		case 1: //JpsiLambda
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiLambda/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiLambda/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda_%s_withiso_%s.root",run,type,version));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiLambda/run%d/jpsilambda%s_withfinalBDT_%s.root",run,type,version),"RECREATE");
			break;

		case 2: //JpsiSigma
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiSigma/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiSigma/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma_%s_withiso_%s.root",run,type,version));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiSigma/run%d/jpsisigma%s_withfinalBDT_%s.root",run,type,version),"RECREATE");
			break;

		case 3: //JpsiXi
			if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/mc/JpsiXi/run%d/%s",run,logFileName));
			if(!gSystem->AccessPathName(TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)))
			{
				genFile.open((TString::Format("logs/mc/JpsiXi/run%d/gen_log.txt",run)).Data());
				genFlag = 1;
			}
			fileIn  = TFile::Open(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi_%s_withiso_%s.root",run,type,version));
			fileOut = new TFile(TString::Format("rootFiles/mcFiles/JpsiXi/run%d/jpsixi%s_withfinalBDT_%s.root",run,type,version),"RECREATE");
			break;
		}
	}
	else // Data
	{
		if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%s.txt",run,logFileName));
		if(flag == 1)
		{
			fileIn      = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withiso_%s.root",run,type,version));
			fileOut     = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda%s_withfinalBDT_%s.root",run,type,version),"RECREATE");
		}
		else if(flag == 2)
		{
			fileIn  = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw.root",run,type));
			fileOut = new TFile(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_withfinalBDT_%s.root",run,type,version),"RECREATE");
		}
	}

	treeIn = (TTree*)fileIn->Get("MyTuple");
	treeOut = (TTree*)treeIn->CloneTree(0);

	cout<<"******************************************"<<endl;
	cout<<"Input file = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	entries_init = treeIn->GetEntries();

	cout<<"Incoming Entries = "<<entries_init;

	#ifdef __CINT__
	gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
	#endif
	TString myMethodList = "";
	//---------------------------------------------------------------

	// This loads the library
	TMVA::Tools::Instance();

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;

	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used

	Float_t log_dtfChi2 = 0., log_lbMinIpChi2 = 0., logAcos_lbDira_ownPV = 0., log_lbFd_ownPV = 0., log_lTau = 0. /*, l_dtfCtauS = 0.*/;
	Float_t log_jpsiMinIpChi2 = 0., log_jpsiMass = 0., jpsi_cosTheta = 0., jpsi_pt = 0.;
	Float_t log_lFdChi2_orivX = 0., logAcos_lDira_orivX = 0., log_lFd_orivX = 0., logAcos_lDira_ownPV = 0., l_dm = 0., log_lMinIpChi2 = 0., l_pt = 0., l_endVertex_Chi2 = 0., l_cosTheta = 0.;
	Float_t p_PIDp = 0., log_pMinIpChi2 = 0., p_ghostProb = 0., log_p_pt = 0., p_probNNp = 0.;
	Float_t pi_PIDK = 0., log_piMinIpChi2 = 0., pi_ghostProb = 0., log_pi_pt = 0., pi_probNNpi = 0.;
	Float_t bdtK = 0.;

	reader->AddVariable( "log_dtfchi2 := log10(Lb_ConsLb_chi2)", &log_dtfChi2 );
	reader->AddVariable( "log_lbminipchi2 := log10(Lb_MINIPCHI2)", &log_lbMinIpChi2 );
	reader->AddVariable( "logacos_lbdira := log10(acos(Lb_DIRA_OWNPV))", &logAcos_lbDira_ownPV );
	reader->AddVariable( "log_lbfd_ownpv := log10(Lb_FD_OWNPV)", &log_lbFd_ownPV );
	reader->AddVariable( "log_ltau := log10(L_TAU)", &log_lTau );
	//reader->AddVariable( "Lb_DTF_CTAUS_L", &l_dtfCtauS );

	reader->AddVariable( "log_jpsiminipchi2 := log10(Jpsi_MINIPCHI2)", &log_jpsiMinIpChi2 );
	reader->AddVariable( "log_jpsimass := log10(Jpsi_M)", &log_jpsiMass );
	reader->AddVariable( "Jpsi_CosTheta", &jpsi_cosTheta );
	reader->AddVariable( "Jpsi_PT", &jpsi_pt );

	reader->AddVariable( "log_lfdchi2 := log10(L_FDCHI2_ORIVX)", &log_lFdChi2_orivX );
	reader->AddVariable( "logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))", &logAcos_lDira_orivX );
	reader->AddVariable( "log_lfd_orivx := log10(L_FD_ORIVX)", &log_lFd_orivX );
	reader->AddVariable( "logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))", &logAcos_lDira_ownPV );
	reader->AddVariable( "L_dm", &l_dm );
	reader->AddVariable( "log_lminipchi2 := log10(L_MINIPCHI2)", &log_lMinIpChi2 );
	reader->AddVariable( "L_PT", &l_pt );
	reader->AddVariable( "L_ENDVERTEX_CHI2", &l_endVertex_Chi2);
	reader->AddVariable( "L_CosTheta", &l_cosTheta );

	reader->AddVariable( "p_PIDp", &p_PIDp );
	reader->AddVariable( "log_pminipchi2 := log10(p_MINIPCHI2)", &log_pMinIpChi2 );
	reader->AddVariable( "p_TRACK_GhostProb", &p_ghostProb );
	reader->AddVariable( "log_p_PT := log10(p_PT)", &log_p_pt);
	reader->AddVariable( "p_ProbNNp", &p_probNNp );

	reader->AddVariable( "pi_TRACK_GhostProb", &pi_ghostProb );
	reader->AddVariable( "pi_PIDK", &pi_PIDK);
	reader->AddVariable( "log_piminipchi2 := log10(pi_MINIPCHI2)", &log_piMinIpChi2);
	reader->AddVariable( "log_pi_PT := log10(pi_PT)", &log_pi_pt );
	reader->AddVariable( "pi_ProbNNpi", &pi_probNNpi );

	reader->AddVariable( TString::Format("BDTkMin_%s",version), &bdtK );

	// --- Book the MVA methods

	TString dir    = "dataset/weights/";
	TString prefix;
	prefix = TString::Format("TMVAClassification-JpsiLambda%s_dataRun%d_%s",type,run,version);

	TString methodName = TString("BDT method");
	TString weightfile = dir + prefix + TString("_") + TString::Format("BDTconf%d",bdtConf) + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile );

	if (!fileIn) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}
	std::cout << "--- TMVAClassificationApp    : Using input file: " << fileIn->GetName() << std::endl;

	// --- Event loop

	// Prepare the event tree
	// - here the variable names have to corresponds to your tree
	// - you can use the same variables as above which is slightly faster,
	//   but of course you can use different ones and copy the values inside the event loop
	//
	std::cout << "--- Select signal sample" << std::endl;
	Float_t chi2array[200];
	Double_t lbMinIpChi2 = 0., lbDira_ownPV = 0., lbFd_ownPV = 0., lTau = 0.;
	Double_t jpsiMinIpChi2 = 0., jpsiMass = 0., jpsiCosTheta = 0., jpsiPT = 0.;
	Double_t lFdChi2_orivX = 0., lDira_orivX = 0., lFd_orivX = 0., lDira_ownPV = 0., lDm = 0., lMinIpChi2 = 0., lPT = 0., lEndVertexChi2 = 0., lCosTheta = 0.;
	Double_t pPIDp = 0., pMinIpChi2 = 0., pGhostProb = 0., pPT = 0., pProbNNp = 0.;
	Double_t piPIDK = 0., piMinIpChi2 = 0., piGhostProb = 0., piPT = 0., piProbNNpi = 0.;
	// Double_t userVar[19], BDT, bmass;
	Double_t myBdtK = 0., BDT = 0.;
	//treeIn->SetBranchAddress("BW", &BW);

	treeIn->SetBranchAddress( "Lb_ConsLb_chi2", chi2array );
	treeIn->SetBranchAddress( "Lb_MINIPCHI2", &lbMinIpChi2 );
	treeIn->SetBranchAddress( "Lb_DIRA_OWNPV", &lbDira_ownPV  );
	treeIn->SetBranchAddress( "Lb_FD_OWNPV", &lbFd_ownPV );
	treeIn->SetBranchAddress( "L_TAU", &lTau);
	//treeIn->SetBranchAddress( "Lb_DTF_CTAUS_L", &l_dtfCtauS );

	treeIn->SetBranchAddress( "Jpsi_MINIPCHI2", &jpsiMinIpChi2 );
	treeIn->SetBranchAddress( "Jpsi_M", &jpsiMass );
	treeIn->SetBranchAddress( "Jpsi_CosTheta", &jpsiCosTheta );
	treeIn->SetBranchAddress( "Jpsi_PT", &jpsiPT );

	treeIn->SetBranchAddress( "L_FDCHI2_ORIVX", &lFdChi2_orivX );
	treeIn->SetBranchAddress( "L_DIRA_ORIVX", &lDira_orivX );
	treeIn->SetBranchAddress( "L_FD_ORIVX", &lFd_orivX );
	treeIn->SetBranchAddress( "L_DIRA_OWNPV", &lDira_ownPV );
	treeIn->SetBranchAddress( "L_dm", &lDm );
	treeIn->SetBranchAddress( "L_MINIPCHI2", &lMinIpChi2 );
	treeIn->SetBranchAddress( "L_PT", &lPT );
	treeIn->SetBranchAddress( "L_ENDVERTEX_CHI2", &lEndVertexChi2 );
	treeIn->SetBranchAddress( "L_CosTheta", &lCosTheta );

	treeIn->SetBranchAddress( "p_PIDp", &pPIDp );
	treeIn->SetBranchAddress( "p_MINIPCHI2", &pMinIpChi2 );
	treeIn->SetBranchAddress( "p_TRACK_GhostProb", &pGhostProb );
	treeIn->SetBranchAddress( "p_PT", &pPT );
	treeIn->SetBranchAddress( "p_ProbNNp", &pProbNNp );

	treeIn->SetBranchAddress( "pi_TRACK_GhostProb", &piGhostProb );
	treeIn->SetBranchAddress( "pi_PIDK", &piPIDK );
	treeIn->SetBranchAddress( "pi_MINIPCHI2", &piMinIpChi2 );
	treeIn->SetBranchAddress( "pi_PT", &piPT );
	treeIn->SetBranchAddress( "pi_ProbNNpi", &piProbNNpi );

	treeIn->SetBranchAddress( TString::Format("BDTkMin_%s",version), &myBdtK );

	treeOut->Branch( "BDT", &BDT, "BDT/D");
	// treeout->Branch( "Lb_ConsLb_chi2", &chi2array[0],"Lb_ConsLb_chi2/F" );
	// treeout->Branch( "Lb_MINIPCHI2", &userVar[1],"Lb_MINIPCHI2/D" );
	// treeout->Branch( "Lb_DIRA_OWNPV", &userVar[2],"Lb_DIRA_OWNPV/D" );
	// treeout->Branch( "Lb_FD_OWNPV", &userVar[3],"Lb_FD_OWNPV/D" );
	// treeout->Branch( "L_TAU", &userVar[4],"L_TAU/D" );
	// treeout->Branch( "Lb_DTF_CTAUS_L", &userVar[5],"Lb_DTF_CTAUS_L/D" );
	//
	// treeout->Branch( "Jpsi_MINIPCHI2", &userVar[6],"Jpsi_MINIPCHI2/D" );
	// treeout->Branch( "Jpsi_M", &userVar[7],"Jpsi_M/D" );
	// treeout->Branch( "Jpsi_CosTheta", &userVar[8],"Jpsi_CosTheta/D" );
	// treeout->Branch( "Jpsi_PT", &userVar[9],"Jpsi_PT/D" );
	//
	// treeout->Branch( "L_FDCHI2_ORIVX", &userVar[10],"L_FDCHI2_ORIVX/D" );
	// treeout->Branch( "L_DIRA_ORIVX", &userVar[11],"L_DIRA_ORIVX/D" );
	// treeout->Branch( "L_FD_ORIVX", &userVar[12],"L_FD_ORIVX/D" );
	// treeout->Branch( "L_DIRA_OWNPV", &userVar[13],"L_DIRA_OWNPV/D" );
	// treeout->Branch( "L_dm", &userVar[14],"L_dm/D" );
	// treeout->Branch( "L_MINIPCHI2", &userVar[15],"L_MINIPCHI2/D" );
	// treeout->Branch( "L_PT", &userVar[16],"L_PT/D" );
	// treeout->Branch( "L_ENDVERTEX_CHI2", &userVar[17],"L_ENDVERTEX_CHI2/D" );
	// treeout->Branch( "L_CosTheta", &userVar[18],"L_CosTheta/D" );
	//
	// treeout->Branch( "p_PIDp", &userVar[19],"p_PIDp/D" );
	// treeout->Branch( "p_MINIPCHI2", &userVar[20],"p_MINIPCHI2/D" );
	// treeout->Branch( "p_TRACK_GhostProb", &userVar[21],"p_TRACK_GhostProb/D" );
	// treeout->Branch( "p_PT", &userVar[22],"p_PT/D" );
	// treeout->Branch( "p_ProbNNp", &userVar[23],"p_ProbNNp/D" );
	//
	// treeout->Branch( "pi_TRACK_GhostProb", &userVar[24],"pi_TRACK_GhostProb/D" );
	// treeout->Branch( "pi_PIDK", &userVar[25], "pi_PIDK/D" );
	// treeout->Branch( "pi_MINIPCHI2", &userVar[26], "pi_MINIPCHI2/D");
	// treeout->Branch( "pi_PT", &userVar[27], "pi_PT/D");
	// treeout->Branch( "pi_ProbNNpi", &userVar[28], "pi_ProbNNpi/D");
	//
	// treeout->Branch( "BDTk", &userVar[29],"BDTk/D" );
	// treeout->Branch( "Lb_DTF_M_JpsiLConstr", &bmass, "Lb_DTF_M_JpsiLConstr/D");

	std::cout << "--- Processing: " << entries_init << " events" << std::endl;

	for (Long64_t ievt = 0; ievt < entries_init; ievt++) {

		if (ievt%50000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

		treeIn->GetEntry(ievt);

		log_dtfChi2 = log10(chi2array[0]);//Lb_ConsLb_chi2
		log_lbMinIpChi2 = log10(lbMinIpChi2);//Lb_MINIPCHI2
		logAcos_lbDira_ownPV = log10(acos(lbDira_ownPV));//Lb_DIRA_OWNPV
		log_lbFd_ownPV = log10(lbFd_ownPV);//Lb_FD_OWNPV
		log_lTau = log10(lTau);//L_TAU
		//	var[5] = userVar[5];//Lb_DTF_CTAUS_L

		log_jpsiMinIpChi2 = log10(jpsiMinIpChi2);//Jpsi_MINIPCHI2
		log_jpsiMass = log10(jpsiMass);//Jpsi_M
		jpsi_cosTheta = jpsiCosTheta;//Jpsi_CosTheta
		jpsi_pt = jpsiPT;//Jpsi_PT

		log_lFdChi2_orivX = log10(lFdChi2_orivX);//L_FDCHI2_ORIVX
		logAcos_lDira_orivX = log10(acos(lDira_orivX));//L_DIRA_ORIVX
		log_lFd_orivX = log10(lFd_orivX);//L_FD_ORIVX
		logAcos_lDira_ownPV = log10(acos(lDira_ownPV));//L_DIRA_OWNPV
		l_dm = lDm;//L_dm
		log_lMinIpChi2 = log10(lMinIpChi2);//L_MINIPCHI2
		l_pt = lPT;//L_PT
		l_endVertex_Chi2 = lEndVertexChi2;//L_ENDVERTEX_CHI2
		l_cosTheta = lCosTheta;//L_CosTheta

		pPIDp = p_PIDp; //p_PIDp
		log_pMinIpChi2 = log10(pMinIpChi2);//p_MINIPCHI2
		p_ghostProb = pGhostProb;//p_TRACK_GhostProb
		log_p_pt = log10(pPT);//p_PT
		p_probNNp = pProbNNp;//p_ProbNNp

		pi_ghostProb = piGhostProb;//pi_TRACK_GhostProb
		pi_PIDK = piPIDK;//pi_PIDK
		log_piMinIpChi2 = log10(piMinIpChi2);//pi_MINIPCHI2
		log_pi_pt = log10(piPT);//pi_PT
		pi_probNNpi = piProbNNpi;//pi_ProbNNpi
		bdtK = myBdtK;//BDTk

		BDT = reader->EvaluateMVA("BDT method");

		treeOut->Fill();
	}
	// Get elapsed time
	sw.Stop();
	std::cout << "--- End of event loop: ";

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	std::cout << "--- Created root file: "<<fileOut->GetName()<<" containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> End of ApplyFinalBDT"; sw.Print();
	if(logFlag) gROOT->ProcessLine(".>");
}
