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

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
using namespace std;

void TMVAClassificationApplication_JpsiLambda(Int_t run = 1, Bool_t isData = true, Int_t mcType = 0, Int_t trackType = 3, const char* version = "v1", Int_t bdtConf = 1, Int_t flag = 1)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
   version = "v1","v2" or "v3"
   flag = 1 when applying on all data, flag = 2 when applying only on signal training sample
 */
{
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

	Float_t log_dtfChi2 = 0., log_lbMinIpChi2 = 0., logAcos_lbDira_ownPV = 0., log_lbFd_ownPV = 0., log_lbDtfCTauL = 0., l_dtfCtauS = 0.;
	Float_t log_jpsiMinIpChi2 = 0., log_jpsiMass = 0., jpsi_cosTheta = 0., jpsi_pt = 0.;
	Float_t log_lFdChi2_orivX = 0., logAcos_lDira_orivX = 0., log_lFd_orivX = 0., logAcos_lDira_ownPV = 0., l_dm = 0., log_lMinIpChi2 = 0., l_pt = 0., l_endVertex_Chi2 = 0., l_cosTheta = 0.;
	Float_t p_PIDp = 0., log_pMinIpChi2 = 0., p_ghostProb = 0., log_p_pt = 0., p_probNNp = 0.;
	Float_t pi_PIDK = 0., log_piMinIpChi2 = 0., pi_ghostProb = 0., log_pi_pt = 0., pi_probNNpi = 0.;
	Float_t bdtK = 0.;

	reader->AddVariable( "log_dtfchi2 := log10(Lb_ConsLb_chi2)", &log_dtfChi2 );
	reader->AddVariable( "log_lbminipchi2 := log10(Lb_MINIPCHI2)", &log_lbMinIpChi2 );
	reader->AddVariable( "logacos_lbdira := log10(acos(Lb_DIRA_OWNPV))", &logAcos_lbDira_ownPV );
	reader->AddVariable( "log_lbfd_ownpv := log10(Lb_FD_OWNPV)", &log_lbFd_ownPV );
	reader->AddVariable( "log_dtfctaul := log10(Lb_DTF_CTAU_L)", &log_lbDtfCTauL );
	reader->AddVariable( "Lb_DTF_CTAUS_L", &l_dtfCtauS );

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

	TString dir    = "weights/";
	TString prefix;
	prefix = TString::Format("TMVAClassification-JpsiLambda_%s_data_%s",type,version);

	// Book method(s)
	// for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	//   //   if (it->second) {
	//   TString methodName = TString("BDT method");
	//   TString weightfile = dir + prefix + TString("_") + TString("BDTconf5") + TString(".weights.xml");
	//   reader->BookMVA( methodName, weightfile );
	//       // }
	// }
	TString methodName = TString("BDT method");
	TString weightfile = dir + prefix + TString("_") + TString::Format("BDTconf%d",bdtConf) + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile );
	// Book output histograms

	// Prepare input tree (this must be replaced by your data source)
	// in this example, there is a toy tree with signal and one with background events
	// we'll later on use only the "signal" events for the test in this example.
	//
	// TFile *input(0);
	// TString fname;
	//
	// switch(flag) {
	// case 1:
	//      fname = "applicationfiles/jpsilambda_LL_forfinalBDTApp.root";
	//      cout<<"TMVA application on all LL"<<endl;
	//      break;
	// case 2:
	//      fname = "applicationfiles/jpsilambda_DD_forfinalBDTApp.root";
	//      cout<<"TMVA application on all DD"<<endl;
	//      break;
	// case 3:
	//      fname = "trainingfiles/jpsilambda_LL_sig.root";
	//      cout<<"TMVA application on LL sig"<<endl;
	//      break;
	// case 4:
	//      fname = "trainingfiles/jpsilambda_DD_sig.root";
	//      cout<<"TMVA application on DD sig"<<endl;
	//      break;
	// case 5:
	//      fname = "trainingfiles/jpsilambda_LL_bkg.root";
	//      cout<<"TMVA application on LL bkg"<<endl;
	//      break;
	// case 6:
	//      fname = "trainingfiles/jpsilambda_DD_bkg.root";
	//      cout<<"TMVA application on DD bkg"<<endl;
	//      break;
	// }
	//input = TFile::Open(fname,"READ");

	if (!fileIn) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}
	std::cout << "--- TMVAClassificationApp    : Using input file: " << fileIn->GetName() << std::endl;

	// TFile *target(0);
	// TTree *treeout(0);
	// TString targetname;
	//
	// switch(flag) {
	// case 1:
	//      targetname = "applicationfiles/TMVApp-JpsiLambda_LL_BDTconf5.root";
	//      break;
	// case 2:
	//      targetname = "applicationfiles/TMVApp-JpsiLambda_DD_BDTconf5.root";
	//      break;
	// case 3:
	//      targetname = "applicationfiles/TMVApp-JpsiLambda_LLsig_BDTconf5.root";
	//      break;
	// case 4:
	//      targetname = "applicationfiles/TMVApp-JpsiLambda_DDsig_BDTconf5.root";
	//      break;
	// case 5:
	//      targetname = "applicationfiles/TMVApp-JpsiLambda_LLbkg_BDTconf5.root";
	//      break;
	// case 6:
	//      targetname = "applicationfiles/TMVApp-JpsiLambda_DDbkg_BDTconf5.root";
	//      break;
	// }
	//
	// target = new TFile(targetname,"RECREATE" );
	//
	// treeout = new TTree("MyTuple","BDT stuff");

	// --- Event loop

	// Prepare the event tree
	// - here the variable names have to corresponds to your tree
	// - you can use the same variables as above which is slightly faster,
	//   but of course you can use different ones and copy the values inside the event loop
	//
	std::cout << "--- Select signal sample" << std::endl;
	Float_t chi2array[200];
	Float_t lbMinIpChi2 = 0., lbDira_ownPV = 0., lbFd_ownPV = 0., lbDtfCTauL = 0.;
	Float_t jpsiMinIpChi2 = 0., jpsiMass = 0.;
	Float_t lFdChi2_orivX = 0., lDira_orivX = 0., lFd_orivX = 0., lDira_ownPV = 0., lMinIpChi2 = 0.;
	Float_t pMinIpChi2 = 0., p_pt = 0.;
	Float_t piMinIpChi2 = 0., pi_pt = 0.;
	// Double_t userVar[19], BDT, bmass;
	Double_t BDT;
	//treeIn->SetBranchAddress("BW", &BW);

	treeIn->SetBranchAddress( "Lb_ConsLb_chi2", chi2array );
	treeIn->SetBranchAddress( "Lb_MINIPCHI2", &lbMinIpChi2 );
	treeIn->SetBranchAddress( "Lb_DIRA_OWNPV", &lbDira_ownPV  );
	treeIn->SetBranchAddress( "Lb_FD_OWNPV", &lbFd_ownPV );
	treeIn->SetBranchAddress( "Lb_DTF_CTAU_L", &lbDtfCTauL);
	treeIn->SetBranchAddress( "Lb_DTF_CTAUS_L", &l_dtfCtauS );

	treeIn->SetBranchAddress( "Jpsi_MINIPCHI2", &jpsiMinIpChi2 );
	treeIn->SetBranchAddress( "Jpsi_M", &jpsiMass );
	treeIn->SetBranchAddress( "Jpsi_CosTheta", &jpsi_cosTheta );
	treeIn->SetBranchAddress( "Jpsi_PT", &jpsi_pt );

	treeIn->SetBranchAddress( "L_FDCHI2_ORIVX", &lFdChi2_orivX );
	treeIn->SetBranchAddress( "L_DIRA_ORIVX", &lDira_orivX );
	treeIn->SetBranchAddress( "L_FD_ORIVX", &lFd_orivX );
	treeIn->SetBranchAddress( "L_DIRA_OWNPV", &lDira_ownPV );
	treeIn->SetBranchAddress( "L_dm", &l_dm );
	treeIn->SetBranchAddress( "L_MINIPCHI2", &lMinIpChi2 );
	treeIn->SetBranchAddress( "L_PT", &l_pt );
	treeIn->SetBranchAddress( "L_ENDVERTEX_CHI2", &l_endVertex_Chi2 );
	treeIn->SetBranchAddress( "L_CosTheta", &l_cosTheta );

	treeIn->SetBranchAddress( "p_PIDp", &p_PIDp );
	treeIn->SetBranchAddress( "p_MINIPCHI2", &pMinIpChi2 );
	treeIn->SetBranchAddress( "p_TRACK_GhostProb", &p_ghostProb );
	treeIn->SetBranchAddress( "p_PT", &p_pt );
	treeIn->SetBranchAddress( "p_ProbNNp", &p_probNNp );

	treeIn->SetBranchAddress( "pi_TRACK_GhostProb", &pi_ghostProb );
	treeIn->SetBranchAddress( "pi_PIDK", &pi_PIDK );
	treeIn->SetBranchAddress( "pi_MINIPCHI2", &piMinIpChi2 );
	treeIn->SetBranchAddress( "pi_PT", &pi_pt );
	treeIn->SetBranchAddress( "pi_ProbNNpi", &pi_probNNpi );

	treeIn->SetBranchAddress( TString::Format("BDTkMin_%s",version), &bdtK );

	treeOut->Branch( "BDT", &BDT, "BDT/D");
	// treeout->Branch( "Lb_ConsLb_chi2", &chi2array[0],"Lb_ConsLb_chi2/F" );
	// treeout->Branch( "Lb_MINIPCHI2", &userVar[1],"Lb_MINIPCHI2/D" );
	// treeout->Branch( "Lb_DIRA_OWNPV", &userVar[2],"Lb_DIRA_OWNPV/D" );
	// treeout->Branch( "Lb_FD_OWNPV", &userVar[3],"Lb_FD_OWNPV/D" );
	// treeout->Branch( "Lb_DTF_CTAU_L", &userVar[4],"Lb_DTF_CTAU_L/D" );
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
	TStopwatch sw;
	sw.Start();
	for (Long64_t ievt = 0; ievt < entries_init; ievt++) {

		if (ievt%50000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

		treeIn->GetEntry(ievt);

		log_dtfChi2 = log10(chi2array[0]);//Lb_ConsLb_chi2
		log_lbMinIpChi2 = log10(lbMinIpChi2);//Lb_MINIPCHI2
		logAcos_lbDira_ownPV = log10(acos(lbDira_ownPV));//Lb_DIRA_OWNPV
		log_lbFd_ownPV = log10(lbFd_ownPV);//Lb_FD_OWNPV
		log_lbDtfCTauL = log10(lbDtfCTauL);//Lb_DTF_CTAU_L
		//	var[5] = userVar[5];//Lb_DTF_CTAUS_L

		log_jpsiMinIpChi2 = log10(jpsiMinIpChi2);//Jpsi_MINIPCHI2
		log_jpsiMass = log10(jpsiMass);//Jpsi_M
		//	var[8] = userVar[8];//Jpsi_CosTheta
		//	var[9] = userVar[9];//Jpsi_PT

		log_lFdChi2_orivX = log10(lFdChi2_orivX);//L_FDCHI2_ORIVX
		logAcos_lDira_orivX = log10(acos(lDira_orivX));//L_DIRA_ORIVX
		log_lFd_orivX = log10(lFd_orivX);//L_FD_ORIVX
		logAcos_lDira_ownPV = log10(acos(lDira_ownPV));//L_DIRA_OWNPV
		//var[14] = userVar[14];//L_dm
		log_lMinIpChi2 = log10(lMinIpChi2);//L_MINIPCHI2
		//var[16] = userVar[16];//L_PT
		//var[17] = userVar[17];//L_ENDVERTEX_CHI2
		//var[18] = userVar[18];//L_CosTheta

		// var[19] = userVar[19];//p_PIDp
		log_pMinIpChi2 = log10(pMinIpChi2);//p_MINIPCHI2
		// var[21] = userVar[21];//p_TRACK_GhostProb
		log_p_pt = log10(p_pt);//p_PT
		// var[23] = userVar[23];//p_ProbNNp

		// var[24] = userVar[24];//pi_TRACK_GhostProb
		// var[25] = userVar[25];//pi_PIDK
		log_piMinIpChi2 = log10(piMinIpChi2);//pi_MINIPCHI2
		log_pi_pt = log10(pi_pt);//pi_PT
		// var[28] = userVar[28];//pi_ProbNNpi
		// var[29] = userVar[29];//BDTk

		BDT = reader->EvaluateMVA("BDT method");

		treeOut->Fill();
	}
	// Get elapsed time
	sw.Stop();
	std::cout << "--- End of event loop: "; sw.Print();

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	std::cout << "--- Created root file: "<<fileOut->GetName()<<" containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
	if(logFlag) gROOT->ProcessLine(".>");
}
