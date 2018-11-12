/**********************************************************************************
* Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
* Package   : TMVA                                                               *
* Exectuable: TMVAClassificationApplication                                      *
*                                                                                *
* This macro provides a simple example on how to use the trained classifiers     *
* within an analysis module                                                      *
**********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
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

void TMVAClassificationApplication_JpsiLambda(Int_t run = 1,Int_t trackType = 3, TString version = "v1", Int_t bdtConf = 1)//flag = 1 for LL, flag = 2 for DD, flag = 3 for LL sig, flag = 4 for DD sig, 5 for LL bkg, 6 for DD bkg
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   trackType = 3 for LL, 5 for DD.
   version = "v1","v2" or "v3"
 */
{
#ifdef __CINT__
	gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif
	TString myMethodList = "";
	//---------------------------------------------------------------

	// This loads the library
	TMVA::Tools::Instance();

	std::cout << std::endl;
	std::cout << "==> Start TMVAClassificationApplication" << std::endl;
	// --------------------------------------------------------------------------------------------------

	// --- Create the Reader object

	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used

	Float_t var[30];

	reader->AddVariable( "log_dtfchi2 := log10(Lb_ConsLb_chi2)", &var[0] );
	reader->AddVariable( "log_lbminipchi2 := log10(Lb_MINIPCHI2)", &var[1] );
	reader->AddVariable( "logacos_lbdira := log10(acos(Lb_DIRA_OWNPV))", &var[2] );
	reader->AddVariable( "log_lbfd_ownpv := log10(Lb_FD_OWNPV)", &var[3] );
	reader->AddVariable( "log_dtfctaul := log10(Lb_DTF_CTAU_L)", &var[4] );
	reader->AddVariable( "Lb_DTF_CTAUS_L", &var[5] );

	reader->AddVariable( "log_jpsiminipchi2 := log10(Jpsi_MINIPCHI2)", &var[6] );
	reader->AddVariable( "log_jpsimass := log10(Jpsi_M)", &var[7] );
	reader->AddVariable( "Jpsi_CosTheta", &var[8] );
	reader->AddVariable( "Jpsi_PT", &var[9] );

	reader->AddVariable( "log_lfdchi2 := log10(L_FDCHI2_ORIVX)", &var[10] );
	reader->AddVariable( "logacos_ldira_orivx := log10(acos(L_DIRA_ORIVX))", &var[11] );
	reader->AddVariable( "log_lfd_orivx := log10(L_FD_ORIVX)", &var[12] );
	reader->AddVariable( "logacos_ldira_ownpv := log10(acos(L_DIRA_OWNPV))", &var[13] );
	reader->AddVariable( "L_dm", &var[14] );
	reader->AddVariable( "log_lminipchi2 := log10(L_MINIPCHI2)", &var[15] );
	reader->AddVariable( "L_PT", &var[16] );
	reader->AddVariable( "L_ENDVERTEX_CHI2", &var[17] );
	reader->AddVariable( "L_CosTheta", &var[18] );

	reader->AddVariable( "p_PIDp", &var[19] );
	reader->AddVariable( "log_pminipchi2 := log10(p_MINIPCHI2)", &var[20] );
	reader->AddVariable( "p_TRACK_GhostProb", &var[21] );
	reader->AddVariable( "log_p_PT := log10(p_PT)", &var[22] );
	reader->AddVariable( "p_ProbNNp", &var[23] );

	reader->AddVariable( "pi_TRACK_GhostProb", &var[24] );
	reader->AddVariable( "pi_PIDK", &var[25]);
	reader->AddVariable( "log_piminipchi2 := log10(pi_MINIPCHI2)", &var[26]);
	reader->AddVariable( "log_pi_PT := log10(pi_PT)", &var[27] );
	reader->AddVariable( "pi_ProbNNpi", &var[28] );

	reader->AddVariable( "BDTk", &var[29] );

	// Spectator variables declared in the training have to be added to the reader, too
	// Float_t spec1,spec2;
	// reader->AddSpectator( "spec1 := var1*2",   &spec1 );
	// reader->AddSpectator( "spec2 := var1*3",   &spec2 );

	// --- Book the MVA methods

	TString dir    = "weights/";
	TString prefix;
	if(flag == 1||flag==3||flag==5)
		prefix = "TMVAClassification-JpsiLambda_LL_data";
	else if(flag == 2||flag==4||flag==6)
		prefix = "TMVAClassification-JpsiLambda_DD_data";

	// Book method(s)
	// for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
	//   //   if (it->second) {
	//   TString methodName = TString("BDT method");
	//   TString weightfile = dir + prefix + TString("_") + TString("BDTconf5") + TString(".weights.xml");
	//   reader->BookMVA( methodName, weightfile );
	//       // }
	// }
	TString methodName = TString("BDT method");
	TString weightfile = dir + prefix + TString("_") + TString("BDTconf5") + TString(".weights.xml");
	reader->BookMVA( methodName, weightfile );
	// Book output histograms
	UInt_t nbin = 100;

	// TH1F  *histBdt(0);

	// histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );

	// Prepare input tree (this must be replaced by your data source)
	// in this example, there is a toy tree with signal and one with background events
	// we'll later on use only the "signal" events for the test in this example.
	//
	TFile *input(0);
	TString fname;

	switch(flag) {
	case 1:
		fname = "applicationfiles/jpsilambda_LL_forfinalBDTApp.root";
		cout<<"TMVA application on all LL"<<endl;
		break;
	case 2:
		fname = "applicationfiles/jpsilambda_DD_forfinalBDTApp.root";
		cout<<"TMVA application on all DD"<<endl;
		break;
	case 3:
		fname = "trainingfiles/jpsilambda_LL_sig.root";
		cout<<"TMVA application on LL sig"<<endl;
		break;
	case 4:
		fname = "trainingfiles/jpsilambda_DD_sig.root";
		cout<<"TMVA application on DD sig"<<endl;
		break;
	case 5:
		fname = "trainingfiles/jpsilambda_LL_bkg.root";
		cout<<"TMVA application on LL bkg"<<endl;
		break;
	case 6:
		fname = "trainingfiles/jpsilambda_DD_bkg.root";
		cout<<"TMVA application on DD bkg"<<endl;
		break;
	}
	input = TFile::Open(fname,"READ");

	if (!input) {
		std::cout << "ERROR: could not open data file" << std::endl;
		exit(1);
	}
	std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;

	TFile *target(0);
	TTree *treeout(0);
	TString targetname;

	switch(flag) {
	case 1:
		targetname = "applicationfiles/TMVApp-JpsiLambda_LL_BDTconf5.root";
		break;
	case 2:
		targetname = "applicationfiles/TMVApp-JpsiLambda_DD_BDTconf5.root";
		break;
	case 3:
		targetname = "applicationfiles/TMVApp-JpsiLambda_LLsig_BDTconf5.root";
		break;
	case 4:
		targetname = "applicationfiles/TMVApp-JpsiLambda_DDsig_BDTconf5.root";
		break;
	case 5:
		targetname = "applicationfiles/TMVApp-JpsiLambda_LLbkg_BDTconf5.root";
		break;
	case 6:
		targetname = "applicationfiles/TMVApp-JpsiLambda_DDbkg_BDTconf5.root";
		break;
	}

	target = new TFile(targetname,"RECREATE" );

	treeout = new TTree("MyTuple","BDT stuff");

	// --- Event loop

	// Prepare the event tree
	// - here the variable names have to corresponds to your tree
	// - you can use the same variables as above which is slightly faster,
	//   but of course you can use different ones and copy the values inside the event loop
	//
	std::cout << "--- Select signal sample" << std::endl;
	TTree* theTree = (TTree*)input->Get("MyTuple");
	Float_t chi2array[200];

	// Double_t userVar[19], BDT, bmass;
	Double_t userVar[30], BDT, bmass;

	Double_t SW,BW;
	//theTree->SetBranchAddress("BW", &BW);

	theTree->SetBranchAddress( "Lb_ConsLb_chi2", chi2array );
	theTree->SetBranchAddress( "Lb_MINIPCHI2", &userVar[1] );
	theTree->SetBranchAddress( "Lb_DIRA_OWNPV", &userVar[2]  );
	theTree->SetBranchAddress( "Lb_FD_OWNPV", &userVar[3] );
	theTree->SetBranchAddress( "Lb_DTF_CTAU_L", &userVar[4] );
	theTree->SetBranchAddress( "Lb_DTF_CTAUS_L", &userVar[5] );

	theTree->SetBranchAddress( "Jpsi_MINIPCHI2", &userVar[6] );
	theTree->SetBranchAddress( "Jpsi_M", &userVar[7] );
	theTree->SetBranchAddress( "Jpsi_CosTheta", &userVar[8] );
	theTree->SetBranchAddress( "Jpsi_PT", &userVar[9] );

	theTree->SetBranchAddress( "L_FDCHI2_ORIVX", &userVar[10] );
	theTree->SetBranchAddress( "L_DIRA_ORIVX", &userVar[11] );
	theTree->SetBranchAddress( "L_FD_ORIVX", &userVar[12] );
	theTree->SetBranchAddress( "L_DIRA_OWNPV", &userVar[13] );
	theTree->SetBranchAddress( "L_dm", &userVar[14] );
	theTree->SetBranchAddress( "L_MINIPCHI2", &userVar[15] );
	theTree->SetBranchAddress( "L_PT", &userVar[16] );
	theTree->SetBranchAddress( "L_ENDVERTEX_CHI2", &userVar[17] );
	theTree->SetBranchAddress( "L_CosTheta", &userVar[18] );

	theTree->SetBranchAddress( "p_PIDp", &userVar[19] );
	theTree->SetBranchAddress( "p_MINIPCHI2", &userVar[20] );
	theTree->SetBranchAddress( "p_TRACK_GhostProb", &userVar[21] );
	theTree->SetBranchAddress( "p_PT", &userVar[22] );
	theTree->SetBranchAddress( "p_ProbNNp", &userVar[23] );

	theTree->SetBranchAddress( "pi_TRACK_GhostProb", &userVar[24] );
	theTree->SetBranchAddress( "pi_PIDK", &userVar[25] );
	theTree->SetBranchAddress( "pi_MINIPCHI2", &userVar[26] );
	theTree->SetBranchAddress( "pi_PT", &userVar[27] );
	theTree->SetBranchAddress( "pi_ProbNNpi", &userVar[28] );

	theTree->SetBranchAddress( "BDTk", &userVar[29] );

	theTree->SetBranchAddress( "Lb_DTF_M_JpsiLConstr",&bmass);

	if(flag!=1 && flag!=2) {
		theTree->SetBranchAddress( "SW",&SW);
		theTree->SetBranchAddress( "BW",&BW);
	}

	treeout->Branch( "BDT", &BDT, "BDT/D");
	treeout->Branch( "Lb_ConsLb_chi2", &chi2array[0],"Lb_ConsLb_chi2/F" );
	treeout->Branch( "Lb_MINIPCHI2", &userVar[1],"Lb_MINIPCHI2/D" );
	treeout->Branch( "Lb_DIRA_OWNPV", &userVar[2],"Lb_DIRA_OWNPV/D" );
	treeout->Branch( "Lb_FD_OWNPV", &userVar[3],"Lb_FD_OWNPV/D" );
	treeout->Branch( "Lb_DTF_CTAU_L", &userVar[4],"Lb_DTF_CTAU_L/D" );
	treeout->Branch( "Lb_DTF_CTAUS_L", &userVar[5],"Lb_DTF_CTAUS_L/D" );

	treeout->Branch( "Jpsi_MINIPCHI2", &userVar[6],"Jpsi_MINIPCHI2/D" );
	treeout->Branch( "Jpsi_M", &userVar[7],"Jpsi_M/D" );
	treeout->Branch( "Jpsi_CosTheta", &userVar[8],"Jpsi_CosTheta/D" );
	treeout->Branch( "Jpsi_PT", &userVar[9],"Jpsi_PT/D" );

	treeout->Branch( "L_FDCHI2_ORIVX", &userVar[10],"L_FDCHI2_ORIVX/D" );
	treeout->Branch( "L_DIRA_ORIVX", &userVar[11],"L_DIRA_ORIVX/D" );
	treeout->Branch( "L_FD_ORIVX", &userVar[12],"L_FD_ORIVX/D" );
	treeout->Branch( "L_DIRA_OWNPV", &userVar[13],"L_DIRA_OWNPV/D" );
	treeout->Branch( "L_dm", &userVar[14],"L_dm/D" );
	treeout->Branch( "L_MINIPCHI2", &userVar[15],"L_MINIPCHI2/D" );
	treeout->Branch( "L_PT", &userVar[16],"L_PT/D" );
	treeout->Branch( "L_ENDVERTEX_CHI2", &userVar[17],"L_ENDVERTEX_CHI2/D" );
	treeout->Branch( "L_CosTheta", &userVar[18],"L_CosTheta/D" );

	treeout->Branch( "p_PIDp", &userVar[19],"p_PIDp/D" );
	treeout->Branch( "p_MINIPCHI2", &userVar[20],"p_MINIPCHI2/D" );
	treeout->Branch( "p_TRACK_GhostProb", &userVar[21],"p_TRACK_GhostProb/D" );
	treeout->Branch( "p_PT", &userVar[22],"p_PT/D" );
	treeout->Branch( "p_ProbNNp", &userVar[23],"p_ProbNNp/D" );

	treeout->Branch( "pi_TRACK_GhostProb", &userVar[24],"pi_TRACK_GhostProb/D" );
	treeout->Branch( "pi_PIDK", &userVar[25], "pi_PIDK/D" );
	treeout->Branch( "pi_MINIPCHI2", &userVar[26], "pi_MINIPCHI2/D");
	treeout->Branch( "pi_PT", &userVar[27], "pi_PT/D");
	treeout->Branch( "pi_ProbNNpi", &userVar[28], "pi_ProbNNpi/D");

	treeout->Branch( "BDTk", &userVar[29],"BDTk/D" );
	treeout->Branch( "Lb_DTF_M_JpsiLConstr", &bmass, "Lb_DTF_M_JpsiLConstr/D");

	if(flag!=1 && flag!=2) {
		treeout->Branch( "SW", &SW,"SW/D" );
		treeout->Branch( "BW", &BW,"BW/D" );
	}

	std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
	TStopwatch sw;
	sw.Start();
	for (Long64_t ievt=0; ievt<theTree->GetEntries(); ievt++) {

		if (ievt%50000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

		theTree->GetEntry(ievt);

		var[0] = log10(chi2array[0]);//Lb_ConsLb_chi2
		var[1] = log10(userVar[1]);//Lb_MINIPCHI2
		var[2] = log10(acos(userVar[2]));//Lb_DIRA_OWNPV
		var[3] = log10(userVar[3]);//Lb_FD_OWNPV
		var[4] = log10(userVar[4]);//Lb_DTF_CTAU_L
		var[5] = userVar[5];//Lb_DTF_CTAUS_L

		var[6] = log10(userVar[6]);//Jpsi_MINIPCHI2
		var[7] = log10(userVar[7]);//Jpsi_M
		var[8] = userVar[8];//Jpsi_CosTheta
		var[9] = userVar[9];//Jpsi_PT

		var[10] = log10(userVar[10]);//L_FDCHI2_ORIVX
		var[11] = log10(acos(userVar[11]));//L_DIRA_ORIVX
		var[12] = log10(userVar[12]);//L_FD_ORIVX
		var[13] = log10(acos(userVar[13]));//L_DIRA_OWNPV
		var[14] = userVar[14];//L_dm
		var[15] = log10(userVar[15]);//L_MINIPCHI2
		var[16] = userVar[16];//L_PT
		var[17] = userVar[17];//L_ENDVERTEX_CHI2
		var[18] = userVar[18];//L_CosTheta

		var[19] = userVar[19];//p_PIDp
		var[20] = log10(userVar[20]);//p_MINIPCHI2
		var[21] = userVar[21];//p_TRACK_GhostProb
		var[22] = log10(userVar[22]);//p_PT
		var[23] = userVar[23];//p_ProbNNp

		var[24] = userVar[24];//pi_TRACK_GhostProb
		var[25] = userVar[25];//pi_PIDK
		var[26] = log10(userVar[26]);//pi_MINIPCHI2
		var[27] = log10(userVar[27]);//pi_PT
		var[28] = userVar[28];//pi_ProbNNpi

		var[29] = userVar[29];//BDTk

		// --- Return the MVA outputs and fill into histograms

		BDT = reader->EvaluateMVA("BDT method");

		//    histBdt ->Fill( reader->EvaluateMVA( "BDT method"           ) );

		treeout->Fill();
	}
	treeout->Write();
	treeout->Print();
	// Get elapsed time
	sw.Stop();
	std::cout << "--- End of event loop: "; sw.Print();

	// --- Write histograms

	//   histBdt    ->Write();

	target->Close();

	std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;

	delete reader;

	std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
}
