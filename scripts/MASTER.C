#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include "Trigger.h"
#include "Sanity.h"
#include "CutOutKs.h"
#include "DoSWeight.h"
#include "TrainIsolation.h"
#include "ApplyIsolation.h"
#include "TrainFinalBDT.h"
#include "ApplyFinalBDT.h"
#include "OptimizeFinalBDT.h"
#include "CutFinalBDT.h"
#include <iostream>
#include <vector>

using namespace std;

void MASTER(Int_t run = 1, Int_t config = 1, Int_t block = 1, Bool_t doFixed = false, Bool_t isData = true)
//MCType = 1 for JpsiLambda MC
//MCType = 2 for JpsiSigma MC
//MCType = 3 for JpsiXi MC
//MCType = 4 for Bu_JpsiX MC
//MCType = 5 for Bd_JpsiX MC
//MCType = 6 for Lambda*(1405)MC
//MCType = 7 for Lambda*(1520)MC
//MCType = 8 for Lambda*(1600)MC
//MCType = 9 for chiC1 MC
{
	TStopwatch sw;
	sw.Start();

	gSystem->Exec("date");

	gROOT->ProcessLine(".L CollateFiles.C+");
	gROOT->ProcessLine(".L Trigger.C+");
	gROOT->ProcessLine(".L Sanity.C+");
	gROOT->ProcessLine(".L CutOutKs.C+");
	gROOT->ProcessLine(".L DoSWeight.C+");
	gROOT->ProcessLine(".L TrainIsolation.C+");
	gROOT->ProcessLine(".L ApplyIsolation.C+");
	gROOT->ProcessLine(".L TrainFinalBDT.C+");
	gROOT->ProcessLine(".L ApplyFinalBDT.C+");
	gROOT->ProcessLine(".L OptimizeFinalBDT.C+");
	//  gROOT->ProcessLine(".L CutFinalBDT.C+");

	gSystem->Load("CollateFiles_C.so");
	gSystem->Load("Trigger_C.so");
	gSystem->Load("Sanity_C.so");
	gSystem->Load("CutOutKs_C.so");
	gSystem->Load("DoSWeight_C.so");
	gSystem->Load("TrainIsolation_C.so");
	gSystem->Load("ApplyIsolation_C.so");
	gSystem->Load("TrainFinalBDT_C.so");
	gSystem->Load("ApplyFinalBDT_C.so");
	gSystem->Load("OptimizeFinalBDT_C.so");
	//  gSystem->Load("CutFinalBDT_C.so");

	cout<<"poop"<<endl;

	Int_t mcType               = 0;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType            = 3;// 3 for LL, 5 for DD.
	Int_t isoConf              = 1;// config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf         = 1;// config for final BDT. Currently 1 or 2 supported
	// Int_t runArray[2]          = {1,2};
	// Int_t isoConfArray[2]      = {1,2};
	// Int_t finalBDTconfArray[2] = {1,2};
	// Int_t mcTypeArray[1]       = {6};
	Int_t year                 = 2018;

	Bool_t testing            = false;// when true, analysis will only run over a subset of data
	Bool_t loose              = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
	// Bool_t isData             = true;// Data or MC?
	Bool_t isoFlag            = true;// when true, isolation will be used in final BDT.
	Bool_t logFlag            = true;// set to false only while testing.

	Float_t bdtCut             = 0.0;
	Float_t bdtCut_ZeroTracks  = 0.0;

	Double_t myChi2_nonZero    = 0.0, myChi2_Zero = 0.0;
	const char* isoVersion     = ""; // which version of isolation BDT? changes input variables used in isolation. v0,1 supported.

	TString FOM                = "Punzi";//Sig for S/sqrt(S+B), Punzi for Punzi FOM with a = 3
	// TString isoVersionArray[2] = {"v0","v1"};

	std::vector <Double_t> bdtCuts;

	switch(config)
	{
	case 1:
	{
		isoVersion   = "v0";
		isoConf      = 1;
		finalBDTconf = 1;
		break;
	}
	case 2:
	{
		isoVersion   = "v0";
		isoConf      = 1;
		finalBDTconf = 2;
		break;
	}
	case 3:
	{
		isoVersion   = "v0";
		isoConf      = 2;
		finalBDTconf = 1;
		break;
	}
	case 4:
	{
		isoVersion   = "v0";
		isoConf      = 2;
		finalBDTconf = 2;
		break;
	}
	case 5:
	{
		isoVersion   = "v1";
		isoConf      = 1;
		finalBDTconf = 1;
		break;
	}
	case 6:
	{
		isoVersion   = "v1";
		isoConf      = 1;
		finalBDTconf = 2;
		break;
	}
	case 7:
	{
		isoVersion   = "v1";
		isoConf      = 2;
		finalBDTconf = 1;
		break;
	}
	case 8:
	{
		isoVersion   = "v1";
		isoConf      = 2;
		finalBDTconf = 2;
		break;
	}
	}
	for(mcType = 10; mcType <=10; mcType++)
	{
		if(block == 1)
		{
			cout<<"$$$$$$$$$$$ Processing MC Type "<<mcType<< " $$$$$$$$$$$$$$"<<endl;

			cout<<"$$$$$$$$$$$ Processing Run "<<run<<" $$$$$$$$$$$"<<endl;
			//Trigger Cut
			cout<<"***Trigger***"<<endl;
			Trigger(run, year, isData, mcType, testing, loose, logFlag);

			//Sanity Cuts
			cout<<"***Sanity***"<<endl;
			Sanity(run, year, isData, mcType, logFlag);

			//Cut Out Ks0
			cout<<"***CutOutKs***"<<endl;
			CutOutKs(run, year, isData, mcType, trackType, logFlag);

			// sWeight data
			if(isData)
			{
				//sWeight nonZeroTracks data
				cout<<"****DoSWeight run "<<run<<" nonZeroTracks"<<endl;

				myChi2_nonZero = DoSWeight(run, trackType, logFlag, false);

				if(myChi2_nonZero > 2.0)
				{
					cout<<"####Fit failed for nonZeroTracks sWeighting. Exiting####"<<endl;
					exit(1);
				}
				//sWeight ZeroTracks data
				cout<<"****DoSWeight run "<<run<<" ZeroTracks"<<endl;

				myChi2_Zero = DoSWeight(run, trackType, logFlag, true);

				if(myChi2_Zero > 2.0)
				{
					cout<<"####Fit failed for ZeroTracks sWeighting. Exiting####"<<endl;
					exit(1);
				}
			}
		}
		else if(block == 2)
		{
			//****ALL THIS STUFF IN THIS BLOCK IS INDEPENDENT OF ISOLATION. ONLY NEEDS TO BE EXECUTED TWICE. ONCE FOR EACH BDT CONF*******

			if(isData)
			{
				//Train Final BDT on data sans isolation
				cout<<"***TrainFinalBDT ZeroTracks run "<<run<<"finalBDTconf "<<finalBDTconf<<" ***"<<endl;

				//                                                isoFlag
				TrainFinalBDT(run, trackType, isoVersion, isoConf, false, finalBDTconf, logFlag);
			}

			//Apply final BDT on zeroTracks data/MC
			cout<<"***ApplyFinalBDT all ZeroTracks run "<<run<<" finalBDTconf "<<
			        finalBDTconf<<" ***"<<endl;
			ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, true, logFlag);

			if(isData)
			{
				//Apply final BDT on ZeroTracks sWeighted data
				cout<<"***ApplyFinalBDT sWeight ZeroTracks run "<<run<<" isoVersion "<<
				        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
				        finalBDTconf<<"***"<<endl;
				ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, true, logFlag);
			}
			// ****************************************************************************************************************************
		}
		else if(block == 3)
		{
			if(isData)
			{
				// Train Isolation BDT on data. 2 configs train separately,
				cout<<"***TrainIsolation run "<<run<<" isoVersion "<<isoVersion<<" isoConf "<<isoConf<<"***"<<endl;
				TrainIsolation(run, trackType, isoVersion, isoConf, logFlag);
			}

			// Apply isolation BDT on all data/MC
			cout<<"***ApplyIsolation all data/MC run "<<run<<" isoVersion "<<isoVersion<<
			        " isoConf "<<isoConf<<"***"<<endl;
			ApplyIsolation(run, isData, mcType, trackType, 1, isoVersion, isoConf, logFlag);

			if(isData)
			{
				//Apply isolation BDT on sWeighted data
				cout<<"***ApplyIsolation sWeight data run "<<run<<" isoVersion "<<isoVersion<<
				        " isoConf "<<isoConf<<"***"<<endl;
				ApplyIsolation(run, isData, mcType, trackType, 2, isoVersion, isoConf, logFlag);
			}
		}
		else if(block == 4)
		{
			if(isData && isoFlag)
			{
				//Train Final BDT on data w/ isolation. 2 configs. Train separately.
				cout<<"***TrainFinalBDT nonZeroTracks run "<<run<<" isoVersion "<<isoVersion<<
				        " isoConf "<<isoConf<<" ***"<<endl;

				//                                                isoFlag
				TrainFinalBDT(run, trackType, isoVersion, isoConf, true, finalBDTconf, logFlag);
			}

			//Apply final BDT on nonZeroTracks data/MC
			cout<<"***ApplyFinalBDT all nonZeroTracks run "<<run<<"isoVersion "<<
			        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
			        finalBDTconf<<" ***"<<endl;
			ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, false, logFlag);

			if(isData)
			{
				// Apply final BDT on nonZeroTracks sWeighted data
				cout<<"***ApplyFinalBDT sWeight nonZeroTracks run "<<run<<" isoVersion "<<
				        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
				        finalBDTconf<<" ***"<<endl;
				ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, false, logFlag);

			}
		}
	}        //end loop on mcTypes

	// Optimize final BDT cut based on some FoM
	// cout<<"***OptimizeFinalBDT run "<<run<<" isoVersion "<<
	//         isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
	//         finalBDTconf<<" FOM "<<FOM<<" ***"<<endl;
	//
	// bdtCuts = OptimizeFinalBDT(run, isoVersion, isoConf,
	//                            finalBDTconf, isoFlag, logFlag, FOM,
	//                            "sigma");

	// bdtCut = bdtCuts[0];
	// bdtCut_ZeroTracks = bdtCuts[1];
	// // Make cut on final BDT at optimum point for both nonZeroTracks and ZeroTracks, and combine the cut files.
	// cout<<"***CutFinalBDT run "<<run<<" isoVersion "<<
	//         isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
	//         finalBDTconf<<" FOM "<<FOM<<" nonZeroCut "<<bdtCut<<" ZeroCut "<<bdtCut_ZeroTracks<<" ***"<<endl;
	// CutFinalBDT(run, isData, mcType, isoVersion, isoConf,
	//             finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
	//             logFlag, FOM, "sigma");
	// CutFinalBDT(run, !isData, mcType, isoVersion, isoConf,
	//             finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
	//             logFlag, FOM, "sigma");
	sw.Stop();
	cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
