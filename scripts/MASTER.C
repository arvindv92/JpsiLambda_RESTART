#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include "DoSWeight.h"
#include "TrainIsolation.h"
#include "ApplyIsolation.h"
#include "TrainFinalBDT.h"
#include "ApplyFinalBDT.h"
#include "OptimizeFinalBDT.h"
#include <iostream>

using namespace std;

void MASTER()
{
	TStopwatch sw;
	sw.Start();

	gSystem->Exec("date");
	// gROOT->ProcessLine(".L TrainIsolation.C+");
	// gROOT->ProcessLine(".L ApplyIsolation.C+");
	gROOT->ProcessLine(".L TrainFinalBDT.C+");
	gROOT->ProcessLine(".L ApplyFinalBDT.C+");
	gROOT->ProcessLine(".L OptimizeFinalBDT.C+");

	// gSystem->Load("TrainIsolation_C.so");
	// gSystem->Load("ApplyIsolation_C.so");
	gSystem->Load("TrainFinalBDT_C.so");
	gSystem->Load("ApplyFinalBDT_C.so");
	gSystem->Load("OptimizeFinalBDT_C.so");

	Int_t run              = 2;// Run 1 or Run 2?
	Int_t mcType           = 0;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;// 3 for LL, 5 for DD.
	Int_t isoConf          = 1;// config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf     = 1;// config for final BDT. Currently 1 or 2 supported
	// Int_t massLow          = 5400;
	// Int_t massHigh         = 5700;

	Bool_t testing         = false;// when true, analysis will only run over a subset of data
	Bool_t loose           = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Bool_t isData          = true;// Data or MC?
	Bool_t isoFlag         = true;// when true, isolation will be used in final BDT.
	Bool_t logFlag         = true;// set to false only while testing.
	Bool_t newFlag         = false;

	const char* isoVersion = ""; // which version of isolation BDT? changes input variables used in isolation. v0,1,2,3,4 supported.

	// TString inFileName  = (trackType == 3) ? ("jpsilambda_LL_forIsoTraining.root") : ("jpsilambda_DD_forIsoTraining.root");
	// TString outFileName = (trackType == 3) ? ("jpsilambda_LL_forIsoTraining_massCut.root") : ("jpsilambda_DD_forIsoTraining_massCut.root");
	Int_t runArray[2]          = {1,2};
	TString isoVersionArray[2] = {"v0","v1"};
	Int_t isoConfArray[2]      = {1,2};
	Int_t finalBDTconfArray[2] = {1,2};
	Bool_t newFlagArray[2]     = {false,true};

	for(Int_t i = 0; i<=1; i++)
	{
		run = runArray[i];
		//Data- Run 1
		//Trigger Cut
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L Trigger.C+");
		// gROOT->ProcessLine(Form("Trigger(%d, %d, %d, %d, %d, %d)",
		//                         run, isData, mcType, testing, loose, logFlag));
		// gROOT->Reset();

		// //Sanity Cuts
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L Sanity.C+");
		// gROOT->ProcessLine(Form("Sanity(%d, %d, %d, %d)", run, isData, mcType, logFlag));
		// gROOT->Reset();
		//
		//Cut Out Ks0
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L CutOutKs.C+");
		// gROOT->ProcessLine(Form("CutOutKs(%d, %d, %d, %d, %d)",
		//                         run, isData, mcType, trackType, logFlag));
		// gROOT->Reset();


		//sWeight data
		if(isData)
		{
			//sWeight nonZeroTracks data
			//gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			//	gROOT->ProcessLine(".L DoSWeight.C+");
			// gROOT->ProcessLine(Form("DoSWeight(%d, %d, %d, false)",
			// run, trackType, logFlag));

			//DoSWeight(run, trackType, logFlag, false);
			//gROOT->Reset();

			//sWeight ZeroTracks data
			//gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			// gROOT->ProcessLine(Form("DoSWeight(%d, %d, %d, true)",
			//                         run, trackType, logFlag));
			//DoSWeight(run, trackType, logFlag, true);
			// gROOT->Reset();
		}

		for (Int_t j = 0; j <= 1; j++)//loop over isolation BDT versions.
		{
			isoVersion = (isoVersionArray[j]);
			//Train Isolation BDT on data
			// cout<<"***TrainIsolation run "<<run<<"isoVersion "<<isoVersion<<"***"<<endl;
			// TrainIsolation(run, trackType, isoVersion, logFlag);

			for (Int_t k = 0; k <= 1; k++)//loop over isolation BDT configurations
			{
				isoConf    = isoConfArray[k];

				// cout<<"***ApplyIsolation all data run "<<run<<"isoVersion "<<isoVersion<<
				//         "isoConf "<<isoConf<<"***"<<endl;

				//Apply isolation BDT on all data/MC
				// ApplyIsolation(run, isData, mcType, trackType, 1, isoVersion, isoConf, logFlag);

				if(isData)
				{
					// cout<<"***ApplyIsolation sWeight data run "<<run<<"isoVersion "<<isoVersion<<
					//         "isoConf "<<isoConf<<"***"<<endl;
					//Apply isolation BDT on sWeighted data
					// ApplyIsolation(run, isData, mcType, trackType, 2, isoVersion, isoConf, logFlag);

					for(Int_t l = 1; l<=1; l++)
					{
						newFlag = newFlagArray[l];

						cout<<"***TrainFinalBDT nonZeroTracks run "<<run<<"isoVersion "<<isoVersion<<
						        "isoConf "<<isoConf<<"newFlag "<<newFlag<<"***"<<endl;
						//Train Final BDT on data w/ isolation (if isoFlag is true)
						TrainFinalBDT(run, trackType, isoVersion, isoConf, isoFlag, logFlag, newFlag);

						cout<<"***TrainFinalBDT ZeroTracks run "<<run<<"isoVersion "<<isoVersion<<
						        "isoConf "<<isoConf<<"newFlag "<<newFlag<<"***"<<endl;
						//Train Final BDT on data sans isolation (if isoFlag is true)
						TrainFinalBDT(run, trackType, isoVersion, isoConf, !isoFlag, logFlag, newFlag);

						for(Int_t m = 0; m <= 1; m++) //Loop over final BDT configurations.
						{
							finalBDTconf = finalBDTconfArray[m];

							cout<<"***ApplyFinalBDT all nonZeroTracks run "<<run<<"isoVersion "<<
							        isoVersion<<"isoConf "<<isoConf<<"newFlag "<<newFlag<<"finalBDTconf "<<
							        finalBDTconf<<"***"<<endl;

							//Apply final BDT on nonZeroTracks data/MC
							ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, false, logFlag, newFlag);

							cout<<"***ApplyFinalBDT all ZeroTracks run "<<run<<"isoVersion "<<
							        isoVersion<<"isoConf "<<isoConf<<"newFlag "<<newFlag<<"finalBDTconf "<<
							        finalBDTconf<<"***"<<endl;

							//Apply final BDT on zeroTracks data/MC
							ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, true, logFlag, newFlag);

							if(isData)
							{
								cout<<"***ApplyFinalBDT sWeight nonZeroTracks run "<<run<<"isoVersion "<<
								        isoVersion<<"isoConf "<<isoConf<<"newFlag "<<newFlag<<"finalBDTconf "<<
								        finalBDTconf<<"***"<<endl;
								//Apply final BDT on nonZeroTracks sWeighted data/MC
								ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, false, logFlag, newFlag);

								cout<<"***ApplyFinalBDT sWeight ZeroTracks run "<<run<<"isoVersion "<<
								        isoVersion<<"isoConf "<<isoConf<<"newFlag "<<newFlag<<"finalBDTconf "<<
								        finalBDTconf<<"***"<<endl;
								//Apply final BDT on ZeroTracks sWeighted data/MC
								ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, true, logFlag, newFlag);

								cout<<"***OptimizeFinalBDT run "<<run<<"isoVersion "<<
								        isoVersion<<"isoConf "<<isoConf<<"newFlag "<<newFlag<<"finalBDTconf "<<
								        finalBDTconf<<"***"<<endl;
								OptimizeFinalBDT(run, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag, newFlag);
							}
						}//end loop on finalBDTconf
					}
				}

			}//end loop on isoConf
		}//end loop on isoVersion
	}

	sw.Stop();
	cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
