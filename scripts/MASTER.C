#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include <iostream>

using namespace std;

void MASTER()
{
	TStopwatch sw;
	sw.Start();

	Bool_t testing         = false; // when true, analysis will only run over a subset of data
	Bool_t loose           = true;  // when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run              = 1;     // Run 1 or Run 2?
	Bool_t isData          = true;  // Data or MC?
	Int_t mcType           = 0;     // 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;     // 3 for LL, 5 for DD.
	Int_t massLow          = 5400;
	Int_t massHigh         = 5700;
	Int_t isoConf          = 1;     // config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf     = 1;     // config for final BDT. Currently 1 or 2 supported
	Bool_t isoFlag         = true;  // when true, isolation will be used in final BDT.
	const char* isoVersion = "";    // which version of isolation BDT? changes input variables used in isolation. v0,1,2,3,4 supported.
	Bool_t logFlag         = true;  // set to false only while testing.

	TString inFileName  = (trackType == 3) ? ("jpsilambda_LL_forIsoTraining.root") : ("jpsilambda_DD_forIsoTraining.root");
	TString outFileName = (trackType == 3) ? ("jpsilambda_LL_forIsoTraining_massCut.root") : ("jpsilambda_DD_forIsoTraining_massCut.root");

	TString isoVersionArray[4] = {"v0", "v1", "v2", "v3"};
	Int_t isoConfArray[2]      = {1,2};
	Int_t finalBDTconfArray[2] = {1,2};

	//Data- Run 1
	//Trigger Cut
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L Trigger.C++");
	// gROOT->ProcessLine(TString::Format("Trigger(%d, %d, %d, %d, %d, %d)",run, isData, mcType, testing, loose, logFlag).Data());
	// gROOT->Reset();
	//
	// //Sanity Cuts
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L Sanity.C++");
	// gROOT->ProcessLine(TString::Format("Sanity(%d, %d, %d, %d)", run, isData, mcType, logFlag).Data());
	// gROOT->Reset();
	//
	// //Cut Out Ks0
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L CutOutKs.C++");
	// gROOT->ProcessLine(TString::Format("CutOutKs(%d, %d, %d, %d, %d)", run, isData, mcType, trackType, logFlag).Data()); //LL
	// gROOT->Reset();
	//

	//sWeight data
	// if(isData)
	// {
	//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	//      gROOT->ProcessLine(".L DosPlot.C++");
	//      gROOT->ProcessLine(TString::Format("DosPlot(%d, %d, %d)", run, trackType, logFlag).Data());
	//      gROOT->Reset();
	// }

	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L massCut.C+");
	// gROOT->ProcessLine(TString::Format("massCut(%d, %d, %d, %d, \"%s\", \"%s\", %d)", run, trackType, massLow, massHigh, inFileName.Data(), outFileName.Data(), logFlag).Data()); //include other versions later
	// gROOT->Reset();

	for (Int_t i = 0; i <= 3; i++)//loop over isolation BDT versions.
	{
		isoVersion = (isoVersionArray[i]).Data();
		//Train Isolation BDT on data
		if(isData)
		{
			gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			gROOT->ProcessLine(".L TrainIsolation.C+");
			gROOT->ProcessLine(TString::Format("TrainIsolation(%d, %d, \"%s\", %d)", run, trackType, isoVersion, logFlag).Data());
			gROOT->Reset();
		}

		for (Int_t j = 0; j <= 1; j++)//loop over isolation BDT configurations
		{
			isoConf    = isoConfArray[j];
			//Apply isolation BDT on all data/MC
			gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			gROOT->ProcessLine(".L ApplyIsolation.C+");
			gROOT->ProcessLine(TString::Format("ApplyIsolation(%d, %d, %d, %d, 1 , \"%s\", %d, %d)", run, isData, mcType, trackType, isoVersion, isoConf, logFlag).Data());
			gROOT->Reset();

			if(isData)
			{
				//Apply isolation BDT on sWeighted data/MC
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(TString::Format("ApplyIsolation(%d, %d, %d, %d, 2 , \"%s\", %d, %d)", run, isData, mcType, trackType, isoVersion, isoConf, logFlag).Data());
				gROOT->Reset();

				//Train Final BDT on data w/ isolation (if isoFlag is true)
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(".L TrainFinalBDT.C+");
				gROOT->ProcessLine(TString::Format("TrainFinalBDT(%d, %d, \"%s\", %d, %d, %d)", run, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());
				gROOT->Reset();

				//Train Final BDT on data sans isolation (if isoFlag is true)
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(TString::Format("TrainFinalBDT(%d, %d, \"%s\", %d, %d, %d)", run, trackType, isoVersion, isoConf, !isoFlag, logFlag).Data());
				gROOT->Reset();
			}

			for(Int_t k = 0; k <= 1; k++)   //Loop over final BDT configurations.
			{
				finalBDTconf = finalBDTconfArray[k];

				//Apply final BDT on nonZeroTracks data/MC
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(".L ApplyFinalBDT.C+");
				gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 1, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());//Apply final conf 1 on all data/MC
				gROOT->Reset();

				//Apply final BDT on zeroTracks data/MC
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				//gROOT->ProcessLine(".L ApplyFinalBDT.C++");
				gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 1, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());//Apply final conf 1 on all data/MC
				gROOT->Reset();

				if(isData)
				{
					//Apply final BDT on nonZeroTracks sWeighted data/MC
					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
					//gROOT->ProcessLine(".L ApplyFinalBDT.C++");
					gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 2, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data()); //Apply final conf 1 only on weighted data
					gROOT->Reset();

					//Apply final BDT on zeroTracks sWeighted data/MC
					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
					//gROOT->ProcessLine(".L ApplyFinalBDT.C++");
					gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 2, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data()); //Apply final conf 1 only on weighted data
					gROOT->Reset();

					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
					gROOT->ProcessLine(".L OptimizeFinalBDT.C+");
					gROOT->ProcessLine(TString::Format("OptimizeFinalBDT(%d, %d, \"%s\", %d, %d, %d, %d)", run, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());//finalBDTconf1
					gROOT->Reset();
				}
			}   //end loop on finalBDTconf
		}//end loop on isoConf
	}//end loop on isoVersion
	sw.Stop();
	cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
