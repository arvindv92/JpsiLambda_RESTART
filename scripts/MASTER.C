#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>

using namespace std;

void MASTER()
{
	Bool_t testing         = false; // when true, analysis will only run over a subset of data
	Bool_t loose           = true;  // when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run              = 1;     // Run 1 or Run 2?
	Bool_t isData          = true;  // Data or MC?
	Int_t mcType           = 0;     // 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;     // 3 for LL, 5 for DD.
	Int_t isoConf          = 1;     // config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf     = 1;     // config for final BDT. Currently 1 or 2 supported
	Bool_t isoFlag         = true;  // when true, isolation will be used in final BDT.
	const char* isoVersion = "";    // which version of isolation BDT? changes input variables used in isolation. v0,1,2,3,4 supported.
	Bool_t logFlag         = true;  // set to false only while testing.

	TString isoVersionArray[4] = {"v0", "v1", "v2", "v3"};
	Int_t isoConfArray[2]      = {1,2};
	Int_t finalBDTconfArray[2] = {1,2};

	//Data- Run 1
	//Trigger Cut
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L Trigger.C++");
	gROOT->ProcessLine(TString::Format("Trigger(%d, %d, %d, %d, %d, %d)",run, isData, mcType, testing, loose, logFlag).Data());
	gROOT->Reset();

	//Sanity Cuts
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L Sanity.C++");
	gROOT->ProcessLine(TString::Format("Sanity(%d, %d, %d, %d)", run, isData, mcType, logFlag).Data());
	gROOT->Reset();

	//Cut Out Ks0
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L CutOutKs.C++");
	gROOT->ProcessLine(TString::Format("CutOutKs(%d, %d, %d, %d, %d)", run, isData, mcType, trackType, logFlag).Data()); //LL
	gROOT->Reset();

	//sWeight data
	if(isData)
	{
		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L DosPlot.C++");
		gROOT->ProcessLine(TString::Format("DosPlot(%d, %d, %d)", run, trackType, logFlag).Data());
		gROOT->Reset();
	}


	for (Int_t i = 0; i <= 3; i++)
	{
		isoVersion = (isoVersionArray[i]).Data();

		//Train Isolation BDT on data
		if(isData)
		{
			gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			gROOT->ProcessLine(".L TrainIsolation.C++");
			gROOT->ProcessLine(TString::Format("TrainIsolation(%d, %d, \"%s\", %d)", run, trackType, isoVersion, logFlag).Data()); //include other versions later
			gROOT->Reset();
		}

		for (Int_t j = 0; j <= 1; j++)
		{
			isoConf    = isoConfArray[j];

			//Apply isolation BDT on all data/MC
			gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
			gROOT->ProcessLine(".L ApplyIsolation.C++");
			gROOT->ProcessLine(TString::Format("ApplyIsolation(%d, %d, %d, %d, 1 , \"%s\", %d, %d)", run, isData, mcType, trackType, isoVersion, isoConf, logFlag).Data()); //apply on all data/MC
			gROOT->Reset();

			if(isData)
			{
				//Apply isolation BDT on sWeighted data/MC
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(".L ApplyIsolation.C++");
				gROOT->ProcessLine(TString::Format("ApplyIsolation(%d, %d, %d, %d, 2 , \"%s\", %d, %d)", run, isData, mcType, trackType, isoVersion, isoConf, logFlag).Data()); //apply only on sWeighted data
				gROOT->Reset();

				//Train Final BDT on data
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(".L TrainFinalBDT.C++");
				gROOT->ProcessLine(TString::Format("TrainFinalBDT(%d, %d, \"%s\", %d, %d, %d)", run, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());
				gROOT->Reset();
			}

			for(Int_t k = 0; k <= 1; k++)
			{
				finalBDTconf = finalBDTconfArray[k];

				//Apply final BDT on nonZeroTracks data/MC
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(".L ApplyFinalBDT.C++");
				gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 1, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());//Apply final conf 1 on all data/MC
				gROOT->Reset();

				//Apply final BDT on zeroTracks data/MC
				gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
				gROOT->ProcessLine(".L ApplyFinalBDT.C++");
				gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 1, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());//Apply final conf 1 on all data/MC
				gROOT->Reset();

				if(isData)
				{
					//Apply final BDT on nonZeroTracks sWeighted data/MC
					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
					//gROOT->ProcessLine(".L ApplyFinalBDT.C++");
					gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 2, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());     //Apply final conf 1 only on weighted data
					gROOT->Reset();

					//Apply final BDT on zeroTracks sWeighted data/MC
					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
					//gROOT->ProcessLine(".L ApplyFinalBDT.C++");
					gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, %d, 2, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag).Data());     //Apply final conf 1 only on weighted data
					gROOT->Reset();
				}
			}//end loop on finalBDTconf
		}//end loop on isoConf
	}//end loop on isoVersion
	/*if(isData)
	   {
	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L TrainFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("TrainFinalBDT(%d, %d, \"%s\", %d, %d, %d)", run, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());
	   gROOT->Reset();
	   }

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 1, 1, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());//Apply final conf 1 on all data/MC
	   gROOT->Reset();

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   // gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 5, 1, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());//Apply final conf 5 on all data/MC
	   gROOT->Reset();

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 1, 1, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());//Apply final conf 1 on all data/MC
	   gROOT->Reset();

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   // gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 5, 1, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());//Apply final conf 5 on all data/MC
	   gROOT->Reset();

	   if(isData)
	   {
	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   //gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 1, 2, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());     //Apply final conf 1 only on weighted data
	   gROOT->Reset();

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   // gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 5, 2, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());     //Apply final conf 5 only on weighted data
	   gROOT->Reset();

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   //gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 1, 2, %d, true, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());     //Apply final conf 1 only on weighted data
	   gROOT->Reset();

	   gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	   // gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	   gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 5, 2, %d, false, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag, logFlag).Data());     //Apply final conf 5 only on weighted data
	   gROOT->Reset();
	 */
	if(isData)
	{
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L OptimizeFinalBDT.C++");
		// gROOT->ProcessLine(TString::Format("OptimizeFinalBDT(%d, %d, \"%s\", %d, 1, %d)", run, trackType, isoVersion, isoConf, isoFlag).Data());//finalBDTconf1
		// gROOT->Reset();
		//
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// // gROOT->ProcessLine(".L OptimizeFinalBDT.C++");
		// gROOT->ProcessLine(TString::Format("OptimizeFinalBDT(%d, %d, \"%s\", %d, 2, %d)", run, trackType, isoVersion, isoConf, isoFlag).Data());//finalBDTconf5
		// gROOT->Reset();
	}


//Data- Run 2
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L Trigger.C++");
// gROOT->ProcessLine("Trigger(2, true, 0, true, true)");
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L Sanity.C++");
// gROOT->ProcessLine("Sanity(2, true, 0)");
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L CutOutKs.C++");
// gROOT->ProcessLine("CutOutKs(2, true, 0, 3)"); //LL
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L DosPlot.C++");
// gROOT->ProcessLine("DosPlot(2, 3)");
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L TrainIsolation.C++");
// gROOT->ProcessLine("TrainIsolation(2, 3, \"v4\")"); //include other versions later
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L ApplyIsolation.C++");
// gROOT->ProcessLine("ApplyIsolation(2, true, 0, 3, 1 , \"v4\", 1)"); //include other versions and conf later
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L ApplyIsolation.C++");
// gROOT->ProcessLine("ApplyIsolation(2, true, 0, 3, 2 , \"v4\", 1)"); //include other versions and conf later
// gROOT->Reset();

// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L TrainFinalBDT.C++");
// gROOT->ProcessLine("TrainFinalBDT(2, 3, \"v4\")");
// gROOT->Reset();
//
// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
// gROOT->ProcessLine(".L ApplyFinalBDT.C++");
// gROOT->ProcessLine("ApplyFinalBDT(2, true, 0, 3, \"v4\", 1, 1)");
// gROOT->Reset();
}
