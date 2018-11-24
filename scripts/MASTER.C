#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include <iostream>

using namespace std;

void MASTER()
{
	Bool_t testing = true;         // when true, analysis will only run over a subset of data
	Bool_t loose = true;           // when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Int_t run = 1;                 // Run 1 or Run 2?
	Bool_t isData = true;          // Data or MC?
	Int_t mcType = 0;              // 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType = 3;            // 3 for LL, 5 for DD.
	Int_t isoConf = 5;                   // config for isolation BDT. Currently 1 or 5 supported
	// Int_t finalBDTconf = 1;        // config for final BDT. Currently 1 or 5 supported
	Bool_t isoFlag = true;         // when true, isolation will be used in final BDT.
	const char* isoVersion = "v1"; // which version of isolation BDT? changes input variables used in isolation. v0,1,2,3,4 supported.


	//Data- Run 1
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L Trigger.C++");
	// gROOT->ProcessLine(TString::Format("Trigger(%d, %d, %d, %d, %d)",run, isData, mcType, testing, loose).Data());
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L Sanity.C++");
	// gROOT->ProcessLine(TString::Format("Sanity(%d, %d, %d)",run,isData,mcType).Data());
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L CutOutKs.C++");
	// gROOT->ProcessLine(TString::Format("CutOutKs(%d, %d, %d, %d)",run,isData,mcType,trackType).Data()); //LL
	// gROOT->Reset();

	// if(isData)
	// {
	//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	//      gROOT->ProcessLine(".L DosPlot.C++");
	//      gROOT->ProcessLine(TString::Format("DosPlot(%d, %d)",run,trackType).Data());
	//      gROOT->Reset();
	//
	//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	//      gROOT->ProcessLine(".L TrainIsolation.C++");
	//      gROOT->ProcessLine(TString::Format("TrainIsolation(%d, %d, \"%s\")",run,trackType,isoVersion).Data()); //include other versions later
	//      gROOT->Reset();
	// }
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L ApplyIsolation.C++");
	// gROOT->ProcessLine(TString::Format("ApplyIsolation(%d, %d, %d, %d, 1 , \"%s\", %d)", run, isData, mcType, trackType, isoVersion, isoConf).Data()); //apply on all data/MC
	// gROOT->Reset();
	//
	// if(isData)
	// {
	//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	//      gROOT->ProcessLine(".L ApplyIsolation.C++");
	//      gROOT->ProcessLine(TString::Format("ApplyIsolation(%d, %d, %d, %d, 2 , \"%s\", %d)", run, isData, mcType, trackType, isoVersion, isoConf).Data()); //apply only on sWeighted data
	//      gROOT->Reset();
	//
	//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	//      gROOT->ProcessLine(".L TrainFinalBDT.C++");
	//      gROOT->ProcessLine(TString::Format("TrainFinalBDT(%d, %d, \"%s\", %d, %d)", run, trackType, isoVersion, isoConf, isoFlag).Data());
	//      gROOT->Reset();
	// }
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	// gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 1, 1, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag).Data());//Apply final conf 1 on all data/MC
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// // gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	// gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 5, 1, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag).Data());//Apply final conf 5 on all data/MC
	// gROOT->Reset();
	//
	if(isData)
	{
		//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		//      gROOT->ProcessLine(".L ApplyFinalBDT.C++");
		//      gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 1, 2, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag).Data());//Apply final conf 1 only on weighted data
		//      gROOT->Reset();
		//
		//      gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		//      // gROOT->ProcessLine(".L ApplyFinalBDT.C++");
		//      gROOT->ProcessLine(TString::Format("ApplyFinalBDT(%d, %d, %d, %d, \"%s\", %d, 5, 2, %d)", run, isData, mcType, trackType, isoVersion, isoConf, isoFlag).Data());//Apply final conf 5 only on weighted data
		//      gROOT->Reset();

		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		gROOT->ProcessLine(".L OptimizeFinalBDT.C++");
		gROOT->ProcessLine(TString::Format("OptimizeFinalBDT(%d, %d, \"%s\", %d, 1, %d)", run, trackType, isoVersion, isoConf, isoFlag).Data());//finalBDTconf1
		gROOT->Reset();

		gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L OptimizeFinalBDT.C++");
		gROOT->ProcessLine(TString::Format("OptimizeFinalBDT(%d, %d, \"%s\", %d, 5, %d)", run, trackType, isoVersion, isoConf, isoFlag).Data());//finalBDTconf5
		gROOT->Reset();
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
