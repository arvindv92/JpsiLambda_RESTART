#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

using namespace std;

void MASTER()
{
	//Data- Run 1
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L Trigger.C++");
	// gROOT->ProcessLine("Trigger(1, true, 0, true, true)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L Sanity.C++");
	// gROOT->ProcessLine("Sanity(1, true, 0)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L CutOutKs.C++");
	// gROOT->ProcessLine("CutOutKs(1, true, 0, 3)"); //LL
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L DosPlot.C++");
	// gROOT->ProcessLine("DosPlot(1, 3)");
	// gROOT->Reset();

	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TrainIsolation.C++");
	// gROOT->ProcessLine("TrainIsolation(1, 3, \"v1\")"); //include other versions later
	// gROOT->Reset();
	//
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L ApplyIsolation.C++");
	gROOT->ProcessLine("ApplyIsolation(1, true, 0, 3, 1 , \"v4\", 5)"); //include other versions and conf later
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L ApplyIsolation.C++");
	gROOT->ProcessLine("ApplyIsolation(1, true, 0, 3, 2 , \"v4\", 5)"); //include other versions and conf later
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TrainFinalBDT.C++");
	gROOT->ProcessLine("TrainFinalBDT(1, 3, \"v4\", 5, true)");
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	gROOT->ProcessLine("ApplyFinalBDT(1, true, 0, 3, \"v4\", 5, 1, 1, true)");//Apply on all data
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	gROOT->ProcessLine("ApplyFinalBDT(1, true, 0, 3, \"v4\", 5, 5, 1, true)");//Apply on all data
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	gROOT->ProcessLine("ApplyFinalBDT(1, true, 0, 3, \"v4\", 5, 1, 2, true)");//Apply only on weighted data
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L ApplyFinalBDT.C++");
	gROOT->ProcessLine("ApplyFinalBDT(1, true, 0, 3, \"v4\", 5, 5, 2, true)");//Apply only on weighted data
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L optimizeFinalBDT.C++");
	gROOT->ProcessLine("optimizeFinalBDT(1, 3, \"v4\", 5, 1, true)");
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L optimizeFinalBDT.C++");
	gROOT->ProcessLine("optimizeFinalBDT(1, 3, \"v4\", 5, 5, true)");
	gROOT->Reset();

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
