#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

using namespace std;

void MASTER()
{
	//Data- Run 1
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L trigger.C++");
	// gROOT->ProcessLine("trigger(1, true, 0, true, true)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L sanity.C++");
	// gROOT->ProcessLine("sanity(1, true, 0)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L cutOutKs.C++");
	// gROOT->ProcessLine("cutOutKs(1, true, 0, 3)"); //LL
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L DosPlot.C++");
	// gROOT->ProcessLine("DosPlot(1, 3)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationIso.C++");
	// gROOT->ProcessLine("TMVAClassificationIso(1, 3, \"v1\")"); //include other versions later
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationApplicationIso.C++");
	// gROOT->ProcessLine("TMVAClassificationApplicationIso(1, true, 0, 3, 1 , \"v1\", 1)"); //include other versions and conf later
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationApplicationIso.C++");
	// gROOT->ProcessLine("TMVAClassificationApplicationIso(1, true, 0, 3, 2 , \"v1\", 1)"); //include other versions and conf later
	// gROOT->Reset();
	//
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TMVAClassification-JpsiLambda.C++");
	gROOT->ProcessLine("TMVAClassification_JpsiLambda(1, 3, \"v1\")");
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TMVAClassificationApplication-JpsiLambda.C++");
	gROOT->ProcessLine("TMVAClassificationApplication_JpsiLambda(1, true, 0, 3, \"v1\", 1, 1)");
	gROOT->Reset();

	//Data- Run 2
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L trigger.C++");
	// gROOT->ProcessLine("trigger(2, true, 0, true, true)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L sanity.C++");
	// gROOT->ProcessLine("sanity(2, true, 0)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L cutOutKs.C++");
	// gROOT->ProcessLine("cutOutKs(2, true, 0, 3)"); //LL
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L DosPlot.C++");
	// gROOT->ProcessLine("DosPlot(2, 3)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationIso.C++");
	// gROOT->ProcessLine("TMVAClassificationIso(2, 3, \"v1\")"); //include other versions later
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationApplicationIso.C++");
	// gROOT->ProcessLine("TMVAClassificationApplicationIso(2, true, 0, 3, 1 , \"v1\", 1)"); //include other versions and conf later
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationApplicationIso.C++");
	// gROOT->ProcessLine("TMVAClassificationApplicationIso(2, true, 0, 3, 2 , \"v1\", 1)"); //include other versions and conf later
	// gROOT->Reset();

	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassification-JpsiLambda.C++");
	// gROOT->ProcessLine("TMVAClassification_JpsiLambda(2, 3, \"v1\")");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L TMVAClassificationApplication-JpsiLambda.C++");
	// gROOT->ProcessLine("TMVAClassificationApplication_JpsiLambda(2, true, 0, 3, \"v1\", 1, 1)");
	// gROOT->Reset();
}
