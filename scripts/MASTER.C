#include <TROOT.h>
#include <TSystem.h>
#include <iostream>

using namespace std;

void MASTER()
{
	//Data
	// gROOT->ProcessLine(".L trigger.C++");
	// gROOT->ProcessLine("trigger(1, true, 0, true, true)");
	// gROOT->Reset();

	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L sanity.C++");
	// gROOT->ProcessLine("sanity(1, true, 0)");
	// gROOT->Reset();
	//
	// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	// gROOT->ProcessLine(".L cutOutKs.C++");
	// gROOT->ProcessLine("cutOutKs(1, true, 0, 3)"); //LL
	// gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L Lb_sPlot.C++");
	gROOT->ProcessLine("Lb_sPlot(1, 3)");
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TMVAClassificationIso.C++");
	gROOT->ProcessLine("TMVAClassificationIso(1, 3, \"v1\")"); //include other versions later
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TMVAClassificationApplicationIso.C++");
	gROOT->ProcessLine("TMVAClassificationIso(1, true, 0, 3, 1 , \"v1\", 1"); //include other versions and conf later
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TMVAClassification-JpsiLambda.C++");
	gROOT->ProcessLine(".L TMVAClassification_JpsiLambda(1, 3, \"v1\"");
	gROOT->Reset();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
	gROOT->ProcessLine(".L TMVAClassificationApplication-JpsiLambda.C++");
	gROOT->ProcessLine(".L TMVAClassificationApplication_JpsiLambda(1, true, 0, 3, \"v1\", 1, 1)");
	gROOT->Reset();
}
