#include "TROOT.h"
#include "TSystem.h"
#include <iostream>
#include <vector>
#include "Fitscript.h"
#include "TGraphErrors.h"
using namespace std;
void RatioOfYields()
{
	gSystem->Exec("date");
	gROOT->ProcessLine(".L Fitscript.C+");
	gSystem->Load("Fitscript_C.so");

	Float_t ratios[10];
	Float_t errs_y[10];
	Float_t errs_x[10];
	Float_t bdt[10];
	std::vector <Float_t> retVals;
	for(Int_t i=0; i<=9; i++)
	{
		Float_t bdtCut = i*0.1;
		bdt[i] = bdtCut;

		retVals = Fitscript(1, 1, 1, "v0", 4700, 1, bdtCut);
		ratios[i] = retVals[0];
		errs_y[i] = retVals[1];
		errs_x[i] = 0.1;
	}
	TGraphErrors *gr = new TGraphErrors(10,bdt,ratios,errs_x,errs_y);
	new TCanvas();
	gr->Draw("ALP");
}
