#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TSystem.h"

#include "RooWorkspace.h"
#include "RooAbsData.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"

#include "iostream"

using namespace RooFit;
using namespace RooStats;
using namespace std;

void GetFeldmanCousins()
{
	gSystem->Load("RooHypatia2_cpp.so");

	TFile *file = TFile::Open("tempModel.root");
// get the workspace out of the file
	RooWorkspace* w = (RooWorkspace*) file->Get("w");
	if(!w) {
		cout <<"workspace not found" << endl;
		return;
	}

	// get the modelConfig out of the file
	ModelConfig* mc = (ModelConfig*) w->obj("ModelConfig");

	// get the modelConfig out of the file
	RooAbsData* data = w->data("combData");

	// make sure ingredients are found
	if(!data || !mc) {
		w->Print();
		cout << "data or ModelConfig was not found" <<endl;
		return;
	}

	/////////////////////////////////////////////
	// create and use the FeldmanCousins tool
	// to find and plot the 95% confidence interval
	// on the parameter of interest as specified
	// in the model config
	FeldmanCousins fc(*data,*mc);
	fc.SetConfidenceLevel(0.95); // 95% interval
	fc.AdditionalNToysFactor(0.1); // to speed up the result
	fc.UseAdaptiveSampling(true); // speed it up a bit
	fc.SetNBins(10); // set how many points per parameter of interest to scan
	fc.CreateConfBelt(true); // save the information in the belt for plotting

	// Since this tool needs to throw toy MC the PDF needs to be
	// extended or the tool needs to know how many entries in a dataset
	// per pseudo experiment.
	// In the 'number counting form' where the entries in the dataset
	// are counts, and not values of discriminating variables, the
	// datasets typically only have one entry and the PDF is not
	// extended.
	if(!mc->GetPdf()->canBeExtended()) {
		if(data->numEntries()==1)
			fc.FluctuateNumDataEntries(false);
		else
			cout <<"Not sure what to do about this model" <<endl;
	}

	// We can use PROOF to speed things along in parallel
	//  ProofConfig pc(*w, 1, "workers=4", kFALSE);
	//  ToyMCSampler*  toymcsampler = (ToyMCSampler*) fc.GetTestStatSampler();
	//  toymcsampler->SetProofConfig(&pc); // enable proof


	// Now get the interval
	PointSetInterval* interval = fc.GetInterval();
	ConfidenceBelt* belt = fc.GetConfidenceBelt();

	// print out the iterval on the first Parameter of Interest
	RooRealVar* firstPOI = (RooRealVar*) mc->GetParametersOfInterest()->first();
	cout << "\n95% interval on " <<firstPOI->GetName()<<" is : ["<<
	        interval->LowerLimit(*firstPOI) << ", "<<
	        interval->UpperLimit(*firstPOI) <<"] "<<endl;

	//////////////////////////////////////////////
	// No nice plots yet, so plot the belt by hand

	// Ask the calculator which points were scanned
	RooDataSet* parameterScan = (RooDataSet*) fc.GetPointsToScan();
	RooArgSet* tmpPoint;

	// make a histogram of parameter vs. threshold
	TH1F* histOfThresholds = new TH1F("histOfThresholds","",
	                                  parameterScan->numEntries(),
	                                  firstPOI->getMin(),
	                                  firstPOI->getMax());

	// loop through the points that were tested and ask confidence belt
	// what the upper/lower thresholds were.
	// For FeldmanCousins, the lower cut off is always 0
	for(Int_t i=0; i<parameterScan->numEntries(); ++i) {
		tmpPoint = (RooArgSet*) parameterScan->get(i)->clone("temp");
		double arMax = belt->GetAcceptanceRegionMax(*tmpPoint);
		double arMin = belt->GetAcceptanceRegionMax(*tmpPoint);
		double poiVal = tmpPoint->getRealValue(firstPOI->GetName());
		histOfThresholds->Fill(poiVal,arMax);
	}
	histOfThresholds->SetMinimum(0);
	histOfThresholds->Draw();
}
