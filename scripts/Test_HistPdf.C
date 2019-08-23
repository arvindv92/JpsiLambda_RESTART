#include "TFile.h"
#include "TTree.h"
#include "RooRealVar.h"
#include "RooHistPdf.h"
#include "TH1.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "RooPlot.h"
using namespace std;
using namespace RooFit;

#define Open TFile::Open

void Test_HistPdf()
{
	RooRealVar *myVar = new RooRealVar("Lb_DTF_M_JpsiLConstr","",5000,6000);
	TFile *fileIn = Open("../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run2/jpsisigma_cutoutks_LL_nonZeroTracks.root","READ");
	TTree *treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->AddFriend("MyTuple","../rootFiles/mcFiles/JpsiLambda/JpsiSigma/run2/jpsisigma_LL_FinalBDT1_iso2_v1.root");
	treeIn->Draw("Lb_DTF_M_JpsiLConstr>>h(250,5000,6000)","BDT1 > 0.475");

	TH1D *h = (TH1D*)gDirectory->Get("h");

	RooDataHist *ds = new RooDataHist("ds","ds",*myVar,h);
	cout<<"Printing RooDataHist"<<endl;
	ds->Print();

	RooHistPdf *myPDF = new RooHistPdf("myPDF","",*myVar,*ds);

	RooPlot *myFrame = myVar->frame();
	ds->plotOn(myFrame, LineColor(kBlue));
	myPDF->plotOn(myFrame, LineColor(kRed));

	new TCanvas();
	myFrame->Draw();
}
