#include "TFile.h"
#include "TTree.h"

void CompareMcData(Int_t run = 1, Int_t mcType = 0, Int_t trackType = 3)
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	TFile *fileIn_data = nullptr, *fileIn_mc = nullptr;
	TTree *treeIn_data = nullptr, *treeIn_mc = nullptr;

	const char *type = (trackType == 3) ? ("LL") : ("DD");

	const char *folder = "", *part = "";
	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Bu_JpsiX";
		part = "bu_jpsix";
		break;
	}
	case 5:
	{
		folder = "Bd_JpsiX";
		part = "bd_jpsix";
		break;
	}
	}

	fileIn_data = TFile::Open(Form("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw.root",run,type),"READ");
	treeIn_data = (TTree*)fileIn_data->Get("MyTuple");

	fileIn_mc   = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s_triggered.root",folder,run,part,type),"READ");
	treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");



}
