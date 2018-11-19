#include <TFile.h>
#include <TTree.h>
#include <iostream>

using namespace std;
void optimizeFinalBDT(Int_t run = 1, Int_t trackType = 3, const char* version = "v1")
{
	TFile *fileIn = nullptr;
	TTree *treeIn = nullptr;
	const char *logFileName = "", *type = "";

	if(trackType == 3)
	{
		cout<<"Processing LL"<<endl;
		type = "LL";
	}
	else if(trackType == 5)
	{
		cout<<"Processing DD"<<endl;
		type = "DD";
	}

	fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_withfinalBDT_%s.root",run,type,version),"READ");
	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->Draw("BDT>>hsig(90,-0.4,0.5)","SW*(Lb_DTF_M_JpsiLConstr < )");
}
