#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>

using namespace std;

Double_t getSumOfWeights(Double_t mybdt, TTree *tree, Int_t flag, Int_t bdtConf);

void optimizeFinalBDT(Int_t run = 1, Int_t trackType = 3, const char* version = "v1", Int_t bdtConf = 1, Bool_t isoFlag = true)
{
	TFile *fileIn = nullptr;
	TTree *treeIn = nullptr;
	TH1D *hsig = nullptr, *hbkg = nullptr;

	Int_t sig = 0, bkg = 0, siginit = 0, bkginit = 0;
	Double_t BDT = 0., BDT_max = 0., FOM = 0., FOM_max = 0.;
	Double_t eff_sig = 0., eff_sig_max = 0., eff_bkg = 0., eff_bkg_max = 0.;
	Bool_t logFlag = false;
	TString t = "";
	const char *logFileName = "", *type = "";

	if(!isoFlag) version = "noIso";

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

	logFileName = (trackType == 3) ? (TString::Format("optimizeFinalBDT%d_LL_%s_log.txt",bdtConf,version)) : (TString::Format("optimizeFinalBDT%d_DD_%s_log.txt",bdtConf,version));

	if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%s.txt",run,logFileName));

	cout<<"******************************************"<<endl;
	cout<<"Starting optimizeFinalBDT "<<endl;
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_withFinalBDT%d_%s.root",run,type,bdtConf,version),"READ");
	treeIn = (TTree*)fileIn->Get("MyTuple");

	treeIn->Draw(TString::Format("BDT%d>>hsig0(90,-0.4,0.5)",bdtConf),"SW","goff");
	treeIn->Draw(TString::Format("BDT%d>>hsig(90,-0.4,0.5)",bdtConf),"SW*(Lb_DTF_M_JpsiLConstr < 5700)","goff");//Should there be a lower bound as well?
	treeIn->Draw(TString::Format("BDT%d>>hsig1(90,-0.4,0.5)",bdtConf),"SW*(Lb_DTF_M_JpsiLConstr < 5700 && Lb_DTF_M_JpsiLConstr > 5300)","goff");
	treeIn->Draw(TString::Format("BDT%d>>hbkg(90,-0.4,0.5)",bdtConf),"BW*(Lb_DTF_M_JpsiLConstr > 5700)","goff");

	hsig = (TH1D*)gDirectory->Get("hsig0");//play around with this
	hbkg = (TH1D*)gDirectory->Get("hbkg");

	siginit =  getSumOfWeights(-0.5,treeIn,1,bdtConf);
	bkginit =  getSumOfWeights(-0.5,treeIn,2,bdtConf);

	cout<<"siginit = "<<siginit<<endl;
	cout<<"bkginit = "<<bkginit<<endl;

	for(Int_t i = 1; i < 90; i++) {

		BDT = hsig->GetBinCenter(i);

		t = Form("BDT > %f",BDT);

		// sig = getSumOfWeights(BDT, sigtree,1);
		// bkg = getSumOfWeights(BDT, bkgtree,2);

		sig = hsig->Integral(i,90);
		bkg = hbkg->Integral(i,90);

		cout<<"SIG = "<<sig<<" BKG = "<<bkg;
		eff_sig = (Double_t) sig/siginit;
		eff_bkg = (Double_t) bkg/bkginit;

		if(sig + bkg > 0)
			FOM = (Double_t)eff_sig/sqrt(sig+bkg);
		else{
			cout<<"NO"<<endl;
			continue;
		}
		if(FOM > FOM_max) {
			FOM_max = FOM;
			BDT_max = BDT;
			eff_sig_max = eff_sig;
			eff_bkg_max = eff_bkg;
		}

		cout<<"For BDT = "<<BDT<<" FOM = "<<FOM<<" sig_eff = "<<eff_sig*100<<"% bkg_eff = "<<eff_bkg*100<<"%"<<endl;
	}
	cout<<"MAXIMUM FOM = "<<FOM_max<<" at BDT = "<<BDT_max<<" with sig_eff = "<<eff_sig_max*100<<"% and bkg_eff = "<<eff_bkg_max*100<<"%"<<endl;
}
Double_t getSumOfWeights(Double_t mybdt, TTree *tree, Int_t flag, Int_t bdtConf)
{
	Double_t sum = 0., mybdt1 = 0., bmass = 0.;
	Float_t myweight = 0.;
	// TTree *treecut = (TTree*)tree->CopyTree(Form("BDT > %f",mybdt));

	Int_t nentries = tree->GetEntries();

	if(flag == 1)
		tree->SetBranchAddress("SW",&myweight);
	else if(flag == 2)
		tree->SetBranchAddress("BW",&myweight);
	tree->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&bmass);
	tree->SetBranchAddress(TString::Format("BDT%d",bdtConf),&mybdt1);

	for(Int_t i=0; i<nentries; i++) {
		tree->GetEntry(i);

		if(mybdt1 > mybdt)
		{
			if((flag == 2 && bmass > 5700) || flag == 1)
			{
				sum += myweight;
			}
		}
	}
	return sum;
}
