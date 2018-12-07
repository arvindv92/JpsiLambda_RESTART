#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TList.h"
#include "TIterator.h"
#include <iostream>

using namespace std;

Double_t getSumOfWeights(Double_t mybdt, TTree *tree, Int_t flag, Int_t bdtConf, Int_t nEntries);

void OptimizeFinalBDT(Int_t run = 1, Int_t trackType = 3, const char* isoVersion = "v1", Int_t isoConf = 1, Int_t bdtConf = 1, Bool_t isoFlag = true, Bool_t logFlag = false)
{
	TFile *fileIn = nullptr, *fileIn_zeroTracks = nullptr;
	TTree *treeIn = nullptr, *treeIn_zeroTracks = nullptr, *myTree = nullptr;
	TH1D *hsig = nullptr, *hbkg = nullptr;

	Int_t nEntries = 0;
	TString t = "";
	const char *type = "";
	TString logFileName = "", friendFileName = "", friendFileName_zeroTracks = "";

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
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	logFileName = TString::Format("optimizeFinalBDT%d_%s_iso%d_%s.txt",bdtConf,type,isoConf,isoVersion);

	if(logFlag) gROOT->ProcessLine(TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName.Data()).Data());
	cout<<"logFile = "<<TString::Format(".> logs/data/JpsiLambda/run%d/%s",run,logFileName.Data()).Data()<<endl;
	cout<<"******************************************"<<endl;
	cout<<"Starting optimizeFinalBDT "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	if(isoFlag)
	{
		fileIn                    = TFile  ::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_iso%d_%s_300.root",run,type,isoConf,isoVersion),"READ");
		treeIn                    = (TTree*)fileIn->Get("MyTuple");
		fileIn_zeroTracks         = TFile  ::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_withsw_ZeroTracks.root",run,type),"READ");
		treeIn_zeroTracks         = (TTree*)fileIn_zeroTracks->Get("MyTuple");
		friendFileName            = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion);
		friendFileName_zeroTracks = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_zeroTracks%ssig_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion);
		// fileIn                 = TFile  ::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_FinalBDT%d_iso%d_%s.root",run,type,bdtConf,isoConf,isoVersion),"READ");
	}
	else
	{
		fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_cutoutks_%s.root",run,type));
		treeIn = (TTree*)fileIn->Get("MyTuple");
		friendFileName = TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%s_FinalBDT%d_noIso.root",run,type,bdtConf);
		// fileIn = TFile::Open(TString::Format("rootFiles/dataFiles/JpsiLambda/run%d/jpsilambda_%ssig_FinalBDT%d_noIso.root",run,type,bdtConf),"READ");
	}
	treeIn->AddFriend("MyTuple",friendFileName);
	treeIn_zeroTracks->AddFriend("MyTuple",friendFileName_zeroTracks);

	Int_t ctr = 0;
	TList *list = new TList();
	list->Add(treeIn);
	list->Add(treeIn_zeroTracks);

	TIterator *iter = (TIterator*)list->MakeIterator();
	while((myTree = (TTree*)iter->Next()))
	{
		Double_t BDT = 0., BDT_max = 0., FOM = 0., FOM_max = 0., sig = 0., bkg = 0., siginit = 0., bkginit = 0.;
		Double_t eff_sig = 0., eff_sig_max = 0., eff_bkg = 0., eff_bkg_max = 0.;

		if(ctr == 0)
		{
			cout<<"****************Optimizing for nonZeroTracks****************"<<endl;
		}
		else if(ctr == 1)
		{
			cout<<"****************Optimizing for ZeroTracks****************"<<endl;
		}
		nEntries = myTree->GetEntries();

		myTree->Draw(TString::Format("BDT%d>>hsig0(200,-1.0,1.0)",bdtConf),"SW","goff");
		// myTree->Draw(TString::Format("BDT%d>>hsig(90,-0.4,0.5)",bdtConf),"SW*(Lb_DTF_M_JpsiLConstr < 5700)","goff");//Should there be a lower bound as well?
		// myTree->Draw(TString::Format("BDT%d>>hsig1(90,-0.4,0.5)",bdtConf),"SW*(Lb_DTF_M_JpsiLConstr < 5700 && Lb_DTF_M_JpsiLConstr > 5300)","goff");
		myTree->Draw(TString::Format("BDT%d>>hbkg(200,-1.0,1.0)",bdtConf),"BW*(Lb_DTF_M_JpsiLConstr > 5700)","goff");

		hsig = (TH1D*)gDirectory->Get("hsig0");//play around with this
		hbkg = (TH1D*)gDirectory->Get("hbkg");

		siginit =  (Int_t)getSumOfWeights(-0.5,myTree,1,bdtConf,nEntries);//rounding off happening here
		bkginit =  (Int_t)getSumOfWeights(-0.5,myTree,2,bdtConf,nEntries);

		cout<<"siginit = "<<siginit<<endl;
		cout<<"bkginit = "<<bkginit<<endl;

		for(Int_t i = 1; i < 200; i++) {

			BDT = hsig->GetBinCenter(i);

			t = Form("BDT > %f",BDT);

			// sig = (Int_t)getSumOfWeights(BDT,myTree,1,bdtConf,nEntries);
			// bkg = (Int_t)getSumOfWeights(BDT,myTree,2,bdtConf,nEntries);

			sig = hsig->Integral(i,200);
			bkg = hbkg->Integral(i,200);

			cout<<"SIG = "<<sig<<" BKG = "<<bkg;
			eff_sig = (Double_t) sig/siginit;
			eff_bkg = (Double_t) bkg/bkginit;

			if(sig + bkg > 0)
				FOM = (Double_t)eff_sig/sqrt(sig+bkg);
			else
			{
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
		ctr++;
	}
	if(logFlag) gROOT->ProcessLine(".>");
}
Double_t getSumOfWeights(Double_t mybdt, TTree *tree, Int_t flag, Int_t bdtConf, Int_t nEntries)
{
	Double_t sum = 0., mybdt1 = 0., bmass = 0.;
	Float_t myweight = 0.;
	// TTree *treecut = (TTree*)tree->CopyTree(Form("BDT > %f",mybdt));

	if(flag == 1)
	{
		tree->SetBranchAddress("SW",&myweight);
	}
	else if(flag == 2)
	{
		tree->SetBranchAddress("BW",&myweight);
	}

	tree->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&bmass);
	tree->SetBranchAddress(TString::Format("BDT%d",bdtConf),&mybdt1);

	for(Int_t i=0; i<nEntries; i++)
	{
		tree->GetEntry(i);

		if(mybdt1 > mybdt)
		{
			if((flag == 2 && bmass > 5700) || flag == 1)
			{
				sum += myweight;
			}
		}
	}
	cout<<"mybdt = "<<mybdt<<" sum = "<<sum<<endl;
	return sum;
}
