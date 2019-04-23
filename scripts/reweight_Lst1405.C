#include "MVmodel.C"
#include "BONN.C"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TEntryList.h"
#include "iostream"

using namespace std;
void reweight_Lst1405(Int_t run = 1)
{
	TFile *filein, *mcfile_MV, *mcfile_BONN;
	TTree *treein, *treein_gen, *treeout_MV, *treeout_BONN;

	Double_t Lst_PE, Lst_PX, Lst_PY, Lst_PZ;
	Double_t Jpsi_PE, Jpsi_PX, Jpsi_PY, Jpsi_PZ;
	Double_t Lambda_PE, Lambda_PX, Lambda_PY,Lambda_PZ;
	Double_t mass_gen, JpsiLmass, MVweight, BONNweight, num, den;

	Int_t nentries, massbin_MV, massbin_BONN;
	Int_t low_MV, high_MV, low_BONN, high_BONN;//limits for Lst mass in which the models work
	Int_t nbins_MV, nbins_BONN, bkgcat;

	ULong64_t evtno_rec, evtno_gen;
	UInt_t runno_rec, runno_gen;

	low_MV = 1300;//these values are rough estimates. Look at paper and get precise limits
	high_MV = 2523;
	low_BONN = 1328;
	high_BONN = 2519;

	Int_t binwidth = 1;

	nbins_MV = (Int_t)(high_MV - low_MV)/binwidth;
	nbins_BONN = (Int_t)(high_BONN - low_BONN)/binwidth;

	TH1D *theory_MVmodel = new TH1D("theory_MVmodel","",nbins_MV,low_MV,high_MV);
	TH1D *theory_BONNmodel = new TH1D("theory_BONNmodel","",nbins_BONN,low_BONN,high_BONN);

	// gSystem->Exec("cp ./Lst1405_total.root ./Lst1405_total_MVrw.root");
	// gSystem->Exec("cp ./Lst1405_total.root ./Lst1405_total_BONNrw.root");

	// mcfile = TFile::Open("./Lst1405_total_MVrw.root","UPDATE");
	// treein = (TTree*)mcfile->Get("MCDecayTreeTuple/MCDecayTree");

	filein = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/lst1405.root",run),"READ");
	treein = (TTree*)filein->Get("Lb2JpsiLTree/MyTuple");
	treein_gen = (TTree*)filein->Get("MCTuple/MCDecayTree");

	cout<<"Copying Trees"<<endl;
	mcfile_MV = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/lst1405_MV.root",run),"RECREATE");
	// treeout_MV = (TTree*)treein->CloneTree(0);
	treeout_MV = new TTree("MyTuple","");
	mcfile_BONN = new TFile(Form("../rootFiles/mcFiles/JpsiLambda/Lst1405/run%d/lst1405_BONN.root",run),"RECREATE");
	// treeout_BONN = (TTree*)treein->CloneTree(0);
	treeout_BONN = new TTree("MyTuple","");

	cout<<"Done copying Trees"<<endl;

	treein_gen->Draw(Form("sqrt(pow(Lambda_1405_0_TRUEP_E,2) - pow(Lambda_1405_0_TRUEP_X,2) - pow(Lambda_1405_0_TRUEP_Y,2) - pow(Lambda_1405_0_TRUEP_Z,2) )>>hmv(%d,%d,%d)",nbins_MV,low_MV,high_MV),"","goff");
	treein_gen->Draw(Form("sqrt(pow(Lambda_1405_0_TRUEP_E,2) - pow(Lambda_1405_0_TRUEP_X,2) - pow(Lambda_1405_0_TRUEP_Y,2) - pow(Lambda_1405_0_TRUEP_Z,2) )>>hbonn(%d,%d,%d)",nbins_BONN,low_BONN,high_BONN),"","goff");

	TH1D *hmv = (TH1D*)gDirectory->Get("hmv");
	hmv->Scale(1.0/hmv->Integral());

	TH1D *hbonn = (TH1D*)gDirectory->Get("hbonn");
	hbonn->Scale(1.0/hbonn->Integral());

	for(Int_t i=0; i<nbins_MV; i++)
	{
		theory_MVmodel->Fill((i+low_MV)*1.0,(Double_t)MVmodel(i+low_MV));
	}
	for(Int_t i=0; i<nbins_BONN; i++)
	{
		theory_BONNmodel->Fill((i+low_BONN)*1.0,(Double_t)invariantmassdistribution_sig0pi0((i+low_BONN)*1.0/1000));//takes input in GeV
	}

	theory_MVmodel->Scale(1.0/theory_MVmodel->Integral());
	theory_BONNmodel->Scale(1.0/theory_BONNmodel->Integral());

	theory_MVmodel->SetLineColor(kMagenta+2);
	hmv->SetLineColor(kRed);
	hbonn->SetLineColor(kRed);

	theory_BONNmodel->Draw("HIST");
	theory_MVmodel->Draw("sameHIST");
	hmv->Draw("same");

	TH1D *wts_MV = (TH1D*)theory_MVmodel->Clone("wts_MV");
	TH1D *wts_BONN = (TH1D*)theory_BONNmodel->Clone("wts_BONN");

	wts_MV->Divide(hmv);
	wts_BONN->Divide(hbonn);

	new TCanvas();
	wts_MV->Draw();

	new TCanvas();
	wts_BONN->Draw();

	nentries = treein->GetEntries();

	treein_gen->SetBranchAddress("Lambda_1405_0_TRUEP_E",&Lst_PE);
	treein_gen->SetBranchAddress("Lambda_1405_0_TRUEP_X",&Lst_PX);
	treein_gen->SetBranchAddress("Lambda_1405_0_TRUEP_Y",&Lst_PY);
	treein_gen->SetBranchAddress("Lambda_1405_0_TRUEP_Z",&Lst_PZ);

	treein_gen->SetBranchAddress("J_psi_1S_TRUEP_E",&Jpsi_PE);
	treein_gen->SetBranchAddress("J_psi_1S_TRUEP_X",&Jpsi_PX);
	treein_gen->SetBranchAddress("J_psi_1S_TRUEP_Y",&Jpsi_PY);
	treein_gen->SetBranchAddress("J_psi_1S_TRUEP_Z",&Jpsi_PZ);

	treein_gen->SetBranchAddress("Lambda0_TRUEP_E",&Lambda_PE);
	treein_gen->SetBranchAddress("Lambda0_TRUEP_X",&Lambda_PX);
	treein_gen->SetBranchAddress("Lambda0_TRUEP_Y",&Lambda_PY);
	treein_gen->SetBranchAddress("Lambda0_TRUEP_Z",&Lambda_PZ);

	treein_gen->SetBranchAddress("eventNumber",&evtno_gen);
	treein_gen->SetBranchAddress("runNumber",&runno_gen);

	treein->SetBranchAddress("eventNumber",&evtno_rec);
	treein->SetBranchAddress("runNumber",&runno_rec);
	treein->SetBranchAddress("Lb_BKGCAT",&bkgcat);

	TBranch *MVwtbranch = treeout_MV->Branch("MVweight",&MVweight,"MVweight/D");
	TBranch *BONNwtbranch = treeout_BONN->Branch("BONNweight",&BONNweight,"BONNweight/D");
	// TBranch *JpsiLMass_MV = treeout_MV->Branch("JpsiLmass",&JpsiLmass,"JpsiLmass/D");
	// TBranch *JpsiLMass_BONN = treeout_BONN->Branch("JpsiLmass",&JpsiLmass,"JpsiLmass/D");

	Int_t ctr = 0;
	for(Int_t i=0; i<nentries; i++)
	{
		if(i==100) break;

		if(i%10000 == 0)
			cout<<i<<endl;

		treein->GetEntry(i);
		cout<<"bkgcat = "<<bkgcat<<endl;
		if(bkgcat==50)
		{
			// treein_gen->Draw("sqrt(pow(Lambda_1405_0_TRUEP_E,2) - pow(Lambda_1405_0_TRUEP_X,2) - pow(Lambda_1405_0_TRUEP_Y,2) - pow(Lambda_1405_0_TRUEP_Z,2) )>>h1",Form("eventNumber==%llu",evtno_rec),"goff");
			// TH1F *h1 = (TH1F*)gDirectory->Get("h1");

			cout<<"poop1"<<endl;
			treein_gen->Draw(">>elist",Form("eventNumber == %llu && runNumber == %u",evtno_rec,runno_rec),"entrylist");
			cout<<"poop2"<<endl;
			TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
			cout<<"poop3"<<endl;
			treein_gen->SetEntryList(elist);

			cout<<"list N = "<<elist->GetN()<<endl;

			Int_t entryNum = treein_gen->GetEntryNumber(0);
			cout<<"entryNum = "<<entryNum<<endl;

			treein_gen->GetEntry(entryNum);

			cout<<"evtno_rec = "<<evtno_rec<<" evtno_gen = "<<evtno_gen<<endl;
			cout<<"runno_rec = "<<runno_rec<<" runno_gen = "<<runno_gen<<endl;
			// mass_gen = h1->GetBinCenter(51);
			mass_gen = sqrt(Lst_PE*Lst_PE - Lst_PX*Lst_PX - Lst_PY*Lst_PY - Lst_PZ*Lst_PZ);
			cout<<"mass_gen = "<<mass_gen<<endl;
			// JpsiLmass = sqrt(pow(Jpsi_PE+Lambda_PE,2) - pow(Jpsi_PX+Lambda_PX,2) - pow(Jpsi_PY+Lambda_PY,2) - pow(Jpsi_PZ+Lambda_PZ,2));

			massbin_MV = wts_MV->FindBin(mass_gen);
			massbin_BONN = wts_BONN->FindBin(mass_gen);
			MVweight = wts_MV->GetBinContent(massbin_MV);
			BONNweight = wts_BONN->GetBinContent(massbin_BONN);

			ctr++;
		}
		else
		{
			MVweight = 1.0;
			BONNweight = 1.0;
		}
		MVwtbranch->Fill();
		BONNwtbranch->Fill();
		// JpsiLMass_MV->Fill();
	}
	cout<<"ctr = "<<ctr<<endl;
	mcfile_MV->cd();
	treeout_MV->Write("",TObject::kOverwrite);
	cout<<"MV DONE"<<endl;

	mcfile_BONN->cd();
	treeout_BONN->Write("",TObject::kOverwrite);
	cout<<"BONN DONE"<<endl;
}
