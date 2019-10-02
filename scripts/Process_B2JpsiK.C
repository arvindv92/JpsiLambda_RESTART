#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TFileCollection.h"
#include "TCollection.h"
#include "TROOT.h"
#include <iostream>

using namespace std;

void Process_B2JpsiK(Int_t run = 1, Bool_t testing = false, Bool_t logFlag = true)
//Input files: Raw data from DaVinci
//Output file: rootFiles/dataFiles/B2JpsiK/run%d/jpsik.root
{
	// ROOT::EnableImplicitMT();
	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");

	if(logFlag) gSystem->RedirectOutput(Form("logs/data/B2JpsiK/run%d/Process_Cuts.txt",run),"w");

	TStopwatch sw;
	sw.Start();

	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"********************************"<<endl;
	cout<<"==> Starting Process_B2JpsiK: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"********************************"<<endl;

	TChain *h1 = new TChain("B2JpsiKTree/MyTuple");
	TString dir = "/data1/avenkate/B+toJpsiK+/RealData_AllKaons";

	if(run == 1)
	{
		gSystem->Exec(Form("ls %s/2011_MagDown/1158_*/jpsik.root > fileList/B2JpsiK_run1Files.txt",dir.Data()));
	}
	else if(run == 2)
	{
		gSystem->Exec(Form("ls %s/2016_MagDown/1170_*/jpsik.root > fileList/B2JpsiK_run2Files.txt",dir.Data()));
	}
	TFileCollection fc(Form("run%d",run),"",Form("fileList/B2JpsiK_run%dFiles.txt",run));
	TCollection *fc_list = (TCollection*)fc.GetList();

	if(testing)
	{
		h1->AddFileInfoList(fc_list,50);
	}
	else
	{
		h1->AddFileInfoList(fc_list);
	}

	cout<<"DONE ATTACHING ROOT FILES"<<endl;

	cout<<"I am making the following cuts. Sit tight"<<endl;
	TCut b_dtf_cut       = "B_DTF_VCHI2NDOF_JpsiConstr > 0 && B_DTF_VCHI2NDOF_JpsiConstr < 20";
	TCut b_minipchi2_cut = "B_MINIPCHI2 < 25";
	TCut b_trig_cut      = "(Jpsi_Hlt1DiMuonHighMassDecision_TOS==1||Jpsi_Hlt1TrackMuonDecision_TOS==1||Jpsi_Hlt1TrackAllL0Decision_TOS==1)&&(Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS==1)";
	TCut b_dira_cut      = "B_DIRA_OWNPV > 0.999";
	TCut b_mass_cut      = "B_DTF_M_JpsiConstr > 5200 && B_DTF_M_JpsiConstr < 5360";

	b_dtf_cut.Print();
	b_minipchi2_cut.Print();
	b_trig_cut.Print();
	b_dira_cut.Print();
	b_mass_cut.Print();

	TTree *treeIn = (TTree*)h1;

	TFile *fileOut = new TFile(Form("rootFiles/dataFiles/B2JpsiK/run%d/jpsik.root",run),"RECREATE");
	TTree *treeOut = (TTree*)treeIn->CloneTree(0);

	Long64_t entries_init = 0, entries_final = 0;

	entries_init = treeIn->GetEntries();

	Double_t B_DTF_VCHI2 = 0.0;
	Double_t B_MINIPCHI2 = 0.0;
	Bool_t Jpsi_Hlt1DiMuon = 0;
	Bool_t Jpsi_Hlt1TrackMuon = 0;
	Bool_t Jpsi_Hlt1TrackAllL0 = 0;
	Bool_t Jpsi_Hlt2DiMuon = 0;
	Double_t B_DIRA_OWNPV = 0.0;
	Double_t B_TAU = 0.0;
	Double_t B_DTF_M_JpsiConstr = 0.0;

	treeIn->SetBranchAddress("B_DTF_VCHI2NDOF_JpsiConstr",&B_DTF_VCHI2);
	treeIn->SetBranchAddress("B_MINIPCHI2",&B_MINIPCHI2);
	treeIn->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_TOS",&Jpsi_Hlt1DiMuon);
	treeIn->SetBranchAddress("Jpsi_Hlt1TrackMuonDecision_TOS",&Jpsi_Hlt1TrackMuon);
	treeIn->SetBranchAddress("Jpsi_Hlt1TrackAllL0Decision_TOS",&Jpsi_Hlt1TrackAllL0);
	treeIn->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS",&Jpsi_Hlt2DiMuon);
	treeIn->SetBranchAddress("B_DIRA_OWNPV",&B_DIRA_OWNPV);
	treeIn->SetBranchAddress("B_TAU",&B_TAU);
	treeIn->SetBranchAddress("B_DTF_M_JpsiConstr",&B_DTF_M_JpsiConstr);

	cout<<"Incoming Entries = "<<entries_init<<endl;

	cout<<"Starting Loop"<<endl;
	for (Long64_t i=0; i<entries_init; i++)
	{
		if(i%100000 == 0)
			cout<<i<<endl;

		treeIn->GetEntry(i);
		if(B_DTF_M_JpsiConstr > 5200 && B_DTF_M_JpsiConstr < 5360)//CHANGE: EARLIER MASS CUT WAS 5000 - 5600.
		{
			if(B_DTF_VCHI2 > 0 && B_DTF_VCHI2 < 20)
			{
				if(B_MINIPCHI2 < 25)
				{
					if((Jpsi_Hlt1DiMuon==1||Jpsi_Hlt1TrackMuon==1||Jpsi_Hlt1TrackAllL0==1)&&(Jpsi_Hlt2DiMuon==1))
					{
						if(B_DIRA_OWNPV > 0.999)
						{
							if(B_TAU > 0)
							{
								treeOut->Fill();
							}
						}
					}
				}
			}
		}
	}
	cout<<"Done with Loop"<<endl;
	entries_final = treeOut->GetEntries();
	cout<<"Outgoing entries = "<<entries_final<<endl;

	fileOut->Write();
	fileOut->Close();

	sw.Stop();

	cout<<"*****Done with Process_B2JpsiK*****"<<endl;
	sw.Print();

	if(logFlag) gSystem->RedirectOutput(0);

}
