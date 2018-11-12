#include <TFile.h>
#include <TTree.h>
#include <iostream>

using namespace std;
void makeFinalBDTTuples(Int_t flag = 1)//flag = 1 for signal training files, flag = 2 for background training files, flag = 3 for total file on which BDT is applied
{
	TFile *filein_LL, *filein_DD, *fileout_LL, *fileout_DD;
	TTree *treein_LL, *treein_DD, *treeout_LL, *treeout_DD;

	if(flag == 1 || flag == 2) {
		if(flag == 1)
			cout<<"Make signal training files for final BDT"<<endl;
		if(flag == 2)
			cout<<"Make background training files for final BDT"<<endl;

		filein_LL = TFile::Open("trainingfiles/jpsilambda_LL_forTMVAisoTraining_temp.root","READ");
		treein_LL = (TTree*)filein_LL->Get("MyTuple");

		filein_DD = TFile::Open("trainingfiles/jpsilambda_DD_forTMVAisoTraining_temp.root","READ");
		treein_DD = (TTree*)filein_DD->Get("MyTuple");
	}
	else if(flag == 3) {
		cout<<"Make TMVA application files for final BDT"<<endl;

		filein_LL = TFile::Open("../isolation/output/jpsilambda_LL_withiso.root","READ");
		treein_LL = (TTree*)filein_LL->Get("MyTuple");

		filein_DD = TFile::Open("../isolation/output/jpsilambda_DD_withiso.root","READ");
		treein_DD = (TTree*)filein_DD->Get("MyTuple");
	}

	treein_LL->SetBranchStatus("*",0);
	treein_LL->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	treein_LL->SetBranchStatus("Lb_ConsLb_chi2",1);
	treein_LL->SetBranchStatus("Lb_MINIPCHI2",1);
	treein_LL->SetBranchStatus("Lb_DIRA_OWNPV",1);
	treein_LL->SetBranchStatus("Lb_FD_OWNPV",1);
	treein_LL->SetBranchStatus("Lb_DTF_CTAUS_L",1);
	treein_LL->SetBranchStatus("Lb_DTF_CTAU_L",1);

	treein_LL->SetBranchStatus("Jpsi_MINIPCHI2",1);
	treein_LL->SetBranchStatus("Jpsi_M",1);
	treein_LL->SetBranchStatus("Jpsi_CosTheta",1);
	treein_LL->SetBranchStatus("Jpsi_PT",1);

	treein_LL->SetBranchStatus("L_FDCHI2_ORIVX",1);
	treein_LL->SetBranchStatus("L_DIRA_ORIVX",1);
	treein_LL->SetBranchStatus("L_FD_ORIVX",1);
	treein_LL->SetBranchStatus("L_DIRA_OWNPV",1);
	treein_LL->SetBranchStatus("L_dm",1);
	treein_LL->SetBranchStatus("L_MINIPCHI2",1);
	treein_LL->SetBranchStatus("L_PT",1);
	treein_LL->SetBranchStatus("L_ENDVERTEX_CHI2",1);
	treein_LL->SetBranchStatus("L_CosTheta",1);

	treein_LL->SetBranchStatus("muplus_MINIPCHI2",1);
	treein_LL->SetBranchStatus("muminus_MINIPCHI2",1);

	treein_LL->SetBranchStatus("p_PIDp",1);
	treein_LL->SetBranchStatus("p_MINIPCHI2",1);
	treein_LL->SetBranchStatus("p_TRACK_GhostProb",1);
	treein_LL->SetBranchStatus("p_PT",1);
	treein_LL->SetBranchStatus("p_ProbNNp",1);

	treein_LL->SetBranchStatus("pi_PIDK",1);
	treein_LL->SetBranchStatus("pi_MINIPCHI2",1);
	treein_LL->SetBranchStatus("pi_TRACK_GhostProb",1);
	treein_LL->SetBranchStatus("pi_PT",1);
	treein_LL->SetBranchStatus("pi_ProbNNpi",1);

	treein_LL->SetBranchStatus("BDTk",1);

	if(flag!=3) {
		treein_LL->SetBranchStatus("SW",1);
		treein_LL->SetBranchStatus("BW",1);
	}
	treein_DD->SetBranchStatus("*",0);
	treein_DD->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	treein_DD->SetBranchStatus("Lb_ConsLb_chi2",1);
	treein_DD->SetBranchStatus("Lb_MINIPCHI2",1);
	treein_DD->SetBranchStatus("Lb_DIRA_OWNPV",1);
	treein_DD->SetBranchStatus("Lb_FD_OWNPV",1);
	treein_DD->SetBranchStatus("Lb_DTF_CTAUS_L",1);
	treein_DD->SetBranchStatus("Lb_DTF_CTAU_L",1);

	treein_DD->SetBranchStatus("Jpsi_MINIPCHI2",1);
	treein_DD->SetBranchStatus("Jpsi_M",1);
	treein_DD->SetBranchStatus("Jpsi_CosTheta",1);
	treein_DD->SetBranchStatus("Jpsi_PT",1);

	treein_DD->SetBranchStatus("L_FDCHI2_ORIVX",1);
	treein_DD->SetBranchStatus("L_DIRA_ORIVX",1);
	treein_DD->SetBranchStatus("L_FD_ORIVX",1);
	treein_DD->SetBranchStatus("L_DIRA_OWNPV",1);
	treein_DD->SetBranchStatus("L_dm",1);
	treein_DD->SetBranchStatus("L_MINIPCHI2",1);
	treein_DD->SetBranchStatus("L_PT",1);
	treein_DD->SetBranchStatus("L_ENDVERTEX_CHI2",1);
	treein_DD->SetBranchStatus("L_CosTheta",1);

	treein_DD->SetBranchStatus("muplus_MINIPCHI2",1);
	treein_DD->SetBranchStatus("muminus_MINIPCHI2",1);

	treein_DD->SetBranchStatus("p_PIDp",1);
	treein_DD->SetBranchStatus("p_MINIPCHI2",1);
	treein_DD->SetBranchStatus("p_TRACK_GhostProb",1);
	treein_DD->SetBranchStatus("p_PT",1);
	treein_DD->SetBranchStatus("p_ProbNNp",1);

	treein_DD->SetBranchStatus("pi_PIDK",1);
	treein_DD->SetBranchStatus("pi_MINIPCHI2",1);
	treein_DD->SetBranchStatus("pi_TRACK_GhostProb",1);
	treein_DD->SetBranchStatus("pi_PT",1);
	treein_DD->SetBranchStatus("pi_ProbNNpi",1);

	treein_DD->SetBranchStatus("BDTk",1);

	if(flag!=3) {
		treein_DD->SetBranchStatus("SW",1);
		treein_DD->SetBranchStatus("BW",1);
	}

	if(flag==2) {
		fileout_LL = new TFile("trainingfiles/jpsilambda_LL_bkg.root","recreate");
		treeout_LL = (TTree*)treein_LL->CopyTree("Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000 && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 7000");

		fileout_DD = new TFile("trainingfiles/jpsilambda_DD_bkg.root","recreate");
		treeout_DD = (TTree*)treein_DD->CopyTree("Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000 && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 7000");
	}

	if(flag==1) {
		fileout_LL = new TFile("trainingfiles/jpsilambda_LL_sig.root","recreate");
		treeout_LL = (TTree*)treein_LL->CopyTree("Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000");

		fileout_DD = new TFile("trainingfiles/jpsilambda_DD_sig.root","recreate");
		treeout_DD = (TTree*)treein_DD->CopyTree("Lb_MINIPCHI2 < 5000 && Jpsi_MINIPCHI2 < 5000 && pi_PIDK < 5 && p_PIDp > 5 && pi_MINIPCHI2 < 5000 && L_MINIPCHI2 < 10000");
	}

	else if(flag==3) {

		fileout_LL = new TFile("applicationfiles/jpsilambda_LL_forfinalBDTApp_etacut.root","recreate");
		treeout_LL = (TTree*)treein_LL->CopyTree("");

		fileout_DD = new TFile("applicationfiles/jpsilambda_DD_forfinalBDTApp_etacut.root","recreate");
		treeout_DD = (TTree*)treein_DD->CopyTree("");
	}

	fileout_LL->Write();

	fileout_DD->Write();

}
