#include "TFile.h"
#include "TTree.h"

void AliasFile(Int_t run = 1, Int_t mcType = 1)
//mcType = 1 for JpsiLambda MC
//mcType = 2 for JpsiSigma MC
//mcType = 6 for JpsiLst(1405) MC
//mcType = 7 for JpsiLst(1520) MC
//mcType = 8 for JpsiLst(1600) MC
//mcType = 9 for chiC1 Lambda MC
{
	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;

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
	case 6:
	{
		folder = "Lst1405";
		part = "lst1405";
		break;
	}
	case 7:
	{
		folder = "Lst1520";
		part = "lst152-";
		break;
	}
	case 8:
	{
		folder = "Lst1600";
		part = "lst1600";
		break;
	}
	case 9:
	{
		folder = "chiC1";
		part = "chic1";
		break;
	}
	}
	fileIn = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/%s.root",folder,run,part),"READ");
	treeIn = (TTree*)fileIn->Get("MCTuple/MCDecayTree");

	if(mcType != 3)
	{
		treeIn->SetAlias("Lb_PT","Lambda_b0_TRUEPT");
		treeIn->SetAlias("Lb_P","sqrt(Lambda_b0_TRUEPT**2+Lambda_b0_TRUEP_Z**2)");
		treeIn->SetAlias("Lb_ETA","-log(tan(0.5*atan(Lambda_b0_TRUEPT/Lambda_b0_TRUEP_Z)))");
	}
	else if(mcType == 3)
	{
		treeIn->SetAlias("Lb_PT","Xi_bminus_TRUEPT");
		treeIn->SetAlias("Lb_P","sqrt(Xi_bminus_TRUEPT**2+Xi_bminus_TRUEP_Z**2)");
		treeIn->SetAlias("Lb_ETA","-log(tan(0.5*atan(Xi_bminus_TRUEPT/Xi_bminus_TRUEP_Z)))");
	}

	treeIn->SetAlias("Jpsi_PT","J_psi_1S_TRUEPT");
	treeIn->SetAlias("L_PT","Lambda0_TRUEPT");
	treeIn->SetAlias("p_PT","pplus_TRUEPT");
	treeIn->SetAlias("pi_PT","piminus_TRUEPT");

	treeIn->SetAlias("Jpsi_P","sqrt(J_psi_1S_TRUEPT**2+J_psi_1S_TRUEP_Z**2)");
	treeIn->SetAlias("L_P","sqrt(Lambda0_TRUEPT**2+Lambda0_TRUEP_Z**2)");
	treeIn->SetAlias("p_P","sqrt(pplus_TRUEPT**2+pplus_TRUEP_Z**2)");
	treeIn->SetAlias("pi_P","sqrt(piminus_TRUEPT**2+piminus_TRUEP_Z**2)");

	treeIn->SetAlias("Jpsi_ETA","-log(tan(0.5*atan(J_psi_1S_TRUEPT/J_psi_1S_TRUEP_Z)))");
	treeIn->SetAlias("L_ETA","-log(tan(0.5*atan(Lambda0_TRUEPT/Lambda0_TRUEP_Z)))");
	treeIn->SetAlias("p_ETA","-log(tan(0.5*atan(pplus_TRUEPT/pplus_TRUEP_Z)))");
	treeIn->SetAlias("pi_ETA","-log(tan(0.5*atan(piminus_TRUEPT/piminus_TRUEP_Z)))");


	fileOut = TFile::Open(Form("../rootFiles/mcFiles/JpsiLambda/%s/run%d/RW/%s_aliased.root",folder,run,part),"RECREATE");
	treeOut = (TTree*)treeIn->CopyTree("");

	fileOut->Write();
	fileOut->Close();

	fileIn->Close();
}
