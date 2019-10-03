#include "TFile.h"
#include "TTree.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include <iostream>

using namespace std;
void AliasFile(Int_t run = 1, Int_t mcType = 1, Int_t logFlag = true)
/*
   mcType = 1 for Lb -> J/psi Lambda MC
   mcType = 2 for Lb -> J/psi Sigma MC        (reco'd JpsiLambda)
   mcType = 3 for Xib -> Jpsi Xi MC           (reco'd JpsiLambda)
   mcType = 4 for Lb -> J/psi Lambda*(1405)MC (reco'd JpsiLambda)
   mcType = 5 for Lb -> J/psi Lambda*(1520)MC (reco'd JpsiLambda)
   mcType = 6 for Lb -> J/psi Lambda*(1600)MC (reco'd JpsiLambda)
   mcType = 7 forB0 -> Jpsi Ks MC             (reco'd JpsiLambda)
   mcType = 8 for Xib0 -> J/psi Lambda MC
   mcType = 9 for Xib- -> J/psi Xi- MC        (reco'd J/psi Xi-)
 */
{
	TStopwatch sw;
	sw.Start();

	const char *folder = "", *part = "";
	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part   = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part   = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part   = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part   = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Lst1405";
		part   = "lst1405";
		break;
	}
	case 5:
	{
		folder = "Lst1520";
		part   = "lst1520";
		break;
	}
	case 6:
	{
		folder = "Lst1600";
		part   = "lst1600";
		break;
	}
	case 7:
	{
		folder = "JpsiKs";
		part   = "jpsiks";
		break;
	}
	case 8:
	{
		folder = "Xib0";
		part   = "xib0";
		break;
	}
	case 9:
	{
		folder = "JpsiXi";
		part   = "jpsixi";
		break;
	}
	default:
	{
		cout<<"$$$ MC Type doesn't match any of the allowed cases. Exiting! $$$"<<endl;
		exit(1);
	}
	}

	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING");
	if(mcType != 9)
	{
		if(logFlag) gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/Alias.txt",folder,run),"w");
	}
	else
	{
		if(logFlag) gSystem->RedirectOutput(Form("logs/mc/JpsiXi/run%d/Alias.txt",run),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting AliasFile: "<<endl;
	cout<<"WD = "<<gSystem->pwd()<<endl;
	gSystem->Exec("date");
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr;
	TTree *treeIn = nullptr;

	if(mcType!=9)
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/%s/run%d/%s.root",folder,run,part),"UPDATE");
	else
		fileIn = TFile::Open(Form("rootFiles/mcFiles/JpsiXi/run%d/%s.root",run,part),"UPDATE");

	treeIn = (TTree*)fileIn->Get("MCTuple/MCDecayTree");

	if(!(mcType == 3 || mcType == 9 || mcType == 8))
	{
		treeIn->SetAlias("Lb_PT","Lambda_b0_TRUEPT");
		treeIn->SetAlias("Lb_P","sqrt(Lambda_b0_TRUEPT**2+Lambda_b0_TRUEP_Z**2)");
		treeIn->SetAlias("Lb_ETA","-log(tan(0.5*atan(Lambda_b0_TRUEPT/Lambda_b0_TRUEP_Z)))");
	}
	else if(mcType == 3 || mcType == 9) //Xib -> J/psi Xi (reco'd J/psi Lambda & J/psi Xi)
	{
		treeIn->SetAlias("Lb_PT","Xi_bminus_TRUEPT");
		treeIn->SetAlias("Lb_P","sqrt(Xi_bminus_TRUEPT**2+Xi_bminus_TRUEP_Z**2)");
		treeIn->SetAlias("Lb_ETA","-log(tan(0.5*atan(Xi_bminus_TRUEPT/Xi_bminus_TRUEP_Z)))");
	}
	else if(mcType == 8) // Xib0 -> J/psi Lambda
	{
		treeIn->SetAlias("Lb_PT","Xi_b0_TRUEPT");
		treeIn->SetAlias("Lb_P","sqrt(Xi_b0_TRUEPT**2+Xi_b0_TRUEP_Z**2)");
		treeIn->SetAlias("Lb_ETA","-log(tan(0.5*atan(Xi_b0_TRUEPT/Xi_b0_TRUEP_Z)))");
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

	fileIn->cd("MCTuple");
	treeIn->Write("",TObject::kOverwrite);

	fileIn->Close();

	sw.Stop();
	cout << "==> AliasFile is done! Huzzah!: "; sw.Print();
	if(logFlag) gSystem->RedirectOutput(0);
}
