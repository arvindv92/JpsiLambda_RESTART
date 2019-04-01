#include "Fitscript_simul.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

#define Open TFile::Open
// std::vector <Float_t> Fitscript_simul(Int_t bdtConf_nonZero, Int_t bdtConf_Zero,
//                                       Int_t isoConf[i], const char *isoVersion,
//                                       Float_t bdtCut_nonZero, Float_t bdtCut_Zero[i])
void Fitscript_simul(const char *option, Int_t myLow, Int_t myHigh)
{

	gSystem->RedirectOutput(Form("../logs/data/JpsiLambda/UpperLimit/Fit_HypatiaSig_ExpBkg_%d_%d.txt",
	                             myLow,myHigh),"w");

	gSystem->Exec("date");

	gSystem->Load("RooHypatia2_cpp.so");

	Float_t bdtCut_nonZero[2] = {0.0,0.0};

	Float_t bdtCut_Zero[2] = {0.0,0.0};

	Int_t bdtConf_nonZero[2] = {0,0};

	Int_t bdtConf_Zero[2] = {0,0};

	Int_t isoConf[2] = {0,0};

	Bool_t calcUL = true;

	const char *isoVersion[2] = {"",""};

	if(!strncmp(option,"best",4))
	{
		isoVersion[0] = "v1";
		isoVersion[1] = "v0";

		isoConf[0] = 1;
		isoConf[1] = 2;

		bdtConf_nonZero[0] = 2;
		bdtConf_nonZero[1] = 1;

		bdtConf_Zero[0] = 1;
		bdtConf_Zero[1] = 1;

		bdtCut_nonZero[0] = 0.375;
		bdtCut_nonZero[1] = 0.535;

		bdtCut_Zero[0] = 0.285;
		bdtCut_Zero[1] = 0.415;
	}

	std::vector <Float_t> yields;

	Double_t xibnorm_LL[2]         = {0.0,0.0};
	Double_t xibnorm_LL_staterr[2] = {0.0,0.0};
	Double_t xibnorm_LL_systerr[2] = {0.0,0.0};
	Double_t xibnorm_LL_err[2]     = {0.0,0.0};

	Double_t eff_Lambda_rec[2]     = {0.0,0.0};
	Double_t eff_Sigma_rec[2]     = {0.0,0.0};
	Double_t eff_Lambda_rec_err[2] = {0.0,0.0};
	Double_t eff_Sigma_rec_err[2]  = {0.0,0.0};

	Double_t eff_Lambda_gen[2]     = {0.0,0.0};
	Double_t eff_Sigma_gen[2]     = {0.0,0.0};
	Double_t eff_Lambda_gen_err[2] = {0.0,0.0};
	Double_t eff_Sigma_gen_err[2]  = {0.0,0.0};

	Double_t eff_Lambda[2] = {0.0,0.0};
	Double_t eff_Lambda_err[2] = {0.0,0.0};

	Double_t eff_Sigma[2] = {0.0,0.0};
	Double_t eff_Sigma_err[2] = {0.0,0.0};

	Int_t nGen_Lambda[2]       = {0,0};
	Int_t nGen_Sigma[2]        = {0,0};

	Int_t lst1405flag = 1;
	Int_t lst1520flag = 1;
	Int_t lst1810flag = 0;
	Int_t xibflag     = 1;
	Int_t sigmaflag   = 1;
	Int_t Lst1405_rwtype = 2;

	Int_t low      = myLow, high = myHigh; //Define range in which fit is performed
	Int_t binwidth = 4; //Final fit is a binned fit with this binning
	Int_t nbins    = (Int_t)(high-low)/binwidth;

	TFile *filein_lst1405 = nullptr;
	TTree *treein_lst1405 = (TTree*)malloc(sizeof(*treein_lst1405));
	TString lst1405wtexp = "";

	RooWorkspace w("w");

	//*********MASTER VARIABLE******************************************
	w.factory(Form("bMass[%d,%d]",low,high));
	w.var("bMass")->setBins((Int_t)(high-low)/binwidth);
	//***********************E******************************************

	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	//*****Run script to get Xib normalization*********
	for(Int_t run = 1; run<=2; run++) {
		Int_t i = run-1;
		gSystem->Exec(Form("python -c \'from GetXibNorm import GetNorm;"
		                   " GetNorm(%d, \"%s\", %d, %d, %d, %f, %f)\'",
		                   run,isoVersion[i],isoConf[i],bdtConf_nonZero[i],
		                   bdtConf_Zero[i],bdtCut_nonZero[i],bdtCut_Zero[i]));

		ifstream infile(Form("../logs/mc/JpsiXi/run%d/xibNorm_log.txt",run));

		infile>>xibnorm_LL[i];
		infile>>xibnorm_LL_staterr[i];

		xibnorm_LL_systerr[i] = xibnorm_LL[i]*0.0; //Assigning 2% systematic error for now

		xibnorm_LL_err[i] = sqrt(pow(xibnorm_LL_staterr[i],2) + pow(xibnorm_LL_systerr[i],2));//Combining stat and syst in quadrature


		xibnorm_LL[i] = xibnorm_LL[i]*2; //ACCOUNT FOR XIB0
		xibnorm_LL_err[i] = xibnorm_LL_err[i] * 1.414;//ACCOUNT FOR XIB0

		cout<<"The LL Xib normalization is "<<xibnorm_LL[i]<<" +/- "
		    <<xibnorm_LL_err[i]<<endl;
	}
	//************************************************
	//****Get J/psi Lambda efficiencies from MC*******
	const char* lambdaMCPath = "../rootFiles/mcFiles/JpsiLambda/JpsiLambda";

	for(Int_t run = 1; run<=2; run++) {
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_Lambda = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_nonZeroTracks.root",
		                                           lambdaMCPath,run));
		TTree *mcTreeIn_nonZero_Lambda = (TTree*)mcFileIn_nonZero_Lambda->Get("MyTuple");

		TFile *mcFileIn_Zero_Lambda    = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_ZeroTracks.root",
		                                           lambdaMCPath,run));
		TTree *mcTreeIn_Zero_Lambda    = (TTree*)mcFileIn_Zero_Lambda->Get("MyTuple");

		mcTreeIn_nonZero_Lambda->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_LL_FinalBDT%d_iso%d_%s.root",
		                                                  lambdaMCPath,run,bdtConf_nonZero[i],
		                                                  isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_Lambda->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_zeroTracksLL_FinalBDT%d.root",
		                                               lambdaMCPath,run,bdtConf_Zero[i]));
		fstream genFile_Lambda;
		genFile_Lambda.open((Form("../logs/mc/JpsiLambda/JpsiLambda/run%d/gen_log.txt",run)));

		genFile_Lambda>>nGen_Lambda[i];

		fstream genEffFile_Lambda;
		genEffFile_Lambda.open(Form("../logs/mc/JpsiLambda/JpsiLambda/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_Lambda>>eff_Lambda_gen[i];
		genEffFile_Lambda>>eff_Lambda_gen_err[i];

		cout<<"Lambda Generator Effs = "<<eff_Lambda_gen[i]*100<<" % +/- "<<eff_Lambda_gen_err[i]*100<<" %"<<endl;

		Int_t num_Lambda = mcTreeIn_nonZero_Lambda->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i]))
		                   + mcTreeIn_Zero_Lambda->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));

		eff_Lambda_rec[i]     = num_Lambda*1.0/nGen_Lambda[i];
		// eff_Lambda_rec_err = eff_Lambda_rec*0.02;

		cout<<"Lambda Recons. Effs = "<<eff_Lambda_rec[i]*100<<" %"<<endl;

		eff_Lambda[i] = eff_Lambda_rec[i]*eff_Lambda_gen[i];
		eff_Lambda_err[i] = eff_Lambda[i]*0.02;// nominally assigning 2% uncertainty for now
		// eff_Lambda_err = sqrt(pow(eff_Lambda_gen_err,2) + pow(eff_Lambda_rec_err,2)); //Combining errors in quadrature

		cout<<"Jpsi Lambda Eff = "<<eff_Lambda[i]*100<<" % +/- "<<eff_Lambda_err[i]*100<<" %"<<endl;
	}
	//************************************************
	//****Get J/psi Sigma efficiencies from MC*******
	const char *sigmaMCPath = "../rootFiles/mcFiles/JpsiLambda/JpsiSigma";

	for(Int_t run = 1; run<=2; run++) {
		Int_t i = run-1;
		TFile *mcFileIn_nonZero_Sigma = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_nonZeroTracks.root",
		                                          sigmaMCPath,run));
		TTree *mcTreeIn_nonZero_Sigma = (TTree*)mcFileIn_nonZero_Sigma->Get("MyTuple");

		TFile *mcFileIn_Zero_Sigma    = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_ZeroTracks.root",
		                                          sigmaMCPath,run));
		TTree *mcTreeIn_Zero_Sigma    = (TTree*)mcFileIn_Zero_Sigma->Get("MyTuple");

		mcTreeIn_nonZero_Sigma->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_LL_FinalBDT%d_iso%d_%s.root",
		                                                 sigmaMCPath,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		mcTreeIn_Zero_Sigma->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_zeroTracksLL_FinalBDT%d.root",
		                                              sigmaMCPath,run,bdtConf_Zero[i]));
		fstream genFile_Sigma;
		genFile_Sigma.open((Form("../logs/mc/JpsiLambda/JpsiSigma/run%d/gen_log.txt",run)));

		genFile_Sigma>>nGen_Sigma[i];

		fstream genEffFile_Sigma;
		genEffFile_Sigma.open(Form("../logs/mc/JpsiLambda/JpsiSigma/run%d/Generator_Effs_Combined.txt",run));

		genEffFile_Sigma>>eff_Sigma_gen[i];
		genEffFile_Sigma>>eff_Sigma_gen_err[i];

		cout<<"Sigma Generator Effs = "<<eff_Sigma_gen[i]*100<<" % +/- "<<eff_Sigma_gen_err[i]*100<<" %"<<endl;

		Int_t num_Sigma = mcTreeIn_nonZero_Sigma->GetEntries(Form("BDT%d > %f", bdtConf_nonZero[i],bdtCut_nonZero[i])) +
		                  mcTreeIn_Zero_Sigma->GetEntries(Form("BDT%d > %f", bdtConf_Zero[i],bdtCut_Zero[i]));

		eff_Sigma_rec[i]     = num_Sigma*1.0/nGen_Sigma[i];
		// eff_Sigma_rec_err = eff_Sigma_rec*0.02;

		cout<<"Sigma Recons. Effs = "<<eff_Sigma_rec[i]*100<<" %"<<endl;

		eff_Sigma[i] = eff_Sigma_gen[i] * eff_Sigma_rec[i];
		eff_Sigma_err[i] = eff_Sigma[i]* 0.02;// nominally assigning 2% uncertainty for now
		// eff_Sigma_err = sqrt(pow(eff_Sigma_gen_err[i],2) + pow(eff_Sigma_rec_err,2)); //Combining errors in quadrature

		cout<<"Jpsi Sigma Eff = "<<eff_Sigma[i]*100<<" % +/- "<<eff_Sigma_err[i]*100<<" %"<<endl;
	}
	//*******************************************************************

	//*********Get shape from Lambda*(1405) background*******************
	//  RooRealVar lstmass("JpsiLmass","JpsiL mass",4200.,5500.);
	//  RooRealVar MVweight("MVweight","MVweight",0.,10.);

	if(Lst1405_rwtype == 0)
	{
		filein_lst1405 = Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1405/Lst1405_total_MVrw.root");
		treein_lst1405 = (TTree*)filein_lst1405->Get("MCDecayTree");
		treein_lst1405->SetBranchStatus("*",0);
		treein_lst1405->SetBranchStatus("JpsiLmass",1);
		lst1405wtexp = "";
	}
	else if(Lst1405_rwtype == 1)
	{
		filein_lst1405 = Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1405/Lst1405_total_MVrw.root");
		lst1405wtexp = "MVweight";
		treein_lst1405 = (TTree*)filein_lst1405->Get("MCDecayTree");
		treein_lst1405->SetBranchStatus("*",0);
		treein_lst1405->SetBranchStatus("JpsiLmass",1);
		treein_lst1405->SetBranchStatus("MVweight",1);
	}
	else if(Lst1405_rwtype == 2)
	{
		filein_lst1405 = Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1405/Lst1405_total_BONNrw.root");
		lst1405wtexp = "BONNweight";
		treein_lst1405 = (TTree*)filein_lst1405->Get("MCDecayTree");
		treein_lst1405->SetBranchStatus("*",0);
		treein_lst1405->SetBranchStatus("JpsiLmass",1);
		treein_lst1405->SetBranchStatus("BONNweight",1);
	}

	treein_lst1405->Draw(Form("JpsiLmass>>hlst1405(%d,%d,%d)",nbins,low,high),lst1405wtexp,"goff");

	TH1D *hlst1405 = (TH1D*)gDirectory->Get("hlst1405");
	TH1D *hlst1405_smooth = (TH1D*)hlst1405->Clone("hlst1405_smooth");
	hlst1405_smooth->Smooth(2);//TODO TWEAK THIS
	RooDataHist *ds_lst1405 = new RooDataHist("ds_lst1405","ds_lst1405",
	                                          *(w.var("bMass")),hlst1405);
	RooDataHist *ds_lst1405_smooth = new RooDataHist("ds_lst1405_smooth","ds_lst1405_smooth",
	                                                 *(w.var("bMass")),hlst1405_smooth);

	ds_lst1405->Print();
	//  RooKeysPdf lstshape("lstshape","lstshape",lstmass,ds_lst);
	RooHistPdf lst1405shape("lst1405shape","lst1405shape",*(w.var("bMass")),
	                        *ds_lst1405,0);
	RooHistPdf LST1405("LST1405","LST1405",
	                   *(w.var("bMass")),*ds_lst1405_smooth,0);
	RooPlot *framelst1405 = (w.var("bMass"))->frame();
	framelst1405->SetTitle("J/#psi #Lambda^{*} (1405)");
	ds_lst1405->plotOn(framelst1405,Name("lst1405data"));
	lst1405shape.plotOn(framelst1405,Name("lst1405fit"),LineColor(kBlue));
	LST1405.plotOn(framelst1405,Name("lst1405fitsmooth"),LineColor(kRed),LineStyle(kDashed));
	// new TCanvas();
	// framelst1405->Draw();
	w.import(LST1405);

	//*******************************************************************

	//*********Get shape from Lambda*(1520) background*******************
	//  RooRealVar lstmass("JpsiLmass","JpsiL mass",4200.,5500.);
	//  RooRealVar MVweight("MVweight","MVweight",0.,10.);
	TFile *filein_lst1520;
	TTree *treein_lst1520 = (TTree*)malloc(sizeof(*treein_lst1520));

	filein_lst1520 = Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1520/Lst1520_total_jpsil.root");
	treein_lst1520 = (TTree*)filein_lst1520->Get("MCDecayTree");
	treein_lst1520->SetBranchStatus("*",0);
	treein_lst1520->SetBranchStatus("JpsiLmass",1);

	//  treein_lst->SetAlias("Lb_DTF_M_JpsiLConstr","JpsiLmass");

	//  RooDataSet ds_lst("ds_lst","ds_lst",treein_lst,RooArgSet(lstmass,MVweight),"MVweight");
	treein_lst1520->Draw(Form("JpsiLmass>>hlst1520(%d,%d,%d)",nbins,low,high),"","goff");

	TH1D *hlst1520        = (TH1D*)gDirectory->Get("hlst1520");
	TH1D *hlst1520_smooth = (TH1D*)hlst1520->Clone("hlst1520_smooth");
	hlst1520_smooth->Smooth(2);//TODO TWEAK THIS

	//  RooDataHist *ds_lst1 = ds_lst.binnedClone();
	RooDataHist *ds_lst1520 = new RooDataHist("ds_lst1520","ds_lst1520",*(w.var("bMass")),hlst1520);
	RooDataHist *ds_lst1520_smooth = new RooDataHist("ds_lst1520_smooth","ds_lst1520_smooth",*(w.var("bMass")),hlst1520_smooth);

	ds_lst1520->Print();

	//  RooKeysPdf lstshape("lstshape","lstshape",lstmass,ds_lst);
	RooHistPdf lst1520shape("lst1520shape","lst1520shape",*(w.var("bMass")),*ds_lst1520,0);
	RooHistPdf LST1520("LST1520","LST1520",*(w.var("bMass")),*ds_lst1520_smooth,0);

	RooPlot *framelst1520 = (w.var("bMass"))->frame();
	framelst1520->SetTitle("J/#psi Lambda^{*} (1520)");
	ds_lst1520->plotOn(framelst1520,Name("lst1520data"));
	lst1520shape.plotOn(framelst1520,Name("lst1520fit"),LineColor(kBlue));
	LST1520.plotOn(framelst1520,Name("lst1520fitsmooth"),LineColor(kRed),LineStyle(kDashed));
	// new TCanvas();
	// framelst1520->Draw();

	w.import(LST1520);
	//*******************************************************************

	RooHistPdf* XIB[2];

	//******************Get shape from Xib background********************
	// RooRealVar xibmass("Lb_DTF_M_JpsiLConstr","xibmass",5200.,5740.);
	const char* xibPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/JpsiXi";

	for(Int_t run = 1; run<=2; run++) {
		Int_t i = run-1;
		TFile *filein_xi_nonZero = Open(Form("%s/run%d/jpsixi_cutoutks_LL_nonZeroTracks.root",xibPath,run));
		TTree *treein_xi_nonZero = (TTree*)filein_xi_nonZero->Get("MyTuple");

		TFile *filein_xi_Zero = Open(Form("%s/run%d/jpsixi_cutoutks_LL_ZeroTracks.root",xibPath,run));
		TTree *treein_xi_Zero = (TTree*)filein_xi_Zero->Get("MyTuple");

		treein_xi_nonZero->AddFriend("MyTuple",Form("%s/run%d/jpsixi_LL_FinalBDT%d_iso%d_%s.root",
		                                            xibPath,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		treein_xi_Zero->AddFriend("MyTuple",Form("%s/run%d/jpsixi_zeroTracksLL_FinalBDT%d.root",
		                                         xibPath,run,bdtConf_Zero[i]));
		// treein_xi->SetBranchStatus("*",0);
		// treein_xi->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
		//RooDataSet ds_xi("ds_xi","ds_xi",treein_xi,ximass);
		treein_xi_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_nonZero%d(%d,%d,%d)",run,nbins,low,high),
		                        Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");//TRUTH MATCHING HERE
		treein_xi_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_Zero%d(%d,%d,%d)",run,nbins,low,high),
		                     Form("Lb_BKGCAT==40 && BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");//TRUTH MATCHING HERE

		TH1D *hxib_nonZero = (TH1D*)gDirectory->Get(Form("hxib_nonZero%d",run));
		TH1D *hxib_Zero    = (TH1D*)gDirectory->Get(Form("hxib_Zero%d",run));
		TH1D *hxib         = new TH1D(Form("hxib%d",run),"",nbins,low,high);
		hxib->Add(hxib_nonZero,hxib_Zero);

		TH1D *hxib_smooth  = (TH1D*)hxib->Clone(Form("hxib_smooth%d",run));
		hxib_smooth->Smooth(4);//TODO TWEAK THIS

		RooDataHist *ds_xib        = new RooDataHist(Form("ds_xib%d",run),Form("ds_xib%d",run),*(w.var("bMass")),hxib);
		RooDataHist *ds_xib_smooth = new RooDataHist(Form("ds_xib_smooth%d",run),Form("ds_xib_smooth%d",run),*(w.var("bMass")),hxib_smooth);

		ds_xib_smooth->Print();

		//  RooKeysPdf xibshape("xibshape","xibshape",xibmass,ds_xib,RooKeysPdf::NoMirror);
		RooHistPdf xibshape(Form("xibshape%d",run),Form("xibshape%d",run),*(w.var("bMass")),*ds_xib,0);
		XIB[i] = new RooHistPdf(Form("XIB%d",run),Form("XIB%d",run),*(w.var("bMass")),*ds_xib_smooth,0);

		// RooPlot *framexib = (w.var("bMass"))->frame();
		// framexib->SetTitle("J/#psi #Xi");
		// ds_xib->plotOn(framexib,Name("xibdata"));
		// xibshape.plotOn(framexib,Name("xibfit"),LineColor(kBlue));
		// XIB[i].plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));
		//  xibshape1.plotOn(framexib,Name("xibfit1"),LineColor(kRed));
		// new TCanvas();
		// framexib->Draw();

		w.import(*(XIB[i]));
	}
	//*******************************************************************

	RooHistPdf* SIG[2];
	//******************Get shape from Jpsi Sigma signal*****************
	//  RooRealVar sigmamass("Lb_DTF_M_JpsiLConstr","sigmamass",5200.,5740.);

	const char* sigmaPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/JpsiSigma";

	for(Int_t run = 1; run<=2; run++) {
		Int_t i = run-1;
		TFile *filein_sigma_nonZero = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_nonZeroTracks.root",sigmaPath,run));
		TTree *treein_sigma_nonZero = (TTree*)filein_sigma_nonZero->Get("MyTuple");

		TFile *filein_sigma_Zero = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_ZeroTracks.root",sigmaPath,run));
		TTree *treein_sigma_Zero = (TTree*)filein_sigma_Zero->Get("MyTuple");

		treein_sigma_nonZero->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_LL_FinalBDT%d_iso%d_%s.root",
		                                               sigmaPath,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		treein_sigma_Zero->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_zeroTracksLL_FinalBDT%d.root",
		                                            sigmaPath,run,bdtConf_Zero[i]));
		// treein_sigma->SetBranchStatus("*",0);
		// treein_sigma->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

		//RooDataSet ds_sigma("ds_sigma","ds_sigma",treein_sigma,sigmamass);
		treein_sigma_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hsigma_nonZero%d(%d,%d,%d)",run,nbins,low,high),
		                           Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");//Not TRUTH MATCHING HERE!

		treein_sigma_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hsigma_Zero%d(%d,%d,%d)",run,nbins,low,high),
		                        Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");//Not TRUTH MATCHING HERE!

		TH1D *hsigma_nonZero = (TH1D*)gDirectory->Get(Form("hsigma_nonZero%d",run));
		TH1D *hsigma_Zero    = (TH1D*)gDirectory->Get(Form("hsigma_Zero%d",run));
		TH1D *hsigma         = new TH1D(Form("hsigma%d",run),"",nbins,low,high);

		hsigma->Add(hsigma_nonZero,hsigma_Zero);
		TH1D *hsigma_smooth  = (TH1D*)hsigma->Clone(Form("hsigma_smooth%d",run));
		hsigma_smooth->Smooth(2);

		RooDataHist *ds_sigma = new RooDataHist(Form("ds_sigma%d",run),Form("ds_sigma%d",run),*(w.var("bMass")),hsigma);
		RooDataHist *ds_sigma_smooth = new RooDataHist(Form("ds_sigma_smooth%d",run),Form("ds_sigma_smooth%d",run),*(w.var("bMass")),hsigma_smooth);

		ds_sigma->Print();

		//  RooKeysPdf sigmashape("sigmashape","sigmashape",sigmamass,ds_sigma,RooKeysPdf::NoMirror);
		//  RooKeysPdf sigmashape("sigmashape","sigmashape",Lb_DTF_M_JpsiLConstr,ds_sigma,RooKeysPdf::NoMirror);
		//sigmashape1("sigmashape1","sigmashape1",sigmamass,ds_sigma,RooKeysPdf::MirrorBoth,2) ;
		RooHistPdf sigmashape(Form("sigmashape%d",run),Form("sigmashape%d",run),*(w.var("bMass")),*ds_sigma,0);
		SIG[i] = new RooHistPdf(Form("SIG%d",run),Form("SIG%d",run),*(w.var("bMass")),*ds_sigma_smooth,0);

		// RooPlot *framesigma = (w.var("bMass"))->frame();
		// framesigma->SetTitle("J/#psi #Sigma");
		// ds_sigma->plotOn(framesigma,Name("sigmadata"));
		// sigmashape.plotOn(framesigma,Name("sigmafit"),LineColor(kBlue));
		// SIG.plotOn(framesigma,Name("sigmafitsmooth"),LineColor(kRed),LineStyle(kDashed));
		//  sigmashape1.plotOn(framesigma,Name("sigmafit1"),LineColor(kRed));

		// new TCanvas();
		// framesigma->Draw();
		w.import(*(SIG[i]));
	}
	//*******************************************************************

	//*********Double Crystal Ball signal shape for Lambda_b0************

	// w.factory("RooCBShape::Lb1_Run1(bMass,mean_Run1[5619.6,5619,5621],"
	//           "sigma_Run1[10.,0.,20.], alpha1_Run1[1.028, 0.8,1.3], 10.0)" );
	// w.factory("RooCBShape::Lb2_Run1(bMass,mean_Run1,sigma_Run1,alpha2_Run1[-1.097,-1.4,-0.7], 10.0)");
	// w.factory("SUM::Lb_Run1(0.5*Lb1_Run1 , 0.5*Lb2_Run1)");
	//
	// w.factory("RooCBShape::Lb1_Run2(bMass,mean_Run2[5619.6,5619,5621],"
	//           "sigma_Run2[10.,0.,20.], alpha1_Run2[1.028, 0.8,1.3], 10.0)" );
	// w.factory("RooCBShape::Lb2_Run2(bMass,mean_Run2,sigma_Run2,alpha2_Run2[-1.097,-1.4,-0.7], 10.0)");
	// w.factory("SUM::Lb_Run2(0.5*Lb1_Run2 , 0.5*Lb2_Run2)");
	//*******************************************************************

	//*********Hypatia signal shape for Lambda_b0************************
	// w.factory("RooHypatia2::Lb_Run1(bMass,lambda_Run1[-2.0,-3.0,0.0],0,0,"
	//           "sigma_Run1[10.,0.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[1.7,1.0,3.0],"
	//           "2 ,a2_Run1[1.8,1.0,3.0], 2)");
	//
	// w.factory("RooHypatia2::Lb_Run2(bMass,lambda_Run2[-2.0,-3.0,0.0],0,0,"
	//           "sigma_Run2[10.,0.,20.], mean_Run2[5619.6,5619,5621], a1_Run2[1.5,1.0,3.0],"
	//           "2 ,a2_Run2[1.5,1.0,3.0], 2)");

	w.factory("RooHypatia2::Lb_Run1(bMass,lambda_Run1[-2.0,-3.0,0.0],0,0,"
	          "sigma_Run1[10.,0.,20.], mean_Run1[5619.6,5619,5621], a1_Run1[1.7,1.0,3.0],"
	          "2 ,a2_Run1[1.8,1.0,3.0], 2)");

	w.factory("RooHypatia2::Lb_Run2(bMass,lambda_Run2[-2.0,-5.0,0.0],0,0,"
	          "sigma_Run2[10.,0.,20.], mean_Run2[5619.6,5619,5621], a1_Run2[1.5,1.0,5.0],"
	          "2 ,a2_Run2[1.5,1.0,3.0], 2)");

	//*********Exponential shape for continuum backgruond****************
	// w.factory("Exponential::Bkg_Run1(bMass,tau_Run1[-0.0007,-0.01,-0.0000001])");
	// w.factory("Exponential::Bkg_Run2(bMass,tau_Run2[-0.0007,-0.01,-0.0000001])");

	w.factory("Exponential::Bkg_Run1(bMass,tau_Run1[-0.001,-0.01,-0.0000001])");
	w.factory("Exponential::Bkg_Run2(bMass,tau_Run2[-0.001,-0.01,-0.0000001])");

	// w.factory("Chebychev::Bkg_Run1(bMass, {a0_Run1[-1.44,-1.6,-1.0], a1_Run1[0.59,0.0,1.0], a2_Run1[-0.13,-0.5,0.0]})");
	// w.factory("Chebychev::Bkg_Run2(bMass, {a0_Run2[-1.46,-1.6,-1.0], a1_Run2[0.64,0.0,1.0], a2_Run2[-0.19,-0.5,0.0]})");

	// w.factory("Chebychev::Bkg_Run1(bMass, {c0_Run1[-0.5,-2.0,2.0], c1_Run1[0.5,-1.0,1.0]})");
	// w.factory("Chebychev::Bkg_Run2(bMass, {c0_Run2[-0.5,-2.0,2.0], c1_Run2[-0.5,-1.0,1.0]})");

	//*******************************************************************

	//*********Gaussian Lump for misc. Lambda*'s ************************
	w.factory("Gaussian::lstLump_Run1(bMass,miscLstMean_Run1[5009.,4980.,5050.],"
	          "miscLstSigma_Run1[60.,20.,100.])");

	w.factory("Gaussian::lstLump_Run2(bMass,miscLstMean_Run2[5009.,4980.,5050.],"
	          "miscLstSigma_Run2[60.,20.,100.])");
	//*******************************************************************

	RooDataHist *ds[2];
	Int_t nentries[2];
	TH1D *myhist[2];
	//*********Input Data************************************************
	const char *dataPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda";

	for(Int_t run = 1; run<=2; run++) {
		cout<<"Importing Run "<<run<<" data"<<endl;
		Int_t i = run-1;

		TFile *filein_nonZero = Open(Form("%s/run%d/"
		                                  "jpsilambda_cutoutks_LL_nonZeroTracks.root",dataPath,run),"READ");
		TTree *treein_nonZero = (TTree*)filein_nonZero->Get("MyTuple");

		TFile *filein_Zero = Open(Form("%s/run%d/"
		                               "jpsilambda_cutoutks_LL_ZeroTracks.root",dataPath,run),"READ");
		TTree *treein_Zero = (TTree*)filein_Zero->Get("MyTuple");

		treein_nonZero->AddFriend("MyTuple",Form("%s/run%d/"
		                                         "jpsilambda_LL_FinalBDT%d_iso%d_%s.root",
		                                         dataPath,run,bdtConf_nonZero[i],isoConf[i],isoVersion[i]));
		treein_Zero->AddFriend("MyTuple",Form("%s/run%d/"
		                                      "jpsilambda_zeroTracksLL_FinalBDT%d.root",
		                                      dataPath,run,bdtConf_Zero[i]));

		// Int_t nentries = treein_nonZero->GetEntries() + treein_Zero->GetEntries();

		treein_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_nonzero%d(%d,%d,%d)",run,nbins,low,high),
		                     Form("BDT%d > %f",bdtConf_nonZero[i],bdtCut_nonZero[i]),"goff");
		treein_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_zero%d(%d,%d,%d)",run,nbins,low,high),
		                  Form("BDT%d > %f",bdtConf_Zero[i],bdtCut_Zero[i]),"goff");
		//  treein->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist1(%d,%d,%d)",nbins,low,high),"BDT1 > 0.065","goff");

		TH1D *myhist_nonzero = (TH1D*)gDirectory->Get(Form("myhist_nonzero%d",run));
		TH1D *myhist_zero = (TH1D*)gDirectory->Get(Form("myhist_zero%d",run));

		myhist[i] = new TH1D(Form("myhist%d",run),"",nbins,low,high);

		myhist[i]->Add(myhist_zero,myhist_nonzero);

		nentries[i] = myhist[i]->Integral();
		cout<<"nentries = "<<nentries[i]<<endl;

		ds[i] = new RooDataHist(Form("ds%d",run),Form("ds%d",run),*(w.var("bMass")),myhist[i]);
		//  RooDataSet ds("ds","ds",treein,Lb_DTF_M_JpsiLConstr);
		cout<<"Done making RooDataHist"<<endl;
		(ds[i])->Print();

		w.import(*(ds[i]));
	}
	//*******************************************************************

	//************************MAKE COMBINED MODEL************************
	// w.factory("R[-6.0,-10.0,0.0]"); //R is the parameter of interest  This is shared b/w Run1 and Run2
	// w.var("R")->setError(0.1);

	w.factory("R[200,0,10000]"); //R*10^5 is the parameter of interest  This is shared b/w Run1 and Run2
	// w.var("R")->setError(0.1);

	// w.factory(Form("sigmaEff1[%f,0.00001,0.005]",eff_Sigma[0]));
	// w.factory(Form("lambdaEff1[%f,0.00001,0.005]",eff_Lambda[0]));
	// w.factory(Form("sigmaEff2[%f,0.00001,0.005]",eff_Sigma[1]));
	// w.factory(Form("lambdaEff2[%f,0.00001,0.005]",eff_Lambda[1]));
	//
	// w.factory(Form("Gaussian::sigmaEff_constraint1(gsigmaEff1[%f,0.00001,0.005],sigmaEff1,%f)",eff_Sigma[0],eff_Sigma_err[0]));
	// w.factory(Form("Gaussian::lambdaEff_constraint1(glambdaEff1[%f,0.00001,0.005],lambdaEff1,%f)",eff_Lambda[0],eff_Lambda_err[0]));
	//
	// w.factory(Form("Gaussian::sigmaEff_constraint2(gsigmaEff2[%f,0.00001,0.005],sigmaEff2,%f)",eff_Sigma[1],eff_Sigma_err[1]));
	// w.factory(Form("Gaussian::lambdaEff_constraint2(glambdaEff2[%f,0.00001,0.005],lambdaEff2,%f)",eff_Lambda[1],eff_Lambda_err[1]));

	w.factory(Form("sigmaEff1[%f,0,1]",eff_Sigma[0]));
	w.factory(Form("lambdaEff1[%f,0,1]",eff_Lambda[0]));
	w.factory(Form("sigmaEff2[%f,0,1]",eff_Sigma[1]));
	w.factory(Form("lambdaEff2[%f,0,1]",eff_Lambda[1]));

	w.factory(Form("Gaussian::sigmaEff_constraint1(gsigmaEff1[%f,0,1],sigmaEff1,%f)",eff_Sigma[0],eff_Sigma_err[0]));
	w.factory(Form("Gaussian::lambdaEff_constraint1(glambdaEff1[%f,0,1],lambdaEff1,%f)",eff_Lambda[0],eff_Lambda_err[0]));

	w.factory(Form("Gaussian::sigmaEff_constraint2(gsigmaEff2[%f,0,1],sigmaEff2,%f)",eff_Sigma[1],eff_Sigma_err[1]));
	w.factory(Form("Gaussian::lambdaEff_constraint2(glambdaEff2[%f,0,1],lambdaEff2,%f)",eff_Lambda[1],eff_Lambda_err[1]));

	w.var("gsigmaEff1")->setConstant();
	w.var("glambdaEff1")->setConstant();

	w.var("gsigmaEff2")->setConstant();
	w.var("glambdaEff2")->setConstant();

	w.factory("nLb_Run1[0,5000]");
	w.factory("nLb_Run2[0,17000]");

	//What should the limits on nXib be?
	Double_t xibCentral_run1 = xibnorm_LL[0];
	Double_t xibErr_run1 = xibnorm_LL_err[0];
	Double_t xibLow_run1 = 0;//xibCentral_run1 - 2*xibErr_run1;
	Double_t xibHigh_run1 = 200;//xibCentral_run1 + 2*xibErr_run1;

	Double_t xibCentral_run2 = xibnorm_LL[1];
	Double_t xibErr_run2 = xibnorm_LL_err[1];
	Double_t xibLow_run2 = 0;//xibCentral_run2 - 2*xibErr_run2;
	Double_t xibHigh_run2 = 200;//xibCentral_run2 + 2*xibErr_run2;

	w.factory(Form("nXib1[%f,%f,%f]",xibCentral_run1,xibLow_run1,xibHigh_run1));
	w.factory(Form("nXib2[%f,%f,%f]",xibCentral_run2,xibLow_run2,xibHigh_run2));

	w.factory(Form("Gaussian::nXib_constraint1(gnXib1[%f,%f,%f],nXib1,%f)",xibCentral_run1,xibLow_run1,xibHigh_run1,xibErr_run1));
	w.factory(Form("Gaussian::nXib_constraint2(gnXib2[%f,%f,%f],nXib2,%f)",xibCentral_run2,xibLow_run2,xibHigh_run2,xibErr_run2));

	w.var("gnXib1")->setConstant();
	w.var("gnXib2")->setConstant();

	// w.factory("expr::nSigma1('pow(10,R)*nLb_Run1*sigmaEff1/(lambdaEff1)',R,nLb_Run1,sigmaEff1,lambdaEff1)");
	// w.factory("expr::nSigma2('pow(10,R)*nLb_Run2*sigmaEff2/(lambdaEff2)',R,nLb_Run2,sigmaEff2,lambdaEff2)");

	w.factory("expr::nSigma1('pow(10,-5)*R*nLb_Run1*sigmaEff1/(lambdaEff1)',R,nLb_Run1,sigmaEff1,lambdaEff1)");
	w.factory("expr::nSigma2('pow(10,-5)*R*nLb_Run2*sigmaEff2/(lambdaEff2)',R,nLb_Run2,sigmaEff2,lambdaEff2)");

	// w.factory("expr::nSigma('pow(10,R)*nLb',R,nLb)");//TAKING EFFICIENCIES OUT OF PICTURE FOR NOW. TESTING.

	// w.factory(Form("SUM:model1(nSigma1*SIG1 , nLb_Run1*Lb_Run1 , nXib1*XIB1 , n1405_Run1[1500,1,5000]*LST1405 , n1520_Run1[3000,1,5000]*LST1520 , nMiscLst_Run1[2500,1,5000]*lstLump_Run1 , nBkg_Run1[2000,1,%d]*Bkg_Run1)",nentries[0]));
	// w.factory(Form("SUM:model2(nSigma2*SIG2 , nLb_Run2*Lb_Run2 , nXib2*XIB2 , n1405_Run2[1500,1,5000]*LST1405 , n1520_Run2[3000,1,5000]*LST1520 , nMiscLst_Run2[2500,1,5000]*lstLump_Run2 , nBkg_Run2[2000,1,%d]*Bkg_Run2)",nentries[1]));

	w.factory(Form("SUM:model1(nSigma1*SIG1 , nLb_Run1*Lb_Run1 , nXib1*XIB1 , n1405_Run1[1500,50,10000]*LST1405 , n1520_Run1[3000,1,10000]*LST1520 , nMiscLst_Run1[2500,1,5000]*lstLump_Run1 , nBkg_Run1[2000,1,%d]*Bkg_Run1)",nentries[0]));
	w.factory(Form("SUM:model2(nSigma2*SIG2 , nLb_Run2*Lb_Run2 , nXib2*XIB2 , n1405_Run2[1500,50,10000]*LST1405 , n1520_Run2[3000,1,10000]*LST1520 , nMiscLst_Run2[2500,1,5000]*lstLump_Run2 , nBkg_Run2[2000,1,%d]*Bkg_Run2)",nentries[1]));

	if(!lst1520flag)
	{
		w.var("n1520_Run1")->setVal(0);
		w.var("n1520_Run1")->setConstant();

		w.var("n1520_Run2")->setVal(0);
		w.var("n1520_Run2")->setConstant();
	}
	w.factory("PROD::model_const1(model1,sigmaEff_constraint1,lambdaEff_constraint1,nXib_constraint1)");//Multiply model by constraint terms
	w.factory("PROD::model_const2(model2,sigmaEff_constraint2,lambdaEff_constraint2,nXib_constraint2)");//Multiply model by constraint terms

	RooAbsPdf* model_const1 = w.pdf("model_const1"); // get the model
	RooAbsPdf* model_const2 = w.pdf("model_const2"); // get the model

	w.defineSet("poi","R"); //parameters of interest

	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,lambda_Run1,a1_Run1,a2_Run1,nBkg_Run1,tau_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,lambda_Run2,a1_Run2,a2_Run2,nBkg_Run2,tau_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff1,lambdaEff1,nXib1,n1405_Run2,n1520_Run2");

	w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,lambda_Run1,a1_Run1,a2_Run1,nBkg_Run1,tau_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,lambda_Run2,a1_Run2,a2_Run2,nBkg_Run2,tau_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff1,lambdaEff1,nXib1,n1405_Run2,n1520_Run2");

	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,alpha1_Run1,alpha2_Run1,nBkg_Run1,tau_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,alpha1_Run2,alpha2_Run2,nBkg_Run2,tau_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff2,lambdaEff2,nXib2,n1405_Run2,n1520_Run2");

	// w.defineSet("nuisParams","nLb_Run1,mean_Run1,sigma_Run1,alpha1_Run1,alpha2_Run1,nBkg_Run1,a0_Run1,a1_Run1,nMiscLst_Run1,miscLstMean_Run1,miscLstSigma_Run1,sigmaEff1,lambdaEff1,nXib1,n1405_Run1,n1520_Run1");// define set of nuisance parameters
	// w.extendSet("nuisParams","nLb_Run2,mean_Run2,sigma_Run2,alpha1_Run2,alpha2_Run2,nBkg_Run2,a0_Run2,a1_Run2,nMiscLst_Run2,miscLstMean_Run2,miscLstSigma_Run2,sigmaEff2,lambdaEff2,nXib2,n1405_Run2,n1520_Run2");

	w.defineSet("globObs","gsigmaEff1,glambdaEff1,gnXib1,gsigmaEff2,glambdaEff2,gnXib2");//define set of global observables
	//*******************************************************************

	//***********************MAKE COMBINED DATASET************************************
	RooCategory sample("sample","sample");
	sample.defineType("run1");
	sample.defineType("run2");

	// Construct combined dataset in (x,sample)
	RooDataHist combData("combData","combined data",*(w.var("bMass")),Index(sample),Import("run1",*myhist[0]),Import("run2",*myhist[1]));

	// Construct a simultaneous pdf using category sample as index
	RooSimultaneous simPdf("simPdf","simultaneous pdf",sample);

	// Associate model with the physics state and model_ctl with the control state
	simPdf.addPdf(*model_const1,"run1");
	simPdf.addPdf(*model_const2,"run2");

	w.import(combData);
	w.import(simPdf);
	// w.import(sample);

	w.defineSet("obs","bMass,sample");

	cout<<"Done importing combData and simPdf"<<endl;
	w.Print();
	//***********************MAKE NLL************************************
	// RooAbsReal* nll = w.pdf("model_const")->createNLL(ds);
	//
	// new TCanvas();
	// RooPlot* frame = w.var("R")->frame();
	//
	// nll->plotOn(frame,ShiftToZero()); //the ShiftToZero option puts the minimum at 0 on the y-axis
	// frame->Draw();
	//
	// RooAbsReal* pll = nll->createProfile(*w.set("poi"));
	// // RooAbsReal* pll = nll->createProfile(*w.var("nLb"));
	// pll->plotOn(frame,LineColor(kRed));
	// frame->Draw();
	// frame->SetAxisRange(0,10,"Y");
	//
	// RooMinuit mini(*nll);
	// // mini.setVerbose(1);
	// mini.setStrategy(2);
	// mini.minos(*w.set("poi")); //could give no arg list and it will calculate minos errors for all the parameters
	// // mini.minos(*w.var("nLb")); //could give no arg list and it will calculate minos errors for all the parameters
	// // mini.improve();
	//
	// RooFitResult* res = mini.save("myResult","My Result");
	// cout<<"STATUS = "<<res->status()<<endl; //should be 0 for success!
	//
	// cout<<"MIN "<<w.var("R")->getVal()<<" ERR LO"<<w.var("R")->getErrorLo()<<" ERR HI"<<w.var("R")->getErrorHi()<<endl;
	//
	// RooFormulaVar p_mu("p_mu","p_{#mu} using asymptotic formulae","TMath::Prob(2*@0,1.)",RooArgList(*pll));
	//
	// TCanvas *c2 = new TCanvas();
	// RooPlot* frame1 = w.var("R")->frame();
	// p_mu.plotOn(frame1,LineColor(kGreen));
	// // frame1->Draw();
	//
	// c2->cd();
	//
	// Double_t limit = p_mu.findRoot(*(w.var("R")),-5.0,0.0,0.05);
	// cout<<"UL ="<<limit<<endl;
	//
	// RooOneSidedProfileLL opll("opll","opll",*pll); //pll is the pointer to your profilell function made earlier
	// opll.plotOn(frame1,LineColor(kMagenta));
	//
	// RooFormulaVar p_mu_q("p_mu_q","p_{#mu} from one-sided pll asymptotic formulae","ROOT::Math::normal_cdf_c(sqrt(2*@0),1.)",RooArgList(opll));
	// p_mu_q.plotOn(frame1,LineColor(kCyan));
	//
	// frame->GetYaxis()->SetRangeUser(-0.1,1.1);
	// frame1->Draw();
	//
	// TLine l; l.SetLineStyle(2);
	// l.DrawLine(-10,0.05,0,0.05);
	//
	// Double_t limit1 = p_mu_q.findRoot(*(w.var("R")),-5.0,0.0,0.05);
	// cout<<"One sided UL ="<<limit1<<endl;
	//
	// yields.push_back(0.0);
	// return yields;

	//************************DO THE FIT***********************
	RooFitResult *res = simPdf.fitTo(combData,Minos(*w.set("poi")),Extended(), Save(), Hesse(false), Strategy(2));
	//*******************************************************************

	//*********************PLOTTING STUFF*********************************************
	TCanvas* c_run1 = new TCanvas("Run1","Run1", 1200, 800);

	RooPlot *frame_run1 = new RooPlot(*(w.var("bMass")),low,high,100);
	frame_run1->SetTitle("Run1 Fit");
	// frame_run1->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");

	combData.plotOn(frame_run1,Name("data_Run1"),Cut("sample==sample::run1"),DataError(RooAbsData::Poisson));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Name("fit_Run1"));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("Lb_Run1"))),Name("lb_Run1"),LineColor(kMagenta+2));
	// simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("Lb1_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	// simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("Lb2_Run1"))),LineStyle(kDotted),LineColor(kMagenta));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("Bkg_Run1"))),LineColor(kRed),Name("bkg_Run1"));
	//  if(xibflag!=0)
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("XIB1"))),LineColor(kGreen),Name("xib_Run1"));
	if(lst1405flag!=0)
		simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("LST1405"))),LineColor(kGreen+2),LineStyle(kDashed),Name("lst1405_Run1"));
	if(lst1520flag!=0)
		simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("LST1520"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1520_Run1"));
	//simPdf.plotOn(frame_run1,Components(lst1810shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1810"));
	// if(sigmaflag!=0)
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("SIG1"))),LineColor(kBlack),Name("sig_Run1"));
	simPdf.plotOn(frame_run1,Slice(sample,"run1"),ProjWData(sample,combData),Components(*(w.pdf("lstLump_Run1"))),LineColor(kRed+2),LineStyle(kDashed),Name("misclst_Run1"));

	frame_run1->GetYaxis()->SetRangeUser(0,60);

	// Double_t chiSquare1 = frame_run1->chiSquare("fit_run1","data_Run1");
	// cout<<"chi square1/dof = "<<chiSquare1<<endl;
	RooArgSet *allpar_run1 = simPdf.getParameters(*(ds[0]));
	RooArgSet *floatpar_run1 = (RooArgSet*)allpar_run1->selectByAttrib("Constant",kFALSE);
	floatpar_run1->Print();
	int floatpars_run1 = (floatpar_run1->selectByAttrib("Constant",kFALSE))->getSize() - 1;//-1 because sample also gets included in this list
	cout<<"run1 float pars = "<<floatpars_run1<<endl;
	Double_t chi2_run1 = frame_run1->chiSquare("fit_Run1","data_Run1",floatpars_run1);
	cout<<"chi square2/dof = "<<chi2_run1<<endl;

	Int_t fit_ndof_run1 = nbins - floatpars_run1;
	cout<<"chi square2 = "<<chi2_run1*fit_ndof_run1<<endl;

	///////////
	TPad *pad1 = new TPad("pad1","pad1",0.0,0.2,1.0,1.0);
	TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.2);

	pad1->SetGridx();
	pad1->SetGridy();
	pad2->SetGridx();
	pad2->SetGridy();

	pad1->SetBottomMargin(0.04);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->SetBorderMode(0);
	pad1->SetBorderMode(0);
	c_run1->SetBorderMode(0);
	pad2->Draw();
	pad1->Draw();
	pad1->cd();
	gPad->SetTopMargin(0.06);
	pad1->Update();

	frame_run1->Draw();

	TLatex l_run1;
	l_run1.SetTextSize(0.025);
	l_run1.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run1));

	c_run1->Modified();

	auto legend_run1 = new TLegend(0.1,0.7,0.3,0.9);
	legend_run1->SetTextSize(0.025);
	legend_run1->AddEntry("data_Run1","Data","lp");
	legend_run1->AddEntry("fit_Run1","Total Fit","l");
	legend_run1->AddEntry("lb_Run1","J/#psi #Lambda shape","l");
	legend_run1->AddEntry("sig_Run1","J/#psi #Sigma shape","l");
	legend_run1->AddEntry("xib_Run1","J/#psi #Xi shape","l");
	if(lst1405flag)
	{
		legend_run1->AddEntry("lst1405_Run1","J/#psi #Lambda(1405) shape","l");
	}
	if(lst1520flag)
	{
		legend_run1->AddEntry("lst1520_Run1","J/#psi #Lambda(1520) shape","l");
	}
	// if(lst1810flag)
	//      legend_run1->AddEntry("lst1810","J/#psi #Lambda(1810) shape","l");
	legend_run1->AddEntry("misclst_Run1","misc. J/#psi #Lambda* shapes","l");
	legend_run1->AddEntry("bkg_Run1","Comb. Bkg. shape","l");
	legend_run1->Draw("same");

	c_run1->Update();

	// Pull distribution
	RooPlot *frame_run1x2 = (w.var("bMass"))->frame();
	// RooPlot *framex2 = new RooPlot(*(w.var("bMass")), low,5800,100);
	RooHist* hpull_run1 = frame_run1->pullHist("data_Run1","fit_Run1");
	frame_run1x2->addPlotable(hpull_run1,"P");
	hpull_run1->SetLineColor(kBlack);
	hpull_run1->SetMarkerColor(kBlack);
	frame_run1x2->SetTitle(0);
	frame_run1x2->GetYaxis()->SetTitle("Pull");
	frame_run1x2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame_run1x2->GetYaxis()->SetTitleSize(0.15);
	frame_run1x2->GetYaxis()->SetLabelSize(0.15);
	frame_run1x2->GetXaxis()->SetTitleSize(0.15);
	frame_run1x2->GetXaxis()->SetLabelSize(0.15);
	frame_run1x2->GetYaxis()->CenterTitle();
	frame_run1x2->GetYaxis()->SetTitleOffset(0.25);
	frame_run1x2->GetXaxis()->SetTitleOffset(0.75);
	frame_run1x2->GetYaxis()->SetNdivisions(505);
	frame_run1x2->GetYaxis()->SetRangeUser(-5.0,5.0);
	pad2->cd();
	frame_run1x2->Draw();

	c_run1->cd();
	// pad1->cd();

	cout<<"Run1 Pull Mean Y = "<<hpull_run1->GetMean(2)<<endl;
	cout<<"Run1 Pull RMS Y = "<<hpull_run1->GetRMS(2)<<endl;

	TCanvas* c_run2 = new TCanvas("Run2","Run2", 1200, 800);

	RooPlot *frame_run2 = new RooPlot(*(w.var("bMass")),low,high,100);

	frame_run2->SetTitle("Run2 Fit");
	frame_run2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	combData.plotOn(frame_run2,Name("data_Run2"),Cut("sample==sample::run2"),DataError(RooAbsData::Poisson));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Name("fit_Run2"));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("Lb_Run2"))),Name("lb_Run2"),LineColor(kMagenta+2));
	// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("Lb1_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	// simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("Lb2_Run2"))),LineStyle(kDotted),LineColor(kMagenta));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("Bkg_Run2"))),LineColor(kRed),Name("bkg_Run2"));
	//  if(xibflag!=0)
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("XIB2"))),LineColor(kGreen),Name("xib_Run2"));
	if(lst1405flag!=0)
		simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("LST1405"))),LineColor(kGreen+2),LineStyle(kDashed),Name("lst1405_Run2"));
	if(lst1520flag!=0)
		simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("LST1520"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1520_Run2"));
	//simPdf.plotOn(frame_run2,Components(lst1810shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1810"));
	// if(sigmaflag!=0)
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("SIG2"))),LineColor(kBlack),Name("sig_Run2"));
	simPdf.plotOn(frame_run2,Slice(sample,"run2"),ProjWData(sample,combData),Components(*(w.pdf("lstLump_Run2"))),LineColor(kRed+2),LineStyle(kDashed),Name("misclst_Run2"));

	frame_run2->GetYaxis()->SetRangeUser(0,160);

	// Double_t chiSquare1 = frame_run2->chiSquare("fit_run2","data_run2");
	// cout<<"chi square1/dof = "<<chiSquare1<<endl;
	RooArgSet *floatpar_run2 = simPdf.getParameters(*(ds[1]));
	floatpar_run2->Print();
	int floatpars_run2 = (floatpar_run2->selectByAttrib("Constant",kFALSE))->getSize() -1;//-1 because sample also gets included in this list
	cout<<"run2 float pars = "<<floatpars_run2<<endl;
	Double_t chi2_run2 = frame_run2->chiSquare("fit_Run2","data_Run2",floatpars_run2);
	cout<<"chi square2/dof = "<<chi2_run2<<endl;

	Int_t fit_ndof_run2 = nbins - floatpars_run2;
	cout<<"chi square2 = "<<chi2_run2*fit_ndof_run2<<endl;

	///////////
	TPad *pad3 = new TPad("pad3","pad3",0.0,0.2,1.0,1.0);
	TPad *pad4 = new TPad("pad4","pad4",0.0,0.0,1.0,0.2);

	pad3->SetGridx();
	pad3->SetGridy();
	pad4->SetGridx();
	pad4->SetGridy();

	pad3->SetBottomMargin(0.04);
	pad4->SetTopMargin(0);
	pad4->SetBottomMargin(0.3);
	pad4->SetBorderMode(0);
	pad3->SetBorderMode(0);
	c_run1->SetBorderMode(0);
	pad4->Draw();
	pad3->Draw();
	pad3->cd();
	gPad->SetTopMargin(0.06);
	pad3->Update();

	frame_run2->Draw();
	// c_run2->cd();

	auto legend_run2 = new TLegend(0.1,0.7,0.3,0.9);
	legend_run2->SetTextSize(0.025);
	legend_run2->AddEntry("data_Run2","Data","lp");
	legend_run2->AddEntry("fit_Run2","Total Fit","l");
	legend_run2->AddEntry("lb_Run2","J/#psi #Lambda shape","l");
	legend_run2->AddEntry("sig_Run2","J/#psi #Sigma shape","l");
	legend_run2->AddEntry("xib_Run2","J/#psi #Xi shape","l");
	if(lst1405flag)
	{
		legend_run2->AddEntry("lst1405_Run2","J/#psi #Lambda(1405) shape","l");
	}
	if(lst1520flag)
	{
		legend_run2->AddEntry("lst1520_Run2","J/#psi #Lambda(1520) shape","l");
	}
	// if(lst1810flag)
	//      legend_run2->AddEntry("lst1810","J/#psi #Lambda(1810) shape","l");
	legend_run2->AddEntry("misclst_Run2","misc. J/#psi #Lambda* shapes","l");
	legend_run2->AddEntry("bkg_Run2","Comb. Bkg. shape","l");
	legend_run2->Draw("same");

	TLatex l_run2;
	l_run2.SetTextSize(0.025);
	l_run2.DrawLatexNDC(0.8,0.6,Form("#chi^{2}/ndf = %f",chi2_run2));

	c_run2->Update();

	// Pull distribution
	RooPlot *frame_run2x2 = (w.var("bMass"))->frame();
	// RooPlot *framex2 = new RooPlot(*(w.var("bMass")), low,5800,100);
	RooHist* hpull_run2 = frame_run2->pullHist("data_Run2","fit_Run2");
	frame_run2x2->addPlotable(hpull_run2,"P");
	hpull_run2->SetLineColor(kBlack);
	hpull_run2->SetMarkerColor(kBlack);
	frame_run2x2->SetTitle(0);
	frame_run2x2->GetYaxis()->SetTitle("Pull");

	frame_run2x2->GetXaxis()->SetTitle("m[J/#psi #Lambda] (MeV)");
	frame_run2x2->GetYaxis()->SetTitleSize(0.15);
	frame_run2x2->GetYaxis()->SetLabelSize(0.15);
	frame_run2x2->GetXaxis()->SetTitleSize(0.15);
	frame_run2x2->GetXaxis()->SetLabelSize(0.15);
	frame_run2x2->GetYaxis()->CenterTitle();
	frame_run2x2->GetYaxis()->SetTitleOffset(0.25);
	frame_run2x2->GetXaxis()->SetTitleOffset(0.75);
	frame_run2x2->GetYaxis()->SetNdivisions(505);
	frame_run2x2->GetYaxis()->SetRangeUser(-5.0,5.0);
	pad4->cd();
	frame_run2x2->Draw();

	c_run2->cd();

	Double_t chi2_global = chi2_run1*fit_ndof_run1 + chi2_run2*fit_ndof_run2;
	Int_t ndof_global = fit_ndof_run1 + fit_ndof_run2;

	Double_t chi2_ndof_global = chi2_global/ndof_global;

	cout<<"****************************"<<endl;
	cout<<"Global Fit chi2/dof = "<<chi2_ndof_global<<endl;
	cout<<"****************************"<<endl;

	if(calcUL) {
		//*************ROOSTATS MODEL CONFIG*********************************

		cout<<"Starting Model Config"<<endl;
		RooStats::ModelConfig mc("ModelConfig",&w);
		mc.SetPdf(*(w.pdf("simPdf")));
		mc.SetParametersOfInterest(*w.set("poi"));
		mc.SetObservables(*w.set("obs"));
		mc.SetGlobalObservables(*w.set("globObs"));
		mc.SetNuisanceParameters(*w.set("nuisParams"));

		// import model in the workspace
		w.import(mc);
		//*******************************************************************

		//***************ASIMOV DATASET*************************************
		cout<<"Starting Asimov Dataset"<<endl;
		RooDataSet *asimovData = (RooDataSet*)RooStats::AsymptoticCalculator::GenerateAsimovData(*(mc.GetPdf()),*(mc.GetObservables()));
		RooPlot *frame_asim = new RooPlot(*(w.var("bMass")), low,5800,100);

		TCanvas *asim = new TCanvas();
		asimovData->plotOn(frame_asim);
		frame_asim->Draw();

		mc.SetSnapshot(*(w.var("R")));

		RooStats::ModelConfig *bkgOnlyModel = mc.Clone();
		bkgOnlyModel->SetName("bkgOnlyModel");
		Double_t oldval = w.var("R")->getVal();
		// w.var("R")->setVal(-10);//What value should I set for the background only hypothesis? log of 0 is -inf

		w.var("R")->setVal(0);//What value should I set for the background only hypothesis? log of 0 is -inf
		bkgOnlyModel->SetSnapshot(*(w.var("R")));

		w.var("R")->setVal(oldval);

		w.import(*bkgOnlyModel);
		// w.Write();
		//*******************************************************************

		//**************Fit like Matt****************************************

/*	cout<<"Starting Fit like Matt"<<endl;
        RooAbsReal* nll = w.pdf("simPdf")->createNLL(combData);
        RooMinimizer min_nll(*nll);
        min_nll.setErrorLevel(0.5);
        min_nll.setStrategy(2);

        Double_t initnll = nll->getVal();

        RooArgSet *pars = w.pdf("simPdf")->getParameters( combData );
        // RooArgSet *initpars = (RooArgSet*)pars->snapshot();
        // cout<<"YYYY"<<pars<<endl;
        // print "YYYY", initpars

        // TIterator *it = initpars->createIterator();
        // RooRealVar *par = (RooRealVar*)it->Next();
        // cout<<"YYYY All par list"<<endl;
        // while(par)
        // {
        //      cout<<"YYYY"<<par<<par->GetName()<<par->getVal()<<par->getError()<<endl;
        //      par = (RooRealVar*)it->Next();
        // }
        min_nll.minimize("Minuit","Migrad");

        min_nll.minos(*(w.var("R")));

        RooFitResult *res1 = min_nll.save();

        Double_t finalnll = nll->getVal();*/

		Double_t origval = w.var("R")->getVal();
		Double_t errhi   = w.var("R")->getErrorHi();
		Double_t errlo   = w.var("R")->getErrorLo();

		//**************Upper Limit Calculation*****************************

		// cout<<"Starting Upper Limit Calculation"<<endl;

		//**************AsymptoticCalculator********************************
		// RooStats::AsymptoticCalculator hc(combData, *bkgOnlyModel, mc, true);
		// hc.SetOneSided(true);
		// RooStats::HypoTestInverter calc(hc);
		//
		// calc.SetConfidenceLevel(0.95);
		// // calc.SetFixedScan(100,-7, -1);
		// calc.SetFixedScan(20,0, 1000);
		// calc.UseCLs(true);
		// // calc.SetVerbose(true);
		// HypoTestInverterResult *r = calc.GetInterval();
		//
		// cout<<"origval = "<<origval<<" error high = "
		//     <<errhi<<" error low = "<<errlo<<endl;
		// cout<<"UL = "<<r->UpperLimit()<<" +/- "<<r->UpperLimitEstimatedError()<<endl;
		//
		// new TCanvas();
		// RooStats::HypoTestInverterPlot plot("hti_plot", "hti_plot", r);
		// plot.Draw("CLb 2CL");

		//******************************************************************

		//**************FrequentistCalculator*******************************
		// RooStats::FrequentistCalculator hc(combData, *bkgOnlyModel, mc);
		// TestStatSampler *toymcs = hc.GetTestStatSampler();
		//
		// RooStats::ProfileLikelihoodTestStat teststat(*(mc.GetPdf()));
		// teststat.SetOneSided(true);
		//
		// teststat.SetPrintLevel(1);
		// teststat.SetStrategy(0);
		//
		// toymcs->SetTestStatistic(teststat);
		// toymcs->SetGenerateBinned(false);
		// toymcs->SetUseMultiGen(false);
	}
	//Save Canvases
	c_run1->SaveAs(Form("../plots/data/JpsiLambda/run1/Fit_HypatiaSig_ExpBkg_%d_%d.pdf",low,high));
	c_run2->SaveAs(Form("../plots/data/JpsiLambda/run2/Fit_HypatiaSig_ExpBkg_%d_%d.pdf",low,high));

	// write the workspace in the file
	TString fileName = Form("../rootFiles/dataFiles/JpsiLambda/ModelConfigs/MyModel_HypatiaSig_ExpBkg_%d_%d.root",low,high);
	w.writeToFile(fileName,true);
	cout << "workspace written to file " << fileName << endl;

	gSystem->RedirectOutput(0);

	// yields.push_back(0.0);
	// return yields;
}
