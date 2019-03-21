#include "Fitscript_new.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

#define Open TFile::Open
std::vector <Float_t> Fitscript_new(Int_t run, Int_t finalBDTConf_nonZero, Int_t finalBDTConf_Zero,
                                    Int_t isoConf, const char *isoVersion,
                                    Int_t mylow, Int_t constraintflag,
                                    Float_t bdtCut_nonZero, Float_t bdtCut_Zero)
{
	// gROOT->ProcessLine(".L RooOneSidedProfileLL.cxx+");
	// gROOT->ProcessLine(".L RooOneSidedProfileLL_cxx.so");
	// gSystem->Load("RooOneSidedProfileLL_cxx.so");

	std::vector <Float_t> yields;

	Double_t xibnorm_LL     = 0.0;
	Double_t xibnorm_LL_err = 0.0;

	Double_t eff_Lambda_rec     = 0.0;
	Double_t eff_Sigma_rec     = 0.0;
	Double_t eff_Lambda_rec_err = 0.0;
	Double_t eff_Sigma_rec_err  = 0.0;

	Double_t eff_Lambda_gen     = 0.0;
	Double_t eff_Sigma_gen     = 0.0;
	Double_t eff_Lambda_gen_err = 0.0;
	Double_t eff_Sigma_gen_err  = 0.0;

	Double_t eff_Lambda = 0.0;
	Double_t eff_Lambda_err = 0.0;

	Double_t eff_Sigma = 0.0;
	Double_t eff_Sigma_err = 0.0;

	Int_t nGen_Lambda       = 0;
	Int_t nGen_Sigma        = 0;

	Int_t lst1405flag = 1;
	Int_t lst1520flag = 1;
	Int_t lst1810flag = 0;
	Int_t xibflag     = 1;
	Int_t sigmaflag   = 1;
	Int_t Lst1405_rwtype = 2;

	Int_t low      = mylow, high = 7000; //Define range in which fit is performed
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

	RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

	//*****Run script to get Xib normalization*********
	gSystem->Exec(Form("python -c \'from GetXibNorm import GetNorm;"
	                   " GetNorm(%d, \"%s\", %d, %d, %d, %f, %f)\'",
	                   run,isoVersion,isoConf,finalBDTConf_nonZero,
	                   finalBDTConf_Zero,bdtCut_nonZero,bdtCut_Zero));

	ifstream infile(Form("../logs/mc/JpsiXi/run%d/xibNorm_log.txt",run));

	infile>>xibnorm_LL;
	infile>>xibnorm_LL_err;

	xibnorm_LL = xibnorm_LL*2; //ACCOUNT FOR XIB0

	xibnorm_LL_err = xibnorm_LL_err * 1.414;//ACCOUNT FOR XIB0
	cout<<"The LL Xib normalization is "<<xibnorm_LL<<" +/- "
	    <<xibnorm_LL_err<<endl;
	//************************************************
	//****Get J/psi Lambda efficiencies from MC*******
	const char* lambdaMCPath = "../rootFiles/mcFiles/JpsiLambda/JpsiLambda";

	TFile *mcFileIn_nonZero_Lambda = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_nonZeroTracks.root",
	                                           lambdaMCPath,run));
	TTree *mcTreeIn_nonZero_Lambda = (TTree*)mcFileIn_nonZero_Lambda->Get("MyTuple");

	TFile *mcFileIn_Zero_Lambda    = Open(Form("%s/run%d/jpsilambda_cutoutks_LL_ZeroTracks.root",
	                                           lambdaMCPath,run));
	TTree *mcTreeIn_Zero_Lambda    = (TTree*)mcFileIn_Zero_Lambda->Get("MyTuple");

	mcTreeIn_nonZero_Lambda->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_LL_FinalBDT%d_iso%d_%s.root",
	                                                  lambdaMCPath,run,finalBDTConf_nonZero,
	                                                  isoConf,isoVersion));
	mcTreeIn_Zero_Lambda->AddFriend("MyTuple",Form("%s/run%d/jpsilambda_zeroTracksLL_FinalBDT%d.root",
	                                               lambdaMCPath,run,finalBDTConf_Zero));
	fstream genFile_Lambda;
	genFile_Lambda.open((Form("../logs/mc/JpsiLambda/JpsiLambda/run%d/gen_log.txt",run)));

	genFile_Lambda>>nGen_Lambda;

	fstream genEffFile_Lambda;
	genEffFile_Lambda.open(Form("../logs/mc/JpsiLambda/JpsiLambda/run%d/Generator_Effs_Combined.txt",run));

	genEffFile_Lambda>>eff_Lambda_gen;
	genEffFile_Lambda>>eff_Lambda_gen_err;

	cout<<"Lambda Generator Effs = "<<eff_Lambda_gen*100<<" % +/- "<<eff_Lambda_gen_err*100<<" %"<<endl;

	Int_t num_Lambda = mcTreeIn_nonZero_Lambda->GetEntries(Form("BDT%d > %f", finalBDTConf_nonZero,bdtCut_nonZero))
	                   + mcTreeIn_Zero_Lambda->GetEntries(Form("BDT%d > %f", finalBDTConf_Zero,bdtCut_Zero));

	eff_Lambda_rec     = num_Lambda*1.0/nGen_Lambda;
	// eff_Lambda_rec_err = eff_Lambda_rec*0.02;

	cout<<"Lambda Recons. Effs = "<<eff_Lambda_rec*100<<" %"<<endl;

	eff_Lambda = eff_Lambda_rec*eff_Lambda_gen;
	eff_Lambda_err = eff_Lambda*0.02;// nominally assigning 2% uncertainty for now
	// eff_Lambda_err = sqrt(pow(eff_Lambda_gen_err,2) + pow(eff_Lambda_rec_err,2)); //Combining errors in quadrature

	cout<<"Jpsi Lambda Eff = "<<eff_Lambda*100<<" % +/- "<<eff_Lambda_err*100<<" %"<<endl;
	//************************************************
	//****Get J/psi Sigma efficiencies from MC*******
	const char *sigmaMCPath = "../rootFiles/mcFiles/JpsiLambda/JpsiSigma";
	TFile *mcFileIn_nonZero_Sigma = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_nonZeroTracks.root",
	                                          sigmaMCPath,run));
	TTree *mcTreeIn_nonZero_Sigma = (TTree*)mcFileIn_nonZero_Sigma->Get("MyTuple");

	TFile *mcFileIn_Zero_Sigma    = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_ZeroTracks.root",
	                                          sigmaMCPath,run));
	TTree *mcTreeIn_Zero_Sigma    = (TTree*)mcFileIn_Zero_Sigma->Get("MyTuple");

	mcTreeIn_nonZero_Sigma->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_LL_FinalBDT%d_iso%d_%s.root",
	                                                 sigmaMCPath,run,finalBDTConf_nonZero,isoConf,isoVersion));
	mcTreeIn_Zero_Sigma->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_zeroTracksLL_FinalBDT%d.root",
	                                              sigmaMCPath,run,finalBDTConf_Zero));
	fstream genFile_Sigma;
	genFile_Sigma.open((Form("../logs/mc/JpsiLambda/JpsiSigma/run%d/gen_log.txt",run)));

	genFile_Sigma>>nGen_Sigma;

	fstream genEffFile_Sigma;
	genEffFile_Sigma.open(Form("../logs/mc/JpsiLambda/JpsiSigma/run%d/Generator_Effs_Combined.txt",run));

	genEffFile_Sigma>>eff_Sigma_gen;
	genEffFile_Sigma>>eff_Sigma_gen_err;

	cout<<"Sigma Generator Effs = "<<eff_Sigma_gen*100<<" % +/- "<<eff_Sigma_gen_err*100<<" %"<<endl;

	Int_t num_Sigma = mcTreeIn_nonZero_Sigma->GetEntries(Form("BDT%d > %f", finalBDTConf_nonZero,bdtCut_nonZero)) +
	                  mcTreeIn_Zero_Sigma->GetEntries(Form("BDT%d > %f", finalBDTConf_Zero,bdtCut_Zero));

	eff_Sigma_rec     = num_Sigma*1.0/nGen_Sigma;
	// eff_Sigma_rec_err = eff_Sigma_rec*0.02;

	cout<<"Sigma Recons. Effs = "<<eff_Sigma_rec*100<<" %"<<endl;

	eff_Sigma = eff_Sigma_gen * eff_Sigma_rec;
	eff_Sigma_err = eff_Sigma* 0.02;// nominally assigning 2% uncertainty for now
	// eff_Sigma_err = sqrt(pow(eff_Sigma_gen_err,2) + pow(eff_Sigma_rec_err,2)); //Combining errors in quadrature

	cout<<"Jpsi Sigma Eff = "<<eff_Sigma*100<<" % +/- "<<eff_Sigma_err*100<<" %"<<endl;

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

	//******************Get shape from Xib background********************
	// RooRealVar xibmass("Lb_DTF_M_JpsiLConstr","xibmass",5200.,5740.);
	const char* xibPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/JpsiXi";
	TFile *filein_xi_nonZero = Open(Form("%s/run%d/jpsixi_cutoutks_LL_nonZeroTracks.root",xibPath,run));
	TTree *treein_xi_nonZero = (TTree*)filein_xi_nonZero->Get("MyTuple");

	TFile *filein_xi_Zero = Open(Form("%s/run%d/jpsixi_cutoutks_LL_ZeroTracks.root",xibPath,run));
	TTree *treein_xi_Zero = (TTree*)filein_xi_Zero->Get("MyTuple");

	treein_xi_nonZero->AddFriend("MyTuple",Form("%s/run%d/jpsixi_LL_FinalBDT%d_iso%d_%s.root",
	                                            xibPath,run,finalBDTConf_nonZero,isoConf,isoVersion));
	treein_xi_Zero->AddFriend("MyTuple",Form("%s/run%d/jpsixi_zeroTracksLL_FinalBDT%d.root",
	                                         xibPath,run,finalBDTConf_Zero));
	// treein_xi->SetBranchStatus("*",0);
	// treein_xi->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	//RooDataSet ds_xi("ds_xi","ds_xi",treein_xi,ximass);
	treein_xi_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_nonZero(%d,%d,%d)",nbins,low,high),
	                        Form("Lb_BKGCAT==40 && BDT%d > %f",finalBDTConf_nonZero,bdtCut_nonZero),"goff");//TRUTH MATCHING HERE
	treein_xi_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib_Zero(%d,%d,%d)",nbins,low,high),
	                     Form("Lb_BKGCAT==40 && BDT%d > %f",finalBDTConf_Zero,bdtCut_Zero),"goff");//TRUTH MATCHING HERE

	TH1D *hxib_nonZero = (TH1D*)gDirectory->Get("hxib_nonZero");
	TH1D *hxib_Zero    = (TH1D*)gDirectory->Get("hxib_Zero");
	TH1D *hxib         = new TH1D("hxib","",nbins,low,high);
	hxib->Add(hxib_nonZero,hxib_Zero);

	TH1D *hxib_smooth  = (TH1D*)hxib->Clone("hxib_smooth");
	hxib_smooth->Smooth(4);//TODO TWEAK THIS

	RooDataHist *ds_xib        = new RooDataHist("ds_xib","ds_xib",*(w.var("bMass")),hxib);
	RooDataHist *ds_xib_smooth = new RooDataHist("ds_xib_smooth","ds_xib_smooth",*(w.var("bMass")),hxib_smooth);

	ds_xib_smooth->Print();

	//  RooKeysPdf xibshape("xibshape","xibshape",xibmass,ds_xib,RooKeysPdf::NoMirror);
	RooHistPdf xibshape("xibshape","xibshape",*(w.var("bMass")),*ds_xib,0);
	RooHistPdf XIB("XIB","XIB",*(w.var("bMass")),*ds_xib_smooth,0);

	RooPlot *framexib = (w.var("bMass"))->frame();
	framexib->SetTitle("J/#psi #Xi");
	ds_xib->plotOn(framexib,Name("xibdata"));
	xibshape.plotOn(framexib,Name("xibfit"),LineColor(kBlue));
	XIB.plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));
	//  xibshape1.plotOn(framexib,Name("xibfit1"),LineColor(kRed));
	// new TCanvas();
	// framexib->Draw();

	w.import(XIB);
	//*******************************************************************

	//******************Get shape from Jpsi Sigma signal*****************
	//  RooRealVar sigmamass("Lb_DTF_M_JpsiLConstr","sigmamass",5200.,5740.);

	const char* sigmaPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/mcFiles/JpsiLambda/JpsiSigma";

	TFile *filein_sigma_nonZero = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_nonZeroTracks.root",sigmaPath,run));
	TTree *treein_sigma_nonZero = (TTree*)filein_sigma_nonZero->Get("MyTuple");

	TFile *filein_sigma_Zero = Open(Form("%s/run%d/jpsisigma_cutoutks_LL_ZeroTracks.root",sigmaPath,run));
	TTree *treein_sigma_Zero = (TTree*)filein_sigma_Zero->Get("MyTuple");

	treein_sigma_nonZero->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_LL_FinalBDT%d_iso%d_%s.root",
	                                               sigmaPath,run,finalBDTConf_nonZero,isoConf,isoVersion));
	treein_sigma_Zero->AddFriend("MyTuple",Form("%s/run%d/jpsisigma_zeroTracksLL_FinalBDT%d.root",
	                                            sigmaPath,run,finalBDTConf_Zero));
	// treein_sigma->SetBranchStatus("*",0);
	// treein_sigma->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

	//RooDataSet ds_sigma("ds_sigma","ds_sigma",treein_sigma,sigmamass);
	treein_sigma_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hsigma_nonZero(%d,%d,%d)",nbins,low,high),
	                           Form("BDT%d > %f",finalBDTConf_nonZero,bdtCut_nonZero),"goff");//Not TRUTH MATCHING HERE!

	treein_sigma_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>hsigma_Zero(%d,%d,%d)",nbins,low,high),
	                        Form("BDT%d > %f",finalBDTConf_Zero,bdtCut_Zero),"goff");//Not TRUTH MATCHING HERE!

	TH1D *hsigma_nonZero = (TH1D*)gDirectory->Get("hsigma_nonZero");
	TH1D *hsigma_Zero    = (TH1D*)gDirectory->Get("hsigma_Zero");
	TH1D *hsigma         = new TH1D("hsigma","",nbins,low,high);

	hsigma->Add(hsigma_nonZero,hsigma_Zero);
	TH1D *hsigma_smooth  = (TH1D*)hsigma->Clone("hsigma_smooth");
	hsigma_smooth->Smooth(2);

	RooDataHist *ds_sigma = new RooDataHist("ds_sigma","ds_sigma",*(w.var("bMass")),hsigma);
	RooDataHist *ds_sigma_smooth = new RooDataHist("ds_sigma_smooth","ds_sigma_smooth",*(w.var("bMass")),hsigma_smooth);

	ds_sigma->Print();

	//  RooKeysPdf sigmashape("sigmashape","sigmashape",sigmamass,ds_sigma,RooKeysPdf::NoMirror);
	//  RooKeysPdf sigmashape("sigmashape","sigmashape",Lb_DTF_M_JpsiLConstr,ds_sigma,RooKeysPdf::NoMirror);
	//sigmashape1("sigmashape1","sigmashape1",sigmamass,ds_sigma,RooKeysPdf::MirrorBoth,2) ;
	RooHistPdf sigmashape("sigmashape","sigmashape",*(w.var("bMass")),*ds_sigma,0);
	RooHistPdf SIG("SIG","SIG",*(w.var("bMass")),*ds_sigma_smooth,0);

	RooPlot *framesigma = (w.var("bMass"))->frame();
	framesigma->SetTitle("J/#psi #Sigma");
	ds_sigma->plotOn(framesigma,Name("sigmadata"));
	sigmashape.plotOn(framesigma,Name("sigmafit"),LineColor(kBlue));
	SIG.plotOn(framesigma,Name("sigmafitsmooth"),LineColor(kRed),LineStyle(kDashed));
	//  sigmashape1.plotOn(framesigma,Name("sigmafit1"),LineColor(kRed));

	// new TCanvas();
	// framesigma->Draw();
	w.import(SIG);
	//*******************************************************************

	//*********Double Crystal Ball signal shape for Lambda_b0************
	//TODO ADD CONSTRAINTS
	w.factory("RooCBShape::Lb1(bMass,mean[5619.6,5619,5621],"
	          "sigma[10.,0.,20.], alpha1[1.028, 0.8,1.3], 10.0)" );
	// w.factory("RooCBShape::Lb1(bMass,mean[5619.86],"
	//           "sigma[6.641], alpha1[1.0529], 10.0)" );
	w.factory("RooCBShape::Lb2(bMass,mean,sigma,alpha2[-1.097,-1.4,-0.7], 10.0)");
	w.factory("SUM::Lb(0.5*Lb1 , 0.5*Lb2)");
	//*******************************************************************

	//*********Exponential shape for continuum backgruond****************
	w.factory("Exponential::Bkg(bMass,tau[-0.0007,-0.01,-0.0000001])");
	//*******************************************************************

	//*********Gaussian Lump for misc. Lambda*'s ************************
	w.factory("Gaussian::lstLump(bMass,miscLstMean[5009.,4980.,5050.],"
	          "miscLstSigma[60.,20.,100.])");
	//*******************************************************************

	//*********Input Data************************************************
	const char *dataPath = "/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/JpsiLambda";

	TFile *filein_nonZero = Open(Form("%s/run%d/"
	                                  "jpsilambda_cutoutks_LL_nonZeroTracks.root",dataPath,run),"READ");
	TTree *treein_nonZero = (TTree*)filein_nonZero->Get("MyTuple");

	TFile *filein_Zero = Open(Form("%s/run%d/"
	                               "jpsilambda_cutoutks_LL_ZeroTracks.root",dataPath,run),"READ");
	TTree *treein_Zero = (TTree*)filein_Zero->Get("MyTuple");

	treein_nonZero->AddFriend("MyTuple",Form("%s/run%d/"
	                                         "jpsilambda_LL_FinalBDT%d_iso%d_%s.root",
	                                         dataPath,run,finalBDTConf_nonZero,isoConf,isoVersion));
	treein_Zero->AddFriend("MyTuple",Form("%s/run%d/"
	                                      "jpsilambda_zeroTracksLL_FinalBDT%d.root",
	                                      dataPath,run,finalBDTConf_Zero));

	// Int_t nentries = treein_nonZero->GetEntries() + treein_Zero->GetEntries();

	treein_nonZero->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_nonzero(%d,%d,%d)",nbins,low,high),
	                     Form("BDT%d > %f",finalBDTConf_nonZero,bdtCut_nonZero),"goff");
	treein_Zero->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist_zero(%d,%d,%d)",nbins,low,high),
	                  Form("BDT%d > %f",finalBDTConf_Zero,bdtCut_Zero),"goff");
	//  treein->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist1(%d,%d,%d)",nbins,low,high),"BDT1 > 0.065","goff");

	TH1D *myhist_nonzero = (TH1D*)gDirectory->Get("myhist_nonzero");
	TH1D *myhist_zero = (TH1D*)gDirectory->Get("myhist_zero");

	TH1D *myhist = new TH1D("myhist","",nbins,low,high);

	myhist->Add(myhist_zero,myhist_nonzero);

	Int_t nentries = myhist->Integral();
	cout<<"nentries = "<<nentries<<endl;

	RooDataHist ds("ds","ds",*(w.var("bMass")),myhist);
	//  RooDataSet ds("ds","ds",treein,Lb_DTF_M_JpsiLConstr);

	ds.Print();

	w.import(ds);
	//*******************************************************************

	//************************MAKE COMBINED MODEL************************
	w.factory("logR[-6.0,-10.0,0.0]"); //logR is the parameter of interest
	w.var("logR")->setError(0.1);
	w.factory(Form("sigmaEff[%f,0.00001,0.005]",eff_Sigma));
	w.factory(Form("lambdaEff[%f,0.00001,0.005]",eff_Lambda));
	w.factory(Form("Gaussian::sigmaEff_constraint(gsigmaEff[%f,0.00001,0.005],sigmaEff,%f)",eff_Sigma,eff_Sigma_err));
	w.factory(Form("Gaussian::lambdaEff_constraint(glambdaEff[%f,0.00001,0.005],lambdaEff,%f)",eff_Lambda,eff_Lambda_err));

	w.var("gsigmaEff")->setConstant();
	w.var("glambdaEff")->setConstant();

	w.factory("nLb[0,17000]");
	w.factory(Form("nXib[%f,%f,%f]",xibnorm_LL,xibnorm_LL-(xibnorm_LL_err*2),xibnorm_LL+(xibnorm_LL_err*2)));
	w.factory(Form("Gaussian::nXib_constraint(gnXib[%f,%f,%f],nXib,%f)",xibnorm_LL,xibnorm_LL-(xibnorm_LL_err*2),xibnorm_LL+(xibnorm_LL_err*2),xibnorm_LL_err));

	w.var("gnXib")->setConstant();

	w.factory("expr::nSigma('pow(10,logR)*nLb*sigmaEff/(lambdaEff)',logR,nLb,sigmaEff,lambdaEff)");
	// w.factory("expr::nSigma('pow(10,logR)*nLb',logR,nLb)");//TAKING EFFICIENCIES OUT OF PICTURE FOR NOW. TESTING.

	w.factory(Form("SUM:model(nSigma*SIG , nLb*Lb , nXib*XIB , n1405[1500,1,5000]*LST1405 , n1520[3000,1,5000]*LST1520 , nMiscLst[2500,1,5000]*lstLump , nBkg[2000,1,%d]*Bkg)",nentries));
	if(!lst1520flag)
	{
		w.var("n1520")->setVal(0);
		w.var("n1520")->setConstant();
	}
	w.factory("PROD::model_const(model,sigmaEff_constraint,lambdaEff_constraint,nXib_constraint)");//Multiply model by constraint terms

	RooAbsPdf* model_const = w.pdf("model_const"); // get the model
	w.defineSet("poi","logR"); //parameters of interest
	w.defineSet("obs","bMass");//not sure if bMass should be observable or total number of events?
	w.defineSet("nuisParams","nLb,mean,sigma,alpha1,alpha2,nBkg,tau,nMiscLst,miscLstMean,miscLstSigma,sigmaEff,lambdaEff,nXib,n1405,n1520");// define set of nuisance parameters
	w.defineSet("globObs","gsigmaEff,glambdaEff,gnXib");//define set of global observables
	//*******************************************************************

	//***********************MAKE NLL************************************
	// RooAbsReal* nll = w.pdf("model_const")->createNLL(ds);
	//
	// new TCanvas();
	// RooPlot* frame = w.var("logR")->frame();
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
	// cout<<"MIN "<<w.var("logR")->getVal()<<" ERR LO"<<w.var("logR")->getErrorLo()<<" ERR HI"<<w.var("logR")->getErrorHi()<<endl;
	//
	// RooFormulaVar p_mu("p_mu","p_{#mu} using asymptotic formulae","TMath::Prob(2*@0,1.)",RooArgList(*pll));
	//
	// TCanvas *c2 = new TCanvas();
	// RooPlot* frame1 = w.var("logR")->frame();
	// p_mu.plotOn(frame1,LineColor(kGreen));
	// // frame1->Draw();
	//
	// c2->cd();
	//
	// Double_t limit = p_mu.findRoot(*(w.var("logR")),-5.0,0.0,0.05);
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
	// Double_t limit1 = p_mu_q.findRoot(*(w.var("logR")),-5.0,0.0,0.05);
	// cout<<"One sided UL ="<<limit1<<endl;
	//
	// yields.push_back(0.0);
	// return yields;


	//************************DO THE FIT************************
/*	RooFitResult *res = model_const->fitTo(ds,Minos(*w.set("poi")),Extended(), Save(), Hesse(false), Strategy(2));
        //*******************************************************************

        //*********************PLOTTING STUFF*********************************************



        //  RooPlot *frame = Lb_DTF_M_JpsiLConstr.frame();
        RooPlot *frame = new RooPlot(*(w.var("bMass")), low,5800,100);
        //  frame->GetYaxis()->SetRangeUser(0,100);
        ds.plotOn(frame,Name("data"));
        model_const->plotOn(frame,Name("fit"));
        model_const->plotOn(frame,Components(*(w.pdf("Lb"))),Name("sig"),LineColor(kMagenta+2));
        model_const->plotOn(frame,Components(*(w.pdf("Lb1"))),LineStyle(kDotted),LineColor(kMagenta));
        model_const->plotOn(frame,Components(*(w.pdf("Lb2"))),LineStyle(kDotted),LineColor(kMagenta));
        model_const->plotOn(frame,Components(*(w.pdf("Bkg"))),LineColor(kRed),Name("bkg"));
        //  if(xibflag!=0)
        model_const->plotOn(frame,Components(*(w.pdf("XIB"))),LineColor(kGreen),Name("xib"));
        if(lst1405flag!=0)
                model_const->plotOn(frame,Components(*(w.pdf("LST1405"))),LineColor(kGreen+2),LineStyle(kDashed),Name("lst1405"));
        if(lst1520flag!=0)
                model_const->plotOn(frame,Components(*(w.pdf("LST1520"))),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1520"));
        //model_const->plotOn(frame,Components(lst1810shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1810"));
        // if(sigmaflag!=0)
        model_const->plotOn(frame,Components(*(w.pdf("SIG"))),LineColor(kBlack),Name("sigma"));
        model_const->plotOn(frame,Components(*(w.pdf("lstLump"))),LineColor(kRed+2),LineStyle(kDashed),Name("misclst"));

        Double_t chiSquare1 = frame->chiSquare("fit","data");
        cout<<"chi square1/dof = "<<chiSquare1<<endl;
        RooArgSet *floatpar = model_const->getParameters(ds);
        floatpar->Print();
        int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
        cout<<"float pars = "<<floatpars<<endl;
        Double_t chi2 = frame->chiSquare("fit","data",floatpars);
        cout<<"chi square2/dof = "<<chi2<<endl;

        Int_t fit_ndof = nbins - floatpars;
        cout<<"chi square2 = "<<chi2*fit_ndof<<endl;

        TCanvas* c = new TCanvas("Ksfit","Ksfit", 1200, 800);

        // ///////////
        TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
        TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
        //  pad1->SetBottomMargin(0);
        pad2->SetTopMargin(0);
        pad2->SetBottomMargin(0.5);
        pad2->SetBorderMode(0);
        pad1->SetBorderMode(0);
        c->SetBorderMode(0);
        pad2->Draw();
        pad1->Draw();
        pad1->cd();

        gPad->SetTopMargin(0.06);
        pad1->Update();

        frame->Draw();

        //  TPaveText *pt = (TPaveText*)pad1->GetListOfPrimitives()->FindObject("model_paramBox");

        // pt->AddText(Form("#chi^{2}/dof = %f",chi2));

        c->Modified();
        // Pull distribution
        //  RooPlot *framex2 = Lb_DTF_M_JpsiLConstr.frame();
        RooPlot *framex2 = new RooPlot(*(w.var("bMass")), low,5800,100);
        RooHist* hpull = frame->pullHist("data","fit");
        framex2->addPlotable(hpull,"P");
        hpull->SetLineColor(kBlack);
        hpull->SetMarkerColor(kBlack);
        framex2->SetTitle(0);
        framex2->GetYaxis()->SetTitle("Pull");
        framex2->GetYaxis()->SetTitleSize(0.15);
        framex2->GetYaxis()->SetLabelSize(0.15);
        framex2->GetXaxis()->SetTitleSize(0.2);
        framex2->GetXaxis()->SetLabelSize(0.15);
        framex2->GetYaxis()->CenterTitle();
        framex2->GetYaxis()->SetTitleOffset(0.25);
        framex2->GetXaxis()->SetTitleOffset(1.1);
        framex2->GetYaxis()->SetNdivisions(505);
        framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
        pad2->cd();
        framex2->Draw();

        c->cd();
        pad1->cd();

        cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
        cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;

        auto legend = new TLegend(0.5,0.7,0.9,0.9);
        legend->AddEntry("data","Data","lp");
        legend->AddEntry("fit","Total Fit","l");
        legend->AddEntry("sigma","J/#psi #Lambda shape","l");
        legend->AddEntry("sigma","J/#psi #Sigma shape","l");
        legend->AddEntry("xib","J/#psi #Xi shape","l");
        if(lst1405flag)
        {
                legend->AddEntry("lst1405","J/#psi #Lambda(1405) shape","l");
        }
        if(lst1520flag)
        {
                legend->AddEntry("lst1520","J/#psi #Lambda(1520) shape","l");
        }
        // if(lst1810flag)
        //      legend->AddEntry("lst1810","J/#psi #Lambda(1810) shape","l");
        legend->AddEntry("misclst","misc. J/#psi #Lambda* shapes","l");
        legend->AddEntry("bkg","Comb. Bkg. shape","l");
        legend->Draw("same");

        c->Update();
 */
	//*************ROOSTATS MODEL CONFIG*********************************
	RooStats::ModelConfig mc("ModelConfig",&w);
	mc.SetPdf(*model_const);
	mc.SetParametersOfInterest(*w.set("poi"));
	mc.SetObservables(*w.set("obs"));
	mc.SetGlobalObservables(*w.set("globObs"));
	mc.SetNuisanceParameters(*w.set("nuisParams"));

	// import model in the workspace
	w.import(mc);
	//*******************************************************************

	//***************ASIMOV DATASET*************************************
	RooDataSet *asimovData = (RooDataSet*)RooStats::AsymptoticCalculator::GenerateAsimovData(*(mc.GetPdf()),*(mc.GetObservables()));
	RooPlot *frame_asim = new RooPlot(*(w.var("bMass")), low,5800,100);

	TCanvas *asim = new TCanvas();
	asimovData->plotOn(frame_asim);
	frame_asim->Draw();

	mc.SetSnapshot(*(w.var("logR")));

	RooStats::ModelConfig *bkgOnlyModel = mc.Clone();
	bkgOnlyModel->SetName("bkgOnlyModel");
	Double_t oldval = w.var("logR")->getVal();
	w.var("logR")->setVal(-10);//What value should I set for the background only hypothesis? log of 0 is -inf
	bkgOnlyModel->SetSnapshot(*(w.var("logR")));

	w.var("logR")->setVal(oldval);

	w.import(*bkgOnlyModel);
	// w.Write();
	//*******************************************************************

	//**************Fit like Matt****************************************
	RooAbsReal* nll = w.pdf("model_const")->createNLL(ds);
	RooMinimizer min_nll(*nll);
	min_nll.setErrorLevel(0.5);
	min_nll.setStrategy(2);

	Double_t initnll = nll->getVal();

	RooArgSet *pars = w.pdf("model_const")->getParameters( ds );
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

	min_nll.minos(*(w.var("logR")));

	RooFitResult *res = min_nll.save();

	Double_t finalnll = nll->getVal();

	Double_t origval = w.var("logR")->getVal();
	Double_t errhi   = w.var("logR")->getErrorHi();
	Double_t errlo   = w.var("logR")->getErrorLo();

	//**************Upper Limit Calculation****************************************
	RooStats::AsymptoticCalculator hc(ds, *bkgOnlyModel, mc, true);
	hc.SetOneSided(true);
	RooStats::HypoTestInverter calc(hc);

	calc.SetConfidenceLevel(0.95);
	calc.SetFixedScan(100,-7, -1);
	calc.UseCLs(true);
	calc.SetVerbose(true);
	HypoTestInverterResult *r = calc.GetInterval();

	cout<<"origval = "<<origval<<" error high = "
	    <<errhi<<" error low = "<<errlo<<endl;
	cout<<"UL = "<<r->UpperLimit()<<" +/- "<<r->UpperLimitEstimatedError()<<endl;

	new TCanvas();
	RooStats::HypoTestInverterPlot plot("hti_plot", "hti_plot", r);
	plot.Draw("CLb 2CL");


	// write the workspace in the file
	// TString fileName = "MyModelConfig.root";
	// w.writeToFile(fileName,true);
	// cout << "workspace written to file " << fileName << endl;



	yields.push_back(0.0);
	return yields;

}
