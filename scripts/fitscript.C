#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TLegend.h>
#include <RooGaussian.h>
#include <RooCBShape.h>
#include <RooPlot.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
//#include <RooChebychev.h>
#include <RooExponential.h>
#include <RooAddPdf.h>
#include <RooProdPdf.h>
#include <RooKeysPdf.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <RooHist.h>
#include <TAxis.h>
#include <TH1D.h>
#include <TString.h>
#include <RooDataHist.h>
#include <TError.h>
#include <RooFitResult.h>
#include <RooHistPdf.h>
#include <RooDataHist.h>
#include <RooConstVar.h>
using namespace RooFit;
using namespace std;

void fitscript(Int_t lst1405flag = 1, Int_t lst1520flag = 1, Int_t lst1810flag = 1, Int_t xibflag = 1, Int_t sigmaflag = 1, Int_t type = 1, Int_t Lst1405_rwtype = 2, Int_t mylow = 4700, Int_t constraintflag = 1)
{/*Set LstXXXXflag = 1 to include the shapes for these decays in the final fit
	Set xibflag = 1 to include Xib background shape.
	Set sigmaflag = 1 to include J/psi Sigma signla shape.
	Set type = 1 to fit LL data, type = 2 to fit DD data.
	Set Lst1405_rwtype = 0 to use raw Lst(1405) shape, Lst1405_rwtype = 1 for MV and 2 for BONN reweighting schemes respectively
	mylow determines the lower limit of the fit range
	Set constraintflag = 0 to fit with Gaussian constraint only on Xib yield, constraintflag = 1 includes(on top of Xib constraint), a Gaussian constraint on the Lb mean */
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	Int_t low = mylow, high = 7000; //Define range in which fit is performed

	Double_t xibminus_norm_LL = 42.5; //37.5;//comes from fit to xib data
	Double_t xibminus_norm_LL_stat = 8.2; //7.7;//statistical error, omes from fit to xib data
	Double_t xibminus_norm_LL_syst = 0.02*xibminus_norm_LL; //conservative 2% systematic, based on f(Xib-) analysis //=0.85
	Double_t eff_xib_jpsixi_LL = 0.030734; //add errors for these effs
	Double_t eff_xib_jpsil_LL = 0.024539;

	Double_t xibminus_norm_LL_err = sqrt(pow(xibminus_norm_LL_stat,2) + pow(xibminus_norm_LL_syst,2) ); //adding stat. and syst. errors in quadrature // = 8.244

	Double_t xibnorm_LL = (xibminus_norm_LL/eff_xib_jpsixi_LL)*eff_xib_jpsil_LL*2.0; //2.0 to account for xib0 as well  //=67.8667
	Double_t xibnorm_LL_err = (xibminus_norm_LL_err/eff_xib_jpsixi_LL)*eff_xib_jpsil_LL*sqrt(2.0); //sqrt because I'm adding the errors in quadrature // = 9.309

	cout<<"The LL Xib normalization is "<<xibnorm_LL<<" +/- "<<xibnorm_LL_err<<endl;

	RooRealVar Lb_DTF_M_JpsiLConstr("Lb_DTF_M_JpsiLConstr","Fit to J/#psi #Lambda inv. mass",low, high,"MeV"); //MASTER VARIABLE

	Int_t binwidth = 4; //Final fit is a binned fit with this binning

	Int_t nbins = (Int_t)(high-low)/binwidth;
	Lb_DTF_M_JpsiLConstr.setBins(nbins);

	//*********Double Gaussian signal shape for Lambda_b0**************************************************
	/*  //  RooRealVar mean("mean","Gaussian Mean",5619.0,5500.0,5700.0);
	    RooRealVar mean("mean","Gaussian Mean",5619.0,5618.0,5620.0,"MeV");//pretty tight
	    RooRealVar sigma1("sigma1","Gaussian sigma1",7.5,1.0,50.0,"MeV");
	    RooRealVar sigma2("sigma2","Gaussian sigma2",10.0,1.0,50.0,"MeV");

	    RooGaussian sig1("sig1","Gaussian signal1",Lb_DTF_M_JpsiLConstr,mean,sigma1);//Gaussian 1
	    RooGaussian sig2("sig2","Gaussian signal2",Lb_DTF_M_JpsiLConstr,mean,sigma2);//Gaussian 2

	    RooRealVar frac1("frac1","Fraction of sig1 in signal",0.3,0.1,1.0);//Relative fraction of Gaussian 1 in total Gaussian

	    RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);//total double gaussian pdf*/
	//*****************************************************************************************************

	//*********Double Crystal Ball signal shape for Lambda_b0 (pars extracted from simulation)**************************************************
	RooRealVar mean("mean","Crystal Ball Mean",5619.6,5618,5622);
	RooRealVar alpha1("alpha1","alpha1",1.028, 0.8,1.3);
	RooRealVar alpha2("alpha2","alpha2",-1.097,-1.4,-0.8);
	RooRealVar sigma("sigma","Crystal Ball sigma",7.069);//,6.4,7.8);
	RooRealVar CBn("CBn","Crystal Ball n",10.0);//fixed in both sim fit and here

	//RooGaussian mean_constraint("mean_constraint","CB mean constraint",mean,RooFit::RooConst(5620.09),RooFit::RooConst(0.107));
	RooGaussian meanpdg_constraint("mean_constraint","CB mean constraint",mean,RooFit::RooConst(5619.6),RooFit::RooConst(0.17));
	// RooGaussian alpha1_constraint("alpha1_constraint","CB alpha1 constraint",mean,RooFit::RooConst(1.028),RooFit::RooConst(0.041));
	// RooGaussian alpha2_constraint("alpha2_constraint","CB alpha2 constraint",mean,RooFit::RooConst(-1.097),RooFit::RooConst(0.046));
	// RooGaussian sigma_constraint("sigma_constraint","CB sigma constraint",mean,RooFit::RooConst(7.069),RooFit::RooConst(0.132));

	RooRealVar fs("fs","Scale factor for width",1.0,0.5,2.0);
	RooFormulaVar sigma1("sigma1","sigma1","sigma*fs",RooArgList(sigma,fs));

	RooCBShape sig1("sig1","Crystal Ball 1",Lb_DTF_M_JpsiLConstr,mean,sigma1,alpha1,CBn);
	RooCBShape sig2("sig2","Crystal Ball 2",Lb_DTF_M_JpsiLConstr,mean,sigma1,alpha2,CBn);

	RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5);//fixed in both sim fit and here
	RooAddPdf sig("sig","Total CB signal",RooArgList(sig1,sig2),frac1);
	//*********Exponential shape for continuum backgruond**************************************************
	RooRealVar tau("tau","tau",-0.0007,-0.01,-0.0000001);
	RooExponential bkg("bkg","Exponential bkg",Lb_DTF_M_JpsiLConstr,tau);

	//*****************************************************************************************************

	//********Gaussian Lump for misc. Lambda*'s ***********************************************************
	//  RooRealVar misclstmass("misclstmass","misclst MASS",4500.,5500.);
	RooRealVar misclstmean("misclstmean","Misc. Lst. Gaussian Mean",5000.,4950.,5050.,"MeV");
	RooRealVar misclstsigma("misclstsigma","Misc. Lst. Gaussian Sigma",60.,20.,100.,"MeV");
	RooGaussian misclstshape("misclstshape","Misc #Lambda^{*}",Lb_DTF_M_JpsiLConstr,misclstmean,misclstsigma);

	//*********Get shape from Lambda*(1405) background******************************************************
	//  RooRealVar lstmass("JpsiLmass","JpsiL mass",4200.,5500.);
	//  RooRealVar MVweight("MVweight","MVweight",0.,10.);
	TFile *filein_lst1405;
	TTree *treein_lst1405 = (TTree*)malloc(sizeof(*treein_lst1405));
	TString lst1405wtexp = "";

	if(Lst1405_rwtype == 0) {
		filein_lst1405 = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1405/Lst1405_total_MVrw.root");
		treein_lst1405 = (TTree*)filein_lst1405->Get("MCDecayTree");
		treein_lst1405->SetBranchStatus("*",0);
		treein_lst1405->SetBranchStatus("JpsiLmass",1);
		lst1405wtexp = "";

	}
	else if(Lst1405_rwtype == 1) {
		filein_lst1405 = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1405/Lst1405_total_MVrw.root");
		lst1405wtexp = "MVweight";
		treein_lst1405 = (TTree*)filein_lst1405->Get("MCDecayTree");
		treein_lst1405->SetBranchStatus("*",0);
		treein_lst1405->SetBranchStatus("JpsiLmass",1);
		treein_lst1405->SetBranchStatus("MVweight",1);
	}
	else if(Lst1405_rwtype == 2) {
		filein_lst1405 = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1405/Lst1405_total_BONNrw.root");
		lst1405wtexp = "BONNweight";
		treein_lst1405 = (TTree*)filein_lst1405->Get("MCDecayTree");
		treein_lst1405->SetBranchStatus("*",0);
		treein_lst1405->SetBranchStatus("JpsiLmass",1);
		treein_lst1405->SetBranchStatus("BONNweight",1);
	}

	//  treein_lst->SetAlias("Lb_DTF_M_JpsiLConstr","JpsiLmass");

	//  RooDataSet ds_lst("ds_lst","ds_lst",treein_lst,RooArgSet(lstmass,MVweight),"MVweight");
	treein_lst1405->Draw(Form("JpsiLmass>>hlst1405(%d,%d,%d)",nbins,low,high),lst1405wtexp,"goff");

	TH1D *hlst1405 = (TH1D*)gDirectory->Get("hlst1405");
	TH1D *hlst1405_smooth = (TH1D*)hlst1405->Clone("hlst1405_smooth");
	hlst1405_smooth->Smooth(2);

	//  RooDataHist *ds_lst1 = ds_lst.binnedClone();
	RooDataHist *ds_lst1405 = new RooDataHist("ds_lst1405","ds_lst1405",Lb_DTF_M_JpsiLConstr,hlst1405);
	RooDataHist *ds_lst1405_smooth = new RooDataHist("ds_lst1405_smooth","ds_lst1405_smooth",Lb_DTF_M_JpsiLConstr,hlst1405_smooth);

	ds_lst1405->Print();

	//  RooKeysPdf lstshape("lstshape","lstshape",lstmass,ds_lst);
	RooHistPdf lst1405shape("lst1405shape","lst1405shape",Lb_DTF_M_JpsiLConstr,*ds_lst1405,0);
	RooHistPdf lst1405shape_smooth("lst1405shape_smooth","lst1405shape_smooth",Lb_DTF_M_JpsiLConstr,*ds_lst1405_smooth,0);

	RooPlot *framelst1405 = Lb_DTF_M_JpsiLConstr.frame();
	framelst1405->SetTitle("J/#psi Lambda^{*} (1405)");
	ds_lst1405->plotOn(framelst1405,Name("lst1405data"));
	lst1405shape.plotOn(framelst1405,Name("lst1405fit"),LineColor(kBlue));
	lst1405shape_smooth.plotOn(framelst1405,Name("lst1405fitsmooth"),LineColor(kRed),LineStyle(kDashed));
	new TCanvas();
	framelst1405->Draw();
	//*****************************************************************************************************

	//*********Get shape from Lambda*(1520) background******************************************************
	//  RooRealVar lstmass("JpsiLmass","JpsiL mass",4200.,5500.);
	//  RooRealVar MVweight("MVweight","MVweight",0.,10.);
	TFile *filein_lst1520;
	TTree *treein_lst1520 = (TTree*)malloc(sizeof(*treein_lst1520));

	filein_lst1520 = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1520/Lst1520_total_jpsil.root");
	treein_lst1520 = (TTree*)filein_lst1520->Get("MCDecayTree");
	treein_lst1520->SetBranchStatus("*",0);
	treein_lst1520->SetBranchStatus("JpsiLmass",1);

	//  treein_lst->SetAlias("Lb_DTF_M_JpsiLConstr","JpsiLmass");

	//  RooDataSet ds_lst("ds_lst","ds_lst",treein_lst,RooArgSet(lstmass,MVweight),"MVweight");
	treein_lst1520->Draw(Form("JpsiLmass>>hlst1520(%d,%d,%d)",nbins,low,high),"","goff");

	TH1D *hlst1520 = (TH1D*)gDirectory->Get("hlst1520");
	TH1D *hlst1520_smooth = (TH1D*)hlst1520->Clone("hlst1520_smooth");
	hlst1520_smooth->Smooth(2);

	//  RooDataHist *ds_lst1 = ds_lst.binnedClone();
	RooDataHist *ds_lst1520 = new RooDataHist("ds_lst1520","ds_lst1520",Lb_DTF_M_JpsiLConstr,hlst1520);
	RooDataHist *ds_lst1520_smooth = new RooDataHist("ds_lst1520_smooth","ds_lst1520_smooth",Lb_DTF_M_JpsiLConstr,hlst1520_smooth);

	ds_lst1520->Print();

	//  RooKeysPdf lstshape("lstshape","lstshape",lstmass,ds_lst);
	RooHistPdf lst1520shape("lst1520shape","lst1520shape",Lb_DTF_M_JpsiLConstr,*ds_lst1520,0);
	RooHistPdf lst1520shape_smooth("lst1520shape_smooth","lst1520shape_smooth",Lb_DTF_M_JpsiLConstr,*ds_lst1520_smooth,0);

	RooPlot *framelst1520 = Lb_DTF_M_JpsiLConstr.frame();
	framelst1520->SetTitle("J/#psi Lambda^{*} (1520)");
	ds_lst1520->plotOn(framelst1520,Name("lst1520data"));
	lst1520shape.plotOn(framelst1520,Name("lst1520fit"),LineColor(kBlue));
	lst1520shape_smooth.plotOn(framelst1520,Name("lst1520fitsmooth"),LineColor(kRed),LineStyle(kDashed));
	new TCanvas();
	framelst1520->Draw();
	//*****************************************************************************************************

	//*********Get shape from Lambda*(1810) background******************************************************
	//  RooRealVar lstmass("JpsiLmass","JpsiL mass",4200.,5500.);
	//  RooRealVar MVweight("MVweight","MVweight",0.,10.);
	TFile *filein_lst1810;
	TTree *treein_lst1810 = (TTree*)malloc(sizeof(*treein_lst1810));

	filein_lst1810 = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/Lst1810/Lst1810_total_jpsil.root");
	treein_lst1810 = (TTree*)filein_lst1810->Get("MCDecayTree");
	treein_lst1810->SetBranchStatus("*",0);
	treein_lst1810->SetBranchStatus("JpsiLmass",1);

	//  treein_lst->SetAlias("Lb_DTF_M_JpsiLConstr","JpsiLmass");

	//  RooDataSet ds_lst("ds_lst","ds_lst",treein_lst,RooArgSet(lstmass,MVweight),"MVweight");
	treein_lst1810->Draw(Form("JpsiLmass>>hlst1810(%d,%d,%d)",nbins,low,high),"","goff");

	TH1D *hlst1810 = (TH1D*)gDirectory->Get("hlst1810");
	TH1D *hlst1810_smooth = (TH1D*)hlst1810->Clone("hlst1810_smooth");
	hlst1810_smooth->Smooth(2);

	//  RooDataHist *ds_lst1 = ds_lst.binnedClone();
	RooDataHist *ds_lst1810 = new RooDataHist("ds_lst1810","ds_lst1810",Lb_DTF_M_JpsiLConstr,hlst1810);
	RooDataHist *ds_lst1810_smooth = new RooDataHist("ds_lst1810_smooth","ds_lst1810_smooth",Lb_DTF_M_JpsiLConstr,hlst1810_smooth);

	ds_lst1810->Print();

	//  RooKeysPdf lstshape("lstshape","lstshape",lstmass,ds_lst);
	RooHistPdf lst1810shape("lst1810shape","lst1810shape",Lb_DTF_M_JpsiLConstr,*ds_lst1810,0);
	RooHistPdf lst1810shape_smooth("lst1810shape_smooth","lst1810shape_smooth",Lb_DTF_M_JpsiLConstr,*ds_lst1810_smooth,0);

	RooPlot *framelst1810 = Lb_DTF_M_JpsiLConstr.frame();
	framelst1810->SetTitle("J/#psi Lambda^{*} (1810)");
	ds_lst1810->plotOn(framelst1810,Name("lst1810data"));
	lst1810shape.plotOn(framelst1810,Name("lst1810fit"),LineColor(kBlue));
	lst1810shape_smooth.plotOn(framelst1810,Name("lst1810fitsmooth"),LineColor(kRed),LineStyle(kDashed));
	new TCanvas();
	framelst1810->Draw();
	//*****************************************************************************************************

	//******************Get shape from Xib background******************************************************
	// RooRealVar xibmass("Lb_DTF_M_JpsiLConstr","xibmass",5200.,5740.);

	TFile *filein_xib = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/JpsiXi/total_run1/jpsixi_DD.root");
	TTree *treein_xib = (TTree*)filein_xib->Get("MyTuple");

	treein_xib->SetBranchStatus("*",0);
	treein_xib->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

	//  RooDataSet ds_xib("ds_xib","ds_xib",treein_xib,xibmass);
	treein_xib->Draw(Form("Lb_DTF_M_JpsiLConstr>>hxib(%d,%d,%d)",nbins,low,high),"","goff");

	TH1D *hxib = (TH1D*)gDirectory->Get("hxib");
	TH1D *hxib_smooth = (TH1D*)hxib->Clone("hxib_smooth");
	hxib_smooth->Smooth(2);

	RooDataHist *ds_xib = new RooDataHist("ds_xib","ds_xib",Lb_DTF_M_JpsiLConstr,hxib);
	RooDataHist *ds_xib_smooth = new RooDataHist("ds_xib_smooth","ds_xib_smooth",Lb_DTF_M_JpsiLConstr,hxib_smooth);

	ds_xib_smooth->Print();

	//  RooKeysPdf xibshape("xibshape","xibshape",xibmass,ds_xib,RooKeysPdf::NoMirror);
	RooHistPdf xibshape("xibshape","xibshape",Lb_DTF_M_JpsiLConstr,*ds_xib,0);
	RooHistPdf xibshape_smooth("xibshape_smooth","xibshape_smooth",Lb_DTF_M_JpsiLConstr,*ds_xib_smooth,0);

	RooPlot *framexib = Lb_DTF_M_JpsiLConstr.frame();
	framexib->SetTitle("J/#psi #Xi");
	ds_xib->plotOn(framexib,Name("xibdata"));
	xibshape.plotOn(framexib,Name("xibfit"),LineColor(kBlue));
	xibshape_smooth.plotOn(framexib,Name("xibfitsmooth"),LineColor(kRed),LineStyle(kDashed));
	//  xibshape1.plotOn(framexib,Name("xibfit1"),LineColor(kRed));
	new TCanvas();
	framexib->Draw();
	//****************************************************************************************************

	//******************Get shape from Jpsi Sigma signal***********************************************
	//  RooRealVar sigmamass("Lb_DTF_M_JpsiLConstr","sigmamass",5200.,5740.);

	TFile *filein_sigma = TFile::Open("/data1/avenkate/JpsiLambda_restart/mc/JpsiSigma/total_run1/jpsisigma_DD.root");//WHY AM I USING THE DD SHAPE HERE, TO FIT LL?
	TTree *treein_sigma = (TTree*)filein_sigma->Get("MyTuple");

	treein_sigma->SetBranchStatus("*",0);
	treein_sigma->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);

	//RooDataSet ds_sigma("ds_sigma","ds_sigma",treein_sigma,sigmamass);
	treein_sigma->Draw(Form("Lb_DTF_M_JpsiLConstr>>hsigma(%d,%d,%d)",nbins,low,high),"","goff");

	TH1D *hsigma = (TH1D*)gDirectory->Get("hsigma");
	TH1D *hsigma_smooth = (TH1D*)hsigma->Clone("hsigma_smooth");
	hsigma_smooth->Smooth(2);

	RooDataHist *ds_sigma = new RooDataHist("ds_sigma","ds_sigma",Lb_DTF_M_JpsiLConstr,hsigma);
	RooDataHist *ds_sigma_smooth = new RooDataHist("ds_sigma_smooth","ds_sigma_smooth",Lb_DTF_M_JpsiLConstr,hsigma_smooth);

	ds_sigma->Print();

	//  RooKeysPdf sigmashape("sigmashape","sigmashape",sigmamass,ds_sigma,RooKeysPdf::NoMirror);
	//  RooKeysPdf sigmashape("sigmashape","sigmashape",Lb_DTF_M_JpsiLConstr,ds_sigma,RooKeysPdf::NoMirror);
	//sigmashape1("sigmashape1","sigmashape1",sigmamass,ds_sigma,RooKeysPdf::MirrorBoth,2) ;
	RooHistPdf sigmashape("sigmashape","sigmashape",Lb_DTF_M_JpsiLConstr,*ds_sigma,0);
	RooHistPdf sigmashape_smooth("sigmashape_smooth","sigmashape_smooth",Lb_DTF_M_JpsiLConstr,*ds_sigma_smooth,0);

	RooPlot *framesigma = Lb_DTF_M_JpsiLConstr.frame();
	framesigma->SetTitle("J/#psi #Sigma");
	ds_sigma->plotOn(framesigma,Name("sigmadata"));
	sigmashape.plotOn(framesigma,Name("sigmafit"),LineColor(kBlue));
	sigmashape_smooth.plotOn(framesigma,Name("sigmafitsmooth"),LineColor(kRed),LineStyle(kDashed));
	//  sigmashape1.plotOn(framesigma,Name("sigmafit1"),LineColor(kRed));

	new TCanvas();
	framesigma->Draw();
	//*****************************************************************************************************

	//*********Input Data*********************************************************************************
	TFile *filein = TFile::Open("../rootFiles/dataFiles/JpsiLambda/run1/jpsilambdaLite_LL_withFinalBDT1_iso1_v0.root","READ");
	// else if(type == 2)
	//   filein = TFile::Open("../jpsilambda_DD.root","READ");

	TTree *treein = (TTree*)filein->Get("MyTuple");

	treein->SetBranchStatus("*",0);
	treein->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
	treein->SetBranchStatus("BDT1",1);

	Int_t nentries = treein->GetEntries();
	cout<<"nentries = "<<nentries<<endl;

	treein->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist(%d,%d,%d)",nbins,low,high),"","goff");
	treein->Draw(Form("Lb_DTF_M_JpsiLConstr>>myhist1(%d,%d,%d)",nbins,low,high),"BDT1 > 0.065","goff");

	TH1D *myhist = (TH1D*)gDirectory->Get("myhist");

	RooDataHist ds("ds","ds",Lb_DTF_M_JpsiLConstr,myhist); //Fit to binned histogram
	//  RooDataSet ds("ds","ds",treein,Lb_DTF_M_JpsiLConstr);

	ds.Print();
	//*****************************************************************************************************

	cout<<"stage1"<<endl;

	//*************YIELDS*********************************************************************************
	RooRealVar *nsig = (RooRealVar*)malloc(sizeof(*nsig));
	RooRealVar *nxib = (RooRealVar*)malloc(sizeof(*nxib));

	if(type == 1) {
		nsig = new RooRealVar("nsig","nsig",1,5000);
		nxib = new RooRealVar("nxib","nxib",xibnorm_LL,xibnorm_LL-xibnorm_LL_err,xibnorm_LL+xibnorm_LL_err);
		//    nxib = new RooRealVar("nxib","nxib",64);
	}
	if(type == 2) {
		nsig = new RooRealVar("nsig","nsig",0,20000);
		nxib = new RooRealVar("nxib","nxib",1114,1054,1174);
	}
	// else if(type == 2) {
	//   nsig = RooRealVar("nsig","nsig",0,20000);
	//   nxib = RooRealVar("nxib","nxib", 1114,1054,1174);
	// }
	//********Gaussian Constraint on xib yield****************
	RooGaussian nxib_constraint("nxib_constraint","Xib yield constraint",*nxib,RooFit::RooConst(xibnorm_LL),RooFit::RooConst(xibnorm_LL_err)); //find a way to put the variables in here rather than hard coding it
	//********************************************************

	RooRealVar nlst1405("nlst1405","nlst1405",1,5000);//Lst(1405) yield
	RooRealVar nlst1520("nlst1520","nlst1520",1,5000);//Lst(1520) yiel
	RooRealVar nlst1810("nlst1810","nlst1810",1,5000);//Lst(1520) yield
	RooRealVar nmisclst("nmisclst","nmisclst",1,10000);//Misc. lst yield
	RooRealVar nsigma("nsigma","nsigma",1,5000);//J/psi Sigma yield
	RooRealVar nbkg("nbkg","nbkg",1,nentries);//Comb. bkg yield
	//*****************************************************************************************************

	cout<<"stage1.1"<<endl;
	if(lst1405flag == 0) {
		nlst1405.setVal(0.0);
		nlst1405.setConstant();
	}
	if(lst1520flag == 0) {
		nlst1520.setVal(0.0);
		nlst1520.setConstant();
	}
	if(lst1810flag == 0) {
		nlst1810.setVal(0.0);
		nlst1810.setConstant();
	}
	if(xibflag == 0) {
		nxib->setVal(0.0);
		nxib->setConstant();
	}
	if(sigmaflag == 0) {
		nsigma.setVal(0.0);
		nsigma.setConstant();
	}
	cout<<"stage1.2"<<endl;

	// if(lstflag == 0 && xibflag == 0 && sigmaflag == 0)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg),RooArgList(nsig,nbkg));
	// else if(lstflag == 1 && xibflag == 0 && sigmaflag == 0)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg,lstshape),RooArgList(nsig,nbkg,nlst));
	// else if(lstflag == 0 && xibflag == 1 && sigmaflag == 0)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg,xibshape),RooArgList(nsig,nbkg,nxib));
	// else if(lstflag == 0 && xibflag == 0 && sigmaflag == 1)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg,sigmashape),RooArgList(nsig,nbkg,nsigma));
	// else if(lstflag == 0 && xibflag == 1 && sigmaflag == 1)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg,xibshape,sigmashape),RooArgList(nsig,nbkg,nxib,nsigma));
	// else if(lstflag == 1 && xibflag == 0 && sigmaflag == 1)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg,sigmashape,lstshape),RooArgList(nsig,nbkg,nsigma,nlst));
	// else if(lstflag == 1 && xibflag == 1 && sigmaflag == 0)
	//   model = RooAddPdf("model","model",RooArgList(sig,bkg,xibshape,lstshape),RooArgList(nsig,nbkg,nxib,nlst));
	// else if(lstflag == 1 && xibflag == 1 && sigmaflag == 1)
	//  RooAddPdf model("model","model",RooArgList(sig,bkg,xibshape,sigmashape,lstshape,misclstshape),RooArgList(nsig,nbkg,nxib,nsigma,nlst,nmisclst));
	// if(sigmaflag == 0)
	//   RooAddPdf model("model","model",RooArgList(sig,bkg,lstshape_smooth,misclstshape,xibshape_smooth,sigmashape_smooth),RooArgList(nsig,nbkg,nlst,nmisclst,nxib,nsigma));
	// else
	RooAddPdf model("model","model",RooArgList(sig,bkg),RooArgList(*nsig,nbkg));

	//RooProdPdf modelc("modelc","model with Gaussian constraint(s)",RooArgSet(model,nxib_constraint,mean_constraint));

	cout<<"stage2"<<endl;
	//*************DO FIT**********************************************************
	RooFitResult *res;

	if(constraintflag==1) {
		res = model.fitTo(ds,ExternalConstraints(RooArgSet(nxib_constraint,meanpdg_constraint)),Extended(), RooFit::Save(true), Strategy(2));
	}
	else {
		//res = modelc.fitTo(ds,Constrain(RooArgSet(nxib_constraint,mean_constraint)),Extended(), RooFit::Save(true), Strategy(2));
		res = model.fitTo(ds,ExternalConstraints(RooArgSet(nxib_constraint)),Extended(), RooFit::Save(true), Strategy(2));
	}

	//*****************************************************************************
	cout<<"stage3"<<endl;

	//*********************PLOTTING STUFF*********************************************
	//  RooPlot *frame = Lb_DTF_M_JpsiLConstr.frame();
	RooPlot *frame = new RooPlot(Lb_DTF_M_JpsiLConstr, low,5800,100);
	//  frame->GetYaxis()->SetRangeUser(0,100);
	ds.plotOn(frame,Name("data"));
	model.plotOn(frame,Name("fit"));
	model.plotOn(frame,Components(sig),Name("sig"),LineColor(kMagenta+2));
	model.plotOn(frame,Components(sig1),LineStyle(kDotted),LineColor(kMagenta));
	model.plotOn(frame,Components(sig2),LineStyle(kDotted),LineColor(kMagenta));
	model.plotOn(frame,Components(bkg),LineColor(kRed),Name("bkg"));
	//  if(xibflag!=0)
	model.plotOn(frame,Components(xibshape_smooth),LineColor(kGreen),Name("xib"));
	// if(lstflag!=0)
	model.plotOn(frame,Components(lst1405shape),LineColor(kGreen+2),LineStyle(kDashed),Name("lst1405"));
	model.plotOn(frame,Components(lst1520shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1520"));
	model.plotOn(frame,Components(lst1810shape),LineColor(kBlue+2),LineStyle(kDashed),Name("lst1810"));
	// if(sigmaflag!=0)
	model.plotOn(frame,Components(sigmashape_smooth),LineColor(kBlack),Name("sigma"));
	model.plotOn(frame,Components(misclstshape),LineColor(kRed+2),LineStyle(kDashed),Name("misclst"));

	Double_t chiSquare1 = frame->chiSquare("fit","data");
	cout<<"chi square1/dof = "<<chiSquare1<<endl;
	RooArgSet *floatpar = model.getParameters(ds);
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
	RooPlot *framex2 = new RooPlot(Lb_DTF_M_JpsiLConstr, low,5800,100);
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
	//Check pull mean and RMS. Ideally the pull mean should be 0 and the RMS should be 1
	cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
	cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;
	nsig->Print();
	nbkg.Print();
	//  nxib.Print();
	// Lb_DTF_M_JpsiLConstr.setRange("signal_window",5520.0,5720.0);
	// Lb_DTF_M_JpsiLConstr.setRange("signal_window1",5570.0,5670.0);

	// RooAbsReal* myint = sig.createIntegral(Lb_DTF_M_JpsiLConstr,NormSet(Lb_DTF_M_JpsiLConstr),Range("signal_window1"));

	// myint->Print();

	// for(UInt_t i=0; i<res->numStatusHistory(); i++) {
	//   if(res->statusCodeHistory(i)%10 != 0) {
	//     Error("JpsiLambda",Form("Fit Failed:%s",res->statusLabelHistory(i)));
	//     exit(EXIT_FAILURE);
	//   }
	// }

	auto legend = new TLegend(0.5,0.7,0.9,0.9);
	legend->AddEntry("data","Data","lp");
	legend->AddEntry("fit","Total Fit","l");
	legend->AddEntry("sig","J/#psi #Lambda shape","l");
	legend->AddEntry("sigma","J/#psi #Sigma shape","l");
	legend->AddEntry("xib","J/#psi #Xi shape","l");
	legend->AddEntry("lst1405","J/#psi #Lambda(1405) shape","l");
	legend->AddEntry("lst1520","J/#psi #Lambda(1520) shape","l");
	legend->AddEntry("lst1810","J/#psi #Lambda(1810) shape","l");
	legend->AddEntry("misclst","misc. J/#psi #Lambda* shapes","l");
	legend->AddEntry("bkg","Comb. Bkg. shape","l");
	legend->Draw("same");
}
