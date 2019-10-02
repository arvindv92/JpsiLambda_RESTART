#include "OptimizeFinalBDT.h"

void addGraphics(TH1F *h, TString Xtitle, TString Ytitle, int iCol){
	h->SetXTitle(Xtitle);
	//h->SetFillColor(30);
	//int bw = h->GetBinWidth(1);
	h->SetYTitle(Ytitle);
	h->SetStats(kFALSE);
	h->SetMinimum(0.1);
	h->SetMaximum(1.3*h->GetMaximum());
	// h->SetTitleSize(0.1);
	h->SetLineColor(iCol);
	// h->SetMarkerColor(iCol);
	// h->SetMarkerSize(0.7);
	// h->SetMarkerStyle(20);
	// h->GetXaxis()->SetTitleOffset(1.0);
	// h->GetYaxis()->SetTitleOffset(1.15);
	// h->GetXaxis()->SetTitleSize(0.055);
	// h->GetYaxis()->SetTitleSize(0.055);
	// h->GetXaxis()->SetLabelSize(0.045);
	// h->GetYaxis()->SetLabelSize(0.045);
	// h->SetNdivisions(404,"X");
	// h->SetNdivisions(505,"Y");
	// h->SetLineWidth(2);
	h->SetTitle("");

}

std::vector <Double_t> OptimizeFinalBDT(Int_t run, const char* isoVersion, Int_t isoConf,
                                        Int_t bdtConf, Bool_t isoFlag,
                                        Bool_t logFlag, const char* FOM,
                                        const char *part, Bool_t simFlag)

{
	gROOT->ProcessLine(".x lhcbStyle.C");
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	if(!simFlag)//optimize BDT trained on data
	{
		if(logFlag && isoFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalBDT%d_iso%d_%s_%s_%s.txt",
			                             run,bdtConf,isoConf,isoVersion, FOM, part),"w");
		}
		if(logFlag && !isoFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalBDT%d_noIso_%s_%s.txt",
			                             run,bdtConf,FOM, part),"w");
		}
	}
	else//optimize BDT trained on MC
	{
		if(logFlag && isoFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalMCBDT%d_iso%d_%s_%s_%s.txt",
			                             run,bdtConf,isoConf,isoVersion, FOM, part),"w");
		}
		if(logFlag && !isoFlag)
		{
			gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalMCBDT%d_noIso_%s_%s.txt",
			                             run,bdtConf,FOM, part),"w");
		}
	}
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"******************************************"<<endl;
	cout<<"==> Starting OptimizeFinalBDT: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	const char *type = "LL";

	TFile *fileIn(0), *fileIn_zeroTracks(0);
	TTree *treeIn(0), *treeIn_zeroTracks(0), *myTree(0);

	TFile *fileIn_mc(0), *fileIn_zeroTracks_mc(0);
	TTree *treeIn_mc(0), *treeIn_zeroTracks_mc(0);

	TString friendFileName = "";
	TString friendFileName_zeroTracks = "", t = "";

	TString rootFolder      = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);
	TString mcFolder        = Form("rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d",run);
	TString mcFolder_Lambda = Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d",run);

	TH1D *hsig(0);
	TH1D *hPass = new TH1D("hPass","pass",1,0,1);
	TH1D *hAll  = new TH1D("hAll","all",1,0,1);

	Int_t nEntries = 0, N = 0;
	Int_t bkgcat   = 0;
	Int_t bkginit  = 0;

	std::vector <Double_t> cuts;

	Float_t width = 0.0; // width of signal region
	Float_t bdtArray_nonZero[99]    = {0.};
	Float_t fomArray_nonZero[99]    = {0.};
	Float_t fomErrArray_nonZero[99] = {0.};

	Float_t bdtArray_Zero[99]    = {0.};
	Float_t fomArray_Zero[99]    = {0.};
	Float_t fomErrArray_Zero[99] = {0.};

	Float_t bdtArray_noIso[99]    = {0.};
	Float_t fomArray_noIso[99]    = {0.};
	Float_t fomErrArray_noIso[99] = {0.};
	Float_t bdtErrArray[99]       = {0.};

	Double_t BDT   = 0., BDT_max   = 0.;
	Double_t myFOM = 0., myFOM_max = 0.;

	Double_t eff_sig_wt    = 0., eff_sig_wt_max    = 0.;
	Double_t eff_sig_wt_TM = 0., eff_sig_wt_TM_max = 0.;
	Double_t eff_sig       = 0., eff_sig_max       = 0.;
	Double_t eff_sig_TM    = 0., eff_sig_TM_max    = 0.;
	Double_t eff_bkg       = 0., eff_bkg_max       = 0.;


	Double_t sig           = 0., bkg               = 0., siginit = 0.;
	Double_t err_sig = 0., err_bkg = 0., err_FOM = 0.;
	Double_t err_siginit = 0., err_bkginit = 0., err_eff = 0.;
	Double_t err_FOM_max = 0.;
	Double_t gbwt = 0.;
	Float_t tauwt = 0.;
	Double_t mcBDT = 0., sumwt_num = 0., sumwt_den = 0.;

	Double_t sigTot = 0., bkgTot = 0., fomTot = 0.;
	Double_t err_sigTot_sq = 0., err_bkgTot_sq = 0.;
	Double_t err_sigTot = 0., err_bkgTot = 0., err_fomTot = 0.;
	Double_t sig_max = 0., err_sig_max = 0., bkg_max = 0., err_bkg_max = 0;
	Double_t tot_siginit = 0., tot_bkginit = 0.;

	TGraphAsymmErrors *gr = new TGraphAsymmErrors(1);

	Float_t siginit_nonZero = 0., err_siginit_nonZero = 0.;
	Float_t siginit_Zero    = 0., err_siginit_Zero    = 0.;

	ifstream inFile_nonZero(Form("logs/data/JpsiLambda/run%d/sigYield_nonZeroTracks.txt",run));
	inFile_nonZero>>siginit_nonZero;
	inFile_nonZero>>err_siginit_nonZero;
	inFile_nonZero.close();

	ifstream inFile_Zero(Form("logs/data/JpsiLambda/run%d/sigYield_ZeroTracks.txt",run));
	inFile_Zero>>siginit_Zero;
	inFile_Zero>>err_siginit_Zero;
	inFile_Zero.close();

	if(isoFlag)//Using isolation
	{
		fileIn = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks.root",
		                          rootFolder.Data(),type),"READ");
		if (!fileIn)
		{
			cout << "ERROR: could not open data file" << endl;
			exit(1);
		}
		treeIn = (TTree*)fileIn->Get("MyTuple");

		fileIn_zeroTracks = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_ZeroTracks.root",
		                                     rootFolder.Data(),type),"READ");
		if (!fileIn_zeroTracks)
		{
			cout << "ERROR: could not open zeroTracks data file" << endl;
			exit(1);
		}
		treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");

		if(!simFlag)
		{
			friendFileName = Form("%s/jpsilambda_%s_FinalBDT%d_iso%d_%s.root",
			                      rootFolder.Data(),type,bdtConf,isoConf,isoVersion);
			friendFileName_zeroTracks = Form("%s/jpsilambda_zeroTracks%s_FinalBDT%d.root",
			                                 rootFolder.Data(),type,bdtConf);
		}
		else
		{
			friendFileName = Form("%s/jpsilambda_%s_MCFinalBDT%d_MCiso%d_%s.root",
			                      rootFolder.Data(),type,bdtConf,isoConf,isoVersion);
			friendFileName_zeroTracks = Form("%s/jpsilambda_zeroTracks%s_MCFinalBDT%d.root",
			                                 rootFolder.Data(),type,bdtConf);
		}

		treeIn->AddFriend("MyTuple",friendFileName);
		treeIn_zeroTracks->AddFriend("MyTuple",friendFileName_zeroTracks);

		if(TString(part) == "Sigma")
		{
			width = 260.0;
			fileIn_mc = TFile::Open(Form("%s/jpsisigma_cutoutks_%s_nonZeroTracks.root",
			                             mcFolder.Data(),type),"READ");
			treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

			fileIn_zeroTracks_mc = TFile::Open(Form("%s/jpsisigma_cutoutks_%s_ZeroTracks.root",
			                                        mcFolder.Data(),type),"READ");
			treeIn_zeroTracks_mc = (TTree*)fileIn_zeroTracks_mc->Get("MyTuple");

			if(!simFlag)
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsisigma_%s_FinalBDT%d_iso%d_%s.root",
				                                    mcFolder.Data(),type,bdtConf,isoConf,isoVersion));

				treeIn_zeroTracks_mc->AddFriend("MyTuple",Form("%s/jpsisigma_zeroTracks%s_FinalBDT%d.root",
				                                               mcFolder.Data(),type,bdtConf));
			}
			else
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsisigma_%s_MCFinalBDT%d_MCiso%d_%s.root",
				                                    mcFolder.Data(),type,bdtConf,isoConf,isoVersion));

				treeIn_zeroTracks_mc->AddFriend("MyTuple",Form("%s/jpsisigma_zeroTracks%s_MCFinalBDT%d.root",
				                                               mcFolder.Data(),type,bdtConf));
			}
		}
		else if(TString(part) == "Lambda")
		{
			width = 28.0;
			fileIn_mc = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks.root",
			                             mcFolder_Lambda.Data(),type),"READ");
			treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

			fileIn_zeroTracks_mc = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_ZeroTracks.root",
			                                        mcFolder_Lambda.Data(),type),"READ");
			treeIn_zeroTracks_mc = (TTree*)fileIn_zeroTracks_mc->Get("MyTuple");

			if(!simFlag)
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsilambda_%s_FinalBDT%d_iso%d_%s.root",
				                                    mcFolder_Lambda.Data(),type,bdtConf,isoConf,isoVersion));

				treeIn_zeroTracks_mc->AddFriend("MyTuple",Form("%s/jpsilambda_zeroTracks%s_FinalBDT%d.root",
				                                               mcFolder_Lambda.Data(),type,bdtConf));
			}
			else
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsilambda_%s_MCFinalBDT%d_MCiso%d_%s.root",
				                                    mcFolder_Lambda.Data(),type,bdtConf,isoConf,isoVersion));

				treeIn_zeroTracks_mc->AddFriend("MyTuple",Form("%s/jpsilambda_zeroTracks%s_MCFinalBDT%d.root",
				                                               mcFolder_Lambda.Data(),type,bdtConf));
			}
		}
	}
	else //No Isolation
	{
		fileIn = TFile::Open(Form("%s/jpsilambda_cutoutks_%s.root",
		                          rootFolder.Data(),type));
		if (!fileIn)
		{
			cout << "ERROR: could not open data file" << endl;
			exit(1);
		}
		treeIn = (TTree*)fileIn->Get("MyTuple");

		if(!simFlag)
		{
			friendFileName = Form("%s/jpsilambda_%s_FinalBDT%d_noIso.root",
			                      rootFolder.Data(),type,bdtConf);
		}
		else
		{
			friendFileName = Form("%s/jpsilambda_%s_MCFinalBDT%d_noIso.root",
			                      rootFolder.Data(),type,bdtConf);
		}
		treeIn->AddFriend("MyTuple",friendFileName);

		if(TString(part) == "Sigma")
		{
			width = 260.0;
			fileIn_mc = TFile::Open(Form("%s/jpsisigma_cutoutks_%s.root",
			                             mcFolder.Data(),type),"READ");
			treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

			if(!simFlag)
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsisigma_%s_FinalBDT%d_noIso.root",
				                                    mcFolder.Data(),type,bdtConf));
			}
			else
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsisigma_%s_MCFinalBDT%d_noIso.root",
				                                    mcFolder.Data(),type,bdtConf));
			}
		}
		else if(TString(part) == "Lambda")
		{
			width = 28.0;
			fileIn_mc = TFile::Open(Form("%s/jpsilambda_cutoutks_%s.root",
			                             mcFolder_Lambda.Data(),type),"READ");
			treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

			if(!simFlag)
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsilambda_%s_FinalBDT%d_noIso.root",
				                                    mcFolder_Lambda.Data(),type,bdtConf));
			}
			else
			{
				treeIn_mc->AddFriend("MyTuple",Form("%s/jpsilambda_%s_MCFinalBDT%d_noIso.root",
				                                    mcFolder_Lambda.Data(),type,bdtConf));
			}
		}
	}

	Int_t ctr = 0;
	TList *list = new TList();
	list->Add(treeIn);

	if(isoFlag)
	{
		list->Add(treeIn_zeroTracks);
	}

	TIterator *iter = (TIterator*)list->MakeIterator();
	TTree *mcTree = nullptr;

	while((myTree = (TTree*)iter->Next()))
	{
		BDT           = 0., BDT_max           = 0.;
		myFOM         = 0., myFOM_max         = 0.;
		eff_sig_wt    = 0., eff_sig_wt_max    = 0.;
		eff_sig_wt_TM = 0., eff_sig_wt_TM_max = 0.;
		eff_sig       = 0., eff_sig_max       = 0.;
		eff_sig_TM    = 0., eff_sig_TM_max    = 0.;
		eff_bkg       = 0., eff_bkg_max       = 0.;
		sig           = 0., bkg               = 0., siginit = 0.;
		gbwt = 0.;
		tauwt = 0.;
		mcBDT = 0., sumwt_num = 0., sumwt_den = 0.;
		bkgcat = 0;
		bkginit = 0;

		if(!isoFlag)
		{
			cout<<"****************Optimizing for noIso****************"<<endl;
			mcTree = treeIn_mc;
			// if(run == 1)
			// {
			siginit     = siginit_nonZero+siginit_Zero;        // 7386; //hardcoded. Dont like this.
			err_siginit = sqrt(pow(err_siginit_nonZero,2)+pow(err_siginit_Zero,2));        //153;
			// }
			// else if(run == 2)
			// {
			//      siginit     = 26050;
			//      err_siginit = 248;
			// }
		}

		if(ctr == 0 && isoFlag)
		{
			cout<<"****************Optimizing for nonZeroTracks****************"<<endl;
			mcTree = treeIn_mc;
			// if(run == 1)
			// {
			siginit     = siginit_nonZero;        //6939;
			err_siginit = err_siginit_nonZero;        //151;
			// }
			// else if(run == 2)
			// {
			//      siginit     = 24898;
			//      err_siginit = 245;
			// }
		}

		else if(ctr == 1 && isoFlag)
		{
			cout<<"****************Optimizing for ZeroTracks****************"<<endl;
			mcTree = treeIn_zeroTracks_mc;
			// if(run == 1)
			// {
			siginit     = siginit_Zero;        //447;
			err_siginit = err_siginit_Zero;        //26;
			// }
			// else if(run == 2)
			// {
			//      siginit     = 1151;
			//      err_siginit = 39;
			// }
		}
		mcTree->SetBranchAddress("GB_WT",&gbwt);
		mcTree->SetBranchAddress("wt_tau",&tauwt);
		mcTree->SetBranchAddress("Lb_BKGCAT",&bkgcat);
		mcTree->SetBranchAddress(Form("BDT%d",bdtConf),&mcBDT);

		Int_t nEntries_mc = mcTree->GetEntries();

		nEntries = myTree->GetEntries();
		cout<<"nEntries= "<<nEntries<<endl;

		myTree->Draw(Form("BDT%d>>hsig(100,0.0,1.0)",bdtConf),"","goff");

		hsig = (TH1D*)gDirectory->Get("hsig");

		bkginit = myTree->GetEntries("!(Lb_DTF_M_JpsiLConstr>5770 && Lb_DTF_M_JpsiLConstr<5810) && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 6000")*width/260;        //scaling background to width of signal region
		err_bkginit = sqrt((double)myTree->GetEntries("!(Lb_DTF_M_JpsiLConstr>5770 && Lb_DTF_M_JpsiLConstr<5810) && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 6000"))*width/260;

		cout<<"siginit = "<<siginit<<" +/- "<<err_siginit<<endl;
		cout<<"bkginit = "<<bkginit<<" +/- "<<err_bkginit<<endl;

		Int_t denom = mcTree->GetEntries("(Lb_BKGCAT==0||Lb_BKGCAT==50)");
		for(Int_t j = 0; j<nEntries_mc; j++)
		{
			mcTree->GetEntry(j);
			if(bkgcat==0||bkgcat==50)
			{
				sumwt_den += gbwt*tauwt;
			}
		}
		for(Int_t i = 1; i < 100; i++)
		{
			sumwt_num = 0.;
			BDT = hsig->GetBinCenter(i);

			for(Int_t j = 0; j<nEntries_mc; j++)
			{
				mcTree->GetEntry(j);
				if(mcBDT < BDT || !(bkgcat==0||bkgcat==50))
				{
					continue;
				}
				else
				{
					sumwt_num += gbwt*tauwt;
				}
			}

			if(ctr == 0 && isoFlag) bdtArray_nonZero[i-1] = BDT;
			else if(ctr == 1 && isoFlag) bdtArray_Zero[i-1] = BDT;
			else if(!isoFlag) bdtArray_noIso[i-1] = BDT;

			eff_sig_TM    = mcTree->GetEntries(Form("(Lb_BKGCAT == 0||Lb_BKGCAT == 50) && BDT%d > %f",bdtConf,BDT))*1.0/denom;
			eff_sig       = mcTree->GetEntries(Form("BDT%d > %f",bdtConf,BDT))*1.0/denom;
			eff_sig_wt_TM = sumwt_num/sumwt_den;

			hPass->SetBinContent(1,sumwt_num);
			hAll->SetBinContent(1,sumwt_den);

			gr->BayesDivide(hPass,hAll);

			err_eff = gr->GetErrorY(0); //errors are actually asymmetric. approximating with sym. error for now.

			sig = siginit*eff_sig_wt_TM;
			if(eff_sig_wt_TM > 0)
				err_sig = sig*sqrt((double) (pow(err_siginit/siginit,2) + pow(err_eff/eff_sig_wt_TM,2)));
			else
				err_sig = 0.;
			Int_t nBkg = myTree->GetEntries(Form("!(Lb_DTF_M_JpsiLConstr>5770 && Lb_DTF_M_JpsiLConstr<5810) && Lb_DTF_M_JpsiLConstr > 5700 && Lb_DTF_M_JpsiLConstr < 6000 && BDT%d > %f",bdtConf,BDT));
			bkg = nBkg*width/260;        //assumes flat background
			err_bkg = sqrt((double)nBkg)*width/260;

			cout<<"SIG = "<<sig<<" BKG = "<<bkg;
			eff_bkg = (Double_t) bkg/bkginit;

			if(TString(part) == "Sigma" && TString(FOM) == "Sig")
			{
				N = myTree->GetEntries(Form("Lb_DTF_M_JpsiLConstr > 5360 && Lb_DTF_M_JpsiLConstr < 5600 && BDT%d > %f",bdtConf,BDT));
			}
			else if(TString(part) == "Lambda" && TString(FOM) == "Sig")
			{
				N = myTree->GetEntries(Form("Lb_DTF_M_JpsiLConstr > 5606 && Lb_DTF_M_JpsiLConstr < 5634 && BDT%d > %f",bdtConf,BDT));
			}
			if(TString(FOM) == "Sig")
			{
				if(N > 0)
				{
					myFOM = (Double_t)sig/sqrt((double)N);
				}
				else
				{
					cout<<"NO"<<endl;
					continue;
				}
			}
			if(TString(FOM) == "Punzi")
			{
				myFOM = (Double_t)sig/(sqrt(bkg) + 1.5);        //a  = 3 chosen

				if(sig > 0 && bkg > 0)
					err_FOM = myFOM * sqrt( (double) ( pow(err_sig/sig,2) + pow( (err_bkg/(2*sqrt(bkg)*(1.5 + sqrt(bkg)))),2 ) )  );
				else if(!(bkg > 0) && sig > 0)
					err_FOM = myFOM * (err_sig/sig);
				else
					err_FOM = 0.;
			}

			if(ctr == 0 && isoFlag)
			{
				fomArray_nonZero[i-1] = myFOM;
				fomErrArray_nonZero[i-1] = err_FOM;
			}
			else if(ctr == 1 && isoFlag)
			{
				fomArray_Zero[i-1] = myFOM;
				fomErrArray_Zero[i-1] = err_FOM;
			}
			else if(!isoFlag)
			{
				fomArray_noIso[i-1] = myFOM;
				fomErrArray_noIso[i-1] = err_FOM;
			}

			if(myFOM > myFOM_max)
			{
				myFOM_max         = myFOM;
				BDT_max           = BDT;
				eff_sig_TM_max    = eff_sig_TM;
				eff_sig_wt_TM_max = eff_sig_wt_TM;
				eff_bkg_max       = eff_bkg;
				err_FOM_max       = err_FOM;

				sig_max = sig;
				bkg_max = bkg;
				err_sig_max = err_sig;
				err_bkg_max = err_bkg;
			}
			cout<<"For BDT = "<<BDT<<" FOM = "<<myFOM<<" +/- "<<err_FOM<<" sig_eff = "
			    <<eff_sig_TM*100<<"% sig_eff_wt = "<<eff_sig_wt_TM*100<<"% bkg_eff = "<<eff_bkg*100<<"%"<<endl;
		}
		cout<<"MAXIMUM FOM = "<<myFOM_max<<" +/- "<<err_FOM_max<<" at BDT = "<<BDT_max<<" with sig_eff = "
		    <<eff_sig_TM_max*100<<"% sig_eff_wt = "<<eff_sig_wt_TM_max*100<<"% and bkg_eff = "<<eff_bkg_max*100<<"%"<<endl;
		cuts.push_back(BDT_max);

		if(isoFlag)
		{
			sigTot += sig_max;//eff_sig_wt_TM_max*siginit;
			bkgTot += bkg_max;//eff_bkg_max*bkginit;

			tot_siginit += siginit;
			tot_bkginit += bkginit;

			err_sigTot_sq += err_sig_max*err_sig_max;
			err_bkgTot_sq += err_bkg_max*err_bkg_max;
		}
		ctr++;
	}
	TLatex *myLatex = new TLatex();
	myLatex->SetTextFont(42);
	myLatex->SetTextColor(1);
	myLatex->SetTextAlign(12);
	myLatex->SetNDC(kTRUE);

	myLatex->SetTextSize(0.065);

	if(isoFlag)
	{
		fomTot = sigTot/(sqrt((double)bkgTot)+1.5);
		err_sigTot = sqrt(err_sigTot_sq);
		err_bkgTot = sqrt(err_bkgTot_sq);

		if(sigTot > 0 && bkgTot > 0)
			err_fomTot = fomTot * sqrt( (double) ( pow(err_sigTot/sigTot,2) + pow( (err_bkgTot/(2*sqrt(bkgTot)*(1.5 + sqrt(bkgTot)))),2 ) )  );
		else if(!(bkgTot > 0) && sigTot > 0)
			err_fomTot = fomTot * (err_sigTot/sigTot);
		else
			err_fomTot = 0.;

		cout<<"$$$$FOM TOTAL = "<<fomTot<<" +/- "<<err_fomTot<<endl;
		cout<<"$$$$SIG EFF TOTAL = "<<(sigTot/tot_siginit)*100<<endl;
		cout<<"$$$$BKG EFF TOTAL = "<<(bkgTot/tot_bkginit)*100<<endl;
		TCanvas *c1 = new TCanvas();
		TGraphErrors *gr_nonZero = new TGraphErrors(99,bdtArray_nonZero,fomArray_nonZero,bdtErrArray,fomErrArray_nonZero);
		gr_nonZero->SetTitle("");
		gr_nonZero->SetLineWidth(2);
		gr_nonZero->SetMarkerSize(0.7);
		gr_nonZero->SetMarkerStyle(20);
		gr_nonZero->GetXaxis()->SetTitle("nonZeroTracks BDT cut");
		gr_nonZero->GetYaxis()->SetTitle("S/(#sqrt{B} + 1.5)");
		// gr_nonZero->SetName("");
		gr_nonZero->Draw("ALP");
		myLatex->DrawLatex(0.18,0.85,Form("LHCb Run %d",run));

		TCanvas *c2 = new TCanvas();
		TGraphErrors *gr_Zero = new TGraphErrors(99,bdtArray_Zero,fomArray_Zero,bdtErrArray,fomErrArray_Zero);
		gr_Zero->SetTitle("");
		gr_Zero->SetLineWidth(2);
		gr_Zero->SetMarkerSize(0.7);
		gr_Zero->SetMarkerStyle(20);
		gr_Zero->GetXaxis()->SetTitle("ZeroTracks BDT cut");
		gr_Zero->GetYaxis()->SetTitle("S/(#sqrt{B} + 1.5)");
		gr_Zero->Draw("ALP");
		myLatex->DrawLatex(0.18,0.85,Form("LHCb Run %d",run));

		TFile *fileOut = new TFile(Form("/data1/avenkate/JpsiLambda_RESTART/"
		                                "rootFiles/dataFiles/JpsiLambda/run%d/FOM"
		                                "_bdtConf%d_iso%d_%s.root",
		                                run,bdtConf,isoConf,isoVersion),"RECREATE");
		gr_nonZero->Write();
		gr_Zero->Write();
		fileOut->Close();

		c1->SaveAs(Form("plots/ANA/FOM_run%d_nonZeroTracks_bdtConf%d_iso%d_%s.pdf",run,bdtConf,isoConf,isoVersion));
		c2->SaveAs(Form("plots/ANA/FOM_run%d_ZeroTracks_bdtConf%d.pdf",run,bdtConf));
	}
	else
	{
		TCanvas *c1 = new TCanvas();
		TGraphErrors *gr_noIso = new TGraphErrors(99,bdtArray_noIso,fomArray_noIso,bdtErrArray,fomArray_noIso);
		gr_noIso->SetTitle("");
		gr_noIso->SetLineWidth(2);
		gr_noIso->SetMarkerSize(0.7);
		gr_noIso->SetMarkerStyle(20);
		gr_noIso->GetXaxis()->SetTitle("noIso BDT cut");
		gr_noIso->GetYaxis()->SetTitle("S/(#sqrt{B} + 1.5)");
		// gr_noIso->SetName("");
		gr_noIso->Draw("ALP");
		myLatex->DrawLatex(0.18,0.85,Form("LHCb Run %d",run));

		TFile *fileOut = new TFile(Form("/data1/avenkate/JpsiLambda_RESTART/"
		                                "rootFiles/dataFiles/JpsiLambda/run%d/FOM"
		                                "_bdtConf%d_noIso.root",
		                                run,bdtConf),"RECREATE");
		gr_noIso->Write();
		fileOut->Close();
		c1->SaveAs(Form("plots/ANA/FOM_run%d_nonZeroTracks_bdtConf%d_noIso.pdf",run,bdtConf));
	}

	if(isoFlag)
	{
		ofstream fileOut(Form("logs/data/JpsiLambda/run%d/bestBDTCuts_BDT%d_iso%d_%s.txt",run,bdtConf,isoConf,isoVersion));
		fileOut<<cuts[0]<<endl;
		fileOut<<cuts[1]<<endl;
		fileOut.close();
	}
	else
	{
		ofstream fileOut(Form("logs/data/JpsiLambda/run%d/bestBDTCuts_noIso_BDT%d.txt",run,bdtConf));
		fileOut<<cuts[0]<<endl;
		fileOut.close();
	}

	if(logFlag) gSystem->RedirectOutput(0);
	return cuts;
}
