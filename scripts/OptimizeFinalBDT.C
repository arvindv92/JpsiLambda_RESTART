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
                                        const char *part)

{
	gROOT->ProcessLine(".x lhcbStyle.C");
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	if(logFlag && isoFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalBDT%d_iso%d_%s_%s_%s_noPID.txt",
		                             run,bdtConf,isoConf,isoVersion, FOM, part),"w");
	}
	if(logFlag && !isoFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalBDT%d_noIso_%s_%s_noPID.txt",
		                             run,bdtConf,FOM, part),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"==> Starting OptimizeFinalBDT: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	const char *type = "LL";

	TFile *fileIn    = nullptr, *fileIn_zeroTracks    = nullptr;
	TFile *fileIn_mc = nullptr, *fileIn_zeroTracks_mc = nullptr;

	TTree *treeIn    = nullptr, *treeIn_zeroTracks    = nullptr, *myTree = nullptr;
	TTree *treeIn_mc = nullptr, *treeIn_zeroTracks_mc = nullptr;

	TString friendFileName            = "";
	TString friendFileName_zeroTracks = "", t = "";

	TH1D *hsig             = nullptr;
	Int_t nEntries         = 0, N = 0;
	Float_t width          = 0.0;
	const char *rootFolder = "";
	std::vector <Double_t> cuts;

	Float_t bdtArray_nonZero[99];
	Float_t fomArray_nonZero[99];

	Float_t bdtArray_Zero[99];
	Float_t fomArray_Zero[99];

	rootFolder = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);

	if(isoFlag)//Using isolation
	{
		fileIn = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks_noPID.root",
		                          rootFolder,type),"READ");
		if (!fileIn)
		{
			cout << "ERROR: could not open data file" << endl;
			exit(1);
		}
		treeIn = (TTree*)fileIn->Get("MyTuple");

		fileIn_zeroTracks = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_ZeroTracks_noPID.root",
		                                     rootFolder,type),"READ");
		if (!fileIn_zeroTracks)
		{
			cout << "ERROR: could not open zeroTracks data file" << endl;
			exit(1);
		}
		treeIn_zeroTracks = (TTree*)fileIn_zeroTracks->Get("MyTuple");

		friendFileName = Form("%s/jpsilambda_%s_FinalBDT%d_iso%d_%s_noPID.root",
		                      rootFolder,type,bdtConf,isoConf,isoVersion);
		friendFileName_zeroTracks = Form("%s/jpsilambda_zeroTracks%s_FinalBDT%d_noPID.root",
		                                 rootFolder,type,bdtConf);

		treeIn->AddFriend("MyTuple",friendFileName);
		treeIn_zeroTracks->AddFriend("MyTuple",friendFileName_zeroTracks);

		if(TString(part) == "sigma")
		{
			width = 260.0;
			fileIn_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_cutoutks_%s_nonZeroTracks_noPID.root",
			                             run,type),"READ");
			treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

			treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_%s_FinalBDT%d_iso%d_%s_noPID.root",
			                                    run,type,bdtConf,isoConf,isoVersion));

			fileIn_zeroTracks_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_cutoutks_%s_ZeroTracks_noPID.root",
			                                        run,type),"READ");
			treeIn_zeroTracks_mc = (TTree*)fileIn_zeroTracks_mc->Get("MyTuple");
			treeIn_zeroTracks_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/JpsiSigma/run%d/jpsisigma_zeroTracks%s_FinalBDT%d_noPID.root",
			                                               run,type,bdtConf));
		}
		else if(TString(part) == "lambda")
		{
			width = 28.0;
			fileIn_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_cutoutks_%s_nonZeroTracks_noPID.root",
			                             run,type),"READ");
			treeIn_mc = (TTree*)fileIn_mc->Get("MyTuple");

			treeIn_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_%s_FinalBDT%d_iso%d_%s_noPID.root",
			                                    run,type,bdtConf,isoConf,isoVersion));

			fileIn_zeroTracks_mc = TFile::Open(Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_cutoutks_%s_ZeroTracks_noPID.root",
			                                        run,type),"READ");
			treeIn_zeroTracks_mc = (TTree*)fileIn_zeroTracks_mc->Get("MyTuple");
			treeIn_zeroTracks_mc->AddFriend("MyTuple",Form("rootFiles/mcFiles/JpsiLambda/JpsiLambda/run%d/jpsilambda_zeroTracks%s_FinalBDT%d_noPID.root",
			                                               run,type,bdtConf));
		}
	}
	else //No Isolation
	{
		fileIn = TFile::Open(Form("%s/jpsilambda_%s_withsw_noPID.root",
		                          rootFolder,type));
		if (!fileIn)
		{
			cout << "ERROR: could not open data file" << endl;
			exit(1);
		}
		treeIn = (TTree*)fileIn->Get("MyTuple");

		friendFileName = Form("%s/jpsilambda_%ssig_FinalBDT%d_noIso_noPID.root",
		                      rootFolder,type,bdtConf);
		treeIn->AddFriend("MyTuple",friendFileName);
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
		Double_t BDT           = 0., BDT_max           = 0.;
		Double_t myFOM         = 0., myFOM_max         = 0.;
		Double_t eff_sig_wt    = 0., eff_sig_wt_max    = 0.;
		Double_t eff_sig_wt_TM = 0., eff_sig_wt_TM_max = 0.;
		Double_t eff_sig       = 0., eff_sig_max       = 0.;
		Double_t eff_sig_TM    = 0., eff_sig_TM_max    = 0.;
		Double_t eff_bkg       = 0., eff_bkg_max       = 0.;
		Double_t sig           = 0., bkg               = 0., siginit = 0.;
		Double_t gbwt = 0.;
		Float_t tauwt = 0.;
		Double_t mcBDT = 0., sumwt_num = 0., sumwt_den = 0.;
		Int_t bkgcat = 0;
		Int_t bkginit = 0;

		if(ctr == 0)
		{
			cout<<"****************Optimizing for nonZeroTracks****************"<<endl;
			mcTree = treeIn_mc;
			if(run == 1)
				siginit = 7385; //6569;
			if(run == 2)
				siginit = 25081; //24741;
		}
		else if(ctr == 1)
		{
			cout<<"****************Optimizing for ZeroTracks****************"<<endl;
			mcTree = treeIn_zeroTracks_mc;
			if(run == 1)
				siginit = 439; //407;
			if(run == 2)
				siginit = 1169; //1195;
		}
		mcTree->SetBranchAddress("gb_wts",&gbwt);
		mcTree->SetBranchAddress("wt_tau",&tauwt);
		mcTree->SetBranchAddress("Lb_BKGCAT",&bkgcat);
		mcTree->SetBranchAddress(Form("BDT%d",bdtConf),&mcBDT);

		Int_t nEntries_mc = mcTree->GetEntries();

		nEntries = myTree->GetEntries();
		cout<<"nEntries= "<<nEntries<<endl;

		myTree->Draw(Form("BDT%d>>hsig0(100,0.0,1.0)",bdtConf),"","goff");

		hsig = (TH1D*)gDirectory->Get("hsig0");  //play around with this

		bkginit = myTree->GetEntries("Lb_DTF_M_JpsiLConstr > 5700")*width/300;//scaling background to width of signal region

		cout<<"siginit = "<<siginit<<endl;
		cout<<"bkginit = "<<bkginit<<endl;

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

			if(ctr == 0) bdtArray_nonZero[i-1] = BDT;
			else if(ctr == 1) bdtArray_Zero[i-1] = BDT;

			eff_sig_TM    = mcTree->GetEntries(Form("(Lb_BKGCAT == 0||Lb_BKGCAT == 50) && BDT%d > %f",bdtConf,BDT))*1.0/denom;
			eff_sig       = mcTree->GetEntries(Form("BDT%d > %f",bdtConf,BDT))*1.0/denom;
			eff_sig_wt_TM = sumwt_num/sumwt_den;

			sig = siginit*eff_sig_wt_TM;
			bkg = myTree->GetEntries(Form("Lb_DTF_M_JpsiLConstr > 5700 && BDT%d > %f",bdtConf,BDT))*width/300;//assumes flat background

			cout<<"SIG = "<<sig<<" BKG = "<<bkg;
			eff_bkg = (Double_t) bkg/bkginit;

			if(TString(part) == "sigma" && TString(FOM) == "Sig")
			{
				N = myTree->GetEntries(Form("Lb_DTF_M_JpsiLConstr > 5360 && Lb_DTF_M_JpsiLConstr < 5600 && BDT%d > %f",bdtConf,BDT));
			}
			else if(TString(part) == "lambda" && TString(FOM) == "Sig")
			{
				N = myTree->GetEntries(Form("Lb_DTF_M_JpsiLConstr > 5606 && Lb_DTF_M_JpsiLConstr < 5634&& BDT%d > %f",bdtConf,BDT));
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
				myFOM = (Double_t)sig/(sqrt(bkg) + 1.5);//a  = 3 chosen
			}

			if(ctr == 0) fomArray_nonZero[i-1] = myFOM;
			else if(ctr == 1) fomArray_Zero[i-1] = myFOM;

			if(myFOM > myFOM_max)
			{
				myFOM_max         = myFOM;
				BDT_max           = BDT;
				eff_sig_TM_max    = eff_sig_TM;
				eff_sig_wt_TM_max = eff_sig_wt_TM;
				eff_bkg_max       = eff_bkg;
			}
			cout<<"For BDT = "<<BDT<<" FOM = "<<myFOM<<" sig_eff = "
			    <<eff_sig_TM*100<<"% sig_eff_wt = "<<eff_sig_wt_TM*100<<"% bkg_eff = "<<eff_bkg*100<<"%"<<endl;
		}
		cout<<"MAXIMUM FOM = "<<myFOM_max<<" at BDT = "<<BDT_max<<" with sig_eff = "
		    <<eff_sig_TM_max*100<<"% sig_eff_wt = "<<eff_sig_wt_TM_max*100<<"% and bkg_eff = "<<eff_bkg_max*100<<"%"<<endl;
		cuts.push_back(BDT_max);
		ctr++;
	}
	TLatex *myLatex = new TLatex();
	myLatex->SetTextFont(42);
	myLatex->SetTextColor(1);
	myLatex->SetTextAlign(12);
	myLatex->SetNDC(kTRUE);

	myLatex->SetTextSize(0.065);

	TCanvas *c1 = new TCanvas();
	TGraph *gr_nonZero = new TGraph(99,bdtArray_nonZero,fomArray_nonZero);
	gr_nonZero->SetTitle("");
	gr_nonZero->SetLineWidth(2);
	gr_nonZero->SetMarkerSize(0.7);
	gr_nonZero->SetMarkerStyle(20);
	gr_nonZero->GetXaxis()->SetTitle("nonZeroTracks BDT cut");
	gr_nonZero->GetYaxis()->SetTitle("S/(#sqrt{B} + 1.5)");
	// gr_nonZero->SetName("");
	gr_nonZero->Draw("AC*");
	myLatex->DrawLatex(0.18,0.85,Form("LHCb Run %d",run));

	TCanvas *c2 = new TCanvas();
	TGraph *gr_Zero = new TGraph(99,bdtArray_Zero,fomArray_Zero);
	gr_Zero->SetTitle("");
	gr_Zero->SetLineWidth(2);
	gr_Zero->SetMarkerSize(0.7);
	gr_Zero->SetMarkerStyle(20);
	gr_Zero->GetXaxis()->SetTitle("ZeroTracks BDT cut");
	gr_Zero->GetYaxis()->SetTitle("S/(#sqrt{B} + 1.5)");
	gr_Zero->Draw("AC*");
	myLatex->DrawLatex(0.18,0.85,Form("LHCb Run %d",run));

	TFile *fileOut = new TFile(Form("/data1/avenkate/JpsiLambda_RESTART/"
	                                "rootFiles/dataFiles/JpsiLambda/run%d/FOM"
	                                "_bdtConf%d_iso%d_%s_noPID.root",
	                                run,bdtConf,isoConf,isoVersion),"RECREATE");
	gr_nonZero->Write();
	gr_Zero->Write();
	fileOut->Close();

	c1->SaveAs(Form("plots/ANA/FOM_run%d_nonZeroTracks_bdtConf%d_iso%d_%s_noPID.pdf",run,bdtConf,isoConf,isoVersion));
	c2->SaveAs(Form("plots/ANA/FOM_run%d_ZeroTracks_bdtConf%d_noPID.pdf",run,bdtConf));

	if(logFlag) gSystem->RedirectOutput(0);
	return cuts;
}
