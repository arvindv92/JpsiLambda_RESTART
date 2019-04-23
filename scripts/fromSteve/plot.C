#include "CodeInput.h"
#include <iostream>
#include <fstream>
using std::ifstream;

using namespace std;
using namespace RooFit;

float mass,err_mass;

// Systematic shift of relative tracking efficiency
//  0 = no shift  (nominal)
// 1 = +1 sigtma
// 1 = -1 sigma
float shift = 0;

// Do efficiency calculation and compute f_Xib / f_Lb
bool doEfficiency_And_Rate = false;

int do_splot = 0;
// do_splot = 0 ==> No sPlot
// do_splot = 1 ==> Do Lambda_b --> J/psi Lambda0 (DD)
// do_splot = 2 ==> Do Lambda_b --> J/psi Lambda0 (LL)
// do_splot = 11 ==> Do Xi_b --> J/psi Xi-


double geteff(TH1F* hn, TH1F* hd)
{

	cout << "hn bin content 10: " << hn->GetBinContent(10) << endl;
	cout << "hd bin content 10: " << hd->GetBinContent(10) << endl;

	TH1F *he = (TH1F*)hn->Clone("he");
	he->Divide(hd);
	double num = 0;
	double den = 0;
	double e, c;
	for(int i=0; i<hn->GetNbinsX(); i++)
	{
		e = he->GetBinContent(i+1);
		c = hd->GetBinContent(i+1);
		cout << i << " " << e << " " << c << endl;
		if(c>0)
		{
			num = num + e*(1.0-e)*c;
			den = den + c;
		}
	}
	if(num <= 0 || den <= 0)
	{
		cout << "Shoot, numeratory or denominator <= 0, ??? " << num << " " << den << endl;
		return 999.0;
	}
	num = sqrt(num);
	delete he;
	return num/den;
}


void makeSomePlots()
{

	TCanvas *cplot = new TCanvas("cplot","Plots",1200,800);
	cplot->Divide(2,2);
	TH1F *hL_dataDD = new TH1F("hL_dataDD","Mass",80,1115.7-20,1115.7+20);
	TH1F *hL_simDD = new TH1F("hL_simDD","Mass",80,1115.7-20,1115.7+20);
	TH1F *hL_dataLL = new TH1F("hL_dataLL","Mass",80,1115.7-20,1115.7+20);
	TH1F *hL_simLL = new TH1F("hL_simLL","Mass",80,1115.7-20,1115.7+20);
	TH1F *hXim_data = new TH1F("hXim_data","Mass",35,1305,1340);
	TH1F *hXim_sim = new TH1F("hXim_sim","Mass",35,1305,1340);
	cplot->cd(1);
	chainl->Draw("L_M>>hL_simDD",MCLambdabDD_Cut);
	hL_simDD->SetMinimum(0);
	hL_simDD->Draw();
	cplot->cd(2);
	treeDD->Draw("L_M>>hL_dataDD",LambdabDD_Cut);
	hL_dataDD->SetMinimum(0);
	hL_dataDD->Draw();

	cplot->cd(3);
	chainl->Draw("L_M>>hL_simLL",MCLambdabLL_Cut);
	hL_simLL->SetMinimum(0);
	hL_simLL->Draw();
	cplot->cd(4);
	treeLL->Draw("L_M>>hL_dataLL",LambdabLL_Cut);
	hL_dataLL->SetMinimum(0);
	hL_dataLL->Draw();

	/*
	   cplot->cd(3);
	   chainx->Draw("Xi_M-L_M+1115.7>>hXim_sim",MCXib_Cut);
	   hXim_sim->SetMinimum(0);
	   cplot->cd(4);
	   treeXi->Draw("Xi_M-L_M+1115.7>>hXim_data",Xib_Cut);
	   hXim_data->SetMinimum(0);
	   hXim_data->Draw();
	 */
}

void plot()
{
	year = "Run1";
	makePDF = false;
	altFitFunction = false;

	bool makeLMassPlots = false; // def = false
	if(makeLMassPlots)
	{
		LbMassWin = 30;
		mcFit = true;
	}

	mcFit = false;
	mcOpt = 5;
	setConditions();
	//cout << applyPrPiKinWeight << endl;
	//return;
	modes[0] = "#Lambda_{b}^{0}#rightarrow#J/#psi#Lambda";
	modes[1] = "#Xi_{b}^{0}#rightarrow#J/#psi#Xi^{#font[122]{-}}";

	// Apply systematic changes, as needed.
	if(doSys) applySys();

	if(mcOpt==1) preScale = 1;

	nbinLb = (mHiLb - mLowLb) / bwLb;
	nbinXib = (mHiXib - mLowXib) / bwXib;

	// Get Histogram weight files
	getWeightHistos();

	// Get data files
	getDataFiles();
	getLbMCFilesForPlot();
	//chainxt->Draw("Xi_bminus_TRUETAU");
	//chainlt->Draw("Lambda_b0_TRUETAU");
	//return;

	// Get cuts
	getCuts();

	if(makeLMassPlots)
	{
		makeSomePlots();
		return;
	}


	// Fill Lb mass histograms
	fillLbMass();

	// Get mass fit functions
	getLbMassFitFunctions();

	// Fit and Plot result
	doLbFitAndPlot();

	//return;

	// Fill Lb mass histograms
	fillXibMass();

	// Get Xib mass fit functions
	getXibMassFitFunctions();

	// Fit and Plot result
	doXibFitAndPlot();



	//return;

	if(doEfficiency_And_Rate)
	{

		//getLbMCFilesForPlot();

		TCanvas *c2 = new TCanvas("c2","",600,600);

		TH1F *h1mcr = new TH1F("h1mcr","#Lambda_{b}^{0} pT",100,0,50);
		TH1F *h1mcg = new TH1F("h1mcg","#Lambda_{b}^{0} pT",100,0,50);
		h1mcr->Sumw2();
		h1mcg->Sumw2();

		double LbNum, LbDen;
		double uLbNum, uLbDen;
		double XibNumL, XibNumD, XibNum, XibDen;
		double uXibNumL, uXibNumD, uXibNum, uXibDen;
		double w, wl, wLb, wL, wPi, wJpsi, wt, wEta;
		double jpsiP, LbP, XibP;
		double wptXib, wetaXib;
		bool useTruth;

		if(mcOpt>=2)
		{
			int ib;
			double etaLb;
			double eta, etaPi;
			//double Lb_P, Lb_PT, Lb_TAU, Lb_TRUETAU;
			//int Lb_BKGCAT;

			TFile *fout = new TFile("/data1/sblusk/junk.root","RECREATE");
			newtree = chainl->CopyTree(MCLambdabDD_Cut);
			cout << "Lb newtree = " << newtree << ", nentry = " << newtree->GetEntriesFast() << endl;
			setLbBranchAddresses();


			Long64_t nentries = newtree->GetEntriesFast();
			TProfile *ha = new TProfile("ha","weight vs p",50,0,300,0,5);
			cout << "Number of reco'd events in MC = " << nentries << endl;
			Long64_t ientry;
			//nentries = 1000;
			// Loop over reconstructed MC tuple
			for (Long64_t i=0; i<nentries; i++)
			{
				ientry = newtree->LoadTree(i);
				if (ientry < 0) break;
				newtree->GetEntry(i);
				if(i%10000 == 0) cout << "At entry = " << i << endl;
				wLb = 1.0; wL = 1.0; wPi = 1.0; wJpsi=1.0, wl = 1.0;


				useTruth = true;
				if(Lb_TRUETAU==0 || Jpsi_TRUEPT==0 || L_TRUEPT==0 || pi_TRUEPT==0 || p_TRUEPT==0) useTruth = false;

				// Lifetime weight
				if(year == "Run1")
				{
					if(Lb_TRUETAU>0)
					{
						wl = exp(-1000*Lb_TRUETAU/1.470)/exp(-1000*Lb_TRUETAU/1.425);
					}
					else
					{
						wl = exp(-1000*Lb_TAU/1.470)/exp(-1000*Lb_TAU/1.425);
					}
				}

				// Pion weight
				if(applyPrPiKinWeight)
				{
					if(mcOpt==4)
					{
						if(applyPiKinWeight && !useTruth) wPi = get2DKinWt(whistPi,pi_PT,pi_PZ,0);
						if(applyPiKinWeight && useTruth) wPi = get2DKinWt(whistPi,pi_TRUEPT,pi_TRUEP_Z,0);
					}
					else if(mcOpt==5)
					{
						if(!useTruth)
						{
							pPr = sqrt(p_PT**2+p_PZ**2);
							pPi = sqrt(pi_PT**2+pi_PZ**2);
						}
						else
						{
							pPr = sqrt(p_TRUEPT**2+p_TRUEP_Z**2);
							pPi = sqrt(pi_TRUEPT**2+pi_TRUEP_Z**2);
						}
						asym = (pPr-pPi) / (pPr+pPi);
						wPi = get1DKinWt(whistPrPi,asym);
					}
				}

				if(!useTruth)
				{
					jpsiP = Jpsi_P;
					LbP = Lb_P;
				}
				else
				{
					jpsiP = sqrt(Jpsi_TRUEPT**2+Jpsi_TRUEP_Z**2);
					LbP = sqrt(Lb_TRUEPT**2+Lb_TRUEP_Z**2);
				}

				if(applyLbKinWeight && !useTruth) wLb = get2DKinWt(whistLb,Lb_PT,Lb_P,1);
				if(applyLbKinWeight && useTruth) wLb = get2DKinWt(whistLb,Lb_TRUEPT,Lb_TRUEP_Z,0);
				if(applyLKinWeight && !useTruth) wL = get2DKinWt(whistL,L_PT,L_P,1);
				if(applyLKinWeight && useTruth) wL = get2DKinWt(whistL,L_TRUEPT,L_TRUEP_Z,0);
				if(applyJpsiLWeight) wJpsi = (JpsiWeightPar[0]+JpsiWeightPar[1]*(jpsiP/LbP)+JpsiWeightPar[2]*(jpsiP/LbP)**2);


				//if(eventNumber==247539 ||  eventNumber==247623 || eventNumber==206593){
				if(i<200) cout << eventNumber << " " << useTruth << " " << eventNumber << " " << wl << " " << wLb << " " << wL << " " << wPi << " " << wJpsi << endl;
				//  cout << L_TRUEPT << " " << L_TRUEP_Z << " " << pi_TRUEPT << " " << pi_TRUEP_Z << " " << Jpsi_TRUEPT << " " << Jpsi_TRUEP_Z << endl;
				// }

				//cout << useTruth << " " << eventNumber << " " << wl << " " << wLb << " " << wL << " " << wPi << " " << wJpsi << endl;

				h1mcr->Fill(0.001*Lb_PT,wl*wLb*wL*wPi*wJpsi);
			}

			// Loop over generated MC tuple
			nentries = chainlt->GetEntriesFast();
			cout << "Number of gen'd events in MC = " << nentries << endl;
			//nentries = 1000;
			for (Long64_t i=0; i<nentries; i++)
			{
				if(preScale > 1 && i%preScale != 0) continue;
				ientry = chainlt->LoadTree(i);
				if (ientry < 0) break;
				chainlt->GetEntry(i);

				if(i%10000 == 0) cout << "At entry = " << i << endl;

				if(chargeCut>0 && J_psi_1S_MC_MOTHER_ID<0) continue;
				if(chargeCut<0 && J_psi_1S_MC_MOTHER_ID>0) continue;

				wLb = 1.0; wL = 1.0; wPi = 1.0; wJpsi=1.0, wl = 1.0;

				// Lifetime weight
				if(year == "Run1")
				{
					wl = exp(-1000*Lambda_b0_TRUETAU/1.470)/exp(-1000*Lambda_b0_TRUETAU/1.425);
				}

				/*
				   if(i>10000) cout << J_psi_1S_TRUEP_Z << " " << J_psi_1S_TRUEPT << " " << Lambda_b0_TRUEP_Z << " " << Lambda_b0_TRUEPT
				                 << " " << Lambda0_TRUEP_Z << " " << Lambda0_TRUEPT << " " << piminus_TRUEPT << " "
				                 << piminus_TRUEP_Z << " " << pplus_TRUEPT << " " << pplus_TRUEP_Z << endl;
				 */

				// Pion weight
				if(applyPrPiKinWeight)
				{
					if(mcOpt==4)
					{
						if(applyPiKinWeight)
						{
							wPi = get2DKinWt(whistPi,piminus_TRUEPT,piminus_TRUEP_Z,0);
						}
					}
					else if(mcOpt==5)
					{
						pPr = sqrt(pplus_TRUEPT**2+pplus_TRUEP_Z**2);
						pPi = sqrt(piminus_TRUEPT**2+piminus_TRUEP_Z**2);
						asym = (pPr-pPi) / (pPr+pPi);
						wPi = get1DKinWt(whistPrPi,asym);
					}
				}

				jpsiP = sqrt(J_psi_1S_TRUEPT**2+J_psi_1S_TRUEP_Z**2);
				LbP = sqrt(Lambda_b0_TRUEPT**2+Lambda_b0_TRUEP_Z**2);

				if(applyLbKinWeight) wLb = get2DKinWt(whistLb,Lambda_b0_TRUEPT,Lambda_b0_TRUEP_Z,0);
				if(applyLKinWeight) wL = get2DKinWt(whistL,Lambda0_TRUEPT,Lambda0_TRUEP_Z,0);
				if(applyJpsiLWeight) wJpsi = (JpsiWeightPar[0]+JpsiWeightPar[1]*(jpsiP/LbP)+JpsiWeightPar[2]*(jpsiP/LbP)**2);

				//if(eventNumber==247539 ||  eventNumber==247623 || eventNumber==206593){
				//    cout << eventNumber << " " << wl << " " << wLb << " " << wL << " " << wPi << " " << wJpsi << endl;
				//  cout << Lambda0_TRUEPT << " " << Lambda0_TRUEP_Z << " " << piminus_TRUEPT << " " << piminus_TRUEP_Z << " " << J_psi_1S_TRUEPT << " " << J_psi_1S_TRUEP_Z << endl;
				//}

				h1mcg->Fill(0.001*Lambda_b0_TRUEPT,wl*wLb*wL*wPi*wJpsi*preScale);
			}
			LbNum = h1mcr->Integral();
			LbDen = h1mcg->Integral();
		}
		else  //if mcOpt < 2
		{
			chainl->Draw("0.001*Lb_PT>>h1mcr",MCLambdabDD_Cut);
			chainlt->Draw("0.001*Lambda_b0_TRUEPT>>h1mcg");
			LbNum = h1mcr->Integral();
			LbDen = h1mcg->Integral();
		}

		//cout << chainlt->GetEntries() << endl;
		double effLambdab = LbNum / LbDen;
		double err_effLambdab = sqrt(effLambdab*(1-effLambdab) / LbDen);

		TCanvas *cc = new TCanvas("cc","",800,400);
		cc->Divide(2,1);
		cc->cd(1);
		h1mcr->Draw();
		cc->cd(2);
		h1mcg->Draw();
		cc->Update();

		//double erreff = geteff(h1mcr,h1mcg);

		cout << "=========================================================="<<endl;
		cout << " Reco Lb = " << h1mcr->GetEntries() << endl;
		cout << " Reco Lb (weighted) = " << LbNum << endl;
		cout << " Weighted Lambda_b efficiency = " << effLambdab << " +- " << err_effLambdab << endl;
		//cout << " Weighted Lambda_b efficiency (using event weights for error) = " << effLambdab << " +- " << erreff << endl;
		cout << " Unweighted Lambda_b efficiency = " << h1mcr->GetEntries()/h1mcg->GetEntries()/preScale << endl;
		cout << "=========================================================="<<endl;

		//return;


		TH1F *bk1 = new TH1F("bk1","Xib BKGCAT",30,0,150);
		TH1F *m1 = new TH1F("m1","Mass",50,5700,5900);
		TH1F *m2 = new TH1F("m2","Mass",50,5700,5900);
		TH1F *m3 = new TH1F("m3","Mass",50,5700,5900);
		TH1F *m4 = new TH1F("m4","Mass",50,5700,5900);

		TLatex *myLatex = new TLatex();
		myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);

		TCanvas *c4 = new TCanvas("c4","",800,800);
		c4->Divide(2,2);
		c4->cd(1);
		chainx->Draw("Xib_BKGCAT>>bk1",Xib_Cut);
		addGraphics(bk1,"#Xi_{b} BKGCAT","",1);
		c4->cd(2);
		chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m1",Xib_Cut&&"Xib_BKGCAT<20");
		addGraphics(m1,"#Xi_{b} mass","",1);
		m1->Draw();
		myLatex->DrawLatex(0.18,0.85,"BKGCAT == 0");
		c4->cd(3);
		chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m2",Xib_Cut&&"Xib_BKGCAT==60");
		//chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m3",Xib_Cut&&"(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5795)<25)");
		addGraphics(m2,"#Xi_{b} mass","",1);
		//addGraphics(m3,"#Xi_{b} mass","",1);
		//m3->SetFillColor(4);
		m2->Draw();
		myLatex->DrawLatex(0.18,0.85,"BKGCAT == 60 (Ghost)");
		//m3->Draw("same");
		c4->cd(4);
		chainx->Draw("Xib_DTF_M_JpsiXiLConstr>>m4",Xib_Cut&&"!((Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5795)<25)||(Xib_BKGCAT<30))");
		addGraphics(m4,"#Xi_{b} mass","",1);
		m4->Draw();
		myLatex->DrawLatex(0.18,0.85,"BKGCAT == 50 (Low mass BG)");
		c4->Update();

		//return;

		TH1F *h2mcg = new TH1F("h2mcg","#Xi_{b}^{-} pT",100,0,50);
		TH1F *h2mcr = new TH1F("h2mcr","#Xi_{b}^{-} pT",100,0,50);
		TH1F *h2mcrL = new TH1F("h2mcrL","#Xi_{b}^{-} pT",100,0,50);
		TH1F *h2mcrD = new TH1F("h2mcrD","#Xi_{b}^{-} pT",100,0,50);
		h2mcr->Sumw2();
		h2mcrL->Sumw2();
		h2mcrD->Sumw2();
		h2mcg->Sumw2();


		TCanvas *c3 = new TCanvas("c3","",600,600);
		c3->Divide(2,2);
		c3->cd(1);


		//ofstream myfile1;
		//myfile1.open("test.dat",ios::out | ios::app);
		double wtot;

		if(mcOpt>=2)
		{
			int ib;
			double etaXib, ewt, etaBach, wXi;
			newtree = chainx->CopyTree(MCXib_Cut);
			cout << "Xib newtree = " << newtree << ", nentry = " << newtree->GetEntriesFast() << endl;
			setXibBranchAddresses();

			Long64_t nentries = newtree->GetEntriesFast();
			Long64_t ientry;
			//nentries = 1000;
			// Loop over reconstructed MC tuple
			for (Long64_t i=0; i<nentries; i++)
			{
				ientry = newtree->LoadTree(i);
				if (ientry < 0) break;
				newtree->GetEntry(i);
				if(i%10000 == 0) cout << "At entry = " << i << endl;
				wLb = 1.0; wL = 1.0; wPi = 1.0; wJpsi=1.0; wEta=1.0; wt = 1.0; wXi = 1.0;

				useTruth = true;
				if(Xib_TRUETAU==0 || Jpsi_TRUEPT==0 || L_TRUEPT==0 || pi_TRUEPT==0 || p_TRUEPT==0 || Xib_TRUEP_Z<=0 || Xib_TRUEPT<=0) useTruth = false;

				// Lifetime weight, for systematics
				if(Xib_TRUETAU>0)
				{
					wl = exp(-1000*Xib_TRUETAU/(1.566+XibLifetimeWeight*0.04))/exp(-1000*Xib_TRUETAU/1.566);
				}else{ wl = exp(-1000*Xib_TAU/(1.566+XibLifetimeWeight*0.04))/exp(-1000*Xib_TAU/1.566);}

				// Pion weight
				if(applyPrPiKinWeight)
				{
					if(mcOpt==4)
					{
						if(applyPiKinWeight && !useTruth) wPi = get2DKinWt(whistPi,pi_PT,pi_PZ,0);
						if(applyPiKinWeight && useTruth) wPi = get2DKinWt(whistPi,pi_TRUEPT,pi_TRUEP_Z,0);
					}
					else if(mcOpt==5)
					{
						if(!useTruth)
						{
							pPr = sqrt(p_PT**2+p_PZ**2);
							pPi = sqrt(pi_PT**2+pi_PZ**2);
						}
						else
						{
							pPr = sqrt(p_TRUEPT**2+p_TRUEP_Z**2);
							pPi = sqrt(pi_TRUEPT**2+pi_TRUEP_Z**2);
						}
						asym = (pPr-pPi) / (pPr+pPi);
						wPi = get1DKinWt(whistPrPi,asym);
					}
				}

				if(!useTruth)
				{
					jpsiP = Jpsi_P;
					XibP = Xib_P;
					etaXib = -log(tan(0.5*asin(Xib_PT/Xib_P)));
				}
				else
				{
					jpsiP = sqrt(Jpsi_TRUEPT**2+Jpsi_TRUEP_Z**2);
					XibP = sqrt(Xib_TRUEPT**2+Xib_TRUEP_Z**2);
					if(XibP<=100) XibP = Xib_P;
					if(jpsiP<=100) jpsiP = Jpsi_P;
					etaXib = -log(tan(0.5*atan(Xib_TRUEPT/Xib_TRUEP_Z)));
				}

				if(applyLbKinWeight && !useTruth) wLb = get2DKinWt(whistLb,Xib_PT,Xib_P,1);
				if(applyLbKinWeight && useTruth) wLb = get2DKinWt(whistLb,Xib_TRUEPT,Xib_TRUEP_Z,0);
				if(applyLKinWeight && !useTruth) wL = get2DKinWt(whistL,L_PT,L_P,1);
				if(applyLKinWeight && useTruth) wL = get2DKinWt(whistL,L_TRUEPT,L_TRUEP_Z,0);
				if(applyJpsiXiWeight && jpsiP>10 && XibP>10)
					wJpsi = (JpsiXiWeightPar[0]+JpsiXiWeightPar[1]*(jpsiP/XibP)+JpsiXiWeightPar[2]*(jpsiP/XibP)**2);
				if(applyEtaWeight) wEta = XibWtPar[0]+XibWtPar[1]*(etaXib-3.5);
				if(applypTWeight) wEta = XibpTWtPar[0]+0.0001*XibpTWtPar[1]*(Xib_PT-10000.0);
				if(applyXiWeight && !useTruth) wXi = (XiWeightPar[0]+XiWeightPar[1]*(Xi_PT/1000.0)+XiWeightPar[2]*(Xi_PT)**2);
				if(applyXiWeight && useTruth) wXi = (XiWeightPar[0]+XiWeightPar[1]*(Xi_TRUEPT/1000.0)+XiWeightPar[2]*(Xi_TRUEPT)**2);
				// Tracking correction
				if(BachPi_TRACK_Type==3)
				{
					wt = 1.00 + shift * 0.05;
					etaBach = -log(tan(0.5*asin(BachPi_PT/BachPi_P)));
					if(etaBach>1.9 && etaBach<4.9 && BachPi_P>5000 && BachPi_P<200000)
					{
						if(year=="Run1") ib = trackEff->FindBin(0.001*BachPi_P,etaBach); // p in GeV
						if(year=="Run2") ib = trackEff->FindBin(BachPi_P,etaBach); // p in MeV
						wt = trackEff->GetBinContent(ib);
						ewt = trackEff->GetBinError(ib);
						wt = wt + shift*ewt;
						if(wt == 0) cout << "Track with zero eff: " << BachPi_P << " " << etaBach << endl;
					}
				}
				else
				{
					// Downstream
					wt = 1.0;
				}
				if(wt<0 || wt>5) wt = 1.0;

				wtot = wl*wLb*wL*wPi*wJpsi*wt*wEta*wXi;
				if(i<100)
				{
					//if(eventNumber==800046 ||eventNumber==800176 ||eventNumber==800341 ||eventNumber==800409 ){        //myfile1 << useTruth << " " << wl << " " << wL << " " << wLb << " " << wPi << " " << wt << " " << wEta << " " << wJpsi << " " << jpsiP << " " << XibP << endl;
					cout << wtot << " " << Xib_PT << " " << wl << " " << wL << " " << wLb << " " << wPi << " " << wt << " " << wEta << " " << wJpsi << " " << wXi << endl;
				}

				if(wtot<=0 || wtot>50)
					cout << "Problem: " << wtot << " " << Xib_PT << " " << wl << " " << wL << " " << wLb << " " << wPi << " " << wt << " " << wEta << " " << wJpsi << " " << wXi << endl;
				if(BachPi_TRACK_Type==3) h2mcrL->Fill(0.001*Xib_PT,wtot); //wl*wLb*wL*wPi*wJpsi*wt*wEta*wXi);
				if(BachPi_TRACK_Type==5) h2mcrD->Fill(0.001*Xib_PT,wtot); //wl*wLb*wL*wPi*wJpsi*wt*wEta*wXi);
				if(BachPi_TRACK_Type==3) h2mcr->Fill(0.001*Xib_PT,wtot); //wl*wLb*wL*wPi*wJpsi*wt*wEta*wXi);
				if(BachPi_TRACK_Type==5) h2mcr->Fill(0.001*Xib_PT,wtot); //wl*wLb*wL*wPi*wJpsi*wt*wEta*wXi);
				//cout << i << " " << h2mcr->Integral() << endl;
			}

			// Loop over generated MC tuple
			nentries = chainxt->GetEntriesFast();
			cout << "Number of gen'd events in MC = " << nentries << endl;

			for (Long64_t i=0; i<nentries; i++)
			{
				if(preScale > 1 && i%preScale != 0) continue;
				ientry = chainxt->LoadTree(i);
				if (ientry < 0) break;
				chainxt->GetEntry(i);
				if(i%10000 == 0) cout << "At entry = " << i << endl;
				if(chargeCut>0 && J_psi_1S_MC_MOTHER_ID<0) continue;
				if(chargeCut<0 && J_psi_1S_MC_MOTHER_ID>0) continue;
				wLb = 1.0; wL = 1.0; wPi = 1.0; wJpsi=1.0; wEta=1.0; wXi=1.0;

				// Lifetime weight
				wl = exp(-1000*Xi_bminus_TRUETAU/(1.566+XibLifetimeWeight*0.04))/exp(-1000*Xi_bminus_TRUETAU/1.566);

				// Pion weight
				if(applyPrPiKinWeight)
				{
					if(mcOpt==4)
					{
						if(applyPiKinWeight) wPi = get2DKinWt(whistPi,piminus_TRUEPT,piminus_TRUEP_Z,0);
					}
					else if(mcOpt==5)
					{
						pPr = sqrt(pplus_TRUEPT**2+pplus_TRUEP_Z**2);
						pPi = sqrt(piminus_TRUEPT**2+piminus_TRUEP_Z**2);
						asym = (pPr-pPi) / (pPr+pPi);
						wPi = get1DKinWt(whistPrPi,asym);
					}
				}

				jpsiP = sqrt(J_psi_1S_TRUEPT**2+J_psi_1S_TRUEP_Z**2);
				XibP = sqrt(Xi_bminus_TRUEPT**2+Xi_bminus_TRUEP_Z**2);

				if(applyLbKinWeight) wLb = get2DKinWt(whistLb,Xi_bminus_TRUEPT,Xi_bminus_TRUEP_Z,0);
				if(applyLKinWeight) wL = get2DKinWt(whistL,Lambda0_TRUEPT,Lambda0_TRUEP_Z,0);
				if(applyJpsiXiWeight) wJpsi = (JpsiXiWeightPar[0]+JpsiXiWeightPar[1]*(jpsiP/XibP)+JpsiXiWeightPar[2]*(jpsiP/XibP)**2);
				if(applyXiWeight) wXi = (XiWeightPar[0]+XiWeightPar[1]*(Ximinus_TRUEPT/1000.0)+XiWeightPar[2]*(Ximinus_TRUEPT)**2);

				etaXib = -log(tan(0.5*atan(Xi_bminus_TRUEPT/Xi_bminus_TRUEP_Z)));
				if(applyEtaWeight) wEta = XibWtPar[0]+XibWtPar[1]*(etaXib-3.5);
				if(applypTWeight) wEta = XibpTWtPar[0]+0.0001*XibpTWtPar[1]*(Xi_bminus_TRUEPT-10000.0);
				if(eventNumber==800046 ||eventNumber==800176 ||eventNumber==800341 ||eventNumber==800409 )
				{
					//myfile1 << useTruth << " " << wl << " " << wL << " " << wLb << " " << wPi << " " << wt << " " << wEta << " " << wJpsi << " " << jpsiP << " " << XibP << endl;
					cout << useTruth << " " << eventNumber << " " << wl << " " << wL << " " << wLb << " " << wPi << " " << wt << " " << wEta << " " << wJpsi << endl;
				}

				h2mcg->Fill(0.001*Xi_bminus_TRUEPT,wl*wLb*wL*wPi*wJpsi*wEta*wXi*preScale);
			}
		}
		else
		{
			chainx->Draw("0.001*Xib_PT>>h2mcr",MCXib_Cut);
			c3->cd(2);
			chainx->Draw("0.001*Xib_PT>>h2mcrL",MCXib_Cut&&XiL_Cut);
			c3->cd(3);
			chainx->Draw("0.001*Xib_PT>>h2mcrD",MCXib_Cut&&XiD_Cut);
			chainxt->Draw("0.001*Xi_bminus_TRUEPT>>h2mcg");
		}



		XibNum = h2mcr->Integral();
		XibNumL = h2mcrL->Integral();
		XibNumD = h2mcrD->Integral();
		XibDen = h2mcg->Integral();

		uLbNum = h1mcr->GetEntries();
		uLbDen = h1mcg->GetEntries()*preScale;
		uXibNum = h2mcr->GetEntries();
		uXibNumL = h2mcrL->GetEntries();
		uXibNumD = h2mcrD->GetEntries();
		uXibDen = h2mcg->GetEntries()*preScale;

		double ueffLb = uLbNum/uLbDen;
		double ueffXib = uXibNum/uXibDen;
		double ueffXibL = uXibNumL/uXibDen;
		double ueffXibD = uXibNumD/uXibDen;

		c3->cd(1);
		h2mcg->Draw();
		c3->cd(2);
		h2mcr->Draw();
		c3->cd(3);
		h2mcrL->Draw();
		c3->cd(4);
		h2mcrD->Draw();
		c3->Update();

		cout << "=========================================================="<<endl;
		cout << " Weighted Lambda_b   efficiency = " << effLambdab << " +- " << err_effLambdab << endl;
		cout << "=========================================================="<<endl;
		cout << " Unweighted Lambda_b efficiency = " << ueffLb << endl;
		cout << " Unweighted Xi_b     efficiency  (L+D)= " << ueffXib << endl;
		cout << " Unweighted Xi_b     efficiency  (L)= " << ueffXibL << endl;
		cout << " Unweighted Xi_b     efficiency  (D)= " << ueffXibD << endl;
		cout << "=========================================================="<<endl;

		cout << chainxt->GetEntries() << endl;
		double effXib = XibNum/XibDen;
		double effXibL = XibNumL/XibDen;
		double effXibD = XibNumD/XibDen;
		double err_effXib = sqrt(effXib*(1-effXib)/ XibDen);
		double err_effXibL = sqrt(effXibL*(1-effXibL)/ XibDen);
		double err_effXibD = sqrt(effXibD*(1-effXibD)/ XibDen);

		/*
		   double errXib_Long  = geteff(h2mcrL,h2mcg);
		   double errXib_Down  = geteff(h2mcrD,h2mcg);
		   cout << "*********************************************************************"<<endl;
		   cout << "Error in Xib Eff " << errXib_Long << " " << errXib_Down << endl;
		   cout << "*********************************************************************"<<endl;
		 */

		//double corr = (2.0/3.0)*(tauLb/tauXib)*(genEffLb/genEffXib);
		//double err_corr = corr*sqrt((err_tauLb/tauLb)**2 +(err_tauXib/tauXib)**2 + (err_genEffLb/genEffLb)**2 + (err_genEffXib/genEffXib)**2);
		double corr = (2.0/3.0)*(tauLb/tauXib);//*(genEffLb/genEffXib);
		double err_corr = corr*sqrt((err_tauLb/tauLb)**2 +(err_tauXib/tauXib)**2 + (err_genEffXib/genEffXib)**2);

		if(year=="Run1")
		{
			lumJpsiL0 = lumJpsiL0_Run1;
			lumJpsiXi = lumJpsiXi_Run1;
		}
		else if(year=="Run2")
		{
			lumJpsiL0 = lumJpsiL0_2016;
			lumJpsiXi = lumJpsiXi_2016;
		}

		double lumRatio = lumJpsiXi / lumJpsiL0;

		double yieldRatio = yieldXib / nsigLb;
		double effRatio = effXib / effLambdab;
		double result = yieldRatio / (effRatio * relacc);
		double err_yieldRatio = yieldRatio*sqrt((err_yieldXib/yieldXib)**2 + (err_nsigLb/nsigLb)**2);
		double err_result = err_yieldRatio / (effRatio * relacc);

		result = result / lumRatio;
		err_result = err_result / lumRatio;

		double frag =  result * corr;
		double err_frag = frag*sqrt((err_result/result)**2 + (err_corr/corr)**2);
		double err_effRatio = effRatio*sqrt( (err_effXib/effXib)**2 + (err_effLambdab/effLambdab)**2);

		cout << "Long + Downstream Pion from Xib" << endl;
		cout << "-------------------------------" << endl;
		cout << "Yield Xib =     " << yieldXib << " +- " << err_yieldXib << endl;
		cout << "Yield Lb =     " << nsigLb << " +- " << err_nsigLb << endl;
		cout << "Yield ratio =     " << yieldRatio << " +- " << err_yieldRatio << endl;
		cout << "Xib efficiency = " << effXib << "+-" << err_effXib << endl;
		cout << "Eff   ratio =     " << effRatio << " +- " << err_effRatio << endl;
		cout << "Corr yield ratio: " << result << " +- " << err_result << endl;
		cout << "f(Xib)/f(Lb)    : " << frag << " +- " << err_frag << "+-" << 0.3*frag << endl;

		double yieldRatioL = yieldXibL / nsigLb;
		double effRatioL = effXibL / effLambdab;
		double resultL = yieldRatioL / (effRatioL * relacc);
		double err_yieldRatioL = yieldRatioL*sqrt((err_yieldXibL/yieldXibL)**2 + (err_nsigLb/nsigLb)**2);
		double err_resultL = err_yieldRatioL / (effRatioL * relacc);
		double err_effRatioL = effRatioL*sqrt( (err_effXibL/effXibL)**2 + (err_effLambdab/effLambdab)**2);

		resultL = resultL / lumRatio;
		err_resultL = err_resultL / lumRatio;

		cout << endl;

		double fragL =  resultL * corr;
		double err_fragL = fragL*sqrt((err_resultL/resultL)**2 + (err_corr/corr)**2);

		cout << "Long Pion from Xib" << endl;
		cout << "------------------" << endl;
		cout << "Yield Xib =     " << yieldXibL << " +- " << err_yieldXibL << endl;
		cout << "Yield Lb =     " << nsigLb << " +- " << err_nsigLb << endl;
		cout << "Yield ratio = " << yieldRatioL << " +- " << err_yieldRatioL << endl;
		cout << "Xib efficiency = " << effXibL << "+-" << err_effXibL << endl;
		cout << "Eff   ratio = " << effRatioL << " +- " << err_effRatioL << endl;
		cout << "Corr yield ratio: " << resultL << " +- " << err_resultL << endl;
		cout << "f(Xib)/f(Lb)    : " << fragL << " +- " << err_fragL << "+-" << 0.3*fragL << endl;
		cout << endl;



		double yieldRatioD = yieldXibD / nsigLb;
		double effRatioD = effXibD / effLambdab;
		double resultD = yieldRatioD / (effRatioD * relacc);
		double err_yieldRatioD = yieldRatioD*sqrt((err_yieldXibD/yieldXibD)**2 + (err_nsigLb/nsigLb)**2);
		double err_resultD = err_yieldRatioD / (effRatioD * relacc);

		resultD = resultD / lumRatio;
		err_resultD = err_resultD / lumRatio;

		double fragD =  resultD * corr;
		double err_fragD = fragD*sqrt((err_resultD/resultD)**2 + (err_corr/corr)**2);
		double err_effRatioD = effRatioD*sqrt( (err_effXibD/effXibD)**2 + (err_effLambdab/effLambdab)**2);

		cout << "Downstream Pion from Xib" << endl;
		cout << "------------------------" << endl;
		cout << "Yield Xib =     " << yieldXibD << " +- " << err_yieldXibD << endl;
		cout << "Yield Lb =     " << nsigLb << " +- " << err_nsigLb << endl;
		cout << "Yield ratio = " << yieldRatioD << " +- " << err_yieldRatioD << endl;
		cout << "Xib efficiency = " << effXibD << "+-" << err_effXibD << endl;
		cout << "Eff   ratio = " << effRatioD << " +- " << err_effRatioD << endl;
		cout << "Corr yield ratio: " << resultD << " +- " << err_resultD << endl;
		cout << "f(Xib)/f(Lb)    : " << fragD << " +- " << err_fragD << "+-" << 0.3*fragD << endl;

		cout << endl;
		cout << "Fitted masses" << endl;
		cout << "   Lambda_b mass(DD)   : " << fitmass_Lb_DD << " +- " << err_fitmass_Lb_DD << endl;
		cout << "   Lambda_b mass(LL)   : " << fitmass_Lb_LL << " +- " << err_fitmass_Lb_LL << endl;
		cout << "   Xi_b mass(L+D)  : " << fitmass_Xib << " +- " << err_fitmass_Xib << endl;
		cout << "   Xi_b mass(L)    : " << fitmassL_Xib << " +- " << err_fitmassL_Xib << endl;
		cout << "   Xi_b mass(D)    : " << fitmassD_Xib << " +- " << err_fitmassD_Xib << endl;

		cout << endl;

		cout << "****************************************************************************" << endl;
		cout << "*** Prescale factor of:    " << preScale <<  "  is applied and corrected for" << endl;
		cout << "****************************************************************************" << endl;


		cout << "****************************************************************************" << endl;
		cout << "Raw yields of events from MC" << endl;
		cout << "Lb Rec & Gen << " << h1mcr->Integral() << " " << h1mcg->Integral() << endl;
		cout << "Xib Rec & Gen (DDL) << " << h2mcrL->Integral() << " " << h2mcg->Integral() << endl;
		cout << "Xib Rec & Gen (DDD) << " << h2mcrD->Integral() << " " << h2mcg->Integral() << endl;
		cout << "****************************************************************************" << endl;

	}


	if(do_splot == 0) return;

	TString sLb_Trig = Lb_Trig;
	TString tag = "DD";
	if(do_splot == 2) tag = "LL";

	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	if(sLb_Trig.Contains("L0Global")>0 && sLb_Trig.Contains("L0DiMuon")>0)
	{
		tag = tag + "_L0JPsiTOSorGlobalTIS";
		cout << "Trigger is (JpsiL0DiMuon TOS || Lb L0 Global TIS)" << endl;
	}
	else if(sLb_Trig.Contains("L0DiMuon")>0 && sLb_Trig.Contains("L0Global")==0  )
	{
		tag = tag + "_L0JPsiTOS";
		cout << "Trigger is (JpsiL0DiMuon TOS)" << endl;
	}
	else if(sLb_Trig.Contains("L0DiMuon")==0 && sLb_Trig.Contains("L0Global")>0 )
	{
		tag = tag + "_L0GlobalTIS";
		cout << "Trigger is (Lb L0 Global TIS)" << endl;
	}
	else
	{
		cout << "Shoot, something wrong here!!" << endl;
	}
	cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;

	if(do_splot == 1 || do_splot==2) {


		TString fileName = "Lb2JpsiL0_"+year+"_"+tag+".root";
		TFile *fout = new TFile(fileName,"RECREATE");

		if(do_splot==1) TTree* treeNew1 = treeDD->CopyTree(LambdabDD_Cut);
		if(do_splot==2) TTree* treeNew1 = treeLL->CopyTree(LambdabLL_Cut);

		//--------------------------------------------
		// Variables that are needed in the sWeighting
		//--------------------------------------------
		RooRealVar* L_ENDVERTEX_Z = new RooRealVar("L_ENDVERTEX_Z","Lambda Z",-1.0e10,1.0e10);
		RooRealVar* L_TAU = new RooRealVar("L_TAU","Lambda decay time",-1.0e10,1.0e10);
		RooRealVar* L_P = new RooRealVar("L_P","Lambda mom",-1.0e10,1.0e10);
		RooRealVar* L_PT = new RooRealVar("L_PT","Lambda mom",-1.0e10,1.0e10);
		RooRealVar* Lb_p = new RooRealVar("Lb_P","Lambda mom",-1.0e10,1.0e10);
		RooRealVar* Lb_pT = new RooRealVar("Lb_PT","Lambda mom",-1.0e10,1.0e10);
		RooRealVar* Lb_OWNPV_Z = new RooRealVar("Lb_OWNPV_Z","Lambdab origin Z",-1.0e10,1.0e10);
		RooRealVar* muplus_pT = new RooRealVar("muplus_PT","mu+ pT",-1.0e10,1.0e10);
		RooRealVar* muminus_pT = new RooRealVar("muminus_PT","mu- pT",-1.0e10,1.0e10);
		RooRealVar* Jpsi_pT = new RooRealVar("Jpsi_PT","Jpsi pT",-1.0e10,1.0e10);
		RooRealVar* Jpsi_p = new RooRealVar("Jpsi_P","Jpsi mom",-1.0e10,1.0e10);
		RooRealVar* L_FD_OWNPV = new RooRealVar("L_FD_OWNPV","Lambda FD",-1.0e10,1.0e10);
		RooRealVar* p_IP_OWNPV = new RooRealVar("p_IP_OWNPV","proton IP",-1.0e10,1.0e10);
		RooRealVar* pi_IP_OWNPV = new RooRealVar("pi_IP_OWNPV","pion IP",-1.0e10,1.0e10);
		RooRealVar* p_IPCHI2_OWNPV = new RooRealVar("p_IPCHI2_OWNPV","proton IP",-1.0e10,1.0e10);
		RooRealVar* pi_IPCHI2_OWNPV = new RooRealVar("pi_IPCHI2_OWNPV","pion IP",-1.0e10,1.0e10);
		RooRealVar* pi_ORIVX_X = new RooRealVar("pi_ORIVX_X","pion X prod",-1.0e10,1.0e10);
		RooRealVar* pi_ORIVX_Y = new RooRealVar("pi_ORIVX_Y","pion Y prod",-1.0e10,1.0e10);
		RooRealVar* pi_ORIVX_Z = new RooRealVar("pi_ORIVX_Z","pion Z prod",-1.0e10,1.0e10);
		RooRealVar* pi_PX = new RooRealVar("pi_PX","pion PX",-1.0e10,1.0e10);
		RooRealVar* pi_PY = new RooRealVar("pi_PY","pion PY",-1.0e10,1.0e10);
		RooRealVar* pi_PZ = new RooRealVar("pi_PZ","pion PZ",-1.0e10,1.0e10);
		RooRealVar* p_PX = new RooRealVar("p_PX","proton PX",-1.0e10,1.0e10);
		RooRealVar* p_PY = new RooRealVar("p_PY","proton PY",-1.0e10,1.0e10);
		RooRealVar* p_PZ = new RooRealVar("p_PZ","proton PZ",-1.0e10,1.0e10);
		RooRealVar* pi_ID = new RooRealVar("pi_ID","pion PZ",-1.0e10,1.0e10);
		RooRealVar* runNumber = new RooRealVar("runNumber","run #",-1.0e10,1.0e10);

		RooArgSet args;
		args.add(*mvar1);
		args.add(*L_TAU);
		args.add(*L_P);
		args.add(*L_PT);
		args.add(*L_ENDVERTEX_Z);
		args.add(*Lb_p);
		args.add(*Lb_pT);
		args.add(*Lb_OWNPV_Z);
		args.add(*Jpsi_pT);
		args.add(*muplus_pT);
		args.add(*muminus_pT);
		args.add(*Jpsi_p);
		args.add(*L_FD_OWNPV);
		args.add(*p_IP_OWNPV);
		args.add(*pi_IP_OWNPV);
		args.add(*p_IPCHI2_OWNPV);
		args.add(*pi_IPCHI2_OWNPV);
		args.add(*pi_ORIVX_X);
		args.add(*pi_ORIVX_Y);
		args.add(*pi_ORIVX_Z);
		args.add(*pi_PX);
		args.add(*pi_PY);
		args.add(*pi_PZ);
		args.add(*pi_ID);
		args.add(*p_PX);
		args.add(*p_PY);
		args.add(*p_PZ);
		args.add(*runNumber);

		RooDataSet* data_red = new RooDataSet("data_red","dataset with reduced vars",treeNew1,args);

		TString wsig = "nsig1_sw";

		RooArgSet *myVars1;
		if(do_splot==1)
		{
			pdf1->fitTo(*rdh1,Hesse(kTRUE),Strategy(2));
			myVars1 = pdf1->getParameters(*rdh1);
			myVars1->Print("v");
			myVars1->setAttribAll("Constant",kTRUE);
			nsig1->setConstant(kFALSE);
			ncomb1->setConstant(kFALSE);
			pdf1->fitTo(*rdh1,Extended(),Hesse(kTRUE));
		}
		else
		{
			pdf1->fitTo(*rdh1LL,Hesse(kTRUE),Strategy(2));
			myVars1 = pdf1->getParameters(*rdh1LL);
			myVars1->Print("v");
			myVars1->setAttribAll("Constant",kTRUE);
			nsig1->setConstant(kFALSE);
			ncomb1->setConstant(kFALSE);
			pdf1->fitTo(*rdh1LL,Extended(),Hesse(kTRUE));
		}

		//--------------
		// Get sWeights
		//--------------
		RooStats::SPlot* sDataA1 = new RooStats::SPlot("sDataA1","An SPlot", *data_red, pdf1, RooArgList(*nsig1,*ncomb1) );

		data_red->Print("v");
		gStyle->SetOptStat(011000);
		//--------------------------
		// create weighted data set
		//--------------------------
		RooDataSet * dataw_signal = new RooDataSet(data_red->GetName(),data_red->GetTitle(),data_red,*data_red->get(),0,wsig);
		dataw_signal->Print("v");
		//-------------------------------------------
		// Project out weighted data onto histograms
		//-------------------------------------------

		TTree *newtree2 = dataw_signal->tree();
		newtree2->SetName("TreeWithWeight");
		newtree2->Write();


		fout->Write();

	}
	else
	{
		//"Lb2JpsiL0_"+year+"_"+tag+".root";

		TString fileName = "Xib2JpsiXi_"+year+"_"+tag;
		if(useNewStrip) fileName = fileName + "_NewStrip"; ;
		if(applyXiPiTrackingAcceptanceCut) fileName = fileName+"_withTrackingAcc";
		fileName = fileName +".root";


		TFile *fout = new TFile(fileName,"RECREATE");
		TCut longCut = Xib_Cut && XiL_Cut;

		TTree* treeNew1 = treeXi->CopyTree(longCut);

		//--------------------------------------------
		// Variables that are needed in the sWeighting
		//--------------------------------------------
		RooRealVar* L_TAU = new RooRealVar("L_TAU","Lambda decay time",-1.0e10,1.0e10);
		RooRealVar* L_P = new RooRealVar("L_P","Lambda mom",-1.0e10,1.0e10);
		RooRealVar* L_PT = new RooRealVar("L_PT","Lambda mom",-1.0e10,1.0e10);
		RooRealVar* Jpsi_pT = new RooRealVar("Jpsi_PT","Jpsi pT",-1.0e10,1.0e10);
		RooRealVar* muplus_pT = new RooRealVar("muplus_PT","mu+ pT",-1.0e10,1.0e10);
		RooRealVar* muminus_pT = new RooRealVar("muminus_PT","mu- pT",-1.0e10,1.0e10);
		RooRealVar* Jpsi_p = new RooRealVar("Jpsi_P","Jpsi mom",-1.0e10,1.0e10);
		RooRealVar* L_FD_OWNPV = new RooRealVar("L_FD_OWNPV","Lambda FD",-1.0e10,1.0e10);
		RooRealVar* L_FD_ORIVX = new RooRealVar("L_FD_ORIVX","Lambda FD",-1.0e10,1.0e10);

		RooRealVar* Xi_TAU = new RooRealVar("Xi_TAU","Xi decay time",-1.0e10,1.0e10);
		RooRealVar* Xi_P = new RooRealVar("Xi_P","Xi mom",-1.0e10,1.0e10);
		RooRealVar* Xi_PT = new RooRealVar("Xi_PT","Xi mom",-1.0e10,1.0e10);
		RooRealVar* Xi_ENDVERTEX_CHI2 = new RooRealVar("Xi_ENDVERTEX_CHI2","Xi chisq_vtx",-1.0e10,1.0e10);
		RooRealVar* L_ENDVERTEX_Z = new RooRealVar("L_ENDVERTEX_Z","Lambda Z",-1.0e10,1.0e10);
		RooRealVar* Xi_DIRA_ORIVX = new RooRealVar("Xi_DIRA_ORIVX","Xi DIRA",-1.0e10,1.0e10);
		//RooRealVar* Xi_FDCHI2_ORIVX = new RooRealVar("Xi_DIRA_ORIVX","Xi DIRA",-1.0e10,1.0e10);


		RooRealVar* Xib_p = new RooRealVar("Xib_P","Xib mom",-1.0e10,1.0e10);
		RooRealVar* Xib_pT = new RooRealVar("Xib_PT","Xib pT",-1.0e10,1.0e10);
		RooRealVar* Xi_FD_OWNPV = new RooRealVar("Xi_FD_OWNPV","Xi FD",-1.0e10,1.0e10);
		RooRealVar* Xib_OWNPV_Z = new RooRealVar("Xib_OWNPV_Z","Xib origin Z",-1.0e10,1.0e10);


		RooRealVar* BachPi_pT = new RooRealVar("BachPi_PT","BachPi PT",-1.0e10,1.0e10);
		RooRealVar* BachPi_p = new RooRealVar("BachPi_P","BachPi P",-1.0e10,1.0e10);
		RooRealVar* BachPi_IPCHI2_OWNPV = new RooRealVar("BachPi_IPCHI2_OWNPV","BachPi chisqip",-1.0e10,1.0e10);

		RooRealVar* pi_PX = new RooRealVar("pi_PX","pion PX",-1.0e10,1.0e10);
		RooRealVar* pi_PY = new RooRealVar("pi_PY","pion PY",-1.0e10,1.0e10);
		RooRealVar* pi_PZ = new RooRealVar("pi_PZ","pion PZ",-1.0e10,1.0e10);
		RooRealVar* p_PX = new RooRealVar("p_PX","proton PX",-1.0e10,1.0e10);
		RooRealVar* p_PY = new RooRealVar("p_PY","proton PY",-1.0e10,1.0e10);
		RooRealVar* p_PZ = new RooRealVar("p_PZ","proton PZ",-1.0e10,1.0e10);


		RooArgSet args;
		args.add(*mvar2);
		args.add(*L_TAU);
		args.add(*L_P);
		args.add(*L_PT);
		args.add(*Jpsi_pT);
		args.add(*muplus_pT);
		args.add(*muminus_pT);
		args.add(*Jpsi_p);
		args.add(*L_FD_OWNPV);
		args.add(*L_FD_ORIVX);
		args.add(*Xi_TAU);
		args.add(*Xi_P);
		args.add(*Xi_PT);
		args.add(*Xib_p);
		args.add(*Xib_pT);
		args.add(*Xib_OWNPV_Z);
		args.add(*Xi_FD_OWNPV);
		args.add(*L_ENDVERTEX_Z);
		args.add(*Xi_ENDVERTEX_CHI2);
		args.add(*Xi_DIRA_ORIVX);
		args.add(*BachPi_pT);
		args.add(*BachPi_p);
		args.add(*BachPi_IPCHI2_OWNPV);
		args.add(*pi_PX);
		args.add(*pi_PY);
		args.add(*pi_PZ);
		args.add(*p_PX);
		args.add(*p_PY);
		args.add(*p_PZ);

		RooDataSet* data_red = new RooDataSet("data_red","dataset with reduced vars",treeNew1,args);

		TString wsig = "nsig2_sw";

		RooArgSet *myVars1;
		pdf2->fitTo(*rdh2L,Hesse(kTRUE),Strategy(2));
		myVars1 = pdf2->getParameters(*rdh2L);
		myVars1->Print("v");
		myVars1->setAttribAll("Constant",kTRUE);
		nsig2->setConstant(kFALSE);
		ncomb2->setConstant(kFALSE);
		pdf2->fitTo(*rdh2L,Extended(),Hesse(kTRUE),Strategy(2));
		//--------------
		// Get sWeights
		//--------------
		cout << "Going to call SPlot now"<< endl;
		RooStats::SPlot* sDataA1 = new RooStats::SPlot("sDataA1","An SPlot", *data_red, pdf2, RooArgList(*nsig2,*ncomb2) );
		cout << "Done calling SPlot now"<< endl;

		data_red->Print("v");
		gStyle->SetOptStat(011000);
		//--------------------------
		// create weighted data set
		//--------------------------
		RooDataSet * dataw_signal = new RooDataSet(data_red->GetName(),data_red->GetTitle(),data_red,*data_red->get(),0,wsig);
		dataw_signal->Print("v");
		//-------------------------------------------
		// Project out weighted data onto histograms
		//-------------------------------------------

		TTree *newtree2 = dataw_signal->tree();
		newtree2->SetName("TreeWithWeight");
		newtree2->Write();


		fout->Write();


	}





	return;

}
