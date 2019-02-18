#include "OptimizeFinalBDT.h"

Double_t getSumOfWeights(Double_t mybdt = 0., TTree* tree = nullptr,
                         Int_t flag = 1, Int_t bdtConf = 1, Int_t nEntries = 0);

std::vector <Double_t> OptimizeFinalBDT(Int_t run, Int_t trackType, const char* isoVersion,
                                        Int_t isoConf, Int_t bdtConf, Bool_t isoFlag,
                                        Bool_t logFlag, Bool_t newFlag, TString FOM)
{
	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");

	std::vector <Double_t> cuts = {-1.0};
	const char *mynew = "", *type = "";
	if(newFlag) mynew = "_new";

	type = (trackType == 3) ? "LL" : "DD";

	if(logFlag && isoFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalBDT%d_%s_iso%d_%s%s.txt",
		                             run,bdtConf,type,isoConf,isoVersion,mynew),"w");
	}
	if(logFlag && !isoFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/OptimizeFinalBDT%d_%s_noIso%s.txt",
		                             run,bdtConf,type,mynew),"w");
	}

	cout<<"******************************************"<<endl;
	cout<<"==> Starting OptimizeFinalBDT: "<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr, *fileIn_zeroTracks = nullptr;
	TTree *treeIn = nullptr, *treeIn_zeroTracks = nullptr, *myTree = nullptr;

	TString friendFileName = "";
	TString friendFileName_zeroTracks = "", t = "";

	TH1D *hsig = nullptr, *hbkg = nullptr;
	Int_t nEntries = 0;
	const char *rootFolder = "";

	rootFolder = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);

	if(isoFlag)
	{
		fileIn                    = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks.root",
		                                             rootFolder,type),"READ");
		if (!fileIn)
		{
			cout << "ERROR: could not open data file" << endl;
			exit(1);
		}
		treeIn                    = (TTree*)fileIn->Get("MyTuple");

		fileIn_zeroTracks         = TFile::Open(Form("%s/jpsilambda_%s_withsw_ZeroTracks.root",
		                                             rootFolder,type),"READ");
		if (!fileIn_zeroTracks)
		{
			cout << "ERROR: could not open zeroTracks data file" << endl;
			exit(1);
		}
		treeIn_zeroTracks         = (TTree*)fileIn_zeroTracks->Get("MyTuple");

		friendFileName            = Form("%s/jpsilambda_%ssig_FinalBDT%d_iso%d_%s%s.root",
		                                 rootFolder,type,bdtConf,isoConf,isoVersion,mynew);
		friendFileName_zeroTracks = Form("%s/jpsilambda_zeroTracks%ssig_FinalBDT%d_iso%d_%s%s.root",
		                                 rootFolder,type,bdtConf,isoConf,isoVersion,mynew);

		treeIn_zeroTracks->AddFriend("MyTuple",friendFileName_zeroTracks);
	}
	else
	{
		fileIn = TFile::Open(Form("%s/jpsilambda_%s_withsw.root",
		                          rootFolder,type));
		if (!fileIn)
		{
			cout << "ERROR: could not open data file" << endl;
			exit(1);
		}
		treeIn = (TTree*)fileIn->Get("MyTuple");

		friendFileName = Form("%s/jpsilambda_%ssig_FinalBDT%d_noIso.root",
		                      rootFolder,type,bdtConf);
	}
	treeIn->AddFriend("MyTuple",friendFileName);

	Int_t ctr = 0;
	TList *list = new TList();
	list->Add(treeIn);

	if(isoFlag)
	{
		list->Add(treeIn_zeroTracks);
	}

	TIterator *iter = (TIterator*)list->MakeIterator();
	while((myTree = (TTree*)iter->Next()))
	{
		Double_t BDT = 0., BDT_max = 0.;
		Double_t myFOM = 0., myFOM_max = 0.;
		Double_t eff_sig = 0., eff_sig_max = 0.;
		Double_t eff_bkg = 0., eff_bkg_max = 0.;
		Double_t sig = 0., bkg = 0., siginit = 0.;
		Int_t bkginit = 0;

		if(ctr == 0)
		{
			cout<<"****************Optimizing for nonZeroTracks****************"<<endl;
		}
		else if(ctr == 1)
		{
			cout<<"****************Optimizing for ZeroTracks****************"<<endl;
		}
		nEntries = myTree->GetEntries();
		cout<<"nEntries= "<<nEntries<<endl;

		myTree->Draw(Form("BDT%d>>hsig0(200,-1.0,1.0)",bdtConf),"SW","goff");
		// myTree->Draw(Form("BDT%d>>hsig(90,-0.4,0.5)",bdtConf),"SW*(Lb_DTF_M_JpsiLConstr < 5700)","goff");//Should there be a lower bound as well?
		// myTree->Draw(Form("BDT%d>>hsig1(90,-0.4,0.5)",bdtConf),"SW*(Lb_DTF_M_JpsiLConstr < 5700 && Lb_DTF_M_JpsiLConstr > 5300)","goff");
		// myTree->Draw(Form("BDT%d>>hbkg(200,-1.0,1.0)",bdtConf),"BW*(Lb_DTF_M_JpsiLConstr > 5700)","goff");

		hsig = (TH1D*)gDirectory->Get("hsig0");  //play around with this
		//      hbkg = (TH1D*)gDirectory->Get("hbkg");

		TFile *tempFile = new TFile(Form("rootFiles/dataFiles/JpsiLambda/run%d/tempFile.root",run),"RECREATE");
		TTree *tempTree = (TTree*)myTree->CopyTree("");

		bkginit = myTree->GetEntries("Lb_DTF_M_JpsiLConstr > 5700")*26/30;
		siginit =  (Int_t)getSumOfWeights(-0.5,tempTree,1,bdtConf,nEntries);//rounding off happening here
		//siginit = 6575;
		//      bkginit =  (Int_t)getSumOfWeights(-0.5,myTree,2,bdtConf,nEntries);
		cout<<"siginit = "<<siginit<<endl;
		cout<<"bkginit = "<<bkginit<<endl;


		for(Int_t i = 1; i < 200; i++) {

			BDT = hsig->GetBinCenter(i);

			t = Form("BDT > %f",BDT);

			// sig = (Int_t)getSumOfWeights(BDT,myTree,1,bdtConf,nEntries);
			// bkg = (Int_t)getSumOfWeights(BDT,myTree,2,bdtConf,nEntries);

			sig = hsig->Integral(i,200);
			//	bkg = hbkg->Integral(i,200);
			bkg = myTree->GetEntries(Form("Lb_DTF_M_JpsiLConstr > 5700 && BDT%d > %f",bdtConf,BDT))*26/30;
			cout<<"SIG = "<<sig<<" BKG = "<<bkg;
			eff_sig = (Double_t) sig/siginit;
			eff_bkg = (Double_t) bkg/bkginit;

			if(sig + bkg > 0)
			{
				if(FOM == "sig")
				{
					myFOM = (Double_t)eff_sig/sqrt(sig+bkg);
				}
				if(FOM == "punzi")
				{
					myFOM = (Double_t)eff_sig/(sqrt(bkg) + 1.5);//a  = 3 chosen
				}
			}

			else
			{
				cout<<"NO"<<endl;
				continue;
			}
			if(myFOM > myFOM_max) {
				myFOM_max = myFOM;
				BDT_max = BDT;
				eff_sig_max = eff_sig;
				eff_bkg_max = eff_bkg;
			}

			cout<<"For BDT = "<<BDT<<" FOM = "<<myFOM<<" sig_eff = "
			    <<eff_sig*100<<"% bkg_eff = "<<eff_bkg*100<<"%"<<endl;
		}
		cout<<"MAXIMUM FOM = "<<myFOM_max<<" at BDT = "<<BDT_max<<" with sig_eff = "
		    <<eff_sig_max*100<<"% and bkg_eff = "<<eff_bkg_max*100<<"%"<<endl;
		cuts.push_back(BDT_max);
		ctr++;
	}
	if(logFlag) gSystem->RedirectOutput(0);
	return cuts;
}
Double_t getSumOfWeights(Double_t mybdt, TTree* tree, Int_t flag,
                         Int_t bdtConf, Int_t nEntries)
{
	Double_t sum = 0., mybdt1 = 0., bmass = 0.;
	Float_t myweight = 0.;
	// TTree *treecut = (TTree*)tree->CopyTree(Form("BDT > %f",mybdt));

	if(flag == 1)
	{
		tree->SetBranchAddress("SW",&myweight);
	}
	else if(flag == 2)
	{
		tree->SetBranchAddress("BW",&myweight);
	}

	tree->SetBranchAddress("Lb_DTF_M_JpsiLConstr",&bmass);
	tree->SetBranchAddress(Form("BDT%d",bdtConf),&mybdt1);

	for(Int_t i=0; i<nEntries; i++)
	{
		tree->GetEntry(i);

		if(mybdt1 > mybdt)
		{
			if((flag == 2 && bmass > 5700) || flag == 1)
			{
				sum += myweight;
			}
		}
	}
	cout<<"mybdt = "<<mybdt<<" sum = "<<sum<<endl;
	return sum;
}
