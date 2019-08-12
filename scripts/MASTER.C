#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
// #include "Trigger.h"
// #include "Sanity.h"
// #include "CutOutKs.h"
// #include "DoSWeight.h"
// #include "DoSWeight_Sanity.h"
// #include "TrainIsolation.h"
// #include "ApplyIsolation.h"
// #include "TrainFinalBDT.h"
// #include "ApplyFinalBDT.h"
// #include "OptimizeFinalBDT.h"
// #include "CutFinalBDT.h"
#include <iostream>
#include <vector>

using namespace std;

void MASTER(Int_t run = 1, Int_t year = 2015, Int_t config = 1, Int_t block = 1, Bool_t isData = true)
//MCType = 1 for JpsiLambda MC
//MCType = 2 for JpsiSigma MC
//MCType = 3 for JpsiXi MC (reco'd JpsiLambda)
//MCType = 4 for Bu_JpsiX MC
//MCType = 5 for Bd_JpsiX MC
//MCType = 6 for Lambda*(1405)MC
//MCType = 7 for Lambda*(1520)MC
//MCType = 8 for Lambda*(1600)MC
//MCType = 9 for chiC1 MC
//MCType = 10 for JpsiKs MC
//MCType = 11 for Xib0 MC
{
	TStopwatch sw;
	sw.Start();

	gSystem->Exec("date");

	// gROOT->ProcessLine(".L CollateFiles.C");
	// gROOT->ProcessLine(".L Trigger.C");
	// gROOT->ProcessLine(".L Sanity.C");
	// gROOT->ProcessLine(".L CutOutKs.C");
	// gROOT->ProcessLine(".L DoSWeight.C");
	// gROOT->ProcessLine(".L DoSWeight_Sanity.C");
	// gROOT->ProcessLine(".L TrainIsolation.C");
	// gROOT->ProcessLine(".L ApplyIsolation.C");
	// gROOT->ProcessLine(".L TrainFinalBDT.C");
	// gROOT->ProcessLine(".L ApplyFinalBDT.C");
	// gROOT->ProcessLine(".L OptimizeFinalBDT.C");
	// gROOT->ProcessLine(".L CutFinalBDT.C+");

	// gSystem->Load("CollateFiles_C.so");
	// gSystem->Load("Trigger_C.so");
	// gSystem->Load("Sanity_C.so");
	// gSystem->Load("CutOutKs_C.so");
	// gSystem->Load("DoSWeight_C.so");
	// gSystem->Load("DoSWeight_Sanity_C.so");
	// gSystem->Load("TrainIsolation_C.so");
	// gSystem->Load("ApplyIsolation_C.so");
	// gSystem->Load("TrainFinalBDT_C.so");
	// gSystem->Load("ApplyFinalBDT_C.so");
	// gSystem->Load("OptimizeFinalBDT_C.so");
	// gSystem->Load("CutFinalBDT_C.so");

	cout<<"poop"<<endl;

	Int_t mcType        = 0;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType     = 3;// 3 for LL, 5 for DD.
	Int_t isoConf       = 1;// config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf  = 1;// config for final BDT. Currently 1 or 2 supported

	Bool_t testing  = false;// when true, analysis will only run over a subset of data
	Bool_t loose    = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Bool_t isoFlag  = true;// when true, isolation will be used in final BDT.
	Bool_t logFlag  = true;// set to false only while testing.
	Bool_t simFlag  = false;

	Float_t bdtCut            = 0.0;
	Float_t bdtCut_ZeroTracks = 0.0;

	Double_t myChi2_nonZero = 0.0, myChi2_Zero = 0.0;
	Double_t myChi2_Sanity  = 0.0;

	const char* isoVersion = ""; // which version of isolation BDT? changes input variables used in isolation. v0,1 supported.

	const char* FOM = "Punzi";//Sig for S/sqrt(S+B), Punzi for Punzi FOM with a = 3

	std::vector <Double_t> bdtCuts;
	Int_t len = 11;

	if(isData)
	{
		len = 1;
	}

	switch(config)
	{
	case 1:
	{
		isoVersion   = "v0";
		isoConf      = 1;
		finalBDTconf = 1;
		break;
	}
	case 2:
	{
		isoVersion   = "v0";
		isoConf      = 1;
		finalBDTconf = 2;
		break;
	}
	case 3:
	{
		isoVersion   = "v0";
		isoConf      = 2;
		finalBDTconf = 1;
		break;
	}
	case 4:
	{
		isoVersion   = "v0";
		isoConf      = 2;
		finalBDTconf = 2;
		break;
	}
	case 5:
	{
		isoVersion   = "v1";
		isoConf      = 1;
		finalBDTconf = 1;
		break;
	}
	case 6:
	{
		isoVersion   = "v1";
		isoConf      = 1;
		finalBDTconf = 2;
		break;
	}
	case 7:
	{
		isoVersion   = "v1";
		isoConf      = 2;
		finalBDTconf = 1;
		break;
	}
	case 8:
	{
		isoVersion   = "v1";
		isoConf      = 2;
		finalBDTconf = 2;
		break;
	}
	}
	for(mcType = 1; mcType <=len; mcType++)
	{
		if(isData) mcType = 0;
		if(mcType != 4 && mcType != 5 && mcType != 9 && mcType !=11)
		{
			cout<<"$$$$$$$$$$$ Processing MC Type "<<mcType<< " $$$$$$$$$$$$$$"<<endl;
			cout<<"$$$$$$$$$$$ Processing Run "<<run<<" $$$$$$$$$$$"<<endl;
			if(block == 0)
			{
				if(!isData) //MC
				{
					cout<<"****Collating MC PIDGEN files"<<endl;
					gSystem->Exec(Form("root -l -b -q \'CollateFiles.C(%d,%d,%d,%d)\'",run, year, isData, mcType));

					if(mcType!=10)
					{
						cout<<"****Applying GB weights on reco MC"<<endl;
						gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
						gSystem->Exec(Form("python ApplyGBWeights.py %d %d 0",run,mcType));
						cout<<"****Applying GB weights on generated MC"<<endl;
						gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
						gSystem->Exec(Form("python ApplyGBWeights.py %d %d 1",run,mcType));
					}
				}
				//Trigger Cut
				cout<<"***Trigger***"<<endl;
				// Trigger(run, year, isData, mcType, testing, loose, logFlag);
				gSystem->Exec(Form("root -l -b -q \'Trigger.C(%d,%d,%d,%d,%d,%d,%d)\'",run, year, isData, mcType, testing, loose, logFlag));
				//Sanity Cuts
				cout<<"***Sanity***"<<endl;
				// Sanity(run, year, isData, mcType, logFlag);
				gSystem->Exec(Form("root -l -q -b \'Sanity.C(%d,%d,%d,%d,%d)\'",run, year, isData, mcType, logFlag));
				//Cut Out Ks0
				cout<<"***CutOutKs***"<<endl;
				// CutOutKs(run, year, isData, mcType, trackType, logFlag);
				gSystem->Exec(Form("root -l -q -b \'CutOutKs.C(%d,%d,%d,%d,%d,%d)\'",run, year, isData, mcType, trackType, logFlag));
			}
			else if(block == 1)
			{
				// sWeight data
				if(isData)
				{
					if(run == 2)
					{
						// cout<<"*****Hadding Run 2 sanity files"<<endl;
						gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/"
						            "JpsiLambda/run2/");
						gSystem->Exec("hadd -f jpsilambda_sanity_LL.root jpsilambda_sanity_LL_2015.root jpsilambda_sanity_LL_2016.root "
						              "jpsilambda_sanity_LL_2017.root jpsilambda_sanity_LL_2018.root");
						gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts/");
					}

					//sWeight data after sanity
					cout<<"***sWeight Sanity run"<<run<<endl;
					//		  myChi2_Sanity = DoSWeight_Sanity(run, trackType, logFlag);
					gSystem->Exec(Form("root -l -b -q \'DoSWeight_Sanity.C(%d,%d,%d)\'",run, trackType, logFlag));
					if(run == 2)
					{
						cout<<"*****Hadding Run 2 cutoutks files"<<endl;
						gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/"
						            "JpsiLambda/run2/");
						gSystem->Exec("hadd -f jpsilambda_cutoutks_LL_nonZeroTracks.root jpsilambda_cutoutks_LL_2015_nonZeroTracks.root jpsilambda_cutoutks_LL_2016_nonZeroTracks.root "
						              "jpsilambda_cutoutks_LL_2017_nonZeroTracks.root jpsilambda_cutoutks_LL_2018_nonZeroTracks.root");
						cout<<"*****Done hadding nonZeroTracks files"<<endl;
						gSystem->Exec("hadd -f jpsilambda_cutoutks_LL_ZeroTracks.root jpsilambda_cutoutks_LL_2015_ZeroTracks.root jpsilambda_cutoutks_LL_2016_ZeroTracks.root "
						              "jpsilambda_cutoutks_LL_2017_ZeroTracks.root jpsilambda_cutoutks_LL_2018_ZeroTracks.root");
						cout<<"*****Done hadding ZeroTracks files"<<endl;
						gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts/");
					}

					//sWeight nonZeroTracks data after cutoutks
					cout<<"****DoSWeight run "<<run<<" nonZeroTracks"<<endl;

					//					myChi2_nonZero = DoSWeight(run, trackType, logFlag, false);
					gSystem->Exec(Form("root -l -q -b \'DoSWeight.C(%d,%d,%d,false)\'",run, trackType, logFlag));
					//sWeight ZeroTracks data after cutoutks
					cout<<"****DoSWeight run "<<run<<" ZeroTracks"<<endl;

					//					myChi2_Zero = DoSWeight(run, trackType, logFlag, true);
					gSystem->Exec(Form("root -l -q -b \'DoSWeight.C(%d,%d,%d,true)\'",run, trackType, logFlag));
					// if(myChi2_Sanity > 2.0)
					// {
					//      cout<<"####Fit failed for Sanity sWeighting. Exiting####"<<endl;
					//      exit(1);
					// }
					// if(myChi2_nonZero > 2.0)
					// {
					//      cout<<"####Fit failed for nonZeroTracks sWeighting. Exiting####"<<endl;
					//      exit(1);
					// }
					// if(myChi2_Zero > 2.0)
					// {
					//      cout<<"####Fit failed for ZeroTracks sWeighting. Exiting####"<<endl;
					//      exit(1);
					// }
					cout<<"*****Calculating GB weights files for data/MC correction"<<endl;
					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
					gSystem->Exec(Form("python GB_RW.py %d",run));
				}
			}
			else if(block == 2)
			{
				//****ALL THIS STUFF IN THIS BLOCK IS INDEPENDENT OF ISOLATION. ONLY NEEDS TO BE EXECUTED TWICE. ONCE FOR EACH BDT CONF*******

				if(isData)
				{
					cout<<"*****Hadding nonZero and Zero sPlot files"<<endl;
					gSystem->cd(Form("/data1/avenkate/JpsiLambda_RESTART/rootFiles/dataFiles/"
					                 "JpsiLambda/run%d/",run));
					gSystem->Exec("hadd -f jpsilambda_LL_withsw.root jpsilambda_LL_withsw_nonZeroTracks.root jpsilambda_LL_withsw_ZeroTracks.root ");
					gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts/");

					//Train Final BDT on data sans isolation
					cout<<"***TrainFinalBDT ZeroTracks run "<<run<<"finalBDTconf "<<finalBDTconf<<" ***"<<endl;

					//                                                isoFlag
					//					TrainFinalBDT(run, trackType, isoVersion, isoConf, false, finalBDTconf, logFlag);
					gSystem->Exec(Form("root -l -b -q \'TrainFinalBDT.C(%d,%d,\"%s\",%d,false,%d,%d,%d)\'",run, trackType, isoVersion, isoConf, finalBDTconf, logFlag, simFlag));
				}

				//Apply final BDT on zeroTracks data/MC
				cout<<"***ApplyFinalBDT all ZeroTracks run "<<run<<" finalBDTconf "<<
				        finalBDTconf<<" ***"<<endl;
				//				ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, true, logFlag);
				gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,%d,\"%s\",%d,%d,1,%d,true,%d,%d)\'",run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag, simFlag));
				if(isData)
				{
					//Apply final BDT on ZeroTracks sWeighted data
					cout<<"***ApplyFinalBDT sWeight ZeroTracks run "<<run<<" isoVersion "<<
					        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
					        finalBDTconf<<"***"<<endl;
					//					ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, true, logFlag);
					gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,%d,\"%s\",%d,%d,2,%d,true,%d, %d)\'",run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag, simFlag));
				}
				// ****************************************************************************************************************************
			}
			else if(block == 3)
			{
				if(isData)
				{
					// Train Isolation BDT on nonZeroTracks data. 2 configs train separately,
					cout<<"***TrainIsolation run "<<run<<" isoVersion "<<isoVersion<<" isoConf "<<isoConf<<"***"<<endl;
					//					TrainIsolation(run, trackType, isoVersion, isoConf, logFlag);
					gSystem->Exec(Form("root -l -b -q \'TrainIsolation.C(%d,%d,\"%s\",%d,%d,%d)\'",run, trackType, isoVersion, isoConf, logFlag, simFlag));
				}

				// Apply isolation BDT on nonZeroTracks data/MC
				cout<<"***ApplyIsolation on nonZeroTracks  data/MC run "<<run<<" isoVersion "<<isoVersion<<
				        " isoConf "<<isoConf<<"***"<<endl;
				//				ApplyIsolation(run, isData, mcType, trackType, 1, isoVersion, isoConf, logFlag);
				gSystem->Exec(Form("root -l -q -b \'ApplyIsolation.C(%d,%d,%d,%d,1,\"%s\",%d,%d, %d)\'",run, isData, mcType, trackType, isoVersion, isoConf, logFlag, simFlag));

				if(isData)
				{
					//Apply isolation BDT on nonZeroTracks Weighted data
					cout<<"***ApplyIsolation sWeight data run "<<run<<" isoVersion "<<isoVersion<<
					        " isoConf "<<isoConf<<"***"<<endl;
					//					ApplyIsolation(run, isData, mcType, trackType, 2, isoVersion, isoConf, logFlag);
					gSystem->Exec(Form("root -l -q -b \'ApplyIsolation.C(%d,%d,%d,%d,2,\"%s\",%d,%d,%d)\'",run, isData, mcType, trackType, isoVersion, isoConf, logFlag, simFlag));
				}
			}
			else if(block == 4)
			{
				if(isData && isoFlag)
				{
					//Train Final BDT on data w/ isolation. 2 configs. Train separately.
					cout<<"***TrainFinalBDT nonZeroTracks run "<<run<<" isoVersion "<<isoVersion<<
					        " isoConf "<<isoConf<<" ***"<<endl;

					//                                                isoFlag
					//					TrainFinalBDT(run, trackType, isoVersion, isoConf, true, finalBDTconf, logFlag);
					gSystem->Exec(Form("root -l -b -q \'TrainFinalBDT.C(%d,%d,\"%s\",%d,true,%d,%d,%d)\'",run, trackType, isoVersion, isoConf, finalBDTconf, logFlag, simFlag));
				}

				//Apply final BDT on nonZeroTracks data/MC
				cout<<"***ApplyFinalBDT all nonZeroTracks run "<<run<<"isoVersion "<<
				        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
				        finalBDTconf<<" ***"<<endl;
				//				ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, false, logFlag);
				gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,%d,\"%s\",%d,%d,1,%d,false,%d,%d)\'",run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag, simFlag));
				if(isData)
				{
					// Apply final BDT on nonZeroTracks sWeighted data
					cout<<"***ApplyFinalBDT sWeight nonZeroTracks run "<<run<<" isoVersion "<<
					        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
					        finalBDTconf<<" ***"<<endl;
					//					ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, false, logFlag);
					gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,%d,\"%s\",%d,%d,2,%d,false,%d, %d)\'",run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, isoFlag, logFlag, simFlag));

				}
			}
			else if(block == 5)
			{
				if(isData)
				{
					// Optimize final BDT cut based on some FoM
					cout<<"***OptimizeFinalBDT run "<<run<<" isoVersion "<<
					        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
					        finalBDTconf<<" FOM "<<FOM<<" ***"<<endl;
					gSystem->Exec(Form("root -l -b -q \'OptimizeFinalBDT.C(%d,\"%s\",%d,%d,%d,%d,\"%s\",\"sigma\")\'",run, isoVersion, isoConf,finalBDTconf, isoFlag, logFlag, FOM));
					// bdtCuts = OptimizeFinalBDT(run, isoVersion, isoConf,finalBDTconf, isoFlag, logFlag, FOM,"sigma");
				}
			}
		}
		if(isData) mcType = 1;
	} //end loop on mcTypes

	// bdtCut = bdtCuts[0];
	// bdtCut_ZeroTracks = bdtCuts[1];
	// // Make cut on final BDT at optimum point for both nonZeroTracks and ZeroTracks, and combine the cut files.
	// cout<<"***CutFinalBDT run "<<run<<" isoVersion "<<
	//         isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
	//         finalBDTconf<<" FOM "<<FOM<<" nonZeroCut "<<bdtCut<<" ZeroCut "<<bdtCut_ZeroTracks<<" ***"<<endl;
	// CutFinalBDT(run, isData, mcType, isoVersion, isoConf,
	//             finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
	//             logFlag, FOM, "sigma");
	// CutFinalBDT(run, !isData, mcType, isoVersion, isoConf,
	//             finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
	//             logFlag, FOM, "sigma");
	sw.Stop();
	cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
