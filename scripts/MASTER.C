#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include "DoSWeight.h"
#include "TrainIsolation.h"
#include "ApplyIsolation.h"
#include "TrainFinalBDT.h"
#include "ApplyFinalBDT.h"
#include "OptimizeFinalBDT.h"
#include "CutFinalBDT.h"
#include <iostream>
#include <vector>

using namespace std;

void MASTER()
{
	TStopwatch sw;
	sw.Start();

	gSystem->Exec("date");
	// gROOT->ProcessLine(".L DoSWeight.C+");
	// gROOT->ProcessLine(".L TrainIsolation.C+");
	// gROOT->ProcessLine(".L ApplyIsolation.C+");
	gROOT->ProcessLine(".L TrainFinalBDT.C+");
	gROOT->ProcessLine(".L ApplyFinalBDT.C+");
	gROOT->ProcessLine(".L OptimizeFinalBDT.C+");
	gROOT->ProcessLine(".L CutFinalBDT.C+");

	// gSystem->Load("DoSWeight_C.so");
	// gSystem->Load("TrainIsolation_C.so");
	// gSystem->Load("ApplyIsolation_C.so");
	gSystem->Load("TrainFinalBDT_C.so");
	gSystem->Load("ApplyFinalBDT_C.so");
	gSystem->Load("OptimizeFinalBDT_C.so");
	gSystem->Load("CutFinalBDT_C.so");

	cout<<"poop"<<endl;

	Int_t run              = 1;// Run 1 or Run 2?
	Int_t mcType           = 0;// 0 when running over data. 1, 2, 3 corresponding to J/psi Lambda MC reco'd as J/psi Lambda, J/psi Sigma and J/psi Xi respectively
	Int_t trackType        = 3;// 3 for LL, 5 for DD.
	Int_t isoConf          = 1;// config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf     = 1;// config for final BDT. Currently 1 or 2 supported
	// Int_t massLow          = 5400;
	// Int_t massHigh         = 5700;

	Bool_t testing         = false;// when true, analysis will only run over a subset of data
	Bool_t loose           = true;// when true, analysis will run over data/MC from "loose" stripping line. Only LL
	Bool_t isData          = true;// Data or MC?
	Bool_t isoFlag         = true;// when true, isolation will be used in final BDT.
	Bool_t logFlag         = true;// set to false only while testing.

	Float_t bdtCut         = 0.0;
	Float_t bdtCut_ZeroTracks = 0.0;

	Double_t myChi2_nonZero = 0.0, myChi2_Zero = 0.0;
	const char* isoVersion = ""; // which version of isolation BDT? changes input variables used in isolation. v0,1 supported.

	TString FOM            = "Punzi";//Sig for S/sqrt(S+B), Punzi for Punzi FOM with a = 3
	// TString inFileName  = (trackType == 3) ? ("jpsilambda_LL_forIsoTraining.root") : ("jpsilambda_DD_forIsoTraining.root");
	// TString outFileName = (trackType == 3) ? ("jpsilambda_LL_forIsoTraining_massCut.root") : ("jpsilambda_DD_forIsoTraining_massCut.root");
	Int_t runArray[2]          = {1,2};
	TString isoVersionArray[2] = {"v0","v1"};
	Int_t isoConfArray[2]      = {1,2};
	Int_t finalBDTconfArray[2] = {1,2};

	Int_t mcTypeArray[5]       = {1,2,3,4,5};

	std::vector <Double_t> bdtCuts;

	// for(Int_t h = 0; h <= 4; h++)
	// {
	//      mcType = mcTypeArray[h];
	//      cout<<"$$$$$$$$$$$ Processing MC Type "<<mcType<< " $$$$$$$$$$$$$$"<<endl;
	for(Int_t i = 0; i<=0; i++)
	{
		run = runArray[i];
		//Data- Run 1
		//Trigger Cut
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L Trigger.C+");
		// gROOT->ProcessLine(Form("Trigger(%d, %d, %d, %d, %d, %d)",
		//                         run, isData, mcType, testing, loose, logFlag));
		// gROOT->Reset();

		// //Sanity Cuts
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L Sanity.C+");
		// gROOT->ProcessLine(Form("Sanity(%d, %d, %d, %d)", run, isData, mcType, logFlag));
		// gROOT->Reset();
		//
		//Cut Out Ks0
		// gSystem->cd("/data1/avenkate/JpsiLambda_RESTART/scripts");
		// gROOT->ProcessLine(".L CutOutKs.C+");
		// gROOT->ProcessLine(Form("CutOutKs(%d, %d, %d, %d, %d)",
		//                         run, isData, mcType, trackType, logFlag));
		// gROOT->Reset();


		//sWeight data
		// if(isData)
		// {
		//      //sWeight nonZeroTracks data
		//      cout<<"****DoSWeight run "<<run<<" nonZeroTracks"<<endl;
		//
		//      myChi2_nonZero = DoSWeight(run, trackType, logFlag, false);
		//
		//      if(myChi2_nonZero > 1.5)
		//      {
		//              cout<<"####Fit failed for nonZeroTracks sWeighting. Exiting####"<<endl;
		//              exit(1);
		//      }
		//      //sWeight ZeroTracks data
		//      cout<<"****DoSWeight run "<<run<<" ZeroTracks"<<endl;
		//
		//      myChi2_Zero = DoSWeight(run, trackType, logFlag, true);
		//
		//      if(myChi2_Zero > 1.5)
		//      {
		//              cout<<"####Fit failed for ZeroTracks sWeighting. Exiting####"<<endl;
		//              exit(1);
		//      }
		// }
		for (Int_t j = 0; j <= 1; j++)//loop over isolation BDT versions.
		{
			isoVersion = (isoVersionArray[j]);

			// if(isData)
			// {
			//      // Train Isolation BDT on data.Both configs get trained here.
			//      cout<<"***TrainIsolation run "<<run<<" isoVersion "<<isoVersion<<"***"<<endl;
			//
			//      TrainIsolation(run, trackType, isoVersion, logFlag);
			// }
			for (Int_t k = 0; k <= 1; k++)//loop over isolation BDT configurations
			{
				isoConf = isoConfArray[k];

				// Apply isolation BDT on all data/MC
				// cout<<"***ApplyIsolation all data/MC run "<<run<<" isoVersion "<<isoVersion<<
				//         " isoConf "<<isoConf<<"***"<<endl;
				// ApplyIsolation(run, isData, mcType, trackType, 1, isoVersion, isoConf, logFlag);
				//
				// if(isData)
				// {
				//      //Apply isolation BDT on sWeighted data
				//      cout<<"***ApplyIsolation sWeight data run "<<run<<" isoVersion "<<isoVersion<<
				//              " isoConf "<<isoConf<<"***"<<endl;
				//      ApplyIsolation(run, isData, mcType, trackType, 2, isoVersion, isoConf, logFlag);
				// }
				for(Int_t l = 0; l <= 0; l++)//loop over new Flag versions
				{

					if(isData && isoFlag)
					{
						//Train Final BDT on data w/ isolation
						cout<<"***TrainFinalBDT nonZeroTracks run "<<run<<" isoVersion "<<isoVersion<<
						        " isoConf "<<isoConf<<" ***"<<endl;

						//                                                 isoFlag
						TrainFinalBDT(run, trackType, isoVersion, isoConf, true, logFlag);

						if(j == 0 && k == 0)
						{
							//Train Final BDT on data sans isolation
							cout<<"***TrainFinalBDT ZeroTracks run "<<run<<" ***"<<endl;

							//                                                 isoFlag
							TrainFinalBDT(run, trackType, isoVersion, isoConf, false, logFlag);
						}
					}
					for(Int_t m = 0; m <= 1; m++) //Loop over final BDT configurations.
					{
						finalBDTconf = finalBDTconfArray[m];

						//Apply final BDT on nonZeroTracks data/MC
						cout<<"***ApplyFinalBDT all nonZeroTracks run "<<run<<"isoVersion "<<
						        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
						        finalBDTconf<<" ***"<<endl;
						ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, false, logFlag);

						if(j == 0 && k == 0)
						{
							//Apply final BDT on zeroTracks data/MC
							cout<<"***ApplyFinalBDT all ZeroTracks run "<<run<<" finalBDTconf "<<
							        finalBDTconf<<" ***"<<endl;
							ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 1, isoFlag, true, logFlag);
						}

						if(isData)
						{
							//Apply final BDT on nonZeroTracks sWeighted data/MC
							cout<<"***ApplyFinalBDT sWeight nonZeroTracks run "<<run<<" isoVersion "<<
							        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
							        finalBDTconf<<" ***"<<endl;
							ApplyFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf, finalBDTconf, 2, isoFlag, false, logFlag);

							if(j == 0 && k == 0)
							{
								//Apply final BDT on ZeroTracks sWeighted data/MC
								cout<<"***ApplyFinalBDT sWeight ZeroTracks run "<<run<<" isoVersion "<<
								        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
								        finalBDTconf<<"***"<<endl;
							}
							//Optimize final BDT cut based on some FoM
							cout<<"***OptimizeFinalBDT run "<<run<<" isoVersion "<<
							        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
							        finalBDTconf<<" FOM "<<FOM<<" ***"<<endl;

							bdtCuts = OptimizeFinalBDT(run, trackType, isoVersion, isoConf,
							                           finalBDTconf, isoFlag, logFlag, FOM,
							                           "sigma");
						}

						//		      bdtCut = cutArray[i][j][k][l][m];
						//		      bdtCut_ZeroTracks = zeroCutArray[i][l][m];
						bdtCut = bdtCuts[0];
						bdtCut_ZeroTracks = bdtCuts[1];
						// Make cut on final BDT at optimum point for both nonZeroTracks and ZeroTracks, and combine the cut files.
						cout<<"***CutFinalBDT run "<<run<<" isoVersion "<<
						        isoVersion<<" isoConf "<<isoConf<<" finalBDTconf "<<
						        finalBDTconf<<" FOM "<<FOM<<" nonZeroCut "<<bdtCut<<" ZeroCut "<<bdtCut_ZeroTracks<<" ***"<<endl;
						CutFinalBDT(run, isData, mcType, trackType, isoVersion, isoConf,
						            finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
						            logFlag, FOM, "sigma");
					}//end loop on finalBDTconf
				}//end loop on newFlag
			}//end loop on isoConf
		}//end loop on isoVersion
	}//end loop on runs
	// }//end loop on mcTypes
	sw.Stop();
	cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
