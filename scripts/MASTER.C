#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TStopwatch.h"
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void MASTER(Int_t run = 1)
/*
   mcType = 0 when processing data. non zero for processing MC.
   mcType = 1 for Lb -> J/psi Lambda MC
   mcType = 2 for Lb -> J/psi Sigma MC        (reco'd JpsiLambda)
   mcType = 3 for Xib -> Jpsi Xi MC           (reco'd JpsiLambda)
   mcType = 4 for Lb -> J/psi Lambda*(1405)MC (reco'd JpsiLambda)
   mcType = 5 for Lb -> J/psi Lambda*(1520)MC (reco'd JpsiLambda)
   mcType = 6 for Lb -> J/psi Lambda*(1600)MC (reco'd JpsiLambda)
   mcType = 7 for B0 -> Jpsi Ks MC            (reco'd JpsiLambda)
   mcType = 8 for Xib0 -> J/psi Lambda MC
   mcType = 9 for Xib- -> J/psi Xi- MC        (reco'd J/psi Xi-)
 */
{
	TStopwatch sw;
	sw.Start();

	cout<<"*****"<<endl;
	cout<<"*****"<<endl;
	cout<<"*****"<<endl;
	cout<<"*****"<<endl;
	cout<<"*****"<<endl;
	gSystem->Exec("date"); //print out date and time

	Int_t config = 7; //Double check that this is the optimal config

	Int_t mcType       = 0;
	Int_t isoConf      = 1;// config for isolation BDT. Currently 1 or 2 supported
	Int_t finalBDTconf = 1;// config for final BDT. Currently 1 or 2 supported

	Bool_t logFlag = true; // set to false only while testing.
	Bool_t testing = true;// when true, analysis will only run over a subset of data
	Bool_t isoFlag = true; // when true, isolation will be used in final BDT.
	Bool_t simFlag = false;// set to true to train on simulation

	Double_t myChi2_nonZero = 0.0;
	Double_t myChi2_Zero = 0.0;
	Double_t myChi2_Sanity  = 0.0;

	TString isoV = ""; // which version of isolation BDT? changes input variables used in isolation. v0,1 supported.
	TString FOM  = "Punzi";//Sig for S/sqrt(S+B), Punzi for Punzi FOM with a = 3

	std::vector <Double_t> bdtCuts;

	switch(config)
	{
	case 1:
	{
		isoV   = "v0";
		isoConf      = 1;
		finalBDTconf = 1;
		break;
	}
	case 2:
	{
		isoV   = "v0";
		isoConf      = 1;
		finalBDTconf = 2;
		break;
	}
	case 3:
	{
		isoV   = "v0";
		isoConf      = 2;
		finalBDTconf = 1;
		break;
	}
	case 4:
	{
		isoV   = "v0";
		isoConf      = 2;
		finalBDTconf = 2;
		break;
	}
	case 5:
	{
		isoV   = "v1";
		isoConf      = 1;
		finalBDTconf = 1;
		break;
	}
	case 6:
	{
		isoV   = "v1";
		isoConf      = 1;
		finalBDTconf = 2;
		break;
	}
	case 7:
	{
		isoV   = "v1";
		isoConf      = 2;
		finalBDTconf = 1;
		break;
	}
	case 8:
	{
		isoV   = "v1";
		isoConf      = 2;
		finalBDTconf = 2;
		break;
	}
	}

	cout<<"$$$$$$$$$$$ PROCESSING RUN "<<run<<" $$$$$$$$$$$"<<endl<<endl;

	//**********
	//Process B -> J/psi K data
	cout<<"*** Processing raw data for B+ -> J/psi K+ ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'Process_B2JpsiK.C(%d,%d,%d)\'",run,testing,logFlag));

	cout<<"*** sWeight B+ -> J/psi K+ data ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'DoSWeight_JpsiK.C(%d,%d)\'",run,logFlag));
	//**********

	//**********
	//Process Lb -> J/psi Lambda
	for(Int_t i = 1; i >=0; i--) // Process data first. Then all MC.
	{
		Bool_t isData = i;
		Int_t len = 0; // This is the number of different MC channels that need to be processed.

		if(isData)
		{
			len = 1;
		}
		else
		{
			len = 8;
		}
		for(mcType = 1; mcType <=len; mcType++)
		{
			if(isData) mcType = 0;

			if(isData) cout<<"$$$$$$$$$$$ PROCESSING DATA $$$$$$$$$$$"<<endl<<endl;
			else cout<<"$$$$$$$$$$$ PROCESSING MC TYPE "<<mcType<< " $$$$$$$$$$$$$$"<<endl<<endl;

			if(!isData) //MC
			{
				cout<<"*** Collating MC PIDGEN files ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -b -q \'CollateFiles.C(%d,1969,%d,%d)\'",run, isData, mcType)); //year argument is irrelevant for MC

				//Run script to make aliases of branch names in generated MC
				cout<<"*** Aliasing branch names ***"<<endl;
				gSystem->Exec(Form("root -l -b -q \'AliasFile.C(%d,%d,%d)\'",run,mcType,logFlag));

				//Reweight Lambda*(1405) simulation to match theory models
				if(mcType == 4)
				{
					cout<<"*** Reweighting Lambda*(1405) MC ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -b -q \'reweight_Lst1405.C(%d,%d)\'",run,logFlag));
				}

				//Apply lifetime Reweighting (only to MC with Lb as parent)
				if(mcType == 1 || mcType == 2 || mcType == 4 || mcType == 5 || mcType == 6)
				{
					cout<<"*** Applying Lifetime RW to reco MC ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -b -q \'ApplyTauWeight.C(%d,%d,0,%d)\'",run,mcType,logFlag));

					cout<<"*** Applying Lifetime RW to generated MC ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -b -q \'ApplyTauWeight.C(%d,%d,1,%d)\'",run,mcType,logFlag));
				}

				if(mcType == 1)
				{
					cout<<"*** Calculating GB weights files for data/MC correction. Data should have already been processed for this to work. ***"<<endl<<endl;
					gSystem->Exec(Form("python GB_RW.py %d",run));
				}

				//Apply Data-MC correction Reweighting (only to Lb -> J/psi Lamda, Lb -> J/psi Sigma, Xib- -> J/psi Xi-, Xib0 -> J/psi Lambda)
				if(mcType == 1 || mcType == 2 || mcType == 3 || mcType == 8)
				{
					cout<<"*** Applying GB weights on reco MC ***"<<endl<<endl;
					gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts");
					gSystem->Exec(Form("python ApplyGBWeights.py %d %d 0", run, mcType));

					cout<<"*** Applying GB weights on generated MC ***"<<endl<<endl;
					gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts");
					gSystem->Exec(Form("python ApplyGBWeights.py %d %d 1", run, mcType));
				}

				//Trigger Cut
				cout<<"*** Trigger ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -b -q \'Trigger.C(%d,1969,%d,%d,%d,%d)\'",run, isData, mcType, testing, logFlag));

				//Sanity Cuts
				cout<<"*** Sanity ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -q -b \'Sanity.C(%d,1969,%d,%d,%d)\'",run, isData, mcType, logFlag));

				//Cut Out Ks0 - PID Veto
				cout<<"*** CutOutKs ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -q -b \'CutOutKs.C(%d,1969,%d,%d,%d)\'",run, isData, mcType, logFlag));
			}
			else if(isData)
			{
				if(run == 1)
				{
					//Trigger Cut
					cout<<"*** Trigger 2011 ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -b -q \'Trigger.C(%d,2011,%d,%d,%d,%d)\'",run, isData, mcType, testing, logFlag));

					cout<<"*** Trigger 2012 ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -b -q \'Trigger.C(%d,2012,%d,%d,%d,%d)\'",run, isData, mcType, testing, logFlag));

					cout<<"*** Hadding Run 1 Trigger files ***"<<endl<<endl;
					gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/rootFiles/dataFiles/"
					            "JpsiLambda/run1/");

					gSystem->Exec("hadd -f jpsilambda_triggered.root jpsilambda_triggered_20{11,12}.root");
					gSystem->Exec("rm -f jpsilambda_triggered_20{11,12}.root");

					gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts/");

					//Sanity Cuts
					cout<<"*** Sanity ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -q -b \'Sanity.C(%d,1969,%d,%d,%d)\'",run, isData, mcType, logFlag));

					//Cut Out Ks0 - PID Veto
					cout<<"*** CutOutKs ***"<<endl<<endl;
					gSystem->Exec(Form("root -l -q -b \'CutOutKs.C(%d,1969,%d,%d,%d)\'",run, isData, mcType, logFlag));
				}
				else if(run == 2)
				{
					for(Int_t year = 2015; year<=2018; year++)
					{
						//Trigger Cut
						cout<<"*** Trigger "<<year<<" ***"<<endl<<endl;
						gSystem->Exec(Form("root -l -b -q \'Trigger.C(%d,%d,%d,%d,%d,%d)\'",run, year, isData, mcType, testing, logFlag));

						//Sanity Cuts
						cout<<"*** Sanity "<<year<<" ***"<<endl<<endl;
						gSystem->Exec(Form("root -l -q -b \'Sanity.C(%d,%d,%d,%d,%d)\'",run, year, isData, mcType, logFlag));

						//Cut Out Ks0 - PID Veto
						cout<<"*** CutOutKs "<<year<<" ***"<<endl<<endl;
						gSystem->Exec(Form("root -l -q -b \'CutOutKs.C(%d,%d,%d,%d,%d)\'",run, year, isData, mcType, logFlag));
					}

					cout<<"*** Hadding Run 2 Sanity files ***"<<endl<<endl;

					gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/rootFiles/dataFiles/JpsiLambda/run2/");
					gSystem->Exec("hadd -f jpsilambda_sanity_LL.root jpsilambda_sanity_LL_20{15,16,17,18}.root");
					gSystem->Exec("rm -f jpsilambda_sanity_LL_20{15,16,17,18}.root");

					cout<<"*** Done hadding Sanity files ***"<<endl<<endl;
					cout<<"*** Hadding Run 2 CutOutKs files ***"<<endl<<endl;

					gSystem->Exec("hadd -f jpsilambda_cutoutks_LL_nonZeroTracks.root jpsilambda_cutoutks_LL_20{15,16,17,18}_nonZeroTracks.root");
					gSystem->Exec("rm -f jpsilambda_cutoutks_LL_20{15,16,17,18}_nonZeroTracks.root");
					cout<<"*** Done hadding nonZeroTracks files ***"<<endl<<endl;

					gSystem->Exec("hadd -f jpsilambda_cutoutks_LL_ZeroTracks.root jpsilambda_cutoutks_LL_20{15,16,17,18}_ZeroTracks.root");
					gSystem->Exec("rm -f jpsilambda_cutoutks_LL_20{15,16,17,18}_ZeroTracks.root");
					cout<<"*** Done hadding ZeroTracks files ***"<<endl<<endl;

					gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts/");
				}

				//sWeight data after sanity
				cout<<"*** sWeight Sanity ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -b -q \'DoSWeight_Sanity.C(%d,%d)\'",run, logFlag));

				//sWeight nonZeroTracks data after cutoutks
				cout<<"*** sWeight nonZeroTracks ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -q -b \'DoSWeight.C(%d,%d,false)\'",run, logFlag));

				//sWeight ZeroTracks data after cutoutks
				cout<<"*** sWeight ZeroTracks ***"<<endl<<endl;
				gSystem->Exec(Form("root -l -q -b \'DoSWeight.C(%d,%d,true)\'",run, logFlag));

				ifstream infile1(Form("chi2_sanity.txt"));
				ifstream infile2(Form("chi2_nonZeroTracks.txt"));
				ifstream infile3(Form("chi2_ZeroTracks.txt"));

				infile1>>myChi2_Sanity;
				infile2>>myChi2_nonZero;
				infile3>>myChi2_Zero;

				if(myChi2_Sanity > 2.0)
				{
					cout<<"#### Fit failed for Sanity sWeighting. Exiting ####"<<endl<<endl;
					exit(1);
				}
				if(myChi2_nonZero > 2.0)
				{
					cout<<"#### Fit failed for nonZeroTracks sWeighting. Exiting ####"<<endl<<endl;
					exit(1);
				}
				if(myChi2_Zero > 2.0)
				{
					cout<<"#### Fit failed for ZeroTracks sWeighting. Exiting ####"<<endl<<endl;
					exit(1);
				}
			}

			{ // Code block for training and applying noIso finalBDT on zeroTracks data/MC
				if((isData && !simFlag) || (!isData && simFlag && mcType == 1))
				{
					if(isData)
					{
						cout<<"*** Hadding nonZero and Zero sPlot files ***"<<endl<<endl;

						gSystem->cd(Form("/data1/avenkate/JpsiLambda_TESTING/rootFiles/dataFiles/JpsiLambda/run%d/",run));
						gSystem->Exec("hadd -f jpsilambda_LL_withsw.root jpsilambda_LL_withsw_nonZeroTracks.root jpsilambda_LL_withsw_ZeroTracks.root ");

						cout<<"*** Done hadding nonZero and Zero sPlot files ***"<<endl<<endl;

						gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts/");
					}

					//Train noIso Final BDT on data/MC sans isolation
					cout<<"*** TrainFinalBDT ZeroTracks finalBDTconf "<<finalBDTconf<<" ***"<<endl<<endl;
					Bool_t isoFlag_temp = false;
					gSystem->Exec(Form("root -l -b -q \'TrainFinalBDT.C(%d,\"%s\",%d,%d,%d,%d,%d)\'",
					                   run, isoV.Data(), isoConf, isoFlag_temp, finalBDTconf, logFlag, simFlag));
				}

				//Apply noIso final BDT on zeroTracks data/MC
				cout<<"*** ApplyFinalBDT all ZeroTracks finalBDTconf "<<finalBDTconf<<" ***"<<endl<<endl;

				Int_t flag = 1;
				Bool_t zeroFlag = true; //these 2 flags make it apply on all zeroTracks data/MC

				gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,\"%s\",%d,%d,%d,%d,%d,%d,%d)\'",
				                   run, isData, mcType, isoV.Data(), isoConf, finalBDTconf, flag, isoFlag, zeroFlag, logFlag, simFlag));
				// if(isData)
				// {
				//      //Apply final BDT on ZeroTracks sWeighted data
				//      cout<<"*** ApplyFinalBDT sWeight ZeroTracks isoV "<<isoV<<" isoConf "<<isoConf<<" finalBDTconf "<<finalBDTconf<<" ***"<<endl<<endl;
				//      flag = 2; // this makes it apply only on the sWeighted data sample
				//
				//      gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,\"%s\",%d,%d,%d,%d,%d,%d,%d)\'",
				//                         run, isData, mcType, isoV.Data(), isoConf, finalBDTconf, flag, isoFlag, zeroFlag, logFlag, simFlag));
				// }
			}

			{//Code block for training and applying isolation BDT on nonZeroTracks data.
				if((isData && !simFlag) || (!isData && simFlag && mcType == 1))
				{
					// Train Isolation BDT on nonZeroTracks data/MC. 2 configs train separately,
					cout<<"*** TrainIsolation run "<<run<<" isoV "<<isoV<<" isoConf "<<isoConf<<" ***"<<endl<<endl;

					gSystem->Exec(Form("root -l -b -q \'TrainIsolation.C(%d,\"%s\",%d,%d,%d)\'",
					                   run, isoV.Data(), isoConf, logFlag, simFlag));
				}

				Int_t flag = 1; //apply on all nonZeroTracks data/MC
				// Apply isolation BDT on nonZeroTracks data/MC
				cout<<"*** ApplyIsolation on nonZeroTracks data/MC isoV "<<isoV<<
				        " isoConf "<<isoConf<<"***"<<endl<<endl;

				gSystem->Exec(Form("root -l -q -b \'ApplyIsolation.C(%d,%d,%d,%d,\"%s\",%d,%d, %d)\'",
				                   run, isData, mcType, flag, isoV.Data(), isoConf, logFlag, simFlag));

				if(isData)
				{
					//Apply isolation BDT on nonZeroTracks sWeighted data
					cout<<"*** ApplyIsolation sWeight data isoV "<<isoV<<" isoConf "<<isoConf<<" ***"<<endl<<endl;
					flag = 2; // apply on only sWeighted nonZeroTracks data

					gSystem->Exec(Form("root -l -q -b \'ApplyIsolation.C(%d,%d,%d,%d,\"%s\",%d,%d,%d)\'",
					                   run, isData, mcType, flag, isoV.Data(), isoConf, logFlag, simFlag));
				}
			}

			{// Code block for training and applying finalBDT
				if( isoFlag && ((isData && !simFlag) || (!isData && simFlag && mcType == 1)) )
				{
					//Train Final BDT on data w/ isolation. 2 configs. Train separately.
					cout<<"*** TrainFinalBDT nonZeroTracks isoV "<<isoV<<" isoConf "<<isoConf<<" ***"<<endl<<endl;

					gSystem->Exec(Form("root -l -b -q \'TrainFinalBDT.C(%d,\"%s\",%d,%d,%d,%d,%d)\'",
					                   run, isoV.Data(), isoConf, isoFlag, finalBDTconf, logFlag, simFlag));
				}

				Int_t flag = 1;
				Bool_t zeroFlag = false;

				//Apply final BDT on nonZeroTracks data/MC (if isoFlag is true. otherwise will apply noIso BDT on ALL data)
				cout<<"*** ApplyFinalBDT all nonZeroTracks isoV "<<isoV<<" isoConf "<<
				        isoConf<<" finalBDTconf "<<finalBDTconf<<" ***"<<endl<<endl;

				gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,\"%s\",%d,%d,%d,%d,%d,%d,%d)\'",
				                   run, isData, mcType, isoV.Data(), isoConf, finalBDTconf, flag, isoFlag, zeroFlag, logFlag, simFlag));
				// if(isData)
				// {
				//      flag = 2;
				//      // Apply final BDT on nonZeroTracks sWeighted data (if isoFlag is true. otherwise will apply noIso BDT on ALL sWeighted data)
				//      cout<<"*** ApplyFinalBDT sWeight nonZeroTracks isoV "<<isoV<<" isoConf "<<
				//              isoConf<<" finalBDTconf "<<finalBDTconf<<" ***"<<endl<<endl;
				//
				//      gSystem->Exec(Form("root -l -b -q \'ApplyFinalBDT.C(%d,%d,%d,\"%s\",%d,%d,%d,%d,%d,%d,%d)\'",
				//                         run, isData, mcType, isoV.Data(), isoConf, finalBDTconf, flag, isoFlag, zeroFlag, logFlag, simFlag));
				// }
			}

			if(isData) mcType = 1;
		} //end loop on mcTypes
	} // End loop processing J/psi Lambda data + MC.


	{// Code block for optimizing the finalBDT cut. MC needs to be processed for this. Consider taking this out of the loop.
		// Optimize final BDT cut based on some FoM
		cout<<"***OptimizeFinalBDT run "<<run<<" isoV "<<isoV<<" isoConf "<<
		        isoConf<<" finalBDTconf "<<finalBDTconf<<" FOM "<<FOM<<" ***"<<endl<<endl;

		gSystem->Exec(Form("root -l -b -q \'OptimizeFinalBDT.C(%d,\"%s\",%d,%d,%d,%d,\"%s\",\"Sigma\",%d)\'",
		                   run, isoV.Data(), isoConf,finalBDTconf, isoFlag, logFlag, FOM.Data(), simFlag));
	}

	//**********
	//Process Xib- -> J/psi Xi- MC
	cout<<"*** Processing MC for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;

	cout<<"*** Collate MC for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'CollateFiles_JpsiXi.C(%d,1969,false)\'",run));

	cout<<"*** Aliasing branch names ***"<<endl;
	gSystem->Exec(Form("root -l -b -q \'AliasFile.C(%d,9,%d)\'",run,logFlag));

	cout<<"*** Triggering MC for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'Trigger_JpsiXi.C(%d,2011,false,false,%d)\'",run,logFlag));

	cout<<"*** Applying Cuts on MC for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'Cuts_JpsiXi.C(%d,2011,false,%d)\'",run,logFlag));

	cout<<"*** Applying Tracking Correction on MC for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'ApplyTrackingCorr.C(%d,%d)\'",run,logFlag));

	cout<<"*** Applying Data/MC correction to MC for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	cout<<"*** Applying GB weights on reco MC ***"<<endl<<endl;
	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts");
	gSystem->Exec(Form("python ApplyGBWeights.py %d 9 0", run));

	cout<<"*** Applying GB weights on generated MC ***"<<endl<<endl;
	gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts");
	gSystem->Exec(Form("python ApplyGBWeights.py %d 9 1", run));

	cout<<"*** Fit to MC after all cuts for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'Fit_JpsiXi.C(%d,false,%d)\'",run,logFlag));
	//**********

	//**********
	//Process Xib- -> J/psi Xi- data
	cout<<"*** Processing Data for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;

	if(run == 1)
	{
		for(Int_t year = 2011; year <= 2012; year++)
		{
			cout<<"*** Collating+Triggering Data for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
			gSystem->Exec(Form("root -l -b -q \'Trigger_JpsiXi.C(%d,%d,true,false,%d)\'",run,year,logFlag));

			cout<<"*** Applying Cuts on Data for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
			gSystem->Exec(Form("root -l -b -q \'Cuts_JpsiXi.C(%d,%d,false,%d)\'",run,year,logFlag));
		}

		cout<<"*** Hadding 2011+2012 data files for Run 1 fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
		gSystem->cd(Form("/data1/avenkate/JpsiLambda_TESTING/rootFiles/dataFiles/JpsiXi/run%d/",run));
		gSystem->Exec("hadd -f jpsixi_cut_Run1.root jpsixi_cut_LL_20{11,12}.root");
		gSystem->Exec("rm -f jpsixi_cut_LL_20{11,12}.root");
		gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts/");
	}
	else if(run == 2)
	{
		for(Int_t year = 2015; year <= 2018; year++)
		{
			cout<<"*** Collating+Triggering Data for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
			gSystem->Exec(Form("root -l -b -q \'Trigger_JpsiXi.C(%d,%d,true,false,true,%d)\'",run,year,logFlag));

			cout<<"*** Applying Cuts on Data for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
			gSystem->Exec(Form("root -l -b -q \'Cuts_JpsiXi.C(%d,%d,false,%d)\'",run,year,logFlag));
		}

		cout<<"*** Hadding 2015+2016+2017+2018 data files for Run 2 fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
		gSystem->cd(Form("/data1/avenkate/JpsiLambda_TESTING/rootFiles/dataFiles/JpsiXi/run%d/",run));
		gSystem->Exec("hadd -f jpsixi_cut_Run2.root jpsixi_cut_LL_20{15,16,17,18}.root");
		gSystem->Exec("rm -f jpsixi_cut_LL_20{15,16,17,18}.root");
		gSystem->cd("/data1/avenkate/JpsiLambda_TESTING/scripts/");
	}

	cout<<"*** Fit to Data after all cuts for fully reco'd Xib- -> J/psi Xi- ***"<<endl<<endl;
	gSystem->Exec(Form("root -l -b -q \'Fit_JpsiXi.C(%d,false,%d)\'",run,logFlag));

	//**********

	//**********
	//Calculate Xib0 -> J/psi Lambda BF
	cout<<"*** Calculating B(Xib0 -> J/psi Lambda)/B(Xib0 -> J/psi Xi0) ***"<<endl;
	gSystem->Exec(Form("root -l -b -q \'Xib_JpsiLambda.C(%d)\'",logFlag));
	//**********

	// bdtCut = bdtCuts[0];
	// bdtCut_ZeroTracks = bdtCuts[1];
	// // Make cut on final BDT at optimum point for both nonZeroTracks and ZeroTracks, and combine the cut files.
	// cout<<"***CutFinalBDT run "<<run<<" isoV "<<
	//         isoV<<" isoConf "<<isoConf<<" finalBDTconf "<<
	//         finalBDTconf<<" FOM "<<FOM<<" nonZeroCut "<<bdtCut<<" ZeroCut "<<bdtCut_ZeroTracks<<" ***"<<endl<<endl;
	// CutFinalBDT(run, isData, mcType, isoV, isoConf,
	//             finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
	//             logFlag, FOM, "sigma");
	// CutFinalBDT(run, !isData, mcType, isoV, isoConf,
	//             finalBDTconf, bdtCut, bdtCut_ZeroTracks, isoFlag,
	//             logFlag, FOM, "sigma");
	sw.Stop();
	cout << "==> MASTER is done! Huzzah!: "; sw.Print();
}
