/********************************
   Author : Aravindhan V.
   The purpose of this script is to apply the isolation BDT to data/MC.
   INPUT: nonZeroTracks file from cutOutKs
   OUTPUT: root file containing only isolation BDT output.
   An isolation BDT value is calculated for every added track
   in a given event. The minimum of thes values for a given event is What
   is written out.
 *********************************/

#include "ApplyIsolation.h"
void ApplyIsolation(Int_t run, Bool_t isData, Int_t mcType, Int_t trackType,
                    Int_t flag, const char* isoVersion, Int_t isoConf,
                    Bool_t logFlag)
/* >run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC.
   Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data

   >isData = true for data, false for MC

   >mcType = 0 when running over data. When running over MC,
   mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.

   >trackType = 3 for LL, 5 for DD.

   >flag = 1 when applying on all data,
   flag = 2 when applying only on signal training sample

   >isoVersion = "v0","v1"

   >isoConf = 1 or 2.
 */
{
	TStopwatch sw;
	sw.Start();

	gSystem->cd("/data1/avenkate/JpsiLambda_RESTART");
	const char *folder = "", *part = "", *type = "";
	type        = (trackType == 3) ? ("LL") : ("DD");

	switch(mcType)
	{
	case 0:
	{
		folder = "";
		part = "";
		break;
	}
	case 1:
	{
		folder = "JpsiLambda";
		part = "jpsilambda";
		break;
	}
	case 2:
	{
		folder = "JpsiSigma";
		part = "jpsisigma";
		break;
	}
	case 3:
	{
		folder = "JpsiXi";
		part = "jpsixi";
		break;
	}
	case 4:
	{
		folder = "Bu_JpsiX";
		part = "bu_jpsix";
		break;
	}
	case 5:
	{
		folder = "Bd_JpsiX";
		part = "bd_jpsix";
		break;
	}
	}

	if(!isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/mc/JpsiLambda/%s/run%d/ApplyIsolation_%s_%s_conf%d.txt",
		                             folder,run,type,isoVersion,isoConf),"w");
	}
	if(isData && logFlag)
	{
		gSystem->RedirectOutput(Form("logs/data/JpsiLambda/run%d/ApplyIsolation_%s_%s_conf%d.txt",
		                             run,type,isoVersion,isoConf),"a");
	}
	cout<<"******************************************"<<endl;
	cout<< "==> Start of ApplyIsolation"<<endl;
	gSystem->Exec("date");
	cout<<"WD = "<<gSystem->pwd()<<endl;
	cout<<"******************************************"<<endl;

	TFile *fileIn = nullptr, *fileOut = nullptr;
	TTree *treeIn = nullptr, *treeOut = nullptr;

	TString dir = "", prefix = "", methodName = "", weightFile = "";
	const Int_t Max = 200;

	Int_t entries_init     = 0, nTracks        = 0;
	Float_t PT[Max]        = {0.}, IPCHI2[Max] = {0.}, VCHI2DOF[Max] = {0.};
	Float_t MINIPCHI2[Max] = {0.}, FD[Max]     = {0.}, FDCHI2[Max] = {0.};
	Float_t ipChi2         = 0., vChi2Dof      = 0., log_minIpChi2 = 0.;
	Float_t log_PT         = 0., log_FD        = 0., log_fdChi2 = 0.;
	Double_t BDTk[Max]     = {0.}, BDTkMin     = 1.1;

	const char *rootFolder  = "";

	TMVA::Tools::Instance();// This loads the library
	TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );

	//set up input, output, logging
	if(!isData) // MC
	{
		rootFolder = Form("rootFiles/mcFiles/JpsiLambda/%s/run%d",folder,run);

		fileIn  = TFile::Open(Form("%s/%s_cutoutks_%s_nonZeroTracks.root",
		                           rootFolder,part,type));
		fileOut = new TFile(Form("%s/%s_%s_iso%d_%s.root",
		                         rootFolder,part,type,isoConf,isoVersion),"RECREATE");
	}//end MC block
	else // Data
	{
		rootFolder = Form("rootFiles/dataFiles/JpsiLambda/run%d",run);

		if(flag == 1)
		{
			fileIn  = TFile::Open(Form("%s/jpsilambda_cutoutks_%s_nonZeroTracks.root",
			                           rootFolder,type));
			fileOut = new TFile(Form("%s/jpsilambda_%s_iso%d_%s.root",
			                         rootFolder,type,isoConf,isoVersion),"RECREATE");
		}
		else if(flag == 2)
		{
			fileIn  = TFile::Open(Form("%s/jpsilambda_%s_withsw_nonZeroTracks.root",
			                           rootFolder,type));
			fileOut = new TFile(Form("%s/jpsilambda_%ssig_iso%d_%s.root",
			                         rootFolder,type,isoConf,isoVersion),"RECREATE");
		}
	}
	//end Data block
	//end set up of input, output
	if (!fileIn)
	{
		cout << "ERROR: could not open input file" << endl;
		exit(1);
	}
	cout<<"******************************************"<<endl;
	cout<<"Processing Run "<<run<<" "<<type
	    <<((isData) ? (" Data ") : (" MC type "))<<mcType
	    <<((flag == 2) ? (" sWeighted") : (""))<<endl;
	cout<<"******************************************"<<endl;

	cout<<"******************************************"<<endl;
	cout<<"Input file  = "<<fileIn->GetName()<<endl;
	cout<<"Output file = "<<fileOut->GetName()<<endl;
	cout<<"******************************************"<<endl;

	treeIn = (TTree*)fileIn->Get("MyTuple");
	treeOut = new TTree("MyTuple","");

	entries_init = treeIn->GetEntries();
	cout<<"Incoming entries = "<<entries_init<<endl;

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and
	// type to those given in the weight file(s) used

	reader->AddVariable("IPCHI2", &ipChi2);
	reader->AddVariable("VCHI2DOF", &vChi2Dof);
	reader->AddVariable("log_minIpChi2 := log10(MINIPCHI2)", &log_minIpChi2);

	if(strncmp(isoVersion,"v1",2) == 0)
	{
		reader->AddVariable("log_PT := log10(PT)", &log_PT);
	}
	// else if(strncmp(isoVersion,"v2",2) == 0)
	// {
	//      reader->AddVariable("log_FD := log10(FD)", &log_FD);
	//      reader->AddVariable("log_fdChi2 := log10(FDCHI2)", &log_fdChi2);
	// }
	// else if(strncmp(isoVersion,"v3",2) == 0)
	// {
	//      reader->AddVariable("log_PT := log10(PT)", &log_PT);
	//      reader->AddVariable("log_FD := log10(FD)", &log_FD);
	//      reader->AddVariable("log_fdChi2 := log10(FDCHI2)", &log_fdChi2);
	// }

	dir        = "dataset/weights/";
	prefix     = Form("TMVAClassification300-isok%s_dataRun%d_%s",
	                  type,run,isoVersion);
	methodName = "BDT method";

	weightFile = dir + prefix + TString("_") +
	             Form("isoConf%d_300",isoConf) + TString(".weights.xml");
	reader->BookMVA(methodName, weightFile);

	treeIn->SetBranchAddress("Added_n_Particles", &nTracks);
	treeIn->SetBranchAddress("psi_1S_H_IPCHI2_NEW", IPCHI2);
	treeIn->SetBranchAddress("psi_1S_H_VERTEXCHI2_NEW", VCHI2DOF);
	treeIn->SetBranchAddress("psi_1S_H_MINIPCHI2", MINIPCHI2);

	if(strncmp(isoVersion,"v1",2) == 0 || strncmp(isoVersion,"v3",2) == 0)
	{
		treeIn->SetBranchAddress("Added_H_PT", PT);
	}
	// if(strncmp(isoVersion,"v2",2) == 0 || strncmp(isoVersion,"v3",2) == 0)
	// {
	//      treeIn->SetBranchAddress("psi_1S_H_FD_NEW", FD);
	//      treeIn->SetBranchAddress("psi_1S_H_FDCHI2_NEW", FDCHI2);
	// }

	treeOut->Branch(Form("BDTkMin_%s",isoVersion),&BDTkMin,
	                Form("BDTkMin_%s/D",isoVersion));

	cout << "--- Processing: " << entries_init << " events" <<endl;

	if(strncmp(isoVersion,"v0",2) == 0)
	{
		for (Long64_t ievt = 0; ievt < entries_init; ievt++)
		{
			if (ievt%100000 == 0)
			{
				cout << "--- ... Processing event: "
				     << ievt << endl;
			}

			treeIn->GetEntry(ievt);
			BDTkMin = 1.1;

			//Loop over all the added tracks in a given event
			for(Int_t j = 0; j < nTracks; j++)
			{
				ipChi2   = IPCHI2[j];
				vChi2Dof = VCHI2DOF[j];

				if(MINIPCHI2[j] < 0)
				{
					BDTk[j] = 1.1;
				}
				else
				{
					log_minIpChi2 = log10(MINIPCHI2[j]);
					BDTk[j]       = reader->EvaluateMVA(methodName);
				}
				if(BDTk[j] < BDTkMin)
				{
					BDTkMin = BDTk[j];
				}
			}
			treeOut->Fill();
		}//end of v0 for loop
	}
	else if(strncmp(isoVersion,"v1",2) == 0)
	{
		for (Long64_t ievt = 0; ievt < entries_init; ievt++)
		{
			if (ievt%50000 == 0) cout << "--- ... Processing event: " << ievt << endl;

			treeIn->GetEntry(ievt);
			BDTkMin = 1.1;

			//Loop over all the added tracks in a given event
			for(Int_t j = 0; j < nTracks; j++)
			{
				ipChi2   = IPCHI2[j];
				vChi2Dof = VCHI2DOF[j];

				if((MINIPCHI2[j] < 0) || (PT[j] < 0))
				{
					BDTk[j] = 1.1;
				}
				else
				{
					log_minIpChi2 = log10(MINIPCHI2[j]);
					log_PT        = log10(PT[j]);
					BDTk[j]       = reader->EvaluateMVA(methodName);
				}
				if(BDTk[j] < BDTkMin)
				{
					BDTkMin = BDTk[j];
				}
			}
			treeOut->Fill();
		}//end of v1 for loop
	}
	// else if(strncmp(isoVersion,"v2",2) == 0)
	// {
	//      for (Long64_t ievt = 0; ievt < entries_init; ievt++)
	//      {
	//              if (ievt%50000 == 0) cout << "--- ... Processing event: " << ievt << endl;
	//
	//              treeIn->GetEntry(ievt);
	//              BDTkMin = 1.1;
	//
	//              //Loop over all the added tracks in a given event
	//              for(Int_t j = 0; j < nTracks; j++)
	//              {
	//                      ipChi2   = IPCHI2[j];
	//                      vChi2Dof = VCHI2DOF[j];
	//
	//                      if((MINIPCHI2[j] < 0) || (FD[j] < 0) || (FDCHI2[j] < 0))
	//                      {
	//                              BDTk[j] = 1.1;
	//                      }
	//                      else
	//                      {
	//                              log_minIpChi2 = log10(MINIPCHI2[j]);
	//                              log_FD        = log10(FD[j]);
	//                              log_fdChi2    = log10(FDCHI2[j]);
	//                              BDTk[j]       = reader->EvaluateMVA(methodName);
	//                      }
	//                      if(BDTk[j] < BDTkMin)
	//                      {
	//                              BDTkMin = BDTk[j];
	//                      }
	//              }
	//              treeOut->Fill();
	//      }//end of v2 for loop
	// }
	// else if(strncmp(isoVersion,"v3",2) == 0)
	// {
	//      for (Long64_t ievt = 0; ievt < entries_init; ievt++)
	//      {
	//              if (ievt%50000 == 0) cout << "--- ... Processing event: " << ievt << endl;
	//
	//              treeIn->GetEntry(ievt);
	//              BDTkMin = 1.1;
	//
	//              //Loop over all the added tracks in a given event
	//              for(Int_t j = 0; j < nTracks; j++)
	//              {
	//                      ipChi2   = IPCHI2[j];
	//                      vChi2Dof = VCHI2DOF[j];
	//
	//                      if((MINIPCHI2[j] < 0) || (PT[j] < 0) || (FD[j] < 0) || (FDCHI2[j] < 0))
	//                      {
	//                              BDTk[j] = 1.1;
	//                      }
	//                      else
	//                      {
	//                              log_minIpChi2 = log10(MINIPCHI2[j]);
	//                              log_PT        = log10(PT[j]);
	//                              log_FD        = log10(FD[j]);
	//                              log_fdChi2    = log10(FDCHI2[j]);
	//                              BDTk[j]       = reader->EvaluateMVA(methodName);
	//                      }
	//                      if(BDTk[j] < BDTkMin)
	//                      {
	//                              BDTkMin = BDTk[j];
	//                      }
	//              }
	//              treeOut->Fill();
	//      }//end of v3 for loop
	// }

	fileOut->cd();
	treeOut->Write();
	fileOut->Close();

	cout<< "--- Created root file: "<<fileOut->GetName()
	    <<" containing the MVA output histograms" <<endl;

	delete reader;

	sw.Stop();// Get elapsed time
	cout<< "==> End of ApplyIsolation! Isolated forever!: "; sw.Print();

	if(logFlag) gSystem->RedirectOutput(0); //redirect back to stdout,stderr
//	if(logFlag && !isData) gROOT->ProcessLine(".>");
}
