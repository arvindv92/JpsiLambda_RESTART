#include "TFile.h"
#include "iostream"

using namespace std;
void Systematics()
{
	Float_t N_JpsiLambda_Run1           = 4126.77;
	Float_t N_JpsiLambda_StatErr_Run1   = 64.9574;
	Float_t N_JpsiLambda_RelSyst_Run1   = 0;//0.9; // relative % unc. comes from choice of fit model
	Float_t N_JpsiLambda_SystErr_Run1   = N_JpsiLambda_Run1 * (N_JpsiLambda_RelSyst_Run1/100);

	Float_t N_JpsiLambda_Run2           = 15541.2;
	Float_t N_JpsiLambda_StatErr_Run2   = 126.413;
	Float_t N_JpsiLambda_RelSyst_Run2   = 0;//0.8; // relative % unc. comes from choice of fit model
	Float_t N_JpsiLambda_SystErr_Run2   = N_JpsiLambda_Run2 * (N_JpsiLambda_RelSyst_Run2/100);

	Float_t window_JpsiLambda_Run1 = 0.0253382;
	Float_t window_JpsiSigma_Run1  = 0.96434;
	Float_t window_JpsiXi_Run1     = 0.64327;

	Float_t window_JpsiLambda_Run2 = 0.0217437;
	Float_t window_JpsiSigma_Run2  = 0.979967;
	Float_t window_JpsiXi_Run2     = 0.620282;

	Float_t N_JpsiLambda_Run1_window = N_JpsiLambda_Run1*window_JpsiLambda_Run1;
	Float_t N_JpsiLambda_Run2_window = N_JpsiLambda_Run2*window_JpsiLambda_Run2;

	Float_t N_Xib_data_Run1         = 36.36;
	Float_t N_Xib_data_StatErr_Run1 = 7.75666;
	Float_t N_Xib_data_RelSyst_Run1 = 9.0; // relative % unc. comes from choice of fit model
	Float_t N_Xib_data_SystErr_Run1   = N_Xib_data_Run1 * (N_Xib_data_RelSyst_Run1/100);

	Float_t N_Xib_data_Run2         = 173.281;
	Float_t N_Xib_data_StatErr_Run2 = 16.4207;
	Float_t N_Xib_data_RelSyst_Run2 = 6.0; // relative % unc. comes from choice of fit model
	Float_t N_Xib_data_SystErr_Run2   = N_Xib_data_Run2 * (N_Xib_data_RelSyst_Run2/100);

	Float_t eff_Xib_JpsiLambda_Run1         = 0.0190; //% Weighted.
	Float_t eff_Xib_JpsiLambda_SystErr_Run1 = 0.0009; //%

	Float_t eff_Xib_JpsiXi_Run1         = 0.0683; //% Weighted.
	Float_t eff_Xib_JpsiXi_SystErr_Run1 = 0.0017; //%

	Float_t eff_Xib_JpsiLambda_Run2         = 0.0170; //% Weighted.
	Float_t eff_Xib_JpsiLambda_SystErr_Run2 = 0.0007; //%

	Float_t eff_Xib_JpsiXi_Run2         = 0.0681; //% Weighted.
	Float_t eff_Xib_JpsiXi_SystErr_Run2 = 0.0015; //%

	Float_t trackingErr = 0.05; // rel error. includes 4.5% tracking unc, and 2% uncertainty for material effects from hadronic interactions
	Float_t xiVtxUnc = 0.014; // rel error. Unc from vertexing the Xi-. Comes from Steves ANA

	Float_t XibNorm_Run1         = 2 * N_Xib_data_Run1 * (eff_Xib_JpsiLambda_Run1/eff_Xib_JpsiXi_Run1) * window_JpsiXi_Run1;
	Float_t XibNorm_StatErr_Run1 = 2 * N_Xib_data_StatErr_Run1 * (eff_Xib_JpsiLambda_Run1/eff_Xib_JpsiXi_Run1) * window_JpsiXi_Run1;
	Float_t XibNorm_SystErr_Run1 = XibNorm_Run1 * sqrt( pow(N_Xib_data_SystErr_Run1/N_Xib_data_Run1,2) +
	                                                    pow(eff_Xib_JpsiLambda_SystErr_Run1/eff_Xib_JpsiLambda_Run1,2) +
	                                                    pow(eff_Xib_JpsiXi_SystErr_Run1/eff_Xib_JpsiXi_Run1,2) +
	                                                    pow(trackingErr,2) + pow(xiVtxUnc,2));

	cout<<XibNorm_SystErr_Run1<<endl;

	Float_t XibNorm_Run2         = 2 * N_Xib_data_Run2 * (eff_Xib_JpsiLambda_Run2/eff_Xib_JpsiXi_Run2) * window_JpsiXi_Run2;
	Float_t XibNorm_StatErr_Run2 = 2 * N_Xib_data_StatErr_Run2 * (eff_Xib_JpsiLambda_Run2/eff_Xib_JpsiXi_Run2) * window_JpsiXi_Run2;
	Float_t XibNorm_SystErr_Run2 = XibNorm_Run2 * sqrt( pow(N_Xib_data_SystErr_Run2/N_Xib_data_Run2,2) +
	                                                    pow(eff_Xib_JpsiLambda_SystErr_Run2/eff_Xib_JpsiLambda_Run2,2) +
	                                                    pow(eff_Xib_JpsiXi_SystErr_Run2/eff_Xib_JpsiXi_Run2,2) +
	                                                    pow(trackingErr,2) + pow(xiVtxUnc,2));


	cout<<XibNorm_SystErr_Run2<<endl;
	Float_t Nobs_Run1 = 191;
	Float_t Ncomb_Run1 = 18;

	Float_t Nobs_StatErr_Run1 = sqrt(Nobs_Run1);
	Float_t Ncomb_StatErr_Run1 = sqrt(Ncomb_Run1);
	Float_t N_JpsiLambda_Run1_StatErr_window = N_JpsiLambda_StatErr_Run1*window_JpsiLambda_Run1;
	Float_t N_JpsiLambda_Run1_SystErr_window = N_JpsiLambda_SystErr_Run1*window_JpsiLambda_Run1;

	Float_t Nobs_Run2 = 598;
	Float_t Ncomb_Run2 = 60;

	Float_t Nobs_StatErr_Run2                = sqrt(Nobs_Run2);
	Float_t Ncomb_StatErr_Run2               = sqrt(Ncomb_Run2);
	Float_t N_JpsiLambda_Run2_StatErr_window = N_JpsiLambda_StatErr_Run2*window_JpsiLambda_Run2;
	Float_t N_JpsiLambda_Run2_SystErr_window = N_JpsiLambda_SystErr_Run2*window_JpsiLambda_Run2;

	Float_t N_JpsiSigma_Run1         = (Nobs_Run1 - Ncomb_Run1 - N_JpsiLambda_Run1_window - XibNorm_Run1)/window_JpsiSigma_Run1;

	cout<<"Nobs_Run1 = "<<Nobs_Run1<<endl;
	cout<<"Ncomb_Run1 = "<<Ncomb_Run1<<endl;
	cout<<"N_JpsiLambda_Run1_window = "<<N_JpsiLambda_Run1_window<<endl;
	cout<<"XibNorm_Run1 = "<<XibNorm_Run1<<endl;
	cout<<"XibNorm_SystErr_Run1 = "<<XibNorm_SystErr_Run1<<endl;
	cout<<"window_JpsiSigma_Run1 = "<<window_JpsiSigma_Run1<<endl;

	Float_t N_JpsiSigma_StatErr_Run1 = sqrt(pow(Nobs_StatErr_Run1,2) + pow(Ncomb_StatErr_Run1,2) +
	                                        pow(N_JpsiLambda_Run1_StatErr_window,2) + pow(XibNorm_StatErr_Run1,2)) / window_JpsiSigma_Run1;
	Float_t N_JpsiSigma_SystErr_Run1 = sqrt(pow(N_JpsiLambda_Run1_SystErr_window,2) + pow(XibNorm_SystErr_Run1,2));

	Float_t N_JpsiSigma_Run2         = (Nobs_Run2 - Ncomb_Run2 - N_JpsiLambda_Run2_window - XibNorm_Run2)/window_JpsiSigma_Run2;
	Float_t N_JpsiSigma_StatErr_Run2 = sqrt(pow(Nobs_StatErr_Run2,2) + pow(Ncomb_StatErr_Run2,2) +
	                                        pow(N_JpsiLambda_Run2_StatErr_window,2) + pow(XibNorm_StatErr_Run2,2)) / window_JpsiSigma_Run2;
	Float_t N_JpsiSigma_SystErr_Run2 = sqrt(pow(N_JpsiLambda_Run2_SystErr_window,2) + pow(XibNorm_SystErr_Run2,2));

	Float_t eff_JpsiLambda_Run1 = 0.112307; //% Weighted
	Float_t eff_JpsiLambda_SystErr_Run1 = 0.00134492; //%

	Float_t eff_JpsiSigma_Run1 = 0.0993729; //% Weighted
	Float_t eff_JpsiSigma_SystErr_Run1 = 0.00200336; //%

	Float_t eff_JpsiLambda_Run2 = 0.121395; //% Weighted
	Float_t eff_JpsiLambda_SystErr_Run2 = 0.00087353; //%

	Float_t eff_JpsiSigma_Run2 = 0.117869; //% Weighted
	Float_t eff_JpsiSigma_SystErr_Run2 = 0.00116115; //%

	Float_t eff_ratio_Run1 = (eff_JpsiLambda_Run1/eff_JpsiSigma_Run1);
	Float_t eff_ratio_StatErr_Run1 = 0.0;
	Float_t eff_ratio_SystErr_Run1 = eff_ratio_Run1 * sqrt( pow(eff_JpsiLambda_SystErr_Run1/eff_JpsiLambda_Run1,2) +
	                                                        pow(eff_JpsiSigma_SystErr_Run1/eff_JpsiSigma_Run1,2) );

	Float_t eff_ratio_Run2 = (eff_JpsiLambda_Run2/eff_JpsiSigma_Run2);
	Float_t eff_ratio_StatErr_Run2 = 0.0;
	Float_t eff_ratio_SystErr_Run2 = eff_ratio_Run2 * sqrt( pow(eff_JpsiLambda_SystErr_Run2/eff_JpsiLambda_Run2,2) +
	                                                        pow(eff_JpsiSigma_SystErr_Run2/eff_JpsiSigma_Run2,2) );

	Float_t R_Run1 = (N_JpsiSigma_Run1/N_JpsiLambda_Run1) * (eff_ratio_Run1);
	Float_t R_StatErr_Run1 = R_Run1 * sqrt( pow(N_JpsiSigma_StatErr_Run1/N_JpsiSigma_Run1,2) +
	                                        pow(N_JpsiLambda_StatErr_Run1/N_JpsiLambda_Run1,2) +
	                                        pow(eff_ratio_StatErr_Run1/eff_ratio_Run1,2));
	Float_t R_SystErr_Run1 = R_Run1 * sqrt( pow(N_JpsiSigma_SystErr_Run1/N_JpsiSigma_Run1,2) +
	                                        pow(N_JpsiLambda_SystErr_Run1/N_JpsiLambda_Run1,2) +
	                                        pow(eff_ratio_SystErr_Run1/eff_ratio_Run1,2));
	cout<<"N_JpsiSigma_SystErr_Run1 = "<<N_JpsiSigma_SystErr_Run1<<endl;
	cout<<"N_JpsiLambda_SystErr_Run1 = "<<N_JpsiLambda_SystErr_Run1<<endl;
	cout<<"eff_ratio_SystErr_Run1 = "<<eff_ratio_SystErr_Run1<<endl;

	Float_t R_Run2 = (N_JpsiSigma_Run2/N_JpsiLambda_Run2) * (eff_ratio_Run2);
	Float_t R_StatErr_Run2 = R_Run2 * sqrt( pow(N_JpsiSigma_StatErr_Run2/N_JpsiSigma_Run2,2) +
	                                        pow(N_JpsiLambda_StatErr_Run2/N_JpsiLambda_Run2,2) +
	                                        pow(eff_ratio_StatErr_Run2/eff_ratio_Run2,2));
	Float_t R_SystErr_Run2 = R_Run2 * sqrt( pow(N_JpsiSigma_SystErr_Run2/N_JpsiSigma_Run2,2) +
	                                        pow(N_JpsiLambda_SystErr_Run2/N_JpsiLambda_Run2,2) +
	                                        pow(eff_ratio_SystErr_Run2/eff_ratio_Run2,2));

	cout<<"***************Run 1 Result***************"<<endl;
	cout<<"N_JpsiLambda = "<<N_JpsiLambda_Run1<<endl;
	cout<<"N_JpsiSigma = "<<N_JpsiSigma_Run1<<endl;
	cout<<"eff_ratio_Run1 = "<<eff_ratio_Run1<<endl;
	cout<<"R = "<<R_Run1<<" +/- "<<R_StatErr_Run1<<" +/- "<<R_SystErr_Run1<<endl;
	cout<<"******************************************"<<endl;

	cout<<"***************Run 2 Result***************"<<endl;
	cout<<"N_JpsiLambda = "<<N_JpsiLambda_Run2<<endl;
	cout<<"N_JpsiSigma = "<<N_JpsiSigma_Run2<<endl;
	cout<<"eff_ratio_Run2 = "<<eff_ratio_Run2<<endl;
	cout<<"R = "<<R_Run2<<" +/- "<<R_StatErr_Run2<<" +/- "<<R_SystErr_Run2<<endl;
	cout<<"******************************************"<<endl;
}
