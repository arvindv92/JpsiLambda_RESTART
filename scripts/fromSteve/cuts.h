// add Xib- weight, w(eta)
bool useEtaWeight = false;   // def = true
bool useBachPiWeight = false;   // def = true
bool applyJpsiLWeight = true;
bool applyJpsiXiWeight = true;
bool applyLbKinWeight = true;
bool applyLKinWeight = true;
bool applyPiKinWeight = true;
bool applyPrPiKinWeight = true;
//bool applyLbPiWeight = false;



bool applyXiPiTrackingAcceptanceCut = false;


bool makeLbPiROOTFile = false; // def = false (don't overwrite existing file)

float XibLifetimeWeight = 0;

double tauLb = 1.470;
double err_tauLb = 0.010;

double tauXib = 1.571;
double err_tauXib = 0.041;

double genEffLb = 0.181;
double err_genEffLb = 0.003;
double genEffXib = 0.171;
double err_genEffXib = 0.001;

double lumJpsiL0 = 2962.0;
double lumJpsiXi = 2816.0;



TCut Jpsi_Cut;
TCut Lambda_Cut;
TCut LambdaLL_Cut;
TCut LambdaDD_Cut;
TCut Xi_Cut;
TCut XiL_Cut;
TCut XiD_Cut;

TCut MCXib_Cut;
TCut Xib_Cut;
TCut Lb_Cut;
TCut Lambdab_Cut;
TCut MCLambdab_Cut;
TCut LambdabLL_Cut;
TCut MCLambdabLL_Cut;
TCut LambdabDD_Cut;
TCut MCLambdabDD_Cut;


double mLowLb=5500;
double mHiLb = 5750;
double mLowXib= 5640;
double mHiXib = 5950;
//double mLowXib= 5700;
//double mHiXib = 5900;
float bwLb = 2.0;
float bwXib = 5.0;
int nbinLb, nbinXib;

TCut Lb_Trig;
TCut Xib_Trig;

TCut XibWt;
TCut BachPiWt;

TCut JpsiLWeight;
TCut JpsiXiWeight;


//float XibWtPar[2] = {6.33, 2.78};   // exponential correction
float XibWtPar[2] = {1.0, 0.3};   // linear correction
float BachPiPar[3] = {0.746, 1.30, 2.12};
float JpsiWeightPar[3] = {0.6, 0.6, 0.0};
float JpsiXiWeightPar[3] = {0.0, 4.5, 0.0};

int nHistsLb;
int histLbf;
int nHistsL;
int histLf;
int nHistsp;
int histpf;
int nHistspi;
int histpif;
int nHistsR;
int histRf;
int histJpsif;
int nHistsJpsi;
int nHistsEx;
int histExf;
int nHistsSys;
int histSys;
int nHistsBachPi;
int histBachPif;
int nHistsXi;
int histXif;

int nHists;

int nRatio = 0;
int nHists;
int nHists2D;
TString histNames[200];
TString histExp[200];
TString histExp2[200];
TString hTit[200];
int NB[200];
double LO[200];
double HI[200];

TH1F *hda[200];
TH1F *hmc[200];
TH1F *hr[20];

TH2F *whistLb ;
TH2F *whistL;
TH2F *whistPi ;
TH2F *whistPrPi;

TTree* TreeWithWeight ;
TTree* newtree;


TString signalWeight;
TCut simCut;

TLegend* legend1;


TChain* chainl;

int mcOpt;
bool makePDF;

TString parentPart;
TString lambdaType;
TString tag;

double bottomMarginOffset = 0.13;

///////////////////////////////////////////////////////////////////////////////

void addGraphics(TH1F *h, TString Xtitle, TString Ytitle, int iCol){
  h->SetXTitle(Xtitle);
  //h->SetFillColor(30);
  int bw = h->GetBinWidth(1);
  h->SetYTitle(Ytitle);
  h->SetStats(kFALSE);
  h->SetMinimum(0.1);
  h->SetMaximum(1.3*h->GetMaximum());
  h->SetTitleSize(0.1);
  h->SetLineColor(iCol);
  h->SetMarkerColor(iCol);
  h->SetMarkerSize(0.7);
  h->SetMarkerStyle(20);
  h->GetXaxis()->SetTitleOffset(1.0);  
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetTitleSize(0.055);  
  h->GetYaxis()->SetTitleSize(0.055);
  h->GetXaxis()->SetLabelSize(0.045);  
  h->GetYaxis()->SetLabelSize(0.045);  
  h->SetNdivisions(404,"X");
  h->SetNdivisions(505,"Y");
  h->SetLineWidth(2);
  h->SetTitle("");
  
}

void getCuts(){
  
  //TCut Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Lb_L0Global_TIS>0)";
  //Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Xib_L0Global_TIS>0)";
  Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0)";
  Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0)";
  //Lb_Trig = "( Lb_L0Global_TIS>0>0 && Jpsi_L0DiMuonDecision_TOS==0 )";
  //Xib_Trig = "( Xib_L0Global_TIS>0 && Jpsi_L0DiMuonDecision_TOS==0 )";
  
  // default mass window on Xi = 11  MeV

  Jpsi_Cut = "abs(Jpsi_MM-3098)<40";// && Jpsi_L0DiMuonDecision_TOS>0";
  Jpsi_Cut = Jpsi_Cut && "muplus_IPCHI2_OWNPV>4&&muminus_IPCHI2_OWNPV>4&&muplus_PT>550&&muminus_PT>550&&Jpsi_ENDVERTEX_CHI2<16 && Jpsi_FDCHI2_OWNPV>12";
  Lambda_Cut  = "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100 && L_FD_OWNPV/L_ENDVERTEX_ZERR>5 && abs(L_M-1115.7)<8 && L_TAU>0";
  LambdaLL_Cut  = Lambda_Cut && "p_TRACK_Type==3 && pi_TRACK_Type==3";
  LambdaDD_Cut  = Lambda_Cut && "p_TRACK_Type==5 && pi_TRACK_Type==5";

  //Lambda_Cut = Lambda_Cut && "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100&&L_FD_OWNPV/L_ENDVERTEX_ZERR>5"; 
  //LambdaLL_Cut = LambdaLL_Cut && "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100&&L_FD_OWNPV/L_ENDVERTEX_ZERR>5"; 

  Xi_Cut = "abs(Xi_M-L_M+1115.7-1321.9)<8 && Xi_DIRA_ORIVX>0.995 && Xi_ENDVERTEX_CHI2<25 && Xi_TAU>0"; // Use Long Track for pion from Xi- decay
  //XiL_Cut = "BachPi_TRACK_Type==3 && BachPi_P>5000"; // Use Long Track for pion from Xi- decay
  //XiD_Cut = "BachPi_TRACK_Type==5 && BachPi_P>5000"; // Use DownStream Track for pion from Xi- decay
  XiL_Cut = "BachPi_TRACK_Type==3"; // Use Long Track for pion from Xi- decay
  XiD_Cut = "BachPi_TRACK_Type==5"; // Use DownStream Track for pion from Xi- decay

  // Apply p and eta cut, for comparison to region where tracking corrections are valid
  if(applyXiPiTrackingAcceptanceCut) XiL_Cut = XiL_Cut && "BachPi_P>5000 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))>1.9 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))<4.9";
  

  Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_MINIPCHI2<12 && Xib_TAU*1000>0.2";
  Xib_Cut = Xib_Cut && Xib_Trig && Jpsi_Cut && LambdaDD_Cut && Xi_Cut && (XiL_Cut||XiD_Cut) ;
  float ml = mLowXib;
  float mh = mHiXib;
  TString mCut = Form("Xib_DTF_M_JpsiXiLConstr>%f && Xib_DTF_M_JpsiXiLConstr<%f",ml,mh);
  Xib_Cut = Xib_Cut && mCut;
  MCXib_Cut = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25))";

  Lb_Cut = "(Lb_ENDVERTEX_CHI2/Lb_ENDVERTEX_NDOF < 10) && L_TAU>0 && Lb_MINIPCHI2 < 12 && Lb_TAU*1000>0.2";
  LambdabDD_Cut = Jpsi_Cut && LambdaDD_Cut && Lb_Cut && Lb_Trig;
  LambdabLL_Cut = Jpsi_Cut && LambdaLL_Cut && Lb_Cut && Lb_Trig;

  ml = mLowLb;
  mh = mHiLb;
  mCut = Form("Lb_DTF_M_JpsiLConstr>%f&& Lb_DTF_M_JpsiLConstr<%f",ml,mh);
  LambdabDD_Cut = LambdabDD_Cut && mCut;
  LambdabLL_Cut = LambdabLL_Cut && mCut;

  MCLambdabDD_Cut = LambdabDD_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";
  MCLambdabLL_Cut = LambdabLL_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";

  //TString ss = Form("(1.0*exp(-0.5*((-log(tan(0.5*asin(Xib_PT/Xib_P))))-%4.2f)**2/%4.2f**2))",XibWtPar[0], XibWtPar[1]); // exponential
  TString ss = Form("%4.2f+%4.2f*(-log(tan(0.5*asin(Xib_PT/Xib_P))))",XibWtPar[0], XibWtPar[1]);
  XibWt = ss;//"(1.0*exp(-0.5*((-log(tan(0.5*asin(Xib_PT/Xib_P))))-6.33)**2/2.78**2))";
  ss = Form("(%6.3f+%6.3f*exp(-0.001*BachPi_P/%6.3f))",BachPiPar[0],BachPiPar[1],BachPiPar[2]);
  BachPiWt = ss;
  ss = Form("(%5.3f+%5.3f*(Jpsi_P/Lb_P)+%5.3f*(Jpsi_P/Lb_P)**2)",JpsiWeightPar[0],JpsiWeightPar[1],JpsiWeightPar[2]);
  JpsiLWeight = ss;
  ss = Form("(%5.3f+%5.3f*(Jpsi_P/Xib_P)+%5.3f*(Jpsi_P/Xib_P)**2)",JpsiXiWeightPar[0],JpsiXiWeightPar[1],JpsiXiWeightPar[2]);
  JpsiXiWeight = ss;
 

}

void bookLbBaseHistos(){

  signalWeight = "nsig1_sw";
  
 TString HB = "#Lambda";

 int i1 = 0;
 int i2 = 0;
 
 double maxDecayTime = 1.0;
 double maxZ = 2600;
 if(lambdaType=="LL"){
   maxDecayTime = 0.1;
   maxZ = 800;
 }
 

 histLf = 0;
 // Lambda0
 histNames[i1] = "taulam"; histExp[i1] = "L_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=50; 
 LO[i1]=0.0; HI[i1]=maxDecayTime; i1++;
 histNames[i1] = "plam"; histExp[i1] = "0.001*L_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
 LO[i1]=0.0; HI[i1]=120.0; i1++;
 histNames[i1] = "ptlam"; histExp[i1] = "0.001*L_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=50; 
 LO[i1]=0.0; HI[i1]=10.0; i1++;
 histNames[i1] = "fdlam"; histExp[i1] = "L_FD_OWNPV"; hTit[i1] = HB+" flight distance [mm]"; NB[i1]=26; 
 LO[i1]=0.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "vzlam"; histExp[i1] = "L_ENDVERTEX_Z"; hTit[i1] = HB+" z decay position [mm]"; NB[i1]=26; 
 LO[i1]=-100.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "etalam"; histExp[i1] = "-log(tan(0.5*asin(L_PT/L_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=14; 
 LO[i1]=1.5; HI[i1]=5.0; i1++;
 nHistsL = i1;

 if(parentPart=="Lb"){ 
   HB = "#Lambda_{b}^{0}";
   histLbf = i1;
   histNames[i1] = "plamb"; histExp[i1] = "0.001*Lb_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
   LO[i1]=0.0; HI[i1]=300.0; i1++;
   histNames[i1] = "ptlamb"; histExp[i1] = "0.001*Lb_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=20.0; i1++;
   histNames[i1] = "etalamb"; histExp[i1] = "-log(tan(0.5*asin(Lb_PT/Lb_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=20; 
   LO[i1]=1.5; HI[i1]=6.5; i1++;
 }else{
   signalWeight = "nsig2_sw";
   HB = "#Xi_{b}^{#font[122]{-}}";
   histLbf = i1;
   histNames[i1] = "pxib"; histExp[i1] = "0.001*Xib_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
   LO[i1]=0.0; HI[i1]=300.0; i1++;
   histNames[i1] = "ptxib"; histExp[i1] = "0.001*Xib_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=20.0; i1++;
   histNames[i1] = "etaxib"; histExp[i1] = "-log(tan(0.5*asin(Xib_PT/Xib_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=20; 
   LO[i1]=1.5; HI[i1]=6.5; i1++;
 }
 

 nHistsLb = i1 - histLbf;

 HB = "proton (from #Lambda)";
 histpf = i1;
 histNames[i1] = "pp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2+p_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=120.0; i1++;
 histNames[i1] = "ptp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 nHistsp = i1 - histpf;
 HB = "pion (from #Lambda)";
 histpif = i1;
 histNames[i1] = "ppi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2+pi_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=30; 
 LO[i1]=0.0; HI[i1]=30.0; i1++;
 histNames[i1] = "ptpi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=2.0; i1++;
 nHists = i1;
 nHistspi = i1 - histpif ;

 histRf = i1;
 histNames[i1] = "jpsioverLb"; histExp[i1] = "Jpsi_P/Lb_P"; hTit[i1] = "p_{J/#psi} / p_{#Lambda_{b}}"; NB[i1]=12; 
 LO[i1]=0.25; HI[i1]=1.0; i1++;
 histNames[i1] = "LoverLb"; histExp[i1] = "L_P/Lb_P"; hTit[i1] = "p_{#Lambda} / p_{#Lambda_{b}}"; NB[i1]=25; 
 LO[i1]=0.0; HI[i1]=0.75; i1++;

 nHistsR = i1 - histRf ;
 
 histJpsif = i1;
 histNames[i1] = "jpsipt"; histExp[i1] = "0.001*Jpsi_PT"; hTit[i1] = "J/#psi p_{T} [GeV/#it{c}]"; NB[i1]=36; 
 LO[i1]=0.0; HI[i1]=18.0; i1++;
 nHistsJpsi = i1 - histJpsif ;

 if(lambdaType=="LL"){
   HB = "proton (from #Lambda)";
   histExf = i1;
   histNames[i1] = "pIP"; histExp[i1] = "p_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=4.0; i1++;
   histNames[i1] = "pIPCHI2"; histExp[i1] = "log(p_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "pDOCAz"; histExp[i1] = "abs(p_PY*(pi_ORIVX_X-0.459)-p_PX*(pi_ORIVX_Y+0.014))/sqrt(p_PX**2+p_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=30; LO[i1]=0.0; HI[i1]=3.0; i1++;

   HB = "pion (from #Lambda)";
   histNames[i1] = "piIP"; histExp[i1] = "pi_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=10.0; i1++;
   histNames[i1] = "piIPCHI2"; histExp[i1] = "log(pi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "piDOCAz"; histExp[i1] = "abs(pi_PY*(pi_ORIVX_X-0.459)-pi_PX*(pi_ORIVX_Y+0.014))/sqrt(pi_PX**2+pi_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=32; LO[i1]=0.0; HI[i1]=8.0; i1++;
   nHistsEx = i1 - histExf;


 } 


 // Do these last
 HB = "proton (from #Lambda)";
 histSys = i1;
 histNames[i1] = "pIPCHI2_sys"; histExp[i1] = "p_IPCHI2_OWNPV"; hTit[i1] = HB+" #chi^{2}_{IP}"; NB[i1]=100; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 HB = "pion (from #Lambda)";
 histNames[i1] = "piIPCHI2_sys"; histExp[i1] = "pi_IPCHI2_OWNPV"; hTit[i1] = HB+" #chi^{2}_{IP}"; NB[i1]=100; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 nHistsSys = i1 - histSys;
 

 nHists = i1;

 cout << "# Histograms = " << nHists << endl;
 
 
 // Book the histograms
 for(int i=0; i<nHists; i++){
   hda[i] = new TH1F(histNames[i]+"_da",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i] = new TH1F(histNames[i]+"_mc",hTit[i],NB[i],LO[i],HI[i]);
 }

}

void bookXibBaseHistos(){

  signalWeight = "nsig2_sw";
  
 TString HB = "#Lambda";

 int i1 = 0;
 int i2 = 0;
 
 double maxDecayTime = 0.6;
 double maxZ = 2600;
 if(lambdaType=="LL"){
   maxDecayTime = 0.1;
   maxZ = 800;
 }
 

 histLf = 0;
 // Lambda0
 histNames[i1] = "taulam"; histExp[i1] = "L_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=12; 
 LO[i1]=0.0; HI[i1]=maxDecayTime; i1++;
 histNames[i1] = "plam"; histExp[i1] = "0.001*L_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 histNames[i1] = "ptlam"; histExp[i1] = "0.001*L_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=12; 
 LO[i1]=0.0; HI[i1]=6.0; i1++;
 histNames[i1] = "fdlam"; histExp[i1] = "L_FD_OWNPV"; hTit[i1] = HB+" flight distance [mm]"; NB[i1]=13; 
 LO[i1]=0.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "vzlam"; histExp[i1] = "L_ENDVERTEX_Z"; hTit[i1] = HB+" z decay position [mm]"; NB[i1]=13; 
 LO[i1]=-100.0; HI[i1]=maxZ; i1++;
 histNames[i1] = "etalam"; histExp[i1] = "-log(tan(0.5*asin(L_PT/L_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=14; 
 LO[i1]=1.5; HI[i1]=5.0; i1++;
 nHistsL = i1;

 HB = "#Xi_{b}^{#font[122]{-}}";
 histLbf = i1;
 histNames[i1] = "pxib"; histExp[i1] = "0.001*Xib_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=300.0; i1++;
 histNames[i1] = "ptxib"; histExp[i1] = "0.001*Xib_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=20.0; i1++;
 histNames[i1] = "etaxib"; histExp[i1] = "-log(tan(0.5*asin(Xib_PT/Xib_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=10; 
 LO[i1]=1.5; HI[i1]=6.5; i1++;
 

 nHistsLb = i1 - histLbf;

 HB = "proton (from #Lambda)";
 histpf = i1;
 histNames[i1] = "pp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2+p_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=80.0; i1++;
 histNames[i1] = "ptp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=15; 
 LO[i1]=0.0; HI[i1]=6.0; i1++;
 nHistsp = i1 - histpf;
 HB = "pion (from #Lambda)";
 histpif = i1;
 histNames[i1] = "ppi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2+pi_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=20; 
 LO[i1]=0.0; HI[i1]=20.0; i1++;
 histNames[i1] = "ptpi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=15; 
 LO[i1]=0.0; HI[i1]=1.5; i1++;
 nHistspi = i1 - histpif ;

 histRf = i1;
 histNames[i1] = "jpsioverLb"; histExp[i1] = "Jpsi_P/Xib_P"; hTit[i1] = "p_{J/#psi} / p_{#Xi_{b}}"; NB[i1]=6; 
 LO[i1]=0.25; HI[i1]=1.0; i1++;
 histNames[i1] = "LoverLb"; histExp[i1] = "L_P/Xib_P"; hTit[i1] = "p_{#Lambda} / p_{#Xi_{b}}"; NB[i1]=6; 
 LO[i1]=0.0; HI[i1]=0.75; i1++;

 nHistsR = i1 - histRf ;
 
 histJpsif = i1;
 histNames[i1] = "jpsipt"; histExp[i1] = "0.001*Jpsi_PT"; hTit[i1] = "J/#psi p_{T} [GeV/#it{c}]"; NB[i1]=18; 
 LO[i1]=0.0; HI[i1]=18.0; i1++;
 nHistsJpsi = i1 - histJpsif ;

 if(lambdaType=="LL"){
   HB = "proton (from #Lambda)";
   histExf = i1;
   histNames[i1] = "pIP"; histExp[i1] = "p_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=4.0; i1++;
   histNames[i1] = "pIPCHI2"; histExp[i1] = "log(p_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "pDOCAz"; histExp[i1] = "abs(p_PY*(pi_ORIVX_X-0.459)-p_PX*(pi_ORIVX_Y+0.014))/sqrt(p_PX**2+p_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=30; LO[i1]=0.0; HI[i1]=3.0; i1++;

   HB = "pion (from #Lambda)";
   histNames[i1] = "piIP"; histExp[i1] = "pi_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
   LO[i1]=0.0; HI[i1]=10.0; i1++;
   histNames[i1] = "piIPCHI2"; histExp[i1] = "log(pi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
   LO[i1]=2.0; HI[i1]=11.0; i1++;
   histNames[i1] = "piDOCAz"; histExp[i1] = "abs(pi_PY*(pi_ORIVX_X-0.459)-pi_PX*(pi_ORIVX_Y+0.014))/sqrt(pi_PX**2+pi_PY**2)"; 
   hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=32; LO[i1]=0.0; HI[i1]=8.0; i1++;
   nHistsEx = i1 - histExf;
 }

   
 HB = "#pi^{#font[122]{-}} (from #Xi^{#font[122]{-}})";
 histBachPif = i1; 
 histNames[i1] = "ppib"; histExp[i1] = "0.001*BachPi_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=18; 
 LO[i1]=0.0; HI[i1]=18.0; i1++;
 histNames[i1] = "ptpib"; histExp[i1] = "0.001*BachPi_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=15; 
 LO[i1]=0.0; HI[i1]=1.5; i1++;
 histNames[i1] = "pibIPCHI2"; histExp[i1] = "log(BachPi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=20; 
 LO[i1]=1.0; HI[i1]=10.0; i1++;
 nHistsBachPi = i1 - histBachPif ;


 HB = "#Xi^{#font[122]{-}}";
 histXif = i1; 
 histNames[i1] = "pXi"; histExp[i1] = "0.001*Xi_P"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=10; 
 LO[i1]=0.0; HI[i1]=100.0; i1++;
 histNames[i1] = "ptXi"; histExp[i1] = "0.001*Xi_PT"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=12; 
 LO[i1]=0.0; HI[i1]=6.0; i1++;
 histNames[i1] = "etaXi"; histExp[i1] = "-log(tan(0.5*asin(Xi_PT/Xi_P)))"; hTit[i1] = HB+" #eta"; NB[i1]=14; 
 LO[i1]=1.5; HI[i1]=5.0; i1++;
 nHistsXi = i1 - histXif ;


 nHists = i1;


 cout << "Number of histograms " << nHists << endl;
 


 nHists = i1;
 histSys = nHists;
 
 
 // Book the histograms
 for(int i=0; i<nHists; i++){
   hda[i] = new TH1F(histNames[i]+"_da",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i] = new TH1F(histNames[i]+"_mc",hTit[i],NB[i],LO[i],HI[i]);
 }

}

void fillLbHistos(TTree* tree_data, TTree *tree_sim){
  
  TString hist;
  //simCut = MCLambdabDD_Cut;
  
  for(int i=0; i<nHists; i++){
    hist = histNames[i]+"_da";
    cout << "Filling histogram: " << i << " " << hTit[i] << " 0 " << endl;
    tree_data->Draw(histExp[i]+">>"+hist,signalWeight);
    addGraphics(hda[i],hTit[i],"",1);
    hist = histNames[i]+"_mc";
    cout << "Filling histogram: " << i << " " << hTit[i] << " 1 " << endl;
    tree_sim->Draw(histExp[i]+">>"+hist,simCut);    
    addGraphics(hmc[i],hTit[i],"",2);
    hmc[i]->Sumw2();
    if(i<histSys) hmc[i]->Scale(hda[i]->Integral()/hmc[i]->Integral());
    cout << "Filling histogram: " << i << " " << hTit[i] << " 2 " << endl;
  }

  legend1 = new TLegend(0.70,0.70,0.88,0.88);
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->SetTextSize(0.055);
  legend1->AddEntry(hda[0],"Data","LEP");
  legend1->AddEntry(hmc[0],"Sim","L");
  legend1->Draw();

  
}

TChain* getLbDDFiles(){

  chainl = new TChain("","MyTuple");
  // reco'd
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  TString sLb_Trig = Lb_Trig;
  tag = "DD";

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if(sLb_Trig.Contains("L0Global")>0 && sLb_Trig.Contains("L0DiMuon")>0) {
    tag = tag + "_L0JPsiTOSorGlobalTIS";
    cout << "Trigger is (JpsiL0DiMuon TOS || Lb L0 Global TIS)" << endl;
  }else if(sLb_Trig.Contains("L0DiMuon")>0 && sLb_Trig.Contains("L0Global")==0  ){
    tag = tag + "_L0JPsiTOS";
    cout << "Trigger is (JpsiL0DiMuon TOS)" << endl;
  }else if(sLb_Trig.Contains("L0DiMuon")==0 && sLb_Trig.Contains("L0Global")>0 ){
    tag = tag + "_L0GlobalTIS";
    cout << "Trigger is (Lb L0 Global TIS)" << endl;
  }else{
    cout << "Shoot, something wrong here!!" << endl;
  }    
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  

  TFile *f = new TFile("Lb2JpsiL0_"+tag+".root");
  TreeWithWeight = (TTree*)f->Get("TreeWithWeight");

  cout << "TreeWithWeight = " << TreeWithWeight  << endl;
  
  return chainl;
}

TChain* getLbLLFiles(){

  chainl = new TChain("","MyTuple");
  // reco'd
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_0/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_1/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_2/jpsilambda.root/Lb2JpsiLTree/MyTuple");

  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagDown/922_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2012_MagUp/923_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagDown/920_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiLambda/massdump/mc/JpsiLambda/Pythia8/2011_MagUp/921_3/jpsilambda.root/Lb2JpsiLTree/MyTuple");


  TString sLb_Trig = Lb_Trig;
  tag = "LL";

  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  if(sLb_Trig.Contains("L0Global")>0 && sLb_Trig.Contains("L0DiMuon")>0) {
    tag = tag + "_L0JPsiTOSorGlobalTIS";
    cout << "Trigger is (JpsiL0DiMuon TOS || Lb L0 Global TIS)" << endl;
  }else if(sLb_Trig.Contains("L0DiMuon")>0 && sLb_Trig.Contains("L0Global")==0  ){
    tag = tag + "_L0JPsiTOS";
    cout << "Trigger is (JpsiL0DiMuon TOS)" << endl;
  }else if(sLb_Trig.Contains("L0DiMuon")==0 && sLb_Trig.Contains("L0Global")>0 ){
    tag = tag + "_L0GlobalTIS";
    cout << "Trigger is (Lb L0 Global TIS)" << endl;
  }else{
    cout << "Shoot, something wrong here!!" << endl;
  }    
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  

  TFile *f = new TFile("Lb2JpsiL0_"+tag+".root");
  TreeWithWeight = (TTree*)f->Get("TreeWithWeight");

  cout << "TreeWithWeight = " << TreeWithWeight  << endl;
  
  return chainl;

}

TChain* getXibFiles(){

  chainl = new TChain("","MyTuple");
  // reco'd
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagDown/1001_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagUp/1002_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagDown/1003_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagUp/1004_0/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagDown/1001_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2011_MagUp/1002_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagDown/1003_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");
  chainl->Add("/data1/avenkate/JpsiXi/mc/Pythia8/2012_MagUp/1004_1/jpsixi.root/Xib2JpsiXiTree/MyTuple");

  TString fileName = "Xib2JpsiXi.root";;
  if(applyXiPiTrackingAcceptanceCut) fileName = "Xib2JpsiXi_withTrackingAcc.root";
  TFile *f = new TFile(fileName);
  TreeWithWeight = (TTree*)f->Get("TreeWithWeight");
  cout << "================================="<<endl;
  cout << TreeWithWeight << endl;
  cout << "================================="<<endl;
  
  return chainl;
  
}


void getWeightHistos(){
  TFile *fw = new TFile("fsig_Lb_new.root");
  //TH2F *whistLb = (TH2F*)fw->Get("hratio");
  whistLb = (TH2F*)fw->Get("hratioLb");
  TFile *fw2 = new TFile("fsig_Lb_new2.root");
  whistL = (TH2F*)fw2->Get("hratioL");
  TFile *fw3 = new TFile("fsig_Lb_new3.root");
  whistPi = (TH2F*)fw3->Get("hratioPi");
  whistPrPi = (TH2F*)fw3->Get("hratioPrPi");

}

TTree* getNewTree(TChain *chain){

  int ib;
  double eta, w, wPi, etaPi, pi_PX, pi_PY, pi_PZ, p_PX, p_PY, p_PZ, ptPi, pPi, pPr, asym;
  double w1, w2, w3, pi_PT, p_PT;
  double L_P, L_PT, Lb_P, Lb_PT, L_TAU;
  chain->SetBranchAddress("Lb_P", &Lb_P);
  chain->SetBranchAddress("Lb_PT", &Lb_PT);
  chain->SetBranchAddress("L_P", &L_P);
  chain->SetBranchAddress("L_PT", &L_PT);
  chain->SetBranchAddress("pi_PT", &pi_PT);
  chain->SetBranchAddress("pi_PZ", &pi_PZ);
  chain->SetBranchAddress("p_PT", &p_PT);
  chain->SetBranchAddress("p_PZ", &p_PZ);
  Long64_t nentries = chain->GetEntriesFast();
  cout << nentries << endl;
  Long64_t ientry ;
  //nentries = 100;
  TFile *ff = new TFile("/data1/sblusk/junk2.root","RECREATE");
  double swt, swtPi, swtLb, swtL, swtPi;
  newtree = (TTree*) chain->CloneTree(0);
  newtree->Branch("swtLb", &swtLb);
  newtree->Branch("swtL", &swtL);
  newtree->Branch("swtPi", &swtPi);
  for (Long64_t ii=0;ii<nentries; ii++) {
    ientry = chain->LoadTree(ii);
    if (ientry < 0) break;
    chain->GetEntry(ii);
    if(ii%10000 == 0) cout << "At entry = " << ii << endl;
    swtPi = 1.0;
    if(mcOpt==4){
      swtPi = get2DKinWt(whistPi,pi_PT,pi_PZ,0);
    }else if(mcOpt==5){
      pPr = sqrt(p_PT**2+p_PZ**2);
      pPi = sqrt(pi_PT**2+pi_PZ**2);
      asym = (pPr-pPi) / (pPr+pPi);      
      swtPi = get1DKinWt(whistPrPi,asym);
    }

    swtLb = get2DKinWt(whistLb,Lb_PT,Lb_P,1);
    swtL = get2DKinWt(whistL,L_PT,L_P,1);

    newtree->Fill(); 
  }
  newtree->AutoSave();
  TCut swLb = "(swtLb)";
  TCut swL = "(swtL)";
  TCut swPi = "(swtPi)";
  if(mcOpt>=2 && applyLbKinWeight) simCut = simCut*swLb;
  if(mcOpt>=3 && applyLKinWeight) simCut = simCut*swL;
  if(mcOpt>=4 && applyPrPiKinWeight) simCut = simCut*swPi;  
  if(applyJpsiLWeight) simCut = simCut*JpsiLWeight;
  
  return newtree;
  
}

TTree* getNewTreeXib(TChain* chain){

  int ib;
  double eta, w, wPi, etaPi, pi_PT, pi_PZ, p_PT, p_PZ, ptPi, pPi, pPr, asym;
  double L_P, L_PT, Xib_P, Xib_PT, Xi_P, Xi_PT, L_TAU, BachPi_P;
  chain->SetBranchAddress("Xib_P", &Xib_P);
  //chain->SetBranchAddress("BachPi_P", &BachPi_P);
  chain->SetBranchAddress("Xib_PT", &Xib_PT);
  chain->SetBranchAddress("L_P", &L_P);
  chain->SetBranchAddress("L_PT", &L_PT);
  //chain->SetBranchAddress("Xi_P", &Xi_P);
  //chain->SetBranchAddress("Xi_PT", &Xi_PT);
  chain->SetBranchAddress("pi_PT", &pi_PT);
  chain->SetBranchAddress("pi_PZ", &pi_PZ);
  chain->SetBranchAddress("p_PT", &p_PT);
  chain->SetBranchAddress("p_PZ", &p_PZ);
  Long64_t nentries = chain->GetEntriesFast();
  cout << nentries << endl;
  Long64_t ientry ;
  TFile *ff = new TFile("/data1/sblusk/junk2.root","RECREATE");
  double swt, swtPi, swtLb, swtL, swtPi;
  newtree = (TTree*) chain->CloneTree(0);
  newtree->Branch("swtLb", &swtLb);
  newtree->Branch("swtL", &swtL);
  newtree->Branch("swtPi", &swtPi);
  for (Long64_t i=0;i<nentries; i++) {
    ientry = chain->LoadTree(i);
    if (ientry < 0) break;
    chain->GetEntry(i);
    if(i%10000 == 0) cout << "At entry = " << i << endl;
    swtPi = 1.0;

    if(mcOpt==4){
      swtPi = get2DKinWt(whistPi,pi_PT,pi_PZ,0);
    }else if(mcOpt==5){
      pPr = sqrt(p_PT**2+p_PZ**2);
      pPi = sqrt(pi_PT**2+pi_PZ**2);
      asym = (pPr-pPi) / (pPr+pPi);
      swtPi = get1DKinWt(whistPrPi,asym);
    }

    swtLb = get2DKinWt(whistLb,Xib_PT,Xib_P,1);
    swtL = get2DKinWt(whistL,L_PT,L_P,1);

    newtree->Fill(); 
    //if(i<100) cout << ib << " " << Xib_PT << " " << etaXib << " " << w << endl;
  }
  newtree->AutoSave();
  TCut swLb = "(swtLb)";
  TCut swL = "(swtL)";
  TCut swPi = "(swtPi)";
  simCut = (simCut && XiL_Cut);
  // Weights from Lb samples
  if(mcOpt>=2 && applyLbKinWeight) simCut = simCut*swLb;
  if(mcOpt>=3 && applyLKinWeight) simCut = simCut*swL;
  if(mcOpt>=4 && applyPrPiKinWeight) simCut = simCut*swPi;
  // Extra weight for Xib- decay
  if(useEtaWeight) simCut = simCut*XibWt;
  if(applyJpsiXiWeight) simCut = simCut*JpsiXiWeight ;
  
  return newtree;


}


void makeWeightFile(){

  if(mcOpt <= 3){    
    TFile *fout;
    if(mcOpt==1) fout = new TFile("fsig_Lb_new.root","RECREATE");
    if(mcOpt==2) fout = new TFile("fsig_Lb_new2.root","RECREATE");
    if(mcOpt==3) fout = new TFile("fsig_Lb_new3.root","RECREATE");
    TCanvas *c2 = new TCanvas("c2","",400,400);

    TH2F *h201 = new TH2F("h201","eta vs pT",25,0,25,10,1.5,6.5);
    TH2F *h301 = new TH2F("h301","eta vs pT",25,0,25,10,1.5,6.5);
    TreeWithWeight->Draw("-log(tan(0.5*asin(Lb_PT/Lb_P))):0.001*Lb_PT>>h201","nsig1_sw");
    newtree->Draw("-log(tan(0.5*asin(Lb_PT/Lb_P))):0.001*Lb_PT>>h301",MCLambdabDD_Cut);
    TH2F *hratioLb = (TH2F*)h201->Clone("hratioLb");
    hratioLb->Divide(h301);    
    h201->Write();
    h301->Write();
    hratioLb->Write();

    TH2F *h202 = new TH2F("h202","eta vs pT",25,0,10,7,1.5,5.0);
    TH2F *h302 = new TH2F("h302","eta vs pT",25,0,10,7,1.5,5.0);
    TreeWithWeight->Draw("-log(tan(0.5*asin(L_PT/L_P))):0.001*L_PT>>h202","nsig1_sw");
    newtree->Draw("-log(tan(0.5*asin(L_PT/L_P))):0.001*L_PT>>h302",MCLambdabDD_Cut);
    TH2F *hratioL = (TH2F*)h202->Clone("hratioL");
    hratioL->Divide(h302);    
    h202->Write();
    h302->Write();
    hratioL->Write();

    TH2F *h203 = new TH2F("h203","eta vs pT",20,0,2,8,1.5,5.5);
    TH2F *h303 = new TH2F("h303","eta vs pT",20,0,2,8,1.5,5.5);
    TreeWithWeight->Draw("-log(tan(0.5*atan(sqrt(pi_PX**2+pi_PY**2)/pi_PZ))):0.001*sqrt(pi_PX**2+pi_PY**2)>>h203","nsig1_sw");
    newtree->Draw("-log(tan(0.5*atan(sqrt(pi_PX**2+pi_PY**2)/pi_PZ))):0.001*sqrt(pi_PX**2+pi_PY**2)>>h303",MCLambdabDD_Cut);
    TH2F *hratioPi = (TH2F*)h203->Clone("hratioPi");
    hratioPi->Divide(h303);    
    h203->Write();
    h303->Write();
    hratioPi->Write();

    TH1F *h204 = new TH1F("h204","Momentum asymmetry",100,0.5,1);
    TH1F *h304 = new TH1F("h304","Momentum asymmetry",100,0.5,1);
    TreeWithWeight->Draw("(sqrt(p_PX**2+p_PY**2+p_PZ**2)-sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))/(sqrt(p_PX**2+p_PY**2+p_PZ**2)+sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))>>h204","nsig1_sw");
    newtree->Draw("(sqrt(p_PX**2+p_PY**2+p_PZ**2)-sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))/(sqrt(p_PX**2+p_PY**2+p_PZ**2)+sqrt(pi_PX**2+pi_PY**2+pi_PZ**2))>>h304",MCLambdabDD_Cut);
    TH1F *hratioPrPi = (TH1F*)h204->Clone("hratioPrPi");
    hratioPrPi->Divide(h304);    
    h204->Write();
    h304->Write();
    hratioPrPi->Write();

    fout->Write();
    fout->Close();
  }


}

void makeLbBasePlots(){
  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);


  // Plot results
  TCanvas *c = new TCanvas("c","",1200,800);
  c->Divide(3,2);
  c->cd(1);
  int j,k;
  for(int i=0;i<nHistsL;i++){
    c->cd(i+1);
    c->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histLf + i;
    cout << i << " " << j << endl;    
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }

  TCanvas *ca = new TCanvas("ca","",1200,400);
  ca->Divide(3,1);
  for(int i=0;i<nHistsLb;i++){
    ca->cd(i+1);
    ca->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histLbf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }

  
  TCanvas *cb = new TCanvas("cb","",800,800);
  cb->Divide(2,2);
  for(int i=0;i<nHistsp+nHistspi;i++){
    cb->cd(i+1);
    cb->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histpf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }

  double v1[2];
  TF1 *func = new TF1("func","[0]+[1]*x",0,1);
  TCanvas *cc = new TCanvas("cc","",800,800);
  cc->Divide(2,2);
  for(int i=0;i<nHistsR;i++){
    cc->cd(i+1);
    cc->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histRf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==2) legend1->Draw();

    k = nRatio+i;
    cc->cd(i+3);
    hr[k] = (TH1F*)hda[j]->Clone(Form("hr_%d",i));
    hr[k]->Divide(hmc[j]);
    hr[k]->SetMaximum(2);
    if(i==0) func->SetRange(0.32,0.94);
    if(i==1) func->SetRange(0.08,0.66);
    hr[i]->Fit("func","RLM");
    v1[0] = func->GetParameter(0);
    v1[1] = func->GetParameter(1);
    if(i==0) myLatex->DrawLatex(0.18,0.8,Form("%4.2f+%5.3f p_{J/#psi}",v1[0],v1[1]));
    if(i==1) myLatex->DrawLatex(0.18,0.8,Form("%4.2f+%5.3f p_{#Lambda}",v1[0],v1[1]));    
  }
  nRatio = nRatio + 2;

  TCanvas *cd = new TCanvas("cd","",500,500);
  cd->Divide(1,1);
  for(int i=0;i<nHistsJpsi;i++){
    cd->cd(i+1);
    cd->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histJpsif + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }  

  if(makePDF){
    makePlots(c,"Lambda_from_"+parentPart+tag);
    makePlots(ca,"Lambdab_from_"+parentPart+tag);
    makePlots(cb,"ProtonPionMomenta_from_"+parentPart+tag);
    makePlots(cc,"MomentumRatio_from_"+parentPart+tag);
    makePlots(cd,"Jpsi_from_"+parentPart+tag);
  }
}

void makeLbPlotsLL(){
  int j,k;

  TCanvas *ce = new TCanvas("ce","",800,600);
  ce->Divide(3,2);
  for(int i=0;i<nHistsEx;i++){
    ce->cd(i+1);
    ce->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histExf + i;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }
  
  if(makePDF){
    makePlots(ce,"ExtraLambda_from_"+parentPart+tag);
  } 
}

void makeXibPlots(){

  cout << histBachPif << " " << nHistsBachPi << endl;

  int j,k;
  TCanvas *cf = new TCanvas("cf","BachPi",900,300);
  cf->Divide(3,1);
  for(int i=0;i<nHistsBachPi;i++){
    cf->cd(i+1);
    cf->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histBachPif + i;
    //cout << " j =  " << j << " " << hTit[j] << endl;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }

  TCanvas *cg = new TCanvas("cg","Xi",900,300);
  cg->Divide(3,1);
  for(int i=0;i<nHistsXi;i++){
    cg->cd(i+1);
    cg->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = histXif + i;
    //cout << " j =  " << j << " " << hTit[j] << endl;
    hda[j]->Draw("e");
    hmc[j]->Draw("hist,same");  
    if(i==0) legend1->Draw();
  }
  
  if(makePDF){
    makePlots(cf,"BachPi_from_"+parentPart+tag);
    makePlots(cg,"Xi_from_"+parentPart+tag);
  } 
}


void getSyst(){

  double frd[10];
  double frm[10];
  
  for(int i=0;i<nHistsSys;i++){
    j = histSys + i;
    frd[i] = hda[j]->Integral(9,25)/hda->GetEntries();
    frm[i] = hmc[j]->Integral(9,25)/hmc->GetEntries();
    cout << hTit[j] << endl;
    cout << " --> Data fraction: " << frd[i] << endl;
    cout << " --> MC fraction: " << frm[i] << endl;
  }
  
}

void makePlots(TCanvas *cc, TString filename){

  if( !applyLbKinWeight ) tagCond = "_unweighted";
  
  cc->Print("./Plots/"+filename+tagCond+".pdf");
  cc->Print("./Plots/"+filename+tagCond+".png");
  cc->Print("./Plots/"+filename+tagCond+".eps");
  cc->Print("./Plots/"+filename+tagCond+".C");
  
}

double get2DKinWt(TH2F* h, double pt, double p, int iopt = 0){
  // iopt = 0 --> provide pt, pz
  // iopt = 1 --> provide pt, pmag 

  double eta;
  double w = 1.0;
  if(p<=0 || pt<=0) return 1.0;
  if(iopt==0) eta = -log(tan(0.5*atan(pt/p)));
  if(iopt==1) eta = -log(tan(0.5*asin(pt/p)));
  int ib = h->FindBin(0.001*pt,eta);
  w = h->GetBinContent(ib);
  if(w<=0 || w>5) w = 1.0;

  return w;
}

double get1DKinWt(TH1F* h, double x){

  
  double w = 1.0;
  int ib = h->FindBin(x);
  w = h->GetBinContent(ib);
  if(w<=0 || w>5) w = 1.0;

  return w;
}

