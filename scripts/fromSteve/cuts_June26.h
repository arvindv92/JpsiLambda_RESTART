// add Xib- weight, w(eta)
bool useEtaWeight = true;   // def = true
bool useBachPiWeight = false;   // def = true
bool applyJpsiLWeight = true;
bool applyJpsiXiWeight = true;
bool applyLbKinWeight = true;
bool applyLbPiWeight = true;


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
float JpsiWeightPar[3] = {0.71, 0.45, 0.0};

  

void getCuts(){
  
  //Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Lb_L0Global_TIS>0)";
  //Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Xib_L0Global_TIS>0)";
  Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0)";
  Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0)";
  //Lb_Trig = "( Lb_L0Global_TIS>0>0 && Jpsi_L0DiMuonDecision_TOS==0 )";
  //Xib_Trig = "( Xib_L0Global_TIS>0 && Jpsi_L0DiMuonDecision_TOS==0 )";
  
  // default mass window on Xi = 11  MeV

  Jpsi_Cut = "abs(Jpsi_MM-3098)<40";// && Jpsi_L0DiMuonDecision_TOS>0";
  Jpsi_Cut = Jpsi_Cut && "muplus_IPCHI2_OWNPV>4&&muminus_IPCHI2_OWNPV>4&&muplus_PT>550&&muminus_PT>550&&Jpsi_ENDVERTEX_CHI2<16 && Jpsi_FDCHI2_OWNPV>12";
  Lambda_Cut  = "p_TRACK_Type==5 && pi_TRACK_Type==5 && abs(L_M-1115)<6 && L_TAU>0";
  // Next line here, to align with J/psi Xi- and Jpsi Lambda stripping line cuts
  Lambda_Cut = Lambda_Cut && "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100&&L_FD_OWNPV/L_ENDVERTEX_ZERR>5"; 
  LambdaLL_Cut  = "p_TRACK_Type==3 && pi_TRACK_Type==3 && abs(L_M-1115)<6 && L_TAU>0";
  Xi_Cut = "abs(Xi_M-L_M+1115.68-1321)<11 && Xi_DIRA_ORIVX>0.995 && Xi_ENDVERTEX_CHI2<25 && Xi_TAU>0"; // Use Long Track for pion from Xi- decay
  //XiL_Cut = "BachPi_TRACK_Type==3 && BachPi_P>5000"; // Use Long Track for pion from Xi- decay
  //XiD_Cut = "BachPi_TRACK_Type==5 && BachPi_P>5000"; // Use DownStream Track for pion from Xi- decay
  XiL_Cut = "BachPi_TRACK_Type==3"; // Use Long Track for pion from Xi- decay
  XiD_Cut = "BachPi_TRACK_Type==5"; // Use DownStream Track for pion from Xi- decay

  // Apply p and eta cut, for comparison to region where tracking corrections are valid
  //XiL_Cut = XiL_Cut && "BachPi_P>5000 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))>1.9 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))<4.9";
  

  Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_MINIPCHI2<12 && Xib_TAU*1000>0.2";
  Xib_Cut = Xib_Cut && Xib_Trig && Jpsi_Cut && Lambda_Cut && Xi_Cut && (XiL_Cut||XiD_Cut) ;
  float ml = mLowXib;
  float mh = mHiXib;
  TString mCut = Form("Xib_DTF_M_JpsiXiLConstr>%f && Xib_DTF_M_JpsiXiLConstr<%f",ml,mh);
  Xib_Cut = Xib_Cut && mCut;
  MCXib_Cut = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25))";

  Lb_Cut = "(Lb_ENDVERTEX_CHI2/Lb_ENDVERTEX_NDOF < 10) && L_TAU>0 && Lb_MINIPCHI2 < 12 && Lb_TAU*1000>0.2";
  Lambdab_Cut = Jpsi_Cut && Lambda_Cut && Lb_Cut && Lb_Trig;
  LambdabLL_Cut = Jpsi_Cut && LambdaLL_Cut && Lb_Cut && Lb_Trig;

  ml = mLowLb;
  mh = mHiLb;
  mCut = Form("Lb_DTF_M_JpsiLConstr>%f&& Lb_DTF_M_JpsiLConstr<%f",ml,mh);
  Lambdab_Cut = Lambdab_Cut && mCut;

  MCLambdab_Cut = Lambdab_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";
  MCLambdabLL_Cut = LambdabLL_Cut&&"(Lb_BKGCAT<30||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";

  //TString ss = Form("(1.0*exp(-0.5*((-log(tan(0.5*asin(Xib_PT/Xib_P))))-%4.2f)**2/%4.2f**2))",XibWtPar[0], XibWtPar[1]); // exponential
  TString ss = Form("%4.2f+%4.2f*(-log(tan(0.5*asin(Xib_PT/Xib_P))))",XibWtPar[0], XibWtPar[1]);
  XibWt = ss;//"(1.0*exp(-0.5*((-log(tan(0.5*asin(Xib_PT/Xib_P))))-6.33)**2/2.78**2))";
  ss = Form("(%6.3f+%6.3f*exp(-0.001*BachPi_P/%6.3f))",BachPiPar[0],BachPiPar[1],BachPiPar[2]);
  BachPiWt = ss;
  ss = Form("(%5.3f+%5.3f*(Jpsi_P/Lb_P)+%5.3f*(Jpsi_P/Lb_P)**2)",JpsiWeightPar[0],JpsiWeightPar[1],JpsiWeightPar[2]);
  JpsiLWeight = ss;
  ss = Form("(%5.3f+%5.3f*(Jpsi_P/Xib_P)+%5.3f*(Jpsi_P/Xib_P)**2)",JpsiWeightPar[0],JpsiWeightPar[1],JpsiWeightPar[2]);
  JpsiXiWeight = ss;
  

}

