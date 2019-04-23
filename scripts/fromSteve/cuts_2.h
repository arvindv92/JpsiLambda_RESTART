double tauLb = 1.470;
double err_tauLb = 0.010;

double tauXib = 1.571;
double err_tauXib = 0.041;

double genEffLb = 0.181;
double err_genEffLb = 0.003;
double genEffXib = 0.171;
double err_genEffXib = 0.001;



TCut Jpsi_Cut;
TCut Lambda_Cut;
TCut Xi_Cut;
TCut XiL_Cut;
TCut XiD_Cut;

TCut MCXib_Cut;
TCut Xib_Cut;
TCut Lb_Cut;
TCut Lambdab_Cut;
TCut MCLambdab_Cut;


double mLowLb=5500;
double mHiLb = 5750;
double mLowXib= 5620;
double mHiXib = 5950;
//double mLowXib= 5700;
//double mHiXib = 5900;
float bwLb = 2.0;
float bwXib = 5.0;
int nbinLb, nbinXib;


void getCuts(){
  
  TCut Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Lb_L0Global_TIS>0)";
  TCut Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Xib_L0Global_TIS>0)";
  

  Jpsi_Cut = "abs(Jpsi_MM-3098)<40";// && Jpsi_L0DiMuonDecision_TOS>0";
  Lambda_Cut  = "p_TRACK_Type==5 && pi_TRACK_Type==5 && abs(L_M-1115)<6 && L_TAU>0";
  Xi_Cut = "abs(Xi_M-1321)<11 && Xi_DIRA_ORIVX>0.999 && Xi_ENDVERTEX_CHI2<25 && Xi_TAU>0"; // Use Long Track for pion from Xi- decay
  XiL_Cut = "BachPi_TRACK_Type==3"; // Use Long Track for pion from Xi- decay
  XiD_Cut = "BachPi_TRACK_Type==5"; // Use DownStream Track for pion from Xi- decay

  Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<100) && Xib_MINIPCHI2<12 && Xib_TAU*1000>0.2";
  Xib_Cut = Xib_Cut && Xib_Trig && Jpsi_Cut && Lambda_Cut && (XiL_Cut||XiD_Cut) ;
  float ml = mLowXib;
  float mh = mHiXib;
  TString mCut = Form("Xib_DTF_M_JpsiXiLConstr>%f && Xib_DTF_M_JpsiXiLConstr<%f",ml,mh);
  Xib_Cut = Xib_Cut && mCut;
  MCXib_Cut = Xib_Cut && "Xib_BKGCAT<30";

  Lb_Cut = "(Lb_ENDVERTEX_CHI2/Lb_ENDVERTEX_NDOF < 10) && L_TAU>0 && Lb_MINIPCHI2 < 12 && Lb_TAU*1000>0.2";
  Lambdab_Cut = Jpsi_Cut && Lambda_Cut && Lb_Cut && Lb_Trig;

  ml = mLowLb;
  mh = mHiLb;
  mCut = Form("Lb_DTF_M_JpsiLConstr>%f&& Lb_DTF_M_JpsiLConstr<%f",ml,mh);
  Lambdab_Cut = Lambdab_Cut && mCut;

  MCLambdab_Cut = Lambdab_Cut&&"Lb_BKGCAT<30";
}
