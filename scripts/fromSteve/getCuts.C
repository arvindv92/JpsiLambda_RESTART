TCut Jpsi_Cut;
TCut Lambda_Cut;
TCut LambdaLL_Cut;
TCut LambdaDD_Cut;
TCut Xi_Cut;
TCut XiL_Cut;
TCut XiD_Cut;

TCut MCXib_Cut;
TCut MCXib_Cut_NoWin;
TCut Xib_Cut;
TCut Lb_Cut;
TCut Lambdab_Cut;
TCut MCLambdab_Cut;
TCut LambdabLL_Cut;
TCut MCLambdabLL_Cut;
TCut MCLambdabLL_Cut_NoWin;
TCut LambdabDD_Cut;
TCut MCLambdabDD_Cut;
TCut MCLambdabDD_Cut_NoWin;

float LbMassWin = 8.0;
bool applyXiPiTrackingAcceptanceCut = false;

double mLowXib= 5640; // def fit range 
double mHiXib = 5950; // def fit range 
double mLowLb=5500;
double mHiLb = 5750;

int iPolarity = 0;

void getCuts(){

  Lb_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Jpsi_L0MuonDecision_TOS>0)";
  Xib_Trig = "(Jpsi_L0DiMuonDecision_TOS>0 || Jpsi_L0MuonDecision_TOS>0)";
  // default mass window on Xi = 11  MeV

  Jpsi_Cut = "abs(Jpsi_MM-3098)<40";// && Jpsi_L0DiMuonDecision_TOS>0";
  Jpsi_Cut = Jpsi_Cut && "muplus_IPCHI2_OWNPV>4&&muminus_IPCHI2_OWNPV>4&&muplus_PT>550&&muminus_PT>550&&Jpsi_ENDVERTEX_CHI2<16 && Jpsi_FDCHI2_OWNPV>12";
  Lambda_Cut  = "L_ENDVERTEX_CHI2/L_ENDVERTEX_NDOF<12 && p_PT>500 && pi_PT>100 && L_FD_OWNPV/L_ENDVERTEX_ZERR>5 && L_TAU>0";
  TString lamWin = Form("abs(L_M-1115.7)<%f",LbMassWin);
  TCut cutLamMass = lamWin;
  Lambda_Cut = Lambda_Cut && cutLamMass;
  if(chargeCut>0) Lambda_Cut = Lambda_Cut && "L_ID>0";
  if(chargeCut<0) Lambda_Cut = Lambda_Cut && "L_ID<0";

  
  LambdaLL_Cut  = Lambda_Cut && "p_TRACK_Type==3 && pi_TRACK_Type==3";
  LambdaDD_Cut  = Lambda_Cut && "p_TRACK_Type==5 && pi_TRACK_Type==5";
  
  Xi_Cut = "abs(Xi_M-L_M+1115.7-1321.9)<10 && Xi_DIRA_ORIVX>0.995 && Xi_ENDVERTEX_CHI2<25 && Xi_TAU>0"; // Use Long Track for pion from Xi- decay
  
  XiL_Cut = "BachPi_TRACK_Type==3&&BachPi_IPCHI2_OWNPV>16"; // Use Long Track for pion from Xi- decay
  XiD_Cut = "BachPi_TRACK_Type==5";                         // Use DownStream Track for pion from Xi- decay

  // Apply p and eta cut, for comparison to region where tracking corrections are valid
  if(applyXiPiTrackingAcceptanceCut) XiL_Cut = XiL_Cut && "BachPi_P>5000 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))>1.9 && -log(tan(0.5*asin(BachPi_PT/BachPi_P)))<4.9";  

  // Build up Xib cuts
  Xib_Cut = "(Xib_ENDVERTEX_CHI2/Xib_ENDVERTEX_NDOF<10) && Xib_IPCHI2_OWNPV<12 && Xib_TAU*1000>0.2";
  if(iPolarity<0) Xib_Cut = Xib_Cut && "Polarity<0";
  if(iPolarity>0) Xib_Cut = Xib_Cut && "Polarity>0";
  if(applyEtaCutToHb) Xib_Cut = Xib_Cut && "(-log(tan(0.5*asin(Xib_PT/Xib_P)))>2 && -log(tan(0.5*asin(Xib_PT/Xib_P)))<6)";
  Xib_Cut = Xib_Cut && Xib_Trig && Jpsi_Cut && LambdaDD_Cut && Xi_Cut && (XiL_Cut||XiD_Cut) ;
  float ml = mLowXib;
  float mh = mHiXib;
  TString mCut = Form("Xib_DTF_M_JpsiXiLConstr>%f && Xib_DTF_M_JpsiXiLConstr<%f",ml,mh);
  Xib_Cut = Xib_Cut && mCut;

  // Build up Lb cuts
  Lb_Cut = "(Lb_ENDVERTEX_CHI2/Lb_ENDVERTEX_NDOF < 10) && L_TAU>0 && Lb_MINIPCHI2 < 12 && Lb_TAU*1000>0.2";
  if(applyEtaCutToHb) Lb_Cut = Lb_Cut && "(-log(tan(0.5*asin(Lb_PT/Lb_P)))>2 && -log(tan(0.5*asin(Lb_PT/Lb_P)))<6)";
  if(iPolarity<0) Lb_Cut = Lb_Cut && "Polarity<0";
  if(iPolarity>0) Lb_Cut = Lb_Cut && "Polarity>0";  
  LambdabDD_Cut = Jpsi_Cut && LambdaDD_Cut && Lb_Cut && Lb_Trig;
  LambdabLL_Cut = Jpsi_Cut && LambdaLL_Cut && Lb_Cut && Lb_Trig;
  ml = mLowLb;
  mh = mHiLb;
  mCut = Form("Lb_DTF_M_JpsiLConstr>%f&& Lb_DTF_M_JpsiLConstr<%f",ml,mh);
  LambdabDD_Cut = LambdabDD_Cut && mCut;
  LambdabLL_Cut = LambdabLL_Cut && mCut;

  // MC version of the cuts
  MCXib_Cut_NoWin = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<300)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<0))";
  MCXib_Cut = Xib_Cut && "(Xib_BKGCAT<20||(Xib_BKGCAT==50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25)||(Xib_BKGCAT==60&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<25))";

  MCLambdabDD_Cut_NoWin = LambdabDD_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<300)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<0))";
  MCLambdabLL_Cut_NoWin = LambdabLL_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<300)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<0))";

  MCLambdabDD_Cut = LambdabDD_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";
  MCLambdabLL_Cut = LambdabLL_Cut&&"(Lb_BKGCAT<20||(Lb_BKGCAT==60&&abs(Lb_DTF_M_JpsiLConstr-5620)<25)||(Lb_BKGCAT==50&&abs(Lb_DTF_M_JpsiLConstr-5620)<25))";


}
