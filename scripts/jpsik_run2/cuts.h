//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct 10 19:07:31 2018 by ROOT version 6.12/06
// from TChain B2JpsiKTree/MyTuple/
//////////////////////////////////////////////////////////

#ifndef cuts_h
#define cuts_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class cuts {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxB_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxB_OWNPV_COV = 1;
   static constexpr Int_t kMaxB_TOPPV_COV = 1;
   static constexpr Int_t kMaxJpsi_ENDVERTEX_COV = 1;
   static constexpr Int_t kMaxJpsi_OWNPV_COV = 1;
   static constexpr Int_t kMaxJpsi_TOPPV_COV = 1;
   static constexpr Int_t kMaxJpsi_ORIVX_COV = 1;
   static constexpr Int_t kMaxmuplus_OWNPV_COV = 1;
   static constexpr Int_t kMaxmuplus_TOPPV_COV = 1;
   static constexpr Int_t kMaxmuplus_ORIVX_COV = 1;
   static constexpr Int_t kMaxmuminus_OWNPV_COV = 1;
   static constexpr Int_t kMaxmuminus_TOPPV_COV = 1;
   static constexpr Int_t kMaxmuminus_ORIVX_COV = 1;
   static constexpr Int_t kMaxK_OWNPV_COV = 1;
   static constexpr Int_t kMaxK_TOPPV_COV = 1;
   static constexpr Int_t kMaxK_ORIVX_COV = 1;

   // Declaration of leaf types
   Double_t        B_MINIP;
   Double_t        B_MINIPCHI2;
   Double_t        B_MINIPNEXTBEST;
   Double_t        B_MINIPCHI2NEXTBEST;
   Double_t        B_ENDVERTEX_X;
   Double_t        B_ENDVERTEX_Y;
   Double_t        B_ENDVERTEX_Z;
   Double_t        B_ENDVERTEX_XERR;
   Double_t        B_ENDVERTEX_YERR;
   Double_t        B_ENDVERTEX_ZERR;
   Double_t        B_ENDVERTEX_CHI2;
   Int_t           B_ENDVERTEX_NDOF;
   Float_t         B_ENDVERTEX_COV_[3][3];
   Double_t        B_OWNPV_X;
   Double_t        B_OWNPV_Y;
   Double_t        B_OWNPV_Z;
   Double_t        B_OWNPV_XERR;
   Double_t        B_OWNPV_YERR;
   Double_t        B_OWNPV_ZERR;
   Double_t        B_OWNPV_CHI2;
   Int_t           B_OWNPV_NDOF;
   Float_t         B_OWNPV_COV_[3][3];
   Double_t        B_IP_OWNPV;
   Double_t        B_IPCHI2_OWNPV;
   Double_t        B_FD_OWNPV;
   Double_t        B_FDCHI2_OWNPV;
   Double_t        B_DIRA_OWNPV;
   Double_t        B_TOPPV_X;
   Double_t        B_TOPPV_Y;
   Double_t        B_TOPPV_Z;
   Double_t        B_TOPPV_XERR;
   Double_t        B_TOPPV_YERR;
   Double_t        B_TOPPV_ZERR;
   Double_t        B_TOPPV_CHI2;
   Int_t           B_TOPPV_NDOF;
   Float_t         B_TOPPV_COV_[3][3];
   Double_t        B_IP_TOPPV;
   Double_t        B_IPCHI2_TOPPV;
   Double_t        B_FD_TOPPV;
   Double_t        B_FDCHI2_TOPPV;
   Double_t        B_DIRA_TOPPV;
   Double_t        B_P;
   Double_t        B_PT;
   Double_t        B_PE;
   Double_t        B_PX;
   Double_t        B_PY;
   Double_t        B_PZ;
   Double_t        B_MM;
   Double_t        B_MMERR;
   Double_t        B_M;
   Int_t           B_ID;
   Double_t        B_TAU;
   Double_t        B_TAUERR;
   Double_t        B_TAUCHI2;
   Bool_t          B_L0Global_Dec;
   Bool_t          B_L0Global_TIS;
   Bool_t          B_L0Global_TOS;
   Bool_t          B_Hlt1Global_Dec;
   Bool_t          B_Hlt1Global_TIS;
   Bool_t          B_Hlt1Global_TOS;
   Bool_t          B_Hlt1Phys_Dec;
   Bool_t          B_Hlt1Phys_TIS;
   Bool_t          B_Hlt1Phys_TOS;
   Bool_t          B_Hlt2Global_Dec;
   Bool_t          B_Hlt2Global_TIS;
   Bool_t          B_Hlt2Global_TOS;
   Bool_t          B_Hlt2Phys_Dec;
   Bool_t          B_Hlt2Phys_TIS;
   Bool_t          B_Hlt2Phys_TOS;
   Double_t        B_ETA;
   Double_t        B_PHI;
   Double_t        B_DTF_CTAU;
   Double_t        B_DTF_CTAUERR;
   Double_t        B_DTF_CTAUS;
   Double_t        B_DTF_M_JpsiConstr;
   Double_t        B_DTF_VCHI2NDOF_JpsiConstr;
   Int_t           B_ConsBPV_nPV;
   Float_t         B_ConsBPV_J_psi_1S_M[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_MERR[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_P[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_PERR[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_ctau[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_ctauErr[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_decayLength[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_decayLengthErr[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_0_ID[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_0_PE[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_0_PX[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_0_PY[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_0_PZ[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_ID[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_PE[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_PX[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_PY[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_J_psi_1S_muminus_PZ[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_Kplus_ID[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_Kplus_PE[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_Kplus_PX[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_Kplus_PY[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_Kplus_PZ[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_M[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_MERR[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_P[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_PERR[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_chi2[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_nDOF[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_nIter[100];   //[B_ConsBPV_nPV]
   Float_t         B_ConsBPV_status[100];   //[B_ConsBPV_nPV]
   Bool_t          B_L0MuonDecision_Dec;
   Bool_t          B_L0MuonDecision_TIS;
   Bool_t          B_L0MuonDecision_TOS;
   Bool_t          B_L0DiMuonDecision_Dec;
   Bool_t          B_L0DiMuonDecision_TIS;
   Bool_t          B_L0DiMuonDecision_TOS;
   Bool_t          B_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          B_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          B_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          B_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          B_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          B_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          B_Hlt1TrackAllL0Decision_Dec;
   Bool_t          B_Hlt1TrackAllL0Decision_TIS;
   Bool_t          B_Hlt1TrackAllL0Decision_TOS;
   Bool_t          B_Hlt1TrackMuonDecision_Dec;
   Bool_t          B_Hlt1TrackMuonDecision_TIS;
   Bool_t          B_Hlt1TrackMuonDecision_TOS;
   Bool_t          B_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          B_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          B_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          B_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          B_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          B_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Double_t        Jpsi_CosTheta;
   Double_t        Jpsi_MINIP;
   Double_t        Jpsi_MINIPCHI2;
   Double_t        Jpsi_MINIPNEXTBEST;
   Double_t        Jpsi_MINIPCHI2NEXTBEST;
   Double_t        Jpsi_ENDVERTEX_X;
   Double_t        Jpsi_ENDVERTEX_Y;
   Double_t        Jpsi_ENDVERTEX_Z;
   Double_t        Jpsi_ENDVERTEX_XERR;
   Double_t        Jpsi_ENDVERTEX_YERR;
   Double_t        Jpsi_ENDVERTEX_ZERR;
   Double_t        Jpsi_ENDVERTEX_CHI2;
   Int_t           Jpsi_ENDVERTEX_NDOF;
   Float_t         Jpsi_ENDVERTEX_COV_[3][3];
   Double_t        Jpsi_OWNPV_X;
   Double_t        Jpsi_OWNPV_Y;
   Double_t        Jpsi_OWNPV_Z;
   Double_t        Jpsi_OWNPV_XERR;
   Double_t        Jpsi_OWNPV_YERR;
   Double_t        Jpsi_OWNPV_ZERR;
   Double_t        Jpsi_OWNPV_CHI2;
   Int_t           Jpsi_OWNPV_NDOF;
   Float_t         Jpsi_OWNPV_COV_[3][3];
   Double_t        Jpsi_IP_OWNPV;
   Double_t        Jpsi_IPCHI2_OWNPV;
   Double_t        Jpsi_FD_OWNPV;
   Double_t        Jpsi_FDCHI2_OWNPV;
   Double_t        Jpsi_DIRA_OWNPV;
   Double_t        Jpsi_TOPPV_X;
   Double_t        Jpsi_TOPPV_Y;
   Double_t        Jpsi_TOPPV_Z;
   Double_t        Jpsi_TOPPV_XERR;
   Double_t        Jpsi_TOPPV_YERR;
   Double_t        Jpsi_TOPPV_ZERR;
   Double_t        Jpsi_TOPPV_CHI2;
   Int_t           Jpsi_TOPPV_NDOF;
   Float_t         Jpsi_TOPPV_COV_[3][3];
   Double_t        Jpsi_IP_TOPPV;
   Double_t        Jpsi_IPCHI2_TOPPV;
   Double_t        Jpsi_FD_TOPPV;
   Double_t        Jpsi_FDCHI2_TOPPV;
   Double_t        Jpsi_DIRA_TOPPV;
   Double_t        Jpsi_ORIVX_X;
   Double_t        Jpsi_ORIVX_Y;
   Double_t        Jpsi_ORIVX_Z;
   Double_t        Jpsi_ORIVX_XERR;
   Double_t        Jpsi_ORIVX_YERR;
   Double_t        Jpsi_ORIVX_ZERR;
   Double_t        Jpsi_ORIVX_CHI2;
   Int_t           Jpsi_ORIVX_NDOF;
   Float_t         Jpsi_ORIVX_COV_[3][3];
   Double_t        Jpsi_IP_ORIVX;
   Double_t        Jpsi_IPCHI2_ORIVX;
   Double_t        Jpsi_FD_ORIVX;
   Double_t        Jpsi_FDCHI2_ORIVX;
   Double_t        Jpsi_DIRA_ORIVX;
   Double_t        Jpsi_P;
   Double_t        Jpsi_PT;
   Double_t        Jpsi_PE;
   Double_t        Jpsi_PX;
   Double_t        Jpsi_PY;
   Double_t        Jpsi_PZ;
   Double_t        Jpsi_MM;
   Double_t        Jpsi_MMERR;
   Double_t        Jpsi_M;
   Int_t           Jpsi_ID;
   Double_t        Jpsi_TAU;
   Double_t        Jpsi_TAUERR;
   Double_t        Jpsi_TAUCHI2;
   Bool_t          Jpsi_L0Global_Dec;
   Bool_t          Jpsi_L0Global_TIS;
   Bool_t          Jpsi_L0Global_TOS;
   Bool_t          Jpsi_Hlt1Global_Dec;
   Bool_t          Jpsi_Hlt1Global_TIS;
   Bool_t          Jpsi_Hlt1Global_TOS;
   Bool_t          Jpsi_Hlt1Phys_Dec;
   Bool_t          Jpsi_Hlt1Phys_TIS;
   Bool_t          Jpsi_Hlt1Phys_TOS;
   Bool_t          Jpsi_Hlt2Global_Dec;
   Bool_t          Jpsi_Hlt2Global_TIS;
   Bool_t          Jpsi_Hlt2Global_TOS;
   Bool_t          Jpsi_Hlt2Phys_Dec;
   Bool_t          Jpsi_Hlt2Phys_TIS;
   Bool_t          Jpsi_Hlt2Phys_TOS;
   Double_t        Jpsi_ETA;
   Double_t        Jpsi_PHI;
   Bool_t          Jpsi_L0MuonDecision_Dec;
   Bool_t          Jpsi_L0MuonDecision_TIS;
   Bool_t          Jpsi_L0MuonDecision_TOS;
   Bool_t          Jpsi_L0DiMuonDecision_Dec;
   Bool_t          Jpsi_L0DiMuonDecision_TIS;
   Bool_t          Jpsi_L0DiMuonDecision_TOS;
   Bool_t          Jpsi_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          Jpsi_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          Jpsi_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          Jpsi_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          Jpsi_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          Jpsi_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          Jpsi_Hlt1TrackAllL0Decision_Dec;
   Bool_t          Jpsi_Hlt1TrackAllL0Decision_TIS;
   Bool_t          Jpsi_Hlt1TrackAllL0Decision_TOS;
   Bool_t          Jpsi_Hlt1TrackMuonDecision_Dec;
   Bool_t          Jpsi_Hlt1TrackMuonDecision_TIS;
   Bool_t          Jpsi_Hlt1TrackMuonDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Double_t        muplus_CosTheta;
   Double_t        muplus_MINIP;
   Double_t        muplus_MINIPCHI2;
   Double_t        muplus_MINIPNEXTBEST;
   Double_t        muplus_MINIPCHI2NEXTBEST;
   Double_t        muplus_OWNPV_X;
   Double_t        muplus_OWNPV_Y;
   Double_t        muplus_OWNPV_Z;
   Double_t        muplus_OWNPV_XERR;
   Double_t        muplus_OWNPV_YERR;
   Double_t        muplus_OWNPV_ZERR;
   Double_t        muplus_OWNPV_CHI2;
   Int_t           muplus_OWNPV_NDOF;
   Float_t         muplus_OWNPV_COV_[3][3];
   Double_t        muplus_IP_OWNPV;
   Double_t        muplus_IPCHI2_OWNPV;
   Double_t        muplus_TOPPV_X;
   Double_t        muplus_TOPPV_Y;
   Double_t        muplus_TOPPV_Z;
   Double_t        muplus_TOPPV_XERR;
   Double_t        muplus_TOPPV_YERR;
   Double_t        muplus_TOPPV_ZERR;
   Double_t        muplus_TOPPV_CHI2;
   Int_t           muplus_TOPPV_NDOF;
   Float_t         muplus_TOPPV_COV_[3][3];
   Double_t        muplus_IP_TOPPV;
   Double_t        muplus_IPCHI2_TOPPV;
   Double_t        muplus_ORIVX_X;
   Double_t        muplus_ORIVX_Y;
   Double_t        muplus_ORIVX_Z;
   Double_t        muplus_ORIVX_XERR;
   Double_t        muplus_ORIVX_YERR;
   Double_t        muplus_ORIVX_ZERR;
   Double_t        muplus_ORIVX_CHI2;
   Int_t           muplus_ORIVX_NDOF;
   Float_t         muplus_ORIVX_COV_[3][3];
   Double_t        muplus_IP_ORIVX;
   Double_t        muplus_IPCHI2_ORIVX;
   Double_t        muplus_P;
   Double_t        muplus_PT;
   Double_t        muplus_PE;
   Double_t        muplus_PX;
   Double_t        muplus_PY;
   Double_t        muplus_PZ;
   Double_t        muplus_M;
   Int_t           muplus_ID;
   Double_t        muplus_PIDe;
   Double_t        muplus_PIDmu;
   Double_t        muplus_PIDK;
   Double_t        muplus_PIDp;
   Double_t        muplus_ProbNNe;
   Double_t        muplus_ProbNNk;
   Double_t        muplus_ProbNNp;
   Double_t        muplus_ProbNNpi;
   Double_t        muplus_ProbNNmu;
   Double_t        muplus_ProbNNghost;
   Bool_t          muplus_hasMuon;
   Bool_t          muplus_isMuon;
   Bool_t          muplus_hasRich;
   Bool_t          muplus_UsedRichAerogel;
   Bool_t          muplus_UsedRich1Gas;
   Bool_t          muplus_UsedRich2Gas;
   Bool_t          muplus_RichAboveElThres;
   Bool_t          muplus_RichAboveMuThres;
   Bool_t          muplus_RichAbovePiThres;
   Bool_t          muplus_RichAboveKaThres;
   Bool_t          muplus_RichAbovePrThres;
   Bool_t          muplus_hasCalo;
   Bool_t          muplus_L0Global_Dec;
   Bool_t          muplus_L0Global_TIS;
   Bool_t          muplus_L0Global_TOS;
   Bool_t          muplus_Hlt1Global_Dec;
   Bool_t          muplus_Hlt1Global_TIS;
   Bool_t          muplus_Hlt1Global_TOS;
   Bool_t          muplus_Hlt1Phys_Dec;
   Bool_t          muplus_Hlt1Phys_TIS;
   Bool_t          muplus_Hlt1Phys_TOS;
   Bool_t          muplus_Hlt2Global_Dec;
   Bool_t          muplus_Hlt2Global_TIS;
   Bool_t          muplus_Hlt2Global_TOS;
   Bool_t          muplus_Hlt2Phys_Dec;
   Bool_t          muplus_Hlt2Phys_TIS;
   Bool_t          muplus_Hlt2Phys_TOS;
   Int_t           muplus_TRACK_Type;
   Int_t           muplus_TRACK_Key;
   Double_t        muplus_TRACK_CHI2NDOF;
   Double_t        muplus_TRACK_PCHI2;
   Double_t        muplus_TRACK_MatchCHI2;
   Double_t        muplus_TRACK_GhostProb;
   Double_t        muplus_TRACK_CloneDist;
   Double_t        muplus_TRACK_Likelihood;
   Double_t        muminus_CosTheta;
   Double_t        muminus_MINIP;
   Double_t        muminus_MINIPCHI2;
   Double_t        muminus_MINIPNEXTBEST;
   Double_t        muminus_MINIPCHI2NEXTBEST;
   Double_t        muminus_OWNPV_X;
   Double_t        muminus_OWNPV_Y;
   Double_t        muminus_OWNPV_Z;
   Double_t        muminus_OWNPV_XERR;
   Double_t        muminus_OWNPV_YERR;
   Double_t        muminus_OWNPV_ZERR;
   Double_t        muminus_OWNPV_CHI2;
   Int_t           muminus_OWNPV_NDOF;
   Float_t         muminus_OWNPV_COV_[3][3];
   Double_t        muminus_IP_OWNPV;
   Double_t        muminus_IPCHI2_OWNPV;
   Double_t        muminus_TOPPV_X;
   Double_t        muminus_TOPPV_Y;
   Double_t        muminus_TOPPV_Z;
   Double_t        muminus_TOPPV_XERR;
   Double_t        muminus_TOPPV_YERR;
   Double_t        muminus_TOPPV_ZERR;
   Double_t        muminus_TOPPV_CHI2;
   Int_t           muminus_TOPPV_NDOF;
   Float_t         muminus_TOPPV_COV_[3][3];
   Double_t        muminus_IP_TOPPV;
   Double_t        muminus_IPCHI2_TOPPV;
   Double_t        muminus_ORIVX_X;
   Double_t        muminus_ORIVX_Y;
   Double_t        muminus_ORIVX_Z;
   Double_t        muminus_ORIVX_XERR;
   Double_t        muminus_ORIVX_YERR;
   Double_t        muminus_ORIVX_ZERR;
   Double_t        muminus_ORIVX_CHI2;
   Int_t           muminus_ORIVX_NDOF;
   Float_t         muminus_ORIVX_COV_[3][3];
   Double_t        muminus_IP_ORIVX;
   Double_t        muminus_IPCHI2_ORIVX;
   Double_t        muminus_P;
   Double_t        muminus_PT;
   Double_t        muminus_PE;
   Double_t        muminus_PX;
   Double_t        muminus_PY;
   Double_t        muminus_PZ;
   Double_t        muminus_M;
   Int_t           muminus_ID;
   Double_t        muminus_PIDe;
   Double_t        muminus_PIDmu;
   Double_t        muminus_PIDK;
   Double_t        muminus_PIDp;
   Double_t        muminus_ProbNNe;
   Double_t        muminus_ProbNNk;
   Double_t        muminus_ProbNNp;
   Double_t        muminus_ProbNNpi;
   Double_t        muminus_ProbNNmu;
   Double_t        muminus_ProbNNghost;
   Bool_t          muminus_hasMuon;
   Bool_t          muminus_isMuon;
   Bool_t          muminus_hasRich;
   Bool_t          muminus_UsedRichAerogel;
   Bool_t          muminus_UsedRich1Gas;
   Bool_t          muminus_UsedRich2Gas;
   Bool_t          muminus_RichAboveElThres;
   Bool_t          muminus_RichAboveMuThres;
   Bool_t          muminus_RichAbovePiThres;
   Bool_t          muminus_RichAboveKaThres;
   Bool_t          muminus_RichAbovePrThres;
   Bool_t          muminus_hasCalo;
   Bool_t          muminus_L0Global_Dec;
   Bool_t          muminus_L0Global_TIS;
   Bool_t          muminus_L0Global_TOS;
   Bool_t          muminus_Hlt1Global_Dec;
   Bool_t          muminus_Hlt1Global_TIS;
   Bool_t          muminus_Hlt1Global_TOS;
   Bool_t          muminus_Hlt1Phys_Dec;
   Bool_t          muminus_Hlt1Phys_TIS;
   Bool_t          muminus_Hlt1Phys_TOS;
   Bool_t          muminus_Hlt2Global_Dec;
   Bool_t          muminus_Hlt2Global_TIS;
   Bool_t          muminus_Hlt2Global_TOS;
   Bool_t          muminus_Hlt2Phys_Dec;
   Bool_t          muminus_Hlt2Phys_TIS;
   Bool_t          muminus_Hlt2Phys_TOS;
   Int_t           muminus_TRACK_Type;
   Int_t           muminus_TRACK_Key;
   Double_t        muminus_TRACK_CHI2NDOF;
   Double_t        muminus_TRACK_PCHI2;
   Double_t        muminus_TRACK_MatchCHI2;
   Double_t        muminus_TRACK_GhostProb;
   Double_t        muminus_TRACK_CloneDist;
   Double_t        muminus_TRACK_Likelihood;
   Double_t        K_CosTheta;
   Double_t        K_MINIP;
   Double_t        K_MINIPCHI2;
   Double_t        K_MINIPNEXTBEST;
   Double_t        K_MINIPCHI2NEXTBEST;
   Double_t        K_OWNPV_X;
   Double_t        K_OWNPV_Y;
   Double_t        K_OWNPV_Z;
   Double_t        K_OWNPV_XERR;
   Double_t        K_OWNPV_YERR;
   Double_t        K_OWNPV_ZERR;
   Double_t        K_OWNPV_CHI2;
   Int_t           K_OWNPV_NDOF;
   Float_t         K_OWNPV_COV_[3][3];
   Double_t        K_IP_OWNPV;
   Double_t        K_IPCHI2_OWNPV;
   Double_t        K_TOPPV_X;
   Double_t        K_TOPPV_Y;
   Double_t        K_TOPPV_Z;
   Double_t        K_TOPPV_XERR;
   Double_t        K_TOPPV_YERR;
   Double_t        K_TOPPV_ZERR;
   Double_t        K_TOPPV_CHI2;
   Int_t           K_TOPPV_NDOF;
   Float_t         K_TOPPV_COV_[3][3];
   Double_t        K_IP_TOPPV;
   Double_t        K_IPCHI2_TOPPV;
   Double_t        K_ORIVX_X;
   Double_t        K_ORIVX_Y;
   Double_t        K_ORIVX_Z;
   Double_t        K_ORIVX_XERR;
   Double_t        K_ORIVX_YERR;
   Double_t        K_ORIVX_ZERR;
   Double_t        K_ORIVX_CHI2;
   Int_t           K_ORIVX_NDOF;
   Float_t         K_ORIVX_COV_[3][3];
   Double_t        K_IP_ORIVX;
   Double_t        K_IPCHI2_ORIVX;
   Double_t        K_P;
   Double_t        K_PT;
   Double_t        K_PE;
   Double_t        K_PX;
   Double_t        K_PY;
   Double_t        K_PZ;
   Double_t        K_M;
   Int_t           K_ID;
   Double_t        K_PIDe;
   Double_t        K_PIDmu;
   Double_t        K_PIDK;
   Double_t        K_PIDp;
   Double_t        K_ProbNNe;
   Double_t        K_ProbNNk;
   Double_t        K_ProbNNp;
   Double_t        K_ProbNNpi;
   Double_t        K_ProbNNmu;
   Double_t        K_ProbNNghost;
   Bool_t          K_hasMuon;
   Bool_t          K_isMuon;
   Bool_t          K_hasRich;
   Bool_t          K_UsedRichAerogel;
   Bool_t          K_UsedRich1Gas;
   Bool_t          K_UsedRich2Gas;
   Bool_t          K_RichAboveElThres;
   Bool_t          K_RichAboveMuThres;
   Bool_t          K_RichAbovePiThres;
   Bool_t          K_RichAboveKaThres;
   Bool_t          K_RichAbovePrThres;
   Bool_t          K_hasCalo;
   Bool_t          K_L0Global_Dec;
   Bool_t          K_L0Global_TIS;
   Bool_t          K_L0Global_TOS;
   Bool_t          K_Hlt1Global_Dec;
   Bool_t          K_Hlt1Global_TIS;
   Bool_t          K_Hlt1Global_TOS;
   Bool_t          K_Hlt1Phys_Dec;
   Bool_t          K_Hlt1Phys_TIS;
   Bool_t          K_Hlt1Phys_TOS;
   Bool_t          K_Hlt2Global_Dec;
   Bool_t          K_Hlt2Global_TIS;
   Bool_t          K_Hlt2Global_TOS;
   Bool_t          K_Hlt2Phys_Dec;
   Bool_t          K_Hlt2Phys_TIS;
   Bool_t          K_Hlt2Phys_TOS;
   Int_t           K_TRACK_Type;
   Int_t           K_TRACK_Key;
   Double_t        K_TRACK_CHI2NDOF;
   Double_t        K_TRACK_PCHI2;
   Double_t        K_TRACK_MatchCHI2;
   Double_t        K_TRACK_GhostProb;
   Double_t        K_TRACK_CloneDist;
   Double_t        K_TRACK_Likelihood;
   Double_t        K_ETA;
   Double_t        K_PHI;
   UInt_t          nCandidate;
   ULong64_t       totCandidates;
   ULong64_t       EventInSequence;
   Double_t        nlong;
   Double_t        ntracks;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          BCID;
   Int_t           BCType;
   UInt_t          OdinTCK;
   UInt_t          L0DUTCK;
   UInt_t          HLT1TCK;
   UInt_t          HLT2TCK;
   ULong64_t       GpsTime;
   Short_t         Polarity;
   Int_t           nPV;
   Float_t         PVX[100];   //[nPV]
   Float_t         PVY[100];   //[nPV]
   Float_t         PVZ[100];   //[nPV]
   Float_t         PVXERR[100];   //[nPV]
   Float_t         PVYERR[100];   //[nPV]
   Float_t         PVZERR[100];   //[nPV]
   Float_t         PVCHI2[100];   //[nPV]
   Float_t         PVNDOF[100];   //[nPV]
   Float_t         PVNTRACKS[100];   //[nPV]
   UInt_t          StrippingFullDSTDiMuonJpsi2MuMuDetachedLineDecision;

   // List of branches
   TBranch        *b_B_MINIP;   //!
   TBranch        *b_B_MINIPCHI2;   //!
   TBranch        *b_B_MINIPNEXTBEST;   //!
   TBranch        *b_B_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_B_ENDVERTEX_X;   //!
   TBranch        *b_B_ENDVERTEX_Y;   //!
   TBranch        *b_B_ENDVERTEX_Z;   //!
   TBranch        *b_B_ENDVERTEX_XERR;   //!
   TBranch        *b_B_ENDVERTEX_YERR;   //!
   TBranch        *b_B_ENDVERTEX_ZERR;   //!
   TBranch        *b_B_ENDVERTEX_CHI2;   //!
   TBranch        *b_B_ENDVERTEX_NDOF;   //!
   TBranch        *b_B_ENDVERTEX_COV_;   //!
   TBranch        *b_B_OWNPV_X;   //!
   TBranch        *b_B_OWNPV_Y;   //!
   TBranch        *b_B_OWNPV_Z;   //!
   TBranch        *b_B_OWNPV_XERR;   //!
   TBranch        *b_B_OWNPV_YERR;   //!
   TBranch        *b_B_OWNPV_ZERR;   //!
   TBranch        *b_B_OWNPV_CHI2;   //!
   TBranch        *b_B_OWNPV_NDOF;   //!
   TBranch        *b_B_OWNPV_COV_;   //!
   TBranch        *b_B_IP_OWNPV;   //!
   TBranch        *b_B_IPCHI2_OWNPV;   //!
   TBranch        *b_B_FD_OWNPV;   //!
   TBranch        *b_B_FDCHI2_OWNPV;   //!
   TBranch        *b_B_DIRA_OWNPV;   //!
   TBranch        *b_B_TOPPV_X;   //!
   TBranch        *b_B_TOPPV_Y;   //!
   TBranch        *b_B_TOPPV_Z;   //!
   TBranch        *b_B_TOPPV_XERR;   //!
   TBranch        *b_B_TOPPV_YERR;   //!
   TBranch        *b_B_TOPPV_ZERR;   //!
   TBranch        *b_B_TOPPV_CHI2;   //!
   TBranch        *b_B_TOPPV_NDOF;   //!
   TBranch        *b_B_TOPPV_COV_;   //!
   TBranch        *b_B_IP_TOPPV;   //!
   TBranch        *b_B_IPCHI2_TOPPV;   //!
   TBranch        *b_B_FD_TOPPV;   //!
   TBranch        *b_B_FDCHI2_TOPPV;   //!
   TBranch        *b_B_DIRA_TOPPV;   //!
   TBranch        *b_B_P;   //!
   TBranch        *b_B_PT;   //!
   TBranch        *b_B_PE;   //!
   TBranch        *b_B_PX;   //!
   TBranch        *b_B_PY;   //!
   TBranch        *b_B_PZ;   //!
   TBranch        *b_B_MM;   //!
   TBranch        *b_B_MMERR;   //!
   TBranch        *b_B_M;   //!
   TBranch        *b_B_ID;   //!
   TBranch        *b_B_TAU;   //!
   TBranch        *b_B_TAUERR;   //!
   TBranch        *b_B_TAUCHI2;   //!
   TBranch        *b_B_L0Global_Dec;   //!
   TBranch        *b_B_L0Global_TIS;   //!
   TBranch        *b_B_L0Global_TOS;   //!
   TBranch        *b_B_Hlt1Global_Dec;   //!
   TBranch        *b_B_Hlt1Global_TIS;   //!
   TBranch        *b_B_Hlt1Global_TOS;   //!
   TBranch        *b_B_Hlt1Phys_Dec;   //!
   TBranch        *b_B_Hlt1Phys_TIS;   //!
   TBranch        *b_B_Hlt1Phys_TOS;   //!
   TBranch        *b_B_Hlt2Global_Dec;   //!
   TBranch        *b_B_Hlt2Global_TIS;   //!
   TBranch        *b_B_Hlt2Global_TOS;   //!
   TBranch        *b_B_Hlt2Phys_Dec;   //!
   TBranch        *b_B_Hlt2Phys_TIS;   //!
   TBranch        *b_B_Hlt2Phys_TOS;   //!
   TBranch        *b_B_ETA;   //!
   TBranch        *b_B_PHI;   //!
   TBranch        *b_B_DTF_CTAU;   //!
   TBranch        *b_B_DTF_CTAUERR;   //!
   TBranch        *b_B_DTF_CTAUS;   //!
   TBranch        *b_B_DTF_M_JpsiConstr;   //!
   TBranch        *b_B_DTF_VCHI2NDOF_JpsiConstr;   //!
   TBranch        *b_B_ConsBPV_nPV;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_M;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_MERR;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_P;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_PERR;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_ctau;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_ctauErr;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_decayLength;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_decayLengthErr;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_0_ID;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_0_PE;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_0_PX;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_0_PY;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_0_PZ;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_ID;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_PE;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_PX;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_PY;   //!
   TBranch        *b_B_ConsBPV_J_psi_1S_muminus_PZ;   //!
   TBranch        *b_B_ConsBPV_Kplus_ID;   //!
   TBranch        *b_B_ConsBPV_Kplus_PE;   //!
   TBranch        *b_B_ConsBPV_Kplus_PX;   //!
   TBranch        *b_B_ConsBPV_Kplus_PY;   //!
   TBranch        *b_B_ConsBPV_Kplus_PZ;   //!
   TBranch        *b_B_ConsBPV_M;   //!
   TBranch        *b_B_ConsBPV_MERR;   //!
   TBranch        *b_B_ConsBPV_P;   //!
   TBranch        *b_B_ConsBPV_PERR;   //!
   TBranch        *b_B_ConsBPV_chi2;   //!
   TBranch        *b_B_ConsBPV_nDOF;   //!
   TBranch        *b_B_ConsBPV_nIter;   //!
   TBranch        *b_B_ConsBPV_status;   //!
   TBranch        *b_B_L0MuonDecision_Dec;   //!
   TBranch        *b_B_L0MuonDecision_TIS;   //!
   TBranch        *b_B_L0MuonDecision_TOS;   //!
   TBranch        *b_B_L0DiMuonDecision_Dec;   //!
   TBranch        *b_B_L0DiMuonDecision_TIS;   //!
   TBranch        *b_B_L0DiMuonDecision_TOS;   //!
   TBranch        *b_B_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_B_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_B_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_B_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_B_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_B_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_B_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_B_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_B_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_B_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_B_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_B_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_B_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_B_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_B_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_B_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_B_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_B_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_Jpsi_CosTheta;   //!
   TBranch        *b_Jpsi_MINIP;   //!
   TBranch        *b_Jpsi_MINIPCHI2;   //!
   TBranch        *b_Jpsi_MINIPNEXTBEST;   //!
   TBranch        *b_Jpsi_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_Jpsi_ENDVERTEX_X;   //!
   TBranch        *b_Jpsi_ENDVERTEX_Y;   //!
   TBranch        *b_Jpsi_ENDVERTEX_Z;   //!
   TBranch        *b_Jpsi_ENDVERTEX_XERR;   //!
   TBranch        *b_Jpsi_ENDVERTEX_YERR;   //!
   TBranch        *b_Jpsi_ENDVERTEX_ZERR;   //!
   TBranch        *b_Jpsi_ENDVERTEX_CHI2;   //!
   TBranch        *b_Jpsi_ENDVERTEX_NDOF;   //!
   TBranch        *b_Jpsi_ENDVERTEX_COV_;   //!
   TBranch        *b_Jpsi_OWNPV_X;   //!
   TBranch        *b_Jpsi_OWNPV_Y;   //!
   TBranch        *b_Jpsi_OWNPV_Z;   //!
   TBranch        *b_Jpsi_OWNPV_XERR;   //!
   TBranch        *b_Jpsi_OWNPV_YERR;   //!
   TBranch        *b_Jpsi_OWNPV_ZERR;   //!
   TBranch        *b_Jpsi_OWNPV_CHI2;   //!
   TBranch        *b_Jpsi_OWNPV_NDOF;   //!
   TBranch        *b_Jpsi_OWNPV_COV_;   //!
   TBranch        *b_Jpsi_IP_OWNPV;   //!
   TBranch        *b_Jpsi_IPCHI2_OWNPV;   //!
   TBranch        *b_Jpsi_FD_OWNPV;   //!
   TBranch        *b_Jpsi_FDCHI2_OWNPV;   //!
   TBranch        *b_Jpsi_DIRA_OWNPV;   //!
   TBranch        *b_Jpsi_TOPPV_X;   //!
   TBranch        *b_Jpsi_TOPPV_Y;   //!
   TBranch        *b_Jpsi_TOPPV_Z;   //!
   TBranch        *b_Jpsi_TOPPV_XERR;   //!
   TBranch        *b_Jpsi_TOPPV_YERR;   //!
   TBranch        *b_Jpsi_TOPPV_ZERR;   //!
   TBranch        *b_Jpsi_TOPPV_CHI2;   //!
   TBranch        *b_Jpsi_TOPPV_NDOF;   //!
   TBranch        *b_Jpsi_TOPPV_COV_;   //!
   TBranch        *b_Jpsi_IP_TOPPV;   //!
   TBranch        *b_Jpsi_IPCHI2_TOPPV;   //!
   TBranch        *b_Jpsi_FD_TOPPV;   //!
   TBranch        *b_Jpsi_FDCHI2_TOPPV;   //!
   TBranch        *b_Jpsi_DIRA_TOPPV;   //!
   TBranch        *b_Jpsi_ORIVX_X;   //!
   TBranch        *b_Jpsi_ORIVX_Y;   //!
   TBranch        *b_Jpsi_ORIVX_Z;   //!
   TBranch        *b_Jpsi_ORIVX_XERR;   //!
   TBranch        *b_Jpsi_ORIVX_YERR;   //!
   TBranch        *b_Jpsi_ORIVX_ZERR;   //!
   TBranch        *b_Jpsi_ORIVX_CHI2;   //!
   TBranch        *b_Jpsi_ORIVX_NDOF;   //!
   TBranch        *b_Jpsi_ORIVX_COV_;   //!
   TBranch        *b_Jpsi_IP_ORIVX;   //!
   TBranch        *b_Jpsi_IPCHI2_ORIVX;   //!
   TBranch        *b_Jpsi_FD_ORIVX;   //!
   TBranch        *b_Jpsi_FDCHI2_ORIVX;   //!
   TBranch        *b_Jpsi_DIRA_ORIVX;   //!
   TBranch        *b_Jpsi_P;   //!
   TBranch        *b_Jpsi_PT;   //!
   TBranch        *b_Jpsi_PE;   //!
   TBranch        *b_Jpsi_PX;   //!
   TBranch        *b_Jpsi_PY;   //!
   TBranch        *b_Jpsi_PZ;   //!
   TBranch        *b_Jpsi_MM;   //!
   TBranch        *b_Jpsi_MMERR;   //!
   TBranch        *b_Jpsi_M;   //!
   TBranch        *b_Jpsi_ID;   //!
   TBranch        *b_Jpsi_TAU;   //!
   TBranch        *b_Jpsi_TAUERR;   //!
   TBranch        *b_Jpsi_TAUCHI2;   //!
   TBranch        *b_Jpsi_L0Global_Dec;   //!
   TBranch        *b_Jpsi_L0Global_TIS;   //!
   TBranch        *b_Jpsi_L0Global_TOS;   //!
   TBranch        *b_Jpsi_Hlt1Global_Dec;   //!
   TBranch        *b_Jpsi_Hlt1Global_TIS;   //!
   TBranch        *b_Jpsi_Hlt1Global_TOS;   //!
   TBranch        *b_Jpsi_Hlt1Phys_Dec;   //!
   TBranch        *b_Jpsi_Hlt1Phys_TIS;   //!
   TBranch        *b_Jpsi_Hlt1Phys_TOS;   //!
   TBranch        *b_Jpsi_Hlt2Global_Dec;   //!
   TBranch        *b_Jpsi_Hlt2Global_TIS;   //!
   TBranch        *b_Jpsi_Hlt2Global_TOS;   //!
   TBranch        *b_Jpsi_Hlt2Phys_Dec;   //!
   TBranch        *b_Jpsi_Hlt2Phys_TIS;   //!
   TBranch        *b_Jpsi_Hlt2Phys_TOS;   //!
   TBranch        *b_Jpsi_ETA;   //!
   TBranch        *b_Jpsi_PHI;   //!
   TBranch        *b_Jpsi_L0MuonDecision_Dec;   //!
   TBranch        *b_Jpsi_L0MuonDecision_TIS;   //!
   TBranch        *b_Jpsi_L0MuonDecision_TOS;   //!
   TBranch        *b_Jpsi_L0DiMuonDecision_Dec;   //!
   TBranch        *b_Jpsi_L0DiMuonDecision_TIS;   //!
   TBranch        *b_Jpsi_L0DiMuonDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_Jpsi_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_Jpsi_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_Jpsi_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_muplus_CosTheta;   //!
   TBranch        *b_muplus_MINIP;   //!
   TBranch        *b_muplus_MINIPCHI2;   //!
   TBranch        *b_muplus_MINIPNEXTBEST;   //!
   TBranch        *b_muplus_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_muplus_OWNPV_X;   //!
   TBranch        *b_muplus_OWNPV_Y;   //!
   TBranch        *b_muplus_OWNPV_Z;   //!
   TBranch        *b_muplus_OWNPV_XERR;   //!
   TBranch        *b_muplus_OWNPV_YERR;   //!
   TBranch        *b_muplus_OWNPV_ZERR;   //!
   TBranch        *b_muplus_OWNPV_CHI2;   //!
   TBranch        *b_muplus_OWNPV_NDOF;   //!
   TBranch        *b_muplus_OWNPV_COV_;   //!
   TBranch        *b_muplus_IP_OWNPV;   //!
   TBranch        *b_muplus_IPCHI2_OWNPV;   //!
   TBranch        *b_muplus_TOPPV_X;   //!
   TBranch        *b_muplus_TOPPV_Y;   //!
   TBranch        *b_muplus_TOPPV_Z;   //!
   TBranch        *b_muplus_TOPPV_XERR;   //!
   TBranch        *b_muplus_TOPPV_YERR;   //!
   TBranch        *b_muplus_TOPPV_ZERR;   //!
   TBranch        *b_muplus_TOPPV_CHI2;   //!
   TBranch        *b_muplus_TOPPV_NDOF;   //!
   TBranch        *b_muplus_TOPPV_COV_;   //!
   TBranch        *b_muplus_IP_TOPPV;   //!
   TBranch        *b_muplus_IPCHI2_TOPPV;   //!
   TBranch        *b_muplus_ORIVX_X;   //!
   TBranch        *b_muplus_ORIVX_Y;   //!
   TBranch        *b_muplus_ORIVX_Z;   //!
   TBranch        *b_muplus_ORIVX_XERR;   //!
   TBranch        *b_muplus_ORIVX_YERR;   //!
   TBranch        *b_muplus_ORIVX_ZERR;   //!
   TBranch        *b_muplus_ORIVX_CHI2;   //!
   TBranch        *b_muplus_ORIVX_NDOF;   //!
   TBranch        *b_muplus_ORIVX_COV_;   //!
   TBranch        *b_muplus_IP_ORIVX;   //!
   TBranch        *b_muplus_IPCHI2_ORIVX;   //!
   TBranch        *b_muplus_P;   //!
   TBranch        *b_muplus_PT;   //!
   TBranch        *b_muplus_PE;   //!
   TBranch        *b_muplus_PX;   //!
   TBranch        *b_muplus_PY;   //!
   TBranch        *b_muplus_PZ;   //!
   TBranch        *b_muplus_M;   //!
   TBranch        *b_muplus_ID;   //!
   TBranch        *b_muplus_PIDe;   //!
   TBranch        *b_muplus_PIDmu;   //!
   TBranch        *b_muplus_PIDK;   //!
   TBranch        *b_muplus_PIDp;   //!
   TBranch        *b_muplus_ProbNNe;   //!
   TBranch        *b_muplus_ProbNNk;   //!
   TBranch        *b_muplus_ProbNNp;   //!
   TBranch        *b_muplus_ProbNNpi;   //!
   TBranch        *b_muplus_ProbNNmu;   //!
   TBranch        *b_muplus_ProbNNghost;   //!
   TBranch        *b_muplus_hasMuon;   //!
   TBranch        *b_muplus_isMuon;   //!
   TBranch        *b_muplus_hasRich;   //!
   TBranch        *b_muplus_UsedRichAerogel;   //!
   TBranch        *b_muplus_UsedRich1Gas;   //!
   TBranch        *b_muplus_UsedRich2Gas;   //!
   TBranch        *b_muplus_RichAboveElThres;   //!
   TBranch        *b_muplus_RichAboveMuThres;   //!
   TBranch        *b_muplus_RichAbovePiThres;   //!
   TBranch        *b_muplus_RichAboveKaThres;   //!
   TBranch        *b_muplus_RichAbovePrThres;   //!
   TBranch        *b_muplus_hasCalo;   //!
   TBranch        *b_muplus_L0Global_Dec;   //!
   TBranch        *b_muplus_L0Global_TIS;   //!
   TBranch        *b_muplus_L0Global_TOS;   //!
   TBranch        *b_muplus_Hlt1Global_Dec;   //!
   TBranch        *b_muplus_Hlt1Global_TIS;   //!
   TBranch        *b_muplus_Hlt1Global_TOS;   //!
   TBranch        *b_muplus_Hlt1Phys_Dec;   //!
   TBranch        *b_muplus_Hlt1Phys_TIS;   //!
   TBranch        *b_muplus_Hlt1Phys_TOS;   //!
   TBranch        *b_muplus_Hlt2Global_Dec;   //!
   TBranch        *b_muplus_Hlt2Global_TIS;   //!
   TBranch        *b_muplus_Hlt2Global_TOS;   //!
   TBranch        *b_muplus_Hlt2Phys_Dec;   //!
   TBranch        *b_muplus_Hlt2Phys_TIS;   //!
   TBranch        *b_muplus_Hlt2Phys_TOS;   //!
   TBranch        *b_muplus_TRACK_Type;   //!
   TBranch        *b_muplus_TRACK_Key;   //!
   TBranch        *b_muplus_TRACK_CHI2NDOF;   //!
   TBranch        *b_muplus_TRACK_PCHI2;   //!
   TBranch        *b_muplus_TRACK_MatchCHI2;   //!
   TBranch        *b_muplus_TRACK_GhostProb;   //!
   TBranch        *b_muplus_TRACK_CloneDist;   //!
   TBranch        *b_muplus_TRACK_Likelihood;   //!
   TBranch        *b_muminus_CosTheta;   //!
   TBranch        *b_muminus_MINIP;   //!
   TBranch        *b_muminus_MINIPCHI2;   //!
   TBranch        *b_muminus_MINIPNEXTBEST;   //!
   TBranch        *b_muminus_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_muminus_OWNPV_X;   //!
   TBranch        *b_muminus_OWNPV_Y;   //!
   TBranch        *b_muminus_OWNPV_Z;   //!
   TBranch        *b_muminus_OWNPV_XERR;   //!
   TBranch        *b_muminus_OWNPV_YERR;   //!
   TBranch        *b_muminus_OWNPV_ZERR;   //!
   TBranch        *b_muminus_OWNPV_CHI2;   //!
   TBranch        *b_muminus_OWNPV_NDOF;   //!
   TBranch        *b_muminus_OWNPV_COV_;   //!
   TBranch        *b_muminus_IP_OWNPV;   //!
   TBranch        *b_muminus_IPCHI2_OWNPV;   //!
   TBranch        *b_muminus_TOPPV_X;   //!
   TBranch        *b_muminus_TOPPV_Y;   //!
   TBranch        *b_muminus_TOPPV_Z;   //!
   TBranch        *b_muminus_TOPPV_XERR;   //!
   TBranch        *b_muminus_TOPPV_YERR;   //!
   TBranch        *b_muminus_TOPPV_ZERR;   //!
   TBranch        *b_muminus_TOPPV_CHI2;   //!
   TBranch        *b_muminus_TOPPV_NDOF;   //!
   TBranch        *b_muminus_TOPPV_COV_;   //!
   TBranch        *b_muminus_IP_TOPPV;   //!
   TBranch        *b_muminus_IPCHI2_TOPPV;   //!
   TBranch        *b_muminus_ORIVX_X;   //!
   TBranch        *b_muminus_ORIVX_Y;   //!
   TBranch        *b_muminus_ORIVX_Z;   //!
   TBranch        *b_muminus_ORIVX_XERR;   //!
   TBranch        *b_muminus_ORIVX_YERR;   //!
   TBranch        *b_muminus_ORIVX_ZERR;   //!
   TBranch        *b_muminus_ORIVX_CHI2;   //!
   TBranch        *b_muminus_ORIVX_NDOF;   //!
   TBranch        *b_muminus_ORIVX_COV_;   //!
   TBranch        *b_muminus_IP_ORIVX;   //!
   TBranch        *b_muminus_IPCHI2_ORIVX;   //!
   TBranch        *b_muminus_P;   //!
   TBranch        *b_muminus_PT;   //!
   TBranch        *b_muminus_PE;   //!
   TBranch        *b_muminus_PX;   //!
   TBranch        *b_muminus_PY;   //!
   TBranch        *b_muminus_PZ;   //!
   TBranch        *b_muminus_M;   //!
   TBranch        *b_muminus_ID;   //!
   TBranch        *b_muminus_PIDe;   //!
   TBranch        *b_muminus_PIDmu;   //!
   TBranch        *b_muminus_PIDK;   //!
   TBranch        *b_muminus_PIDp;   //!
   TBranch        *b_muminus_ProbNNe;   //!
   TBranch        *b_muminus_ProbNNk;   //!
   TBranch        *b_muminus_ProbNNp;   //!
   TBranch        *b_muminus_ProbNNpi;   //!
   TBranch        *b_muminus_ProbNNmu;   //!
   TBranch        *b_muminus_ProbNNghost;   //!
   TBranch        *b_muminus_hasMuon;   //!
   TBranch        *b_muminus_isMuon;   //!
   TBranch        *b_muminus_hasRich;   //!
   TBranch        *b_muminus_UsedRichAerogel;   //!
   TBranch        *b_muminus_UsedRich1Gas;   //!
   TBranch        *b_muminus_UsedRich2Gas;   //!
   TBranch        *b_muminus_RichAboveElThres;   //!
   TBranch        *b_muminus_RichAboveMuThres;   //!
   TBranch        *b_muminus_RichAbovePiThres;   //!
   TBranch        *b_muminus_RichAboveKaThres;   //!
   TBranch        *b_muminus_RichAbovePrThres;   //!
   TBranch        *b_muminus_hasCalo;   //!
   TBranch        *b_muminus_L0Global_Dec;   //!
   TBranch        *b_muminus_L0Global_TIS;   //!
   TBranch        *b_muminus_L0Global_TOS;   //!
   TBranch        *b_muminus_Hlt1Global_Dec;   //!
   TBranch        *b_muminus_Hlt1Global_TIS;   //!
   TBranch        *b_muminus_Hlt1Global_TOS;   //!
   TBranch        *b_muminus_Hlt1Phys_Dec;   //!
   TBranch        *b_muminus_Hlt1Phys_TIS;   //!
   TBranch        *b_muminus_Hlt1Phys_TOS;   //!
   TBranch        *b_muminus_Hlt2Global_Dec;   //!
   TBranch        *b_muminus_Hlt2Global_TIS;   //!
   TBranch        *b_muminus_Hlt2Global_TOS;   //!
   TBranch        *b_muminus_Hlt2Phys_Dec;   //!
   TBranch        *b_muminus_Hlt2Phys_TIS;   //!
   TBranch        *b_muminus_Hlt2Phys_TOS;   //!
   TBranch        *b_muminus_TRACK_Type;   //!
   TBranch        *b_muminus_TRACK_Key;   //!
   TBranch        *b_muminus_TRACK_CHI2NDOF;   //!
   TBranch        *b_muminus_TRACK_PCHI2;   //!
   TBranch        *b_muminus_TRACK_MatchCHI2;   //!
   TBranch        *b_muminus_TRACK_GhostProb;   //!
   TBranch        *b_muminus_TRACK_CloneDist;   //!
   TBranch        *b_muminus_TRACK_Likelihood;   //!
   TBranch        *b_K_CosTheta;   //!
   TBranch        *b_K_MINIP;   //!
   TBranch        *b_K_MINIPCHI2;   //!
   TBranch        *b_K_MINIPNEXTBEST;   //!
   TBranch        *b_K_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_K_OWNPV_X;   //!
   TBranch        *b_K_OWNPV_Y;   //!
   TBranch        *b_K_OWNPV_Z;   //!
   TBranch        *b_K_OWNPV_XERR;   //!
   TBranch        *b_K_OWNPV_YERR;   //!
   TBranch        *b_K_OWNPV_ZERR;   //!
   TBranch        *b_K_OWNPV_CHI2;   //!
   TBranch        *b_K_OWNPV_NDOF;   //!
   TBranch        *b_K_OWNPV_COV_;   //!
   TBranch        *b_K_IP_OWNPV;   //!
   TBranch        *b_K_IPCHI2_OWNPV;   //!
   TBranch        *b_K_TOPPV_X;   //!
   TBranch        *b_K_TOPPV_Y;   //!
   TBranch        *b_K_TOPPV_Z;   //!
   TBranch        *b_K_TOPPV_XERR;   //!
   TBranch        *b_K_TOPPV_YERR;   //!
   TBranch        *b_K_TOPPV_ZERR;   //!
   TBranch        *b_K_TOPPV_CHI2;   //!
   TBranch        *b_K_TOPPV_NDOF;   //!
   TBranch        *b_K_TOPPV_COV_;   //!
   TBranch        *b_K_IP_TOPPV;   //!
   TBranch        *b_K_IPCHI2_TOPPV;   //!
   TBranch        *b_K_ORIVX_X;   //!
   TBranch        *b_K_ORIVX_Y;   //!
   TBranch        *b_K_ORIVX_Z;   //!
   TBranch        *b_K_ORIVX_XERR;   //!
   TBranch        *b_K_ORIVX_YERR;   //!
   TBranch        *b_K_ORIVX_ZERR;   //!
   TBranch        *b_K_ORIVX_CHI2;   //!
   TBranch        *b_K_ORIVX_NDOF;   //!
   TBranch        *b_K_ORIVX_COV_;   //!
   TBranch        *b_K_IP_ORIVX;   //!
   TBranch        *b_K_IPCHI2_ORIVX;   //!
   TBranch        *b_K_P;   //!
   TBranch        *b_K_PT;   //!
   TBranch        *b_K_PE;   //!
   TBranch        *b_K_PX;   //!
   TBranch        *b_K_PY;   //!
   TBranch        *b_K_PZ;   //!
   TBranch        *b_K_M;   //!
   TBranch        *b_K_ID;   //!
   TBranch        *b_K_PIDe;   //!
   TBranch        *b_K_PIDmu;   //!
   TBranch        *b_K_PIDK;   //!
   TBranch        *b_K_PIDp;   //!
   TBranch        *b_K_ProbNNe;   //!
   TBranch        *b_K_ProbNNk;   //!
   TBranch        *b_K_ProbNNp;   //!
   TBranch        *b_K_ProbNNpi;   //!
   TBranch        *b_K_ProbNNmu;   //!
   TBranch        *b_K_ProbNNghost;   //!
   TBranch        *b_K_hasMuon;   //!
   TBranch        *b_K_isMuon;   //!
   TBranch        *b_K_hasRich;   //!
   TBranch        *b_K_UsedRichAerogel;   //!
   TBranch        *b_K_UsedRich1Gas;   //!
   TBranch        *b_K_UsedRich2Gas;   //!
   TBranch        *b_K_RichAboveElThres;   //!
   TBranch        *b_K_RichAboveMuThres;   //!
   TBranch        *b_K_RichAbovePiThres;   //!
   TBranch        *b_K_RichAboveKaThres;   //!
   TBranch        *b_K_RichAbovePrThres;   //!
   TBranch        *b_K_hasCalo;   //!
   TBranch        *b_K_L0Global_Dec;   //!
   TBranch        *b_K_L0Global_TIS;   //!
   TBranch        *b_K_L0Global_TOS;   //!
   TBranch        *b_K_Hlt1Global_Dec;   //!
   TBranch        *b_K_Hlt1Global_TIS;   //!
   TBranch        *b_K_Hlt1Global_TOS;   //!
   TBranch        *b_K_Hlt1Phys_Dec;   //!
   TBranch        *b_K_Hlt1Phys_TIS;   //!
   TBranch        *b_K_Hlt1Phys_TOS;   //!
   TBranch        *b_K_Hlt2Global_Dec;   //!
   TBranch        *b_K_Hlt2Global_TIS;   //!
   TBranch        *b_K_Hlt2Global_TOS;   //!
   TBranch        *b_K_Hlt2Phys_Dec;   //!
   TBranch        *b_K_Hlt2Phys_TIS;   //!
   TBranch        *b_K_Hlt2Phys_TOS;   //!
   TBranch        *b_K_TRACK_Type;   //!
   TBranch        *b_K_TRACK_Key;   //!
   TBranch        *b_K_TRACK_CHI2NDOF;   //!
   TBranch        *b_K_TRACK_PCHI2;   //!
   TBranch        *b_K_TRACK_MatchCHI2;   //!
   TBranch        *b_K_TRACK_GhostProb;   //!
   TBranch        *b_K_TRACK_CloneDist;   //!
   TBranch        *b_K_TRACK_Likelihood;   //!
   TBranch        *b_K_ETA;   //!
   TBranch        *b_K_PHI;   //!
   TBranch        *b_nCandidate;   //!
   TBranch        *b_totCandidates;   //!
   TBranch        *b_EventInSequence;   //!
   TBranch        *b_nlong;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_BCID;   //!
   TBranch        *b_BCType;   //!
   TBranch        *b_OdinTCK;   //!
   TBranch        *b_L0DUTCK;   //!
   TBranch        *b_HLT1TCK;   //!
   TBranch        *b_HLT2TCK;   //!
   TBranch        *b_GpsTime;   //!
   TBranch        *b_Polarity;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_PVX;   //!
   TBranch        *b_PVY;   //!
   TBranch        *b_PVZ;   //!
   TBranch        *b_PVXERR;   //!
   TBranch        *b_PVYERR;   //!
   TBranch        *b_PVZERR;   //!
   TBranch        *b_PVCHI2;   //!
   TBranch        *b_PVNDOF;   //!
   TBranch        *b_PVNTRACKS;   //!
   TBranch        *b_StrippingFullDSTDiMuonJpsi2MuMuDetachedLineDecision;   //!

   cuts(TTree *tree=0);
   virtual ~cuts();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef cuts_cxx
cuts::cuts(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("B2JpsiKTree/MyTuple",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("B2JpsiKTree/MyTuple","");
      chain->Add("../2016_MagDown/1170_0/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_1/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_2/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_3/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_4/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_5/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_6/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_7/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_8/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_9/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_10/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_11/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_12/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_13/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_14/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_15/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_16/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_17/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_18/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_19/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_20/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_21/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_22/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_23/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_24/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_25/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_26/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_27/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_28/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_29/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_30/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_31/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_32/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_33/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_34/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_35/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_36/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_37/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_38/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_39/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_40/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_41/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_42/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_43/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_45/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_46/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_47/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_48/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_49/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_50/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_51/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_52/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_53/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_54/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_55/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_56/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_57/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_58/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_59/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_60/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_61/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_62/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_63/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_64/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_65/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_66/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_67/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_68/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_69/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_70/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_71/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_72/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_73/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_74/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_75/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_76/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_77/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_78/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_79/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_80/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_81/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_82/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_83/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_84/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_85/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_86/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_87/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_88/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_89/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_90/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_91/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_92/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_93/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_94/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_95/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_96/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_97/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_98/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_99/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_100/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_101/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_102/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_103/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_104/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_105/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_106/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_107/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_108/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_109/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_110/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_111/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_112/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_113/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_114/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_115/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_116/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_117/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_118/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_119/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_120/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_121/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_122/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_123/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_124/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_125/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_126/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_127/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_128/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_129/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_130/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_131/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_132/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_133/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_134/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_135/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_136/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_137/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_138/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_139/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_140/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_141/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_142/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_143/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_144/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_145/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_146/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_147/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_148/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_149/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_150/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_151/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_152/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_153/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_154/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_155/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_156/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_157/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_158/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_159/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_160/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_161/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_162/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_163/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_164/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_165/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_166/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_167/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_168/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_169/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_170/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_171/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_172/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_173/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_174/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_175/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_176/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_177/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_178/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_179/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_180/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_181/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_182/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_183/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_184/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_185/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_186/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_187/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_188/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_189/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_190/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_191/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_192/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_193/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_194/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_195/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_196/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_197/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_198/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_199/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_200/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_201/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_202/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_203/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_204/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_205/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_206/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_207/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_208/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_209/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_210/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_211/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_212/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_213/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_214/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_215/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_216/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_217/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_218/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_219/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_220/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_221/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_222/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_223/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_224/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_225/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_226/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_227/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_228/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_229/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_230/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_231/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_232/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_233/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_234/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_235/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_236/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_237/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_238/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_239/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_240/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_241/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_242/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_243/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_244/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_245/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_246/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_247/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_248/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_249/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_250/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_251/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_252/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_253/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_254/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_255/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_256/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_257/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_258/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_259/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_260/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_261/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_262/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_263/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_264/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_265/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_266/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_267/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_268/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_270/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_271/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_272/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_273/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_274/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_275/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_276/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_277/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_278/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_279/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_280/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_281/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_282/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_283/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_284/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_285/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_286/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_287/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_288/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_289/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_290/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_291/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_292/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_293/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_294/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_295/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_296/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_297/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_298/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_299/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_300/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_301/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_302/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_303/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_304/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_305/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_306/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_307/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_308/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_309/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_310/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_311/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_312/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_313/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_314/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_315/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_316/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_317/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_318/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_319/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_320/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_321/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_322/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_323/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_324/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_325/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_326/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_327/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_328/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_329/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_330/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_331/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_332/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_333/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_334/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_335/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_336/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_337/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_338/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_339/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_340/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_341/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_342/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_343/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_344/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_345/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_346/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_347/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_348/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_349/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_350/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_351/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_352/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_353/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_354/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_355/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_356/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_357/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_358/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_359/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_360/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_361/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_362/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_363/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_364/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_365/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_366/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_367/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_368/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_369/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_370/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_371/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_372/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_373/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_374/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_375/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_376/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_377/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_378/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_379/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_380/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_381/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_382/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_383/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_384/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_385/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_386/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_387/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_388/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_389/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_390/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_391/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_392/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_393/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_394/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_395/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_396/jpsik.root/B2JpsiKTree/MyTuple");
      chain->Add("../2016_MagDown/1170_397/jpsik.root/B2JpsiKTree/MyTuple");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

cuts::~cuts()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t cuts::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t cuts::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void cuts::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("B_MINIP", &B_MINIP, &b_B_MINIP);
   fChain->SetBranchAddress("B_MINIPCHI2", &B_MINIPCHI2, &b_B_MINIPCHI2);
   fChain->SetBranchAddress("B_MINIPNEXTBEST", &B_MINIPNEXTBEST, &b_B_MINIPNEXTBEST);
   fChain->SetBranchAddress("B_MINIPCHI2NEXTBEST", &B_MINIPCHI2NEXTBEST, &b_B_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("B_ENDVERTEX_X", &B_ENDVERTEX_X, &b_B_ENDVERTEX_X);
   fChain->SetBranchAddress("B_ENDVERTEX_Y", &B_ENDVERTEX_Y, &b_B_ENDVERTEX_Y);
   fChain->SetBranchAddress("B_ENDVERTEX_Z", &B_ENDVERTEX_Z, &b_B_ENDVERTEX_Z);
   fChain->SetBranchAddress("B_ENDVERTEX_XERR", &B_ENDVERTEX_XERR, &b_B_ENDVERTEX_XERR);
   fChain->SetBranchAddress("B_ENDVERTEX_YERR", &B_ENDVERTEX_YERR, &b_B_ENDVERTEX_YERR);
   fChain->SetBranchAddress("B_ENDVERTEX_ZERR", &B_ENDVERTEX_ZERR, &b_B_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("B_ENDVERTEX_CHI2", &B_ENDVERTEX_CHI2, &b_B_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("B_ENDVERTEX_NDOF", &B_ENDVERTEX_NDOF, &b_B_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("B_ENDVERTEX_COV_", B_ENDVERTEX_COV_, &b_B_ENDVERTEX_COV_);
   fChain->SetBranchAddress("B_OWNPV_X", &B_OWNPV_X, &b_B_OWNPV_X);
   fChain->SetBranchAddress("B_OWNPV_Y", &B_OWNPV_Y, &b_B_OWNPV_Y);
   fChain->SetBranchAddress("B_OWNPV_Z", &B_OWNPV_Z, &b_B_OWNPV_Z);
   fChain->SetBranchAddress("B_OWNPV_XERR", &B_OWNPV_XERR, &b_B_OWNPV_XERR);
   fChain->SetBranchAddress("B_OWNPV_YERR", &B_OWNPV_YERR, &b_B_OWNPV_YERR);
   fChain->SetBranchAddress("B_OWNPV_ZERR", &B_OWNPV_ZERR, &b_B_OWNPV_ZERR);
   fChain->SetBranchAddress("B_OWNPV_CHI2", &B_OWNPV_CHI2, &b_B_OWNPV_CHI2);
   fChain->SetBranchAddress("B_OWNPV_NDOF", &B_OWNPV_NDOF, &b_B_OWNPV_NDOF);
   fChain->SetBranchAddress("B_OWNPV_COV_", B_OWNPV_COV_, &b_B_OWNPV_COV_);
   fChain->SetBranchAddress("B_IP_OWNPV", &B_IP_OWNPV, &b_B_IP_OWNPV);
   fChain->SetBranchAddress("B_IPCHI2_OWNPV", &B_IPCHI2_OWNPV, &b_B_IPCHI2_OWNPV);
   fChain->SetBranchAddress("B_FD_OWNPV", &B_FD_OWNPV, &b_B_FD_OWNPV);
   fChain->SetBranchAddress("B_FDCHI2_OWNPV", &B_FDCHI2_OWNPV, &b_B_FDCHI2_OWNPV);
   fChain->SetBranchAddress("B_DIRA_OWNPV", &B_DIRA_OWNPV, &b_B_DIRA_OWNPV);
   fChain->SetBranchAddress("B_TOPPV_X", &B_TOPPV_X, &b_B_TOPPV_X);
   fChain->SetBranchAddress("B_TOPPV_Y", &B_TOPPV_Y, &b_B_TOPPV_Y);
   fChain->SetBranchAddress("B_TOPPV_Z", &B_TOPPV_Z, &b_B_TOPPV_Z);
   fChain->SetBranchAddress("B_TOPPV_XERR", &B_TOPPV_XERR, &b_B_TOPPV_XERR);
   fChain->SetBranchAddress("B_TOPPV_YERR", &B_TOPPV_YERR, &b_B_TOPPV_YERR);
   fChain->SetBranchAddress("B_TOPPV_ZERR", &B_TOPPV_ZERR, &b_B_TOPPV_ZERR);
   fChain->SetBranchAddress("B_TOPPV_CHI2", &B_TOPPV_CHI2, &b_B_TOPPV_CHI2);
   fChain->SetBranchAddress("B_TOPPV_NDOF", &B_TOPPV_NDOF, &b_B_TOPPV_NDOF);
   fChain->SetBranchAddress("B_TOPPV_COV_", B_TOPPV_COV_, &b_B_TOPPV_COV_);
   fChain->SetBranchAddress("B_IP_TOPPV", &B_IP_TOPPV, &b_B_IP_TOPPV);
   fChain->SetBranchAddress("B_IPCHI2_TOPPV", &B_IPCHI2_TOPPV, &b_B_IPCHI2_TOPPV);
   fChain->SetBranchAddress("B_FD_TOPPV", &B_FD_TOPPV, &b_B_FD_TOPPV);
   fChain->SetBranchAddress("B_FDCHI2_TOPPV", &B_FDCHI2_TOPPV, &b_B_FDCHI2_TOPPV);
   fChain->SetBranchAddress("B_DIRA_TOPPV", &B_DIRA_TOPPV, &b_B_DIRA_TOPPV);
   fChain->SetBranchAddress("B_P", &B_P, &b_B_P);
   fChain->SetBranchAddress("B_PT", &B_PT, &b_B_PT);
   fChain->SetBranchAddress("B_PE", &B_PE, &b_B_PE);
   fChain->SetBranchAddress("B_PX", &B_PX, &b_B_PX);
   fChain->SetBranchAddress("B_PY", &B_PY, &b_B_PY);
   fChain->SetBranchAddress("B_PZ", &B_PZ, &b_B_PZ);
   fChain->SetBranchAddress("B_MM", &B_MM, &b_B_MM);
   fChain->SetBranchAddress("B_MMERR", &B_MMERR, &b_B_MMERR);
   fChain->SetBranchAddress("B_M", &B_M, &b_B_M);
   fChain->SetBranchAddress("B_ID", &B_ID, &b_B_ID);
   fChain->SetBranchAddress("B_TAU", &B_TAU, &b_B_TAU);
   fChain->SetBranchAddress("B_TAUERR", &B_TAUERR, &b_B_TAUERR);
   fChain->SetBranchAddress("B_TAUCHI2", &B_TAUCHI2, &b_B_TAUCHI2);
   fChain->SetBranchAddress("B_L0Global_Dec", &B_L0Global_Dec, &b_B_L0Global_Dec);
   fChain->SetBranchAddress("B_L0Global_TIS", &B_L0Global_TIS, &b_B_L0Global_TIS);
   fChain->SetBranchAddress("B_L0Global_TOS", &B_L0Global_TOS, &b_B_L0Global_TOS);
   fChain->SetBranchAddress("B_Hlt1Global_Dec", &B_Hlt1Global_Dec, &b_B_Hlt1Global_Dec);
   fChain->SetBranchAddress("B_Hlt1Global_TIS", &B_Hlt1Global_TIS, &b_B_Hlt1Global_TIS);
   fChain->SetBranchAddress("B_Hlt1Global_TOS", &B_Hlt1Global_TOS, &b_B_Hlt1Global_TOS);
   fChain->SetBranchAddress("B_Hlt1Phys_Dec", &B_Hlt1Phys_Dec, &b_B_Hlt1Phys_Dec);
   fChain->SetBranchAddress("B_Hlt1Phys_TIS", &B_Hlt1Phys_TIS, &b_B_Hlt1Phys_TIS);
   fChain->SetBranchAddress("B_Hlt1Phys_TOS", &B_Hlt1Phys_TOS, &b_B_Hlt1Phys_TOS);
   fChain->SetBranchAddress("B_Hlt2Global_Dec", &B_Hlt2Global_Dec, &b_B_Hlt2Global_Dec);
   fChain->SetBranchAddress("B_Hlt2Global_TIS", &B_Hlt2Global_TIS, &b_B_Hlt2Global_TIS);
   fChain->SetBranchAddress("B_Hlt2Global_TOS", &B_Hlt2Global_TOS, &b_B_Hlt2Global_TOS);
   fChain->SetBranchAddress("B_Hlt2Phys_Dec", &B_Hlt2Phys_Dec, &b_B_Hlt2Phys_Dec);
   fChain->SetBranchAddress("B_Hlt2Phys_TIS", &B_Hlt2Phys_TIS, &b_B_Hlt2Phys_TIS);
   fChain->SetBranchAddress("B_Hlt2Phys_TOS", &B_Hlt2Phys_TOS, &b_B_Hlt2Phys_TOS);
   fChain->SetBranchAddress("B_ETA", &B_ETA, &b_B_ETA);
   fChain->SetBranchAddress("B_PHI", &B_PHI, &b_B_PHI);
   fChain->SetBranchAddress("B_DTF_CTAU", &B_DTF_CTAU, &b_B_DTF_CTAU);
   fChain->SetBranchAddress("B_DTF_CTAUERR", &B_DTF_CTAUERR, &b_B_DTF_CTAUERR);
   fChain->SetBranchAddress("B_DTF_CTAUS", &B_DTF_CTAUS, &b_B_DTF_CTAUS);
   fChain->SetBranchAddress("B_DTF_M_JpsiConstr", &B_DTF_M_JpsiConstr, &b_B_DTF_M_JpsiConstr);
   fChain->SetBranchAddress("B_DTF_VCHI2NDOF_JpsiConstr", &B_DTF_VCHI2NDOF_JpsiConstr, &b_B_DTF_VCHI2NDOF_JpsiConstr);
   fChain->SetBranchAddress("B_ConsBPV_nPV", &B_ConsBPV_nPV, &b_B_ConsBPV_nPV);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_M", B_ConsBPV_J_psi_1S_M, &b_B_ConsBPV_J_psi_1S_M);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_MERR", B_ConsBPV_J_psi_1S_MERR, &b_B_ConsBPV_J_psi_1S_MERR);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_P", B_ConsBPV_J_psi_1S_P, &b_B_ConsBPV_J_psi_1S_P);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_PERR", B_ConsBPV_J_psi_1S_PERR, &b_B_ConsBPV_J_psi_1S_PERR);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_ctau", B_ConsBPV_J_psi_1S_ctau, &b_B_ConsBPV_J_psi_1S_ctau);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_ctauErr", B_ConsBPV_J_psi_1S_ctauErr, &b_B_ConsBPV_J_psi_1S_ctauErr);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_decayLength", B_ConsBPV_J_psi_1S_decayLength, &b_B_ConsBPV_J_psi_1S_decayLength);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_decayLengthErr", B_ConsBPV_J_psi_1S_decayLengthErr, &b_B_ConsBPV_J_psi_1S_decayLengthErr);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_0_ID", B_ConsBPV_J_psi_1S_muminus_0_ID, &b_B_ConsBPV_J_psi_1S_muminus_0_ID);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_0_PE", B_ConsBPV_J_psi_1S_muminus_0_PE, &b_B_ConsBPV_J_psi_1S_muminus_0_PE);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_0_PX", B_ConsBPV_J_psi_1S_muminus_0_PX, &b_B_ConsBPV_J_psi_1S_muminus_0_PX);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_0_PY", B_ConsBPV_J_psi_1S_muminus_0_PY, &b_B_ConsBPV_J_psi_1S_muminus_0_PY);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_0_PZ", B_ConsBPV_J_psi_1S_muminus_0_PZ, &b_B_ConsBPV_J_psi_1S_muminus_0_PZ);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_ID", B_ConsBPV_J_psi_1S_muminus_ID, &b_B_ConsBPV_J_psi_1S_muminus_ID);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_PE", B_ConsBPV_J_psi_1S_muminus_PE, &b_B_ConsBPV_J_psi_1S_muminus_PE);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_PX", B_ConsBPV_J_psi_1S_muminus_PX, &b_B_ConsBPV_J_psi_1S_muminus_PX);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_PY", B_ConsBPV_J_psi_1S_muminus_PY, &b_B_ConsBPV_J_psi_1S_muminus_PY);
   fChain->SetBranchAddress("B_ConsBPV_J_psi_1S_muminus_PZ", B_ConsBPV_J_psi_1S_muminus_PZ, &b_B_ConsBPV_J_psi_1S_muminus_PZ);
   fChain->SetBranchAddress("B_ConsBPV_Kplus_ID", B_ConsBPV_Kplus_ID, &b_B_ConsBPV_Kplus_ID);
   fChain->SetBranchAddress("B_ConsBPV_Kplus_PE", B_ConsBPV_Kplus_PE, &b_B_ConsBPV_Kplus_PE);
   fChain->SetBranchAddress("B_ConsBPV_Kplus_PX", B_ConsBPV_Kplus_PX, &b_B_ConsBPV_Kplus_PX);
   fChain->SetBranchAddress("B_ConsBPV_Kplus_PY", B_ConsBPV_Kplus_PY, &b_B_ConsBPV_Kplus_PY);
   fChain->SetBranchAddress("B_ConsBPV_Kplus_PZ", B_ConsBPV_Kplus_PZ, &b_B_ConsBPV_Kplus_PZ);
   fChain->SetBranchAddress("B_ConsBPV_M", B_ConsBPV_M, &b_B_ConsBPV_M);
   fChain->SetBranchAddress("B_ConsBPV_MERR", B_ConsBPV_MERR, &b_B_ConsBPV_MERR);
   fChain->SetBranchAddress("B_ConsBPV_P", B_ConsBPV_P, &b_B_ConsBPV_P);
   fChain->SetBranchAddress("B_ConsBPV_PERR", B_ConsBPV_PERR, &b_B_ConsBPV_PERR);
   fChain->SetBranchAddress("B_ConsBPV_chi2", B_ConsBPV_chi2, &b_B_ConsBPV_chi2);
   fChain->SetBranchAddress("B_ConsBPV_nDOF", B_ConsBPV_nDOF, &b_B_ConsBPV_nDOF);
   fChain->SetBranchAddress("B_ConsBPV_nIter", B_ConsBPV_nIter, &b_B_ConsBPV_nIter);
   fChain->SetBranchAddress("B_ConsBPV_status", B_ConsBPV_status, &b_B_ConsBPV_status);
   fChain->SetBranchAddress("B_L0MuonDecision_Dec", &B_L0MuonDecision_Dec, &b_B_L0MuonDecision_Dec);
   fChain->SetBranchAddress("B_L0MuonDecision_TIS", &B_L0MuonDecision_TIS, &b_B_L0MuonDecision_TIS);
   fChain->SetBranchAddress("B_L0MuonDecision_TOS", &B_L0MuonDecision_TOS, &b_B_L0MuonDecision_TOS);
   fChain->SetBranchAddress("B_L0DiMuonDecision_Dec", &B_L0DiMuonDecision_Dec, &b_B_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("B_L0DiMuonDecision_TIS", &B_L0DiMuonDecision_TIS, &b_B_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("B_L0DiMuonDecision_TOS", &B_L0DiMuonDecision_TOS, &b_B_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("B_Hlt1SingleMuonNoIPDecision_Dec", &B_Hlt1SingleMuonNoIPDecision_Dec, &b_B_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("B_Hlt1SingleMuonNoIPDecision_TIS", &B_Hlt1SingleMuonNoIPDecision_TIS, &b_B_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("B_Hlt1SingleMuonNoIPDecision_TOS", &B_Hlt1SingleMuonNoIPDecision_TOS, &b_B_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("B_Hlt1DiMuonHighMassDecision_Dec", &B_Hlt1DiMuonHighMassDecision_Dec, &b_B_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("B_Hlt1DiMuonHighMassDecision_TIS", &B_Hlt1DiMuonHighMassDecision_TIS, &b_B_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("B_Hlt1DiMuonHighMassDecision_TOS", &B_Hlt1DiMuonHighMassDecision_TOS, &b_B_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("B_Hlt1TrackAllL0Decision_Dec", &B_Hlt1TrackAllL0Decision_Dec, &b_B_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("B_Hlt1TrackAllL0Decision_TIS", &B_Hlt1TrackAllL0Decision_TIS, &b_B_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("B_Hlt1TrackAllL0Decision_TOS", &B_Hlt1TrackAllL0Decision_TOS, &b_B_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("B_Hlt1TrackMuonDecision_Dec", &B_Hlt1TrackMuonDecision_Dec, &b_B_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("B_Hlt1TrackMuonDecision_TIS", &B_Hlt1TrackMuonDecision_TIS, &b_B_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("B_Hlt1TrackMuonDecision_TOS", &B_Hlt1TrackMuonDecision_TOS, &b_B_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("B_Hlt2DiMuonDetachedHeavyDecision_Dec", &B_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_B_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("B_Hlt2DiMuonDetachedHeavyDecision_TIS", &B_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_B_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("B_Hlt2DiMuonDetachedHeavyDecision_TOS", &B_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_B_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("B_Hlt2DiMuonDetachedJPsiDecision_Dec", &B_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_B_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("B_Hlt2DiMuonDetachedJPsiDecision_TIS", &B_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_B_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("B_Hlt2DiMuonDetachedJPsiDecision_TOS", &B_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_B_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("Jpsi_CosTheta", &Jpsi_CosTheta, &b_Jpsi_CosTheta);
   fChain->SetBranchAddress("Jpsi_MINIP", &Jpsi_MINIP, &b_Jpsi_MINIP);
   fChain->SetBranchAddress("Jpsi_MINIPCHI2", &Jpsi_MINIPCHI2, &b_Jpsi_MINIPCHI2);
   fChain->SetBranchAddress("Jpsi_MINIPNEXTBEST", &Jpsi_MINIPNEXTBEST, &b_Jpsi_MINIPNEXTBEST);
   fChain->SetBranchAddress("Jpsi_MINIPCHI2NEXTBEST", &Jpsi_MINIPCHI2NEXTBEST, &b_Jpsi_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_X", &Jpsi_ENDVERTEX_X, &b_Jpsi_ENDVERTEX_X);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_Y", &Jpsi_ENDVERTEX_Y, &b_Jpsi_ENDVERTEX_Y);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_Z", &Jpsi_ENDVERTEX_Z, &b_Jpsi_ENDVERTEX_Z);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_XERR", &Jpsi_ENDVERTEX_XERR, &b_Jpsi_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_YERR", &Jpsi_ENDVERTEX_YERR, &b_Jpsi_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_ZERR", &Jpsi_ENDVERTEX_ZERR, &b_Jpsi_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_CHI2", &Jpsi_ENDVERTEX_CHI2, &b_Jpsi_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_NDOF", &Jpsi_ENDVERTEX_NDOF, &b_Jpsi_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Jpsi_ENDVERTEX_COV_", Jpsi_ENDVERTEX_COV_, &b_Jpsi_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Jpsi_OWNPV_X", &Jpsi_OWNPV_X, &b_Jpsi_OWNPV_X);
   fChain->SetBranchAddress("Jpsi_OWNPV_Y", &Jpsi_OWNPV_Y, &b_Jpsi_OWNPV_Y);
   fChain->SetBranchAddress("Jpsi_OWNPV_Z", &Jpsi_OWNPV_Z, &b_Jpsi_OWNPV_Z);
   fChain->SetBranchAddress("Jpsi_OWNPV_XERR", &Jpsi_OWNPV_XERR, &b_Jpsi_OWNPV_XERR);
   fChain->SetBranchAddress("Jpsi_OWNPV_YERR", &Jpsi_OWNPV_YERR, &b_Jpsi_OWNPV_YERR);
   fChain->SetBranchAddress("Jpsi_OWNPV_ZERR", &Jpsi_OWNPV_ZERR, &b_Jpsi_OWNPV_ZERR);
   fChain->SetBranchAddress("Jpsi_OWNPV_CHI2", &Jpsi_OWNPV_CHI2, &b_Jpsi_OWNPV_CHI2);
   fChain->SetBranchAddress("Jpsi_OWNPV_NDOF", &Jpsi_OWNPV_NDOF, &b_Jpsi_OWNPV_NDOF);
   fChain->SetBranchAddress("Jpsi_OWNPV_COV_", Jpsi_OWNPV_COV_, &b_Jpsi_OWNPV_COV_);
   fChain->SetBranchAddress("Jpsi_IP_OWNPV", &Jpsi_IP_OWNPV, &b_Jpsi_IP_OWNPV);
   fChain->SetBranchAddress("Jpsi_IPCHI2_OWNPV", &Jpsi_IPCHI2_OWNPV, &b_Jpsi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Jpsi_FD_OWNPV", &Jpsi_FD_OWNPV, &b_Jpsi_FD_OWNPV);
   fChain->SetBranchAddress("Jpsi_FDCHI2_OWNPV", &Jpsi_FDCHI2_OWNPV, &b_Jpsi_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Jpsi_DIRA_OWNPV", &Jpsi_DIRA_OWNPV, &b_Jpsi_DIRA_OWNPV);
   fChain->SetBranchAddress("Jpsi_TOPPV_X", &Jpsi_TOPPV_X, &b_Jpsi_TOPPV_X);
   fChain->SetBranchAddress("Jpsi_TOPPV_Y", &Jpsi_TOPPV_Y, &b_Jpsi_TOPPV_Y);
   fChain->SetBranchAddress("Jpsi_TOPPV_Z", &Jpsi_TOPPV_Z, &b_Jpsi_TOPPV_Z);
   fChain->SetBranchAddress("Jpsi_TOPPV_XERR", &Jpsi_TOPPV_XERR, &b_Jpsi_TOPPV_XERR);
   fChain->SetBranchAddress("Jpsi_TOPPV_YERR", &Jpsi_TOPPV_YERR, &b_Jpsi_TOPPV_YERR);
   fChain->SetBranchAddress("Jpsi_TOPPV_ZERR", &Jpsi_TOPPV_ZERR, &b_Jpsi_TOPPV_ZERR);
   fChain->SetBranchAddress("Jpsi_TOPPV_CHI2", &Jpsi_TOPPV_CHI2, &b_Jpsi_TOPPV_CHI2);
   fChain->SetBranchAddress("Jpsi_TOPPV_NDOF", &Jpsi_TOPPV_NDOF, &b_Jpsi_TOPPV_NDOF);
   fChain->SetBranchAddress("Jpsi_TOPPV_COV_", Jpsi_TOPPV_COV_, &b_Jpsi_TOPPV_COV_);
   fChain->SetBranchAddress("Jpsi_IP_TOPPV", &Jpsi_IP_TOPPV, &b_Jpsi_IP_TOPPV);
   fChain->SetBranchAddress("Jpsi_IPCHI2_TOPPV", &Jpsi_IPCHI2_TOPPV, &b_Jpsi_IPCHI2_TOPPV);
   fChain->SetBranchAddress("Jpsi_FD_TOPPV", &Jpsi_FD_TOPPV, &b_Jpsi_FD_TOPPV);
   fChain->SetBranchAddress("Jpsi_FDCHI2_TOPPV", &Jpsi_FDCHI2_TOPPV, &b_Jpsi_FDCHI2_TOPPV);
   fChain->SetBranchAddress("Jpsi_DIRA_TOPPV", &Jpsi_DIRA_TOPPV, &b_Jpsi_DIRA_TOPPV);
   fChain->SetBranchAddress("Jpsi_ORIVX_X", &Jpsi_ORIVX_X, &b_Jpsi_ORIVX_X);
   fChain->SetBranchAddress("Jpsi_ORIVX_Y", &Jpsi_ORIVX_Y, &b_Jpsi_ORIVX_Y);
   fChain->SetBranchAddress("Jpsi_ORIVX_Z", &Jpsi_ORIVX_Z, &b_Jpsi_ORIVX_Z);
   fChain->SetBranchAddress("Jpsi_ORIVX_XERR", &Jpsi_ORIVX_XERR, &b_Jpsi_ORIVX_XERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_YERR", &Jpsi_ORIVX_YERR, &b_Jpsi_ORIVX_YERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_ZERR", &Jpsi_ORIVX_ZERR, &b_Jpsi_ORIVX_ZERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_CHI2", &Jpsi_ORIVX_CHI2, &b_Jpsi_ORIVX_CHI2);
   fChain->SetBranchAddress("Jpsi_ORIVX_NDOF", &Jpsi_ORIVX_NDOF, &b_Jpsi_ORIVX_NDOF);
   fChain->SetBranchAddress("Jpsi_ORIVX_COV_", Jpsi_ORIVX_COV_, &b_Jpsi_ORIVX_COV_);
   fChain->SetBranchAddress("Jpsi_IP_ORIVX", &Jpsi_IP_ORIVX, &b_Jpsi_IP_ORIVX);
   fChain->SetBranchAddress("Jpsi_IPCHI2_ORIVX", &Jpsi_IPCHI2_ORIVX, &b_Jpsi_IPCHI2_ORIVX);
   fChain->SetBranchAddress("Jpsi_FD_ORIVX", &Jpsi_FD_ORIVX, &b_Jpsi_FD_ORIVX);
   fChain->SetBranchAddress("Jpsi_FDCHI2_ORIVX", &Jpsi_FDCHI2_ORIVX, &b_Jpsi_FDCHI2_ORIVX);
   fChain->SetBranchAddress("Jpsi_DIRA_ORIVX", &Jpsi_DIRA_ORIVX, &b_Jpsi_DIRA_ORIVX);
   fChain->SetBranchAddress("Jpsi_P", &Jpsi_P, &b_Jpsi_P);
   fChain->SetBranchAddress("Jpsi_PT", &Jpsi_PT, &b_Jpsi_PT);
   fChain->SetBranchAddress("Jpsi_PE", &Jpsi_PE, &b_Jpsi_PE);
   fChain->SetBranchAddress("Jpsi_PX", &Jpsi_PX, &b_Jpsi_PX);
   fChain->SetBranchAddress("Jpsi_PY", &Jpsi_PY, &b_Jpsi_PY);
   fChain->SetBranchAddress("Jpsi_PZ", &Jpsi_PZ, &b_Jpsi_PZ);
   fChain->SetBranchAddress("Jpsi_MM", &Jpsi_MM, &b_Jpsi_MM);
   fChain->SetBranchAddress("Jpsi_MMERR", &Jpsi_MMERR, &b_Jpsi_MMERR);
   fChain->SetBranchAddress("Jpsi_M", &Jpsi_M, &b_Jpsi_M);
   fChain->SetBranchAddress("Jpsi_ID", &Jpsi_ID, &b_Jpsi_ID);
   fChain->SetBranchAddress("Jpsi_TAU", &Jpsi_TAU, &b_Jpsi_TAU);
   fChain->SetBranchAddress("Jpsi_TAUERR", &Jpsi_TAUERR, &b_Jpsi_TAUERR);
   fChain->SetBranchAddress("Jpsi_TAUCHI2", &Jpsi_TAUCHI2, &b_Jpsi_TAUCHI2);
   fChain->SetBranchAddress("Jpsi_L0Global_Dec", &Jpsi_L0Global_Dec, &b_Jpsi_L0Global_Dec);
   fChain->SetBranchAddress("Jpsi_L0Global_TIS", &Jpsi_L0Global_TIS, &b_Jpsi_L0Global_TIS);
   fChain->SetBranchAddress("Jpsi_L0Global_TOS", &Jpsi_L0Global_TOS, &b_Jpsi_L0Global_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1Global_Dec", &Jpsi_Hlt1Global_Dec, &b_Jpsi_Hlt1Global_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1Global_TIS", &Jpsi_Hlt1Global_TIS, &b_Jpsi_Hlt1Global_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1Global_TOS", &Jpsi_Hlt1Global_TOS, &b_Jpsi_Hlt1Global_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1Phys_Dec", &Jpsi_Hlt1Phys_Dec, &b_Jpsi_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1Phys_TIS", &Jpsi_Hlt1Phys_TIS, &b_Jpsi_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1Phys_TOS", &Jpsi_Hlt1Phys_TOS, &b_Jpsi_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2Global_Dec", &Jpsi_Hlt2Global_Dec, &b_Jpsi_Hlt2Global_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2Global_TIS", &Jpsi_Hlt2Global_TIS, &b_Jpsi_Hlt2Global_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2Global_TOS", &Jpsi_Hlt2Global_TOS, &b_Jpsi_Hlt2Global_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2Phys_Dec", &Jpsi_Hlt2Phys_Dec, &b_Jpsi_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2Phys_TIS", &Jpsi_Hlt2Phys_TIS, &b_Jpsi_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2Phys_TOS", &Jpsi_Hlt2Phys_TOS, &b_Jpsi_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Jpsi_ETA", &Jpsi_ETA, &b_Jpsi_ETA);
   fChain->SetBranchAddress("Jpsi_PHI", &Jpsi_PHI, &b_Jpsi_PHI);
   fChain->SetBranchAddress("Jpsi_L0MuonDecision_Dec", &Jpsi_L0MuonDecision_Dec, &b_Jpsi_L0MuonDecision_Dec);
   fChain->SetBranchAddress("Jpsi_L0MuonDecision_TIS", &Jpsi_L0MuonDecision_TIS, &b_Jpsi_L0MuonDecision_TIS);
   fChain->SetBranchAddress("Jpsi_L0MuonDecision_TOS", &Jpsi_L0MuonDecision_TOS, &b_Jpsi_L0MuonDecision_TOS);
   fChain->SetBranchAddress("Jpsi_L0DiMuonDecision_Dec", &Jpsi_L0DiMuonDecision_Dec, &b_Jpsi_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("Jpsi_L0DiMuonDecision_TIS", &Jpsi_L0DiMuonDecision_TIS, &b_Jpsi_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("Jpsi_L0DiMuonDecision_TOS", &Jpsi_L0DiMuonDecision_TOS, &b_Jpsi_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1SingleMuonNoIPDecision_Dec", &Jpsi_Hlt1SingleMuonNoIPDecision_Dec, &b_Jpsi_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1SingleMuonNoIPDecision_TIS", &Jpsi_Hlt1SingleMuonNoIPDecision_TIS, &b_Jpsi_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1SingleMuonNoIPDecision_TOS", &Jpsi_Hlt1SingleMuonNoIPDecision_TOS, &b_Jpsi_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_Dec", &Jpsi_Hlt1DiMuonHighMassDecision_Dec, &b_Jpsi_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_TIS", &Jpsi_Hlt1DiMuonHighMassDecision_TIS, &b_Jpsi_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1DiMuonHighMassDecision_TOS", &Jpsi_Hlt1DiMuonHighMassDecision_TOS, &b_Jpsi_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1TrackAllL0Decision_Dec", &Jpsi_Hlt1TrackAllL0Decision_Dec, &b_Jpsi_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1TrackAllL0Decision_TIS", &Jpsi_Hlt1TrackAllL0Decision_TIS, &b_Jpsi_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1TrackAllL0Decision_TOS", &Jpsi_Hlt1TrackAllL0Decision_TOS, &b_Jpsi_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt1TrackMuonDecision_Dec", &Jpsi_Hlt1TrackMuonDecision_Dec, &b_Jpsi_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt1TrackMuonDecision_TIS", &Jpsi_Hlt1TrackMuonDecision_TIS, &b_Jpsi_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt1TrackMuonDecision_TOS", &Jpsi_Hlt1TrackMuonDecision_TOS, &b_Jpsi_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec", &Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS", &Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS", &Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_Jpsi_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec", &Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS", &Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS", &Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_Jpsi_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("muplus_CosTheta", &muplus_CosTheta, &b_muplus_CosTheta);
   fChain->SetBranchAddress("muplus_MINIP", &muplus_MINIP, &b_muplus_MINIP);
   fChain->SetBranchAddress("muplus_MINIPCHI2", &muplus_MINIPCHI2, &b_muplus_MINIPCHI2);
   fChain->SetBranchAddress("muplus_MINIPNEXTBEST", &muplus_MINIPNEXTBEST, &b_muplus_MINIPNEXTBEST);
   fChain->SetBranchAddress("muplus_MINIPCHI2NEXTBEST", &muplus_MINIPCHI2NEXTBEST, &b_muplus_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("muplus_OWNPV_X", &muplus_OWNPV_X, &b_muplus_OWNPV_X);
   fChain->SetBranchAddress("muplus_OWNPV_Y", &muplus_OWNPV_Y, &b_muplus_OWNPV_Y);
   fChain->SetBranchAddress("muplus_OWNPV_Z", &muplus_OWNPV_Z, &b_muplus_OWNPV_Z);
   fChain->SetBranchAddress("muplus_OWNPV_XERR", &muplus_OWNPV_XERR, &b_muplus_OWNPV_XERR);
   fChain->SetBranchAddress("muplus_OWNPV_YERR", &muplus_OWNPV_YERR, &b_muplus_OWNPV_YERR);
   fChain->SetBranchAddress("muplus_OWNPV_ZERR", &muplus_OWNPV_ZERR, &b_muplus_OWNPV_ZERR);
   fChain->SetBranchAddress("muplus_OWNPV_CHI2", &muplus_OWNPV_CHI2, &b_muplus_OWNPV_CHI2);
   fChain->SetBranchAddress("muplus_OWNPV_NDOF", &muplus_OWNPV_NDOF, &b_muplus_OWNPV_NDOF);
   fChain->SetBranchAddress("muplus_OWNPV_COV_", muplus_OWNPV_COV_, &b_muplus_OWNPV_COV_);
   fChain->SetBranchAddress("muplus_IP_OWNPV", &muplus_IP_OWNPV, &b_muplus_IP_OWNPV);
   fChain->SetBranchAddress("muplus_IPCHI2_OWNPV", &muplus_IPCHI2_OWNPV, &b_muplus_IPCHI2_OWNPV);
   fChain->SetBranchAddress("muplus_TOPPV_X", &muplus_TOPPV_X, &b_muplus_TOPPV_X);
   fChain->SetBranchAddress("muplus_TOPPV_Y", &muplus_TOPPV_Y, &b_muplus_TOPPV_Y);
   fChain->SetBranchAddress("muplus_TOPPV_Z", &muplus_TOPPV_Z, &b_muplus_TOPPV_Z);
   fChain->SetBranchAddress("muplus_TOPPV_XERR", &muplus_TOPPV_XERR, &b_muplus_TOPPV_XERR);
   fChain->SetBranchAddress("muplus_TOPPV_YERR", &muplus_TOPPV_YERR, &b_muplus_TOPPV_YERR);
   fChain->SetBranchAddress("muplus_TOPPV_ZERR", &muplus_TOPPV_ZERR, &b_muplus_TOPPV_ZERR);
   fChain->SetBranchAddress("muplus_TOPPV_CHI2", &muplus_TOPPV_CHI2, &b_muplus_TOPPV_CHI2);
   fChain->SetBranchAddress("muplus_TOPPV_NDOF", &muplus_TOPPV_NDOF, &b_muplus_TOPPV_NDOF);
   fChain->SetBranchAddress("muplus_TOPPV_COV_", muplus_TOPPV_COV_, &b_muplus_TOPPV_COV_);
   fChain->SetBranchAddress("muplus_IP_TOPPV", &muplus_IP_TOPPV, &b_muplus_IP_TOPPV);
   fChain->SetBranchAddress("muplus_IPCHI2_TOPPV", &muplus_IPCHI2_TOPPV, &b_muplus_IPCHI2_TOPPV);
   fChain->SetBranchAddress("muplus_ORIVX_X", &muplus_ORIVX_X, &b_muplus_ORIVX_X);
   fChain->SetBranchAddress("muplus_ORIVX_Y", &muplus_ORIVX_Y, &b_muplus_ORIVX_Y);
   fChain->SetBranchAddress("muplus_ORIVX_Z", &muplus_ORIVX_Z, &b_muplus_ORIVX_Z);
   fChain->SetBranchAddress("muplus_ORIVX_XERR", &muplus_ORIVX_XERR, &b_muplus_ORIVX_XERR);
   fChain->SetBranchAddress("muplus_ORIVX_YERR", &muplus_ORIVX_YERR, &b_muplus_ORIVX_YERR);
   fChain->SetBranchAddress("muplus_ORIVX_ZERR", &muplus_ORIVX_ZERR, &b_muplus_ORIVX_ZERR);
   fChain->SetBranchAddress("muplus_ORIVX_CHI2", &muplus_ORIVX_CHI2, &b_muplus_ORIVX_CHI2);
   fChain->SetBranchAddress("muplus_ORIVX_NDOF", &muplus_ORIVX_NDOF, &b_muplus_ORIVX_NDOF);
   fChain->SetBranchAddress("muplus_ORIVX_COV_", muplus_ORIVX_COV_, &b_muplus_ORIVX_COV_);
   fChain->SetBranchAddress("muplus_IP_ORIVX", &muplus_IP_ORIVX, &b_muplus_IP_ORIVX);
   fChain->SetBranchAddress("muplus_IPCHI2_ORIVX", &muplus_IPCHI2_ORIVX, &b_muplus_IPCHI2_ORIVX);
   fChain->SetBranchAddress("muplus_P", &muplus_P, &b_muplus_P);
   fChain->SetBranchAddress("muplus_PT", &muplus_PT, &b_muplus_PT);
   fChain->SetBranchAddress("muplus_PE", &muplus_PE, &b_muplus_PE);
   fChain->SetBranchAddress("muplus_PX", &muplus_PX, &b_muplus_PX);
   fChain->SetBranchAddress("muplus_PY", &muplus_PY, &b_muplus_PY);
   fChain->SetBranchAddress("muplus_PZ", &muplus_PZ, &b_muplus_PZ);
   fChain->SetBranchAddress("muplus_M", &muplus_M, &b_muplus_M);
   fChain->SetBranchAddress("muplus_ID", &muplus_ID, &b_muplus_ID);
   fChain->SetBranchAddress("muplus_PIDe", &muplus_PIDe, &b_muplus_PIDe);
   fChain->SetBranchAddress("muplus_PIDmu", &muplus_PIDmu, &b_muplus_PIDmu);
   fChain->SetBranchAddress("muplus_PIDK", &muplus_PIDK, &b_muplus_PIDK);
   fChain->SetBranchAddress("muplus_PIDp", &muplus_PIDp, &b_muplus_PIDp);
   fChain->SetBranchAddress("muplus_ProbNNe", &muplus_ProbNNe, &b_muplus_ProbNNe);
   fChain->SetBranchAddress("muplus_ProbNNk", &muplus_ProbNNk, &b_muplus_ProbNNk);
   fChain->SetBranchAddress("muplus_ProbNNp", &muplus_ProbNNp, &b_muplus_ProbNNp);
   fChain->SetBranchAddress("muplus_ProbNNpi", &muplus_ProbNNpi, &b_muplus_ProbNNpi);
   fChain->SetBranchAddress("muplus_ProbNNmu", &muplus_ProbNNmu, &b_muplus_ProbNNmu);
   fChain->SetBranchAddress("muplus_ProbNNghost", &muplus_ProbNNghost, &b_muplus_ProbNNghost);
   fChain->SetBranchAddress("muplus_hasMuon", &muplus_hasMuon, &b_muplus_hasMuon);
   fChain->SetBranchAddress("muplus_isMuon", &muplus_isMuon, &b_muplus_isMuon);
   fChain->SetBranchAddress("muplus_hasRich", &muplus_hasRich, &b_muplus_hasRich);
   fChain->SetBranchAddress("muplus_UsedRichAerogel", &muplus_UsedRichAerogel, &b_muplus_UsedRichAerogel);
   fChain->SetBranchAddress("muplus_UsedRich1Gas", &muplus_UsedRich1Gas, &b_muplus_UsedRich1Gas);
   fChain->SetBranchAddress("muplus_UsedRich2Gas", &muplus_UsedRich2Gas, &b_muplus_UsedRich2Gas);
   fChain->SetBranchAddress("muplus_RichAboveElThres", &muplus_RichAboveElThres, &b_muplus_RichAboveElThres);
   fChain->SetBranchAddress("muplus_RichAboveMuThres", &muplus_RichAboveMuThres, &b_muplus_RichAboveMuThres);
   fChain->SetBranchAddress("muplus_RichAbovePiThres", &muplus_RichAbovePiThres, &b_muplus_RichAbovePiThres);
   fChain->SetBranchAddress("muplus_RichAboveKaThres", &muplus_RichAboveKaThres, &b_muplus_RichAboveKaThres);
   fChain->SetBranchAddress("muplus_RichAbovePrThres", &muplus_RichAbovePrThres, &b_muplus_RichAbovePrThres);
   fChain->SetBranchAddress("muplus_hasCalo", &muplus_hasCalo, &b_muplus_hasCalo);
   fChain->SetBranchAddress("muplus_L0Global_Dec", &muplus_L0Global_Dec, &b_muplus_L0Global_Dec);
   fChain->SetBranchAddress("muplus_L0Global_TIS", &muplus_L0Global_TIS, &b_muplus_L0Global_TIS);
   fChain->SetBranchAddress("muplus_L0Global_TOS", &muplus_L0Global_TOS, &b_muplus_L0Global_TOS);
   fChain->SetBranchAddress("muplus_Hlt1Global_Dec", &muplus_Hlt1Global_Dec, &b_muplus_Hlt1Global_Dec);
   fChain->SetBranchAddress("muplus_Hlt1Global_TIS", &muplus_Hlt1Global_TIS, &b_muplus_Hlt1Global_TIS);
   fChain->SetBranchAddress("muplus_Hlt1Global_TOS", &muplus_Hlt1Global_TOS, &b_muplus_Hlt1Global_TOS);
   fChain->SetBranchAddress("muplus_Hlt1Phys_Dec", &muplus_Hlt1Phys_Dec, &b_muplus_Hlt1Phys_Dec);
   fChain->SetBranchAddress("muplus_Hlt1Phys_TIS", &muplus_Hlt1Phys_TIS, &b_muplus_Hlt1Phys_TIS);
   fChain->SetBranchAddress("muplus_Hlt1Phys_TOS", &muplus_Hlt1Phys_TOS, &b_muplus_Hlt1Phys_TOS);
   fChain->SetBranchAddress("muplus_Hlt2Global_Dec", &muplus_Hlt2Global_Dec, &b_muplus_Hlt2Global_Dec);
   fChain->SetBranchAddress("muplus_Hlt2Global_TIS", &muplus_Hlt2Global_TIS, &b_muplus_Hlt2Global_TIS);
   fChain->SetBranchAddress("muplus_Hlt2Global_TOS", &muplus_Hlt2Global_TOS, &b_muplus_Hlt2Global_TOS);
   fChain->SetBranchAddress("muplus_Hlt2Phys_Dec", &muplus_Hlt2Phys_Dec, &b_muplus_Hlt2Phys_Dec);
   fChain->SetBranchAddress("muplus_Hlt2Phys_TIS", &muplus_Hlt2Phys_TIS, &b_muplus_Hlt2Phys_TIS);
   fChain->SetBranchAddress("muplus_Hlt2Phys_TOS", &muplus_Hlt2Phys_TOS, &b_muplus_Hlt2Phys_TOS);
   fChain->SetBranchAddress("muplus_TRACK_Type", &muplus_TRACK_Type, &b_muplus_TRACK_Type);
   fChain->SetBranchAddress("muplus_TRACK_Key", &muplus_TRACK_Key, &b_muplus_TRACK_Key);
   fChain->SetBranchAddress("muplus_TRACK_CHI2NDOF", &muplus_TRACK_CHI2NDOF, &b_muplus_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("muplus_TRACK_PCHI2", &muplus_TRACK_PCHI2, &b_muplus_TRACK_PCHI2);
   fChain->SetBranchAddress("muplus_TRACK_MatchCHI2", &muplus_TRACK_MatchCHI2, &b_muplus_TRACK_MatchCHI2);
   fChain->SetBranchAddress("muplus_TRACK_GhostProb", &muplus_TRACK_GhostProb, &b_muplus_TRACK_GhostProb);
   fChain->SetBranchAddress("muplus_TRACK_CloneDist", &muplus_TRACK_CloneDist, &b_muplus_TRACK_CloneDist);
   fChain->SetBranchAddress("muplus_TRACK_Likelihood", &muplus_TRACK_Likelihood, &b_muplus_TRACK_Likelihood);
   fChain->SetBranchAddress("muminus_CosTheta", &muminus_CosTheta, &b_muminus_CosTheta);
   fChain->SetBranchAddress("muminus_MINIP", &muminus_MINIP, &b_muminus_MINIP);
   fChain->SetBranchAddress("muminus_MINIPCHI2", &muminus_MINIPCHI2, &b_muminus_MINIPCHI2);
   fChain->SetBranchAddress("muminus_MINIPNEXTBEST", &muminus_MINIPNEXTBEST, &b_muminus_MINIPNEXTBEST);
   fChain->SetBranchAddress("muminus_MINIPCHI2NEXTBEST", &muminus_MINIPCHI2NEXTBEST, &b_muminus_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("muminus_OWNPV_X", &muminus_OWNPV_X, &b_muminus_OWNPV_X);
   fChain->SetBranchAddress("muminus_OWNPV_Y", &muminus_OWNPV_Y, &b_muminus_OWNPV_Y);
   fChain->SetBranchAddress("muminus_OWNPV_Z", &muminus_OWNPV_Z, &b_muminus_OWNPV_Z);
   fChain->SetBranchAddress("muminus_OWNPV_XERR", &muminus_OWNPV_XERR, &b_muminus_OWNPV_XERR);
   fChain->SetBranchAddress("muminus_OWNPV_YERR", &muminus_OWNPV_YERR, &b_muminus_OWNPV_YERR);
   fChain->SetBranchAddress("muminus_OWNPV_ZERR", &muminus_OWNPV_ZERR, &b_muminus_OWNPV_ZERR);
   fChain->SetBranchAddress("muminus_OWNPV_CHI2", &muminus_OWNPV_CHI2, &b_muminus_OWNPV_CHI2);
   fChain->SetBranchAddress("muminus_OWNPV_NDOF", &muminus_OWNPV_NDOF, &b_muminus_OWNPV_NDOF);
   fChain->SetBranchAddress("muminus_OWNPV_COV_", muminus_OWNPV_COV_, &b_muminus_OWNPV_COV_);
   fChain->SetBranchAddress("muminus_IP_OWNPV", &muminus_IP_OWNPV, &b_muminus_IP_OWNPV);
   fChain->SetBranchAddress("muminus_IPCHI2_OWNPV", &muminus_IPCHI2_OWNPV, &b_muminus_IPCHI2_OWNPV);
   fChain->SetBranchAddress("muminus_TOPPV_X", &muminus_TOPPV_X, &b_muminus_TOPPV_X);
   fChain->SetBranchAddress("muminus_TOPPV_Y", &muminus_TOPPV_Y, &b_muminus_TOPPV_Y);
   fChain->SetBranchAddress("muminus_TOPPV_Z", &muminus_TOPPV_Z, &b_muminus_TOPPV_Z);
   fChain->SetBranchAddress("muminus_TOPPV_XERR", &muminus_TOPPV_XERR, &b_muminus_TOPPV_XERR);
   fChain->SetBranchAddress("muminus_TOPPV_YERR", &muminus_TOPPV_YERR, &b_muminus_TOPPV_YERR);
   fChain->SetBranchAddress("muminus_TOPPV_ZERR", &muminus_TOPPV_ZERR, &b_muminus_TOPPV_ZERR);
   fChain->SetBranchAddress("muminus_TOPPV_CHI2", &muminus_TOPPV_CHI2, &b_muminus_TOPPV_CHI2);
   fChain->SetBranchAddress("muminus_TOPPV_NDOF", &muminus_TOPPV_NDOF, &b_muminus_TOPPV_NDOF);
   fChain->SetBranchAddress("muminus_TOPPV_COV_", muminus_TOPPV_COV_, &b_muminus_TOPPV_COV_);
   fChain->SetBranchAddress("muminus_IP_TOPPV", &muminus_IP_TOPPV, &b_muminus_IP_TOPPV);
   fChain->SetBranchAddress("muminus_IPCHI2_TOPPV", &muminus_IPCHI2_TOPPV, &b_muminus_IPCHI2_TOPPV);
   fChain->SetBranchAddress("muminus_ORIVX_X", &muminus_ORIVX_X, &b_muminus_ORIVX_X);
   fChain->SetBranchAddress("muminus_ORIVX_Y", &muminus_ORIVX_Y, &b_muminus_ORIVX_Y);
   fChain->SetBranchAddress("muminus_ORIVX_Z", &muminus_ORIVX_Z, &b_muminus_ORIVX_Z);
   fChain->SetBranchAddress("muminus_ORIVX_XERR", &muminus_ORIVX_XERR, &b_muminus_ORIVX_XERR);
   fChain->SetBranchAddress("muminus_ORIVX_YERR", &muminus_ORIVX_YERR, &b_muminus_ORIVX_YERR);
   fChain->SetBranchAddress("muminus_ORIVX_ZERR", &muminus_ORIVX_ZERR, &b_muminus_ORIVX_ZERR);
   fChain->SetBranchAddress("muminus_ORIVX_CHI2", &muminus_ORIVX_CHI2, &b_muminus_ORIVX_CHI2);
   fChain->SetBranchAddress("muminus_ORIVX_NDOF", &muminus_ORIVX_NDOF, &b_muminus_ORIVX_NDOF);
   fChain->SetBranchAddress("muminus_ORIVX_COV_", muminus_ORIVX_COV_, &b_muminus_ORIVX_COV_);
   fChain->SetBranchAddress("muminus_IP_ORIVX", &muminus_IP_ORIVX, &b_muminus_IP_ORIVX);
   fChain->SetBranchAddress("muminus_IPCHI2_ORIVX", &muminus_IPCHI2_ORIVX, &b_muminus_IPCHI2_ORIVX);
   fChain->SetBranchAddress("muminus_P", &muminus_P, &b_muminus_P);
   fChain->SetBranchAddress("muminus_PT", &muminus_PT, &b_muminus_PT);
   fChain->SetBranchAddress("muminus_PE", &muminus_PE, &b_muminus_PE);
   fChain->SetBranchAddress("muminus_PX", &muminus_PX, &b_muminus_PX);
   fChain->SetBranchAddress("muminus_PY", &muminus_PY, &b_muminus_PY);
   fChain->SetBranchAddress("muminus_PZ", &muminus_PZ, &b_muminus_PZ);
   fChain->SetBranchAddress("muminus_M", &muminus_M, &b_muminus_M);
   fChain->SetBranchAddress("muminus_ID", &muminus_ID, &b_muminus_ID);
   fChain->SetBranchAddress("muminus_PIDe", &muminus_PIDe, &b_muminus_PIDe);
   fChain->SetBranchAddress("muminus_PIDmu", &muminus_PIDmu, &b_muminus_PIDmu);
   fChain->SetBranchAddress("muminus_PIDK", &muminus_PIDK, &b_muminus_PIDK);
   fChain->SetBranchAddress("muminus_PIDp", &muminus_PIDp, &b_muminus_PIDp);
   fChain->SetBranchAddress("muminus_ProbNNe", &muminus_ProbNNe, &b_muminus_ProbNNe);
   fChain->SetBranchAddress("muminus_ProbNNk", &muminus_ProbNNk, &b_muminus_ProbNNk);
   fChain->SetBranchAddress("muminus_ProbNNp", &muminus_ProbNNp, &b_muminus_ProbNNp);
   fChain->SetBranchAddress("muminus_ProbNNpi", &muminus_ProbNNpi, &b_muminus_ProbNNpi);
   fChain->SetBranchAddress("muminus_ProbNNmu", &muminus_ProbNNmu, &b_muminus_ProbNNmu);
   fChain->SetBranchAddress("muminus_ProbNNghost", &muminus_ProbNNghost, &b_muminus_ProbNNghost);
   fChain->SetBranchAddress("muminus_hasMuon", &muminus_hasMuon, &b_muminus_hasMuon);
   fChain->SetBranchAddress("muminus_isMuon", &muminus_isMuon, &b_muminus_isMuon);
   fChain->SetBranchAddress("muminus_hasRich", &muminus_hasRich, &b_muminus_hasRich);
   fChain->SetBranchAddress("muminus_UsedRichAerogel", &muminus_UsedRichAerogel, &b_muminus_UsedRichAerogel);
   fChain->SetBranchAddress("muminus_UsedRich1Gas", &muminus_UsedRich1Gas, &b_muminus_UsedRich1Gas);
   fChain->SetBranchAddress("muminus_UsedRich2Gas", &muminus_UsedRich2Gas, &b_muminus_UsedRich2Gas);
   fChain->SetBranchAddress("muminus_RichAboveElThres", &muminus_RichAboveElThres, &b_muminus_RichAboveElThres);
   fChain->SetBranchAddress("muminus_RichAboveMuThres", &muminus_RichAboveMuThres, &b_muminus_RichAboveMuThres);
   fChain->SetBranchAddress("muminus_RichAbovePiThres", &muminus_RichAbovePiThres, &b_muminus_RichAbovePiThres);
   fChain->SetBranchAddress("muminus_RichAboveKaThres", &muminus_RichAboveKaThres, &b_muminus_RichAboveKaThres);
   fChain->SetBranchAddress("muminus_RichAbovePrThres", &muminus_RichAbovePrThres, &b_muminus_RichAbovePrThres);
   fChain->SetBranchAddress("muminus_hasCalo", &muminus_hasCalo, &b_muminus_hasCalo);
   fChain->SetBranchAddress("muminus_L0Global_Dec", &muminus_L0Global_Dec, &b_muminus_L0Global_Dec);
   fChain->SetBranchAddress("muminus_L0Global_TIS", &muminus_L0Global_TIS, &b_muminus_L0Global_TIS);
   fChain->SetBranchAddress("muminus_L0Global_TOS", &muminus_L0Global_TOS, &b_muminus_L0Global_TOS);
   fChain->SetBranchAddress("muminus_Hlt1Global_Dec", &muminus_Hlt1Global_Dec, &b_muminus_Hlt1Global_Dec);
   fChain->SetBranchAddress("muminus_Hlt1Global_TIS", &muminus_Hlt1Global_TIS, &b_muminus_Hlt1Global_TIS);
   fChain->SetBranchAddress("muminus_Hlt1Global_TOS", &muminus_Hlt1Global_TOS, &b_muminus_Hlt1Global_TOS);
   fChain->SetBranchAddress("muminus_Hlt1Phys_Dec", &muminus_Hlt1Phys_Dec, &b_muminus_Hlt1Phys_Dec);
   fChain->SetBranchAddress("muminus_Hlt1Phys_TIS", &muminus_Hlt1Phys_TIS, &b_muminus_Hlt1Phys_TIS);
   fChain->SetBranchAddress("muminus_Hlt1Phys_TOS", &muminus_Hlt1Phys_TOS, &b_muminus_Hlt1Phys_TOS);
   fChain->SetBranchAddress("muminus_Hlt2Global_Dec", &muminus_Hlt2Global_Dec, &b_muminus_Hlt2Global_Dec);
   fChain->SetBranchAddress("muminus_Hlt2Global_TIS", &muminus_Hlt2Global_TIS, &b_muminus_Hlt2Global_TIS);
   fChain->SetBranchAddress("muminus_Hlt2Global_TOS", &muminus_Hlt2Global_TOS, &b_muminus_Hlt2Global_TOS);
   fChain->SetBranchAddress("muminus_Hlt2Phys_Dec", &muminus_Hlt2Phys_Dec, &b_muminus_Hlt2Phys_Dec);
   fChain->SetBranchAddress("muminus_Hlt2Phys_TIS", &muminus_Hlt2Phys_TIS, &b_muminus_Hlt2Phys_TIS);
   fChain->SetBranchAddress("muminus_Hlt2Phys_TOS", &muminus_Hlt2Phys_TOS, &b_muminus_Hlt2Phys_TOS);
   fChain->SetBranchAddress("muminus_TRACK_Type", &muminus_TRACK_Type, &b_muminus_TRACK_Type);
   fChain->SetBranchAddress("muminus_TRACK_Key", &muminus_TRACK_Key, &b_muminus_TRACK_Key);
   fChain->SetBranchAddress("muminus_TRACK_CHI2NDOF", &muminus_TRACK_CHI2NDOF, &b_muminus_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("muminus_TRACK_PCHI2", &muminus_TRACK_PCHI2, &b_muminus_TRACK_PCHI2);
   fChain->SetBranchAddress("muminus_TRACK_MatchCHI2", &muminus_TRACK_MatchCHI2, &b_muminus_TRACK_MatchCHI2);
   fChain->SetBranchAddress("muminus_TRACK_GhostProb", &muminus_TRACK_GhostProb, &b_muminus_TRACK_GhostProb);
   fChain->SetBranchAddress("muminus_TRACK_CloneDist", &muminus_TRACK_CloneDist, &b_muminus_TRACK_CloneDist);
   fChain->SetBranchAddress("muminus_TRACK_Likelihood", &muminus_TRACK_Likelihood, &b_muminus_TRACK_Likelihood);
   fChain->SetBranchAddress("K_CosTheta", &K_CosTheta, &b_K_CosTheta);
   fChain->SetBranchAddress("K_MINIP", &K_MINIP, &b_K_MINIP);
   fChain->SetBranchAddress("K_MINIPCHI2", &K_MINIPCHI2, &b_K_MINIPCHI2);
   fChain->SetBranchAddress("K_MINIPNEXTBEST", &K_MINIPNEXTBEST, &b_K_MINIPNEXTBEST);
   fChain->SetBranchAddress("K_MINIPCHI2NEXTBEST", &K_MINIPCHI2NEXTBEST, &b_K_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("K_OWNPV_X", &K_OWNPV_X, &b_K_OWNPV_X);
   fChain->SetBranchAddress("K_OWNPV_Y", &K_OWNPV_Y, &b_K_OWNPV_Y);
   fChain->SetBranchAddress("K_OWNPV_Z", &K_OWNPV_Z, &b_K_OWNPV_Z);
   fChain->SetBranchAddress("K_OWNPV_XERR", &K_OWNPV_XERR, &b_K_OWNPV_XERR);
   fChain->SetBranchAddress("K_OWNPV_YERR", &K_OWNPV_YERR, &b_K_OWNPV_YERR);
   fChain->SetBranchAddress("K_OWNPV_ZERR", &K_OWNPV_ZERR, &b_K_OWNPV_ZERR);
   fChain->SetBranchAddress("K_OWNPV_CHI2", &K_OWNPV_CHI2, &b_K_OWNPV_CHI2);
   fChain->SetBranchAddress("K_OWNPV_NDOF", &K_OWNPV_NDOF, &b_K_OWNPV_NDOF);
   fChain->SetBranchAddress("K_OWNPV_COV_", K_OWNPV_COV_, &b_K_OWNPV_COV_);
   fChain->SetBranchAddress("K_IP_OWNPV", &K_IP_OWNPV, &b_K_IP_OWNPV);
   fChain->SetBranchAddress("K_IPCHI2_OWNPV", &K_IPCHI2_OWNPV, &b_K_IPCHI2_OWNPV);
   fChain->SetBranchAddress("K_TOPPV_X", &K_TOPPV_X, &b_K_TOPPV_X);
   fChain->SetBranchAddress("K_TOPPV_Y", &K_TOPPV_Y, &b_K_TOPPV_Y);
   fChain->SetBranchAddress("K_TOPPV_Z", &K_TOPPV_Z, &b_K_TOPPV_Z);
   fChain->SetBranchAddress("K_TOPPV_XERR", &K_TOPPV_XERR, &b_K_TOPPV_XERR);
   fChain->SetBranchAddress("K_TOPPV_YERR", &K_TOPPV_YERR, &b_K_TOPPV_YERR);
   fChain->SetBranchAddress("K_TOPPV_ZERR", &K_TOPPV_ZERR, &b_K_TOPPV_ZERR);
   fChain->SetBranchAddress("K_TOPPV_CHI2", &K_TOPPV_CHI2, &b_K_TOPPV_CHI2);
   fChain->SetBranchAddress("K_TOPPV_NDOF", &K_TOPPV_NDOF, &b_K_TOPPV_NDOF);
   fChain->SetBranchAddress("K_TOPPV_COV_", K_TOPPV_COV_, &b_K_TOPPV_COV_);
   fChain->SetBranchAddress("K_IP_TOPPV", &K_IP_TOPPV, &b_K_IP_TOPPV);
   fChain->SetBranchAddress("K_IPCHI2_TOPPV", &K_IPCHI2_TOPPV, &b_K_IPCHI2_TOPPV);
   fChain->SetBranchAddress("K_ORIVX_X", &K_ORIVX_X, &b_K_ORIVX_X);
   fChain->SetBranchAddress("K_ORIVX_Y", &K_ORIVX_Y, &b_K_ORIVX_Y);
   fChain->SetBranchAddress("K_ORIVX_Z", &K_ORIVX_Z, &b_K_ORIVX_Z);
   fChain->SetBranchAddress("K_ORIVX_XERR", &K_ORIVX_XERR, &b_K_ORIVX_XERR);
   fChain->SetBranchAddress("K_ORIVX_YERR", &K_ORIVX_YERR, &b_K_ORIVX_YERR);
   fChain->SetBranchAddress("K_ORIVX_ZERR", &K_ORIVX_ZERR, &b_K_ORIVX_ZERR);
   fChain->SetBranchAddress("K_ORIVX_CHI2", &K_ORIVX_CHI2, &b_K_ORIVX_CHI2);
   fChain->SetBranchAddress("K_ORIVX_NDOF", &K_ORIVX_NDOF, &b_K_ORIVX_NDOF);
   fChain->SetBranchAddress("K_ORIVX_COV_", K_ORIVX_COV_, &b_K_ORIVX_COV_);
   fChain->SetBranchAddress("K_IP_ORIVX", &K_IP_ORIVX, &b_K_IP_ORIVX);
   fChain->SetBranchAddress("K_IPCHI2_ORIVX", &K_IPCHI2_ORIVX, &b_K_IPCHI2_ORIVX);
   fChain->SetBranchAddress("K_P", &K_P, &b_K_P);
   fChain->SetBranchAddress("K_PT", &K_PT, &b_K_PT);
   fChain->SetBranchAddress("K_PE", &K_PE, &b_K_PE);
   fChain->SetBranchAddress("K_PX", &K_PX, &b_K_PX);
   fChain->SetBranchAddress("K_PY", &K_PY, &b_K_PY);
   fChain->SetBranchAddress("K_PZ", &K_PZ, &b_K_PZ);
   fChain->SetBranchAddress("K_M", &K_M, &b_K_M);
   fChain->SetBranchAddress("K_ID", &K_ID, &b_K_ID);
   fChain->SetBranchAddress("K_PIDe", &K_PIDe, &b_K_PIDe);
   fChain->SetBranchAddress("K_PIDmu", &K_PIDmu, &b_K_PIDmu);
   fChain->SetBranchAddress("K_PIDK", &K_PIDK, &b_K_PIDK);
   fChain->SetBranchAddress("K_PIDp", &K_PIDp, &b_K_PIDp);
   fChain->SetBranchAddress("K_ProbNNe", &K_ProbNNe, &b_K_ProbNNe);
   fChain->SetBranchAddress("K_ProbNNk", &K_ProbNNk, &b_K_ProbNNk);
   fChain->SetBranchAddress("K_ProbNNp", &K_ProbNNp, &b_K_ProbNNp);
   fChain->SetBranchAddress("K_ProbNNpi", &K_ProbNNpi, &b_K_ProbNNpi);
   fChain->SetBranchAddress("K_ProbNNmu", &K_ProbNNmu, &b_K_ProbNNmu);
   fChain->SetBranchAddress("K_ProbNNghost", &K_ProbNNghost, &b_K_ProbNNghost);
   fChain->SetBranchAddress("K_hasMuon", &K_hasMuon, &b_K_hasMuon);
   fChain->SetBranchAddress("K_isMuon", &K_isMuon, &b_K_isMuon);
   fChain->SetBranchAddress("K_hasRich", &K_hasRich, &b_K_hasRich);
   fChain->SetBranchAddress("K_UsedRichAerogel", &K_UsedRichAerogel, &b_K_UsedRichAerogel);
   fChain->SetBranchAddress("K_UsedRich1Gas", &K_UsedRich1Gas, &b_K_UsedRich1Gas);
   fChain->SetBranchAddress("K_UsedRich2Gas", &K_UsedRich2Gas, &b_K_UsedRich2Gas);
   fChain->SetBranchAddress("K_RichAboveElThres", &K_RichAboveElThres, &b_K_RichAboveElThres);
   fChain->SetBranchAddress("K_RichAboveMuThres", &K_RichAboveMuThres, &b_K_RichAboveMuThres);
   fChain->SetBranchAddress("K_RichAbovePiThres", &K_RichAbovePiThres, &b_K_RichAbovePiThres);
   fChain->SetBranchAddress("K_RichAboveKaThres", &K_RichAboveKaThres, &b_K_RichAboveKaThres);
   fChain->SetBranchAddress("K_RichAbovePrThres", &K_RichAbovePrThres, &b_K_RichAbovePrThres);
   fChain->SetBranchAddress("K_hasCalo", &K_hasCalo, &b_K_hasCalo);
   fChain->SetBranchAddress("K_L0Global_Dec", &K_L0Global_Dec, &b_K_L0Global_Dec);
   fChain->SetBranchAddress("K_L0Global_TIS", &K_L0Global_TIS, &b_K_L0Global_TIS);
   fChain->SetBranchAddress("K_L0Global_TOS", &K_L0Global_TOS, &b_K_L0Global_TOS);
   fChain->SetBranchAddress("K_Hlt1Global_Dec", &K_Hlt1Global_Dec, &b_K_Hlt1Global_Dec);
   fChain->SetBranchAddress("K_Hlt1Global_TIS", &K_Hlt1Global_TIS, &b_K_Hlt1Global_TIS);
   fChain->SetBranchAddress("K_Hlt1Global_TOS", &K_Hlt1Global_TOS, &b_K_Hlt1Global_TOS);
   fChain->SetBranchAddress("K_Hlt1Phys_Dec", &K_Hlt1Phys_Dec, &b_K_Hlt1Phys_Dec);
   fChain->SetBranchAddress("K_Hlt1Phys_TIS", &K_Hlt1Phys_TIS, &b_K_Hlt1Phys_TIS);
   fChain->SetBranchAddress("K_Hlt1Phys_TOS", &K_Hlt1Phys_TOS, &b_K_Hlt1Phys_TOS);
   fChain->SetBranchAddress("K_Hlt2Global_Dec", &K_Hlt2Global_Dec, &b_K_Hlt2Global_Dec);
   fChain->SetBranchAddress("K_Hlt2Global_TIS", &K_Hlt2Global_TIS, &b_K_Hlt2Global_TIS);
   fChain->SetBranchAddress("K_Hlt2Global_TOS", &K_Hlt2Global_TOS, &b_K_Hlt2Global_TOS);
   fChain->SetBranchAddress("K_Hlt2Phys_Dec", &K_Hlt2Phys_Dec, &b_K_Hlt2Phys_Dec);
   fChain->SetBranchAddress("K_Hlt2Phys_TIS", &K_Hlt2Phys_TIS, &b_K_Hlt2Phys_TIS);
   fChain->SetBranchAddress("K_Hlt2Phys_TOS", &K_Hlt2Phys_TOS, &b_K_Hlt2Phys_TOS);
   fChain->SetBranchAddress("K_TRACK_Type", &K_TRACK_Type, &b_K_TRACK_Type);
   fChain->SetBranchAddress("K_TRACK_Key", &K_TRACK_Key, &b_K_TRACK_Key);
   fChain->SetBranchAddress("K_TRACK_CHI2NDOF", &K_TRACK_CHI2NDOF, &b_K_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("K_TRACK_PCHI2", &K_TRACK_PCHI2, &b_K_TRACK_PCHI2);
   fChain->SetBranchAddress("K_TRACK_MatchCHI2", &K_TRACK_MatchCHI2, &b_K_TRACK_MatchCHI2);
   fChain->SetBranchAddress("K_TRACK_GhostProb", &K_TRACK_GhostProb, &b_K_TRACK_GhostProb);
   fChain->SetBranchAddress("K_TRACK_CloneDist", &K_TRACK_CloneDist, &b_K_TRACK_CloneDist);
   fChain->SetBranchAddress("K_TRACK_Likelihood", &K_TRACK_Likelihood, &b_K_TRACK_Likelihood);
   fChain->SetBranchAddress("K_ETA", &K_ETA, &b_K_ETA);
   fChain->SetBranchAddress("K_PHI", &K_PHI, &b_K_PHI);
   fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
   fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
   fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
   fChain->SetBranchAddress("nlong", &nlong, &b_nlong);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
   fChain->SetBranchAddress("BCType", &BCType, &b_BCType);
   fChain->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
   fChain->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
   fChain->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
   fChain->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
   fChain->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
   fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("PVX", PVX, &b_PVX);
   fChain->SetBranchAddress("PVY", PVY, &b_PVY);
   fChain->SetBranchAddress("PVZ", PVZ, &b_PVZ);
   fChain->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
   fChain->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
   fChain->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
   fChain->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
   fChain->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
   fChain->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
   fChain->SetBranchAddress("StrippingFullDSTDiMuonJpsi2MuMuDetachedLineDecision", &StrippingFullDSTDiMuonJpsi2MuMuDetachedLineDecision, &b_StrippingFullDSTDiMuonJpsi2MuMuDetachedLineDecision);
   Notify();
}

Bool_t cuts::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void cuts::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t cuts::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef cuts_cxx
