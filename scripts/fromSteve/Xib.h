//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Dec 11 09:26:53 2018 by ROOT version 5.34/18
// from TTree MyTuple/MyTuple
// found on file: Xib2JpsiXi_Run2_DD_L0JPsiTOS_NewStrip.root
//////////////////////////////////////////////////////////

#ifndef Xib_h
#define Xib_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxXib_ENDVERTEX_COV = 1;
const Int_t kMaxXib_OWNPV_COV = 1;
const Int_t kMaxJpsi_ENDVERTEX_COV = 1;
const Int_t kMaxJpsi_OWNPV_COV = 1;
const Int_t kMaxJpsi_ORIVX_COV = 1;
const Int_t kMaxmuplus_OWNPV_COV = 1;
const Int_t kMaxmuplus_ORIVX_COV = 1;
const Int_t kMaxmuminus_OWNPV_COV = 1;
const Int_t kMaxmuminus_ORIVX_COV = 1;
const Int_t kMaxXi_ENDVERTEX_COV = 1;
const Int_t kMaxXi_OWNPV_COV = 1;
const Int_t kMaxXi_ORIVX_COV = 1;
const Int_t kMaxL_ENDVERTEX_COV = 1;
const Int_t kMaxL_OWNPV_COV = 1;
const Int_t kMaxL_ORIVX_COV = 1;
const Int_t kMaxp_OWNPV_COV = 1;
const Int_t kMaxp_ORIVX_COV = 1;
const Int_t kMaxpi_OWNPV_COV = 1;
const Int_t kMaxpi_ORIVX_COV = 1;
const Int_t kMaxBachPi_OWNPV_COV = 1;
const Int_t kMaxBachPi_ORIVX_COV = 1;

class Xib {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        Xib_ENDVERTEX_X;
   Double_t        Xib_ENDVERTEX_Y;
   Double_t        Xib_ENDVERTEX_Z;
   Double_t        Xib_ENDVERTEX_XERR;
   Double_t        Xib_ENDVERTEX_YERR;
   Double_t        Xib_ENDVERTEX_ZERR;
   Double_t        Xib_ENDVERTEX_CHI2;
   Int_t           Xib_ENDVERTEX_NDOF;
   Float_t         Xib_ENDVERTEX_COV_[3][3];
   Double_t        Xib_OWNPV_X;
   Double_t        Xib_OWNPV_Y;
   Double_t        Xib_OWNPV_Z;
   Double_t        Xib_OWNPV_XERR;
   Double_t        Xib_OWNPV_YERR;
   Double_t        Xib_OWNPV_ZERR;
   Double_t        Xib_OWNPV_CHI2;
   Int_t           Xib_OWNPV_NDOF;
   Float_t         Xib_OWNPV_COV_[3][3];
   Double_t        Xib_IP_OWNPV;
   Double_t        Xib_IPCHI2_OWNPV;
   Double_t        Xib_FD_OWNPV;
   Double_t        Xib_FDCHI2_OWNPV;
   Double_t        Xib_DIRA_OWNPV;
   Double_t        Xib_P;
   Double_t        Xib_PT;
   Double_t        Xib_PE;
   Double_t        Xib_PX;
   Double_t        Xib_PY;
   Double_t        Xib_PZ;
   Double_t        Xib_MM;
   Double_t        Xib_MMERR;
   Double_t        Xib_M;
   Int_t           Xib_ID;
   Double_t        Xib_TAU;
   Double_t        Xib_TAUERR;
   Double_t        Xib_TAUCHI2;
   Bool_t          Xib_L0Global_Dec;
   Bool_t          Xib_L0Global_TIS;
   Bool_t          Xib_L0Global_TOS;
   Bool_t          Xib_Hlt1Global_Dec;
   Bool_t          Xib_Hlt1Global_TIS;
   Bool_t          Xib_Hlt1Global_TOS;
   Bool_t          Xib_Hlt1Phys_Dec;
   Bool_t          Xib_Hlt1Phys_TIS;
   Bool_t          Xib_Hlt1Phys_TOS;
   Bool_t          Xib_Hlt2Global_Dec;
   Bool_t          Xib_Hlt2Global_TIS;
   Bool_t          Xib_Hlt2Global_TOS;
   Bool_t          Xib_Hlt2Phys_Dec;
   Bool_t          Xib_Hlt2Phys_TIS;
   Bool_t          Xib_Hlt2Phys_TOS;
   Double_t        Xib_ETA;
   Double_t        Xib_PHI;
   Double_t        Xib_DTF_CTAU;
   Double_t        Xib_DTF_CTAUERR;
   Double_t        Xib_DTF_CTAUS;
   Double_t        Xib_DTF_M_JpsiXiLConstr;
   Double_t        Xib_DTF_VCHI2NDOF;
   Int_t           Added_n_Particles;
   Float_t         psi_1S_H_M[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_M_pion[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_M_kaon[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_M_proton[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PT[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_ETA[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PHI[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PX[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PY[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_FDCHI2_OLD[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_FDCHI2_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_FD_OLD[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_FD_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEXCHI2_OLD[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEXCHI2_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_X_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_Y_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_Z_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_X_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_Y_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_Z_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_COV_XX[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_COV_XX[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_COV_XY[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_COV_XY[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_COV_XZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_COV_XZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_COV_YY[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_COV_YY[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_COV_YZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_COV_YZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VERTEX_COV_ZZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_PV_COV_ZZ[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_VertexRefitStatus[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_IPCHI2_OLD[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_IP_OLD[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_IPCHI2_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_IP_NEW[200];   //[Added_n_Particles]
   Float_t         psi_1S_H_MINIPCHI2[200];   //[Added_n_Particles]
   Float_t         Xib_H_OPENING[200];   //[Added_n_Particles]
   Float_t         Added_H_PX[200];   //[Added_n_Particles]
   Float_t         Added_H_PY[200];   //[Added_n_Particles]
   Float_t         Added_H_PZ[200];   //[Added_n_Particles]
   Float_t         Added_H_PT[200];   //[Added_n_Particles]
   Float_t         Added_H_ETA[200];   //[Added_n_Particles]
   Float_t         Added_H_PHI[200];   //[Added_n_Particles]
   Float_t         Added_H_THETA[200];   //[Added_n_Particles]
   Float_t         Added_H_DELTAPHI[200];   //[Added_n_Particles]
   Float_t         Added_H_PID[200];   //[Added_n_Particles]
   Int_t           Added_n_Tracks;
   Float_t         Added_H_PROBNNPID[200];   //[Added_n_Tracks]
   Float_t         Added_H_GHOST[200];   //[Added_n_Tracks]
   Float_t         Added_H_TRACKCHI2[200];   //[Added_n_Tracks]
   Float_t         Added_H_TRACKTYPE[200];   //[Added_n_Tracks]
   Float_t         Added_H_PIDe[200];   //[Added_n_Tracks]
   Float_t         Added_H_PIDK[200];   //[Added_n_Tracks]
   Float_t         Added_H_PIDmu[200];   //[Added_n_Tracks]
   Float_t         Added_H_PIDp[200];   //[Added_n_Tracks]
   Float_t         Added_H_ProbNNpi[200];   //[Added_n_Tracks]
   Float_t         Added_H_ProbNNk[200];   //[Added_n_Tracks]
   Float_t         Added_H_ProbNNp[200];   //[Added_n_Tracks]
   Float_t         Added_H_ProbNNe[200];   //[Added_n_Tracks]
   Float_t         Added_H_ProbNNmu[200];   //[Added_n_Tracks]
   Int_t           Xib_ConsXib_nPV;
   Float_t         Xib_ConsXib_J_psi_1S_M[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_MERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_P[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_PERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_ctau[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_ctauErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_decayLength[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_decayLengthErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_0_ID[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_0_PE[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_0_PX[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_0_PY[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_0_PZ[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_ID[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_PE[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_PX[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_PY[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_J_psi_1S_muminus_PZ[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_M[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_MERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_P[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_PERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_ctau[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_ctauErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_decayLength[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Lambda0_decayLengthErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_M[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_MERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_P[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_PERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_PV_X[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_PV_Y[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_PV_Z[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_PV_key[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_piplus_ID[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_piplus_PE[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_piplus_PX[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_piplus_PY[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_piplus_PZ[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_pplus_ID[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_pplus_PE[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_pplus_PX[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_pplus_PY[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_Lambda0_pplus_PZ[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_M[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_MERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_P[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_PERR[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_ctau[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_ctauErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_decayLength[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_decayLengthErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_piplus_ID[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_piplus_PE[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_piplus_PX[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_piplus_PY[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_Ximinus_piplus_PZ[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_chi2[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_ctau[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_ctauErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_decayLength[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_decayLengthErr[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_nDOF[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_nIter[100];   //[Xib_ConsXib_nPV]
   Float_t         Xib_ConsXib_status[100];   //[Xib_ConsXib_nPV]
   Bool_t          Xib_L0MuonDecision_Dec;
   Bool_t          Xib_L0MuonDecision_TIS;
   Bool_t          Xib_L0MuonDecision_TOS;
   Bool_t          Xib_L0DiMuonDecision_Dec;
   Bool_t          Xib_L0DiMuonDecision_TIS;
   Bool_t          Xib_L0DiMuonDecision_TOS;
   Bool_t          Xib_Hlt1SingleMuonNoIPDecision_Dec;
   Bool_t          Xib_Hlt1SingleMuonNoIPDecision_TIS;
   Bool_t          Xib_Hlt1SingleMuonNoIPDecision_TOS;
   Bool_t          Xib_Hlt1DiMuonHighMassDecision_Dec;
   Bool_t          Xib_Hlt1DiMuonHighMassDecision_TIS;
   Bool_t          Xib_Hlt1DiMuonHighMassDecision_TOS;
   Bool_t          Xib_Hlt1TrackAllL0Decision_Dec;
   Bool_t          Xib_Hlt1TrackAllL0Decision_TIS;
   Bool_t          Xib_Hlt1TrackAllL0Decision_TOS;
   Bool_t          Xib_Hlt1TrackMuonDecision_Dec;
   Bool_t          Xib_Hlt1TrackMuonDecision_TIS;
   Bool_t          Xib_Hlt1TrackMuonDecision_TOS;
   Bool_t          Xib_Hlt2DiMuonDetachedHeavyDecision_Dec;
   Bool_t          Xib_Hlt2DiMuonDetachedHeavyDecision_TIS;
   Bool_t          Xib_Hlt2DiMuonDetachedHeavyDecision_TOS;
   Bool_t          Xib_Hlt2DiMuonDetachedJPsiDecision_Dec;
   Bool_t          Xib_Hlt2DiMuonDetachedJPsiDecision_TIS;
   Bool_t          Xib_Hlt2DiMuonDetachedJPsiDecision_TOS;
   Double_t        Jpsi_CosTheta;
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
   Double_t        Jpsi_ORIVX_X;
   Double_t        Jpsi_ORIVX_Y;
   Double_t        Jpsi_ORIVX_Z;
   Double_t        Jpsi_ORIVX_XERR;
   Double_t        Jpsi_ORIVX_YERR;
   Double_t        Jpsi_ORIVX_ZERR;
   Double_t        Jpsi_ORIVX_CHI2;
   Int_t           Jpsi_ORIVX_NDOF;
   Float_t         Jpsi_ORIVX_COV_[3][3];
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
   Double_t        muplus_ORIVX_X;
   Double_t        muplus_ORIVX_Y;
   Double_t        muplus_ORIVX_Z;
   Double_t        muplus_ORIVX_XERR;
   Double_t        muplus_ORIVX_YERR;
   Double_t        muplus_ORIVX_ZERR;
   Double_t        muplus_ORIVX_CHI2;
   Int_t           muplus_ORIVX_NDOF;
   Float_t         muplus_ORIVX_COV_[3][3];
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
   Double_t        muminus_ORIVX_X;
   Double_t        muminus_ORIVX_Y;
   Double_t        muminus_ORIVX_Z;
   Double_t        muminus_ORIVX_XERR;
   Double_t        muminus_ORIVX_YERR;
   Double_t        muminus_ORIVX_ZERR;
   Double_t        muminus_ORIVX_CHI2;
   Int_t           muminus_ORIVX_NDOF;
   Float_t         muminus_ORIVX_COV_[3][3];
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
   Double_t        Xi_CosTheta;
   Double_t        Xi_ENDVERTEX_X;
   Double_t        Xi_ENDVERTEX_Y;
   Double_t        Xi_ENDVERTEX_Z;
   Double_t        Xi_ENDVERTEX_XERR;
   Double_t        Xi_ENDVERTEX_YERR;
   Double_t        Xi_ENDVERTEX_ZERR;
   Double_t        Xi_ENDVERTEX_CHI2;
   Int_t           Xi_ENDVERTEX_NDOF;
   Float_t         Xi_ENDVERTEX_COV_[3][3];
   Double_t        Xi_OWNPV_X;
   Double_t        Xi_OWNPV_Y;
   Double_t        Xi_OWNPV_Z;
   Double_t        Xi_OWNPV_XERR;
   Double_t        Xi_OWNPV_YERR;
   Double_t        Xi_OWNPV_ZERR;
   Double_t        Xi_OWNPV_CHI2;
   Int_t           Xi_OWNPV_NDOF;
   Float_t         Xi_OWNPV_COV_[3][3];
   Double_t        Xi_IP_OWNPV;
   Double_t        Xi_IPCHI2_OWNPV;
   Double_t        Xi_FD_OWNPV;
   Double_t        Xi_FDCHI2_OWNPV;
   Double_t        Xi_DIRA_OWNPV;
   Double_t        Xi_ORIVX_X;
   Double_t        Xi_ORIVX_Y;
   Double_t        Xi_ORIVX_Z;
   Double_t        Xi_ORIVX_XERR;
   Double_t        Xi_ORIVX_YERR;
   Double_t        Xi_ORIVX_ZERR;
   Double_t        Xi_ORIVX_CHI2;
   Int_t           Xi_ORIVX_NDOF;
   Float_t         Xi_ORIVX_COV_[3][3];
   Double_t        Xi_FD_ORIVX;
   Double_t        Xi_FDCHI2_ORIVX;
   Double_t        Xi_DIRA_ORIVX;
   Double_t        Xi_P;
   Double_t        Xi_PT;
   Double_t        Xi_PE;
   Double_t        Xi_PX;
   Double_t        Xi_PY;
   Double_t        Xi_PZ;
   Double_t        Xi_MM;
   Double_t        Xi_MMERR;
   Double_t        Xi_M;
   Int_t           Xi_ID;
   Double_t        Xi_TAU;
   Double_t        Xi_TAUERR;
   Double_t        Xi_TAUCHI2;
   Bool_t          Xi_L0Global_Dec;
   Bool_t          Xi_L0Global_TIS;
   Bool_t          Xi_L0Global_TOS;
   Bool_t          Xi_Hlt1Global_Dec;
   Bool_t          Xi_Hlt1Global_TIS;
   Bool_t          Xi_Hlt1Global_TOS;
   Bool_t          Xi_Hlt1Phys_Dec;
   Bool_t          Xi_Hlt1Phys_TIS;
   Bool_t          Xi_Hlt1Phys_TOS;
   Bool_t          Xi_Hlt2Global_Dec;
   Bool_t          Xi_Hlt2Global_TIS;
   Bool_t          Xi_Hlt2Global_TOS;
   Bool_t          Xi_Hlt2Phys_Dec;
   Bool_t          Xi_Hlt2Phys_TIS;
   Bool_t          Xi_Hlt2Phys_TOS;
   Double_t        Xi_ETA;
   Double_t        Xi_PHI;
   Double_t        Xi_DTF_L_ARMEN_JpsiXiConstr;
   Double_t        Xi_DTF_L_ARMEN_JpsiXiLConstr;
   Double_t        Xi_DTF_L_LV01_JpsiXiConstr;
   Double_t        Xi_DTF_L_LV01_JpsiXiLConstr;
   Double_t        Xi_DTF_L_WMpipi_JpsiXiConstr;
   Double_t        Xi_DTF_L_WMpipi_JpsiXiLConstr;
   Double_t        L_CosTheta;
   Double_t        L_ENDVERTEX_X;
   Double_t        L_ENDVERTEX_Y;
   Double_t        L_ENDVERTEX_Z;
   Double_t        L_ENDVERTEX_XERR;
   Double_t        L_ENDVERTEX_YERR;
   Double_t        L_ENDVERTEX_ZERR;
   Double_t        L_ENDVERTEX_CHI2;
   Int_t           L_ENDVERTEX_NDOF;
   Float_t         L_ENDVERTEX_COV_[3][3];
   Double_t        L_OWNPV_X;
   Double_t        L_OWNPV_Y;
   Double_t        L_OWNPV_Z;
   Double_t        L_OWNPV_XERR;
   Double_t        L_OWNPV_YERR;
   Double_t        L_OWNPV_ZERR;
   Double_t        L_OWNPV_CHI2;
   Int_t           L_OWNPV_NDOF;
   Float_t         L_OWNPV_COV_[3][3];
   Double_t        L_IP_OWNPV;
   Double_t        L_IPCHI2_OWNPV;
   Double_t        L_FD_OWNPV;
   Double_t        L_FDCHI2_OWNPV;
   Double_t        L_DIRA_OWNPV;
   Double_t        L_ORIVX_X;
   Double_t        L_ORIVX_Y;
   Double_t        L_ORIVX_Z;
   Double_t        L_ORIVX_XERR;
   Double_t        L_ORIVX_YERR;
   Double_t        L_ORIVX_ZERR;
   Double_t        L_ORIVX_CHI2;
   Int_t           L_ORIVX_NDOF;
   Float_t         L_ORIVX_COV_[3][3];
   Double_t        L_FD_ORIVX;
   Double_t        L_FDCHI2_ORIVX;
   Double_t        L_DIRA_ORIVX;
   Double_t        L_P;
   Double_t        L_PT;
   Double_t        L_PE;
   Double_t        L_PX;
   Double_t        L_PY;
   Double_t        L_PZ;
   Double_t        L_MM;
   Double_t        L_MMERR;
   Double_t        L_M;
   Int_t           L_ID;
   Double_t        L_TAU;
   Double_t        L_TAUERR;
   Double_t        L_TAUCHI2;
   Bool_t          L_L0Global_Dec;
   Bool_t          L_L0Global_TIS;
   Bool_t          L_L0Global_TOS;
   Bool_t          L_Hlt1Global_Dec;
   Bool_t          L_Hlt1Global_TIS;
   Bool_t          L_Hlt1Global_TOS;
   Bool_t          L_Hlt1Phys_Dec;
   Bool_t          L_Hlt1Phys_TIS;
   Bool_t          L_Hlt1Phys_TOS;
   Bool_t          L_Hlt2Global_Dec;
   Bool_t          L_Hlt2Global_TIS;
   Bool_t          L_Hlt2Global_TOS;
   Bool_t          L_Hlt2Phys_Dec;
   Bool_t          L_Hlt2Phys_TIS;
   Bool_t          L_Hlt2Phys_TOS;
   Double_t        L_ETA;
   Double_t        L_PHI;
   Double_t        L_ARMEN;
   Double_t        L_BPVLTIME;
   Double_t        L_WMkk;
   Double_t        L_WMkpi;
   Double_t        L_WMpik;
   Double_t        L_WMpipi;
   Double_t        L_WMpp;
   Double_t        L_dm;
   Double_t        L_lv01;
   Double_t        L_lv02;
   Double_t        p_CosTheta;
   Double_t        p_OWNPV_X;
   Double_t        p_OWNPV_Y;
   Double_t        p_OWNPV_Z;
   Double_t        p_OWNPV_XERR;
   Double_t        p_OWNPV_YERR;
   Double_t        p_OWNPV_ZERR;
   Double_t        p_OWNPV_CHI2;
   Int_t           p_OWNPV_NDOF;
   Float_t         p_OWNPV_COV_[3][3];
   Double_t        p_IP_OWNPV;
   Double_t        p_IPCHI2_OWNPV;
   Double_t        p_ORIVX_X;
   Double_t        p_ORIVX_Y;
   Double_t        p_ORIVX_Z;
   Double_t        p_ORIVX_XERR;
   Double_t        p_ORIVX_YERR;
   Double_t        p_ORIVX_ZERR;
   Double_t        p_ORIVX_CHI2;
   Int_t           p_ORIVX_NDOF;
   Float_t         p_ORIVX_COV_[3][3];
   Double_t        p_P;
   Double_t        p_PT;
   Double_t        p_PE;
   Double_t        p_PX;
   Double_t        p_PY;
   Double_t        p_PZ;
   Double_t        p_M;
   Int_t           p_ID;
   Double_t        p_PIDe;
   Double_t        p_PIDmu;
   Double_t        p_PIDK;
   Double_t        p_PIDp;
   Double_t        p_ProbNNe;
   Double_t        p_ProbNNk;
   Double_t        p_ProbNNp;
   Double_t        p_ProbNNpi;
   Double_t        p_ProbNNmu;
   Double_t        p_ProbNNghost;
   Bool_t          p_hasMuon;
   Bool_t          p_isMuon;
   Bool_t          p_hasRich;
   Bool_t          p_UsedRichAerogel;
   Bool_t          p_UsedRich1Gas;
   Bool_t          p_UsedRich2Gas;
   Bool_t          p_RichAboveElThres;
   Bool_t          p_RichAboveMuThres;
   Bool_t          p_RichAbovePiThres;
   Bool_t          p_RichAboveKaThres;
   Bool_t          p_RichAbovePrThres;
   Bool_t          p_hasCalo;
   Bool_t          p_L0Global_Dec;
   Bool_t          p_L0Global_TIS;
   Bool_t          p_L0Global_TOS;
   Bool_t          p_Hlt1Global_Dec;
   Bool_t          p_Hlt1Global_TIS;
   Bool_t          p_Hlt1Global_TOS;
   Bool_t          p_Hlt1Phys_Dec;
   Bool_t          p_Hlt1Phys_TIS;
   Bool_t          p_Hlt1Phys_TOS;
   Bool_t          p_Hlt2Global_Dec;
   Bool_t          p_Hlt2Global_TIS;
   Bool_t          p_Hlt2Global_TOS;
   Bool_t          p_Hlt2Phys_Dec;
   Bool_t          p_Hlt2Phys_TIS;
   Bool_t          p_Hlt2Phys_TOS;
   Int_t           p_TRACK_Type;
   Int_t           p_TRACK_Key;
   Double_t        p_TRACK_CHI2NDOF;
   Double_t        p_TRACK_PCHI2;
   Double_t        p_TRACK_MatchCHI2;
   Double_t        p_TRACK_GhostProb;
   Double_t        p_TRACK_CloneDist;
   Double_t        p_TRACK_Likelihood;
   Double_t        p_ETA;
   Double_t        p_PHI;
   Double_t        pi_CosTheta;
   Double_t        pi_OWNPV_X;
   Double_t        pi_OWNPV_Y;
   Double_t        pi_OWNPV_Z;
   Double_t        pi_OWNPV_XERR;
   Double_t        pi_OWNPV_YERR;
   Double_t        pi_OWNPV_ZERR;
   Double_t        pi_OWNPV_CHI2;
   Int_t           pi_OWNPV_NDOF;
   Float_t         pi_OWNPV_COV_[3][3];
   Double_t        pi_IP_OWNPV;
   Double_t        pi_IPCHI2_OWNPV;
   Double_t        pi_ORIVX_X;
   Double_t        pi_ORIVX_Y;
   Double_t        pi_ORIVX_Z;
   Double_t        pi_ORIVX_XERR;
   Double_t        pi_ORIVX_YERR;
   Double_t        pi_ORIVX_ZERR;
   Double_t        pi_ORIVX_CHI2;
   Int_t           pi_ORIVX_NDOF;
   Float_t         pi_ORIVX_COV_[3][3];
   Double_t        pi_P;
   Double_t        pi_PT;
   Double_t        pi_PE;
   Double_t        pi_PX;
   Double_t        pi_PY;
   Double_t        pi_PZ;
   Double_t        pi_M;
   Int_t           pi_ID;
   Double_t        pi_PIDe;
   Double_t        pi_PIDmu;
   Double_t        pi_PIDK;
   Double_t        pi_PIDp;
   Double_t        pi_ProbNNe;
   Double_t        pi_ProbNNk;
   Double_t        pi_ProbNNp;
   Double_t        pi_ProbNNpi;
   Double_t        pi_ProbNNmu;
   Double_t        pi_ProbNNghost;
   Bool_t          pi_hasMuon;
   Bool_t          pi_isMuon;
   Bool_t          pi_hasRich;
   Bool_t          pi_UsedRichAerogel;
   Bool_t          pi_UsedRich1Gas;
   Bool_t          pi_UsedRich2Gas;
   Bool_t          pi_RichAboveElThres;
   Bool_t          pi_RichAboveMuThres;
   Bool_t          pi_RichAbovePiThres;
   Bool_t          pi_RichAboveKaThres;
   Bool_t          pi_RichAbovePrThres;
   Bool_t          pi_hasCalo;
   Bool_t          pi_L0Global_Dec;
   Bool_t          pi_L0Global_TIS;
   Bool_t          pi_L0Global_TOS;
   Bool_t          pi_Hlt1Global_Dec;
   Bool_t          pi_Hlt1Global_TIS;
   Bool_t          pi_Hlt1Global_TOS;
   Bool_t          pi_Hlt1Phys_Dec;
   Bool_t          pi_Hlt1Phys_TIS;
   Bool_t          pi_Hlt1Phys_TOS;
   Bool_t          pi_Hlt2Global_Dec;
   Bool_t          pi_Hlt2Global_TIS;
   Bool_t          pi_Hlt2Global_TOS;
   Bool_t          pi_Hlt2Phys_Dec;
   Bool_t          pi_Hlt2Phys_TIS;
   Bool_t          pi_Hlt2Phys_TOS;
   Int_t           pi_TRACK_Type;
   Int_t           pi_TRACK_Key;
   Double_t        pi_TRACK_CHI2NDOF;
   Double_t        pi_TRACK_PCHI2;
   Double_t        pi_TRACK_MatchCHI2;
   Double_t        pi_TRACK_GhostProb;
   Double_t        pi_TRACK_CloneDist;
   Double_t        pi_TRACK_Likelihood;
   Double_t        pi_ETA;
   Double_t        pi_PHI;
   Double_t        BachPi_CosTheta;
   Double_t        BachPi_OWNPV_X;
   Double_t        BachPi_OWNPV_Y;
   Double_t        BachPi_OWNPV_Z;
   Double_t        BachPi_OWNPV_XERR;
   Double_t        BachPi_OWNPV_YERR;
   Double_t        BachPi_OWNPV_ZERR;
   Double_t        BachPi_OWNPV_CHI2;
   Int_t           BachPi_OWNPV_NDOF;
   Float_t         BachPi_OWNPV_COV_[3][3];
   Double_t        BachPi_IP_OWNPV;
   Double_t        BachPi_IPCHI2_OWNPV;
   Double_t        BachPi_ORIVX_X;
   Double_t        BachPi_ORIVX_Y;
   Double_t        BachPi_ORIVX_Z;
   Double_t        BachPi_ORIVX_XERR;
   Double_t        BachPi_ORIVX_YERR;
   Double_t        BachPi_ORIVX_ZERR;
   Double_t        BachPi_ORIVX_CHI2;
   Int_t           BachPi_ORIVX_NDOF;
   Float_t         BachPi_ORIVX_COV_[3][3];
   Double_t        BachPi_P;
   Double_t        BachPi_PT;
   Double_t        BachPi_PE;
   Double_t        BachPi_PX;
   Double_t        BachPi_PY;
   Double_t        BachPi_PZ;
   Double_t        BachPi_M;
   Int_t           BachPi_ID;
   Double_t        BachPi_PIDe;
   Double_t        BachPi_PIDmu;
   Double_t        BachPi_PIDK;
   Double_t        BachPi_PIDp;
   Double_t        BachPi_ProbNNe;
   Double_t        BachPi_ProbNNk;
   Double_t        BachPi_ProbNNp;
   Double_t        BachPi_ProbNNpi;
   Double_t        BachPi_ProbNNmu;
   Double_t        BachPi_ProbNNghost;
   Bool_t          BachPi_hasMuon;
   Bool_t          BachPi_isMuon;
   Bool_t          BachPi_hasRich;
   Bool_t          BachPi_UsedRichAerogel;
   Bool_t          BachPi_UsedRich1Gas;
   Bool_t          BachPi_UsedRich2Gas;
   Bool_t          BachPi_RichAboveElThres;
   Bool_t          BachPi_RichAboveMuThres;
   Bool_t          BachPi_RichAbovePiThres;
   Bool_t          BachPi_RichAboveKaThres;
   Bool_t          BachPi_RichAbovePrThres;
   Bool_t          BachPi_hasCalo;
   Bool_t          BachPi_L0Global_Dec;
   Bool_t          BachPi_L0Global_TIS;
   Bool_t          BachPi_L0Global_TOS;
   Bool_t          BachPi_Hlt1Global_Dec;
   Bool_t          BachPi_Hlt1Global_TIS;
   Bool_t          BachPi_Hlt1Global_TOS;
   Bool_t          BachPi_Hlt1Phys_Dec;
   Bool_t          BachPi_Hlt1Phys_TIS;
   Bool_t          BachPi_Hlt1Phys_TOS;
   Bool_t          BachPi_Hlt2Global_Dec;
   Bool_t          BachPi_Hlt2Global_TIS;
   Bool_t          BachPi_Hlt2Global_TOS;
   Bool_t          BachPi_Hlt2Phys_Dec;
   Bool_t          BachPi_Hlt2Phys_TIS;
   Bool_t          BachPi_Hlt2Phys_TOS;
   Int_t           BachPi_TRACK_Type;
   Int_t           BachPi_TRACK_Key;
   Double_t        BachPi_TRACK_CHI2NDOF;
   Double_t        BachPi_TRACK_PCHI2;
   Double_t        BachPi_TRACK_MatchCHI2;
   Double_t        BachPi_TRACK_GhostProb;
   Double_t        BachPi_TRACK_CloneDist;
   Double_t        BachPi_TRACK_Likelihood;
   Double_t        BachPi_ETA;
   Double_t        BachPi_PHI;
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
   TBranch        *b_Xib_ENDVERTEX_X;   //!
   TBranch        *b_Xib_ENDVERTEX_Y;   //!
   TBranch        *b_Xib_ENDVERTEX_Z;   //!
   TBranch        *b_Xib_ENDVERTEX_XERR;   //!
   TBranch        *b_Xib_ENDVERTEX_YERR;   //!
   TBranch        *b_Xib_ENDVERTEX_ZERR;   //!
   TBranch        *b_Xib_ENDVERTEX_CHI2;   //!
   TBranch        *b_Xib_ENDVERTEX_NDOF;   //!
   TBranch        *b_Xib_ENDVERTEX_COV_;   //!
   TBranch        *b_Xib_OWNPV_X;   //!
   TBranch        *b_Xib_OWNPV_Y;   //!
   TBranch        *b_Xib_OWNPV_Z;   //!
   TBranch        *b_Xib_OWNPV_XERR;   //!
   TBranch        *b_Xib_OWNPV_YERR;   //!
   TBranch        *b_Xib_OWNPV_ZERR;   //!
   TBranch        *b_Xib_OWNPV_CHI2;   //!
   TBranch        *b_Xib_OWNPV_NDOF;   //!
   TBranch        *b_Xib_OWNPV_COV_;   //!
   TBranch        *b_Xib_IP_OWNPV;   //!
   TBranch        *b_Xib_IPCHI2_OWNPV;   //!
   TBranch        *b_Xib_FD_OWNPV;   //!
   TBranch        *b_Xib_FDCHI2_OWNPV;   //!
   TBranch        *b_Xib_DIRA_OWNPV;   //!
   TBranch        *b_Xib_P;   //!
   TBranch        *b_Xib_PT;   //!
   TBranch        *b_Xib_PE;   //!
   TBranch        *b_Xib_PX;   //!
   TBranch        *b_Xib_PY;   //!
   TBranch        *b_Xib_PZ;   //!
   TBranch        *b_Xib_MM;   //!
   TBranch        *b_Xib_MMERR;   //!
   TBranch        *b_Xib_M;   //!
   TBranch        *b_Xib_ID;   //!
   TBranch        *b_Xib_TAU;   //!
   TBranch        *b_Xib_TAUERR;   //!
   TBranch        *b_Xib_TAUCHI2;   //!
   TBranch        *b_Xib_L0Global_Dec;   //!
   TBranch        *b_Xib_L0Global_TIS;   //!
   TBranch        *b_Xib_L0Global_TOS;   //!
   TBranch        *b_Xib_Hlt1Global_Dec;   //!
   TBranch        *b_Xib_Hlt1Global_TIS;   //!
   TBranch        *b_Xib_Hlt1Global_TOS;   //!
   TBranch        *b_Xib_Hlt1Phys_Dec;   //!
   TBranch        *b_Xib_Hlt1Phys_TIS;   //!
   TBranch        *b_Xib_Hlt1Phys_TOS;   //!
   TBranch        *b_Xib_Hlt2Global_Dec;   //!
   TBranch        *b_Xib_Hlt2Global_TIS;   //!
   TBranch        *b_Xib_Hlt2Global_TOS;   //!
   TBranch        *b_Xib_Hlt2Phys_Dec;   //!
   TBranch        *b_Xib_Hlt2Phys_TIS;   //!
   TBranch        *b_Xib_Hlt2Phys_TOS;   //!
   TBranch        *b_Xib_ETA;   //!
   TBranch        *b_Xib_PHI;   //!
   TBranch        *b_Xib_DTF_CTAU;   //!
   TBranch        *b_Xib_DTF_CTAUERR;   //!
   TBranch        *b_Xib_DTF_CTAUS;   //!
   TBranch        *b_Xib_DTF_M_JpsiXiLConstr;   //!
   TBranch        *b_Xib_DTF_VCHI2NDOF;   //!
   TBranch        *b_Added_n_Particles;   //!
   TBranch        *b_psi_1S_H_M;   //!
   TBranch        *b_psi_1S_H_M_pion;   //!
   TBranch        *b_psi_1S_H_M_kaon;   //!
   TBranch        *b_psi_1S_H_M_proton;   //!
   TBranch        *b_psi_1S_H_PT;   //!
   TBranch        *b_psi_1S_H_ETA;   //!
   TBranch        *b_psi_1S_H_PHI;   //!
   TBranch        *b_psi_1S_H_PX;   //!
   TBranch        *b_psi_1S_H_PY;   //!
   TBranch        *b_psi_1S_H_PZ;   //!
   TBranch        *b_psi_1S_H_FDCHI2_OLD;   //!
   TBranch        *b_psi_1S_H_FDCHI2_NEW;   //!
   TBranch        *b_psi_1S_H_FD_OLD;   //!
   TBranch        *b_psi_1S_H_FD_NEW;   //!
   TBranch        *b_psi_1S_H_VERTEXCHI2_OLD;   //!
   TBranch        *b_psi_1S_H_VERTEXCHI2_NEW;   //!
   TBranch        *b_psi_1S_H_VERTEX_X_NEW;   //!
   TBranch        *b_psi_1S_H_VERTEX_Y_NEW;   //!
   TBranch        *b_psi_1S_H_VERTEX_Z_NEW;   //!
   TBranch        *b_psi_1S_H_PV_X_NEW;   //!
   TBranch        *b_psi_1S_H_PV_Y_NEW;   //!
   TBranch        *b_psi_1S_H_PV_Z_NEW;   //!
   TBranch        *b_psi_1S_H_VERTEX_COV_XX;   //!
   TBranch        *b_psi_1S_H_PV_COV_XX;   //!
   TBranch        *b_psi_1S_H_VERTEX_COV_XY;   //!
   TBranch        *b_psi_1S_H_PV_COV_XY;   //!
   TBranch        *b_psi_1S_H_VERTEX_COV_XZ;   //!
   TBranch        *b_psi_1S_H_PV_COV_XZ;   //!
   TBranch        *b_psi_1S_H_VERTEX_COV_YY;   //!
   TBranch        *b_psi_1S_H_PV_COV_YY;   //!
   TBranch        *b_psi_1S_H_VERTEX_COV_YZ;   //!
   TBranch        *b_psi_1S_H_PV_COV_YZ;   //!
   TBranch        *b_psi_1S_H_VERTEX_COV_ZZ;   //!
   TBranch        *b_psi_1S_H_PV_COV_ZZ;   //!
   TBranch        *b_psi_1S_H_VertexRefitStatus;   //!
   TBranch        *b_psi_1S_H_IPCHI2_OLD;   //!
   TBranch        *b_psi_1S_H_IP_OLD;   //!
   TBranch        *b_psi_1S_H_IPCHI2_NEW;   //!
   TBranch        *b_psi_1S_H_IP_NEW;   //!
   TBranch        *b_psi_1S_H_MINIPCHI2;   //!
   TBranch        *b_Xib_H_OPENING;   //!
   TBranch        *b_Added_H_PX;   //!
   TBranch        *b_Added_H_PY;   //!
   TBranch        *b_Added_H_PZ;   //!
   TBranch        *b_Added_H_PT;   //!
   TBranch        *b_Added_H_ETA;   //!
   TBranch        *b_Added_H_PHI;   //!
   TBranch        *b_Added_H_THETA;   //!
   TBranch        *b_Added_H_DELTAPHI;   //!
   TBranch        *b_Added_H_PID;   //!
   TBranch        *b_Added_n_Tracks;   //!
   TBranch        *b_Added_H_PROBNNPID;   //!
   TBranch        *b_Added_H_GHOST;   //!
   TBranch        *b_Added_H_TRACKCHI2;   //!
   TBranch        *b_Added_H_TRACKTYPE;   //!
   TBranch        *b_Added_H_PIDe;   //!
   TBranch        *b_Added_H_PIDK;   //!
   TBranch        *b_Added_H_PIDmu;   //!
   TBranch        *b_Added_H_PIDp;   //!
   TBranch        *b_Added_H_ProbNNpi;   //!
   TBranch        *b_Added_H_ProbNNk;   //!
   TBranch        *b_Added_H_ProbNNp;   //!
   TBranch        *b_Added_H_ProbNNe;   //!
   TBranch        *b_Added_H_ProbNNmu;   //!
   TBranch        *b_Xib_ConsXib_nPV;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_M;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_MERR;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_P;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_PERR;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_ctau;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_ctauErr;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_decayLength;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_decayLengthErr;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_0_ID;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_0_PE;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_0_PX;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_0_PY;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_0_PZ;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_ID;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_PE;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_PX;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_PY;   //!
   TBranch        *b_Xib_ConsXib_J_psi_1S_muminus_PZ;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_M;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_MERR;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_P;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_PERR;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_ctau;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_ctauErr;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_decayLength;   //!
   TBranch        *b_Xib_ConsXib_Lambda0_decayLengthErr;   //!
   TBranch        *b_Xib_ConsXib_M;   //!
   TBranch        *b_Xib_ConsXib_MERR;   //!
   TBranch        *b_Xib_ConsXib_P;   //!
   TBranch        *b_Xib_ConsXib_PERR;   //!
   TBranch        *b_Xib_ConsXib_PV_X;   //!
   TBranch        *b_Xib_ConsXib_PV_Y;   //!
   TBranch        *b_Xib_ConsXib_PV_Z;   //!
   TBranch        *b_Xib_ConsXib_PV_key;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_piplus_ID;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_piplus_PE;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_piplus_PX;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_piplus_PY;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_piplus_PZ;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_pplus_ID;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_pplus_PE;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_pplus_PX;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_pplus_PY;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_Lambda0_pplus_PZ;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_M;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_MERR;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_P;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_PERR;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_ctau;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_ctauErr;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_decayLength;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_decayLengthErr;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_piplus_ID;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_piplus_PE;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_piplus_PX;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_piplus_PY;   //!
   TBranch        *b_Xib_ConsXib_Ximinus_piplus_PZ;   //!
   TBranch        *b_Xib_ConsXib_chi2;   //!
   TBranch        *b_Xib_ConsXib_ctau;   //!
   TBranch        *b_Xib_ConsXib_ctauErr;   //!
   TBranch        *b_Xib_ConsXib_decayLength;   //!
   TBranch        *b_Xib_ConsXib_decayLengthErr;   //!
   TBranch        *b_Xib_ConsXib_nDOF;   //!
   TBranch        *b_Xib_ConsXib_nIter;   //!
   TBranch        *b_Xib_ConsXib_status;   //!
   TBranch        *b_Xib_L0MuonDecision_Dec;   //!
   TBranch        *b_Xib_L0MuonDecision_TIS;   //!
   TBranch        *b_Xib_L0MuonDecision_TOS;   //!
   TBranch        *b_Xib_L0DiMuonDecision_Dec;   //!
   TBranch        *b_Xib_L0DiMuonDecision_TIS;   //!
   TBranch        *b_Xib_L0DiMuonDecision_TOS;   //!
   TBranch        *b_Xib_Hlt1SingleMuonNoIPDecision_Dec;   //!
   TBranch        *b_Xib_Hlt1SingleMuonNoIPDecision_TIS;   //!
   TBranch        *b_Xib_Hlt1SingleMuonNoIPDecision_TOS;   //!
   TBranch        *b_Xib_Hlt1DiMuonHighMassDecision_Dec;   //!
   TBranch        *b_Xib_Hlt1DiMuonHighMassDecision_TIS;   //!
   TBranch        *b_Xib_Hlt1DiMuonHighMassDecision_TOS;   //!
   TBranch        *b_Xib_Hlt1TrackAllL0Decision_Dec;   //!
   TBranch        *b_Xib_Hlt1TrackAllL0Decision_TIS;   //!
   TBranch        *b_Xib_Hlt1TrackAllL0Decision_TOS;   //!
   TBranch        *b_Xib_Hlt1TrackMuonDecision_Dec;   //!
   TBranch        *b_Xib_Hlt1TrackMuonDecision_TIS;   //!
   TBranch        *b_Xib_Hlt1TrackMuonDecision_TOS;   //!
   TBranch        *b_Xib_Hlt2DiMuonDetachedHeavyDecision_Dec;   //!
   TBranch        *b_Xib_Hlt2DiMuonDetachedHeavyDecision_TIS;   //!
   TBranch        *b_Xib_Hlt2DiMuonDetachedHeavyDecision_TOS;   //!
   TBranch        *b_Xib_Hlt2DiMuonDetachedJPsiDecision_Dec;   //!
   TBranch        *b_Xib_Hlt2DiMuonDetachedJPsiDecision_TIS;   //!
   TBranch        *b_Xib_Hlt2DiMuonDetachedJPsiDecision_TOS;   //!
   TBranch        *b_Jpsi_CosTheta;   //!
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
   TBranch        *b_Jpsi_ORIVX_X;   //!
   TBranch        *b_Jpsi_ORIVX_Y;   //!
   TBranch        *b_Jpsi_ORIVX_Z;   //!
   TBranch        *b_Jpsi_ORIVX_XERR;   //!
   TBranch        *b_Jpsi_ORIVX_YERR;   //!
   TBranch        *b_Jpsi_ORIVX_ZERR;   //!
   TBranch        *b_Jpsi_ORIVX_CHI2;   //!
   TBranch        *b_Jpsi_ORIVX_NDOF;   //!
   TBranch        *b_Jpsi_ORIVX_COV_;   //!
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
   TBranch        *b_muplus_ORIVX_X;   //!
   TBranch        *b_muplus_ORIVX_Y;   //!
   TBranch        *b_muplus_ORIVX_Z;   //!
   TBranch        *b_muplus_ORIVX_XERR;   //!
   TBranch        *b_muplus_ORIVX_YERR;   //!
   TBranch        *b_muplus_ORIVX_ZERR;   //!
   TBranch        *b_muplus_ORIVX_CHI2;   //!
   TBranch        *b_muplus_ORIVX_NDOF;   //!
   TBranch        *b_muplus_ORIVX_COV_;   //!
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
   TBranch        *b_muminus_ORIVX_X;   //!
   TBranch        *b_muminus_ORIVX_Y;   //!
   TBranch        *b_muminus_ORIVX_Z;   //!
   TBranch        *b_muminus_ORIVX_XERR;   //!
   TBranch        *b_muminus_ORIVX_YERR;   //!
   TBranch        *b_muminus_ORIVX_ZERR;   //!
   TBranch        *b_muminus_ORIVX_CHI2;   //!
   TBranch        *b_muminus_ORIVX_NDOF;   //!
   TBranch        *b_muminus_ORIVX_COV_;   //!
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
   TBranch        *b_Xi_CosTheta;   //!
   TBranch        *b_Xi_ENDVERTEX_X;   //!
   TBranch        *b_Xi_ENDVERTEX_Y;   //!
   TBranch        *b_Xi_ENDVERTEX_Z;   //!
   TBranch        *b_Xi_ENDVERTEX_XERR;   //!
   TBranch        *b_Xi_ENDVERTEX_YERR;   //!
   TBranch        *b_Xi_ENDVERTEX_ZERR;   //!
   TBranch        *b_Xi_ENDVERTEX_CHI2;   //!
   TBranch        *b_Xi_ENDVERTEX_NDOF;   //!
   TBranch        *b_Xi_ENDVERTEX_COV_;   //!
   TBranch        *b_Xi_OWNPV_X;   //!
   TBranch        *b_Xi_OWNPV_Y;   //!
   TBranch        *b_Xi_OWNPV_Z;   //!
   TBranch        *b_Xi_OWNPV_XERR;   //!
   TBranch        *b_Xi_OWNPV_YERR;   //!
   TBranch        *b_Xi_OWNPV_ZERR;   //!
   TBranch        *b_Xi_OWNPV_CHI2;   //!
   TBranch        *b_Xi_OWNPV_NDOF;   //!
   TBranch        *b_Xi_OWNPV_COV_;   //!
   TBranch        *b_Xi_IP_OWNPV;   //!
   TBranch        *b_Xi_IPCHI2_OWNPV;   //!
   TBranch        *b_Xi_FD_OWNPV;   //!
   TBranch        *b_Xi_FDCHI2_OWNPV;   //!
   TBranch        *b_Xi_DIRA_OWNPV;   //!
   TBranch        *b_Xi_ORIVX_X;   //!
   TBranch        *b_Xi_ORIVX_Y;   //!
   TBranch        *b_Xi_ORIVX_Z;   //!
   TBranch        *b_Xi_ORIVX_XERR;   //!
   TBranch        *b_Xi_ORIVX_YERR;   //!
   TBranch        *b_Xi_ORIVX_ZERR;   //!
   TBranch        *b_Xi_ORIVX_CHI2;   //!
   TBranch        *b_Xi_ORIVX_NDOF;   //!
   TBranch        *b_Xi_ORIVX_COV_;   //!
   TBranch        *b_Xi_FD_ORIVX;   //!
   TBranch        *b_Xi_FDCHI2_ORIVX;   //!
   TBranch        *b_Xi_DIRA_ORIVX;   //!
   TBranch        *b_Xi_P;   //!
   TBranch        *b_Xi_PT;   //!
   TBranch        *b_Xi_PE;   //!
   TBranch        *b_Xi_PX;   //!
   TBranch        *b_Xi_PY;   //!
   TBranch        *b_Xi_PZ;   //!
   TBranch        *b_Xi_MM;   //!
   TBranch        *b_Xi_MMERR;   //!
   TBranch        *b_Xi_M;   //!
   TBranch        *b_Xi_ID;   //!
   TBranch        *b_Xi_TAU;   //!
   TBranch        *b_Xi_TAUERR;   //!
   TBranch        *b_Xi_TAUCHI2;   //!
   TBranch        *b_Xi_L0Global_Dec;   //!
   TBranch        *b_Xi_L0Global_TIS;   //!
   TBranch        *b_Xi_L0Global_TOS;   //!
   TBranch        *b_Xi_Hlt1Global_Dec;   //!
   TBranch        *b_Xi_Hlt1Global_TIS;   //!
   TBranch        *b_Xi_Hlt1Global_TOS;   //!
   TBranch        *b_Xi_Hlt1Phys_Dec;   //!
   TBranch        *b_Xi_Hlt1Phys_TIS;   //!
   TBranch        *b_Xi_Hlt1Phys_TOS;   //!
   TBranch        *b_Xi_Hlt2Global_Dec;   //!
   TBranch        *b_Xi_Hlt2Global_TIS;   //!
   TBranch        *b_Xi_Hlt2Global_TOS;   //!
   TBranch        *b_Xi_Hlt2Phys_Dec;   //!
   TBranch        *b_Xi_Hlt2Phys_TIS;   //!
   TBranch        *b_Xi_Hlt2Phys_TOS;   //!
   TBranch        *b_Xi_ETA;   //!
   TBranch        *b_Xi_PHI;   //!
   TBranch        *b_Xi_DTF_L_ARMEN_JpsiXiConstr;   //!
   TBranch        *b_Xi_DTF_L_ARMEN_JpsiXiLConstr;   //!
   TBranch        *b_Xi_DTF_L_LV01_JpsiXiConstr;   //!
   TBranch        *b_Xi_DTF_L_LV01_JpsiXiLConstr;   //!
   TBranch        *b_Xi_DTF_L_WMpipi_JpsiXiConstr;   //!
   TBranch        *b_Xi_DTF_L_WMpipi_JpsiXiLConstr;   //!
   TBranch        *b_L_CosTheta;   //!
   TBranch        *b_L_ENDVERTEX_X;   //!
   TBranch        *b_L_ENDVERTEX_Y;   //!
   TBranch        *b_L_ENDVERTEX_Z;   //!
   TBranch        *b_L_ENDVERTEX_XERR;   //!
   TBranch        *b_L_ENDVERTEX_YERR;   //!
   TBranch        *b_L_ENDVERTEX_ZERR;   //!
   TBranch        *b_L_ENDVERTEX_CHI2;   //!
   TBranch        *b_L_ENDVERTEX_NDOF;   //!
   TBranch        *b_L_ENDVERTEX_COV_;   //!
   TBranch        *b_L_OWNPV_X;   //!
   TBranch        *b_L_OWNPV_Y;   //!
   TBranch        *b_L_OWNPV_Z;   //!
   TBranch        *b_L_OWNPV_XERR;   //!
   TBranch        *b_L_OWNPV_YERR;   //!
   TBranch        *b_L_OWNPV_ZERR;   //!
   TBranch        *b_L_OWNPV_CHI2;   //!
   TBranch        *b_L_OWNPV_NDOF;   //!
   TBranch        *b_L_OWNPV_COV_;   //!
   TBranch        *b_L_IP_OWNPV;   //!
   TBranch        *b_L_IPCHI2_OWNPV;   //!
   TBranch        *b_L_FD_OWNPV;   //!
   TBranch        *b_L_FDCHI2_OWNPV;   //!
   TBranch        *b_L_DIRA_OWNPV;   //!
   TBranch        *b_L_ORIVX_X;   //!
   TBranch        *b_L_ORIVX_Y;   //!
   TBranch        *b_L_ORIVX_Z;   //!
   TBranch        *b_L_ORIVX_XERR;   //!
   TBranch        *b_L_ORIVX_YERR;   //!
   TBranch        *b_L_ORIVX_ZERR;   //!
   TBranch        *b_L_ORIVX_CHI2;   //!
   TBranch        *b_L_ORIVX_NDOF;   //!
   TBranch        *b_L_ORIVX_COV_;   //!
   TBranch        *b_L_FD_ORIVX;   //!
   TBranch        *b_L_FDCHI2_ORIVX;   //!
   TBranch        *b_L_DIRA_ORIVX;   //!
   TBranch        *b_L_P;   //!
   TBranch        *b_L_PT;   //!
   TBranch        *b_L_PE;   //!
   TBranch        *b_L_PX;   //!
   TBranch        *b_L_PY;   //!
   TBranch        *b_L_PZ;   //!
   TBranch        *b_L_MM;   //!
   TBranch        *b_L_MMERR;   //!
   TBranch        *b_L_M;   //!
   TBranch        *b_L_ID;   //!
   TBranch        *b_L_TAU;   //!
   TBranch        *b_L_TAUERR;   //!
   TBranch        *b_L_TAUCHI2;   //!
   TBranch        *b_L_L0Global_Dec;   //!
   TBranch        *b_L_L0Global_TIS;   //!
   TBranch        *b_L_L0Global_TOS;   //!
   TBranch        *b_L_Hlt1Global_Dec;   //!
   TBranch        *b_L_Hlt1Global_TIS;   //!
   TBranch        *b_L_Hlt1Global_TOS;   //!
   TBranch        *b_L_Hlt1Phys_Dec;   //!
   TBranch        *b_L_Hlt1Phys_TIS;   //!
   TBranch        *b_L_Hlt1Phys_TOS;   //!
   TBranch        *b_L_Hlt2Global_Dec;   //!
   TBranch        *b_L_Hlt2Global_TIS;   //!
   TBranch        *b_L_Hlt2Global_TOS;   //!
   TBranch        *b_L_Hlt2Phys_Dec;   //!
   TBranch        *b_L_Hlt2Phys_TIS;   //!
   TBranch        *b_L_Hlt2Phys_TOS;   //!
   TBranch        *b_L_ETA;   //!
   TBranch        *b_L_PHI;   //!
   TBranch        *b_L_ARMEN;   //!
   TBranch        *b_L_BPVLTIME;   //!
   TBranch        *b_L_WMkk;   //!
   TBranch        *b_L_WMkpi;   //!
   TBranch        *b_L_WMpik;   //!
   TBranch        *b_L_WMpipi;   //!
   TBranch        *b_L_WMpp;   //!
   TBranch        *b_L_dm;   //!
   TBranch        *b_L_lv01;   //!
   TBranch        *b_L_lv02;   //!
   TBranch        *b_p_CosTheta;   //!
   TBranch        *b_p_OWNPV_X;   //!
   TBranch        *b_p_OWNPV_Y;   //!
   TBranch        *b_p_OWNPV_Z;   //!
   TBranch        *b_p_OWNPV_XERR;   //!
   TBranch        *b_p_OWNPV_YERR;   //!
   TBranch        *b_p_OWNPV_ZERR;   //!
   TBranch        *b_p_OWNPV_CHI2;   //!
   TBranch        *b_p_OWNPV_NDOF;   //!
   TBranch        *b_p_OWNPV_COV_;   //!
   TBranch        *b_p_IP_OWNPV;   //!
   TBranch        *b_p_IPCHI2_OWNPV;   //!
   TBranch        *b_p_ORIVX_X;   //!
   TBranch        *b_p_ORIVX_Y;   //!
   TBranch        *b_p_ORIVX_Z;   //!
   TBranch        *b_p_ORIVX_XERR;   //!
   TBranch        *b_p_ORIVX_YERR;   //!
   TBranch        *b_p_ORIVX_ZERR;   //!
   TBranch        *b_p_ORIVX_CHI2;   //!
   TBranch        *b_p_ORIVX_NDOF;   //!
   TBranch        *b_p_ORIVX_COV_;   //!
   TBranch        *b_p_P;   //!
   TBranch        *b_p_PT;   //!
   TBranch        *b_p_PE;   //!
   TBranch        *b_p_PX;   //!
   TBranch        *b_p_PY;   //!
   TBranch        *b_p_PZ;   //!
   TBranch        *b_p_M;   //!
   TBranch        *b_p_ID;   //!
   TBranch        *b_p_PIDe;   //!
   TBranch        *b_p_PIDmu;   //!
   TBranch        *b_p_PIDK;   //!
   TBranch        *b_p_PIDp;   //!
   TBranch        *b_p_ProbNNe;   //!
   TBranch        *b_p_ProbNNk;   //!
   TBranch        *b_p_ProbNNp;   //!
   TBranch        *b_p_ProbNNpi;   //!
   TBranch        *b_p_ProbNNmu;   //!
   TBranch        *b_p_ProbNNghost;   //!
   TBranch        *b_p_hasMuon;   //!
   TBranch        *b_p_isMuon;   //!
   TBranch        *b_p_hasRich;   //!
   TBranch        *b_p_UsedRichAerogel;   //!
   TBranch        *b_p_UsedRich1Gas;   //!
   TBranch        *b_p_UsedRich2Gas;   //!
   TBranch        *b_p_RichAboveElThres;   //!
   TBranch        *b_p_RichAboveMuThres;   //!
   TBranch        *b_p_RichAbovePiThres;   //!
   TBranch        *b_p_RichAboveKaThres;   //!
   TBranch        *b_p_RichAbovePrThres;   //!
   TBranch        *b_p_hasCalo;   //!
   TBranch        *b_p_L0Global_Dec;   //!
   TBranch        *b_p_L0Global_TIS;   //!
   TBranch        *b_p_L0Global_TOS;   //!
   TBranch        *b_p_Hlt1Global_Dec;   //!
   TBranch        *b_p_Hlt1Global_TIS;   //!
   TBranch        *b_p_Hlt1Global_TOS;   //!
   TBranch        *b_p_Hlt1Phys_Dec;   //!
   TBranch        *b_p_Hlt1Phys_TIS;   //!
   TBranch        *b_p_Hlt1Phys_TOS;   //!
   TBranch        *b_p_Hlt2Global_Dec;   //!
   TBranch        *b_p_Hlt2Global_TIS;   //!
   TBranch        *b_p_Hlt2Global_TOS;   //!
   TBranch        *b_p_Hlt2Phys_Dec;   //!
   TBranch        *b_p_Hlt2Phys_TIS;   //!
   TBranch        *b_p_Hlt2Phys_TOS;   //!
   TBranch        *b_p_TRACK_Type;   //!
   TBranch        *b_p_TRACK_Key;   //!
   TBranch        *b_p_TRACK_CHI2NDOF;   //!
   TBranch        *b_p_TRACK_PCHI2;   //!
   TBranch        *b_p_TRACK_MatchCHI2;   //!
   TBranch        *b_p_TRACK_GhostProb;   //!
   TBranch        *b_p_TRACK_CloneDist;   //!
   TBranch        *b_p_TRACK_Likelihood;   //!
   TBranch        *b_p_ETA;   //!
   TBranch        *b_p_PHI;   //!
   TBranch        *b_pi_CosTheta;   //!
   TBranch        *b_pi_OWNPV_X;   //!
   TBranch        *b_pi_OWNPV_Y;   //!
   TBranch        *b_pi_OWNPV_Z;   //!
   TBranch        *b_pi_OWNPV_XERR;   //!
   TBranch        *b_pi_OWNPV_YERR;   //!
   TBranch        *b_pi_OWNPV_ZERR;   //!
   TBranch        *b_pi_OWNPV_CHI2;   //!
   TBranch        *b_pi_OWNPV_NDOF;   //!
   TBranch        *b_pi_OWNPV_COV_;   //!
   TBranch        *b_pi_IP_OWNPV;   //!
   TBranch        *b_pi_IPCHI2_OWNPV;   //!
   TBranch        *b_pi_ORIVX_X;   //!
   TBranch        *b_pi_ORIVX_Y;   //!
   TBranch        *b_pi_ORIVX_Z;   //!
   TBranch        *b_pi_ORIVX_XERR;   //!
   TBranch        *b_pi_ORIVX_YERR;   //!
   TBranch        *b_pi_ORIVX_ZERR;   //!
   TBranch        *b_pi_ORIVX_CHI2;   //!
   TBranch        *b_pi_ORIVX_NDOF;   //!
   TBranch        *b_pi_ORIVX_COV_;   //!
   TBranch        *b_pi_P;   //!
   TBranch        *b_pi_PT;   //!
   TBranch        *b_pi_PE;   //!
   TBranch        *b_pi_PX;   //!
   TBranch        *b_pi_PY;   //!
   TBranch        *b_pi_PZ;   //!
   TBranch        *b_pi_M;   //!
   TBranch        *b_pi_ID;   //!
   TBranch        *b_pi_PIDe;   //!
   TBranch        *b_pi_PIDmu;   //!
   TBranch        *b_pi_PIDK;   //!
   TBranch        *b_pi_PIDp;   //!
   TBranch        *b_pi_ProbNNe;   //!
   TBranch        *b_pi_ProbNNk;   //!
   TBranch        *b_pi_ProbNNp;   //!
   TBranch        *b_pi_ProbNNpi;   //!
   TBranch        *b_pi_ProbNNmu;   //!
   TBranch        *b_pi_ProbNNghost;   //!
   TBranch        *b_pi_hasMuon;   //!
   TBranch        *b_pi_isMuon;   //!
   TBranch        *b_pi_hasRich;   //!
   TBranch        *b_pi_UsedRichAerogel;   //!
   TBranch        *b_pi_UsedRich1Gas;   //!
   TBranch        *b_pi_UsedRich2Gas;   //!
   TBranch        *b_pi_RichAboveElThres;   //!
   TBranch        *b_pi_RichAboveMuThres;   //!
   TBranch        *b_pi_RichAbovePiThres;   //!
   TBranch        *b_pi_RichAboveKaThres;   //!
   TBranch        *b_pi_RichAbovePrThres;   //!
   TBranch        *b_pi_hasCalo;   //!
   TBranch        *b_pi_L0Global_Dec;   //!
   TBranch        *b_pi_L0Global_TIS;   //!
   TBranch        *b_pi_L0Global_TOS;   //!
   TBranch        *b_pi_Hlt1Global_Dec;   //!
   TBranch        *b_pi_Hlt1Global_TIS;   //!
   TBranch        *b_pi_Hlt1Global_TOS;   //!
   TBranch        *b_pi_Hlt1Phys_Dec;   //!
   TBranch        *b_pi_Hlt1Phys_TIS;   //!
   TBranch        *b_pi_Hlt1Phys_TOS;   //!
   TBranch        *b_pi_Hlt2Global_Dec;   //!
   TBranch        *b_pi_Hlt2Global_TIS;   //!
   TBranch        *b_pi_Hlt2Global_TOS;   //!
   TBranch        *b_pi_Hlt2Phys_Dec;   //!
   TBranch        *b_pi_Hlt2Phys_TIS;   //!
   TBranch        *b_pi_Hlt2Phys_TOS;   //!
   TBranch        *b_pi_TRACK_Type;   //!
   TBranch        *b_pi_TRACK_Key;   //!
   TBranch        *b_pi_TRACK_CHI2NDOF;   //!
   TBranch        *b_pi_TRACK_PCHI2;   //!
   TBranch        *b_pi_TRACK_MatchCHI2;   //!
   TBranch        *b_pi_TRACK_GhostProb;   //!
   TBranch        *b_pi_TRACK_CloneDist;   //!
   TBranch        *b_pi_TRACK_Likelihood;   //!
   TBranch        *b_pi_ETA;   //!
   TBranch        *b_pi_PHI;   //!
   TBranch        *b_BachPi_CosTheta;   //!
   TBranch        *b_BachPi_OWNPV_X;   //!
   TBranch        *b_BachPi_OWNPV_Y;   //!
   TBranch        *b_BachPi_OWNPV_Z;   //!
   TBranch        *b_BachPi_OWNPV_XERR;   //!
   TBranch        *b_BachPi_OWNPV_YERR;   //!
   TBranch        *b_BachPi_OWNPV_ZERR;   //!
   TBranch        *b_BachPi_OWNPV_CHI2;   //!
   TBranch        *b_BachPi_OWNPV_NDOF;   //!
   TBranch        *b_BachPi_OWNPV_COV_;   //!
   TBranch        *b_BachPi_IP_OWNPV;   //!
   TBranch        *b_BachPi_IPCHI2_OWNPV;   //!
   TBranch        *b_BachPi_ORIVX_X;   //!
   TBranch        *b_BachPi_ORIVX_Y;   //!
   TBranch        *b_BachPi_ORIVX_Z;   //!
   TBranch        *b_BachPi_ORIVX_XERR;   //!
   TBranch        *b_BachPi_ORIVX_YERR;   //!
   TBranch        *b_BachPi_ORIVX_ZERR;   //!
   TBranch        *b_BachPi_ORIVX_CHI2;   //!
   TBranch        *b_BachPi_ORIVX_NDOF;   //!
   TBranch        *b_BachPi_ORIVX_COV_;   //!
   TBranch        *b_BachPi_P;   //!
   TBranch        *b_BachPi_PT;   //!
   TBranch        *b_BachPi_PE;   //!
   TBranch        *b_BachPi_PX;   //!
   TBranch        *b_BachPi_PY;   //!
   TBranch        *b_BachPi_PZ;   //!
   TBranch        *b_BachPi_M;   //!
   TBranch        *b_BachPi_ID;   //!
   TBranch        *b_BachPi_PIDe;   //!
   TBranch        *b_BachPi_PIDmu;   //!
   TBranch        *b_BachPi_PIDK;   //!
   TBranch        *b_BachPi_PIDp;   //!
   TBranch        *b_BachPi_ProbNNe;   //!
   TBranch        *b_BachPi_ProbNNk;   //!
   TBranch        *b_BachPi_ProbNNp;   //!
   TBranch        *b_BachPi_ProbNNpi;   //!
   TBranch        *b_BachPi_ProbNNmu;   //!
   TBranch        *b_BachPi_ProbNNghost;   //!
   TBranch        *b_BachPi_hasMuon;   //!
   TBranch        *b_BachPi_isMuon;   //!
   TBranch        *b_BachPi_hasRich;   //!
   TBranch        *b_BachPi_UsedRichAerogel;   //!
   TBranch        *b_BachPi_UsedRich1Gas;   //!
   TBranch        *b_BachPi_UsedRich2Gas;   //!
   TBranch        *b_BachPi_RichAboveElThres;   //!
   TBranch        *b_BachPi_RichAboveMuThres;   //!
   TBranch        *b_BachPi_RichAbovePiThres;   //!
   TBranch        *b_BachPi_RichAboveKaThres;   //!
   TBranch        *b_BachPi_RichAbovePrThres;   //!
   TBranch        *b_BachPi_hasCalo;   //!
   TBranch        *b_BachPi_L0Global_Dec;   //!
   TBranch        *b_BachPi_L0Global_TIS;   //!
   TBranch        *b_BachPi_L0Global_TOS;   //!
   TBranch        *b_BachPi_Hlt1Global_Dec;   //!
   TBranch        *b_BachPi_Hlt1Global_TIS;   //!
   TBranch        *b_BachPi_Hlt1Global_TOS;   //!
   TBranch        *b_BachPi_Hlt1Phys_Dec;   //!
   TBranch        *b_BachPi_Hlt1Phys_TIS;   //!
   TBranch        *b_BachPi_Hlt1Phys_TOS;   //!
   TBranch        *b_BachPi_Hlt2Global_Dec;   //!
   TBranch        *b_BachPi_Hlt2Global_TIS;   //!
   TBranch        *b_BachPi_Hlt2Global_TOS;   //!
   TBranch        *b_BachPi_Hlt2Phys_Dec;   //!
   TBranch        *b_BachPi_Hlt2Phys_TIS;   //!
   TBranch        *b_BachPi_Hlt2Phys_TOS;   //!
   TBranch        *b_BachPi_TRACK_Type;   //!
   TBranch        *b_BachPi_TRACK_Key;   //!
   TBranch        *b_BachPi_TRACK_CHI2NDOF;   //!
   TBranch        *b_BachPi_TRACK_PCHI2;   //!
   TBranch        *b_BachPi_TRACK_MatchCHI2;   //!
   TBranch        *b_BachPi_TRACK_GhostProb;   //!
   TBranch        *b_BachPi_TRACK_CloneDist;   //!
   TBranch        *b_BachPi_TRACK_Likelihood;   //!
   TBranch        *b_BachPi_ETA;   //!
   TBranch        *b_BachPi_PHI;   //!
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

   Xib(TTree *tree=0);
   virtual ~Xib();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Xib_cxx
Xib::Xib(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Xib2JpsiXi_Run2_DD_L0JPsiTOS_NewStrip.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Xib2JpsiXi_Run2_DD_L0JPsiTOS_NewStrip.root");
      }
      f->GetObject("MyTuple",tree);

   }
   Init(tree);
}

Xib::~Xib()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Xib::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Xib::LoadTree(Long64_t entry)
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

void Xib::Init(TTree *tree)
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

   fChain->SetBranchAddress("Xib_ENDVERTEX_X", &Xib_ENDVERTEX_X, &b_Xib_ENDVERTEX_X);
   fChain->SetBranchAddress("Xib_ENDVERTEX_Y", &Xib_ENDVERTEX_Y, &b_Xib_ENDVERTEX_Y);
   fChain->SetBranchAddress("Xib_ENDVERTEX_Z", &Xib_ENDVERTEX_Z, &b_Xib_ENDVERTEX_Z);
   fChain->SetBranchAddress("Xib_ENDVERTEX_XERR", &Xib_ENDVERTEX_XERR, &b_Xib_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Xib_ENDVERTEX_YERR", &Xib_ENDVERTEX_YERR, &b_Xib_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Xib_ENDVERTEX_ZERR", &Xib_ENDVERTEX_ZERR, &b_Xib_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Xib_ENDVERTEX_CHI2", &Xib_ENDVERTEX_CHI2, &b_Xib_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Xib_ENDVERTEX_NDOF", &Xib_ENDVERTEX_NDOF, &b_Xib_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Xib_ENDVERTEX_COV_", Xib_ENDVERTEX_COV_, &b_Xib_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Xib_OWNPV_X", &Xib_OWNPV_X, &b_Xib_OWNPV_X);
   fChain->SetBranchAddress("Xib_OWNPV_Y", &Xib_OWNPV_Y, &b_Xib_OWNPV_Y);
   fChain->SetBranchAddress("Xib_OWNPV_Z", &Xib_OWNPV_Z, &b_Xib_OWNPV_Z);
   fChain->SetBranchAddress("Xib_OWNPV_XERR", &Xib_OWNPV_XERR, &b_Xib_OWNPV_XERR);
   fChain->SetBranchAddress("Xib_OWNPV_YERR", &Xib_OWNPV_YERR, &b_Xib_OWNPV_YERR);
   fChain->SetBranchAddress("Xib_OWNPV_ZERR", &Xib_OWNPV_ZERR, &b_Xib_OWNPV_ZERR);
   fChain->SetBranchAddress("Xib_OWNPV_CHI2", &Xib_OWNPV_CHI2, &b_Xib_OWNPV_CHI2);
   fChain->SetBranchAddress("Xib_OWNPV_NDOF", &Xib_OWNPV_NDOF, &b_Xib_OWNPV_NDOF);
   fChain->SetBranchAddress("Xib_OWNPV_COV_", Xib_OWNPV_COV_, &b_Xib_OWNPV_COV_);
   fChain->SetBranchAddress("Xib_IP_OWNPV", &Xib_IP_OWNPV, &b_Xib_IP_OWNPV);
   fChain->SetBranchAddress("Xib_IPCHI2_OWNPV", &Xib_IPCHI2_OWNPV, &b_Xib_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Xib_FD_OWNPV", &Xib_FD_OWNPV, &b_Xib_FD_OWNPV);
   fChain->SetBranchAddress("Xib_FDCHI2_OWNPV", &Xib_FDCHI2_OWNPV, &b_Xib_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Xib_DIRA_OWNPV", &Xib_DIRA_OWNPV, &b_Xib_DIRA_OWNPV);
   fChain->SetBranchAddress("Xib_P", &Xib_P, &b_Xib_P);
   fChain->SetBranchAddress("Xib_PT", &Xib_PT, &b_Xib_PT);
   fChain->SetBranchAddress("Xib_PE", &Xib_PE, &b_Xib_PE);
   fChain->SetBranchAddress("Xib_PX", &Xib_PX, &b_Xib_PX);
   fChain->SetBranchAddress("Xib_PY", &Xib_PY, &b_Xib_PY);
   fChain->SetBranchAddress("Xib_PZ", &Xib_PZ, &b_Xib_PZ);
   fChain->SetBranchAddress("Xib_MM", &Xib_MM, &b_Xib_MM);
   fChain->SetBranchAddress("Xib_MMERR", &Xib_MMERR, &b_Xib_MMERR);
   fChain->SetBranchAddress("Xib_M", &Xib_M, &b_Xib_M);
   fChain->SetBranchAddress("Xib_ID", &Xib_ID, &b_Xib_ID);
   fChain->SetBranchAddress("Xib_TAU", &Xib_TAU, &b_Xib_TAU);
   fChain->SetBranchAddress("Xib_TAUERR", &Xib_TAUERR, &b_Xib_TAUERR);
   fChain->SetBranchAddress("Xib_TAUCHI2", &Xib_TAUCHI2, &b_Xib_TAUCHI2);
   fChain->SetBranchAddress("Xib_L0Global_Dec", &Xib_L0Global_Dec, &b_Xib_L0Global_Dec);
   fChain->SetBranchAddress("Xib_L0Global_TIS", &Xib_L0Global_TIS, &b_Xib_L0Global_TIS);
   fChain->SetBranchAddress("Xib_L0Global_TOS", &Xib_L0Global_TOS, &b_Xib_L0Global_TOS);
   fChain->SetBranchAddress("Xib_Hlt1Global_Dec", &Xib_Hlt1Global_Dec, &b_Xib_Hlt1Global_Dec);
   fChain->SetBranchAddress("Xib_Hlt1Global_TIS", &Xib_Hlt1Global_TIS, &b_Xib_Hlt1Global_TIS);
   fChain->SetBranchAddress("Xib_Hlt1Global_TOS", &Xib_Hlt1Global_TOS, &b_Xib_Hlt1Global_TOS);
   fChain->SetBranchAddress("Xib_Hlt1Phys_Dec", &Xib_Hlt1Phys_Dec, &b_Xib_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Xib_Hlt1Phys_TIS", &Xib_Hlt1Phys_TIS, &b_Xib_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Xib_Hlt1Phys_TOS", &Xib_Hlt1Phys_TOS, &b_Xib_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Xib_Hlt2Global_Dec", &Xib_Hlt2Global_Dec, &b_Xib_Hlt2Global_Dec);
   fChain->SetBranchAddress("Xib_Hlt2Global_TIS", &Xib_Hlt2Global_TIS, &b_Xib_Hlt2Global_TIS);
   fChain->SetBranchAddress("Xib_Hlt2Global_TOS", &Xib_Hlt2Global_TOS, &b_Xib_Hlt2Global_TOS);
   fChain->SetBranchAddress("Xib_Hlt2Phys_Dec", &Xib_Hlt2Phys_Dec, &b_Xib_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Xib_Hlt2Phys_TIS", &Xib_Hlt2Phys_TIS, &b_Xib_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Xib_Hlt2Phys_TOS", &Xib_Hlt2Phys_TOS, &b_Xib_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Xib_ETA", &Xib_ETA, &b_Xib_ETA);
   fChain->SetBranchAddress("Xib_PHI", &Xib_PHI, &b_Xib_PHI);
   fChain->SetBranchAddress("Xib_DTF_CTAU", &Xib_DTF_CTAU, &b_Xib_DTF_CTAU);
   fChain->SetBranchAddress("Xib_DTF_CTAUERR", &Xib_DTF_CTAUERR, &b_Xib_DTF_CTAUERR);
   fChain->SetBranchAddress("Xib_DTF_CTAUS", &Xib_DTF_CTAUS, &b_Xib_DTF_CTAUS);
   fChain->SetBranchAddress("Xib_DTF_M_JpsiXiLConstr", &Xib_DTF_M_JpsiXiLConstr, &b_Xib_DTF_M_JpsiXiLConstr);
   fChain->SetBranchAddress("Xib_DTF_VCHI2NDOF", &Xib_DTF_VCHI2NDOF, &b_Xib_DTF_VCHI2NDOF);
   fChain->SetBranchAddress("Added_n_Particles", &Added_n_Particles, &b_Added_n_Particles);
   fChain->SetBranchAddress("psi_1S_H_M", psi_1S_H_M, &b_psi_1S_H_M);
   fChain->SetBranchAddress("psi_1S_H_M_pion", psi_1S_H_M_pion, &b_psi_1S_H_M_pion);
   fChain->SetBranchAddress("psi_1S_H_M_kaon", psi_1S_H_M_kaon, &b_psi_1S_H_M_kaon);
   fChain->SetBranchAddress("psi_1S_H_M_proton", psi_1S_H_M_proton, &b_psi_1S_H_M_proton);
   fChain->SetBranchAddress("psi_1S_H_PT", psi_1S_H_PT, &b_psi_1S_H_PT);
   fChain->SetBranchAddress("psi_1S_H_ETA", psi_1S_H_ETA, &b_psi_1S_H_ETA);
   fChain->SetBranchAddress("psi_1S_H_PHI", psi_1S_H_PHI, &b_psi_1S_H_PHI);
   fChain->SetBranchAddress("psi_1S_H_PX", psi_1S_H_PX, &b_psi_1S_H_PX);
   fChain->SetBranchAddress("psi_1S_H_PY", psi_1S_H_PY, &b_psi_1S_H_PY);
   fChain->SetBranchAddress("psi_1S_H_PZ", psi_1S_H_PZ, &b_psi_1S_H_PZ);
   fChain->SetBranchAddress("psi_1S_H_FDCHI2_OLD", psi_1S_H_FDCHI2_OLD, &b_psi_1S_H_FDCHI2_OLD);
   fChain->SetBranchAddress("psi_1S_H_FDCHI2_NEW", psi_1S_H_FDCHI2_NEW, &b_psi_1S_H_FDCHI2_NEW);
   fChain->SetBranchAddress("psi_1S_H_FD_OLD", psi_1S_H_FD_OLD, &b_psi_1S_H_FD_OLD);
   fChain->SetBranchAddress("psi_1S_H_FD_NEW", psi_1S_H_FD_NEW, &b_psi_1S_H_FD_NEW);
   fChain->SetBranchAddress("psi_1S_H_VERTEXCHI2_OLD", psi_1S_H_VERTEXCHI2_OLD, &b_psi_1S_H_VERTEXCHI2_OLD);
   fChain->SetBranchAddress("psi_1S_H_VERTEXCHI2_NEW", psi_1S_H_VERTEXCHI2_NEW, &b_psi_1S_H_VERTEXCHI2_NEW);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_X_NEW", psi_1S_H_VERTEX_X_NEW, &b_psi_1S_H_VERTEX_X_NEW);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_Y_NEW", psi_1S_H_VERTEX_Y_NEW, &b_psi_1S_H_VERTEX_Y_NEW);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_Z_NEW", psi_1S_H_VERTEX_Z_NEW, &b_psi_1S_H_VERTEX_Z_NEW);
   fChain->SetBranchAddress("psi_1S_H_PV_X_NEW", psi_1S_H_PV_X_NEW, &b_psi_1S_H_PV_X_NEW);
   fChain->SetBranchAddress("psi_1S_H_PV_Y_NEW", psi_1S_H_PV_Y_NEW, &b_psi_1S_H_PV_Y_NEW);
   fChain->SetBranchAddress("psi_1S_H_PV_Z_NEW", psi_1S_H_PV_Z_NEW, &b_psi_1S_H_PV_Z_NEW);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_COV_XX", psi_1S_H_VERTEX_COV_XX, &b_psi_1S_H_VERTEX_COV_XX);
   fChain->SetBranchAddress("psi_1S_H_PV_COV_XX", psi_1S_H_PV_COV_XX, &b_psi_1S_H_PV_COV_XX);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_COV_XY", psi_1S_H_VERTEX_COV_XY, &b_psi_1S_H_VERTEX_COV_XY);
   fChain->SetBranchAddress("psi_1S_H_PV_COV_XY", psi_1S_H_PV_COV_XY, &b_psi_1S_H_PV_COV_XY);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_COV_XZ", psi_1S_H_VERTEX_COV_XZ, &b_psi_1S_H_VERTEX_COV_XZ);
   fChain->SetBranchAddress("psi_1S_H_PV_COV_XZ", psi_1S_H_PV_COV_XZ, &b_psi_1S_H_PV_COV_XZ);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_COV_YY", psi_1S_H_VERTEX_COV_YY, &b_psi_1S_H_VERTEX_COV_YY);
   fChain->SetBranchAddress("psi_1S_H_PV_COV_YY", psi_1S_H_PV_COV_YY, &b_psi_1S_H_PV_COV_YY);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_COV_YZ", psi_1S_H_VERTEX_COV_YZ, &b_psi_1S_H_VERTEX_COV_YZ);
   fChain->SetBranchAddress("psi_1S_H_PV_COV_YZ", psi_1S_H_PV_COV_YZ, &b_psi_1S_H_PV_COV_YZ);
   fChain->SetBranchAddress("psi_1S_H_VERTEX_COV_ZZ", psi_1S_H_VERTEX_COV_ZZ, &b_psi_1S_H_VERTEX_COV_ZZ);
   fChain->SetBranchAddress("psi_1S_H_PV_COV_ZZ", psi_1S_H_PV_COV_ZZ, &b_psi_1S_H_PV_COV_ZZ);
   fChain->SetBranchAddress("psi_1S_H_VertexRefitStatus", psi_1S_H_VertexRefitStatus, &b_psi_1S_H_VertexRefitStatus);
   fChain->SetBranchAddress("psi_1S_H_IPCHI2_OLD", psi_1S_H_IPCHI2_OLD, &b_psi_1S_H_IPCHI2_OLD);
   fChain->SetBranchAddress("psi_1S_H_IP_OLD", psi_1S_H_IP_OLD, &b_psi_1S_H_IP_OLD);
   fChain->SetBranchAddress("psi_1S_H_IPCHI2_NEW", psi_1S_H_IPCHI2_NEW, &b_psi_1S_H_IPCHI2_NEW);
   fChain->SetBranchAddress("psi_1S_H_IP_NEW", psi_1S_H_IP_NEW, &b_psi_1S_H_IP_NEW);
   fChain->SetBranchAddress("psi_1S_H_MINIPCHI2", psi_1S_H_MINIPCHI2, &b_psi_1S_H_MINIPCHI2);
   fChain->SetBranchAddress("Xib_H_OPENING", Xib_H_OPENING, &b_Xib_H_OPENING);
   fChain->SetBranchAddress("Added_H_PX", Added_H_PX, &b_Added_H_PX);
   fChain->SetBranchAddress("Added_H_PY", Added_H_PY, &b_Added_H_PY);
   fChain->SetBranchAddress("Added_H_PZ", Added_H_PZ, &b_Added_H_PZ);
   fChain->SetBranchAddress("Added_H_PT", Added_H_PT, &b_Added_H_PT);
   fChain->SetBranchAddress("Added_H_ETA", Added_H_ETA, &b_Added_H_ETA);
   fChain->SetBranchAddress("Added_H_PHI", Added_H_PHI, &b_Added_H_PHI);
   fChain->SetBranchAddress("Added_H_THETA", Added_H_THETA, &b_Added_H_THETA);
   fChain->SetBranchAddress("Added_H_DELTAPHI", Added_H_DELTAPHI, &b_Added_H_DELTAPHI);
   fChain->SetBranchAddress("Added_H_PID", Added_H_PID, &b_Added_H_PID);
   fChain->SetBranchAddress("Added_n_Tracks", &Added_n_Tracks, &b_Added_n_Tracks);
   fChain->SetBranchAddress("Added_H_PROBNNPID", Added_H_PROBNNPID, &b_Added_H_PROBNNPID);
   fChain->SetBranchAddress("Added_H_GHOST", Added_H_GHOST, &b_Added_H_GHOST);
   fChain->SetBranchAddress("Added_H_TRACKCHI2", Added_H_TRACKCHI2, &b_Added_H_TRACKCHI2);
   fChain->SetBranchAddress("Added_H_TRACKTYPE", Added_H_TRACKTYPE, &b_Added_H_TRACKTYPE);
   fChain->SetBranchAddress("Added_H_PIDe", Added_H_PIDe, &b_Added_H_PIDe);
   fChain->SetBranchAddress("Added_H_PIDK", Added_H_PIDK, &b_Added_H_PIDK);
   fChain->SetBranchAddress("Added_H_PIDmu", Added_H_PIDmu, &b_Added_H_PIDmu);
   fChain->SetBranchAddress("Added_H_PIDp", Added_H_PIDp, &b_Added_H_PIDp);
   fChain->SetBranchAddress("Added_H_ProbNNpi", Added_H_ProbNNpi, &b_Added_H_ProbNNpi);
   fChain->SetBranchAddress("Added_H_ProbNNk", Added_H_ProbNNk, &b_Added_H_ProbNNk);
   fChain->SetBranchAddress("Added_H_ProbNNp", Added_H_ProbNNp, &b_Added_H_ProbNNp);
   fChain->SetBranchAddress("Added_H_ProbNNe", Added_H_ProbNNe, &b_Added_H_ProbNNe);
   fChain->SetBranchAddress("Added_H_ProbNNmu", Added_H_ProbNNmu, &b_Added_H_ProbNNmu);
   fChain->SetBranchAddress("Xib_ConsXib_nPV", &Xib_ConsXib_nPV, &b_Xib_ConsXib_nPV);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_M", Xib_ConsXib_J_psi_1S_M, &b_Xib_ConsXib_J_psi_1S_M);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_MERR", Xib_ConsXib_J_psi_1S_MERR, &b_Xib_ConsXib_J_psi_1S_MERR);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_P", Xib_ConsXib_J_psi_1S_P, &b_Xib_ConsXib_J_psi_1S_P);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_PERR", Xib_ConsXib_J_psi_1S_PERR, &b_Xib_ConsXib_J_psi_1S_PERR);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_ctau", Xib_ConsXib_J_psi_1S_ctau, &b_Xib_ConsXib_J_psi_1S_ctau);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_ctauErr", Xib_ConsXib_J_psi_1S_ctauErr, &b_Xib_ConsXib_J_psi_1S_ctauErr);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_decayLength", Xib_ConsXib_J_psi_1S_decayLength, &b_Xib_ConsXib_J_psi_1S_decayLength);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_decayLengthErr", Xib_ConsXib_J_psi_1S_decayLengthErr, &b_Xib_ConsXib_J_psi_1S_decayLengthErr);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_0_ID", Xib_ConsXib_J_psi_1S_muminus_0_ID, &b_Xib_ConsXib_J_psi_1S_muminus_0_ID);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_0_PE", Xib_ConsXib_J_psi_1S_muminus_0_PE, &b_Xib_ConsXib_J_psi_1S_muminus_0_PE);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_0_PX", Xib_ConsXib_J_psi_1S_muminus_0_PX, &b_Xib_ConsXib_J_psi_1S_muminus_0_PX);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_0_PY", Xib_ConsXib_J_psi_1S_muminus_0_PY, &b_Xib_ConsXib_J_psi_1S_muminus_0_PY);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_0_PZ", Xib_ConsXib_J_psi_1S_muminus_0_PZ, &b_Xib_ConsXib_J_psi_1S_muminus_0_PZ);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_ID", Xib_ConsXib_J_psi_1S_muminus_ID, &b_Xib_ConsXib_J_psi_1S_muminus_ID);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_PE", Xib_ConsXib_J_psi_1S_muminus_PE, &b_Xib_ConsXib_J_psi_1S_muminus_PE);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_PX", Xib_ConsXib_J_psi_1S_muminus_PX, &b_Xib_ConsXib_J_psi_1S_muminus_PX);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_PY", Xib_ConsXib_J_psi_1S_muminus_PY, &b_Xib_ConsXib_J_psi_1S_muminus_PY);
   fChain->SetBranchAddress("Xib_ConsXib_J_psi_1S_muminus_PZ", Xib_ConsXib_J_psi_1S_muminus_PZ, &b_Xib_ConsXib_J_psi_1S_muminus_PZ);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_M", Xib_ConsXib_Lambda0_M, &b_Xib_ConsXib_Lambda0_M);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_MERR", Xib_ConsXib_Lambda0_MERR, &b_Xib_ConsXib_Lambda0_MERR);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_P", Xib_ConsXib_Lambda0_P, &b_Xib_ConsXib_Lambda0_P);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_PERR", Xib_ConsXib_Lambda0_PERR, &b_Xib_ConsXib_Lambda0_PERR);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_ctau", Xib_ConsXib_Lambda0_ctau, &b_Xib_ConsXib_Lambda0_ctau);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_ctauErr", Xib_ConsXib_Lambda0_ctauErr, &b_Xib_ConsXib_Lambda0_ctauErr);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_decayLength", Xib_ConsXib_Lambda0_decayLength, &b_Xib_ConsXib_Lambda0_decayLength);
   fChain->SetBranchAddress("Xib_ConsXib_Lambda0_decayLengthErr", Xib_ConsXib_Lambda0_decayLengthErr, &b_Xib_ConsXib_Lambda0_decayLengthErr);
   fChain->SetBranchAddress("Xib_ConsXib_M", Xib_ConsXib_M, &b_Xib_ConsXib_M);
   fChain->SetBranchAddress("Xib_ConsXib_MERR", Xib_ConsXib_MERR, &b_Xib_ConsXib_MERR);
   fChain->SetBranchAddress("Xib_ConsXib_P", Xib_ConsXib_P, &b_Xib_ConsXib_P);
   fChain->SetBranchAddress("Xib_ConsXib_PERR", Xib_ConsXib_PERR, &b_Xib_ConsXib_PERR);
   fChain->SetBranchAddress("Xib_ConsXib_PV_X", Xib_ConsXib_PV_X, &b_Xib_ConsXib_PV_X);
   fChain->SetBranchAddress("Xib_ConsXib_PV_Y", Xib_ConsXib_PV_Y, &b_Xib_ConsXib_PV_Y);
   fChain->SetBranchAddress("Xib_ConsXib_PV_Z", Xib_ConsXib_PV_Z, &b_Xib_ConsXib_PV_Z);
   fChain->SetBranchAddress("Xib_ConsXib_PV_key", Xib_ConsXib_PV_key, &b_Xib_ConsXib_PV_key);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_piplus_ID", Xib_ConsXib_Ximinus_Lambda0_piplus_ID, &b_Xib_ConsXib_Ximinus_Lambda0_piplus_ID);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_piplus_PE", Xib_ConsXib_Ximinus_Lambda0_piplus_PE, &b_Xib_ConsXib_Ximinus_Lambda0_piplus_PE);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_piplus_PX", Xib_ConsXib_Ximinus_Lambda0_piplus_PX, &b_Xib_ConsXib_Ximinus_Lambda0_piplus_PX);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_piplus_PY", Xib_ConsXib_Ximinus_Lambda0_piplus_PY, &b_Xib_ConsXib_Ximinus_Lambda0_piplus_PY);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_piplus_PZ", Xib_ConsXib_Ximinus_Lambda0_piplus_PZ, &b_Xib_ConsXib_Ximinus_Lambda0_piplus_PZ);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_pplus_ID", Xib_ConsXib_Ximinus_Lambda0_pplus_ID, &b_Xib_ConsXib_Ximinus_Lambda0_pplus_ID);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_pplus_PE", Xib_ConsXib_Ximinus_Lambda0_pplus_PE, &b_Xib_ConsXib_Ximinus_Lambda0_pplus_PE);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_pplus_PX", Xib_ConsXib_Ximinus_Lambda0_pplus_PX, &b_Xib_ConsXib_Ximinus_Lambda0_pplus_PX);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_pplus_PY", Xib_ConsXib_Ximinus_Lambda0_pplus_PY, &b_Xib_ConsXib_Ximinus_Lambda0_pplus_PY);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_Lambda0_pplus_PZ", Xib_ConsXib_Ximinus_Lambda0_pplus_PZ, &b_Xib_ConsXib_Ximinus_Lambda0_pplus_PZ);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_M", Xib_ConsXib_Ximinus_M, &b_Xib_ConsXib_Ximinus_M);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_MERR", Xib_ConsXib_Ximinus_MERR, &b_Xib_ConsXib_Ximinus_MERR);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_P", Xib_ConsXib_Ximinus_P, &b_Xib_ConsXib_Ximinus_P);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_PERR", Xib_ConsXib_Ximinus_PERR, &b_Xib_ConsXib_Ximinus_PERR);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_ctau", Xib_ConsXib_Ximinus_ctau, &b_Xib_ConsXib_Ximinus_ctau);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_ctauErr", Xib_ConsXib_Ximinus_ctauErr, &b_Xib_ConsXib_Ximinus_ctauErr);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_decayLength", Xib_ConsXib_Ximinus_decayLength, &b_Xib_ConsXib_Ximinus_decayLength);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_decayLengthErr", Xib_ConsXib_Ximinus_decayLengthErr, &b_Xib_ConsXib_Ximinus_decayLengthErr);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_piplus_ID", Xib_ConsXib_Ximinus_piplus_ID, &b_Xib_ConsXib_Ximinus_piplus_ID);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_piplus_PE", Xib_ConsXib_Ximinus_piplus_PE, &b_Xib_ConsXib_Ximinus_piplus_PE);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_piplus_PX", Xib_ConsXib_Ximinus_piplus_PX, &b_Xib_ConsXib_Ximinus_piplus_PX);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_piplus_PY", Xib_ConsXib_Ximinus_piplus_PY, &b_Xib_ConsXib_Ximinus_piplus_PY);
   fChain->SetBranchAddress("Xib_ConsXib_Ximinus_piplus_PZ", Xib_ConsXib_Ximinus_piplus_PZ, &b_Xib_ConsXib_Ximinus_piplus_PZ);
   fChain->SetBranchAddress("Xib_ConsXib_chi2", Xib_ConsXib_chi2, &b_Xib_ConsXib_chi2);
   fChain->SetBranchAddress("Xib_ConsXib_ctau", Xib_ConsXib_ctau, &b_Xib_ConsXib_ctau);
   fChain->SetBranchAddress("Xib_ConsXib_ctauErr", Xib_ConsXib_ctauErr, &b_Xib_ConsXib_ctauErr);
   fChain->SetBranchAddress("Xib_ConsXib_decayLength", Xib_ConsXib_decayLength, &b_Xib_ConsXib_decayLength);
   fChain->SetBranchAddress("Xib_ConsXib_decayLengthErr", Xib_ConsXib_decayLengthErr, &b_Xib_ConsXib_decayLengthErr);
   fChain->SetBranchAddress("Xib_ConsXib_nDOF", Xib_ConsXib_nDOF, &b_Xib_ConsXib_nDOF);
   fChain->SetBranchAddress("Xib_ConsXib_nIter", Xib_ConsXib_nIter, &b_Xib_ConsXib_nIter);
   fChain->SetBranchAddress("Xib_ConsXib_status", Xib_ConsXib_status, &b_Xib_ConsXib_status);
   fChain->SetBranchAddress("Xib_L0MuonDecision_Dec", &Xib_L0MuonDecision_Dec, &b_Xib_L0MuonDecision_Dec);
   fChain->SetBranchAddress("Xib_L0MuonDecision_TIS", &Xib_L0MuonDecision_TIS, &b_Xib_L0MuonDecision_TIS);
   fChain->SetBranchAddress("Xib_L0MuonDecision_TOS", &Xib_L0MuonDecision_TOS, &b_Xib_L0MuonDecision_TOS);
   fChain->SetBranchAddress("Xib_L0DiMuonDecision_Dec", &Xib_L0DiMuonDecision_Dec, &b_Xib_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("Xib_L0DiMuonDecision_TIS", &Xib_L0DiMuonDecision_TIS, &b_Xib_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("Xib_L0DiMuonDecision_TOS", &Xib_L0DiMuonDecision_TOS, &b_Xib_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("Xib_Hlt1SingleMuonNoIPDecision_Dec", &Xib_Hlt1SingleMuonNoIPDecision_Dec, &b_Xib_Hlt1SingleMuonNoIPDecision_Dec);
   fChain->SetBranchAddress("Xib_Hlt1SingleMuonNoIPDecision_TIS", &Xib_Hlt1SingleMuonNoIPDecision_TIS, &b_Xib_Hlt1SingleMuonNoIPDecision_TIS);
   fChain->SetBranchAddress("Xib_Hlt1SingleMuonNoIPDecision_TOS", &Xib_Hlt1SingleMuonNoIPDecision_TOS, &b_Xib_Hlt1SingleMuonNoIPDecision_TOS);
   fChain->SetBranchAddress("Xib_Hlt1DiMuonHighMassDecision_Dec", &Xib_Hlt1DiMuonHighMassDecision_Dec, &b_Xib_Hlt1DiMuonHighMassDecision_Dec);
   fChain->SetBranchAddress("Xib_Hlt1DiMuonHighMassDecision_TIS", &Xib_Hlt1DiMuonHighMassDecision_TIS, &b_Xib_Hlt1DiMuonHighMassDecision_TIS);
   fChain->SetBranchAddress("Xib_Hlt1DiMuonHighMassDecision_TOS", &Xib_Hlt1DiMuonHighMassDecision_TOS, &b_Xib_Hlt1DiMuonHighMassDecision_TOS);
   fChain->SetBranchAddress("Xib_Hlt1TrackAllL0Decision_Dec", &Xib_Hlt1TrackAllL0Decision_Dec, &b_Xib_Hlt1TrackAllL0Decision_Dec);
   fChain->SetBranchAddress("Xib_Hlt1TrackAllL0Decision_TIS", &Xib_Hlt1TrackAllL0Decision_TIS, &b_Xib_Hlt1TrackAllL0Decision_TIS);
   fChain->SetBranchAddress("Xib_Hlt1TrackAllL0Decision_TOS", &Xib_Hlt1TrackAllL0Decision_TOS, &b_Xib_Hlt1TrackAllL0Decision_TOS);
   fChain->SetBranchAddress("Xib_Hlt1TrackMuonDecision_Dec", &Xib_Hlt1TrackMuonDecision_Dec, &b_Xib_Hlt1TrackMuonDecision_Dec);
   fChain->SetBranchAddress("Xib_Hlt1TrackMuonDecision_TIS", &Xib_Hlt1TrackMuonDecision_TIS, &b_Xib_Hlt1TrackMuonDecision_TIS);
   fChain->SetBranchAddress("Xib_Hlt1TrackMuonDecision_TOS", &Xib_Hlt1TrackMuonDecision_TOS, &b_Xib_Hlt1TrackMuonDecision_TOS);
   fChain->SetBranchAddress("Xib_Hlt2DiMuonDetachedHeavyDecision_Dec", &Xib_Hlt2DiMuonDetachedHeavyDecision_Dec, &b_Xib_Hlt2DiMuonDetachedHeavyDecision_Dec);
   fChain->SetBranchAddress("Xib_Hlt2DiMuonDetachedHeavyDecision_TIS", &Xib_Hlt2DiMuonDetachedHeavyDecision_TIS, &b_Xib_Hlt2DiMuonDetachedHeavyDecision_TIS);
   fChain->SetBranchAddress("Xib_Hlt2DiMuonDetachedHeavyDecision_TOS", &Xib_Hlt2DiMuonDetachedHeavyDecision_TOS, &b_Xib_Hlt2DiMuonDetachedHeavyDecision_TOS);
   fChain->SetBranchAddress("Xib_Hlt2DiMuonDetachedJPsiDecision_Dec", &Xib_Hlt2DiMuonDetachedJPsiDecision_Dec, &b_Xib_Hlt2DiMuonDetachedJPsiDecision_Dec);
   fChain->SetBranchAddress("Xib_Hlt2DiMuonDetachedJPsiDecision_TIS", &Xib_Hlt2DiMuonDetachedJPsiDecision_TIS, &b_Xib_Hlt2DiMuonDetachedJPsiDecision_TIS);
   fChain->SetBranchAddress("Xib_Hlt2DiMuonDetachedJPsiDecision_TOS", &Xib_Hlt2DiMuonDetachedJPsiDecision_TOS, &b_Xib_Hlt2DiMuonDetachedJPsiDecision_TOS);
   fChain->SetBranchAddress("Jpsi_CosTheta", &Jpsi_CosTheta, &b_Jpsi_CosTheta);
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
   fChain->SetBranchAddress("Jpsi_ORIVX_X", &Jpsi_ORIVX_X, &b_Jpsi_ORIVX_X);
   fChain->SetBranchAddress("Jpsi_ORIVX_Y", &Jpsi_ORIVX_Y, &b_Jpsi_ORIVX_Y);
   fChain->SetBranchAddress("Jpsi_ORIVX_Z", &Jpsi_ORIVX_Z, &b_Jpsi_ORIVX_Z);
   fChain->SetBranchAddress("Jpsi_ORIVX_XERR", &Jpsi_ORIVX_XERR, &b_Jpsi_ORIVX_XERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_YERR", &Jpsi_ORIVX_YERR, &b_Jpsi_ORIVX_YERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_ZERR", &Jpsi_ORIVX_ZERR, &b_Jpsi_ORIVX_ZERR);
   fChain->SetBranchAddress("Jpsi_ORIVX_CHI2", &Jpsi_ORIVX_CHI2, &b_Jpsi_ORIVX_CHI2);
   fChain->SetBranchAddress("Jpsi_ORIVX_NDOF", &Jpsi_ORIVX_NDOF, &b_Jpsi_ORIVX_NDOF);
   fChain->SetBranchAddress("Jpsi_ORIVX_COV_", Jpsi_ORIVX_COV_, &b_Jpsi_ORIVX_COV_);
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
   fChain->SetBranchAddress("muplus_ORIVX_X", &muplus_ORIVX_X, &b_muplus_ORIVX_X);
   fChain->SetBranchAddress("muplus_ORIVX_Y", &muplus_ORIVX_Y, &b_muplus_ORIVX_Y);
   fChain->SetBranchAddress("muplus_ORIVX_Z", &muplus_ORIVX_Z, &b_muplus_ORIVX_Z);
   fChain->SetBranchAddress("muplus_ORIVX_XERR", &muplus_ORIVX_XERR, &b_muplus_ORIVX_XERR);
   fChain->SetBranchAddress("muplus_ORIVX_YERR", &muplus_ORIVX_YERR, &b_muplus_ORIVX_YERR);
   fChain->SetBranchAddress("muplus_ORIVX_ZERR", &muplus_ORIVX_ZERR, &b_muplus_ORIVX_ZERR);
   fChain->SetBranchAddress("muplus_ORIVX_CHI2", &muplus_ORIVX_CHI2, &b_muplus_ORIVX_CHI2);
   fChain->SetBranchAddress("muplus_ORIVX_NDOF", &muplus_ORIVX_NDOF, &b_muplus_ORIVX_NDOF);
   fChain->SetBranchAddress("muplus_ORIVX_COV_", muplus_ORIVX_COV_, &b_muplus_ORIVX_COV_);
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
   fChain->SetBranchAddress("muminus_ORIVX_X", &muminus_ORIVX_X, &b_muminus_ORIVX_X);
   fChain->SetBranchAddress("muminus_ORIVX_Y", &muminus_ORIVX_Y, &b_muminus_ORIVX_Y);
   fChain->SetBranchAddress("muminus_ORIVX_Z", &muminus_ORIVX_Z, &b_muminus_ORIVX_Z);
   fChain->SetBranchAddress("muminus_ORIVX_XERR", &muminus_ORIVX_XERR, &b_muminus_ORIVX_XERR);
   fChain->SetBranchAddress("muminus_ORIVX_YERR", &muminus_ORIVX_YERR, &b_muminus_ORIVX_YERR);
   fChain->SetBranchAddress("muminus_ORIVX_ZERR", &muminus_ORIVX_ZERR, &b_muminus_ORIVX_ZERR);
   fChain->SetBranchAddress("muminus_ORIVX_CHI2", &muminus_ORIVX_CHI2, &b_muminus_ORIVX_CHI2);
   fChain->SetBranchAddress("muminus_ORIVX_NDOF", &muminus_ORIVX_NDOF, &b_muminus_ORIVX_NDOF);
   fChain->SetBranchAddress("muminus_ORIVX_COV_", muminus_ORIVX_COV_, &b_muminus_ORIVX_COV_);
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
   fChain->SetBranchAddress("Xi_CosTheta", &Xi_CosTheta, &b_Xi_CosTheta);
   fChain->SetBranchAddress("Xi_ENDVERTEX_X", &Xi_ENDVERTEX_X, &b_Xi_ENDVERTEX_X);
   fChain->SetBranchAddress("Xi_ENDVERTEX_Y", &Xi_ENDVERTEX_Y, &b_Xi_ENDVERTEX_Y);
   fChain->SetBranchAddress("Xi_ENDVERTEX_Z", &Xi_ENDVERTEX_Z, &b_Xi_ENDVERTEX_Z);
   fChain->SetBranchAddress("Xi_ENDVERTEX_XERR", &Xi_ENDVERTEX_XERR, &b_Xi_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Xi_ENDVERTEX_YERR", &Xi_ENDVERTEX_YERR, &b_Xi_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Xi_ENDVERTEX_ZERR", &Xi_ENDVERTEX_ZERR, &b_Xi_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Xi_ENDVERTEX_CHI2", &Xi_ENDVERTEX_CHI2, &b_Xi_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Xi_ENDVERTEX_NDOF", &Xi_ENDVERTEX_NDOF, &b_Xi_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Xi_ENDVERTEX_COV_", Xi_ENDVERTEX_COV_, &b_Xi_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Xi_OWNPV_X", &Xi_OWNPV_X, &b_Xi_OWNPV_X);
   fChain->SetBranchAddress("Xi_OWNPV_Y", &Xi_OWNPV_Y, &b_Xi_OWNPV_Y);
   fChain->SetBranchAddress("Xi_OWNPV_Z", &Xi_OWNPV_Z, &b_Xi_OWNPV_Z);
   fChain->SetBranchAddress("Xi_OWNPV_XERR", &Xi_OWNPV_XERR, &b_Xi_OWNPV_XERR);
   fChain->SetBranchAddress("Xi_OWNPV_YERR", &Xi_OWNPV_YERR, &b_Xi_OWNPV_YERR);
   fChain->SetBranchAddress("Xi_OWNPV_ZERR", &Xi_OWNPV_ZERR, &b_Xi_OWNPV_ZERR);
   fChain->SetBranchAddress("Xi_OWNPV_CHI2", &Xi_OWNPV_CHI2, &b_Xi_OWNPV_CHI2);
   fChain->SetBranchAddress("Xi_OWNPV_NDOF", &Xi_OWNPV_NDOF, &b_Xi_OWNPV_NDOF);
   fChain->SetBranchAddress("Xi_OWNPV_COV_", Xi_OWNPV_COV_, &b_Xi_OWNPV_COV_);
   fChain->SetBranchAddress("Xi_IP_OWNPV", &Xi_IP_OWNPV, &b_Xi_IP_OWNPV);
   fChain->SetBranchAddress("Xi_IPCHI2_OWNPV", &Xi_IPCHI2_OWNPV, &b_Xi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Xi_FD_OWNPV", &Xi_FD_OWNPV, &b_Xi_FD_OWNPV);
   fChain->SetBranchAddress("Xi_FDCHI2_OWNPV", &Xi_FDCHI2_OWNPV, &b_Xi_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Xi_DIRA_OWNPV", &Xi_DIRA_OWNPV, &b_Xi_DIRA_OWNPV);
   fChain->SetBranchAddress("Xi_ORIVX_X", &Xi_ORIVX_X, &b_Xi_ORIVX_X);
   fChain->SetBranchAddress("Xi_ORIVX_Y", &Xi_ORIVX_Y, &b_Xi_ORIVX_Y);
   fChain->SetBranchAddress("Xi_ORIVX_Z", &Xi_ORIVX_Z, &b_Xi_ORIVX_Z);
   fChain->SetBranchAddress("Xi_ORIVX_XERR", &Xi_ORIVX_XERR, &b_Xi_ORIVX_XERR);
   fChain->SetBranchAddress("Xi_ORIVX_YERR", &Xi_ORIVX_YERR, &b_Xi_ORIVX_YERR);
   fChain->SetBranchAddress("Xi_ORIVX_ZERR", &Xi_ORIVX_ZERR, &b_Xi_ORIVX_ZERR);
   fChain->SetBranchAddress("Xi_ORIVX_CHI2", &Xi_ORIVX_CHI2, &b_Xi_ORIVX_CHI2);
   fChain->SetBranchAddress("Xi_ORIVX_NDOF", &Xi_ORIVX_NDOF, &b_Xi_ORIVX_NDOF);
   fChain->SetBranchAddress("Xi_ORIVX_COV_", Xi_ORIVX_COV_, &b_Xi_ORIVX_COV_);
   fChain->SetBranchAddress("Xi_FD_ORIVX", &Xi_FD_ORIVX, &b_Xi_FD_ORIVX);
   fChain->SetBranchAddress("Xi_FDCHI2_ORIVX", &Xi_FDCHI2_ORIVX, &b_Xi_FDCHI2_ORIVX);
   fChain->SetBranchAddress("Xi_DIRA_ORIVX", &Xi_DIRA_ORIVX, &b_Xi_DIRA_ORIVX);
   fChain->SetBranchAddress("Xi_P", &Xi_P, &b_Xi_P);
   fChain->SetBranchAddress("Xi_PT", &Xi_PT, &b_Xi_PT);
   fChain->SetBranchAddress("Xi_PE", &Xi_PE, &b_Xi_PE);
   fChain->SetBranchAddress("Xi_PX", &Xi_PX, &b_Xi_PX);
   fChain->SetBranchAddress("Xi_PY", &Xi_PY, &b_Xi_PY);
   fChain->SetBranchAddress("Xi_PZ", &Xi_PZ, &b_Xi_PZ);
   fChain->SetBranchAddress("Xi_MM", &Xi_MM, &b_Xi_MM);
   fChain->SetBranchAddress("Xi_MMERR", &Xi_MMERR, &b_Xi_MMERR);
   fChain->SetBranchAddress("Xi_M", &Xi_M, &b_Xi_M);
   fChain->SetBranchAddress("Xi_ID", &Xi_ID, &b_Xi_ID);
   fChain->SetBranchAddress("Xi_TAU", &Xi_TAU, &b_Xi_TAU);
   fChain->SetBranchAddress("Xi_TAUERR", &Xi_TAUERR, &b_Xi_TAUERR);
   fChain->SetBranchAddress("Xi_TAUCHI2", &Xi_TAUCHI2, &b_Xi_TAUCHI2);
   fChain->SetBranchAddress("Xi_L0Global_Dec", &Xi_L0Global_Dec, &b_Xi_L0Global_Dec);
   fChain->SetBranchAddress("Xi_L0Global_TIS", &Xi_L0Global_TIS, &b_Xi_L0Global_TIS);
   fChain->SetBranchAddress("Xi_L0Global_TOS", &Xi_L0Global_TOS, &b_Xi_L0Global_TOS);
   fChain->SetBranchAddress("Xi_Hlt1Global_Dec", &Xi_Hlt1Global_Dec, &b_Xi_Hlt1Global_Dec);
   fChain->SetBranchAddress("Xi_Hlt1Global_TIS", &Xi_Hlt1Global_TIS, &b_Xi_Hlt1Global_TIS);
   fChain->SetBranchAddress("Xi_Hlt1Global_TOS", &Xi_Hlt1Global_TOS, &b_Xi_Hlt1Global_TOS);
   fChain->SetBranchAddress("Xi_Hlt1Phys_Dec", &Xi_Hlt1Phys_Dec, &b_Xi_Hlt1Phys_Dec);
   fChain->SetBranchAddress("Xi_Hlt1Phys_TIS", &Xi_Hlt1Phys_TIS, &b_Xi_Hlt1Phys_TIS);
   fChain->SetBranchAddress("Xi_Hlt1Phys_TOS", &Xi_Hlt1Phys_TOS, &b_Xi_Hlt1Phys_TOS);
   fChain->SetBranchAddress("Xi_Hlt2Global_Dec", &Xi_Hlt2Global_Dec, &b_Xi_Hlt2Global_Dec);
   fChain->SetBranchAddress("Xi_Hlt2Global_TIS", &Xi_Hlt2Global_TIS, &b_Xi_Hlt2Global_TIS);
   fChain->SetBranchAddress("Xi_Hlt2Global_TOS", &Xi_Hlt2Global_TOS, &b_Xi_Hlt2Global_TOS);
   fChain->SetBranchAddress("Xi_Hlt2Phys_Dec", &Xi_Hlt2Phys_Dec, &b_Xi_Hlt2Phys_Dec);
   fChain->SetBranchAddress("Xi_Hlt2Phys_TIS", &Xi_Hlt2Phys_TIS, &b_Xi_Hlt2Phys_TIS);
   fChain->SetBranchAddress("Xi_Hlt2Phys_TOS", &Xi_Hlt2Phys_TOS, &b_Xi_Hlt2Phys_TOS);
   fChain->SetBranchAddress("Xi_ETA", &Xi_ETA, &b_Xi_ETA);
   fChain->SetBranchAddress("Xi_PHI", &Xi_PHI, &b_Xi_PHI);
   fChain->SetBranchAddress("Xi_DTF_L_ARMEN_JpsiXiConstr", &Xi_DTF_L_ARMEN_JpsiXiConstr, &b_Xi_DTF_L_ARMEN_JpsiXiConstr);
   fChain->SetBranchAddress("Xi_DTF_L_ARMEN_JpsiXiLConstr", &Xi_DTF_L_ARMEN_JpsiXiLConstr, &b_Xi_DTF_L_ARMEN_JpsiXiLConstr);
   fChain->SetBranchAddress("Xi_DTF_L_LV01_JpsiXiConstr", &Xi_DTF_L_LV01_JpsiXiConstr, &b_Xi_DTF_L_LV01_JpsiXiConstr);
   fChain->SetBranchAddress("Xi_DTF_L_LV01_JpsiXiLConstr", &Xi_DTF_L_LV01_JpsiXiLConstr, &b_Xi_DTF_L_LV01_JpsiXiLConstr);
   fChain->SetBranchAddress("Xi_DTF_L_WMpipi_JpsiXiConstr", &Xi_DTF_L_WMpipi_JpsiXiConstr, &b_Xi_DTF_L_WMpipi_JpsiXiConstr);
   fChain->SetBranchAddress("Xi_DTF_L_WMpipi_JpsiXiLConstr", &Xi_DTF_L_WMpipi_JpsiXiLConstr, &b_Xi_DTF_L_WMpipi_JpsiXiLConstr);
   fChain->SetBranchAddress("L_CosTheta", &L_CosTheta, &b_L_CosTheta);
   fChain->SetBranchAddress("L_ENDVERTEX_X", &L_ENDVERTEX_X, &b_L_ENDVERTEX_X);
   fChain->SetBranchAddress("L_ENDVERTEX_Y", &L_ENDVERTEX_Y, &b_L_ENDVERTEX_Y);
   fChain->SetBranchAddress("L_ENDVERTEX_Z", &L_ENDVERTEX_Z, &b_L_ENDVERTEX_Z);
   fChain->SetBranchAddress("L_ENDVERTEX_XERR", &L_ENDVERTEX_XERR, &b_L_ENDVERTEX_XERR);
   fChain->SetBranchAddress("L_ENDVERTEX_YERR", &L_ENDVERTEX_YERR, &b_L_ENDVERTEX_YERR);
   fChain->SetBranchAddress("L_ENDVERTEX_ZERR", &L_ENDVERTEX_ZERR, &b_L_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("L_ENDVERTEX_CHI2", &L_ENDVERTEX_CHI2, &b_L_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("L_ENDVERTEX_NDOF", &L_ENDVERTEX_NDOF, &b_L_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("L_ENDVERTEX_COV_", L_ENDVERTEX_COV_, &b_L_ENDVERTEX_COV_);
   fChain->SetBranchAddress("L_OWNPV_X", &L_OWNPV_X, &b_L_OWNPV_X);
   fChain->SetBranchAddress("L_OWNPV_Y", &L_OWNPV_Y, &b_L_OWNPV_Y);
   fChain->SetBranchAddress("L_OWNPV_Z", &L_OWNPV_Z, &b_L_OWNPV_Z);
   fChain->SetBranchAddress("L_OWNPV_XERR", &L_OWNPV_XERR, &b_L_OWNPV_XERR);
   fChain->SetBranchAddress("L_OWNPV_YERR", &L_OWNPV_YERR, &b_L_OWNPV_YERR);
   fChain->SetBranchAddress("L_OWNPV_ZERR", &L_OWNPV_ZERR, &b_L_OWNPV_ZERR);
   fChain->SetBranchAddress("L_OWNPV_CHI2", &L_OWNPV_CHI2, &b_L_OWNPV_CHI2);
   fChain->SetBranchAddress("L_OWNPV_NDOF", &L_OWNPV_NDOF, &b_L_OWNPV_NDOF);
   fChain->SetBranchAddress("L_OWNPV_COV_", L_OWNPV_COV_, &b_L_OWNPV_COV_);
   fChain->SetBranchAddress("L_IP_OWNPV", &L_IP_OWNPV, &b_L_IP_OWNPV);
   fChain->SetBranchAddress("L_IPCHI2_OWNPV", &L_IPCHI2_OWNPV, &b_L_IPCHI2_OWNPV);
   fChain->SetBranchAddress("L_FD_OWNPV", &L_FD_OWNPV, &b_L_FD_OWNPV);
   fChain->SetBranchAddress("L_FDCHI2_OWNPV", &L_FDCHI2_OWNPV, &b_L_FDCHI2_OWNPV);
   fChain->SetBranchAddress("L_DIRA_OWNPV", &L_DIRA_OWNPV, &b_L_DIRA_OWNPV);
   fChain->SetBranchAddress("L_ORIVX_X", &L_ORIVX_X, &b_L_ORIVX_X);
   fChain->SetBranchAddress("L_ORIVX_Y", &L_ORIVX_Y, &b_L_ORIVX_Y);
   fChain->SetBranchAddress("L_ORIVX_Z", &L_ORIVX_Z, &b_L_ORIVX_Z);
   fChain->SetBranchAddress("L_ORIVX_XERR", &L_ORIVX_XERR, &b_L_ORIVX_XERR);
   fChain->SetBranchAddress("L_ORIVX_YERR", &L_ORIVX_YERR, &b_L_ORIVX_YERR);
   fChain->SetBranchAddress("L_ORIVX_ZERR", &L_ORIVX_ZERR, &b_L_ORIVX_ZERR);
   fChain->SetBranchAddress("L_ORIVX_CHI2", &L_ORIVX_CHI2, &b_L_ORIVX_CHI2);
   fChain->SetBranchAddress("L_ORIVX_NDOF", &L_ORIVX_NDOF, &b_L_ORIVX_NDOF);
   fChain->SetBranchAddress("L_ORIVX_COV_", L_ORIVX_COV_, &b_L_ORIVX_COV_);
   fChain->SetBranchAddress("L_FD_ORIVX", &L_FD_ORIVX, &b_L_FD_ORIVX);
   fChain->SetBranchAddress("L_FDCHI2_ORIVX", &L_FDCHI2_ORIVX, &b_L_FDCHI2_ORIVX);
   fChain->SetBranchAddress("L_DIRA_ORIVX", &L_DIRA_ORIVX, &b_L_DIRA_ORIVX);
   fChain->SetBranchAddress("L_P", &L_P, &b_L_P);
   fChain->SetBranchAddress("L_PT", &L_PT, &b_L_PT);
   fChain->SetBranchAddress("L_PE", &L_PE, &b_L_PE);
   fChain->SetBranchAddress("L_PX", &L_PX, &b_L_PX);
   fChain->SetBranchAddress("L_PY", &L_PY, &b_L_PY);
   fChain->SetBranchAddress("L_PZ", &L_PZ, &b_L_PZ);
   fChain->SetBranchAddress("L_MM", &L_MM, &b_L_MM);
   fChain->SetBranchAddress("L_MMERR", &L_MMERR, &b_L_MMERR);
   fChain->SetBranchAddress("L_M", &L_M, &b_L_M);
   fChain->SetBranchAddress("L_ID", &L_ID, &b_L_ID);
   fChain->SetBranchAddress("L_TAU", &L_TAU, &b_L_TAU);
   fChain->SetBranchAddress("L_TAUERR", &L_TAUERR, &b_L_TAUERR);
   fChain->SetBranchAddress("L_TAUCHI2", &L_TAUCHI2, &b_L_TAUCHI2);
   fChain->SetBranchAddress("L_L0Global_Dec", &L_L0Global_Dec, &b_L_L0Global_Dec);
   fChain->SetBranchAddress("L_L0Global_TIS", &L_L0Global_TIS, &b_L_L0Global_TIS);
   fChain->SetBranchAddress("L_L0Global_TOS", &L_L0Global_TOS, &b_L_L0Global_TOS);
   fChain->SetBranchAddress("L_Hlt1Global_Dec", &L_Hlt1Global_Dec, &b_L_Hlt1Global_Dec);
   fChain->SetBranchAddress("L_Hlt1Global_TIS", &L_Hlt1Global_TIS, &b_L_Hlt1Global_TIS);
   fChain->SetBranchAddress("L_Hlt1Global_TOS", &L_Hlt1Global_TOS, &b_L_Hlt1Global_TOS);
   fChain->SetBranchAddress("L_Hlt1Phys_Dec", &L_Hlt1Phys_Dec, &b_L_Hlt1Phys_Dec);
   fChain->SetBranchAddress("L_Hlt1Phys_TIS", &L_Hlt1Phys_TIS, &b_L_Hlt1Phys_TIS);
   fChain->SetBranchAddress("L_Hlt1Phys_TOS", &L_Hlt1Phys_TOS, &b_L_Hlt1Phys_TOS);
   fChain->SetBranchAddress("L_Hlt2Global_Dec", &L_Hlt2Global_Dec, &b_L_Hlt2Global_Dec);
   fChain->SetBranchAddress("L_Hlt2Global_TIS", &L_Hlt2Global_TIS, &b_L_Hlt2Global_TIS);
   fChain->SetBranchAddress("L_Hlt2Global_TOS", &L_Hlt2Global_TOS, &b_L_Hlt2Global_TOS);
   fChain->SetBranchAddress("L_Hlt2Phys_Dec", &L_Hlt2Phys_Dec, &b_L_Hlt2Phys_Dec);
   fChain->SetBranchAddress("L_Hlt2Phys_TIS", &L_Hlt2Phys_TIS, &b_L_Hlt2Phys_TIS);
   fChain->SetBranchAddress("L_Hlt2Phys_TOS", &L_Hlt2Phys_TOS, &b_L_Hlt2Phys_TOS);
   fChain->SetBranchAddress("L_ETA", &L_ETA, &b_L_ETA);
   fChain->SetBranchAddress("L_PHI", &L_PHI, &b_L_PHI);
   fChain->SetBranchAddress("L_ARMEN", &L_ARMEN, &b_L_ARMEN);
   fChain->SetBranchAddress("L_BPVLTIME", &L_BPVLTIME, &b_L_BPVLTIME);
   fChain->SetBranchAddress("L_WMkk", &L_WMkk, &b_L_WMkk);
   fChain->SetBranchAddress("L_WMkpi", &L_WMkpi, &b_L_WMkpi);
   fChain->SetBranchAddress("L_WMpik", &L_WMpik, &b_L_WMpik);
   fChain->SetBranchAddress("L_WMpipi", &L_WMpipi, &b_L_WMpipi);
   fChain->SetBranchAddress("L_WMpp", &L_WMpp, &b_L_WMpp);
   fChain->SetBranchAddress("L_dm", &L_dm, &b_L_dm);
   fChain->SetBranchAddress("L_lv01", &L_lv01, &b_L_lv01);
   fChain->SetBranchAddress("L_lv02", &L_lv02, &b_L_lv02);
   fChain->SetBranchAddress("p_CosTheta", &p_CosTheta, &b_p_CosTheta);
   fChain->SetBranchAddress("p_OWNPV_X", &p_OWNPV_X, &b_p_OWNPV_X);
   fChain->SetBranchAddress("p_OWNPV_Y", &p_OWNPV_Y, &b_p_OWNPV_Y);
   fChain->SetBranchAddress("p_OWNPV_Z", &p_OWNPV_Z, &b_p_OWNPV_Z);
   fChain->SetBranchAddress("p_OWNPV_XERR", &p_OWNPV_XERR, &b_p_OWNPV_XERR);
   fChain->SetBranchAddress("p_OWNPV_YERR", &p_OWNPV_YERR, &b_p_OWNPV_YERR);
   fChain->SetBranchAddress("p_OWNPV_ZERR", &p_OWNPV_ZERR, &b_p_OWNPV_ZERR);
   fChain->SetBranchAddress("p_OWNPV_CHI2", &p_OWNPV_CHI2, &b_p_OWNPV_CHI2);
   fChain->SetBranchAddress("p_OWNPV_NDOF", &p_OWNPV_NDOF, &b_p_OWNPV_NDOF);
   fChain->SetBranchAddress("p_OWNPV_COV_", p_OWNPV_COV_, &b_p_OWNPV_COV_);
   fChain->SetBranchAddress("p_IP_OWNPV", &p_IP_OWNPV, &b_p_IP_OWNPV);
   fChain->SetBranchAddress("p_IPCHI2_OWNPV", &p_IPCHI2_OWNPV, &b_p_IPCHI2_OWNPV);
   fChain->SetBranchAddress("p_ORIVX_X", &p_ORIVX_X, &b_p_ORIVX_X);
   fChain->SetBranchAddress("p_ORIVX_Y", &p_ORIVX_Y, &b_p_ORIVX_Y);
   fChain->SetBranchAddress("p_ORIVX_Z", &p_ORIVX_Z, &b_p_ORIVX_Z);
   fChain->SetBranchAddress("p_ORIVX_XERR", &p_ORIVX_XERR, &b_p_ORIVX_XERR);
   fChain->SetBranchAddress("p_ORIVX_YERR", &p_ORIVX_YERR, &b_p_ORIVX_YERR);
   fChain->SetBranchAddress("p_ORIVX_ZERR", &p_ORIVX_ZERR, &b_p_ORIVX_ZERR);
   fChain->SetBranchAddress("p_ORIVX_CHI2", &p_ORIVX_CHI2, &b_p_ORIVX_CHI2);
   fChain->SetBranchAddress("p_ORIVX_NDOF", &p_ORIVX_NDOF, &b_p_ORIVX_NDOF);
   fChain->SetBranchAddress("p_ORIVX_COV_", p_ORIVX_COV_, &b_p_ORIVX_COV_);
   fChain->SetBranchAddress("p_P", &p_P, &b_p_P);
   fChain->SetBranchAddress("p_PT", &p_PT, &b_p_PT);
   fChain->SetBranchAddress("p_PE", &p_PE, &b_p_PE);
   fChain->SetBranchAddress("p_PX", &p_PX, &b_p_PX);
   fChain->SetBranchAddress("p_PY", &p_PY, &b_p_PY);
   fChain->SetBranchAddress("p_PZ", &p_PZ, &b_p_PZ);
   fChain->SetBranchAddress("p_M", &p_M, &b_p_M);
   fChain->SetBranchAddress("p_ID", &p_ID, &b_p_ID);
   fChain->SetBranchAddress("p_PIDe", &p_PIDe, &b_p_PIDe);
   fChain->SetBranchAddress("p_PIDmu", &p_PIDmu, &b_p_PIDmu);
   fChain->SetBranchAddress("p_PIDK", &p_PIDK, &b_p_PIDK);
   fChain->SetBranchAddress("p_PIDp", &p_PIDp, &b_p_PIDp);
   fChain->SetBranchAddress("p_ProbNNe", &p_ProbNNe, &b_p_ProbNNe);
   fChain->SetBranchAddress("p_ProbNNk", &p_ProbNNk, &b_p_ProbNNk);
   fChain->SetBranchAddress("p_ProbNNp", &p_ProbNNp, &b_p_ProbNNp);
   fChain->SetBranchAddress("p_ProbNNpi", &p_ProbNNpi, &b_p_ProbNNpi);
   fChain->SetBranchAddress("p_ProbNNmu", &p_ProbNNmu, &b_p_ProbNNmu);
   fChain->SetBranchAddress("p_ProbNNghost", &p_ProbNNghost, &b_p_ProbNNghost);
   fChain->SetBranchAddress("p_hasMuon", &p_hasMuon, &b_p_hasMuon);
   fChain->SetBranchAddress("p_isMuon", &p_isMuon, &b_p_isMuon);
   fChain->SetBranchAddress("p_hasRich", &p_hasRich, &b_p_hasRich);
   fChain->SetBranchAddress("p_UsedRichAerogel", &p_UsedRichAerogel, &b_p_UsedRichAerogel);
   fChain->SetBranchAddress("p_UsedRich1Gas", &p_UsedRich1Gas, &b_p_UsedRich1Gas);
   fChain->SetBranchAddress("p_UsedRich2Gas", &p_UsedRich2Gas, &b_p_UsedRich2Gas);
   fChain->SetBranchAddress("p_RichAboveElThres", &p_RichAboveElThres, &b_p_RichAboveElThres);
   fChain->SetBranchAddress("p_RichAboveMuThres", &p_RichAboveMuThres, &b_p_RichAboveMuThres);
   fChain->SetBranchAddress("p_RichAbovePiThres", &p_RichAbovePiThres, &b_p_RichAbovePiThres);
   fChain->SetBranchAddress("p_RichAboveKaThres", &p_RichAboveKaThres, &b_p_RichAboveKaThres);
   fChain->SetBranchAddress("p_RichAbovePrThres", &p_RichAbovePrThres, &b_p_RichAbovePrThres);
   fChain->SetBranchAddress("p_hasCalo", &p_hasCalo, &b_p_hasCalo);
   fChain->SetBranchAddress("p_L0Global_Dec", &p_L0Global_Dec, &b_p_L0Global_Dec);
   fChain->SetBranchAddress("p_L0Global_TIS", &p_L0Global_TIS, &b_p_L0Global_TIS);
   fChain->SetBranchAddress("p_L0Global_TOS", &p_L0Global_TOS, &b_p_L0Global_TOS);
   fChain->SetBranchAddress("p_Hlt1Global_Dec", &p_Hlt1Global_Dec, &b_p_Hlt1Global_Dec);
   fChain->SetBranchAddress("p_Hlt1Global_TIS", &p_Hlt1Global_TIS, &b_p_Hlt1Global_TIS);
   fChain->SetBranchAddress("p_Hlt1Global_TOS", &p_Hlt1Global_TOS, &b_p_Hlt1Global_TOS);
   fChain->SetBranchAddress("p_Hlt1Phys_Dec", &p_Hlt1Phys_Dec, &b_p_Hlt1Phys_Dec);
   fChain->SetBranchAddress("p_Hlt1Phys_TIS", &p_Hlt1Phys_TIS, &b_p_Hlt1Phys_TIS);
   fChain->SetBranchAddress("p_Hlt1Phys_TOS", &p_Hlt1Phys_TOS, &b_p_Hlt1Phys_TOS);
   fChain->SetBranchAddress("p_Hlt2Global_Dec", &p_Hlt2Global_Dec, &b_p_Hlt2Global_Dec);
   fChain->SetBranchAddress("p_Hlt2Global_TIS", &p_Hlt2Global_TIS, &b_p_Hlt2Global_TIS);
   fChain->SetBranchAddress("p_Hlt2Global_TOS", &p_Hlt2Global_TOS, &b_p_Hlt2Global_TOS);
   fChain->SetBranchAddress("p_Hlt2Phys_Dec", &p_Hlt2Phys_Dec, &b_p_Hlt2Phys_Dec);
   fChain->SetBranchAddress("p_Hlt2Phys_TIS", &p_Hlt2Phys_TIS, &b_p_Hlt2Phys_TIS);
   fChain->SetBranchAddress("p_Hlt2Phys_TOS", &p_Hlt2Phys_TOS, &b_p_Hlt2Phys_TOS);
   fChain->SetBranchAddress("p_TRACK_Type", &p_TRACK_Type, &b_p_TRACK_Type);
   fChain->SetBranchAddress("p_TRACK_Key", &p_TRACK_Key, &b_p_TRACK_Key);
   fChain->SetBranchAddress("p_TRACK_CHI2NDOF", &p_TRACK_CHI2NDOF, &b_p_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("p_TRACK_PCHI2", &p_TRACK_PCHI2, &b_p_TRACK_PCHI2);
   fChain->SetBranchAddress("p_TRACK_MatchCHI2", &p_TRACK_MatchCHI2, &b_p_TRACK_MatchCHI2);
   fChain->SetBranchAddress("p_TRACK_GhostProb", &p_TRACK_GhostProb, &b_p_TRACK_GhostProb);
   fChain->SetBranchAddress("p_TRACK_CloneDist", &p_TRACK_CloneDist, &b_p_TRACK_CloneDist);
   fChain->SetBranchAddress("p_TRACK_Likelihood", &p_TRACK_Likelihood, &b_p_TRACK_Likelihood);
   fChain->SetBranchAddress("p_ETA", &p_ETA, &b_p_ETA);
   fChain->SetBranchAddress("p_PHI", &p_PHI, &b_p_PHI);
   fChain->SetBranchAddress("pi_CosTheta", &pi_CosTheta, &b_pi_CosTheta);
   fChain->SetBranchAddress("pi_OWNPV_X", &pi_OWNPV_X, &b_pi_OWNPV_X);
   fChain->SetBranchAddress("pi_OWNPV_Y", &pi_OWNPV_Y, &b_pi_OWNPV_Y);
   fChain->SetBranchAddress("pi_OWNPV_Z", &pi_OWNPV_Z, &b_pi_OWNPV_Z);
   fChain->SetBranchAddress("pi_OWNPV_XERR", &pi_OWNPV_XERR, &b_pi_OWNPV_XERR);
   fChain->SetBranchAddress("pi_OWNPV_YERR", &pi_OWNPV_YERR, &b_pi_OWNPV_YERR);
   fChain->SetBranchAddress("pi_OWNPV_ZERR", &pi_OWNPV_ZERR, &b_pi_OWNPV_ZERR);
   fChain->SetBranchAddress("pi_OWNPV_CHI2", &pi_OWNPV_CHI2, &b_pi_OWNPV_CHI2);
   fChain->SetBranchAddress("pi_OWNPV_NDOF", &pi_OWNPV_NDOF, &b_pi_OWNPV_NDOF);
   fChain->SetBranchAddress("pi_OWNPV_COV_", pi_OWNPV_COV_, &b_pi_OWNPV_COV_);
   fChain->SetBranchAddress("pi_IP_OWNPV", &pi_IP_OWNPV, &b_pi_IP_OWNPV);
   fChain->SetBranchAddress("pi_IPCHI2_OWNPV", &pi_IPCHI2_OWNPV, &b_pi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("pi_ORIVX_X", &pi_ORIVX_X, &b_pi_ORIVX_X);
   fChain->SetBranchAddress("pi_ORIVX_Y", &pi_ORIVX_Y, &b_pi_ORIVX_Y);
   fChain->SetBranchAddress("pi_ORIVX_Z", &pi_ORIVX_Z, &b_pi_ORIVX_Z);
   fChain->SetBranchAddress("pi_ORIVX_XERR", &pi_ORIVX_XERR, &b_pi_ORIVX_XERR);
   fChain->SetBranchAddress("pi_ORIVX_YERR", &pi_ORIVX_YERR, &b_pi_ORIVX_YERR);
   fChain->SetBranchAddress("pi_ORIVX_ZERR", &pi_ORIVX_ZERR, &b_pi_ORIVX_ZERR);
   fChain->SetBranchAddress("pi_ORIVX_CHI2", &pi_ORIVX_CHI2, &b_pi_ORIVX_CHI2);
   fChain->SetBranchAddress("pi_ORIVX_NDOF", &pi_ORIVX_NDOF, &b_pi_ORIVX_NDOF);
   fChain->SetBranchAddress("pi_ORIVX_COV_", pi_ORIVX_COV_, &b_pi_ORIVX_COV_);
   fChain->SetBranchAddress("pi_P", &pi_P, &b_pi_P);
   fChain->SetBranchAddress("pi_PT", &pi_PT, &b_pi_PT);
   fChain->SetBranchAddress("pi_PE", &pi_PE, &b_pi_PE);
   fChain->SetBranchAddress("pi_PX", &pi_PX, &b_pi_PX);
   fChain->SetBranchAddress("pi_PY", &pi_PY, &b_pi_PY);
   fChain->SetBranchAddress("pi_PZ", &pi_PZ, &b_pi_PZ);
   fChain->SetBranchAddress("pi_M", &pi_M, &b_pi_M);
   fChain->SetBranchAddress("pi_ID", &pi_ID, &b_pi_ID);
   fChain->SetBranchAddress("pi_PIDe", &pi_PIDe, &b_pi_PIDe);
   fChain->SetBranchAddress("pi_PIDmu", &pi_PIDmu, &b_pi_PIDmu);
   fChain->SetBranchAddress("pi_PIDK", &pi_PIDK, &b_pi_PIDK);
   fChain->SetBranchAddress("pi_PIDp", &pi_PIDp, &b_pi_PIDp);
   fChain->SetBranchAddress("pi_ProbNNe", &pi_ProbNNe, &b_pi_ProbNNe);
   fChain->SetBranchAddress("pi_ProbNNk", &pi_ProbNNk, &b_pi_ProbNNk);
   fChain->SetBranchAddress("pi_ProbNNp", &pi_ProbNNp, &b_pi_ProbNNp);
   fChain->SetBranchAddress("pi_ProbNNpi", &pi_ProbNNpi, &b_pi_ProbNNpi);
   fChain->SetBranchAddress("pi_ProbNNmu", &pi_ProbNNmu, &b_pi_ProbNNmu);
   fChain->SetBranchAddress("pi_ProbNNghost", &pi_ProbNNghost, &b_pi_ProbNNghost);
   fChain->SetBranchAddress("pi_hasMuon", &pi_hasMuon, &b_pi_hasMuon);
   fChain->SetBranchAddress("pi_isMuon", &pi_isMuon, &b_pi_isMuon);
   fChain->SetBranchAddress("pi_hasRich", &pi_hasRich, &b_pi_hasRich);
   fChain->SetBranchAddress("pi_UsedRichAerogel", &pi_UsedRichAerogel, &b_pi_UsedRichAerogel);
   fChain->SetBranchAddress("pi_UsedRich1Gas", &pi_UsedRich1Gas, &b_pi_UsedRich1Gas);
   fChain->SetBranchAddress("pi_UsedRich2Gas", &pi_UsedRich2Gas, &b_pi_UsedRich2Gas);
   fChain->SetBranchAddress("pi_RichAboveElThres", &pi_RichAboveElThres, &b_pi_RichAboveElThres);
   fChain->SetBranchAddress("pi_RichAboveMuThres", &pi_RichAboveMuThres, &b_pi_RichAboveMuThres);
   fChain->SetBranchAddress("pi_RichAbovePiThres", &pi_RichAbovePiThres, &b_pi_RichAbovePiThres);
   fChain->SetBranchAddress("pi_RichAboveKaThres", &pi_RichAboveKaThres, &b_pi_RichAboveKaThres);
   fChain->SetBranchAddress("pi_RichAbovePrThres", &pi_RichAbovePrThres, &b_pi_RichAbovePrThres);
   fChain->SetBranchAddress("pi_hasCalo", &pi_hasCalo, &b_pi_hasCalo);
   fChain->SetBranchAddress("pi_L0Global_Dec", &pi_L0Global_Dec, &b_pi_L0Global_Dec);
   fChain->SetBranchAddress("pi_L0Global_TIS", &pi_L0Global_TIS, &b_pi_L0Global_TIS);
   fChain->SetBranchAddress("pi_L0Global_TOS", &pi_L0Global_TOS, &b_pi_L0Global_TOS);
   fChain->SetBranchAddress("pi_Hlt1Global_Dec", &pi_Hlt1Global_Dec, &b_pi_Hlt1Global_Dec);
   fChain->SetBranchAddress("pi_Hlt1Global_TIS", &pi_Hlt1Global_TIS, &b_pi_Hlt1Global_TIS);
   fChain->SetBranchAddress("pi_Hlt1Global_TOS", &pi_Hlt1Global_TOS, &b_pi_Hlt1Global_TOS);
   fChain->SetBranchAddress("pi_Hlt1Phys_Dec", &pi_Hlt1Phys_Dec, &b_pi_Hlt1Phys_Dec);
   fChain->SetBranchAddress("pi_Hlt1Phys_TIS", &pi_Hlt1Phys_TIS, &b_pi_Hlt1Phys_TIS);
   fChain->SetBranchAddress("pi_Hlt1Phys_TOS", &pi_Hlt1Phys_TOS, &b_pi_Hlt1Phys_TOS);
   fChain->SetBranchAddress("pi_Hlt2Global_Dec", &pi_Hlt2Global_Dec, &b_pi_Hlt2Global_Dec);
   fChain->SetBranchAddress("pi_Hlt2Global_TIS", &pi_Hlt2Global_TIS, &b_pi_Hlt2Global_TIS);
   fChain->SetBranchAddress("pi_Hlt2Global_TOS", &pi_Hlt2Global_TOS, &b_pi_Hlt2Global_TOS);
   fChain->SetBranchAddress("pi_Hlt2Phys_Dec", &pi_Hlt2Phys_Dec, &b_pi_Hlt2Phys_Dec);
   fChain->SetBranchAddress("pi_Hlt2Phys_TIS", &pi_Hlt2Phys_TIS, &b_pi_Hlt2Phys_TIS);
   fChain->SetBranchAddress("pi_Hlt2Phys_TOS", &pi_Hlt2Phys_TOS, &b_pi_Hlt2Phys_TOS);
   fChain->SetBranchAddress("pi_TRACK_Type", &pi_TRACK_Type, &b_pi_TRACK_Type);
   fChain->SetBranchAddress("pi_TRACK_Key", &pi_TRACK_Key, &b_pi_TRACK_Key);
   fChain->SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_TRACK_CHI2NDOF, &b_pi_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("pi_TRACK_PCHI2", &pi_TRACK_PCHI2, &b_pi_TRACK_PCHI2);
   fChain->SetBranchAddress("pi_TRACK_MatchCHI2", &pi_TRACK_MatchCHI2, &b_pi_TRACK_MatchCHI2);
   fChain->SetBranchAddress("pi_TRACK_GhostProb", &pi_TRACK_GhostProb, &b_pi_TRACK_GhostProb);
   fChain->SetBranchAddress("pi_TRACK_CloneDist", &pi_TRACK_CloneDist, &b_pi_TRACK_CloneDist);
   fChain->SetBranchAddress("pi_TRACK_Likelihood", &pi_TRACK_Likelihood, &b_pi_TRACK_Likelihood);
   fChain->SetBranchAddress("pi_ETA", &pi_ETA, &b_pi_ETA);
   fChain->SetBranchAddress("pi_PHI", &pi_PHI, &b_pi_PHI);
   fChain->SetBranchAddress("BachPi_CosTheta", &BachPi_CosTheta, &b_BachPi_CosTheta);
   fChain->SetBranchAddress("BachPi_OWNPV_X", &BachPi_OWNPV_X, &b_BachPi_OWNPV_X);
   fChain->SetBranchAddress("BachPi_OWNPV_Y", &BachPi_OWNPV_Y, &b_BachPi_OWNPV_Y);
   fChain->SetBranchAddress("BachPi_OWNPV_Z", &BachPi_OWNPV_Z, &b_BachPi_OWNPV_Z);
   fChain->SetBranchAddress("BachPi_OWNPV_XERR", &BachPi_OWNPV_XERR, &b_BachPi_OWNPV_XERR);
   fChain->SetBranchAddress("BachPi_OWNPV_YERR", &BachPi_OWNPV_YERR, &b_BachPi_OWNPV_YERR);
   fChain->SetBranchAddress("BachPi_OWNPV_ZERR", &BachPi_OWNPV_ZERR, &b_BachPi_OWNPV_ZERR);
   fChain->SetBranchAddress("BachPi_OWNPV_CHI2", &BachPi_OWNPV_CHI2, &b_BachPi_OWNPV_CHI2);
   fChain->SetBranchAddress("BachPi_OWNPV_NDOF", &BachPi_OWNPV_NDOF, &b_BachPi_OWNPV_NDOF);
   fChain->SetBranchAddress("BachPi_OWNPV_COV_", BachPi_OWNPV_COV_, &b_BachPi_OWNPV_COV_);
   fChain->SetBranchAddress("BachPi_IP_OWNPV", &BachPi_IP_OWNPV, &b_BachPi_IP_OWNPV);
   fChain->SetBranchAddress("BachPi_IPCHI2_OWNPV", &BachPi_IPCHI2_OWNPV, &b_BachPi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("BachPi_ORIVX_X", &BachPi_ORIVX_X, &b_BachPi_ORIVX_X);
   fChain->SetBranchAddress("BachPi_ORIVX_Y", &BachPi_ORIVX_Y, &b_BachPi_ORIVX_Y);
   fChain->SetBranchAddress("BachPi_ORIVX_Z", &BachPi_ORIVX_Z, &b_BachPi_ORIVX_Z);
   fChain->SetBranchAddress("BachPi_ORIVX_XERR", &BachPi_ORIVX_XERR, &b_BachPi_ORIVX_XERR);
   fChain->SetBranchAddress("BachPi_ORIVX_YERR", &BachPi_ORIVX_YERR, &b_BachPi_ORIVX_YERR);
   fChain->SetBranchAddress("BachPi_ORIVX_ZERR", &BachPi_ORIVX_ZERR, &b_BachPi_ORIVX_ZERR);
   fChain->SetBranchAddress("BachPi_ORIVX_CHI2", &BachPi_ORIVX_CHI2, &b_BachPi_ORIVX_CHI2);
   fChain->SetBranchAddress("BachPi_ORIVX_NDOF", &BachPi_ORIVX_NDOF, &b_BachPi_ORIVX_NDOF);
   fChain->SetBranchAddress("BachPi_ORIVX_COV_", BachPi_ORIVX_COV_, &b_BachPi_ORIVX_COV_);
   fChain->SetBranchAddress("BachPi_P", &BachPi_P, &b_BachPi_P);
   fChain->SetBranchAddress("BachPi_PT", &BachPi_PT, &b_BachPi_PT);
   fChain->SetBranchAddress("BachPi_PE", &BachPi_PE, &b_BachPi_PE);
   fChain->SetBranchAddress("BachPi_PX", &BachPi_PX, &b_BachPi_PX);
   fChain->SetBranchAddress("BachPi_PY", &BachPi_PY, &b_BachPi_PY);
   fChain->SetBranchAddress("BachPi_PZ", &BachPi_PZ, &b_BachPi_PZ);
   fChain->SetBranchAddress("BachPi_M", &BachPi_M, &b_BachPi_M);
   fChain->SetBranchAddress("BachPi_ID", &BachPi_ID, &b_BachPi_ID);
   fChain->SetBranchAddress("BachPi_PIDe", &BachPi_PIDe, &b_BachPi_PIDe);
   fChain->SetBranchAddress("BachPi_PIDmu", &BachPi_PIDmu, &b_BachPi_PIDmu);
   fChain->SetBranchAddress("BachPi_PIDK", &BachPi_PIDK, &b_BachPi_PIDK);
   fChain->SetBranchAddress("BachPi_PIDp", &BachPi_PIDp, &b_BachPi_PIDp);
   fChain->SetBranchAddress("BachPi_ProbNNe", &BachPi_ProbNNe, &b_BachPi_ProbNNe);
   fChain->SetBranchAddress("BachPi_ProbNNk", &BachPi_ProbNNk, &b_BachPi_ProbNNk);
   fChain->SetBranchAddress("BachPi_ProbNNp", &BachPi_ProbNNp, &b_BachPi_ProbNNp);
   fChain->SetBranchAddress("BachPi_ProbNNpi", &BachPi_ProbNNpi, &b_BachPi_ProbNNpi);
   fChain->SetBranchAddress("BachPi_ProbNNmu", &BachPi_ProbNNmu, &b_BachPi_ProbNNmu);
   fChain->SetBranchAddress("BachPi_ProbNNghost", &BachPi_ProbNNghost, &b_BachPi_ProbNNghost);
   fChain->SetBranchAddress("BachPi_hasMuon", &BachPi_hasMuon, &b_BachPi_hasMuon);
   fChain->SetBranchAddress("BachPi_isMuon", &BachPi_isMuon, &b_BachPi_isMuon);
   fChain->SetBranchAddress("BachPi_hasRich", &BachPi_hasRich, &b_BachPi_hasRich);
   fChain->SetBranchAddress("BachPi_UsedRichAerogel", &BachPi_UsedRichAerogel, &b_BachPi_UsedRichAerogel);
   fChain->SetBranchAddress("BachPi_UsedRich1Gas", &BachPi_UsedRich1Gas, &b_BachPi_UsedRich1Gas);
   fChain->SetBranchAddress("BachPi_UsedRich2Gas", &BachPi_UsedRich2Gas, &b_BachPi_UsedRich2Gas);
   fChain->SetBranchAddress("BachPi_RichAboveElThres", &BachPi_RichAboveElThres, &b_BachPi_RichAboveElThres);
   fChain->SetBranchAddress("BachPi_RichAboveMuThres", &BachPi_RichAboveMuThres, &b_BachPi_RichAboveMuThres);
   fChain->SetBranchAddress("BachPi_RichAbovePiThres", &BachPi_RichAbovePiThres, &b_BachPi_RichAbovePiThres);
   fChain->SetBranchAddress("BachPi_RichAboveKaThres", &BachPi_RichAboveKaThres, &b_BachPi_RichAboveKaThres);
   fChain->SetBranchAddress("BachPi_RichAbovePrThres", &BachPi_RichAbovePrThres, &b_BachPi_RichAbovePrThres);
   fChain->SetBranchAddress("BachPi_hasCalo", &BachPi_hasCalo, &b_BachPi_hasCalo);
   fChain->SetBranchAddress("BachPi_L0Global_Dec", &BachPi_L0Global_Dec, &b_BachPi_L0Global_Dec);
   fChain->SetBranchAddress("BachPi_L0Global_TIS", &BachPi_L0Global_TIS, &b_BachPi_L0Global_TIS);
   fChain->SetBranchAddress("BachPi_L0Global_TOS", &BachPi_L0Global_TOS, &b_BachPi_L0Global_TOS);
   fChain->SetBranchAddress("BachPi_Hlt1Global_Dec", &BachPi_Hlt1Global_Dec, &b_BachPi_Hlt1Global_Dec);
   fChain->SetBranchAddress("BachPi_Hlt1Global_TIS", &BachPi_Hlt1Global_TIS, &b_BachPi_Hlt1Global_TIS);
   fChain->SetBranchAddress("BachPi_Hlt1Global_TOS", &BachPi_Hlt1Global_TOS, &b_BachPi_Hlt1Global_TOS);
   fChain->SetBranchAddress("BachPi_Hlt1Phys_Dec", &BachPi_Hlt1Phys_Dec, &b_BachPi_Hlt1Phys_Dec);
   fChain->SetBranchAddress("BachPi_Hlt1Phys_TIS", &BachPi_Hlt1Phys_TIS, &b_BachPi_Hlt1Phys_TIS);
   fChain->SetBranchAddress("BachPi_Hlt1Phys_TOS", &BachPi_Hlt1Phys_TOS, &b_BachPi_Hlt1Phys_TOS);
   fChain->SetBranchAddress("BachPi_Hlt2Global_Dec", &BachPi_Hlt2Global_Dec, &b_BachPi_Hlt2Global_Dec);
   fChain->SetBranchAddress("BachPi_Hlt2Global_TIS", &BachPi_Hlt2Global_TIS, &b_BachPi_Hlt2Global_TIS);
   fChain->SetBranchAddress("BachPi_Hlt2Global_TOS", &BachPi_Hlt2Global_TOS, &b_BachPi_Hlt2Global_TOS);
   fChain->SetBranchAddress("BachPi_Hlt2Phys_Dec", &BachPi_Hlt2Phys_Dec, &b_BachPi_Hlt2Phys_Dec);
   fChain->SetBranchAddress("BachPi_Hlt2Phys_TIS", &BachPi_Hlt2Phys_TIS, &b_BachPi_Hlt2Phys_TIS);
   fChain->SetBranchAddress("BachPi_Hlt2Phys_TOS", &BachPi_Hlt2Phys_TOS, &b_BachPi_Hlt2Phys_TOS);
   fChain->SetBranchAddress("BachPi_TRACK_Type", &BachPi_TRACK_Type, &b_BachPi_TRACK_Type);
   fChain->SetBranchAddress("BachPi_TRACK_Key", &BachPi_TRACK_Key, &b_BachPi_TRACK_Key);
   fChain->SetBranchAddress("BachPi_TRACK_CHI2NDOF", &BachPi_TRACK_CHI2NDOF, &b_BachPi_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("BachPi_TRACK_PCHI2", &BachPi_TRACK_PCHI2, &b_BachPi_TRACK_PCHI2);
   fChain->SetBranchAddress("BachPi_TRACK_MatchCHI2", &BachPi_TRACK_MatchCHI2, &b_BachPi_TRACK_MatchCHI2);
   fChain->SetBranchAddress("BachPi_TRACK_GhostProb", &BachPi_TRACK_GhostProb, &b_BachPi_TRACK_GhostProb);
   fChain->SetBranchAddress("BachPi_TRACK_CloneDist", &BachPi_TRACK_CloneDist, &b_BachPi_TRACK_CloneDist);
   fChain->SetBranchAddress("BachPi_TRACK_Likelihood", &BachPi_TRACK_Likelihood, &b_BachPi_TRACK_Likelihood);
   fChain->SetBranchAddress("BachPi_ETA", &BachPi_ETA, &b_BachPi_ETA);
   fChain->SetBranchAddress("BachPi_PHI", &BachPi_PHI, &b_BachPi_PHI);
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

Bool_t Xib::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Xib::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Xib::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Xib_cxx
