#include <cmath>
#include <vector>
#include <map>
#include <iostream>
#include <TMatrixD.h>
#include "matrix.h"
//---------useful matrix functions defined by me-------------//
complex<double> Trace(matrix in){
  return in.trace();
}
matrix Transpose(matrix in){
  return in.transpose();
}
matrix Chop(matrix in){
  return in.chop();
}
matrix Inverse(matrix in){
  return in.inverse();
}
matrix Re(matrix in){
  return in.re();
}
double Re(complex<double> in){
  return in.real();
}
//----------end matrix functions----------------------------//
//----------declarations------------------------------------//
int SWITCH=3; //3 for Fig. 3 4 for Fig. 4
double PAR[] = {0.0401445, 1.87719, 0.530877, -0.251579, 5.65448, 
		0.132386, -1.56158, 5.95771, 0.697887, 1.26759, 2.87638, -1.65421, 
		1.13512, 0.298039, -3.30635, -1.53112, 0.270024, 
		0.27917, -1.37604, -1.54434};
double Mai[] = {0.49368, 0.49761, 0.13498, 0.13498, 0.13957, 0.13957, 0.54786, 
		0.54786, 0.49368, 0.49761};
double mai[] = {0.93827, 0.93956, 1.11568, 1.19264, 1.18937, 1.19745, 1.11568, 
		1.19264, 1.32171, 1.31486};
double Mai0x[] = {(Mai[0] + Mai[1])/2, (Mai[0] + Mai[1])/
		  2, (Mai[3] + Mai[4] + Mai[5])/
		  3, (Mai[3] + Mai[4] + Mai[5])/
		  3, (Mai[3] + Mai[4] + Mai[5])/
		  3, (Mai[3] + Mai[4] + Mai[5])/3, Mai[6], 
		  Mai[7], (Mai[0] + Mai[1])/2, (Mai[0] + Mai[1])/2};
double mai0x[] = {(mai[0] + mai[1])/2, (mai[0] + mai[1])/2, 
		  mai[2], (mai[3] + mai[4] + mai[5])/
		  3, (mai[3] + mai[4] + mai[5])/
		  3, (mai[3] + mai[4] + mai[5])/3, 
		  mai[6], (mai[3] + mai[4] + mai[5])/
		  3, (mai[9] + mai[8])/2, (mai[9] + mai[8])/2};
const int Channels = 10;
matrix mbn(Channels,Channels);//mbn = Chop[DiagonalMatrix[Table[mai[[i]], {i, 1, Channels}]]];
matrix mmn(Channels,Channels);//mmn = Chop[DiagonalMatrix[Table[Mai[[i]], {i, 1, Channels}]]];
double fmespi = 0.09240, fmesk = 0.1130, fmeseta = 1.30*fmespi;
matrix fmes(Channels,Channels);//fmes = Chop[DiagonalMatrix[Table[{fmesk, fmesk, fmespi, fmespi, fmespi, fmespi, fmeseta,fmeseta, fmesk, fmesk}[[i]], {i, 1, Channels}]]];
matrix eins(Channels,Channels);//eins = IdentityMatrix[Channels];
matrix null(Channels,Channels);//null = 0*eins;
matrix im(Channels,Channels);//im = null;
matrix ib(Channels,Channels);//ib = null;
  
double mu = (0.5)*(Mai[0]*Mai[0] - Mai[1]*Mai[1] + Mai[2]*Mai[2]);
double md = (0.5)*(Mai[1]*Mai[1] - Mai[0]*Mai[0] + Mai[2]*Mai[2]);
double ms = (0.5)*(Mai[0]*Mai[0] + Mai[1]*Mai[1] - Mai[2]*Mai[2]);
matrix Mquark(3,3),l1(3,3),l2(3,3),l3(3,3),l4(3,3),l5(3,3),l6(3,3),l7(3,3),l8(3,3);
double Pi=3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745;
double wd4 = 0/(16*Pi*Pi);

map<string,matrix > fm;
map<string,matrix > fb;

string mlist[] = {"Kminus","Kbar0","Pi0","Pi0","Piminus","Piplus","Eta","Eta","Kplus","K0"};
string blist[] = {"proton","neutron","lambda","sigma0","sigmaplus","sigmaminus","lambda","sigma0","cascademinus","cascadenull"};
complex<double> spurWT(matrix lb,matrix lj,matrix li,matrix la){
  return Trace(Transpose(lb)*((Transpose(lj)*li - li*Transpose(lj))*la - la*(Transpose(lj)*li - li*Transpose(lj))));
}
complex<double> spurg1(matrix lb,matrix lj,matrix li,matrix la){
  return Trace(Transpose(lb)*(Transpose(lj)*(li*la - la*li) - (li*la - la*li)*Transpose(lj))) + Trace(Transpose(lb)*(li*(Transpose(lj)*la - la*Transpose(lj)) - (Transpose(lj)*la - la*Transpose(lj))*li));
}
complex<double> spurg2(matrix lb,matrix lj,matrix li,matrix la){
  return Trace(Transpose(lb)*(Transpose(lj)*(li*la + la*li) - (li*la + la*li)*Transpose(lj))) + Trace(Transpose(lb)*(li*(Transpose(lj)*la + la*Transpose(lj)) - (Transpose(lj)*la + la*Transpose(lj))*li)); 
}
complex<double> spurg3(matrix lb,matrix lj,matrix li,matrix la){
  return Trace(Transpose(lb)*(Transpose(lj)*(li*la + la*li) + (li*la + la*li)*Transpose(lj))) + Trace(Transpose(lb)*(li*(Transpose(lj)*la + la*Transpose(lj)) + (Transpose(lj)*la + la*Transpose(lj))*li)); 
}
complex<double> spurg4(matrix lb,matrix lj,matrix li,matrix la){
  return 2*Trace(Transpose(lb)*la)*Trace(Transpose(lj)*li); 
}
complex<double> spurg5(matrix lb,matrix lj,matrix li,matrix la){
  return Trace(Transpose(lb)*((Transpose(lj)*li - li*Transpose(lj))*la - la*(Transpose(lj)*li - li*Transpose(lj)))); 
}
complex<double> spurg6(matrix lb,matrix lj,matrix li,matrix la){
  return Trace(Transpose(lb)*((Transpose(lj)*li - li*Transpose(lj))*la + la*(Transpose(lj)*li - li*Transpose(lj)))); 
}
complex<double> spurg7(matrix lb,matrix lj,matrix li,matrix la){
  return (0.5)*(Trace(Transpose(lb)*Transpose(lj))*Trace(li*la) - Trace(Transpose(lb)*li)*Trace(Transpose(lj)*la)); 
}
complex<double> spurgnull(matrix lb,matrix lj,matrix li,matrix la,matrix Mquark){
  return Trace(Transpose(lb)*la)*Trace((Transpose(lj)*li + li*Transpose(lj))*Mquark); 
}
matrix helpgDF(matrix lb,matrix lj,matrix li,matrix la,matrix Mquark){
  matrix out1 = (Transpose(lj)*li + li*Transpose(lj))*Mquark;
  matrix out2 = Mquark*(Transpose(lj)*li + li*Transpose(lj));
  matrix out3 = Transpose(lj)*Mquark*li;
  out3*=2;
  matrix out4 = li*Mquark*Transpose(lj);
  out4*=2;
  return out1 + out2 + out3 + out4; 
}
complex<double> spurgD(matrix lb,matrix lj,matrix li,matrix la,matrix Mquark){
  return Trace(Transpose(lb)*(helpgDF(lb, lj, li, la, Mquark)*la + la*helpgDF(lb, lj, li, la, Mquark))); 
}
complex<double> spurgF(matrix lb,matrix lj,matrix li,matrix la,matrix Mquark){
  return Trace(Transpose(lb)*(helpgDF(lb, lj, li, la, Mquark)*la - la*helpgDF(lb, lj, li, la, Mquark))); 
}
matrix spurWTmatrix(Channels,Channels);
matrix spurn(Channels,Channels);
matrix g(Channels,Channels);//g = -(1/4)*Chop[((Inverse[fmes]).spurn.(Inverse[fmes]))];
double fpi = 0.1;
matrix spurg1matrix(Channels,Channels); 
matrix spurg2matrix(Channels,Channels); 
matrix spurg3matrix(Channels,Channels); 
matrix spurg4matrix(Channels,Channels); 
matrix spurg5matrix(Channels,Channels);
matrix spurg6matrix(Channels,Channels);
matrix spurg7matrix(Channels,Channels);
matrix spurgnullmatrix(Channels,Channels);
matrix spurgDmatrix(Channels,Channels);
matrix spurgFmatrix(Channels,Channels);
matrix spurg8matrix(Channels,Channels);
matrix spurg9matrix(Channels,Channels);
matrix spurg10matrix(Channels,Channels);
matrix spurg11matrix(Channels,Channels);
int m0=1;
//----------------------end declarations---------------------//
void fillthings(){
  if(SWITCH == 4){std::copy(Mai0x,Mai0x + 10,Mai); std::copy(mai0x,mai0x+10,mai);}//decide which set to use

  for(int i=0; i<Channels; i++) for(int j=0; j<Channels; j++){
      null[i][j]=0; im[i][j]=0; ib[i][j]=0;
      if(i==j){
	eins[i][j]=1;
	mbn[i][j]=mai[i];
	mmn[i][j]=Mai[i];
	if((i>=0 && i<2)||(i>=8 && i<10)){ fmes[i][j]=fmesk;
	}else if(i>=2 && i<6){ fmes[i][j]=fmespi;
	}else if(i>=6 && i<8){ fmes[i][j]=fmeseta;}
      }else{
	mbn[i][j]=0; mmn[i][j]=0; fmes[i][j]=0; eins[i][j]=0;
      }
    }
  mbn=Chop(mbn);
  mmn=Chop(mmn);
  for(int i=0; i<Mquark.GetNrows(); i++) for(int j=0; j<Mquark.GetNcols(); j++){
      if(i==j){
	if(i==0) Mquark[i][j] = mu;
	if(i==1) Mquark[i][j] = md;
	if(i==2) Mquark[i][j] = ms;
      }else Mquark[i][j]=0;
    }
  l1[0][0]=(1/sqrt(2)); l1[0][1]=0;            l1[0][2]=0;
  l1[1][0]=0;           l1[1][1]=(-1/sqrt(2)); l1[1][2]=0;
  l1[2][0]=0;           l1[2][1]=0;            l1[2][2]=0;

  l2[0][0]=(1/sqrt(6)); l2[0][1]=0;            l2[0][2]=0;
  l2[1][0]=0;           l2[1][1]=(1/sqrt(6));  l2[1][2]=0;
  l2[2][0]=0;           l2[2][1]=0;            l2[2][2]=(-2/sqrt(6));

  l3[0][0]=0;           l3[0][1]=1;            l3[0][2]=0;
  l3[1][0]=0;           l3[1][1]=0;            l3[1][2]=0;
  l3[2][0]=0;           l3[2][1]=0;            l3[2][2]=0;

  l4[0][0]=0;           l4[0][1]=0;            l4[0][2]=0;
  l4[1][0]=1;           l4[1][1]=0;            l4[1][2]=0;
  l4[2][0]=0;           l4[2][1]=0;            l4[2][2]=0;

  l5[0][0]=0;           l5[0][1]=0;            l5[0][2]=1;
  l5[1][0]=0;           l5[1][1]=0;            l5[1][2]=0;
  l5[2][0]=0;           l5[2][1]=0;            l5[2][2]=0;

  l6[0][0]=0;           l6[0][1]=0;            l6[0][2]=0;
  l6[1][0]=0;           l6[1][1]=0;            l6[1][2]=1;
  l6[2][0]=0;           l6[2][1]=0;            l6[2][2]=0;

  l7[0][0]=0;           l7[0][1]=0;            l7[0][2]=0;
  l7[1][0]=0;           l7[1][1]=0;            l7[1][2]=0;
  l7[2][0]=0;           l7[2][1]=1;            l7[2][2]=0;

  l8[0][0]=0;           l8[0][1]=0;            l8[0][2]=0;
  l8[1][0]=0;           l8[1][1]=0;            l8[1][2]=0;
  l8[2][0]=1;           l8[2][1]=0;            l8[2][2]=0;
  fm = {{"Pi0",l1},{"Pi0",l1},{"Piplus",l3},{"Eta",l2},{"Kplus",l5},{"K0",l6},{"Kbar0",l7},{"Kminus",l8},{"Piminus",l4}};
  fb = {{"proton",l5},{"neutron",l6},{"lambda",l2},{"sigma0",l1},{"sigmaplus",l3},{"sigmaminus",l4},{"cascademinus",l8},{"cascadenull",l7}};

  for(int i=0; i<Channels; i++)
    for(int j=0; j<Channels; j++)
      spurWTmatrix[i][j] = spurWT(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
  spurn = spurWTmatrix;
  
  g = fmes.inverse()*spurn*fmes.inverse();
  g *= (-0.25);
  g = Chop(g);
  for(int i=0; i<Channels; i++)
    for(int j=0; j<Channels; j++){
      spurg1matrix[i][j] = spurg1(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurg2matrix[i][j] = spurg2(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurg3matrix[i][j] = spurg3(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurg4matrix[i][j] = spurg4(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurg5matrix[i][j] = spurg5(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurg6matrix[i][j] = spurg6(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurg7matrix[i][j] = spurg7(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]]);
      spurgnullmatrix[i][j] = spurgnull(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]],Mquark);
      spurgDmatrix[i][j] = spurgD(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]],Mquark);
      spurgFmatrix[i][j] = spurgF(fb[blist[i]],fm[mlist[i]],fm[mlist[j]],fb[blist[j]],Mquark);
    }
  spurg8matrix  = spurg1matrix;
  spurg9matrix  = spurg2matrix;
  spurg10matrix = spurg3matrix;
  spurg11matrix = spurg4matrix;
}
matrix Gvec; //set by BONN()
matrix t(Channels,Channels); //set by BONN()
void BONN(complex<double> svar, complex<double> logKN, complex<double> logPiL, complex<double> logPiS, complex<double> logEtaL, complex<double> logEtaS, complex<double> logKXi, complex<double> b1var, complex<double> b2var, complex<double> b3var, complex<double> b4var, complex<double> b5var, complex<double> b6var, complex<double> b7var, complex<double> b8var, complex<double> b9var, complex<double> b10var, complex<double> b11var, complex<double> bnullvar, complex<double> bDvar, complex<double> bFvar){
  fillthings();
  complex<double> s = svar, bnull = bnullvar, bD = bDvar, bF = bFvar, b1 = -b1var, b2 = -b2var, b3 = -b3var, b4 = -b4var, b5 = b5var, b6 = b6var, b7 = b7var, b8 = -b8var, b9 = -b9var, b10 = -b10var, b11 = -b11var;
  complex<double> lg[] = {logKN, logKN, logPiL, logPiS, logPiS, logPiS, logEtaL, logEtaS, logKXi, logKXi};
  complex<double> m5 = mai[0], mp = m5, d = 4, p = sqrt(s), W = sqrt(s), s0 = m5*m5;

  matrix qcms(Channels,Channels);
  for(int i=0; i<Channels; i++) 
    for(int j=0; j<Channels; j++)
      if(i==j)
	qcms[i][j] = (sqrt((pow(svar,2) + pow(mai[i],4) + pow(Mai[i],4)) + ( - 2*pow(mai[i],2)*pow(Mai[i],2)-2*svar*pow(mai[i],2) - 2*svar*pow(Mai[i],2))))/(2*sqrt(svar));//extra parentheses in numerator because of weirdness with order of operations--don't understand why it goes wrong
      else qcms[i][j]=0;
  //qcms = DiagonalMatrix[Table[(Sqrt[ svar^2 + mai[[i]]^4 + Mai[[i]]^4 - 2*mai[[i]]^2*Mai[[i]]^2 -2*svar*mai[[i]]^2 - 2*svar*Mai[[i]]^2])/(2*Sqrt[svar]), {i, 1,Channels}]];
  matrix imb(Channels,Channels);
  for(int i=0; i<Channels; i++) 
    for(int j=0; j<Channels; j++)
      if(i==j)
	imb[i][j] = (1/(16*Pi*Pi))*(-1 + 2*log(mai[i]) - 2*lg[i] + ((pow(Mai[i],2) - pow(mai[i],2) + svar)/(2*svar))*log(pow((Mai[i]/mai[i]),2)) - ((4*qcms[i][i])/sqrt(svar))*ArcTanh((2*qcms[i][i]*sqrt(svar))/(pow((mai[i] + Mai[i]),2) - svar)));
      else imb[i][j]=0;
  //imb = DiagonalMatrix[ Table[(1/(16*Pi^2))*(-1 + 2 Log[mai[[i]]] - 2 lg[[i]] + ((Mai[[i]]^2 - mai[[i]]^2 + svar)/(2*svar))*Log[(Mai[[i]]/mai[[i]])^2] - ((4*qcms[[i, i]])/Sqrt[svar])*ArcTanh[(2*qcms[[i, i]]*Sqrt[svar])/((mai[[i]] + Mai[[i]])^2 - svar)]), {i, 1, Channels}]];
  matrix reqcms = qcms.re();
  
  matrix a = (1/(4*p*p*(d - 1)))*((4*p*p*mmn*mmn - pow(p*p*eins + mmn*mmn - mbn*mbn,2))*imb + (p*p*eins +  mmn*mmn - mbn*mbn)*im + (p*p*eins - mmn*mmn + mbn*mbn)*ib) - (wd4/(18*p*p))*(4*p*p*mmn*mmn - pow(p*p*eins + mmn*mmn - mbn*mbn,2) +  mmn*mmn*(p*p*eins + mmn*mmn - mbn*mbn) +  mbn*mbn*(p*p*eins - mmn*mmn + mbn*mbn));
  matrix b = (1/(4* p*p*(d - 1)))*((d*pow(p*p*eins + mmn*mmn - mbn*mbn,2) -  4*p*p*mmn*mmn)*imb - d*(p*p*eins + mmn*mmn -  mbn*mbn)*im + (d*(3*p*p*eins + mmn*mmn - mbn*mbn) - 4*p*p*eins)*ib) - (wd4/(18*p*p))*(pow(p*p*eins + mmn*mmn - mbn*mbn,2) - 4*p*p*mmn*mmn -  mmn*mmn* (p*p*eins + mmn*mmn - mbn*mbn) +  mbn*mbn* (mmn*mmn - mbn*mbn - p*p*eins));
  matrix G0 = Chop(mbn*imb);
  matrix G1 = Chop((1/(2*p*p))*((p*p*eins - mmn*mmn + mbn*mbn)*imb + im - ib));
  matrix A0 = Chop((1/(4*p*p*(d - 1)))*((pow(p*p*eins - mmn*mmn + mbn*mbn,2) - 4*p*p*mbn*mbn)*imb - (p*p*eins + mmn*mmn - mbn*mbn)*im - (p*p*eins - mmn*mmn + mbn*mbn)*ib) - (wd4/(18* p*p))*(pow(p*p*eins - mmn*mmn + mbn*mbn,2) - 4*p*p*mbn*mbn - mmn*mmn* (p*p*eins + mmn*mmn - mbn*mbn) - mbn*mbn* (p*p*eins - mmn*mmn + mbn*mbn)));
  matrix A = A0;
  matrix B0 = Chop((1/(2*p*p))*mbn*((p*p*eins + mmn*mmn - mbn*mbn)*imb + ib - im));
  matrix B1 =  Chop((1/(4*p*p*p*p*(d -  1)))*(((d*(p*p*eins - mbn*mbn + mmn*mmn) -  2*p*p*eins)*(p*p*eins + mbn*mbn - mmn*mmn) + 4*p*p*mbn*mbn)*imb - (d*(mbn*mbn - mmn*mmn) + (d - 2)*p*p*eins)*im + (d*(mbn*mbn - mmn*mmn) - (d - 2)*p*p*eins)*ib) - (wd4/(18*p*p*p*p))*((mmn*mmn - mbn*mbn - p*p*eins)*(p*p*eins - mmn*mmn +  mbn*mbn) + 4*p*p*mbn*mbn - mmn*mmn* (mbn*mbn - mmn*mmn) +  mbn*mbn* (mbn*mbn - mmn*mmn) +  p*p*(mbn*mbn + mmn*mmn)));
  matrix Ha = (1/(2* p*p))*((p*p*eins + mmn*mmn - mbn*mbn)*a + ((1/d)*mbn*mbn)*ib - ((1/d)*mmn*mmn)*im) - (wd4/(16*p*p))*(pow(mbn,4) - pow(mmn,4));
  matrix Hb = (1/(2* p*p))*((p*p*eins + mmn*mmn - mbn*mbn)*(b - 2*a) + (p*p*eins - ((1/d)*(2*mbn*mbn)))*ib + ((1/d)*(2* mmn*mmn))*im) - (wd4/(8*p*p))*(pow(mmn,4) - pow(mbn,4));
  matrix C0 = Chop(mbn*a);
  matrix D0 = Chop(mbn*b);
  matrix C1 = Chop(a - Ha);
  matrix D1 = Chop(b - Hb);
  matrix E0 = Chop(-1*Ha);
  matrix h0s = Chop(p*p*G1 - mbn*G0 - im);
  matrix h1s = Chop(G0 - mbn*G1);
  
  matrix G1matrix = -1*Chop( Inverse(fmes)*(-2*b1*spurg1matrix - 2*b2*spurg2matrix - 2*b3*spurg3matrix - 2*b4*spurg4matrix)*Inverse(fmes));
  matrix G2matrix = Chop(Inverse( fmes)*(2*b5*spurg5matrix + 2*b6*spurg6matrix + 2*b7*spurg7matrix)*Inverse(fmes));
  matrix G3matrix = -1*Chop((1/m0)*Inverse(fmes)*(-b8*spurg8matrix - b9*spurg9matrix -b10*spurg10matrix - b11*spurg11matrix)*Inverse(fmes));
  matrix G0matrix = Chop(Inverse( fmes)*(4*bnull*spurgnullmatrix + bD*spurgDmatrix + bF*spurgFmatrix)*Inverse(fmes));
  
  matrix pmalqm = Chop(0.5*(s*eins - mbn*mbn + mmn*mmn));
  matrix qmalpm = pmalqm;
  matrix V1p = 2*g - 2*(mbn*G2matrix + G2matrix*mbn) + G3matrix*pmalqm + qmalpm*G3matrix;
  matrix V11 = -mbn*g - g*mbn + G0matrix + 2*(s*G2matrix + mbn*G2matrix*mbn) - mbn*G3matrix*pmalqm - qmalpm*G3matrix*mbn;
  matrix V21 = G1matrix - 2*G2matrix;
  matrix vpt2 = -C1*V21;
  matrix v1t2 = eins - C0*V21;
  matrix T2on1 = V21*Inverse(v1t2 - s*vpt2*Inverse(v1t2)*vpt2);
  matrix T2onp = -V21*Inverse(v1t2)*vpt2*Inverse(v1t2 - s*vpt2*Inverse(v1t2)*vpt2);
  
  matrix vpt1 = -G0*V1p - G1*V11 - A0*V21 - B1*V21*pmalqm;
  matrix v1t1 = eins - G0*V11 - s*G1*V1p + A0*V21*mbn - B0*V21*pmalqm;
  matrix inv1t1 = Inverse(v1t1 - s*vpt1*Inverse(v1t1)*vpt1);
  matrix invpt1 = -Inverse(v1t1)*vpt1*inv1t1;
  matrix pot1t1 = V11 + s*(T2on1*A*V1p) + (2*qmalpm - s*eins)*T2onp*A*V11 - mbn*(T2on1*A*V11) + s*mbn*T2onp*A*V1p + qmalpm*((T2on1*B0 + s*T2onp*B1)*V11 + s*(T2on1*B1 + T2onp*B0)*V1p) + (1/s)*qmalpm*(T2on1*D0 + s*T2onp*D1)*V21*pmalqm + ((2*qmalpm - s*eins)*T2onp - mbn*T2on1)*E0*V21*pmalqm + qmalpm*(s*T2onp*E0*V21 - T2on1*E0*V21*mbn);
  matrix potpt1 = V1p + T2on1*A*V11 + (2*qmalpm - s*eins)*T2onp*A*V1p - mbn*(T2on1*A*V1p - T2onp*A*V11) + qmalpm*((T2on1*B1 + T2onp*B0)*V11 + (T2on1*B0 +s*T2onp*B1)*V1p) + (1/s)*qmalpm*(T2on1*D1 + T2onp*D0)*V21*pmalqm + (T2on1 +mbn*T2onp)*E0*V21*pmalqm + qmalpm*(T2on1*E0*V21 - T2onp*E0*V21*mbn);
  matrix T1on1 = pot1t1*inv1t1 + s*potpt1*invpt1;
  matrix T1onp = pot1t1*invpt1 + potpt1*inv1t1;
  
  matrix qnullcms = sqrt(qcms*qcms + mmn*mmn);
  // T0ON = T1on1 + qnullcms*T2on1*qnullcms - z*qcms*T2on1*qcms;
  // T1ON = T1onp + qnullcms*T2onp*qnullcms - z*qcms*T2onp*qcms;
  matrix T0ON_0 = T1on1 + qnullcms*T2on1*qnullcms - 0*qcms*T2on1*qcms;//T0ON with z->0
  matrix T1ON_0 = T1onp + qnullcms*T2onp*qnullcms - 0*qcms*T2onp*qcms;//T1On with z->0

  matrix Ecms = sqrt(mbn*mbn + qcms*qcms);
  matrix fnullpluson0 = Chop(-(1/(16*Pi*sqrt(s)))*(sqrt(Ecms + mbn)*(2*(T0ON_0 + p*T1ON_0) )*sqrt(Ecms + mbn) +sqrt(Ecms - mbn)*(((double)2/(double)3)*(qcms*T2on1*qcms - p*qcms*T2onp*qcms)/* Coefficient(-T0ON + p*T1ON, z)*/)*sqrt(Ecms - mbn)));
  matrix maimatrix(sizeof(mai)/sizeof(mai[0]),1);
  for(unsigned int i=0; i<(sizeof(mai)/sizeof(mai[0])); i++)
    maimatrix[i][0] = mai[i];
  Gvec = 2*imb*maimatrix;
  //t = Table[fnullpluson0[[i, j]]*(-4 \[Pi] Sqrt[s]/Sqrt[mai[[i]]]/Sqrt[mai[[j]]]), {i, 1, Channels}, {j, 1, Channels}];
  for(int i=0; i<Channels; i++)
    for(int j=0; j<Channels; j++)
      t[i][j] = fnullpluson0[i][j]*(-4*Pi*sqrt(s)/sqrt(mai[i])/sqrt(mai[j]));
}
void invariantmassdistribution_original(){
  double MLb = 5.619;
  double MJpsi = 3.100;
  int NN = 501;
  matrix RES(501,7,0);
  double NORMALIZATION = 10^5;
  double h1 = 1, h4 = 0, h7 = -sqrt(2)/3, h10 = 0;
  double maix[10]; double Maix[10];
  std::copy(mai,mai + 10,maix); std::copy(Mai,Mai + 10,Maix);

  for(int i=0; i<NN; i++){
    complex<double> Wx = complex<double>(1.34 + 0.26*(i-1)/NN, 0); //1.34 + 0.26 (i - 1)/NN;
    // Call FSI from Bonn model
    BONN(Wx*Wx, PAR[0], PAR[1], PAR[2], PAR[3], PAR[4], PAR[5], PAR[6], PAR[7], PAR[8], PAR[9], PAR[10], PAR[11], PAR[12], PAR[13], PAR[14], PAR[15], PAR[16], PAR[17], PAR[18], PAR[19]);
    // Eq. 2,3
    RES[i] =
      {Wx,
       NORMALIZATION*1/pow(2*Pi,3)*Re(sqrt((MLb*MLb - pow(MJpsi + Wx,2))*(MLb*MLb - pow(MJpsi - Wx,2)))*sqrt((Wx*Wx - pow(Maix[0] + maix[0],2))*(Wx*Wx - pow(Maix[0] -maix[0],2))))/(16*MLb*MLb*MLb*Wx)
       *pow(abs((h1 + h1*Gvec[0][0]*t[0][0] + h1*Gvec[1][0]*t[1][0] + h4*Gvec[3][0]*t[3][0] + h4*Gvec[4][0]*t[4][0] + h4*Gvec[5][0]*t[5][0] + h7*Gvec[6][0]*t[6][0] + h10*Gvec[8][0]*t[8][0] + h10*Gvec[9][0]*t[9][0])),2),
       NORMALIZATION*1/pow(2*Pi,3)*Re(sqrt((MLb*MLb - pow(MJpsi + Wx,2))*(MLb*MLb - pow(MJpsi - Wx,2)))*sqrt((Wx*Wx - pow(Maix[1] + maix[1],2))*(Wx*Wx - pow(Maix[1] -maix[1],2))))/(16*MLb*MLb*MLb*Wx)
       *pow(abs((h1 + h1*Gvec[0][0]*t[0][1] + h1*Gvec[1][0]*t[1][1] +h4*Gvec[3][0]*t[3][1] + h4*Gvec[4][0]*t[4][1] +h4*Gvec[5][0]*t[5][1] + h7*Gvec[6][0]*t[6][1] +h10*Gvec[8][0]*t[8][1] + h10*Gvec[9][0]*t[9][1])),2),
       0,
       NORMALIZATION*1/pow(2*Pi,3)*Re( sqrt((MLb*MLb - pow(MJpsi + Wx,2))*(MLb*MLb - pow(MJpsi -Wx,2)))*sqrt((Wx*Wx - pow(Maix[3] +maix[3],2))*(Wx*Wx - pow(Maix[3] - maix[3],2))))/(16*MLb*MLb*MLb*Wx)
       *pow(abs((h4 +h4*Gvec[3][0]*t[3][3] + h4*Gvec[4][0]*t[4][3] +h4*Gvec[5][0]*t[5][3] + h7*Gvec[6][0]*t[6][3] +h10*Gvec[8][0]*t[8][3] + h10*Gvec[9][0]*t[9][3] +h1*Gvec[0][0]*t[0][3] + h1*Gvec[1][0]*t[1][3])),2),
       NORMALIZATION*1/pow(2*Pi,3)*Re( sqrt((MLb*MLb - pow(MJpsi + Wx,2))*(MLb*MLb - pow(MJpsi -Wx,2)))*sqrt((Wx*Wx - pow(Maix[4] +maix[4],2))*(Wx*Wx - pow(Maix[4] - maix[4],2))))/(16*MLb*MLb*MLb*Wx)
       *pow(abs((h4 +h4*Gvec[3][0]*t[3][4] + h4*Gvec[4][0]*t[4][4] +h4*Gvec[5][0]*t[5][4] + h7*Gvec[6][0]*t[6][4] +h10*Gvec[8][0]*t[8][4] + h10*Gvec[9][0]*t[9][4] +h1*Gvec[0][0]*t[0][4] + h1*Gvec[1][0]*t[1][4])),2),
       NORMALIZATION*1/pow(2*Pi,3)*Re( sqrt((MLb*MLb - pow(MJpsi + Wx,2))*(MLb*MLb - pow(MJpsi -Wx,2)))*sqrt((Wx*Wx - pow(Maix[5] +maix[5],2))*(Wx*Wx - pow(Maix[5] - maix[5],2))))/(16*MLb*MLb*MLb*Wx)
       *pow(abs((h4 +h4*Gvec[3][0]*t[3][5] + h4*Gvec[4][0]*t[4][5] +h4*Gvec[5][0]*t[5][5] + h7*Gvec[6][0]*t[6][5] +h10*Gvec[8][0]*t[8][5] + h10*Gvec[9][0]*t[9][5] +h1*Gvec[0][0]*t[0][5] + h1*Gvec[1][0]*t[1][5])),2)};
  }
}
double invariantmassdistribution_sig0pi0(double Wx){//input mass in GeV
  double MLb = 5.619;
  double MJpsi = 3.100;
  double NORMALIZATION = pow(10,5);
  double h1 = 1, h4 = 0, h7 = -sqrt(2)/3, h10 = 0;
  double maix[10]; double Maix[10];
  std::copy(mai,mai + 10,maix); std::copy(Mai,Mai + 10,Maix);

  // Call FSI from Bonn model
  BONN(Wx*Wx, PAR[0], PAR[1], PAR[2], PAR[3], PAR[4], PAR[5], PAR[6], PAR[7], PAR[8], PAR[9], PAR[10], PAR[11], PAR[12], PAR[13], PAR[14], PAR[15], PAR[16], PAR[17], PAR[18], PAR[19]);
  return NORMALIZATION*1/pow(2*Pi,3)*Re( sqrt((MLb*MLb - pow(MJpsi + Wx,2))*(MLb*MLb - pow(MJpsi -Wx,2)))*sqrt((Wx*Wx - pow(Maix[3] +maix[3],2))*(Wx*Wx - pow(Maix[3] - maix[3],2))))/(16*MLb*MLb*MLb*Wx)
    *pow(abs((h4 +h4*Gvec[3][0]*t[3][3] + h4*Gvec[4][0]*t[4][3] +h4*Gvec[5][0]*t[5][3] + h7*Gvec[6][0]*t[6][3] +h10*Gvec[8][0]*t[8][3] + h10*Gvec[9][0]*t[9][3] +h1*Gvec[0][0]*t[0][3] + h1*Gvec[1][0]*t[1][3])),2);//value of Bonn model in arb. units.

}
