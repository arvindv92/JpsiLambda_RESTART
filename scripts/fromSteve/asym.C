void asym(){

  double ylb_Run1[2] = {6827, 6480};
  double ylb_Run2[2] = {7210, 6822};
  double eylb_Run1[2] = {94, 92};
  double eylb_Run2[2] = {99, 97};

  double yxb_Run1[2] = {236, 230};
  double yxb_Run2[2] = {304, 326};
  double eyxb_Run1[2] = {18, 18};
  double eyxb_Run2[2] = {21, 22};
  

  double asym_Lb_Run1 = (ylb_Run1[0] - ylb_Run1[1]) / (ylb_Run1[0] + ylb_Run1[1]) ;
  double asym_Lb_Run2 = (ylb_Run2[0] - ylb_Run2[1]) / (ylb_Run2[0] + ylb_Run2[1]) ;
  double asym_Xb_Run1 = (yxb_Run1[0] - yxb_Run1[1]) / (yxb_Run1[0] + yxb_Run1[1]) ;
  double asym_Xb_Run2 = (yxb_Run2[0] - yxb_Run2[1]) / (yxb_Run2[0] + yxb_Run2[1]) ;
  
  double err_asym_Lb_Run1 = ((2*ylb_Run1[0]*ylb_Run1[1])/(ylb_Run1[0] + ylb_Run1[1])**2);
  err_asym_Lb_Run1 = err_asym_Lb_Run1*sqrt((eylb_Run1[0]/ylb_Run1[0])**2+(eylb_Run1[1]/ylb_Run1[1])**2);

  double err_asym_Lb_Run2 = ((2*ylb_Run2[0]*ylb_Run2[1])/(ylb_Run2[0] + ylb_Run2[1])**2);
  err_asym_Lb_Run2 = err_asym_Lb_Run2*sqrt((eylb_Run2[0]/ylb_Run2[0])**2+(eylb_Run2[1]/ylb_Run2[1])**2);


  double err_asym_Xb_Run1 = ((2*yxb_Run1[0]*yxb_Run1[1])/(yxb_Run1[0] + yxb_Run1[1])**2);
  err_asym_Xb_Run1 = err_asym_Xb_Run1*sqrt((eyxb_Run1[0]/yxb_Run1[0])**2+(eyxb_Run1[1]/yxb_Run1[1])**2);

  double err_asym_Xb_Run2 = ((2*yxb_Run2[0]*yxb_Run2[1])/(yxb_Run2[0] + yxb_Run2[1])**2);
  err_asym_Xb_Run2 = err_asym_Xb_Run2*sqrt((eyxb_Run2[0]/yxb_Run2[0])**2+(eyxb_Run2[1]/yxb_Run2[1])**2);
  

  double Lb_Ave = (asym_Lb_Run1/err_asym_Lb_Run1**2 + asym_Lb_Run2/err_asym_Lb_Run2**2) / (1/err_asym_Lb_Run1**2 + 1/err_asym_Lb_Run2**2);
  double Xb_Ave = (asym_Xb_Run1/err_asym_Xb_Run1**2 + asym_Xb_Run2/err_asym_Xb_Run2**2) / (1/err_asym_Xb_Run1**2 + 1/err_asym_Xb_Run2**2);
  double Lb_errAve = sqrt(1./((1/err_asym_Lb_Run1**2 + 1/err_asym_Lb_Run2**2)));
  double Xb_errAve = sqrt(1./((1/err_asym_Xb_Run1**2 + 1/err_asym_Xb_Run2**2)));

  cout << "Lambda_b" << endl;
  cout << "Run 1 Asym = " << asym_Lb_Run1 << "+-" << err_asym_Lb_Run1 << endl;
  cout << "Run 2 Asym = " << asym_Lb_Run2 << "+-" << err_asym_Lb_Run2 << endl;
  cout << "Ave Asym = " << Lb_Ave << "+-" << Lb_errAve << endl;

  cout << "Xi_b" << endl;
  cout << "Run 1 Asym = " << asym_Xb_Run1 << "+-" << err_asym_Xb_Run1 << endl;
  cout << "Run 2 Asym = " << asym_Xb_Run2 << "+-" << err_asym_Xb_Run2 << endl;
  cout << "Ave Asym = " << Xb_Ave << "+-" << Xb_errAve << endl;  

  double d1 = asym_Xb_Run1 - asym_Lb_Run1;
  double ed1 = sqrt(err_asym_Xb_Run1**2 + err_asym_Lb_Run1**2);
  double d2 = asym_Xb_Run2 - asym_Lb_Run2;
  double ed2 = sqrt(err_asym_Xb_Run2**2 + err_asym_Lb_Run2**2);
  cout << "A(Xib)-A(Lb)" << endl;
  cout << "Run 1 delta Asym = " << d1 << "+-" << ed1 << endl;
  cout << "Run 2 delta Asym = " << d2 << "+-" << ed2 << endl;
  
  double d1c = d1 + 0.024;
  double d2c = d2 + 0.024;
  
  cout << "Absolute asymmetries" << endl;
  cout << "Run 1 Xib- Asym = " << d1c << "+-" << ed1 << endl;
  cout << "Run 2 Xib- Asym = " << d2c << "+-" << ed2 << endl;
  
  //return;
  

  double dXib_Asym = Xb_Ave - Lb_Ave;
  double derr_Xib_Asym = sqrt(Xb_errAve**2 + Lb_errAve**2);
  double Xib_Asym = Xb_Ave - Lb_Ave + 0.024;
  double err_Xib_Asym = sqrt(Xb_errAve**2 + Lb_errAve**2 + 0.014**2);
  
  cout << endl;
  
  cout << "Meas Delta Asym = " << dXib_Asym << " +- " <<  derr_Xib_Asym << endl;
  cout << "Meas Asym = " << Xib_Asym << " +- " <<  err_Xib_Asym << "+-" << sqrt(0.009**2+0.014**2) << endl;
  

}
