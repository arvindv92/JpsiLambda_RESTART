void asym2018_062(){
  
  int sample = 12;

  double asym[2] = {};
  



  TCanvas *c = new TCanvas("c","",800,800);
  c->Divide(2,2);
  c->cd(1);

  if(sample==1){  
    double p[2] ={6952, 6563};
    double ep[2] ={100, 96};
  }else if(sample==2){
    double p[2] ={7493, 7118};
    double ep[2] ={108, 105};
  }else if(sample==11){
    double p[2] ={254, 236};
    double ep[2] ={19, 18};
  }else if(sample==12){
    double p[2] ={312, 342};
    double ep[2] ={22, 22};
  }
  
  cout << p[0] << " " << p[1] << endl;  
  
  double v1, v2, a;
  TH1F *h1 = new TH1F("h1","",300,-0.3,0.3);
  
  for(int i=0; i<20000;i++){
    v1 = p[0] + ep[0]*gRandom->Gaus();
    v2 = p[1] + ep[1]*gRandom->Gaus();
    a = (v1-v2)/(v1+v2);
    h1->Fill(a);
  }
  
  h1->Draw();
  

  double as[2] = {-0.11, 3.44};
  double eaStat[2] = {2.53, 1.61};
  double eaSyst[2] = {1.08, 0.76};

  TH1F *h2 = new TH1F("h2","",400,-20,20);
  TH1F *h3 = new TH1F("h3","",400,-20,20);
  TH1F *h4 = new TH1F("h4","",400,-20,20);
  
  for(int i=0 ; i<50000; i++){
    double v1 = gRandom->Gaus();
    double v2 = gRandom->Gaus();
    double v3 = gRandom->Gaus();
    double v4 = gRandom->Gaus();
    double a1 = as[0] + v1*eaStat[0];
    double a2 = as[1] + v2*eaStat[1];
    double a_ave = (a1/eaStat[0]**2 + a2/eaStat[1]**2) / (1./eaStat[0]**2 + 1./eaStat[1]**2) ;
    h2->Fill(a_ave);
    a1 = as[0] + v1*eaStat[0] + v3*eaSyst[0];
    a2 = as[1] + v2*eaStat[1] + v3*eaSyst[1];
    a_ave = (a1/eaStat[0]**2 + a2/eaStat[1]**2) / (1./eaStat[0]**2 + 1./eaStat[1]**2) ;
    h3->Fill(a_ave);
    a1 = as[0] + v1*eaStat[0] + v3*eaSyst[0];
    a2 = as[1] + v2*eaStat[1] + v4*eaSyst[1];
    a_ave = (a1/eaStat[0]**2 + a2/eaStat[1]**2) / (1./eaStat[0]**2 + 1./eaStat[1]**2) ;
    h4->Fill(a_ave);
  }

  c->cd(2);  
  h2->Draw();
  c->cd(3);  
  h3->Draw();
  c->cd(4);  
  h4->Draw();

  double statErr = h2->GetRMS();
  double systErr = sqrt(h3->GetRMS()**2-h2->GetRMS()**2);
  double mean = h2->GetMean();
  
  cout << "Asym = " << mean << "+-" << statErr << "+-" << systErr << "%" << endl;
  


}


