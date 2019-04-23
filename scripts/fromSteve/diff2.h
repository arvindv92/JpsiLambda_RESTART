364c364
<   h->SetStats(kFALSE);
---
>   //h->SetStats(kFALSE);
558c558
<  LO[i1]=1.5; HI[i1]=6.5; i1++;
---
>  LO[i1]=1.5; HI[i1]=6.5; i1++; 
647a648
>    hmc[i]->Sumw2();
652d652
< 
657,667c657,660
<  TString HB = "#Lambda";
< 
<  int i1 = 0;
<  int i2 = 0;
<  
<  double maxDecayTime = 0.6;
<  double maxZ = 2600;
<  if(lambdaType=="LL"){
<    maxDecayTime = 0.1;
<    maxZ = 800;
<  }
---
>   TString HB = "#Lambda";
>   
>   int i1 = 0;
>   int i2 = 0;
668a662,668
>   double maxDecayTime = 0.6;
>   double maxZ = 2600;
>   if(lambdaType=="LL"){
>     maxDecayTime = 0.1;
>     maxZ = 800;
>   }
>   
694,695d693
<  
< 
727,746d724
<  if(lambdaType=="LL"){
<    HB = "proton (from #Lambda)";
<    histExf = i1;
<    histNames[i1] = "pIP"; histExp[i1] = "p_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
<    LO[i1]=0.0; HI[i1]=4.0; i1++;
<    histNames[i1] = "pIPCHI2"; histExp[i1] = "log(p_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
<    LO[i1]=2.0; HI[i1]=11.0; i1++;
<    histNames[i1] = "pDOCAz"; histExp[i1] = "abs(p_PY*(pi_ORIVX_X-0.459)-p_PX*(pi_ORIVX_Y+0.014))/sqrt(p_PX**2+p_PY**2)"; 
<    hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=30; LO[i1]=0.0; HI[i1]=3.0; i1++;
< 
<    HB = "pion (from #Lambda)";
<    histNames[i1] = "piIP"; histExp[i1] = "pi_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
<    LO[i1]=0.0; HI[i1]=10.0; i1++;
<    histNames[i1] = "piIPCHI2"; histExp[i1] = "log(pi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
<    LO[i1]=2.0; HI[i1]=11.0; i1++;
<    histNames[i1] = "piDOCAz"; histExp[i1] = "abs(pi_PY*(pi_ORIVX_X-0.459)-pi_PX*(pi_ORIVX_Y+0.014))/sqrt(pi_PX**2+pi_PY**2)"; 
<    hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=32; LO[i1]=0.0; HI[i1]=8.0; i1++;
<    nHistsEx = i1 - histExf;
<  }
< 
772,773c750,751
<  histNames[i1] = "XiDira"; histExp[i1] = "1000*acos(Xi_DIRA_ORIVX)"; hTit[i1] = HB+" #thata_{DIRA, ORIVX}"; NB[i1]=30; 
<  LO[i1]=0.0; HI[i1]=15.0; i1++;
---
>  histNames[i1] = "XiDira"; histExp[i1] = "1000*acos(Xi_DIRA_ORIVX)"; hTit[i1] = HB+" #thata_{DIRA, ORIVX}"; NB[i1]=20; 
>  LO[i1]=0.0; HI[i1]=20.0; i1++;
775a754,774
>  if(lambdaType=="LL"){
>    HB = "proton (from #Lambda)";
>    histExf = i1;
>    histNames[i1] = "pIP"; histExp[i1] = "p_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
>    LO[i1]=0.0; HI[i1]=4.0; i1++;
>    histNames[i1] = "pIPCHI2"; histExp[i1] = "log(p_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
>    LO[i1]=2.0; HI[i1]=11.0; i1++;
>    histNames[i1] = "pDOCAz"; histExp[i1] = "abs(p_PY*(pi_ORIVX_X-0.459)-p_PX*(pi_ORIVX_Y+0.014))/sqrt(p_PX**2+p_PY**2)"; 
>    hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=30; LO[i1]=0.0; HI[i1]=3.0; i1++;
> 
>    HB = "pion (from #Lambda)";
>    histNames[i1] = "piIP"; histExp[i1] = "pi_IP_OWNPV"; hTit[i1] = HB+" IP [mm]"; NB[i1]=50; 
>    LO[i1]=0.0; HI[i1]=10.0; i1++;
>    histNames[i1] = "piIPCHI2"; histExp[i1] = "log(pi_IPCHI2_OWNPV)"; hTit[i1] = HB+" log(#chi^{2}_{IP})"; NB[i1]=45; 
>    LO[i1]=2.0; HI[i1]=11.0; i1++;
>    histNames[i1] = "piDOCAz"; histExp[i1] = "abs(pi_PY*(pi_ORIVX_X-0.459)-pi_PX*(pi_ORIVX_Y+0.014))/sqrt(pi_PX**2+pi_PY**2)"; 
>    hTit[i1] = HB+" DOCAz [mm]"; NB[i1]=32; LO[i1]=0.0; HI[i1]=8.0; i1++;
>    nHistsEx = i1 - histExf;
>  }
> 
> 
780d778
<  
782,784d779
< 
<  nHists = i1;
<  histSys = nHists;
790a786
>    hmc[i]->Sumw2();
800a797,800
>   cout << "simCut" << endl;
>   cout << simCut.Print() << endl;
>   cout << "********************************************************************"<<endl;
>   
806,812c806,811
<     if(i==1){
<       hda[i]->Draw();
<       ctest->Update();      
<     }
<     cout << "Integral: " << hda[1]->Integral() << endl;
<     
<     addGraphics(hda[i],hTit[i],"",1);
---
>     //    if(i==1){
>     //hda[i]->Draw();
>     //ctest->Update();      
>       // }
>     cout << "Integral: " << hda[i]->Integral() << endl;
>     //addGraphics(hda[i],hTit[i],"",1);
814,815c813
<     cout << "Filling histogram: " << i << " " << hTit[i] << " 1 " << endl;
<     hmc[i]->Sumw2();
---
>     cout << "Filling histogram: " << i << " " << hTit[i] << " 1 " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
817c815,817
<     addGraphics(hmc[i],hTit[i],"",2);
---
>     cout << "Filling histogram: " << i << " " << hTit[i] << " 1a " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
>     //addGraphics(hmc[i],hTit[i],"",2);
>     cout << "Filling histogram: " << i << " " << hTit[i] << " 1b " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
819c819
<     cout << "Filling histogram: " << i << " " << hTit[i] << " 2 " << endl;
---
>     cout << "Filling histogram: " << i << " " << hTit[i] << " 2 " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
843d842
<     
1282a1282,1303
> 
> void makeLbPlotsLL(){
>   int j,k;
> 
>   TCanvas *ce = new TCanvas("ce","",800,600);
>   ce->Divide(3,2);
>   for(int i=0;i<nHistsEx;i++){
>     ce->cd(i+1);
>     ce->cd(i+1)->SetBottomMargin(bottomMarginOffset);
>     j = histExf + i;
>     hda[j]->Draw("e");
>     hmc[j]->Draw("hist,same");  
>     if(i==0) legend1->Draw();
>   }
>   
>   if(makePDF){
>     makePlots(ce,"ExtraLambda_from_"+parentPart+tag);
>   } 
> }
> 
> 
> 
1289c1310,1311
<   cout << "I am here - 0" << endl;
---
>   cout << "I am here - 0  " << nHistsL << " " << histLf << endl;
>   cout << hda[0] << " " << hmc[0] << endl;
1291a1314,1315
>   
>   
1293,1295c1317,1322
<   TCanvas *c0 = new TCanvas("c0","Lambda0",1200,800);
<   c0->Divide(3,2);
<   c0->cd(1);
---
>   TCanvas *cc0 = new TCanvas("cc0","Lambda0",1200,800);
>   cout << hda[0] << " " << hmc[0] << endl;
>   cc0->Divide(3,2);
>   cout << hda[0] << " " << hmc[0] << endl;
>   cc0->cd(1);
>   cout << hda[0] << " " << hmc[0] << endl;
1296a1324
>   cout << hda[0] << " " << hmc[0] << endl;
1298,1299c1326,1329
<     c0->cd(i+1);
<     c0->cd(i+1)->SetBottomMargin(bottomMarginOffset);
---
>     cout << hda[i] << " " << hmc[i] << endl;
>     cc0->cd(i+1);
>     cout << hda[i] << " " << hmc[i] << endl;
>     cc0->cd(i+1)->SetBottomMargin(bottomMarginOffset);
1304c1334,1336
<     if(i==0) legend1->Draw();
---
>     cout << hda[0] << " " << hmc[0] << endl;
>     //if(i==0) legend1->Draw();
>     //cc0->Update();
1307c1339,1340
< 
---
>   //return;
>   
1322c1355
< 
---
>   ca->Update();
1349c1382
< 
---
>   cb->Update();
1367a1401
>     nRatio++;
1391c1425
<   nRatio = nRatio + 2;
---
>   cc->Update();
1408c1442
< 
---
>   cd->Update();
1421a1456
>     nRatio++;
1431c1466
<   nRatio++;
---
>   cz->Update();
1440c1475
<     makePlots(c0,"Lambda_from_"+parentPart+tag);
---
>     makePlots(cc0,"Lambda_from_"+parentPart+tag);
1455,1473d1489
< void makeLbPlotsLL(){
<   int j,k;
< 
<   TCanvas *ce = new TCanvas("ce","",800,600);
<   ce->Divide(3,2);
<   for(int i=0;i<nHistsEx;i++){
<     ce->cd(i+1);
<     ce->cd(i+1)->SetBottomMargin(bottomMarginOffset);
<     j = histExf + i;
<     hda[j]->Draw("e");
<     hmc[j]->Draw("hist,same");  
<     if(i==0) legend1->Draw();
<   }
<   
<   if(makePDF){
<     makePlots(ce,"ExtraLambda_from_"+parentPart+tag);
<   } 
< }
< 
1477a1494
> 
1518a1536
> 
1526,1527c1544,1545
<     hmc[j]->Draw("hist");  
<     hda[j]->Draw("e,same");
---
>     hda[j]->Draw("e");
>     hmc[j]->Draw("hist,same");  
1530a1549,1550
> 
> 
1551a1572
>     makePlots(cj,"XiVtx_from_"+parentPart+tag);
1896,1897c1917,1922
<   sg1A = new RooRealVar("#sigma_{1A}","",7.56,2,12);
<   sg1scale = new RooRealVar("#sigma_{1B}/#sigma_{1A}","",2.32,0.5,3.5);
---
>   //sg1A_mc = new RooRealVar("#sigma_{1A,mc}","",7.5,2,12);
>   sg1A_mc = new RooRealVar("sg1A_mc","",7.5,2,12);
>   sgscale = new RooRealVar("sf_{mc1}","",1.0,0.9,1.4);//,0.5,1.5);
>   sg1A = new RooFormulaVar("sg1A","","@0*@1",RooArgList(*sg1A_mc,*sgscale));
>   //sg1A = new RooRealVar("#sigma_{1A}","",7.5,2,12);
>   sg1scale = new RooRealVar("#sigma_{1B}/#sigma_{1A}","",1.0,0.5,1.5);
1899,1900c1924,1927
<   alpha1A = new RooRealVar("#alpha_{1A}","",1.5,0.5,3.0);
<   alpha1B = new RooRealVar("#alpha_{1B}","",-2.0,-5.0,-0.2);
---
>   minusOne = new RooRealVar("minusOne","",-1.0);
>   alpha1A = new RooRealVar("#alpha_{1A}","",1.0,0.5,3.0);
>   //alpha1B = new RooFormulaVar("alpha1B","","@0*@1",RooArgList(*minusOne,*alpha1A));
>   alpha1B = new RooRealVar("#alpha_{1B}","",-1.3,-5.0,-0.2);
1902,1903c1929,1940
<   gauss1A = new RooCBShape("gauss1A","Signal",*mvar1,*mg1,*sg1A,*alpha1A,*nCB);
<   gauss1B = new RooCBShape("gauss1B","Signal",*mvar1,*mg1,*sg1B,*alpha1B,*nCB);
---
>   fr1A = new RooRealVar("f_{1A}","fraction",0.50,0.2,0.8);
>   fr1B = new RooRealVar("f_{1B}","fraction",0.3,0.1,0.7);
> 
>   if(!altFitFunction){
>     gauss1A = new RooCBShape("gauss1A","Signal",*mvar1,*mg1,*sg1A,*alpha1A,*nCB);
>     gauss1B = new RooCBShape("gauss1B","Signal",*mvar1,*mg1,*sg1B,*alpha1B,*nCB);
>   }else{ 
>     gauss1A = new RooGaussian("gauss1A","Signal",*mvar1,*mg1,*sg1A);
>     gauss1B = new RooGaussian("gauss1B","Signal",*mvar1,*mg1,*sg1B);
>   }
>   
>   sg1C = new RooRealVar("#sigma_{1C}","",8,5,15);
1906,1907c1943,1944
<   fr1 = new RooRealVar("f_{1A}","fraction",0.5);//,0.3,0.7);
<   sig1 = new RooAddPdf("sig1","Signal",RooArgList(*gauss1A,*gauss1B),RooArgList(*fr1));
---
>   //sig1 = new RooAddPdf("sig1","Signal",RooArgList(*gauss1A,*gauss1B,*gauss1C),RooArgList(*fr1A,*fr1B));
>   sig1 = new RooAddPdf("sig1","Signal",RooArgList(*gauss1A,*gauss1B),RooArgList(*fr1A));
1909c1946,1952
<   comb1 = new RooExponential("comb1","f_{back}",*mvar1,*b1);
---
>   if(!altBackFunction){
>     comb1 = new RooExponential("comb1","f_{back}",*mvar1,*b1);
>   }else{
>     comb1 = new RooChebychev("comb1","f_{back}",*mvar1,*b1);
>   }
>   
>   
1914a1958,1959
> 
> 
1921,1922c1966,1969
<   sg2A = new RooRealVar("#sigma_{2A}","",sg1A->getVal(),3,12);
<   sg2scale = new RooRealVar("#sigma_{2B}/#sigma_{2A}","",sg1scale->getVal(),1.0,3.5);
---
>   sg2A_mc = new RooRealVar("sg2A_mc","",7.5,2,12);
>   sgscale2 = new RooRealVar("sf_{mc2}","",1.0,0.9,1.4);//,0.5,1.5);
>   sg2A = new RooFormulaVar("sg2A","","@0*@1",RooArgList(*sg2A_mc,*sgscale2));
>   sg2scale = new RooRealVar("#sigma_{2B}/#sigma_{2A}","",1.0);//,0.5,1.5);
1926,1929d1972
<   gauss2A = new RooCBShape("gauss2A","Signal",*mvar2,*mg2,*sg2A,*alpha2A,*nCB);
<   gauss2B = new RooCBShape("gauss2B","Signal",*mvar2,*mg2,*sg2B,*alpha2B,*nCB);
<   //gauss2A = new RooGaussian("gauss2A","Signal",*mvar2,*mg2,*sg2A);
<   //gauss2B = new RooGaussian("gauss2B","Signal",*mvar2,*mg2,*sg2B);
1930a1974,1982
> 
>   if(!altFitFunction){
>     gauss2A = new RooCBShape("gauss2A","Signal",*mvar2,*mg2,*sg2A,*alpha2A,*nCB);
>     gauss2B = new RooCBShape("gauss2B","Signal",*mvar2,*mg2,*sg2B,*alpha2B,*nCB);
>   }else{ 
>     gauss2A = new RooGaussian("gauss2A","Signal",*mvar2,*mg2,*sg2A);
>     gauss2B = new RooGaussian("gauss2B","Signal",*mvar2,*mg2,*sg2B);
>   }
> 
1932d1983
<   //sig2 = new RooGaussian("sig2","Signal",*mvar2,*mg2,*sg2A);
1934c1985,1990
<   comb2 = new RooExponential("comb2","f_{back}",*mvar2,*b2);
---
>   if(!altBackFunction){
>     comb2 = new RooExponential("comb2","f_{back}",*mvar2,*b2);
>   }else{
>     comb2 = new RooChebychev("comb2","f_{back}",*mvar2,*b2);
>   }
>   //comb2 = new RooExponential("comb2","f_{back}",*mvar2,*b2);
1942a1999,2006
>   if(!mc){
>     alpha2A->setVal(1.08);
>     alpha2B->setVal(-1.17);
>     alpha2A->setConstant();
>     alpha2B->setConstant();    
>     sg2scale->setVal(0.855);
>     sg2scale->setConstant();
>   }
1945a2010
> void setXibDDLparams(){
1947c2012
< void doLbFitAndPlot(){
---
>   if(altFitFunction) return;
1948a2014,2055
>   if(year=="Run1"){
>     sg2A_mc->setVal(7.1);
>     sg2A_mc->setConstant();
>   }else if(year=="Run2"){
>     sg2A_mc->setVal(7.45);
>     sg2A_mc->setConstant();
>   }
> 
>   if(mc) return;
> 
>   if(year=="Run1"){
>     alpha2A->setVal(0.96);
>     alpha2B->setVal(-1.155);
>     sg2scale->setVal(1.0);
>     fr2->setVal(0.5);
>     alpha2A->setConstant();
>     alpha2B->setConstant();    
>     sg2scale->setConstant();
>     fr2->setConstant();
>   }else if(year=="Run2"){
>     alpha2A->setVal(0.91);
>     alpha2B->setVal(-1.05);
>     sg2scale->setVal(1.0);
>     fr2->setVal(0.5);
>     alpha2A->setConstant();
>     alpha2B->setConstant();    
>     sg2scale->setConstant();
>     fr2->setConstant();
>   }
>   
> 
>   int ib1 = h2L->FindBin(5750);
>   int ib2 = h2L->FindBin(5840);
>   double ns = h2L->Integral(ib1,ib2);
>   double nb = h2L->Integral()-ns;
>   nsig2->setVal(ns);
>   ncomb2->setVal(nb);
> 
> 
> }
> 
> void setXibDDDparams(){
1949a2057,2208
>   if(altFitFunction) return;
> 
>   if(year=="Run1"){
>     sg2A_mc->setVal(8.0);
>     sg2A_mc->setConstant();
>   }else if(year=="Run2"){
>     sg2A_mc->setVal(8.4);
>     sg2A_mc->setConstant();
>   }
> 
>   if(mc) return;
> 
>   if(year=="Run1"){
>     alpha2A->setVal(1.14);
>     alpha2B->setVal(-1.25);
>     alpha2A->setConstant();
>     alpha2B->setConstant();    
>     sg2scale->setVal(1.0);
>     sg2scale->setConstant();
>     fr2->setVal(0.5);
>   }else if(year=="Run2"){
>     alpha2A->setVal(1.08);
>     alpha2B->setVal(-1.154);
>     alpha2A->setConstant();
>     alpha2B->setConstant();    
>     sg2scale->setVal(1.0);
>     sg2scale->setConstant();
>     fr2->setVal(0.5);
>   }
>   
> 
>   int ib1 = h2D->FindBin(5750);
>   int ib2 = h2D->FindBin(5840);
>   double ns = h2D->Integral(ib1,ib2);
>   double nb = h2D->Integral()-ns;
>   nsig2->setVal(ns);
>   ncomb2->setVal(nb);
> 
> 
> }
> 
> 
> void setLbDDparams(){
> 
>   if(altFitFunction) return;
> 
>   sg1scale->setVal(1.0);
>   sg1scale->setConstant();
>   if(year=="Run1"){
>     sg1A_mc->setVal(7.5);
>     sg1A_mc->setConstant();
>   }else if(year=="Run2"){
>     sg1A_mc->setVal(7.5);
>     sg1A_mc->setConstant();
>   }
>   
>   if(mc) return;
> 
>   if(year=="Run1"){
>     alpha1A->setVal(1.02);
>     alpha1B->setVal(-1.31);
>     fr1A->setVal(0.35); //0.35
>     fr1B->setVal(0.28);
>     sg1C->setVal(8.4);
>     alpha1A->setConstant();
>     alpha1B->setConstant();
>     sg1C->setConstant();
>     fr1A->setConstant();
>     fr1B->setConstant();
>   }else if(year=="Run2"){
>     alpha1A->setVal(0.99);
>     alpha1B->setVal(-1.216);
>     fr1A->setVal(0.40);
>     fr1B->setVal(0.28);
>     sg1C->setVal(8.4);
>     alpha1A->setConstant();
>     alpha1B->setConstant();
>     sg1C->setConstant();
>     fr1A->setConstant();
>     fr1B->setConstant();
>   }
>   
> 
>   int ib1 = h1->FindBin(5600);
>   int ib2 = h1->FindBin(5640);
>   double ns = h1->Integral(ib1,ib2);
>   double nb = h1->Integral()-ns;
>   nsig1->setVal(ns);
>   ncomb1->setVal(nb);
> 
> 
> }
> 
> 
> void setLbLLparams(){
> 
>   if(altFitFunction) return;
>   
>   sg1scale->setVal(1.0);
>   sg1scale->setConstant();
> 
>   if(year=="Run1"){
>     //sg1A_mc->setVal(7.5);
>     sg1A_mc->setVal(6.96);
>     sg1A_mc->setConstant();
>   }else if(year=="Run2"){
>     sg1A_mc->setVal(7.47);
>     sg1A_mc->setConstant();
>   }
> 
> 
>   if(mc) return;
> 
>   if(year=="Run1"){
>     alpha1A->setVal(0.82);
>     alpha1B->setVal(-1.06);
>     fr1A->setVal(0.37); //0.37
>     fr1B->setVal(0.23);
>     sg1C->setVal(5.6);
>     alpha1A->setConstant();
>     alpha1B->setConstant();    
>     sg1C->setConstant();
>     fr1A->setConstant();
>     fr1B->setConstant();
>   }else if(year=="Run2"){
>     alpha1A->setVal(1.02);
>     alpha1B->setVal(-0.940);
>     fr1A->setVal(0.59);
>     fr1B->setVal(0.23);
>     sg1C->setVal(5.6);
>     alpha1A->setConstant();
>     alpha1B->setConstant();    
>     sg1C->setConstant();
>     fr1A->setConstant();
>     fr1B->setConstant();
>   }
>   
> 
>   int ib1 = h1LL->FindBin(5600);
>   int ib2 = h1LL->FindBin(5640);
>   double ns = h1LL->Integral(ib1,ib2);
>   double nb = h1LL->Integral()-ns;
>   nsig1->setVal(ns);
>   ncomb1->setVal(nb);
> 
> 
> }
> 
> 
> void doLbFitAndPlot(){
> 
>   setLbDDparams();
1961a2221
>   suff = suff + "_" + year;
1970a2231
>   setLbLLparams();
1980a2242
>   suff = suff + "_" + year;
2111,2122c2373,2385
<   /*
<   pdf2->fitTo(*rdh2,Hesse(kTRUE),Strategy(1));
<   mass = mg2->getVal();
<   err_mass = mg2->getError();
<   TCanvas* cc2 = makeRooPlot(pdf2,rdh2,mvar2,1,nsig2,modes[1],mLowXib,mHiXib,nbinXib);
<   cc2->Draw();
<   cc2->Update();
<   yieldXib = nsig2->getVal();
<   err_yieldXib = nsig2->getError();
<   fitmass_Xib = mg2->getVal();
<   err_fitmass_Xib = mg2->getError();
<   */
---
>   setXibDDDparams();  
>   if(!mc){
>     pdf2->fitTo(*rdh2,Hesse(kTRUE),Strategy(1));
>     mass = mg2->getVal();
>     err_mass = mg2->getError();
>     TCanvas* cc2 = makeRooPlot(pdf2,rdh2,mvar2,1,nsig2,modes[1],mLowXib,mHiXib,nbinXib);
>     cc2->Draw();
>     cc2->Update();
>     yieldXib = nsig2->getVal();
>     err_yieldXib = nsig2->getError();
>     fitmass_Xib = mg2->getVal();
>     err_fitmass_Xib = mg2->getError();
>   }
2123a2387
>   setXibDDDparams();  
2134a2399
>   setXibDDLparams();  
2148c2413,2414
<     TString suff = "XiLorD";
---
>     TString suff = "XiLorD_"+year;
>     /*
2154c2420,2421
<     suff = "XiL";
---
>     */
>     suff = "XiL_"+year;
2159,2160c2426,2427
<     cc2K->Print("./Plots/XibMass_LambdaDD_"+suff+".eps");
<     suff = "XiD";
---
>     cc2L->Print("./Plots/XibMass_LambdaDD_"+suff+".eps");
>     suff = "XiD_"+year;
2225a2493,2505
> 
> void applySys(){
>   makePDF = false;
> 
>   if(halfBinWidth == true){
>     bwLb = bwLb / 2.0;
>     bwXib = bwXib / 2.0;
>   }
> 
>   if(noReweighting) mcOpt = 1;
>   
> 
> }
