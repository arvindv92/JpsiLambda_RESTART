#include "CodeInput.h"
bool makePDF2 = true;

void bookLbBaseHistos2(){

  //signalWeight = "nsig1_sw";  
  TString HB = "#Lambda";

  int i1 = 0;
 int i2 = 0;
 

 // Lambda_b                                                            
 HB = "#Lambda_{b}^{0}";
 histNames[i1] = "taulamb"; histExp[i1] = "1000*Lb_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=80; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 histNames[i1] = "diralb"; histExp[i1] = "abs(Lb_DIRA_OWNPV)"; hTit[i1] = HB+" #theta_{DIRA}"; NB[i1]=100; 
 LO[i1]=0.999; HI[i1]=1.00; i1++;
 // Lambda0
 HB = "#Lambda";
 histNames[i1] = "fdslam"; histExp[i1] = "log((L_ENDVERTEX_Z-L_OWNPV_Z)/L_ENDVERTEX_ZERR)"; hTit[i1] = HB+" FD Significance"; NB[i1]=50; 
 //histNames[i1] = "fdslam"; histExp[i1] = "L_ENDVERTEX_Z-L_OWNPV_Z"; hTit[i1] = HB+" FD Significance"; NB[i1]=50; 
 LO[i1]=0.0; HI[i1]=10.0; i1++;
 HB = "proton (from #Lambda)";
 histNames[i1] = "pp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2+p_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=120.0; i1++;
 histNames[i1] = "ptp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 HB = "pion (from #Lambda)";
 histNames[i1] = "ppi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2+pi_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
 LO[i1]=0.0; HI[i1]=30.0; i1++;
 histNames[i1] = "ptpi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=2.0; i1++;

 nHists = i1;
 //nHists = 2;
 
 cout << "# Histograms = " << nHists << endl;
 
 
 // Book the histograms
 for(int i=0; i<nHists; i++){
   hda[i] = new TH1F(histNames[i]+"_da",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i] = new TH1F(histNames[i]+"_mc",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i]->Sumw2();
 }


}



void fillHistos2(TTree* tree_data, TTree *tree_sim){
  
  TString hist;
  double vint;
  //simCut = MCLambdabDD_Cut;
  TCanvas *ctest = new TCanvas("ctest","",600,600);
  cout << "simCut" << endl;
  cout << simCut.Print() << endl;
  cout << "********************************************************************"<<endl;
  
  


  for(int i=0; i<nHists; i++){
    hist = histNames[i]+"_da";
    cout << "Filling histogram: " << i << " " << hTit[i] << " " << histExp[i] << " 0 " << " " << signalWeight << " " << tree_data << endl;
    tree_data->Draw(histExp[i]+">>"+hist,signalWeight);
    cout << "Integral: " << hda[i]->Integral() << endl;
    addGraphics(hda[i],hTit[i],"",4);
    hist = histNames[i]+"_mc";
    cout << "Filling histogram: " << i << " " << hTit[i] << " 1 " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
    tree_sim->Draw(histExp[i]+">>"+hist,simCut);    
    cout << "Filling histogram: " << i << " " << hTit[i] << " 1a " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
    addGraphics(hmc[i],hTit[i],"",2);
    cout << "Filling histogram: " << i << " " << hTit[i] << " 1b " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
    hmc[i]->Scale(hda[i]->Integral()/hmc[i]->Integral());
    cout << "Filling histogram: " << i << " " << hTit[i] << " 2 " << " " << tree_sim << " " << hist << " " << histExp[i] << " " << hmc[i] << endl;
    hda[i]->Draw();
  }
  
}

void makeBasePlots2(){
  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);


  cout << "============================="<<endl;
  cout << "I am here - 0  " << nHistsL << " " << histLf << endl;
  cout << hda[0] << " " << hmc[0] << endl;
  cout << "============================="<<endl;

  
  legend1 = new TLegend(0.25,0.70,0.88,0.88);
  legend1->SetFillStyle(0);
  legend1->SetBorderSize(0);
  legend1->SetFillStyle(0);
  legend1->SetTextSize(0.065);
  legend1->AddEntry(hda[0],"Sideband Data","L");
  legend1->AddEntry(hmc[0],"Signal - Simulation","L");
  legend1->Draw();
  
  int j;
  // Plot results
  TCanvas* cc0 = new TCanvas("cc0","Lambda_0",1200,1000);
  cc0->Divide(3,3);
  for(int i=0;i<nHists;i++){
    cout << hda[i] << " " << hmc[i] << endl;
    cc0->cd(i+1);
    cout << hda[i] << " " << hmc[i] << endl;
    cc0->cd(i+1)->SetBottomMargin(bottomMarginOffset);
    j = i;
    if(i!=1){
      hda[j]->Draw("hist");
      hmc[j]->Draw("hist,same");  
    }else{
      cc0->cd(i+1)->SetLogy();
      hmc[j]->Draw("hist");
      hda[j]->Draw("hist,same");  
    }
    
    cout << hda[0] << " " << hmc[0] << endl;
    //if(i==0) legend1->Draw();
    //cc0->Update();
  }
  cc0->cd(1);
  legend1->Draw();
  
}

void compSigBackLb(){
  year = "Run1";
  lambdaType = "DD";
  parentPart = "Lb";  
  mcOpt = 5;  //def = 5

  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);

  // Get Cuts
  getCuts();
  simCut = MCLambdabDD_Cut;

  // Get needed files
  getLbDDFiles();
  TreeWithWeight = (TTree*)fsigdata->Get("MyTuple");    // Tuple with all variables

  // Get the histogram weights
  getWeightHistos();  

  // Get the MC tree for plotting MC results
  if(mcOpt>=2){
    getNewTree(chainl);
  }else{
    newtree = (TTree*)chainl;
  }
  
  // print out cuts
  cout << simCut.Print() << endl;  

  // Book histograms
  bookLbBaseHistos2();

  // Fill histograms
  signalWeight  = LambdabDD_Cut&&"abs(Lb_DTF_M_JpsiLConstr-5620)>50&&abs(Lb_DTF_M_JpsiLConstr-5620)<100";
  fillHistos2(TreeWithWeight,newtree);

  // make Plots
  makeBasePlots2();

  if(makePDF2){
    cc0->Print("./Plots/CompSigBackLb.png");
  }


}

void compSigBackXib(){
  year = "Run1";
  lambdaType = "DD";
  parentPart = "Xib";  
  mcOpt = 5;  //def = 5

  TLatex *myLatex = new TLatex();
  myLatex->SetTextFont(42); myLatex->SetTextColor(1); myLatex->SetTextAlign(12); myLatex->SetNDC(kTRUE); myLatex->SetTextSize(0.06);

  // Get Cuts
  getCuts();
  simCut = MCXib_Cut&&XiL_Cut;

  // Get needed files
  getXibFiles();
  TreeWithWeight = (TTree*)fsigdata->Get("MyTuple");    // Tuple with all variables

  // Get the histogram weights
  getWeightHistos();  

  // Get the MC tree for plotting MC results
  if(mcOpt>=2){
    getNewTreeXib(chainl);
  }else{
    newtree = (TTree*)chainl;
  }
  
  // print out cuts
  cout << simCut.Print() << endl;  

  // Book histograms
  bookXibBaseHistos2();

  // Fill histograms
  signalWeight  = Xib_Cut&&XiL_Cut&&"abs(Xib_DTF_M_JpsiXiLConstr-5797)>50&&abs(Xib_DTF_M_JpsiXiLConstr-5797)<150";
  fillHistos2(TreeWithWeight,newtree);

  // make Plots
  makeBasePlots2();

  if(makePDF2){
    cc0->Print("./Plots/CompSigBackXib.png");
  }


}


void bookXibBaseHistos2(){

  //signalWeight = "nsig1_sw";  
  TString HB = "#Lambda";

  int i1 = 0;
 int i2 = 0;
 

 // Lambda_b                                                            
 HB = "#Xi_{b}^{-}";
 histNames[i1] = "taulamb"; histExp[i1] = "1000*Xib_TAU"; hTit[i1] = HB+" decay time"; NB[i1]=80; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 histNames[i1] = "diralb"; histExp[i1] = "abs(Xib_DIRA_OWNPV)"; hTit[i1] = HB+" DIRA"; NB[i1]=100; 
 LO[i1]=0.999; HI[i1]=1.00; i1++;
 // Lambda0
 HB = "#Lambda";
 histNames[i1] = "fdslam"; histExp[i1] = "log((L_ENDVERTEX_Z-L_OWNPV_Z)/L_ENDVERTEX_ZERR)"; hTit[i1] = HB+" FD Significance"; NB[i1]=50; 
 //histNames[i1] = "fdslam"; histExp[i1] = "L_ENDVERTEX_Z-L_OWNPV_Z"; hTit[i1] = HB+" FD Significance"; NB[i1]=50; 
 LO[i1]=0.0; HI[i1]=10.0; i1++;
 HB = "proton (from #Lambda)";
 histNames[i1] = "pp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2+p_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=120.0; i1++;
 histNames[i1] = "ptp"; histExp[i1] = "0.001*sqrt(p_PX**2+p_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=8.0; i1++;
 HB = "pion (from #Lambda)";
 histNames[i1] = "ppi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2+pi_PZ**2)"; hTit[i1] = HB+" momentum [GeV/#it{c}]"; NB[i1]=60; 
 LO[i1]=0.0; HI[i1]=30.0; i1++;
 histNames[i1] = "ptpi"; histExp[i1] = "0.001*sqrt(pi_PX**2+pi_PY**2)"; hTit[i1] = HB+" p_{T} [GeV/#it{c}]"; NB[i1]=40; 
 LO[i1]=0.0; HI[i1]=2.0; i1++;

 nHists = i1;
 //nHists = 2;
 
 cout << "# Histograms = " << nHists << endl;
 
 
 // Book the histograms
 for(int i=0; i<nHists; i++){
   hda[i] = new TH1F(histNames[i]+"_da",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i] = new TH1F(histNames[i]+"_mc",hTit[i],NB[i],LO[i],HI[i]);
   hmc[i]->Sumw2();
 }


}

