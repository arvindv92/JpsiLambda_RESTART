void checkDup(TString mode = "Lb",TString run = "Run1"){

  TString infile;
  UInt_t          runNumber;
  ULong64_t       eventNumber;
  Double_t        Lb_DTF_M_JpsiLConstr;
  Double_t        Xib_DTF_M_JpsiXiLConstr;
  Double_t mass, mLo, mHi;
  

  if(mode=="Xib"){
    infile = "Xib2JpsiXi_"+run+"_DD_L0JPsiTOS_NewStrip.root";
  }else if(mode=="Lb"){
    infile = "Lb2JpsiL0_"+run+"_DD_L0JPsiTOS.root";
  }
  TString treeName = "MyTuple";
  
  
  
  cout << " -------------------------------------------------------------------------" << endl;
  //cout << "Running tuning of input MC variables, Decay.= " << decay << endl;
  cout << "Input File: " << infile << endl;
  cout << "----> Tree: " << treeName << endl;
  cout << " -------------------------------------------------------------------------" << endl;

  TFile *f = new TFile(infile);
  
  TTree *oldtree = (TTree*)f->Get(treeName);
  oldtree->SetBranchAddress("runNumber", &runNumber);  
  oldtree->SetBranchAddress("eventNumber", &eventNumber);  

  if(mode=="Xib"){
    oldtree->SetBranchAddress("Xib_DTF_M_JpsiXiLConstr", &mass);  
    mLo = 5640; mHi = 5950;
  }else if(mode=="Lb"){
    oldtree->SetBranchAddress("Lb_DTF_M_JpsiLConstr", &mass);  
    mLo = 5500; mHi = 5750;
  }


  Long64_t nentries = oldtree->GetEntriesFast();
  cout << "Looping over " << nentries << " entries" << endl;
   
   //nentries = 1000;
  TH1F* hdup = new TH1F("hdup","Duplicates",10,0.0,10.0);
  
  int ndup = 0;
  ULong64_t lastEv = 0;
  UInt_t    lastRun = 0;

   double v;
   for (Long64_t i=0;i<nentries; i++) {
      oldtree->GetEntry(i);
      if(i%10000 == 0) cout << "At entry = " << i << endl;
      if(mass<mLo || mass>mHi) continue;
      
      if(runNumber==lastRun && eventNumber==lastEv) ndup++;

      if(!(runNumber==lastRun && eventNumber==lastEv)){
        //if(!(runNumber==lastRun)){        
        lastEv = eventNumber;
        lastRun = runNumber;
        hdup->Fill(ndup);
        ndup = 1;
      }
   }

   hdup->Draw();
   
}
