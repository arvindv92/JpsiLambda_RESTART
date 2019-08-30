#include <TFile.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TCut.h>
#include <iostream>

using namespace std;
void run()
{
   TChain *h1 = new TChain("B2JpsiKTree/MyTuple");

    for(Int_t i=0; i<=398;i++)
    {
      h1->Add(TString::Format("../2016_MagDown/1170_%d/jpsik.root",i),-1);
    }

  cout<<"DONE ATTACHING ROOT FILES"<<endl;

  h1->MakeClass("cuts");
}
