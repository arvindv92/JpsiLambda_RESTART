#include <TROOT.h>
#include <TSystem.h>

void master()
{
  // gROOT->ProcessLine(".L trigger.C++");
  // gROOT->ProcessLine("trigger t");
  // gROOT->ProcessLine("t.Loop()");

  //  gROOT->Reset();

  gROOT->ProcessLine(".L cuts.C++");
  gROOT->ProcessLine("cuts t");
  gROOT->ProcessLine("t.Loop()");

  // gROOT->Reset();

  // gROOT->ProcessLine(".L separate.C++");
  // gROOT->ProcessLine("separate()");

}
