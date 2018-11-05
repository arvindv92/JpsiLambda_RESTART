#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooStats/SPlot.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooAddPdf.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include "TCanvas.h"
#include "RooAbsPdf.h"
#include "RooFit.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"
#include "TCut.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TMath.h"
// use this order for safety on library loading
using namespace RooFit ;
using namespace RooStats ;
using namespace std;

// see below for implementation
void AddModel(RooWorkspace*);
void AddData(RooWorkspace*);
void DoSPlot(RooWorkspace*);
void MakePlots(RooWorkspace*);

void Lb_splot_LL(Int_t run = 1, Int_t trackType = 3)
/*
   run = 1/2 for Run 1/2 data/MC. Run 1 = 2011,2012 for both data and MC. Run 2 = 2015,2016 for MC, 2015,2016,2017,2018 for data
   isData = 1 for data, 0 for MC
   mcType = 0 when running over data. When running over MC, mcType = 1 for JpsiLambda, 2 for JpsiSigma, 3 for JpsiXi.
   trackType = 3 for LL, 5 for DD.
 */
{
    // Create a new workspace to manage the project.
    RooWorkspace* wspace = new RooWorkspace("myWS");

    // add the signal and background models to the workspace.
    // Inside this function you will find a discription our model.
    AddModel(wspace);

    // add some toy data to the workspace
    AddData(wspace);

    // inspect the workspace if you wish
    //  wspace->Print();

    // do sPlot.
    //This wil make a new dataset with sWeights added for every event.
    DoSPlot(wspace);

    // Make some plots showing the discriminating variable and
    // the control variable after unfolding.
    MakePlots(wspace);

    // cleanup
    delete wspace;
}

//____________________________________
void AddModel(RooWorkspace* ws){

    cout<<"Starting AddModel()"<<endl;
    // Make models for signal and background

    // set range of observable
    Double_t lowRange = 5200.0, highRange = 6000.0;//highRange changed from 6800

    // make a RooRealVar for the observables
    RooRealVar Lb_DTF_M_JpsiLConstr("Lb_DTF_M_JpsiLConstr", "M_{inv}", lowRange, highRange,"MeV");//discriminating variable

    Int_t nbins = (Int_t)(highRange-lowRange)/4;
    Lb_DTF_M_JpsiLConstr.setBins(nbins);

    // SIGNAL MODEL
    cout << "Making signal model" << endl;
    RooRealVar mean("mean","Gaussian Mean",5619.42,5619.0,5620.0);
    RooRealVar sigma1("sigma1","Gaussian sigma1",7.5,1.0,20.0);
    RooRealVar sigma2("sigma2","Gaussian sigma2",10.0,1.0,20.0);

    RooGaussian sig1("sig1","Gaussian signal1",Lb_DTF_M_JpsiLConstr,mean,sigma1);
    RooGaussian sig2("sig2","Gaussian signal2",Lb_DTF_M_JpsiLConstr,mean,sigma2);

    RooRealVar frac1("frac1","Fraction of sig1 in signal",0.5,0.1,0.9);

    RooAddPdf sig("sig","Gaussian signal",RooArgList(sig1,sig2),frac1);

    // BACKGROUND MODEL
    cout << "make background model" << endl;
    RooRealVar c0("c0","c0",-0.1,-1.0,-1.0);
    RooRealVar c1("c1","c1",0.1,-1.0,1.0);

    RooChebychev bkg("bkg","bkg",Lb_DTF_M_JpsiLConstr,RooArgList(c0,c1));

    // COMBINED MODEL
    cout << "Making full model" << endl;

    RooRealVar sigYield("sigYield","fitted yield for sig",5000 ,0, 15000) ;
    RooRealVar bkgYield("bkgYield","fitted yield for bkg",100000 ,0, 1.1e+06) ;//adjust this based on nentries

    RooAddPdf model("model","signal+background models", RooArgList(sig, bkg), RooArgList(sigYield,bkgYield));

    // interesting for debugging and visualizing the model
    //model.graphVizTree("fullModel.dot");

    cout << "Importing model into workspace" << endl;

    ws->import(model);
    ws->Print();

    cout<<"AddModel() complete"<<endl;
}

void AddData(RooWorkspace* ws){

    cout<<"Starting AddData()"<<endl;
    // get what we need out of the workspace to make toy data

    RooRealVar Lb_DTF_M_JpsiLConstr = *(ws->var("Lb_DTF_M_JpsiLConstr"));

    const char* cuts = "";
    //Import the data
    TFile *f = new TFile("../jpsilambda_LL_forTMVAisoApp.root","READ");

    TTree* tree = (TTree*)f->Get("MyTuple");

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("Lb_DTF_M_JpsiLConstr",1);
    tree->SetBranchStatus("pi_PT",1);

    cout<<tree->GetEntries()<<endl;
    cout<<"poop4"<<endl;

    RooArgSet s(Lb_DTF_M_JpsiLConstr,pi_PT);

    RooDataSet data("data","dataset with x",tree,s) ;

    // import data into workspace
    cout<<"poop5"<<endl;
    ws->import(data, Rename("data"));
    ws->Print();
    cout<<"poop6"<<endl;

}

//____________________________________
void DoSPlot(RooWorkspace* ws){

    cout << "Calculate sWeights" << endl;

    // get what we need out of the workspace to do the fit
    RooAbsPdf* model = ws->pdf("model");
    RooRealVar* sigYield = ws->var("sigYield");
    RooRealVar* bkgYield = ws->var("bkgYield");
    RooDataSet* data = (RooDataSet*) ws->data("data");

    // fit the model to the data.
    model->fitTo(*data, Extended());

    // The sPlot technique requires that we fix the parameters
    // of the model that are not yields after doing the fit.
    RooRealVar* sigma1 = ws->var("sigma1");
    RooRealVar* sigma2 = ws->var("sigma2");

    RooRealVar* mean = ws->var("mean");
    RooRealVar* c0 = ws->var("c0");
    RooRealVar* c1 = ws->var("c1");

    // RooRealVar* c2 = ws->var("c2");
    // RooRealVar* c3 = ws->var("c3");

    mean->setConstant(kTRUE);
    sigma1->setConstant(kTRUE);
    sigma2->setConstant(kTRUE);

    c0->setConstant(kTRUE);
    c1->setConstant(kTRUE);
    // c2->setConstant(kTRUE);
    // c3->setConstant(kTRUE);

    RooMsgService::instance().setSilentMode(true);


    // Now we use the SPlot class to add SWeights to our data set
    // based on our model and our yield variables
    RooStats::SPlot* sData = new RooStats::SPlot("sData","An SPlot",*data, model, RooArgList(*sigYield,*bkgYield) );


    // Check that our weights have the desired properties

    cout << "Check SWeights:" << endl;


    cout << endl <<  "Yield of signal is "
    << sigYield->getVal() << ".  From sWeights it is "
    << sData->GetYieldFromSWeight("sigYield") << endl;


    cout << "Yield of background is "
    << bkgYield->getVal() << ".  From sWeights it is "
    << sData->GetYieldFromSWeight("bkgYield") << endl
    << endl;

    for(Int_t i=0; i < 10; i++)
    {
        cout << "z Weight   " << sData->GetSWeight(i,"sigYield")
        << "   qcd Weight   " << sData->GetSWeight(i,"bkgYield")
        << "  Total Weight   " << sData->GetSumOfEventSWeight(i)
        << endl;
    }

    cout << endl;

    // import this new dataset with sWeights
    cout << "import new dataset with sWeights" << endl;
    ws->import(*data, Rename("dataWithSWeights"));
    cout <<" after splot import data set looks like"<<endl;
    ws->Print();

    Int_t dsentries = data->numEntries();
    cout<<"No. of entries in dataset = "<<dsentries<<endl;

    cout<<"Writing output tree with unbinned sweighted data"<<endl;
    TFile *fileout = new TFile("jpsilambda_LL_withsw.root","RECREATE");
    TTree *treeout = new TTree("MyTuple","SWeighted Stuff");

    Double_t pi_PT, B_MASS, SIGW, BKGW;

    treeout->Branch("B_MASS",&B_MASS,"B_MASS/D");
    treeout->Branch("pi_PT",&pi_PT,"pi_PT/D");

    treeout->Branch("SIGW",&SIGW,"SIGW/D");
    treeout->Branch("BKGW",&BKGW,"BKGW/D");

    for(Int_t i=0; i<dsentries; i++)
    {
        B_MASS = (data->get(i))->getRealValue("Lb_DTF_M_JpsiLConstr");
        pi_PT = (data->get(i))->getRealValue("pi_PT");

        SIGW = (data->get(i))->getRealValue("sigYield_sw");
        BKGW = (data->get(i))->getRealValue("bkgYield_sw");

        treeout->Fill();
    }
    fileout->Write();
    fileout->Close();
    cout<<"Done Writing output tree with unbinned sweighted data"<<endl;
}

void MakePlots(RooWorkspace* ws){

    // Here we make plots of the discriminating variable (Lb_DTF_M_JpsiLConstr) after the fit
    // and of the control variable (isolation) after unfolding with sPlot.
    cout << "make plots" << endl;

    // make our canvas
    TCanvas* cdata = new TCanvas("sPlot","sPlot demo", 1200, 800);

    ///////////
    TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
    TPad *pad2 = new TPad("pad2","pad2",0.0,0.0,1.0,0.3);
    pad1->SetBottomMargin(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.5);
    pad2->SetBorderMode(0);
    pad1->SetBorderMode(0);
    cdata->SetBorderMode(0);
    pad2->Draw();
    pad1->Draw();
    pad1->cd();

    gPad->SetTopMargin(0.06);
    pad1->Update();
    ////////////

    // get what we need out of the workspace
    RooAbsPdf* model = ws->pdf("model");
    RooAbsPdf* sig1 = ws->pdf("sig1");
    RooAbsPdf* sig2 = ws->pdf("sig2");
    RooAbsPdf* sig = ws->pdf("sig");
    RooAbsPdf* bkg = ws->pdf("bkg");

    RooRealVar* pi_PT = ws->var("pi_PT");
    RooRealVar* Lb_DTF_M_JpsiLConstr = ws->var("Lb_DTF_M_JpsiLConstr");

    cout<<"poop9"<<endl;
    //  RooFormulaVar sum_pt("sum_pt","@0+@1",RooArgSet(*pi_PT,*Pi_PT));

    // RooFormulaVar minipchi2("minipchi2","TMath::Min(@0,@1)",RooArgSet(*K_MINIPCHI2,*Pi_MINIPCHI2));
    ws->Print();
    // note, we get the dataset with sWeights
    RooDataSet* data = (RooDataSet*) ws->data("dataWithSWeights");

    // this shouldn't be necessary, need to fix something with workspace
    // do this to set parameters back to their fitted values.
    model->fitTo(*data);

    //plot Lb_DTF_M_JpsiLConstr for data with full model and individual componenets overlayed

    RooPlot* frame = Lb_DTF_M_JpsiLConstr->frame() ;
    data->plotOn(frame,Name("Hist"),DataError(RooAbsData::SumW2) ) ;
    model->plotOn(frame,Name("curvetot") ) ;
    model->plotOn(frame,Components(*sig1),LineStyle(kDotted),LineColor(kMagenta));
    model->plotOn(frame,Components(*sig2),LineStyle(kDotted),LineColor(kMagenta));
    model->plotOn(frame,Components(*sig),LineStyle(kDashed), LineColor(kRed)) ;
    model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen)) ;

    frame->SetTitle("Fit of model to discriminating variable");
    frame->Draw() ;

    ////Fit chi2

    Double_t chiSquare1 = frame->chiSquare("curvetot","Hist");
    cout<<"chi square1 = "<<chiSquare1<<endl;
    RooArgSet *floatpar = model->getParameters(data);
    floatpar->Print();
    int floatpars = (floatpar->selectByAttrib("Constant",kFALSE))->getSize();
    Double_t chi2 = frame->chiSquare("curvetot","Hist",floatpars);
    cout<<"chi square2 = "<<chi2<<endl;


    // Pull distribution
    RooPlot *framex2 = Lb_DTF_M_JpsiLConstr->frame();
    RooHist* hpull = frame->pullHist("Hist","curvetot");
    framex2->addPlotable(hpull,"P");
    hpull->SetLineColor(kBlack);
    hpull->SetMarkerColor(kBlack);
    framex2->SetTitle(0);
    framex2->GetYaxis()->SetTitle("Pull");
    framex2->GetYaxis()->SetTitleSize(0.15);
    framex2->GetYaxis()->SetLabelSize(0.15);
    framex2->GetXaxis()->SetTitleSize(0.2);
    framex2->GetXaxis()->SetLabelSize(0.15);
    framex2->GetYaxis()->CenterTitle();
    framex2->GetYaxis()->SetTitleOffset(0.45);
    framex2->GetXaxis()->SetTitleOffset(1.1);
    framex2->GetYaxis()->SetNdivisions(505);
    framex2->GetYaxis()->SetRangeUser(-8.8,8.8);
    pad2->cd();
    framex2->Draw();

    cdata->cd();

    cout<<"Pull Mean Y = "<<hpull->GetMean(2)<<endl;
    cout<<"Pull RMS Y = "<<hpull->GetRMS(2)<<endl;


    // Now use the sWeights to show isolation distribution for Z and QCD.
    // The SPlot class can make this easier, but here we demonstrait in more
    // detail how the sWeights are used.  The SPlot class should make this
    // very easy and needs some more development.

    // Plot isolation for Z component.
    // Do this by plotting all events weighted by the sWeight for the Z component.
    // The SPlot class adds a new variable that has the name of the corresponding
    // yield + "_sw".

    TCanvas *mycan1 = new TCanvas();
    RooRealVar* sigYield_sw = ws->var("sigYield_sw");
    sigYield_sw->setRange(-1.0,1.5);
    RooPlot* frame_sigsw = sigYield_sw->frame();
    frame_sigsw->SetTitle("Signal sweights");

    // create weightfed data set
    RooDataSet * dataw_z = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"sigYield_sw") ;

    RooDataSet * dataw_qcd = new RooDataSet(data->GetName(),data->GetTitle(),data,*data->get(),0,"bkgYield_sw") ;
    data->plotOn(frame_sigsw) ;
    frame_sigsw->Draw();

    // cdata->cd(3);
    // RooPlot* frame_sigsw_x = Lb_DTF_M_JpsiLConstr->frame();
    // frame_sigsw_x->SetTitle("sWeights vs bmass");
    // data->plotOnXY(frame_sigsw_x,YVar(*sigYield_sw),MarkerColor(kBlue));
    // frame_sigsw_x->Draw();

    TCanvas *can2 = new TCanvas();
    can2->Divide(1,2);

    can2->cd(1);

    RooPlot* frame2_1 = pi_PT->frame() ;

    dataw_z->plotOn(frame2_1, DataError(RooAbsData::SumW2) ) ;

    frame2_1->SetTitle("pi_PT distribution for Sig");
    frame2_1->Draw() ;

    can2->cd(2);

    RooPlot* frame2_2 = pi_PT->frame() ;
    dataw_qcd->plotOn(frame2_2,DataError(RooAbsData::SumW2) ) ;

    frame2_2->SetTitle("pi_PT distribution for Bkg");
    frame2_2->Draw() ;

    // TCanvas *can3 = new TCanvas();
    // can3->Divide(1,2);

    // can3->cd(1);
    // RooRealVar* sumpt1 = (RooRealVar*) dataw_z->addColumn(sum_pt);
    // sumpt1->setRange(0.0,10000.0);
    // RooPlot* frame3_1 = sumpt1->frame() ;
    // dataw_z->plotOn(frame3_1, DataError(RooAbsData::SumW2) ) ;
    // frame3_1->SetTitle("sum_PT distribution for Sig");
    // frame3_1->Draw();
    // can3->cd(2);
    // RooRealVar* sumpt2 = (RooRealVar*) dataw_qcd->addColumn(sum_pt);
    // sumpt2->setRange(0.0,10000.0);
    // RooPlot* frame3_2 = sumpt2->frame() ;
    // dataw_qcd->plotOn(frame3_2,DataError(RooAbsData::SumW2) ) ;
    // frame3_2->SetTitle("sum_PT distribution for Bkg");
    // frame3_2->Draw();

    // TCanvas *can4 = new TCanvas();
    // can4->cd();
    // can4->Divide(1,2);

    // can4->cd(1);
    // RooRealVar* minipchi2_1 = (RooRealVar*) dataw_z->addColumn(minipchi2);
    // minipchi2_1->setRange(0.0,200.0);
    // RooPlot* frame4_1 = minipchi2_1->frame() ;
    // dataw_z->plotOn(frame4_1, DataError(RooAbsData::SumW2) ) ;
    // frame4_1->SetTitle("MINIPCHI2 distribution for Sig");
    // frame4_1->Draw();

    // can4->cd(2);
    // RooRealVar* minipchi2_2 = (RooRealVar*) dataw_qcd->addColumn(minipchi2);
    // minipchi2_2->setRange(0.0,200.0);
    // RooPlot* frame4_2 = minipchi2_2->frame() ;
    // dataw_qcd->plotOn(frame4_2,DataError(RooAbsData::SumW2) ) ;
    // frame4_2->SetTitle("MINIPCHI2 distribution for Bkg");
    // frame4_2->Draw();

    // TCanvas *can5 = new TCanvas();
    // can5->cd();
    // can5->Divide(1,2);

    // can5->cd(1);
    // RooPlot* frame5_1 = B_ENDVERTEX_CHI2->frame() ;
    // dataw_z->plotOn(frame5_1, DataError(RooAbsData::SumW2) ) ;
    // frame5_1->SetTitle("B_ENDVERTEX_CHI2 distribution for Sig");
    // frame5_1->Draw();
    // can5->cd(2);
    // RooPlot* frame5_2 = B_ENDVERTEX_CHI2->frame() ;
    // dataw_qcd->plotOn(frame5_2,DataError(RooAbsData::SumW2) ) ;
    // frame5_2->SetTitle("B_ENDVERTEX_CHI2 distribution for Bkg");
    // frame5_2->Draw();

    // TCanvas *can6 = new TCanvas();
    // can6->cd();
    // can6->Divide(1,2);

    // can6->cd(1);
    // RooPlot* frame6_1 = K_IPCHI2_ORIVX->frame() ;
    // dataw_z->plotOn(frame6_1, DataError(RooAbsData::SumW2) ) ;
    // frame6_1->SetTitle("K_IPCHI2_ORIVX distribution for Sig");
    // frame6_1->Draw();
    // can6->cd(2);
    // RooPlot* frame6_2 = K_IPCHI2_ORIVX->frame() ;
    // dataw_qcd->plotOn(frame6_2,DataError(RooAbsData::SumW2) ) ;
    // frame6_2->SetTitle("K_IPCHI2_ORIVX distribution for Bkg");
    // frame6_2->Draw();

    // TFile *f1 = new TFile("save.root","RECREATE");
    // f1->cd();
    // TH1D *sigminipchi2 = (TH1D*)dataw_z->createHistogram("minipchi2",100);
    // TH1D *sigvchi2 = (TH1D*)dataw_z->createHistogram("B_ENDVERTEX_CHI2",100);
    // TH1D *sigkipchi2orivx = (TH1D*)dataw_z->createHistogram("K_IPCHI2_ORIVX",100);
    // TH1D *sigkpt = (TH1D*)dataw_z->createHistogram("pi_PT",200);
    // TH1D *sigsumpt = (TH1D*)dataw_z->createHistogram("sum_pt",100);

    // sigminipchi2->SetName("sigminipchi2");
    // sigvchi2->SetName("sigvchi2");
    // sigfd->SetName("sigfd");
    // sigdira->SetName("sigdira");
    // sigpt->SetName("sigpt");
    // sigsumpt->SetName("sigsumpt");
    // sigkpidk->SetName("sigkpidk");
    // sigminipchi2->Write();
    // sigvchi2->Write();
    // sigfd->Write();
    // sigdira->Write();
    // sigpt->Write();
    // sigsumpt->Write();
    // sigkpidk->Write();
    // f1->Close();
    //cdata->SaveAs("SPlot.gif");

}
