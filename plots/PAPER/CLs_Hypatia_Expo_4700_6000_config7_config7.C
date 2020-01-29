void CLs_Hypatia_Expo_4700_6000_config7_config7()
{
//=========Macro generated from canvas: Asymptotic_Scan/
//=========  (Mon Nov 18 09:43:40 2019) by ROOT version 6.12/06
   TCanvas *Asymptotic_Scan = new TCanvas("Asymptotic_Scan", "",10,32,700,520);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   Asymptotic_Scan->ToggleEventStatus();
   Asymptotic_Scan->Range(-0.8650617,-0.2226187,5.313951,1.169587);
   Asymptotic_Scan->SetFillColor(0);
   Asymptotic_Scan->SetBorderMode(0);
   Asymptotic_Scan->SetBorderSize(2);
   Asymptotic_Scan->SetTickx(1);
   Asymptotic_Scan->SetTicky(1);
   Asymptotic_Scan->SetLeftMargin(0.14);
   Asymptotic_Scan->SetRightMargin(0.05);
   Asymptotic_Scan->SetTopMargin(0.05);
   Asymptotic_Scan->SetBottomMargin(0.16);
   Asymptotic_Scan->SetFrameLineWidth(2);
   Asymptotic_Scan->SetFrameBorderMode(0);
   Asymptotic_Scan->SetFrameLineWidth(2);
   Asymptotic_Scan->SetFrameBorderMode(0);
   
   Double_t CLs_observed_fx1001[20] = {
   0,
   0.2631579,
   0.5263158,
   0.7894737,
   1.052632,
   1.315789,
   1.578947,
   1.842105,
   2.105263,
   2.368421,
   2.631579,
   2.894737,
   3.157895,
   3.421053,
   3.684211,
   3.947368,
   4.210526,
   4.473684,
   4.736842,
   5};
   Double_t CLs_observed_fy1001[20] = {
   0.9999925,
   0.7520516,
   0.5569476,
   0.4038145,
   0.2881196,
   0.2048232,
   0.1411612,
   0.09526804,
   0.06281666,
   0.0403935,
   0.02543887,
   0.01566569,
   0.009490475,
   0.005598148,
   0.00323037,
   0.001824035,
   0.001008962,
   0.000545652,
   0.0002888652,
   0.0001491543};
   Double_t CLs_observed_fex1001[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t CLs_observed_fey1001[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   TGraphErrors *gre = new TGraphErrors(20,CLs_observed_fx1001,CLs_observed_fy1001,CLs_observed_fex1001,CLs_observed_fey1001);
   gre->SetName("CLs_observed");
   gre->SetTitle("Observed CLs");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_CLs_observed1001 = new TH1F("Graph_CLs_observed1001","",100,0,5.5);
   Graph_CLs_observed1001->SetMinimum(0.0001342389);
   Graph_CLs_observed1001->SetMaximum(1.099977);
   Graph_CLs_observed1001->SetDirectory(0);
   Graph_CLs_observed1001->SetStats(0);
   Graph_CLs_observed1001->SetLineWidth(2);
   Graph_CLs_observed1001->SetMarkerStyle(20);
   Graph_CLs_observed1001->GetXaxis()->SetTitle("R#upoint 10^{3}");
   Graph_CLs_observed1001->GetXaxis()->SetRange(1,91);
   Graph_CLs_observed1001->GetXaxis()->SetNdivisions(505);
   Graph_CLs_observed1001->GetXaxis()->SetLabelFont(132);
   Graph_CLs_observed1001->GetXaxis()->SetLabelOffset(0.01);
   Graph_CLs_observed1001->GetXaxis()->SetLabelSize(0.06);
   Graph_CLs_observed1001->GetXaxis()->SetTitleSize(0.072);
   Graph_CLs_observed1001->GetXaxis()->SetTitleOffset(0.95);
   Graph_CLs_observed1001->GetXaxis()->SetTitleFont(132);
   Graph_CLs_observed1001->GetYaxis()->SetTitle("p-value");
   Graph_CLs_observed1001->GetYaxis()->SetLabelFont(132);
   Graph_CLs_observed1001->GetYaxis()->SetLabelOffset(0.01);
   Graph_CLs_observed1001->GetYaxis()->SetLabelSize(0.06);
   Graph_CLs_observed1001->GetYaxis()->SetTitleSize(0.072);
   Graph_CLs_observed1001->GetYaxis()->SetTitleOffset(0.95);
   Graph_CLs_observed1001->GetYaxis()->SetTitleFont(132);
   Graph_CLs_observed1001->GetZaxis()->SetLabelFont(132);
   Graph_CLs_observed1001->GetZaxis()->SetLabelSize(0.06);
   Graph_CLs_observed1001->GetZaxis()->SetTitleSize(0.072);
   Graph_CLs_observed1001->GetZaxis()->SetTitleOffset(1.2);
   Graph_CLs_observed1001->GetZaxis()->SetTitleFont(132);
   gre->SetHistogram(Graph_CLs_observed1001);
   
   gre->Draw("apl");
   
   TMultiGraph *multigraph = new TMultiGraph();
   multigraph->SetName("HTI_Result_Plot_expected");
   multigraph->SetTitle("Expected ");
   
   Double_t _fx3001[20] = {
   0,
   0.2631579,
   0.5263158,
   0.7894737,
   1.052632,
   1.315789,
   1.578947,
   1.842105,
   2.105263,
   2.368421,
   2.631579,
   2.894737,
   3.157895,
   3.421053,
   3.684211,
   3.947368,
   4.210526,
   4.473684,
   4.736842,
   5};
   Double_t _fy3001[20] = {
   0.9999925,
   0.8722031,
   0.7505977,
   0.6341695,
   0.5297403,
   0.438006,
   0.3533042,
   0.2795796,
   0.2163853,
   0.1634808,
   0.1216473,
   0.08834271,
   0.06252871,
   0.04339713,
   0.02946654,
   0.01956435,
   0.01277483,
   0.008124584,
   0.005055586,
   0.003077424};
   Double_t _felx3001[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3001[20] = {
   1.471978e-05,
   0.1973477,
   0.3009871,
   0.3420541,
   0.3411987,
   0.3168948,
   0.2783747,
   0.234278,
   0.1897627,
   0.1483126,
   0.1131128,
   0.0836822,
   0.06006236,
   0.04211734,
   0.02881789,
   0.01924342,
   0.01261828,
   0.008050509,
   0.005021313,
   0.003061923};
   Double_t _fehx3001[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3001[20] = {
   6.939238e-06,
   0.1173623,
   0.2253349,
   0.3238887,
   0.4064641,
   0.4723098,
   0.5247673,
   0.5604135,
   0.5792225,
   0.581387,
   0.5684556,
   0.5423113,
   0.5050551,
   0.4601328,
   0.4102166,
   0.3579887,
   0.3065622,
   0.2567669,
   0.2107044,
   0.1694001};
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(20,_fx3001,_fy3001,_felx3001,_fehx3001,_fely3001,_fehy3001);
   grae->SetName("");
   grae->SetTitle("Expected CLs #pm 2 #sigma");

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph3001 = new TH1F("Graph_Graph3001","Expected CLs #pm 2 #sigma",100,0,5.5);
   Graph_Graph3001->SetMinimum(1.395106e-05);
   Graph_Graph3001->SetMaximum(1.099998);
   Graph_Graph3001->SetDirectory(0);
   Graph_Graph3001->SetStats(0);
   Graph_Graph3001->SetLineWidth(2);
   Graph_Graph3001->SetMarkerStyle(20);
   Graph_Graph3001->GetXaxis()->SetNdivisions(505);
   Graph_Graph3001->GetXaxis()->SetLabelFont(132);
   Graph_Graph3001->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph3001->GetXaxis()->SetLabelSize(0.06);
   Graph_Graph3001->GetXaxis()->SetTitleSize(0.072);
   Graph_Graph3001->GetXaxis()->SetTitleOffset(0.95);
   Graph_Graph3001->GetXaxis()->SetTitleFont(132);
   Graph_Graph3001->GetYaxis()->SetLabelFont(132);
   Graph_Graph3001->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph3001->GetYaxis()->SetLabelSize(0.06);
   Graph_Graph3001->GetYaxis()->SetTitleSize(0.072);
   Graph_Graph3001->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph3001->GetYaxis()->SetTitleFont(132);
   Graph_Graph3001->GetZaxis()->SetLabelFont(132);
   Graph_Graph3001->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph3001->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph3001->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph3001->GetZaxis()->SetTitleFont(132);
   grae->SetHistogram(Graph_Graph3001);
   
   multigraph->Add(grae,"3");
   
   Double_t _fx3002[20] = {
   0,
   0.2631579,
   0.5263158,
   0.7894737,
   1.052632,
   1.315789,
   1.578947,
   1.842105,
   2.105263,
   2.368421,
   2.631579,
   2.894737,
   3.157895,
   3.421053,
   3.684211,
   3.947368,
   4.210526,
   4.473684,
   4.736842,
   5};
   Double_t _fy3002[20] = {
   0.9999925,
   0.8722031,
   0.7505977,
   0.6341695,
   0.5297403,
   0.438006,
   0.3533042,
   0.2795796,
   0.2163853,
   0.1634808,
   0.1216473,
   0.08834271,
   0.06252871,
   0.04339713,
   0.02946654,
   0.01956435,
   0.01277483,
   0.008124584,
   0.005055586,
   0.003077424};
   Double_t _felx3002[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fely3002[20] = {
   6.795409e-06,
   0.09788731,
   0.1595255,
   0.1930251,
   0.2037518,
   0.1991085,
   0.1836586,
   0.1616833,
   0.1365302,
   0.1108815,
   0.08749315,
   0.06676513,
   0.04928376,
   0.03542748,
   0.02477994,
   0.01687228,
   0.01125251,
   0.007287938,
   0.004605922,
   0.002841119};
   Double_t _fehx3002[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t _fehy3002[20] = {
   4.768118e-06,
   0.07782891,
   0.1437177,
   0.1977215,
   0.2367751,
   0.2618149,
   0.2749992,
   0.2762162,
   0.266936,
   0.2489422,
   0.2252754,
   0.1976904,
   0.1682905,
   0.1395191,
   0.1126396,
   0.08860027,
   0.06818367,
   0.05106,
   0.03732559,
   0.0266416};
   grae = new TGraphAsymmErrors(20,_fx3002,_fy3002,_felx3002,_fehx3002,_fely3002,_fehy3002);
   grae->SetName("");
   grae->SetTitle("Expected CLs #pm 1 #sigma");

   ci = TColor::GetColor("#00ff00");
   grae->SetFillColor(ci);
   grae->SetFillStyle(1000);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   
   TH1F *Graph_Graph3002 = new TH1F("Graph_Graph3002","Expected CLs #pm 1 #sigma",100,0,5.5);
   Graph_Graph3002->SetMinimum(0.0002126751);
   Graph_Graph3002->SetMaximum(1.099973);
   Graph_Graph3002->SetDirectory(0);
   Graph_Graph3002->SetStats(0);
   Graph_Graph3002->SetLineWidth(2);
   Graph_Graph3002->SetMarkerStyle(20);
   Graph_Graph3002->GetXaxis()->SetNdivisions(505);
   Graph_Graph3002->GetXaxis()->SetLabelFont(132);
   Graph_Graph3002->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph3002->GetXaxis()->SetLabelSize(0.06);
   Graph_Graph3002->GetXaxis()->SetTitleSize(0.072);
   Graph_Graph3002->GetXaxis()->SetTitleOffset(0.95);
   Graph_Graph3002->GetXaxis()->SetTitleFont(132);
   Graph_Graph3002->GetYaxis()->SetLabelFont(132);
   Graph_Graph3002->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph3002->GetYaxis()->SetLabelSize(0.06);
   Graph_Graph3002->GetYaxis()->SetTitleSize(0.072);
   Graph_Graph3002->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph3002->GetYaxis()->SetTitleFont(132);
   Graph_Graph3002->GetZaxis()->SetLabelFont(132);
   Graph_Graph3002->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph3002->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph3002->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph3002->GetZaxis()->SetTitleFont(132);
   grae->SetHistogram(Graph_Graph3002);
   
   multigraph->Add(grae,"3");
   
   Double_t _fx1[20] = {
   0,
   0.2631579,
   0.5263158,
   0.7894737,
   1.052632,
   1.315789,
   1.578947,
   1.842105,
   2.105263,
   2.368421,
   2.631579,
   2.894737,
   3.157895,
   3.421053,
   3.684211,
   3.947368,
   4.210526,
   4.473684,
   4.736842,
   5};
   Double_t _fy1[20] = {
   0.9999925,
   0.8722031,
   0.7505977,
   0.6341695,
   0.5297403,
   0.438006,
   0.3533042,
   0.2795796,
   0.2163853,
   0.1634808,
   0.1216473,
   0.08834271,
   0.06252871,
   0.04339713,
   0.02946654,
   0.01956435,
   0.01277483,
   0.008124584,
   0.005055586,
   0.003077424};
   TGraph *graph = new TGraph(20,_fx1,_fy1);
   graph->SetName("");
   graph->SetTitle("Expected CLs - Median");
   graph->SetFillStyle(1000);
   graph->SetLineStyle(2);
   graph->SetLineWidth(2);
   graph->SetMarkerStyle(20);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Expected CLs - Median",100,0,5.5);
   Graph_Graph1->SetMinimum(0.002769682);
   Graph_Graph1->SetMaximum(1.099684);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineWidth(2);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetNdivisions(505);
   Graph_Graph1->GetXaxis()->SetLabelFont(132);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.06);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.072);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.95);
   Graph_Graph1->GetXaxis()->SetTitleFont(132);
   Graph_Graph1->GetYaxis()->SetLabelFont(132);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.06);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.072);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph1->GetYaxis()->SetTitleFont(132);
   Graph_Graph1->GetZaxis()->SetLabelFont(132);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph1->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph1->GetZaxis()->SetTitleFont(132);
   graph->SetHistogram(Graph_Graph1);
   
   multigraph->Add(graph,"L");
   multigraph->Draw("");
   TLine *line = new TLine(0,0.05,5,0.05);

   ci = TColor::GetColor("#ff0000");
   line->SetLineColor(ci);
   line->SetLineWidth(2);
   line->Draw();
   
   Double_t CLs_observed_fx1002[20] = {
   0,
   0.2631579,
   0.5263158,
   0.7894737,
   1.052632,
   1.315789,
   1.578947,
   1.842105,
   2.105263,
   2.368421,
   2.631579,
   2.894737,
   3.157895,
   3.421053,
   3.684211,
   3.947368,
   4.210526,
   4.473684,
   4.736842,
   5};
   Double_t CLs_observed_fy1002[20] = {
   0.9999925,
   0.7520516,
   0.5569476,
   0.4038145,
   0.2881196,
   0.2048232,
   0.1411612,
   0.09526804,
   0.06281666,
   0.0403935,
   0.02543887,
   0.01566569,
   0.009490475,
   0.005598148,
   0.00323037,
   0.001824035,
   0.001008962,
   0.000545652,
   0.0002888652,
   0.0001491543};
   Double_t CLs_observed_fex1002[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   Double_t CLs_observed_fey1002[20] = {
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0,
   0};
   gre = new TGraphErrors(20,CLs_observed_fx1002,CLs_observed_fy1002,CLs_observed_fex1002,CLs_observed_fey1002);
   gre->SetName("CLs_observed");
   gre->SetTitle("Observed CLs");
   gre->SetFillStyle(1000);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   
   TH1F *Graph_Graph_CLs_observed10011002 = new TH1F("Graph_Graph_CLs_observed10011002","",100,0,5.5);
   Graph_Graph_CLs_observed10011002->SetMinimum(0.0001342389);
   Graph_Graph_CLs_observed10011002->SetMaximum(1.099977);
   Graph_Graph_CLs_observed10011002->SetDirectory(0);
   Graph_Graph_CLs_observed10011002->SetStats(0);
   Graph_Graph_CLs_observed10011002->SetLineWidth(2);
   Graph_Graph_CLs_observed10011002->SetMarkerStyle(20);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetTitle("R#upoint 10^{3}");
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetRange(1,91);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetNdivisions(505);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetLabelFont(132);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetLabelSize(0.06);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetTitleSize(0.072);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetTitleOffset(0.95);
   Graph_Graph_CLs_observed10011002->GetXaxis()->SetTitleFont(132);
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetTitle("p-value");
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetLabelFont(132);
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetLabelSize(0.06);
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetTitleSize(0.072);
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetTitleOffset(0.95);
   Graph_Graph_CLs_observed10011002->GetYaxis()->SetTitleFont(132);
   Graph_Graph_CLs_observed10011002->GetZaxis()->SetLabelFont(132);
   Graph_Graph_CLs_observed10011002->GetZaxis()->SetLabelSize(0.06);
   Graph_Graph_CLs_observed10011002->GetZaxis()->SetTitleSize(0.072);
   Graph_Graph_CLs_observed10011002->GetZaxis()->SetTitleOffset(1.2);
   Graph_Graph_CLs_observed10011002->GetZaxis()->SetTitleFont(132);
   gre->SetHistogram(Graph_Graph_CLs_observed10011002);
   
   gre->Draw("pl");
   
   TLegend *leg = new TLegend(0.6217765,0.5991561,0.9226361,0.9177215,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(132);
   leg->SetTextSize(0.04008439);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(2);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("CLs_observed","Observed CLs","PEL");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("","Expected CLs - Median","L");
   entry->SetLineColor(1);
   entry->SetLineStyle(2);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("","Expected CLs #pm 1 #sigma","F");

   ci = TColor::GetColor("#00ff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   entry=leg->AddEntry("","Expected CLs #pm 2 #sigma","F");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1000);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(132);
   leg->Draw();
   
   TH1F *CLs_observed_copy__1 = new TH1F("CLs_observed_copy__1","",100,0,5.5);
   CLs_observed_copy__1->SetMinimum(0.0001342389);
   CLs_observed_copy__1->SetMaximum(1.099977);
   CLs_observed_copy__1->SetDirectory(0);
   CLs_observed_copy__1->SetStats(0);
   CLs_observed_copy__1->SetLineWidth(2);
   CLs_observed_copy__1->SetMarkerStyle(20);
   CLs_observed_copy__1->GetXaxis()->SetTitle("R");
   CLs_observed_copy__1->GetXaxis()->SetNdivisions(505);
   CLs_observed_copy__1->GetXaxis()->SetLabelFont(132);
   CLs_observed_copy__1->GetXaxis()->SetLabelOffset(0.01);
   CLs_observed_copy__1->GetXaxis()->SetLabelSize(0.06);
   CLs_observed_copy__1->GetXaxis()->SetTitleSize(0.072);
   CLs_observed_copy__1->GetXaxis()->SetTitleOffset(0.95);
   CLs_observed_copy__1->GetXaxis()->SetTitleFont(132);
   CLs_observed_copy__1->GetYaxis()->SetTitle("p value");
   CLs_observed_copy__1->GetYaxis()->SetLabelFont(132);
   CLs_observed_copy__1->GetYaxis()->SetLabelOffset(0.01);
   CLs_observed_copy__1->GetYaxis()->SetLabelSize(0.06);
   CLs_observed_copy__1->GetYaxis()->SetTitleSize(0.072);
   CLs_observed_copy__1->GetYaxis()->SetTitleOffset(0.95);
   CLs_observed_copy__1->GetYaxis()->SetTitleFont(132);
   CLs_observed_copy__1->GetZaxis()->SetLabelFont(132);
   CLs_observed_copy__1->GetZaxis()->SetLabelSize(0.06);
   CLs_observed_copy__1->GetZaxis()->SetTitleSize(0.072);
   CLs_observed_copy__1->GetZaxis()->SetTitleOffset(1.2);
   CLs_observed_copy__1->GetZaxis()->SetTitleFont(132);
   CLs_observed_copy__1->Draw("sameaxis");
   
   TPaveText *pt = new TPaveText(0.3796562,0.8122363,0.491404,0.9240506,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetLineWidth(2);
   pt->SetTextAlign(12);
   pt->SetTextFont(132);
   TText *pt_LaTex = pt->AddText("LHCb");
   pt->Draw();
   Asymptotic_Scan->Modified();
   Asymptotic_Scan->cd();
   Asymptotic_Scan->SetSelected(Asymptotic_Scan);
}
