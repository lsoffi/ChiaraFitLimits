{
//=========Macro generated from canvas: c/c
//=========  (Fri Aug 23 18:36:14 2013) by ROOT version5.32/00
   TCanvas *c = new TCanvas("c", "c",0,22,700,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c->Range(1.954125,-6.728395,2.747377,5.617284);
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetLogx();
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetLeftMargin(0.18);
   c->SetRightMargin(0.15);
   c->SetTopMargin(0.05);
   c->SetBottomMargin(0.14);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5);
   grae->SetName("Graph");
   grae->SetTitle("Graph");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#00ff00");
   grae->SetFillColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,150,0.03968254);
   grae->SetPointError(0,0,0,1,1);
   grae->SetPoint(1,200,0.03968254);
   grae->SetPointError(1,0,0,1,1);
   grae->SetPoint(2,250,0.03968254);
   grae->SetPointError(2,0,0,1,1);
   grae->SetPoint(3,300,0.03968254);
   grae->SetPointError(3,0,0,1,1);
   grae->SetPoint(4,400,0.03968254);
   grae->SetPointError(4,0,0,1,1);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,125,425);
   Graph_Graph1->SetMinimum(-5);
   Graph_Graph1->SetMaximum(5);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineStyle(0);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetTitle("m_{H}");
   Graph_Graph1->GetXaxis()->SetNdivisions(505);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("median((N_{true}^{FWHM}-N_{fit}^{FWHM})/#Delta N_{fit}^{FWHM})");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.3);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1);
   
   grae->Draw("a3");
   
   grae = new TGraphAsymmErrors(5);
   grae->SetName("Graph");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,150,0.03968254);
   grae->SetPointError(0,0,0,0,0);
   grae->SetPoint(1,200,0.03968254);
   grae->SetPointError(1,0,0,0,0);
   grae->SetPoint(2,250,0.03968254);
   grae->SetPointError(2,0,0,0,0);
   grae->SetPoint(3,300,0.03968254);
   grae->SetPointError(3,0,0,0,0);
   grae->SetPoint(4,400,0.03968254);
   grae->SetPointError(4,0,0,0,0);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,125,425);
   Graph_Graph2->SetMinimum(0);
   Graph_Graph2->SetMaximum(1.139683);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   Graph_Graph2->SetLineStyle(0);
   Graph_Graph2->SetMarkerStyle(20);
   Graph_Graph2->GetXaxis()->SetNdivisions(505);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetYaxis()->SetTitleOffset(1.3);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph2->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph2);
   
   grae->Draw("3");
   
   TGraphErrors *gre = new TGraphErrors(5);
   gre->SetName("Graph");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);

   ci = TColor::GetColor("#0000ff");
   gre->SetLineColor(ci);
   gre->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   gre->SetMarkerColor(ci);
   gre->SetMarkerStyle(20);
   gre->SetPoint(0,150,0.03968254);
   gre->SetPointError(0,0,0);
   gre->SetPoint(1,200,0.03968254);
   gre->SetPointError(1,0,0);
   gre->SetPoint(2,250,0.03968254);
   gre->SetPointError(2,0,0);
   gre->SetPoint(3,300,0.03968254);
   gre->SetPointError(3,0,0);
   gre->SetPoint(4,400,0.03968254);
   gre->SetPointError(4,0,0);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,125,425);
   Graph_Graph1->SetMinimum(0);
   Graph_Graph1->SetMaximum(1.139683);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineStyle(0);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetNdivisions(505);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.08);
   Graph_Graph1->GetXaxis()->SetTitleOffset(0.8);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("median((N_{true}^{FWHM}-N_{fit}^{FWHM})/#Delta N_{fit}^{FWHM})");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.08);
   Graph_Graph1->GetYaxis()->SetTitleOffset(0.5);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelOffset(0.007);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetZaxis()->SetTitleSize(0.06);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_Graph1);
   
   gre->Draw("pxl ");
   
   TLegend *leg = new TLegend(0.5890805,0.7377622,0.8074713,0.9143357,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("Graph","95%","f");

   ci = TColor::GetColor("#00ff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph","68%","f");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry=leg->AddEntry("Graph","mean","pl");

   ci = TColor::GetColor("#0000ff");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);

   ci = TColor::GetColor("#0000ff");
   entry->SetMarkerColor(ci);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1);
   leg->Draw();
   
   TPaveText *pt = new TPaveText(0.2298851,0.7744755,0.3994253,0.9248252,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextFont(42);
   TText *text = pt->AddText("Gen fun: 2exp");
   text = pt->AddText("Fit fun: 2exp");
   text = pt->AddText("Category: 0");
   pt->Draw();
   
   pt = new TPaveText(0.175,0.953,0.6,0.975,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextAlign(11);
   pt->SetTextFont(42);
   pt->SetTextSize(0.038);
   text = pt->AddText("CMS Simulation");
   pt->Draw();
   
   pt = new TPaveText(0.65,0.953,0.87,0.975,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextAlign(31);
   pt->SetTextFont(42);
   pt->SetTextSize(0.038);
   text = pt->AddText("#sqrt{s} = 8 TeV");
   pt->Draw();
   TLine *line = new TLine(0.1,0.5,1,0.5);
   line->Draw();
   line = new TLine(125.2364,0.005611675,422.11,0.005611675);
   line->SetLineStyle(7);
   line->SetLineWidth(2);
   line->Draw();
   c->Modified();
   c->cd();
   c->SetSelected(c);
}
