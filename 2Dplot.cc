#include <iostream> 
#include <algorithm>
#include <vector>

#include <fstream>
#include <TROOT.h>

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphSmooth.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>
#include "TGraph2D.h"
#include <TH1F.h>
#include <TString.h>
#include <TCut.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TFile.h>
#include "TH2F.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataHist.h"
#include "RooHistFunc.h"


TPaveText* get_labelCMS( int legendQuadrant = 0 ,  bool sim=false) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for CMS label. Using 2." << std::endl;
    legendQuadrant = 2;
  }

  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.83;
    x2 = 0.8;
    y2 = 0.87;
  } else if( legendQuadrant==2 ) {
    x1 =  0.25;
    y1 = 0.83;
    x2 =  0.42;
    y2 = 0.87;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.2;
    x2 = 0.42;
    y2 = 0.24;
  } else if( legendQuadrant==0 ) {
    x1 = 0.175;
    y1 = 0.953;
    x2 = 0.6;
    y2 = 0.975;
  }

  
  TPaveText* cmslabel = new TPaveText( x1, y1, x2, y2, "brNDC" );
  cmslabel->SetFillColor(kWhite);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetBorderSize(0.);
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
 
    std::string leftText;
   
     
    if (sim)  leftText = "CMS Simulation"; //cwr ->remove 2011
    else {
     leftText = "CMS Preliminary, 19.5 fb^{-1}";
    }
    cmslabel->AddText(leftText.c_str());
    return cmslabel;

}




TPaveText* get_labelSqrt( int legendQuadrant ) {

  if( legendQuadrant!=0 && legendQuadrant!=1 && legendQuadrant!=2 && legendQuadrant!=3 ) {
    std::cout << "WARNING! Legend quadrant '" << legendQuadrant << "' not yet implemented for Sqrt label. Using 2." << std::endl;
    legendQuadrant = 2;
  }


  float x1, y1, x2, y2;
  if( legendQuadrant==1 ) {
    x1 = 0.63;
    y1 = 0.78;
    x2 = 0.8;
    y2 = 0.82;
  } else if( legendQuadrant==2 ) {
    x1 = 0.25;
    y1 = 0.78;
    x2 = 0.42;
    y2 = 0.82;
  } else if( legendQuadrant==3 ) {
    x1 = 0.25;
    y1 = 0.16;
    x2 = 0.42;
    y2 = 0.2;
  } else if( legendQuadrant==0 ) {
    x1 = 0.65;
    y1 = 0.953;
    x2 = 0.87;
    y2 = 0.975;
  }


  TPaveText* label_sqrt = new TPaveText(x1,y1,x2,y2, "brNDC");
  label_sqrt->SetFillColor(kWhite);
  label_sqrt->SetTextSize(0.038);
  label_sqrt->SetTextFont(42);
  label_sqrt->SetBorderSize(0.);
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->AddText("#sqrt{s} = 8 TeV");
  return label_sqrt;

}



double getMax(double* array){
  double temp = 0.;
  for(int i = 0;i<5; i++) if(fabs(array[i])> fabs(temp)) temp = array[i];
  return fabs(temp);
}



void make2Dplot(int cat){
 int NMAX=10;

 ifstream file;
 if(cat!=4) file.open(TString::Format("limitsResults_CH%d.txt",cat));
 if (cat==4)file.open("limitsResults_COMB.txt");
 if (cat==5)file.open("limitsResults_COMBOnlyEB.txt");
 int i = 0;
 TGraph2D* gr= new TGraph2D();
 std::vector<double> vec_limit;
 while (!file.eof()) {
    	  double mass;
	  double limit;
	  double up68;
	  double dn68;
	  double width;
	  file >> mass >> limit >> up68>>dn68>>width;
	  
	  if(width>1) width=width/100*mass;	  
	  std::cout<<mass<<" "<<limit<<" "<<width<<std::endl;
	  if(width>100)continue;
	  vec_limit.push_back(limit);
	  gr->SetPoint(i, mass, width, limit);
	  i++;
 }

 TCanvas* c = new TCanvas("c", "c", 1);
 c->cd();
 TPaveText* label_cms = get_labelCMS(0, false);
 TPaveText* label_sqrt = get_labelSqrt(0);  



  gr->Draw("COLZ");
  gr->GetXaxis()->SetTitle("m_{H} [GeV]");
  gr->GetYaxis()->SetTitle("#Gamma_{H} [GeV]");
  //gr->GetZaxis()->SetTitle("#sigma #times BR (H #rightarrow #gamma #gamma)_{95%CL} (pb)");
  gr->GetZaxis()->SetTitleSize(0.04);
  gr->GetZaxis()->SetTitleOffset(1.8);
  
  
  if(cat==4 || cat==5)gr->GetZaxis()->SetRangeUser(0.0007, 0.04);
  else gr->GetZaxis()->SetRangeUser(0.001, 0.08);
  
  gr->GetXaxis()->SetTitleSize(0.04);
  gr->GetXaxis()->SetTitleOffset(1.8);
  gr->GetYaxis()->SetTitleSize(0.04);
  gr->GetYaxis()->SetTitleOffset(1.8);

  double cont[5];
  /*if(cat==2)for(int i=0; i<5;i++)cont[i]=0.03+i*0.02;
  if( cat==3)for(int i=0; i<5;i++)cont[i]=0.03+i*0.02;
  if(cat==0 )for(int i=0; i<5;i++)cont[i]=0.01+i*0.01;
  if(cat==4 )for(int i=0; i<5;i++)cont[i]=0.005+i*0.005;
  if(cat==1 )for(int i=0; i<5;i++)cont[i]=0.01+i*0.01;*/
  if(cat==4 || cat==5){
  cont[0]= 0.001;
  cont[1] = 0.005;
  cont[2] = 0.008;
  cont[3] = 0.01;
  cont[4] = 0.02;
  }else {//if(cat==0
    cont[0]= 0.005;
    cont[1] = 0.01;
    cont[2] = 0.02;
    cont[3] = 0.04;
    cont[4] = 0.05;
 
  }
  TLegend* leg = new TLegend(0.2,0.6744755,0.494253,0.8948252,"", "brNDC");
  if(cat!=4)leg->SetHeader(TString::Format("Cat: %d",cat));
  if(cat==4) leg->SetHeader(TString::Format(" All Events",cat));
  if(cat==5) leg->SetHeader(TString::Format(" Events in EB",cat));
  leg->SetBorderSize(0.);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  
  TList* list1 = (TList*)gr->GetContourList(cont[0]);
  TGraph* gr1 = (TGraph*)list1->First();
  gr1->SetLineColor(kBlack);  
  gr1->SetLineStyle(1);
  gr1->SetLineWidth(2);
  gr1->Draw("same");
  leg->AddEntry(gr1, TString::Format("%.3f pb",cont[0]), "L");
   
  TList* list2 = (TList*)gr->GetContourList(cont[1]);
  TGraph* gr2 = (TGraph*)list2->First();
  gr2->SetLineColor(kBlack);
  gr2->SetLineStyle(3);
  gr2->SetLineWidth(2);
  gr2->Draw("same");
  leg->AddEntry(gr2, TString::Format("%.3f pb",cont[1]), "L");


  TList* list3 = (TList*)gr->GetContourList(cont[2]);
  list3->Print("V");
  TGraph* gr3 = (TGraph*)list3->First();
  gr3->SetLineColor(kBlack);
  gr3->SetLineStyle(6);
  gr3->SetLineWidth(2);
   gr3->Draw("same");
  leg->AddEntry(gr3,TString::Format("%.3f pb",cont[2]), "L");


  TList* list4 = (TList*)gr->GetContourList(cont[3]);
  list4->Print("V");
  TGraph* gr4 = (TGraph*)list4->First();
  gr4->SetLineColor(kBlack);
  gr4->SetLineStyle(7);
  gr4->SetLineWidth(2);
  gr4->Draw("same");
  leg->AddEntry(gr4,TString::Format("%.3f pb",cont[3]), "L");



  TList* list5 = (TList*)gr->GetContourList(cont[4]);
  list5->Print("V");
  TGraph* gr5 = (TGraph*)list5->First();
  gr5->SetLineColor(kBlack);
  gr5->SetLineStyle(10);
  gr5->SetLineWidth(2);
  gr5->Draw("same");
  leg->AddEntry(gr5,TString::Format("%.3f pb",cont[4]),"L");
  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  
  leg->Draw("same");
 // h2->RedrawAxis();
  if(cat<4){
    c->SaveAs(TString::Format("plots/2Dlimit_CH%d.png", cat));
    c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_CH%d.png", cat));
    c->SaveAs(TString::Format("plots/2Dlimit_CH%d.pdf", cat));
    c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_CH%d.pdf", cat));
    c->SetLogz();
    c->SaveAs(TString::Format("plots/2Dlimit_CH%d_LOG.png", cat));
    c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_CH%d_LOG.png", cat));
    c->SaveAs(TString::Format("plots/2Dlimit_CH%d_LOG.pdf", cat));
    c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_CH%d_LOG.pdf", cat));
  }else if(cat==4){
    c->SaveAs("plots/2Dlimit_COMB.png");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMB.png");
    c->SaveAs("plots/2Dlimit_COMB.pdf");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMB.pdf");
    c->SetLogz();
    c->SaveAs("plots/2Dlimit_COMB_LOG.png");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMB_LOG.png");
    c->SaveAs("plots/2Dlimit_COMB_LOG.pdf");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMB_LOG.pdf");

 } else if(cat==5){
    c->SaveAs("plots/2Dlimit_COMBOnlyEB.png");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMBOnlyEB.png");
    c->SaveAs("plots/2Dlimit_COMBOnlyEB.pdf");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMBOnlyEB.pdf");
    c->SetLogz();
    c->SaveAs("plots/2Dlimit_COMBOnlyEB_LOG.png");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMBOnlyEB_LOG.png");
    c->SaveAs("plots/2Dlimit_COMBOnlyEB_LOG.pdf");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/2Dlimit_COMBOnlyEB_LOG.pdf");

  }
}

 

void makeRatioBkgModels(int w){


 ifstream file01;
 file01.open(TString::Format("limitsResults_COMB_Expol-DiJet_w01.txt"));
 int i = 0;
 TGraph* gr01= new TGraph();
 ifstream file15;
 file15.open(TString::Format("limitsResults_COMB_Expol-DiJet_w15.txt"));
 int j = 0;
 TGraph* gr15= new TGraph();
 while (!file01.eof()) {
    	  double mass;
	  double limit_expol;
	  double limit_dijet;
	 
	  file01 >> mass >> limit_expol >> limit_dijet ;
	  std::cout<<mass<<" "<<limit_expol<<"  "<<limit_dijet<<"   "<<(limit_expol-limit_dijet)/limit_expol <<std::endl;
	  gr01->SetPoint(i, mass,(limit_expol-limit_dijet)/limit_expol );
	  i++;
 }
while (!file15.eof()) {
    	  double mass;
	  double limit_expol;
	  double limit_dijet;
	

	  file15 >> mass >> limit_expol >> limit_dijet;
	  std::cout<<mass<<" "<<limit_expol<<"  "<<limit_dijet<<"  "<<(limit_expol-limit_dijet)/limit_expol <<std::endl;
	  gr15->SetPoint(j, mass,(limit_expol-limit_dijet)/limit_expol );
	  j++;
 }

 TCanvas* c = new TCanvas("c", "c", 1);
 c->cd();
 TPaveText* label_cms = get_labelCMS(0, false);
 TPaveText* label_sqrt = get_labelSqrt(0);  
 TMultiGraph* mg = new TMultiGraph();
 mg->Add(gr01);
 mg->Add(gr15);
 gr15->SetMarkerColor(kRed);

  
  mg->Draw("AP");
 
  mg->GetXaxis()->SetTitle("m_{H} [GeV]");
  mg->GetYaxis()->SetTitle("#Delta #epsilon");
  //gr->GetZaxis()->SetTitle("#sigma #times BR (H #rightarrow #gamma #gamma)_{95%CL} (pb)");
 
  
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitleOffset(1.8);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleOffset(1.8);
 
  TLegend* leg = new TLegend(0.598851,0.5744755,0.794253,0.8248252,"", "brNDC");
  leg->SetHeader(TString::Format("Cat: All Events "));
  leg->SetBorderSize(0.);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  leg->AddEntry(gr01, "Width = 0.1 GeV", "P");
  leg->AddEntry(gr15, "Width = 15 GeV", "P");
 
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  leg->Draw("same");


    c->SaveAs("plots/limitsResults_COMB_Expol-DiJet_w%d.png");
    c->SaveAs("plots/limitsResults_COMB_Expol-DiJet_w%d.pdf");

}



 

void makeRatioEbEe(){


  ifstream fileeb;
  fileeb.open(TString::Format("limitsResults_COMBOnlyEB.txt"));
  ifstream file;
  file.open(TString::Format("limitsResults_COMB.txt"));
  TGraphAsymmErrors* gr01= new TGraphAsymmErrors();
  std::vector<double> lim01eb;
  std::vector<double> lim01; 
  std::vector<double> mass01;
  std::vector<double> lim2eb;
  std::vector<double> lim2;
  std::vector<double> mass2;
  std::vector<double> lim5eb;
  std::vector<double> lim5;
  std::vector<double> mass5;
  std::vector<double> lim7eb;
  std::vector<double> lim7;
  std::vector<double> mass7;
  std::vector<double> lim10eb;
  std::vector<double> lim10;
  std::vector<double> mass10;

  std::vector<double> lim01ebUP68;
  std::vector<double> lim01ebDN68;
  std::vector<double> lim01UP68;
  std::vector<double> lim01DN68;
  std::vector<double> lim2ebUP68;
  std::vector<double> lim2ebDN68;
  std::vector<double> lim2UP68;
  std::vector<double> lim2DN68;
  std::vector<double> lim5ebUP68;
  std::vector<double> lim5ebDN68;
  std::vector<double> lim5UP68;
  std::vector<double> lim5DN68;
  std::vector<double> lim7ebUP68;
  std::vector<double> lim7ebDN68;
  std::vector<double> lim7UP68;
  std::vector<double> lim7DN68;
  std::vector<double> lim10ebUP68;
  std::vector<double> lim10ebDN68;
  std::vector<double> lim10UP68;
  std::vector<double> lim10DN68;


  TGraphAsymmErrors* gr2= new TGraphAsymmErrors();  
  TGraphAsymmErrors* gr5= new TGraphAsymmErrors();  
  TGraphAsymmErrors* gr7= new TGraphAsymmErrors();
  TGraphAsymmErrors* gr10= new TGraphAsymmErrors();
  
  while (!fileeb.eof()) {
    	  double mass;
	  double limit_eb;
	  double limit_ebUP68;
	  double limit_ebDN68;
	  double width;
	  fileeb >> mass >> limit_eb>>limit_ebUP68>>limit_ebDN68>>width;
	  if (mass==600 || mass==650 || mass==850 || mass== 750 || mass== 800)continue;
	  std::cout<<mass<<" "<<limit_eb<<"  "<<width <<std::endl; 
    
	  if(width<1){
	    lim01eb.push_back(limit_eb);
	    lim01ebUP68.push_back(limit_ebUP68);
	    lim01ebDN68.push_back(limit_ebDN68);
	    mass01.push_back(mass);
	  }else if(width==2){
	    lim2eb.push_back(limit_eb);
	    lim2ebUP68.push_back(limit_ebUP68);
	    lim2ebDN68.push_back(limit_ebDN68);
	    mass2.push_back(mass);
	  }else if(width==5){
	    lim5eb.push_back(limit_eb);
	    lim5ebUP68.push_back(limit_ebUP68);
	    lim5ebDN68.push_back(limit_ebDN68);	   
	    mass5.push_back(mass);
	  }else if(width==7){
	    lim7eb.push_back(limit_eb);
	    lim7ebUP68.push_back(limit_ebUP68);
	    lim7ebDN68.push_back(limit_ebDN68);	   
	    mass7.push_back(mass);
	  }else if(width==10){
	    lim10eb.push_back(limit_eb);
	    lim10ebUP68.push_back(limit_ebUP68);
	    lim10ebDN68.push_back(limit_ebDN68);	   
	    mass10.push_back(mass);
	  }
	  
 }
  while (!file.eof()) {
    	  double mass;
	  double limit;
	  double limit_UP68;
	  double limit_DN68;
	  double width;
	  file >> mass >> limit>>limit_UP68>>limit_DN68>>width;
	  if (mass==600 || mass==650 || mass==850 || mass== 750 || mass== 800)continue;
	  std::cout<<mass<<" "<<limit<<"  "<<width <<std::endl;
	  if(width<1){
	    lim01.push_back(limit);
	    lim01UP68.push_back(limit_UP68);
	    lim01DN68.push_back(limit_DN68);	   

	  }
	  if(width==2){
	    lim2.push_back(limit);
	    lim2UP68.push_back(limit_UP68);
	    lim2DN68.push_back(limit_DN68);	   

	  }
	  if(width==5){
	    lim5.push_back(limit);
	    lim5UP68.push_back(limit_UP68);
	    lim5DN68.push_back(limit_DN68);	   

	  }
	  if(width==7){
	    lim7.push_back(limit);
	    lim7UP68.push_back(limit_UP68);
	    lim7DN68.push_back(limit_DN68);	   

	  }
	  if(width==10){
	    lim10.push_back(limit);
	    lim10UP68.push_back(limit_UP68);
	    lim10DN68.push_back(limit_DN68);	   
	   
	  }
	  
 }
  for(int i = 0; i< lim01.size();i++){
    double errDN=sqrt(pow(lim01eb[i], -2)*lim01DN68[i]*lim01DN68[i]+pow(lim01[i]/(lim01eb[i]*lim01eb[i]), -2)*lim01ebDN68[i]*lim01ebDN68[i]);
    double errUP=sqrt(pow(lim01eb[i], -2)*lim01UP68[i]*lim01UP68[i]+pow(lim01[i]/(lim01eb[i]*lim01eb[i]), -2)*lim01ebUP68[i]*lim01ebUP68[i]);
    errDN=errUP=0;
    gr01->SetPoint(i, mass01[i],(lim01eb[i]-lim01[i])/lim01eb[i] );
    gr01->SetPointError(i,0, 0, errDN, errUP);
  }
  for(int i = 0; i< lim2.size();i++){
    double errDN=sqrt(pow(lim2eb[i], -2)*lim2DN68[i]*lim2DN68[i]+pow(lim2[i]/(lim2eb[i]*lim2eb[i]), -2)*lim2ebDN68[i]*lim2ebDN68[i]);
    double errUP=sqrt(pow(lim2eb[i], -2)*lim2UP68[i]*lim2UP68[i]+pow(lim2[i]/(lim2eb[i]*lim2eb[i]), -2)*lim2ebUP68[i]*lim2ebUP68[i]);
    errDN=errUP=0;
    gr2->SetPoint(i, mass2[i],(lim2eb[i]-lim2[i])/lim2eb[i] );
    gr2->SetPointError(i,0, 0, errDN, errUP);
  }
  for(int i = 0; i< lim5.size();i++){
   double errDN=sqrt(pow(lim5eb[i], -2)*lim5DN68[i]*lim5DN68[i]+pow(lim5[i]/(lim5eb[i]*lim5eb[i]), -2)*lim5ebDN68[i]*lim5ebDN68[i]);
    double errUP=sqrt(pow(lim5eb[i], -2)*lim5UP68[i]*lim5UP68[i]+pow(lim5[i]/(lim5eb[i]*lim5eb[i]), -2)*lim5ebUP68[i]*lim5ebUP68[i]);
    errDN=errUP=0;
    gr5->SetPoint(i, mass5[i],(lim5eb[i]-lim5[i])/lim5eb[i] );
    gr5->SetPointError(i,0, 0, errDN, errUP);
  }
  for(int i = 0; i< lim7.size();i++){
   double errDN=sqrt(pow(lim7eb[i], -2)*lim7DN68[i]*lim7DN68[i]+pow(lim7[i]/(lim7eb[i]*lim7eb[i]), -2)*lim7ebDN68[i]*lim7ebDN68[i]);
    double errUP=sqrt(pow(lim7eb[i], -2)*lim7UP68[i]*lim7UP68[i]+pow(lim7[i]/(lim7eb[i]*lim7eb[i]), -2)*lim7ebUP68[i]*lim7ebUP68[i]);
    errDN=errUP=0;
    gr7->SetPoint(i, mass7[i],(lim7eb[i]-lim7[i])/lim7eb[i] );
    gr7->SetPointError(i,0, 0, errDN, errUP);
  }
  for(int i = 0; i< lim10.size();i++){
    double errDN=sqrt(pow(lim10eb[i], -2)*lim10DN68[i]*lim10DN68[i]+pow(lim10[i]/(lim10eb[i]*lim10eb[i]), -2)*lim10ebDN68[i]*lim10ebDN68[i]);
    double errUP=sqrt(pow(lim10eb[i], -2)*lim10UP68[i]*lim10UP68[i]+pow(lim10[i]/(lim10eb[i]*lim10eb[i]), -2)*lim10ebUP68[i]*lim10ebUP68[i]);
    errDN=errUP=0;
    gr10->SetPoint(i, mass10[i],(lim10eb[i]-lim10[i])/lim10eb[i]);
    gr10->SetPointError(i,0, 0, errDN, errUP);

  }

  TGraphSmooth sgr1;
  TGraph* gr01_sm = sgr1.Approx(gr01);
  TGraphSmooth sgr5;
  TGraph* gr5_sm = sgr5.Approx(gr5);
  TGraphSmooth sgr10;
  TGraph* gr10_sm = sgr10.Approx(gr10);

 TCanvas* c = new TCanvas("c", "c", 1);
 c->cd();
 TPaveText* label_cms = get_labelCMS(0, false);
 TPaveText* label_sqrt = get_labelSqrt(0);  
 TMultiGraph* mg = new TMultiGraph();
 mg->Add(gr01);
 //mg->Add(gr2);
 mg->Add(gr5);
 //mg->Add(gr7);
 mg->Add(gr10);
 gr10->SetMarkerColor(kAzure+3);
  gr7->SetMarkerColor(kAzure+2);
 gr5->SetMarkerColor(kAzure+1);
  gr2->SetMarkerColor(kViolet+1);
 gr01->SetMarkerColor(kViolet-9);
 gr10->SetLineColor(kAzure+3);
 gr7->SetLineColor(kAzure+2);
 gr5->SetLineColor(kAzure+1);
 gr2->SetLineColor(kViolet+1);
 gr01->SetLineColor(kViolet-9);

 gr10->SetLineWidth(2);
 //gr7->SetLineWidth(2);
 gr5->SetLineWidth(2);
 //gr2->SetLineWidth(2);
 gr01->SetLineWidth(2);

  
  mg->Draw("APL");
 
  mg->GetXaxis()->SetTitle("m_{H} [GeV]");
  mg->GetYaxis()->SetTitle("UL Improvement Including EE ");
  mg->GetYaxis()->SetRangeUser(-0.5, 0.5);
  //gr->GetZaxis()->SetTitle("#sigma #times BR (H #rightarrow #gamma #gamma)_{95%CL} (pb)");
 
  
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitleOffset(1.8);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleOffset(1.8);
 
  TLegend* leg = new TLegend(0.298851,0.2044755,0.54253,0.428252,"", "brNDC");
  // leg->SetHeader(TString::Format("Cat: All Events "));
  leg->SetBorderSize(0.);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  
  leg->AddEntry(gr10, "Width = 10% m_{H} [GeV]", "P");
  //leg->AddEntry(gr7, "Width = 7% m_{H} [GeV]", "P");
  leg->AddEntry(gr5, "Width = 5% m_{H} [GeV]", "P");
  //leg->AddEntry(gr2, "Width = 2% m_{H} [GeV]", "P");
  leg->AddEntry(gr01, "Width = 0.1 GeV", "P");
 
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  leg->Draw("same");


    c->SaveAs("plots/limitsResults_COMB_EB-EE.png");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/limitsResults_COMB_EB-EE.png");
    c->SaveAs("plots/limitsResults_COMB_EB-EE.pdf");
    c->SaveAs("/afs/cern.ch/user/s/soffi/www/Limits/limitsResults_COMB_EB-EE.pdf");
 
}






 

void makeRatioSyst(){


 ifstream file01;
 file01.open(TString::Format("limitsResults_COMB_w0.10.txt"));
 ifstream file01syst;
 file01syst.open(TString::Format("limitsResults_COMBWSyst_w0.10.txt"));
 TGraph* gr01= new TGraph();

 
 ifstream file5;
 file5.open(TString::Format("limitsResults_COMB_w5.00.txt"));
 ifstream file5syst;
 file5syst.open(TString::Format("limitsResults_COMBWSyst_w5.00.txt"));
 TGraph* gr5= new TGraph();
 


 ifstream file10;
 file10.open(TString::Format("limitsResults_COMB_w10.00.txt"));
 ifstream file10syst;
 file10syst.open(TString::Format("limitsResults_COMBWSyst_w10.00.txt"));
 TGraph* gr10= new TGraph();
 

 double limit01[7];
 double limit_syst01[7];
 double limit5[7];
 double limit_syst5[7];
 double limit10[7];
 double limit_syst10[7];
 double mass;
 double width;
 double limit_UP68;
 double limit_DN68;

 for(int i = 0;i<7; i++){
   std::cout<< i<< std::endl;
   file01 >> mass >> limit01[i]>>limit_UP68>>limit_DN68>>width;
   std::cout<<mass<<" "<<limit01[i]<<std::endl;   
 }
 for(int j = 0;j<7; j++) {
   std::cout<< j<< std::endl;
   file01syst >> mass >> limit_syst01[j] >>limit_UP68>>limit_DN68>> width;
   std::cout<<mass<<" "<<limit_syst01[j]<<std::endl;
   gr01->SetPoint(j, mass,(limit_syst01[j]-limit01[j])/limit_syst01[j] );   
 }
 for(int k = 0;k<7; k++){
std::cout<< k<< std::endl;
   file5 >> mass >> limit5[k]>>limit_UP68>>limit_DN68>>width;
   std::cout<<mass<<" "<<limit5[k]<<std::endl;   
 }
 for(int r = 0;r<7; r++) {
   file5syst >> mass >> limit_syst5[r] >> limit_UP68>>limit_DN68>>width;
   std::cout<<mass<<" "<<limit_syst5[r]<<std::endl;
   gr5->SetPoint(r, mass,(limit_syst5[r]-limit5[r])/limit_syst5[r] );   
 }
 for(int k = 0;k<7; k++){
   file10 >> mass >> limit10[k]>>limit_UP68>>limit_DN68>>width;
   std::cout<<mass<<" "<<limit10[k]<<std::endl;   
 }
 for(int r = 0;r<7; r++) {
   file10syst >> mass >> limit_syst10[r] >>limit_UP68>>limit_DN68>> width;
   std::cout<<mass<<" "<<limit_syst10[r]<<std::endl;
   gr10->SetPoint(r, mass,(limit_syst10[r]-limit10[r])/limit_syst10[r] );   
 }



 TCanvas* c = new TCanvas("c", "c", 1);
 c->cd();
 TPaveText* label_cms = get_labelCMS(0, false);
 TPaveText* label_sqrt = get_labelSqrt(0);  
 TMultiGraph* mg = new TMultiGraph();
 mg->Add(gr01);
 mg->Add(gr5);
 mg->Add(gr10);


 gr10->SetMarkerColor(kAzure+3);
 gr5->SetMarkerColor(kAzure+1);
 gr01->SetMarkerColor(kViolet-9);
 gr10->SetLineColor(kAzure+3);
 gr5->SetLineColor(kAzure+1);
 gr01->SetLineColor(kViolet-9);
 gr10->SetLineWidth(2);
 gr5->SetLineWidth(2);
 gr01->SetLineWidth(2);

  
  mg->Draw("APL");
 
  mg->GetXaxis()->SetTitle("m_{H} [GeV]");
  mg->GetYaxis()->SetTitle("UL Worsening Including Systematics");
  mg->GetYaxis()->SetRangeUser(-0.5, 0.5);
  
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitleOffset(1.8);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleOffset(1.8);
 
  TLegend* leg = new TLegend(0.238851,0.654755,0.52253,0.938252,"", "brNDC");
  leg->SetHeader(TString::Format("All Events "));
  leg->SetBorderSize(0.);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  
  leg->AddEntry(gr10, "Width = 10% m_{H} [GeV]", "P");
  leg->AddEntry(gr5, "Width = 5% m_{H} [GeV]", "P");
  leg->AddEntry(gr01, "Width = 0.1 GeV", "P");

  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  leg->Draw("same");


    c->SaveAs("plots/limitsResults_COMB_Syst.png");
    c->SaveAs("plots/limitsResults_COMB_Syst.pdf");

}





void makeplotSignalYield(){
  double lumi = 19620.;
  double mass[5] = {150., 200., 250., 300., 400.};
  double cat0[5] = {57.2743/lumi,73.98/lumi, 86.1266/lumi,99.099/lumi,115.707/lumi};
  double cat1[5] = {69.246/lumi,76.9358/lumi,75.4159/lumi,72.8432/lumi,69.0506/lumi };
  double cat2[5] = {27.697/lumi,33.0972/lumi,34.1901/lumi,37.6118/lumi,36.0121/lumi };
  double cat3[5] = {54.2063/lumi,60.2946/lumi,59.8451/lumi,60.5352/lumi,53.4554/lumi };
  double all[5] = {187.921/lumi, 217.376/lumi, 226.456/lumi,237.32/lumi,241.692/lumi};

  TGraph* y_all = new TGraph(5, mass, all);
  TGraph* y_cat0 = new TGraph(5, mass, cat0);
  TGraph* y_cat1 = new TGraph(5, mass, cat1);
  TGraph* y_cat2 = new TGraph(5, mass, cat2);
  TGraph* y_cat3 = new TGraph(5, mass, cat3);

  
  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();
  TPaveText* label_cms = get_labelCMS(0, false);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  TMultiGraph* mg = new TMultiGraph();
  mg->Add(y_cat0);
  mg->Add(y_cat1);
  mg->Add(y_cat2);
  mg->Add(y_cat3);
  mg->Add(y_all);
  
  y_all->SetMarkerColor(kRed+4);
  y_cat0->SetMarkerColor(kRed+2);
  y_cat1->SetMarkerColor(kOrange+10);
  y_cat2->SetMarkerColor(kOrange+1);
  y_cat3->SetMarkerColor(kOrange-2);
 
  /* mg->Draw("APL");
 
  mg->GetXaxis()->SetTitle("m_{H} [GeV]");
  mg->GetYaxis()->SetTitle("Signal Yield in 1 Pb");
  mg->GetYaxis()->SetRangeUser(0.0001, 0.03);
  //gr->GetZaxis()->SetTitle("#sigma #times BR (H #rightarrow #gamma #gamma)_{95%CL} (pb)");
 
  
  mg->GetXaxis()->SetTitleSize(0.04);
  mg->GetXaxis()->SetTitleOffset(1.8);
  mg->GetYaxis()->SetTitleSize(0.04);
  mg->GetYaxis()->SetTitleOffset(1.8);*/
 
  TLegend* leg = new TLegend(0.498851,0.6044755,0.84253,0.928252,"", "brNDC");
  // leg->SetHeader(TString::Format("Cat: All Events "));
  leg->SetBorderSize(0.);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(42);
  
  leg->AddEntry(y_all, "All  Events", "P");
  leg->AddEntry(y_cat0, "Cat 0", "P");
  leg->AddEntry(y_cat1, "Cat 1", "P");
  leg->AddEntry(y_cat2, "Cat 2", "P");
  leg->AddEntry(y_cat3, "Cat 3", "P");
 
  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  leg->Draw("same");

  /* c->SaveAs("plots/signalYield_cats.png");
     c->SaveAs("plots/signalYield_cats.pdf");*/
  double yield[5];
  
  double massBin[7]={125., 175., 225., 275., 325.,375., 425.};
  RooRealVar* PhotonsMass = new RooRealVar("PhotonsMass", "PhotonsMass",300,130, 950);
  //try to create the roohist func
  TH1F* h_0 = new TH1F("h_0", "h_0", 14, 130, 825);
  TH1F* h_1 = new TH1F("h_1", "h_1", 14, 130, 825);
  TH1F* h_2 = new TH1F("h_2", "h_2", 14, 130, 825);
  TH1F* h_3 = new TH1F("h_3", "h_3", 14, 130, 825);

  TH1F* h_all = new TH1F("h_all", "h_all", 15, 130, 875);

  for(int i=0;i<5;i++){
    h_0->SetBinContent(h_0->FindBin(mass[i]),cat0[i]);
    h_1->SetBinContent(h_1->FindBin(mass[i]),cat1[i]);
    h_2->SetBinContent(h_2->FindBin(mass[i]),cat2[i]);
    h_3->SetBinContent(h_3->FindBin(mass[i]),cat3[i]);
    h_all->SetBinContent(h_all->FindBin(mass[i]),all[i]);

    std::cout<<mass[i]<<std::endl;
    
  }
  h_0->SetBinContent(h_0->FindBin(350),(cat0[3]+cat0[4])/2);
  h_1->SetBinContent(h_1->FindBin(350),(cat1[3]+cat1[4])/2); 
  h_2->SetBinContent(h_2->FindBin(350),(cat2[3]+cat2[4])/2);
  h_3->SetBinContent(h_3->FindBin(350),(cat3[3]+cat3[4])/2);
  h_all->SetBinContent(h_all->FindBin(350),(all[3]+all[4])/2);
  double highmass[9]={450., 500., 550., 600., 650., 700., 750., 800., 850};
  for(int i = 0; i<9;i++){
   h_0->SetBinContent(h_0->FindBin(highmass[i]),(cat0[4]));
   h_1->SetBinContent(h_1->FindBin(highmass[i]),(cat1[4]));
   h_2->SetBinContent(h_2->FindBin(highmass[i]),(cat2[4]));
   h_3->SetBinContent(h_3->FindBin(highmass[i]),(cat3[4]));
   h_all->SetBinContent(h_all->FindBin(highmass[i]),(all[4]));
  }
  // h_all->GetYaxis()->SetRangeUser(0., 500.);

  
  // h_all->Draw();
  RooDataHist* rooData_0 = new RooDataHist("rooData_0", "roData_0",*PhotonsMass,h_0);
  RooHistFunc* rooFunc_0 = new RooHistFunc("rooFunc_0", "rooFunc_0", *PhotonsMass,*rooData_0,3);
  RooDataHist* rooData_1 = new RooDataHist("rooData_1", "roData_1",*PhotonsMass,h_1);
  RooHistFunc* rooFunc_1 = new RooHistFunc("rooFunc_1", "rooFunc_1", *PhotonsMass,*rooData_1,3);
  RooDataHist* rooData_2 = new RooDataHist("rooData_2", "roData_2",*PhotonsMass,h_2);
  RooHistFunc* rooFunc_2 = new RooHistFunc("rooFunc_2", "rooFunc_2", *PhotonsMass,*rooData_2,3);
  RooDataHist* rooData_3 = new RooDataHist("rooData_3", "roData_3",*PhotonsMass,h_3);
  RooHistFunc* rooFunc_3 = new RooHistFunc("rooFunc_3", "rooFunc_3", *PhotonsMass,*rooData_3,3);
  RooDataHist* rooData_all = new RooDataHist("rooData_all", "roData_all",*PhotonsMass,h_all);
  RooHistFunc* rooFunc_all = new RooHistFunc("rooFunc_all", "rooFunc_all", *PhotonsMass,*rooData_all,3);

  RooPlot* plot = PhotonsMass->frame();
  plot->GetYaxis()->SetRangeUser(0.0001, 0.03);
  plot->GetYaxis()->SetTitle("Interpolated Signal Yield in 1 Pb");
  plot->GetXaxis()->SetTitle("m_{H} [GeV]");
  plot->GetXaxis()->SetTitleSize(0.04);
  plot->GetXaxis()->SetTitleOffset(1.8);
  plot->GetYaxis()->SetTitleSize(0.04);
  plot->GetYaxis()->SetTitleOffset(1.8);
 
  // rooData_all->plotOn(plot);
  rooFunc_all->plotOn(plot,RooFit::LineColor(kRed+4));
  rooFunc_0->plotOn(plot, RooFit::LineColor(kRed+2));
  rooFunc_1->plotOn(plot, RooFit::LineColor(kOrange+10));
  rooFunc_2->plotOn(plot, RooFit::LineColor(kOrange+1));
  rooFunc_3->plotOn(plot, RooFit::LineColor(kOrange+2));

  plot->Draw();

  label_cms->Draw("same");
  label_sqrt->Draw("same");
 
  leg->Draw("same");

 
  c->SaveAs("plots/signalYield_roofcn.png");
  c->SaveAs("plots/signalYield_roofcn.pdf");

}

