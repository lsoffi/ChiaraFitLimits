#include"TPaveText.h"
#include "TChain.h"
#include "TH1F.h"
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include <stdio.h>
#include "TCanvas.h"
#include "TPad.h"
#include "TH1F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"
#include "TLatex.h"

using namespace std;
  

TString indirname = "./REDUCED/";
TString savedirname = "./PAS_PLOT";

double lumi = 4860.;

double fracQCD = 0.672;
double fracG = 0.328;

void protectZeroBin(TH1* hin) {
  hin->Print("v");

   for(int i = 1; i<=hin->GetNbinsX()+1; i++){
     
       if(hin->GetBinContent(i) <=0)  { 
	 double area = hin->GetBinWidth(i);
	 double bin_content = 3./(4*area);
	 double bin_error = sqrt(bin_content);
             hin->SetBinContent(i, bin_content);
             hin->SetBinError(i,  bin_error);
       }
     }
   }




std::string get_sqrtText() {

   char label_sqrt_text[150];
  
    sprintf( label_sqrt_text, "#sqrt{s} = 7 TeV");
    std::string returnString(label_sqrt_text);

  return returnString;

}


TPaveText* get_labelCMS( int legendQuadrant = 0 , std::string year="2011", bool sim=false, std::string run = "ALL") {

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
  if( legendQuadrant==0 ) cmslabel->SetTextAlign(11);
  cmslabel->SetTextSize(0.038);
  cmslabel->SetTextFont(42);
 
    std::string leftText;
   
    if(year == "2010")  leftText = "CMS Preliminary 2010, 34 pb^{-1}";
    if (sim)  leftText = "CMS Simulation"; //cwr ->remove 2011
    else {
      if(year == "2011" && run == "RUN2011A")  leftText = "CMS Preliminary RUN2011A 2.034 fb^{-1}";
      if(year == "2011" && run == "RUN2011B")  leftText = "CMS Preliminary 2011, 2.516 fb^{-1}";
      if(year == "2011" && run == "ALL")  leftText = "CMS 4.9 fb^{-1}"; //cwr ->remove 2011
      if(year == "none" && run == "ALL")  leftText = "CMS Data"; //cwr ->remove 2011
      if(year == "May2011")leftText = "CMS Preliminary 2011, 858.4 pb^{-1}";

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
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->AddText("#sqrt{s} = 7 TeV");
  return label_sqrt;

}


void hdensity( TH1F* hout){

   
  for(int i = 1; i<=hout->GetNbinsX(); i++){
   
      double area = hout->GetBinWidth(i);
      
      hout->SetBinContent(i,(double)hout->GetBinContent(i)/area);
      hout->SetBinError(i,(double)hout->GetBinError(i)/area);

      // std::cout << "i: " << i << " j: " << j <<
	//	" val: " << hout->GetBinContent(i,j) << " area: " << area <<std::endl; 

  }

}



void makePlotSignal_SameCat(){
  TFile* file;
  
  file = new TFile("/afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/AnalysisScripts/scriptoutput_HighMass/CMS-HGG.root", "READ");

  std::vector<std::string> cat;

  cat.push_back("0");
  cat.push_back("1");
  cat.push_back("2");
  cat.push_back("3");


  TH1F* h_ggh_150[4];
  TH1F* h_ggh_200[4];
  TH1F* h_ggh_250[4];
  TH1F* h_ggh_300[4];
  TH1F* h_ggh_400[4];
		 
  TH1F* h_vbf_150[4];
  TH1F* h_vbf_200[4];
  TH1F* h_vbf_250[4];
  TH1F* h_vbf_300[4];
  TH1F* h_vbf_400[4];
  
  TH1F* h_wzh_150[4];
  TH1F* h_wzh_200[4];
  TH1F* h_wzh_250[4];
  TH1F* h_wzh_300[4];
  TH1F* h_wzh_400[4];
	     	 
  TH1F* h_tth_150[4];
  TH1F* h_tth_200[4];
  TH1F* h_tth_250[4];
  TH1F* h_tth_300[4];
  TH1F* h_tth_400[4];
  

  TH1F* h_sum_150[4];
  TH1F* h_sum_200[4];
  TH1F* h_sum_250[4];
  TH1F* h_sum_300[4];
  TH1F* h_sum_400[4];

  TCanvas* c = new TCanvas("c", "c", 1); 

  TLegend* legSameCat = new TLegend(0.5, 0.65, 0.82, 0.9);
  legSameCat->SetFillColor(kWhite);
  legSameCat->SetBorderSize(0);

  TPaveText* label_cms = get_labelCMS(0, "2013", true, "");
  TPaveText* label_sqrt = get_labelSqrt(0);


  for(int k = 0;k<cat.size();k++){
                                     
  h_ggh_150[k] = (TH1F*) file->Get(("th1f_sig_ggh_mass_m150_cat"+cat[k]).c_str());
  h_ggh_200[k] = (TH1F*) file->Get(("th1f_sig_ggh_mass_m200_cat"+cat[k]).c_str());
  h_ggh_250[k] = (TH1F*) file->Get(("th1f_sig_ggh_mass_m250_cat"+cat[k]).c_str());
  h_ggh_300[k] = (TH1F*) file->Get(("th1f_sig_ggh_mass_m300_cat"+cat[k]).c_str());
  h_ggh_400[k] = (TH1F*) file->Get(("th1f_sig_ggh_mass_m400_cat"+cat[k]).c_str());
 		 
  h_vbf_150[k] = (TH1F*) file->Get(("th1f_sig_vbf_mass_m150_cat"+cat[k]).c_str());
  h_vbf_200[k] = (TH1F*) file->Get(("th1f_sig_vbf_mass_m200_cat"+cat[k]).c_str());
  h_vbf_250[k] = (TH1F*) file->Get(("th1f_sig_vbf_mass_m250_cat"+cat[k]).c_str());
  h_vbf_300[k] = (TH1F*) file->Get(("th1f_sig_vbf_mass_m300_cat"+cat[k]).c_str());
  h_vbf_400[k] = (TH1F*) file->Get(("th1f_sig_vbf_mass_m400_cat"+cat[k]).c_str());
 	     	    				    
  h_wzh_150[k] = (TH1F*) file->Get(("th1f_sig_wzh_mass_m150_cat"+cat[k]).c_str());
  h_wzh_200[k] = (TH1F*) file->Get(("th1f_sig_wzh_mass_m200_cat"+cat[k]).c_str());
  h_wzh_250[k] = (TH1F*) file->Get(("th1f_sig_wzh_mass_m250_cat"+cat[k]).c_str());
  h_wzh_300[k] = (TH1F*) file->Get(("th1f_sig_wzh_mass_m300_cat"+cat[k]).c_str());
  h_wzh_400[k] = (TH1F*) file->Get(("th1f_sig_wzh_mass_m400_cat"+cat[k]).c_str());
 	     	    				    
  h_tth_150[k] = (TH1F*) file->Get(("th1f_sig_tth_mass_m150_cat"+cat[k]).c_str());
  h_tth_200[k] = (TH1F*) file->Get(("th1f_sig_tth_mass_m200_cat"+cat[k]).c_str());
  h_tth_250[k] = (TH1F*) file->Get(("th1f_sig_tth_mass_m250_cat"+cat[k]).c_str());
  h_tth_300[k] = (TH1F*) file->Get(("th1f_sig_tth_mass_m300_cat"+cat[k]).c_str());
  h_tth_400[k] = (TH1F*) file->Get(("th1f_sig_tth_mass_m400_cat"+cat[k]).c_str());
 
  
  h_ggh_150[k]->SetLineColor(kAzure-2);
  h_ggh_200[k]->SetLineColor(kAzure-2);
  h_ggh_250[k]->SetLineColor(kAzure-2);
  h_ggh_300[k]->SetLineColor(kAzure-2);
  h_ggh_400[k]->SetLineColor(kAzure-2);
 
  h_vbf_150[k]->SetLineColor(kPink-2);
  h_vbf_200[k]->SetLineColor(kPink-2);
  h_vbf_250[k]->SetLineColor(kPink-2);
  h_vbf_300[k]->SetLineColor(kPink-2);
  h_vbf_400[k]->SetLineColor(kPink-2);

  h_wzh_150[k]->SetLineColor(kGreen-3);
  h_wzh_200[k]->SetLineColor(kGreen-3);
  h_wzh_250[k]->SetLineColor(kGreen-3);
  h_wzh_300[k]->SetLineColor(kGreen-3);
  h_wzh_400[k]->SetLineColor(kGreen-3);

  h_tth_150[k]->SetLineColor(kOrange+7);
  h_tth_200[k]->SetLineColor(kOrange+7);
  h_tth_250[k]->SetLineColor(kOrange+7);
  h_tth_300[k]->SetLineColor(kOrange+7);
  h_tth_400[k]->SetLineColor(kOrange+7);

  h_ggh_150[k]->SetLineWidth(2);
  h_ggh_200[k]->SetLineWidth(2);
  h_ggh_250[k]->SetLineWidth(2);
  h_ggh_300[k]->SetLineWidth(2);
  h_ggh_400[k]->SetLineWidth(2);
 
  h_vbf_150[k]->SetLineWidth(2);
  h_vbf_200[k]->SetLineWidth(2);
  h_vbf_250[k]->SetLineWidth(2);
  h_vbf_300[k]->SetLineWidth(2);
  h_vbf_400[k]->SetLineWidth(2);

  h_wzh_150[k]->SetLineWidth(2);
  h_wzh_200[k]->SetLineWidth(2);
  h_wzh_250[k]->SetLineWidth(2);
  h_wzh_300[k]->SetLineWidth(2);
  h_wzh_400[k]->SetLineWidth(2);

  h_tth_150[k]->SetLineWidth(2);
  h_tth_200[k]->SetLineWidth(2);
  h_tth_250[k]->SetLineWidth(2);
  h_tth_300[k]->SetLineWidth(2);
  h_tth_400[k]->SetLineWidth(2);


  

  h_sum_150[k] = (TH1F*)  h_ggh_150[k]->Clone();
  h_sum_150[k]->Add(h_vbf_150[k]);
  h_sum_150[k]->Add(h_wzh_150[k]);
  h_sum_150[k]->Add(h_tth_150[k]);

  h_sum_200[k] = (TH1F*)  h_ggh_200[k]->Clone();
  h_sum_200[k]->Add(h_vbf_200[k]);
  h_sum_200[k]->Add(h_wzh_200[k]);
  h_sum_200[k]->Add(h_tth_200[k]);

  h_sum_250[k] = (TH1F*)  h_ggh_250[k]->Clone();
  h_sum_250[k]->Add(h_vbf_250[k]);
  h_sum_250[k]->Add(h_wzh_250[k]);
  h_sum_250[k]->Add(h_tth_250[k]);
  
  h_sum_300[k] = (TH1F*)  h_ggh_300[k]->Clone();
  h_sum_300[k]->Add(h_vbf_300[k]);
  h_sum_300[k]->Add(h_wzh_300[k]);
  h_sum_300[k]->Add(h_tth_300[k]);

  h_sum_400[k] = (TH1F*)  h_ggh_400[k]->Clone();
  h_sum_400[k]->Add(h_vbf_400[k]);
  h_sum_400[k]->Add(h_wzh_400[k]);
  h_sum_400[k]->Add(h_tth_400[k]);



  h_sum_150[k]->SetLineWidth(2);
  h_sum_200[k]->SetLineWidth(2);
  h_sum_250[k]->SetLineWidth(2);
  h_sum_300[k]->SetLineWidth(2);
  h_sum_400[k]->SetLineWidth(2);

  h_sum_150[k]->SetLineColor(kBlack);
  h_sum_200[k]->SetLineColor(kBlack);
  h_sum_250[k]->SetLineColor(kBlack);
  h_sum_300[k]->SetLineColor(kBlack);
  h_sum_400[k]->SetLineColor(kBlack);



 
  legSameCat->Clear();

  legSameCat->AddEntry(h_ggh_150[k],"M = 150 GeV/c^{2}", "");
  legSameCat->AddEntry(h_ggh_200[k],"M = 200 GeV/c^{2}", "");
  legSameCat->AddEntry(h_ggh_250[k],"M = 250 GeV/c^{2}", "");
  legSameCat->AddEntry(h_ggh_300[k],"M = 300 GeV/c^{2}", "");
  legSameCat->AddEntry(h_ggh_400[k],"M = 400 GeV/c^{2}", "");



  //Cat k - alla masses only ggh
  std::cout<<"Cat "<<k<<" - all masses only ggh"<<std::endl;
  c->cd();
  c->SetLogy();
  legSameCat->SetHeader(("ggH prodcution mode - CAT "+cat[k]).c_str());

  h_ggh_150[k]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  h_ggh_150[k]->Draw("hist");
  h_ggh_200[k]->Draw("histsame");
  h_ggh_250[k]->Draw("histsame");
  h_ggh_300[k]->Draw("histsame");
  h_ggh_400[k]->Draw("histsame");
  legSameCat->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c->SaveAs(("img/Signal_ggh_allMasses_cat_"+cat[k]+".png").c_str());
  c->SaveAs(("img/Signal_ggh_allMasses_cat_"+cat[k]+".pdf").c_str());
   
  
  //Cat k - alla masses only vbf
  std::cout<<"Cat "<<k<<" - all masses only vbf"<<std::endl;
  legSameCat->SetHeader(("vbf prodcution mode - CAT "+cat[k]).c_str());

  h_vbf_150[k]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  h_vbf_150[k]->Draw("hist");
  h_vbf_200[k]->Draw("histsame");
  h_vbf_250[k]->Draw("histsame");
  h_vbf_300[k]->Draw("histsame");
  h_vbf_400[k]->Draw("histsame");
  legSameCat->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c->SaveAs(("img/Signal_vbf_allMasses_cat_"+cat[k]+".png").c_str());
  c->SaveAs(("img/Signal_vbf_allMasses_cat_"+cat[k]+".pdf").c_str());

  //Cat k - alla masses only wzh
  std::cout<<"Cat "<<k<<" - all masses only wzh"<<std::endl;
  legSameCat->SetHeader(("wzH prodcution mode - CAT "+cat[k]).c_str());
  h_wzh_150[k]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  h_wzh_150[k]->Draw("hist");
  h_wzh_200[k]->Draw("histsame");
  h_wzh_250[k]->Draw("histsame");
  h_wzh_300[k]->Draw("histsame");
  h_wzh_400[k]->Draw("histsame");
  legSameCat->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c->SaveAs(("img/Signal_wzh_allMasses_cat_"+cat[k]+".png").c_str());
  c->SaveAs(("img/Signal_wzh_allMasses_cat_"+cat[k]+".pdf").c_str());

  //Cat k - alla masses only tth
  std::cout<<"Cat "<<k<<" - all masses only tth"<<std::endl;
  legSameCat->SetHeader(("ttH prodcution mode - CAT "+cat[k]).c_str());
  h_tth_150[k]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  h_tth_150[k]->Draw("hist");
  h_tth_200[k]->Draw("histsame");
  h_tth_250[k]->Draw("histsame");
  h_tth_300[k]->Draw("histsame");
  h_tth_400[k]->Draw("histsame");
  legSameCat->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c->SaveAs(("img/Signal_tth_allMasses_cat_"+cat[k]+".png").c_str());
  c->SaveAs(("img/Signal_tth_allMasses_cat_"+cat[k]+".pdf").c_str());


  std::cout<<"Cat "<<k<<" - all masses all production modes"<<std::endl;
  legSameCat->SetHeader(("All prodcution mode - CAT "+cat[k]).c_str());
  h_sum_150[k]->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
  h_sum_150[k]->Draw("hist");
  h_sum_200[k]->Draw("histsame");
  h_sum_250[k]->Draw("histsame");
  h_sum_300[k]->Draw("histsame");
  h_sum_400[k]->Draw("histsame");
  legSameCat->Draw("same");
  label_cms->Draw("same");
  label_sqrt->Draw("same");

  c->SaveAs(("img/Signal_sum_allMasses_cat_"+cat[k]+".png").c_str());
  c->SaveAs(("img/Signal_sum_allMasses_cat_"+cat[k]+".pdf").c_str());



  }


 


 
    
}


