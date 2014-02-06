#include <iostream> 

#include <TROOT.h>

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>

#include <TH1F.h>
#include <TString.h>
#include <TCut.h>
#include <TF1.h>

//#define GAUSFIT
#define CI
#define REMOVEBADFIT
#define REMOVELOWSIGMA
#define HISTQUANTILES

void AddPoint(TGraphAsymmErrors *graph, int x_var, TH1F *hist, TString genFun, TString fitFun, TString option){
  //  TCanvas c_tmp("c_tmp","c_tmp"); // non attivare questa canvas altrimenti crasha!!!!
  //  c_tmp.cd();

  
  int i_point=graph->GetN();
  //  std::cout << "i_point = " << i_point << std::endl;

  // save the histogram
  hist->Draw();
  //  c.SaveAs("hist-"+genFun+"-"+fitFun+".C");
  int max_bin = hist->GetMaximumBin();

  // get variables 
  float  hist_min=std::max(hist->GetBinLowEdge(1), hist->GetMean()-4*hist->GetRMS());
  float  hist_max=std::min(hist->GetBinLowEdge(hist->GetNbinsX()-1), hist->GetMean()+4*hist->GetRMS());

  if(option.Contains("mean"))  hist->GetXaxis()->SetRangeUser(hist_min,hist_max);
  //  std::cout << "hist range: " << hist_min << "\t" << hist_max << std::endl;
  double y_var = hist->GetMean(); //hist->GetBinCenter(max_bin); //GetMean(); 
  double y_err_min=hist->GetMeanError();
  double y_err_max=hist->GetMeanError();


  if(option.Contains("median"))  
    {
      int nqm=1;
      double xqm[1];
      double yqm[1];
      xqm[0]=0.5;
      hist->GetQuantiles(nqm,yqm,xqm);
      y_var=yqm[0];
      std::cout << "**** median **** " << y_var << std::endl;
    }
  
  if(option.Contains("gaus")){
    //#ifdef GAUSFIT
    hist->Fit("gaus","Q");
    y_var=hist->GetFunction("gaus")->GetParameter(1);
    y_err_min=hist->GetFunction("gaus")->GetParameter(2)/sqrt(hist->Integral());
    y_err_max=y_err_min;
  //#endif
  }

  //#ifdef CI	
  if(option.Contains("68") || option.Contains("95"))
    {
      float band_lim=0;
      
      if (option.Contains("68")) band_lim=0.68;
      else if (option.Contains("95")) band_lim=0.95;
      
#ifndef HISTQUANTILES
      y_var=hist->GetBinCenter(max_bin);
      
      double integral = hist->Integral();
      
      int i_min=max_bin; 
      int i_max=max_bin;
      
      int n_move=0;
      //  std::cout << "------------------------------" << std::endl;
      //  std::cout << "cc = " << max_bin << "\t" << i_min << "\t" << i_max << "\t" << y_var << "\t" << hist->Integral(i_min,i_max)/integral << std::endl;
      while(hist->Integral(i_min,i_max)<band_lim*integral){
	//          std::cout << max_bin << "\t" << i_min << "\t" << i_max << "\t" << hist->Integral(i_min,i_max)/integral << std::endl;
	if(n_move%2==1)  i_min=std::max(--i_min,1); //i_min--; 
	else  i_max=std::min(++i_max,hist->GetNbinsX()-1);
	//    If(i_min < 1) i_min=1;
	//    if(i_max >= hist->GetNbinsX()-1) i_max=hist->GetNbinsX()-1;
	
	
	n_move++;
      }
      if(hist->Integral(i_min,i_max)/integral > band_lim*1.05){
	std::cerr << "[WARNING] " << std::endl;
	
	std::cout << i_min << "\t" << i_max << "\t" << hist->Integral(i_min,i_max)/integral << std::endl;
      }
      y_err_min=fabs(hist->GetBinCenter(i_min)-y_var);
      y_err_max=fabs(hist->GetBinCenter(i_max)-y_var);
      //  std::cout << y_var << "\t" << y_err_min << "\t" << y_err_max << std::endl;
      //#endif
#else
      int nq=3;
      double xq[3];
      double yq[3];
      if (option.Contains("68"))
	{
	  xq[0]=0.5-0.68/2.;
	  xq[1]=0.5;
	  xq[2]=0.5+0.68/2;
	}
      else if (option.Contains("95"))
	{
	  xq[0]=0.5-0.95/2.;
	  xq[1]=0.5;
	  xq[2]=0.5+0.95/2;
	}
      
      hist->GetQuantiles(nq,yq,xq);
      y_var=yq[1];
      y_err_min=TMath::Abs(yq[1]-yq[0]);
      y_err_max=TMath::Abs(yq[2]-yq[1]);
#endif
    }
  
  graph->SetPoint(i_point, x_var, y_var);
  
  
  graph->SetPointEYhigh(i_point, y_err_min);
  graph->SetPointEYlow(i_point,  y_err_max);
  
  
  i_point++;  
  graph->Set(i_point);

  return;
}


TGraphAsymmErrors *MakeGraph(TTree *tree, std::vector<int> mass_vec, TString genFun, TString fitFun, TString varName, TString option="mean"){

  //  TMultiGraph *graph = new TMultiGraph();
  TGraphAsymmErrors *graph = new TGraphAsymmErrors();

//   TCanvas c_graph("c_graph","c_graph");
//   c_graph.cd();

  for(std::vector<int>::const_iterator mass_itr=mass_vec.begin();
      mass_itr!=mass_vec.end();
      mass_itr++){

    //    std::cout << *mass_itr << std::endl;
    TString mass_string;
    mass_string+=(*mass_itr);
    //	mass_string.Print();
    TCut mass_cut("mass == "+mass_string);
    //	mass_cut.Print();
    TCut cut=mass_cut;
    cut+="genFun==\""+genFun+"\"";
    cut+="fitFun==\""+fitFun+"\"";
#ifdef REMOVEBADFIT
    cut+="mu>-15 && mu<15";
#endif
#ifdef REMOBELOWSIGMA
    cut+="sigma_mu>1e-2";
#endif
    cut.Print();


    if(varName.CompareTo("mu")==0) tree->Draw(varName+">>hist(600,-15,15)",cut);
    if(varName.CompareTo("bkgSig1fwhm/bkgTrue1fwhm")==0) tree->Draw(varName+">>hist(300,-5,5)",cut);
    if(varName.CompareTo("bkgErrSig1fwhm/bkgTrue1fwhm")==0) tree->Draw(varName+">>hist(400,0.,2.)",cut);
    if(varName.CompareTo("(bkgSig1fwhm-bkgTrue1fwhm)/bkgErrSig1fwhm")==0) tree->Draw(varName+">>hist(400,-3,3.)",cut);
    if(varName.CompareTo("(bkgSig2fwhm-bkgTrue2fwhm)/bkgErrSig2fwhm")==0) tree->Draw(varName+">>hist(400,-3,3.)",cut);
    if(varName.CompareTo("sigma_mu")==0) tree->Draw(varName+">>hist(400,0,20)",cut);
    if(varName.CompareTo("(mu-muTruth)/sigma_mu")==0) tree->Draw(varName+">>hist(300,-3,3)",cut);
    if(varName.CompareTo("(mu-muTruth)/bkgErrNormSig2fwhm")==0) tree->Draw(varName+">>hist(400,-8,8)",cut);
    if(varName.CompareTo("bkgErrNormSig2fwhm")==0) tree->Draw(varName+">>hist(200,0.,10.)",cut);

    
    TH1F *hist = (TH1F *) gROOT->FindObject("hist");
   
    if(hist==NULL || hist->GetEntries()==0){
      std::cerr << "no entries for cut: " << std::endl;
      cut.Print();
    } else     AddPoint(graph, *mass_itr, hist, genFun, fitFun, option);
    //    graph->SetMinimum(-100);
    //    graph->SetMaximum(100);

    delete hist;
  }
  
  
//       std::cout << "Np = " << i_point << std::endl;
//       }
//       c.cd();

  return graph;
}

