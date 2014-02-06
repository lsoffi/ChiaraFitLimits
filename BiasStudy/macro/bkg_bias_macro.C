#define SIGMA
  //#define PULL
  //#define FABIANPULL
#include<map>

void bkg_bias_macro(
		    int mu_gen=0,
		    TString fileName="biasCheck.root"
		    )
{
  gROOT->ProcessLine(".L macro/bkg_bias_fun.C+");
  
  //  TString fileName="biasCheck-mu_100.root";
  
  //  TLatex latex_
  
  TString muGen;
//   int mu_min;
//   int mu_max;

  //  TString muGen="mu_1/";      
  TString muGen="mu_";      
  muGen+=mu_gen;
  muGen+="/";      

//   if(mu_gen==0){
//     //   fileName="biasCheck.root";
//     mu_min=-0.5;
//     mu_max=0.5;
//   }
//   if(mu_gen==1){
//     //    fileName="biasCheck-sigGen_1.root";	
//     //  muGen="mu_1/";      
//     mu_min=0;
//     mu_max=2;
//   }

//   if(mu_gen==3){
//     //    fileName="biasCheck-sigGen_1.root";	
//     //  muGen="mu_1/";      
//     mu_min=1;
//     mu_max=5;
//   }

//   if(mu_gen==10){
//     //    fileName="biasCheck-sigGen_10.root";	
//     //    muGen="mu_10/";      
//     mu_min=5;
//     mu_max=15;
//   }  

  
  TFile f(fileName);
  f.cd();


  TTree *tree = (TTree *) f.Get("muTree");

//   TString variable[5] = { "mu", "sigma_mu" , "(mu-muTruth)/sigma_mu" , "(mu-muTruth)/bkgErrNormSig2fwhm" , "bkgErrNormSig2fwhm" };
//   TString variableLabel[5] = { "#mu", "#sigma_{#mu}" , "#frac{#mu-#mu_{Truth}}{#sigma_{#mu}}" , "#frac{#mu-#mu_{Truth}}{bkgErr}" , "bkgErrNormSig2fwhm" };
//   TString variableName[5] = { "mu", "sigma_mu" , "pull" , "fabianPull" , "bkgErrNormSig2fwhm" };

  TString variable[5] = { "bkgSig1fwhm/bkgTrue1fwhm", "bkgErrSig1fwhm/bkgTrue1fwhm" , "(bkgSig1fwhm-bkgTrue1fwhm)/bkgErrSig1fwhm" , "(bkgSig2fwhm-bkgTrue2fwhm)/bkgErrSig2fwhm" , "bkgErrNormSig2fwhm" };
  TString variableLabel[5] = { "normBkg", "normBkgErr" , "#Delta_{BkgFit}/#sigma_{BkgFit}" , "#Delta_{BkgFit}/#sigma_{BkgFit} (2FWHM)" , "bkgErrNormSig2fwhm" };
  TString variableName[5] = { "normBkg", "normBkgErr" , "fabianPull" , "fabianPull2fwhm" , "bkgErrNormSig2fwhm" };

//   float limitLow[5] = { -1. ,3., -1., -1. , 0.};
//   float limitHigh[5] = { 1., 7., 1., 1. , 2.} ;
//   float limitBandLow[5] = { -15 , 0. , -4., -10. , 0. };
//   float limitBandHigh[5] = { 15 , 12. , 4. , 10. , 5. };

//   float limitLow[5] = { -1. ,0., -1., -1. , 0.};
//   float limitHigh[5] = { 1., 3., 1., 1. , 1.} ;
//   float limitBandLow[5] = { -15 , 0. , -4., -7. , 0. };
//   float limitBandHigh[5] = { 15 , 3. , 4. , 7. , 3. };


  float limitLow[5] = { -1. ,0., -0.8, -0.8. , 0.};
  float limitHigh[5] = { 1., 3., 0.8, 0.8 , 1.} ;
  float limitBandLow[5] = { -15 , 0. , -2.5, -2.5 , 0. };
  float limitBandHigh[5] = { 15 , 3. , 2.5 , 2.5 , 3. };

  std::vector<TString> genFun_vec;
  std::vector<TString> fitFun_vec;
  genFun_vec.push_back("1exp");
  genFun_vec.push_back("1pow");
  genFun_vec.push_back("2pol");

  
  fitFun_vec.push_back("1exp");
  //  fitFun_vec.push_back("2exp"); 
  fitFun_vec.push_back("1pow"); 
  //  fitFun_vec.push_back("2pow"); 
  //fitFun_vec.push_back("1pol"); 
  fitFun_vec.push_back("2pol"); 
  fitFun_vec.push_back("3pol"); 
  fitFun_vec.push_back("4pol");

//     fitFun_vec.push_back("pow2_2");  
//     fitFun_vec.push_back("exp2_2"); 
//     fitFun_vec.push_back("exp3_2"); 
//     fitFun_vec.push_back("exp3_3");
  //fitFun_vec.push_back("4pol");
  
//   fitFun_vec.push_back("pow2");
//   fitFun_vec.push_back("exp2");
//   fitFun_vec.push_back("exp3");


  std::vector<int> mass_vec;
  mass_vec.push_back(110);
  mass_vec.push_back(115);
  mass_vec.push_back(120);
  mass_vec.push_back(125);
  mass_vec.push_back(130);
  mass_vec.push_back(140);
  // prendo la lista delle funzioni di fit e di gen

  //  TDirectory dir("ciao","ciao");
  //  dir.cd();
  // 
  TCanvas c("canvass","c");
  //  c.Divide(2,1);
//   c.cd(0);
   //  pad1->Draw();
  //  pad2->Draw();
  //  c.cd(0);
  //  TH1D *prj_x=NULL;



  for (int ivar=2;ivar<3;ivar++)
    {
      gSystem->mkdir("img/"+muGen+variableName[ivar],true);
      
      for(unsigned int i_gen=0; i_gen < genFun_vec.size(); i_gen++){
	
	std::cout << "============================== " << genFun_vec[i_gen] << std::endl;
	
	for(unsigned int i_fit=0; i_fit < fitFun_vec.size(); i_fit++){
	  


      
	  std::cout << "------------------------------ " << fitFun_vec[i_fit] << std::endl;
	  TPad pad1("pad1","pad1",0,0.5,1,1);
	  TPad pad2("pad2","pad2",0,0,1,0.5);

	  pad1.SetLeftMargin(pad1.GetLeftMargin()-0.04);
	  pad2.SetLeftMargin(pad2.GetLeftMargin()-0.04);
	  pad1.SetRightMargin(pad1.GetRightMargin()+0.09);
	  pad2.SetRightMargin(pad2.GetRightMargin()+0.09);

	  c.cd();
	  pad1.Draw();

	  c.cd();
	  pad2.Draw();

	  pad1.cd();

	  TGraphAsymmErrors *mu_graph   = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], variable[ivar], "median");
	  TGraphAsymmErrors *mu_68_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], variable[ivar], "68");
	  TGraphAsymmErrors *mu_95_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], variable[ivar], "95");
      
	  mu_graph->GetXaxis()->SetTitle("m_{H}");
	  mu_graph->SetMarkerStyle(20);
	  mu_graph->SetMarkerSize(0.8);
	  mu_graph->SetLineWidth(2);
	  mu_graph->GetXaxis()->SetLabelSize(0.05);
	  mu_graph->GetYaxis()->SetLabelSize(0.05);
	  mu_graph->GetXaxis()->SetTitleSize(0.08);
	  mu_graph->GetYaxis()->SetTitleSize(0.08);
	  mu_graph->GetXaxis()->SetTitleOffset(0.8);
	  mu_graph->GetYaxis()->SetTitleOffset(0.5);
	  mu_graph->GetYaxis()->SetTitle(variableLabel[ivar]);
	  mu_graph->GetYaxis()->SetRangeUser(limitLow[ivar],limitHigh[ivar]);
      
	  mu_68_band->SetFillColor(kBlue-10);
	  mu_68_band->SetFillStyle(1);
	  mu_68_band->SetDrawOption("3");
	  mu_68_band->GetYaxis()->SetRangeUser(limitBandLow[ivar],limitBandHigh[ivar]);
	  
	  mu_95_band->SetFillColor(kOrange-1);
	  mu_95_band->SetFillStyle(1);
	  mu_95_band->SetDrawOption("3");
	  mu_95_band->GetYaxis()->SetRangeUser(limitBandLow[ivar],limitBandHigh[ivar]);
	  mu_95_band->GetXaxis()->SetTitle("m_{H}");
	  mu_95_band->GetXaxis()->SetLabelSize(0.05);
	  mu_95_band->GetYaxis()->SetLabelSize(0.05);
	  mu_95_band->GetXaxis()->SetTitleSize(0.08);
	  mu_95_band->GetYaxis()->SetTitleSize(0.08);
	  mu_95_band->GetXaxis()->SetTitleOffset(0.8);
	  mu_95_band->GetYaxis()->SetTitleOffset(0.5);
	  mu_95_band->GetYaxis()->SetTitle(variableLabel[ivar]);

	  TLegend* leg= new TLegend(0.79,0.5,0.97,0.8);
	  leg->SetBorderSize(1);
	  leg->SetFillStyle(1001);
	  leg->SetFillColor(0);	  
	  leg->SetTextSize(0.075);

	  char line[300];

	  leg->AddEntry(mu_95_band,"95%","f");
	  leg->AddEntry(mu_68_band,"68%","f");
	  leg->AddEntry(mu_graph, "mean","pl");
	  //      Pad1.cd();
	  //      c.cd(0);

// 	  TLatex* lat=new TLatex();

// 	  //      lat->SetNDC();
// 	  lat->SetTextSize(0.33);
// 	  lat->SetLineWidth(2);
// 	  lat->SetTextColor(1);

// 	  float yhi = 0.1;
// 	  float xmin= 0.1;
// 	  float ypass= 0.05;
	  
// 	  char tmpStr[300];
	  
	  //	  pad1.cd();
	  mu_95_band->Draw("a3");
	  mu_68_band->Draw("3 same");
	  mu_graph->Draw("pXl same");

//  	  sprintf(tmpStr,"Gen fun: ");
//  	  sprintf(line,"%s%s", tmpStr, genFun_vec[i_gen].Data());

//   	  std::cout << "++++++++++ " << line << std::endl;
//   	  lat->DrawLatex(xmin, yhi, line);
    
//  	  //       sprintf(tmpStr,"95%");
//  	  sprintf(line,"%s", "95\%");
//   	  std::cout << "++++++++++ " << line << std::endl;
//   	  lat->DrawLatex(xmin, yhi-ypass, line);
    
//  	  //       sprintf(tmpStr,"68%");
//  	  sprintf(line,"%s", "68\%");
//   	  std::cout << "++++++++++ " << line << std::endl;
//   	  lat->DrawLatex(xmin, yhi-2*ypass, line);
//  	  //      mu_graph->Draw("Pl same");

	  leg->Draw();      
	  TPaveText *pt = new TPaveText(0.82,0.8,0.99,0.95,"NDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(0);	  
	  pt->SetTextSize(0.07);
	  TString genFunName= genFun_vec[i_gen];
	  TString fitFunName= fitFun_vec[i_fit];

 	  if (genFunName.CompareTo("2pol")==0)
 	    genFunName="1pol";
 	  if (genFunName.CompareTo("3pol")==0)
 	    genFunName="2pol";
 	  if (genFunName.CompareTo("4pol")==0)
 	    genFunName="3pol";
	  

 	  if (fitFunName.CompareTo("2pol")==0)
 	    fitFunName="1pol";
 	  if (fitFunName.CompareTo("3pol")==0)
 	    fitFunName="2pol";
 	  if (fitFunName.CompareTo("4pol")==0)
 	    fitFunName="3pol";
	  

 	  sprintf(line,"Gen fun: %s", genFunName.Data() );
	  pt->AddText(line);
  	  sprintf(line,"Fit fun: %s",  fitFunName.Data() );
 	  pt->AddText(line);
  	  pt->Draw();
	  //      lat->Paint();

	  //      pad1->Draw();
	  
	  pad2.cd();
	  //	  leg->Clear();
	  //	  leg->AddEntry(mu_graph, "median value", "pl");
	  //	  leg->Draw();
	  mu_graph->Draw("aPl");
	  if (ivar>1 && ivar<4)
	    {
	      TLine l1(107,-0.2,143,-0.2);
	      l1.SetLineColor(2);
	      l1.Draw();
	      TLine l2(107,0.2,143,0.2);
	      l2.SetLineColor(2);
	      l2.Draw();
	    }
	  mu_graph->Draw("PlSAME");

	  //      pad1.Draw();

	  pad1.cd();
	  c.SaveAs("./img/"+muGen+variableName[ivar]+"/"+variableName[ivar]+"_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");
	  c.SaveAs("./img/"+muGen+variableName[ivar]+"/"+variableName[ivar]+"_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".C");
	  c.SaveAs("./img/"+muGen+variableName[ivar]+"/"+variableName[ivar]+"_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".eps");
	  
	}
      }
    }
  /*
#ifdef SIGMA
      TGraphAsymmErrors *sigma_mu_graph = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], "sigma_mu","mean");
      TGraphAsymmErrors *sigma_mu_68_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], "sigma_mu", "68");
      TGraphAsymmErrors *sigma_mu_95_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], "sigma_mu", "95");


      std::cout << "Plotting" << std::endl;
      //sigma_mu_graph->Draw("AP");
      sigma_mu_graph->GetXaxis()->SetTitle("H mass");
      sigma_mu_graph->GetYaxis()->SetTitle("#sigma_{#mu}");
      sigma_mu_graph->GetYaxis()->SetRangeUser(3.,7.);
      sigma_mu_68_band->SetFillColor(kYellow);
      sigma_mu_68_band->SetFillStyle(1);
      sigma_mu_68_band->SetDrawOption("3");
      sigma_mu_68_band->GetYaxis()->SetRangeUser(0.,10.);
      
      sigma_mu_95_band->SetFillColor(4);
      sigma_mu_95_band->SetFillStyle(1);
      sigma_mu_95_band->SetDrawOption("3");
      sigma_mu_95_band->GetYaxis()->SetRangeUser(0,15.);
      sigma_mu_95_band->GetXaxis()->SetTitle("H mass");
      sigma_mu_95_band->GetYaxis()->SetTitle("#sigma_{#mu}");

      //      std::cout << "-------" << std::endl;
      //      pad1.Draw();
      pad1.cd();
      sigma_mu_95_band->Draw("a3");
      sigma_mu_68_band->Draw("3 same");
      sigma_mu_graph->Draw("Pl same");


      pad2.cd();
      sigma_mu_graph->Draw("aPl");

      //std::cout << "-------============" << std::endl;

      c.SaveAs("./img/"+muGen+"sigma_mu/sigma_mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");
      c.SaveAs("./img/"+muGen+"sigma_mu/sigma_mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".C");
      c.SaveAs("./img/"+muGen+"sigma_mu/sigma_mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".eps");
#endif

#ifdef PULL      
      TString pullString="(mu-muTruth)/sigma_mu";
      TGraphAsymmErrors *pull_graph = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], pullString, "median");
      TGraphAsymmErrors *pull_68_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], pullString, "68");
      TGraphAsymmErrors *pull_95_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], pullString, "95");
      
      
      pull_graph->Draw("AP");
      pull_graph->GetXaxis()->SetTitle("H mass");
      pull_graph->GetYaxis()->SetTitle("pull");
      pull_graph->GetYaxis()->SetRangeUser(-1,1);
      pull_68_band->SetFillColor(kYellow);
      pull_68_band->SetFillStyle(1);
      pull_68_band->SetDrawOption("3");
      pull_68_band->GetYaxis()->SetRangeUser(-2,2);
      
      pull_95_band->SetFillColor(4);
      pull_95_band->SetFillStyle(1);
      pull_95_band->SetDrawOption("3");
      pull_95_band->GetYaxis()->SetRangeUser(-4,4);
      pull_95_band->GetXaxis()->SetTitle("H mass");
      pull_95_band->GetYaxis()->SetTitle("pull ((#mu - #mu_{truth})/#sigma_{#mu})");

      pad1.cd();
      pull_95_band->Draw("a3");
      pull_68_band->Draw("3 same");
      pull_graph->Draw("Pl same");
      //      leg.Draw();
      
      pad2.cd();
      pull_graph->Draw("aPl");

      c.SaveAs("./img/"+muGen+"pull/pull_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");
      c.SaveAs("./img/"+muGen+"pull/pull_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".C");
      c.SaveAs("./img/"+muGen+"pull/pull_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".eps");

      //  leg.Clear(); // devo fare il clear prima di distruggere gli oggetti a cui e' riferito      

#endif

#ifdef FABIANPULL      
      TString pullString="(mu-muTruth)/bkgErrNormSig2fwhm";
      TGraphAsymmErrors *pull_graph = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], pullString, "median");
      TGraphAsymmErrors *pull_68_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], pullString, "68");
      TGraphAsymmErrors *pull_95_band = MakeGraph(tree, mass_vec, genFun_vec[i_gen], fitFun_vec[i_fit], pullString, "95");
      
      
      pull_graph->Draw("AP");
      pull_graph->GetXaxis()->SetTitle("H mass");
      pull_graph->GetYaxis()->SetTitle("pull");
      pull_graph->GetYaxis()->SetRangeUser(-1,1);
      pull_68_band->SetFillColor(kYellow);
      pull_68_band->SetFillStyle(1);
      pull_68_band->SetDrawOption("3");
      pull_68_band->GetYaxis()->SetRangeUser(-5,5);
      
      pull_95_band->SetFillColor(4);
      pull_95_band->SetFillStyle(1);
      pull_95_band->SetDrawOption("3");
      pull_95_band->GetYaxis()->SetRangeUser(-12,12);
      pull_95_band->GetXaxis()->SetTitle("H mass");
      pull_95_band->GetYaxis()->SetTitle("pull ((#mu - #mu_{truth})/#sigma_{#mu})");

      pad1.cd();
      pull_95_band->Draw("a3");
      pull_68_band->Draw("3 same");
      pull_graph->Draw("Pl same");
      //      leg.Draw();
      
      pad2.cd();
      pull_graph->Draw("aPl");

      c.SaveAs("./img/"+muGen+"fabianpull/pull_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");
      c.SaveAs("./img/"+muGen+"fabianpull/pull_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".C");
      c.SaveAs("./img/"+muGen+"fabianpull/pull_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".eps");

      //  leg.Clear(); // devo fare il clear prima di distruggere gli oggetti a cui e' riferito      

#endif




//       delete mu_graph;
//       delete mu_68_band;
//       delete mu_95_band;
//       delete sigma_mu_graph;
//       delete sigma_mu_68_band;
//       delete sigma_mu_95_band;
//       delete pull_graph;
//       delete pull_68_band;
//       delete pull_95_band;

//      delete pad1;
//      delete pad2;
  //  leg.Clear();
  
  
}


#ifdef shervin

std::map<TString, int> genFun_map;
  std::map<TString, int> fitFun_map;

  char fitFun_[30];
  char genFun_[30];
  tree->SetBranchAddress("fitFun", fitFun_);
  tree->SetBranchAddress("genFun", genFun_);

  Long64_t nEntries = tree->GetEntries();
  nEntries=20;
  for(Long64_t jEntry=0; jEntry < nEntries; jEntry++){
    tree->GetEntry(jEntry);
    //std::cout << fitFun_ << std::endl;
    fitFun_map[fitFun_]=1;
    genFun_map[genFun_]=1;

  }    
    
  for(std::map<TString, int>::const_iterator itr=fitFun_vec.begin();
      itr! fitFun_vec.end(); itr++){
    fitFun_vec.push_back(itr->first);
  }
				  
  for(std::map<TString, int>::const_iterator itr=genFun_vec.begin();
      itr! genFun_vec.end(); itr++){
    genFun_vec.push_back(itr->first);
  }


#endif

#ifdef oldPlot
      TGraphAsymmErrors mu_graph;
      TGraphAsymmErrors sigma_mu_graph;
      int i_point=0;

      for(std::vector<int>::const_iterator mass_itr=mass_vec.begin();
	  mass_itr!=mass_vec.end();
	  mass_itr++){

	std::cout << *mass_itr << std::endl;
       	TString mass_string;
	mass_string+=(*mass_itr);
	//	mass_string.Print();
	TCut mass_cut("mass == "+mass_string);
	//	mass_cut.Print();
	TCut cut=mass_cut;
	cut+="genFun==\""+genFun_vec[i_gen]+"\"";
	cut+="fitFun==\""+fitFun_vec[i_fit]+"\"";
	//	cut.Print();



	tree->Draw("mu>>mu_hist(100,-10,10)",cut);
	TH1F *mu_hist = (TH1F *) gROOT->FindObject("mu_hist");
	mu_hist->Draw();
	c.SaveAs("mu_hist-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");
	double mu = mu_hist->GetMean();
	
	int max_bin = mu_hist->GetMaximumBin();
	double mu_max = mu; //mu_hist->GetBinContent(max_bin);

	double integral = mu_hist->Integral();
	
	int i_min=max_bin;
	int i_max=max_bin;

	int n_move=0;
	std::cout << "------------------------------" << std::endl;
	std::cout << "cc = " << max_bin << "\t" << i_min << "\t" << i_max << "\t" << mu_max << "\t" << mu_hist->Integral(i_min,i_max)/integral << std::endl;

	while(mu_hist->Integral(i_min,i_max)<0.68*integral){
	  if(n_move%2)  i_min--; //=std::max(i_min--,1);
	  else  i_max++; //=std::min(i_max++,mu_hist->GetNbinsX()-1);
	  if(i_min < 1) i_min=1;
	  if(i_max >= mu_hist->GetNbinsX()-1) i_max=mu_hist->GetNbinsX()-1;

	  //	  std::cout << max_bin << "\t" << i_min << "\t" << i_max << "\t" << mu_hist->Integral(i_min,i_max)/integral << std::endl;
	  n_move++;
	}
	if(mu_hist->Integral(i_min,i_max)/integral > 0.70) std::cerr << "[WARNING] " << std::endl;

	std::cout << i_min << "\t" << i_max << "\t" << mu_hist->Integral(i_min,i_max)/integral << std::endl;

	mu_graph.SetPoint(i_point, *mass_itr, mu_max);

	mu_graph.SetPointEYhigh(i_point, mu_hist->GetBinContent(i_max));
	mu_graph.SetPointEYlow(i_point, mu_hist->GetBinContent(i_min));

	mu_graph.SetPointEYhigh(i_point, mu_hist->GetMeanError());
	mu_graph.SetPointEYlow(i_point, mu_hist->GetMeanError());
	
	i_point++;
	//   Long64_t nEntries = tree->GetEntries();
	//   nEntries=20;
	//   for(Long64_t jEntry=0; jEntry < nEntries; jEntry++){
	//     tree->GetEntry(jEntry);

#ifdef ciao  
	tree->Draw("mu:mass>>mu_vs_mass_hist","genFun==\""+genFun_vec[i_gen]+"\" && fitFun==\""+fitFun_vec[i_fit]+"\"");
	TH2F *mu_hist= (TH2F *) gROOT->FindObject("mu_vs_mass_hist");
	prj_x = mu_hist->ProfileX();
	prj_x->Draw();
	c.SaveAs("mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");

	tree->Draw("mu:mass>>mu_vs_mass_hist","genFun==\""+genFun_vec[i_gen]+"\" && fitFun==\""+fitFun_vec[i_fit]+"\"");
	TH2F *mu_hist= (TH2F *) gROOT->FindObject("mu_vs_mass_hist");
	prj_x = mu_hist->ProfileX();
	prj_x->Draw();
	c.SaveAs("mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");

	delete hist;
	delete prj_x;
#endif

	mu_graph.Set(i_point);
	std::cout << "Np = " << i_point << std::endl;
      }
      c.cd();
      mu_graph.Draw("AP");
      mu_graph.GetXaxis()->SetTitle("H mass");
      mu_graph.GetYaxis()->SetTitle("#mu");
      c.SaveAs("mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");


#ifdef ciao
      i_point=0;
      tree->Draw("sigma_mu>>sigma_mu_hist",cut);
      TH1F *sigma_mu_hist = (TH1F *) gROOT->FindObject("sigma_mu_hist");
      //	sigma_mu_hist->Draw();
      //	c.SaveAs("sigma_mu_hist-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");
      double sigma_mu = sigma_mu_hist->GetMean();
	
      max_bin = sigma_mu_hist->GetMaximumBin();
      double sigma_mu_max = sigma_mu; //mu_hist->GetBinContent(max_bin);
	
      integral = sigma_mu_hist->Integral();
      
      i_min=max_bin;
      i_max=max_bin;
      
      n_move=0;
      std::cout << "------------------------------" << std::endl;
      std::cout << "cc = " << max_bin << "\t" << i_min << "\t" << i_max << "\t" << sigma_mu_max << "\t" << sigma_mu_hist->Integral(i_min,i_max)/integral << std::endl;
      
      while(sigma_mu_hist->Integral(i_min,i_max)<0.68*integral){
	if(n_move%2)  i_min--; //=std::max(i_min--,1);
	else  i_max++; //=std::min(i_max++,sigma_mu_hist->GetNbinsX()-1);
	if(i_min < 1) i_min=1;
	if(i_max >= sigma_mu_hist->GetNbinsX()-1) i_max=sigma_mu_hist->GetNbinsX()-1;
	
	//	  std::cout << max_bin << "\t" << i_min << "\t" << i_max << "\t" << sigma_mu_hist->Integral(i_min,i_max)/integral << std::endl;
	n_move++;
      }
      if(sigma_mu_hist->Integral(i_min,i_max)/integral > 0.70) std::cerr << "[WARNING] " << std::endl;

      std::cout << i_min << "\t" << i_max << "\t" << sigma_mu_hist->Integral(i_min,i_max)/integral << std::endl;

      sigma_mu_graph.SetPoint(i_point, *mass_itr, sigma_mu_max);

      sigma_mu_graph.SetPointEYhigh(i_point, sigma_mu_hist->GetBinContent(i_max));
      sigma_mu_graph.SetPointEYlow(i_point, sigma_mu_hist->GetBinContent(i_min));

      sigma_mu_graph.SetPointEYhigh(i_point, sigma_mu_hist->GetRMS());
      sigma_mu_graph.SetPointEYlow(i_point, sigma_mu_hist->GetRMS());
	
      i_point++;

      sigma_mu_graph.Set(i_point);
      std::cout << "Np = " << i_point << std::endl;
    }
    c.cd();
    sigma_mu_graph.Draw("AP");
    c.SaveAs("sigma_mu_vs_mass-"+genFun_vec[i_gen]+"-"+fitFun_vec[i_fit]+".png");

#endif
#endif
      */
}
