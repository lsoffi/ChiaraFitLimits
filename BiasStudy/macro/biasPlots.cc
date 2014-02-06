#include <iostream> 
#include <iostream> 

#include <TROOT.h>

#include <TTree.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TMultiGraph.h>

#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCut.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TLegend.h>
#include <TFile.h>



#include <string>
#include <sstream>
#include <cmath>
#include <iostream>

std::string toPrecision(double num, int n) {
    
    if(num == 0) {
      return "0";
    }

    double d = std::ceil(std::log10(num < 0 ? -num : num));
    int power = n - (int)d;
    double magnitude = std::pow(10., power);
    long shifted = ::round(num*magnitude);

    std::ostringstream oss;
    oss << shifted/magnitude;
    return oss.str();
}

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



double* getPar(double* array, double* par){
  double temp = 0.;
  double sum;
  double sum2;
  double max;
  double mean;
  double rms;

  for(int i = 0;i<12; i++){
    if(fabs(array[i])> fabs(temp)) temp = array[i];
    sum+= array[i];
    sum2+= (array[i]*array[i]);
  }
  mean= sum/12;
  max=temp;
  rms =  sqrt(sum2/(144)-mean*mean);
  par[0]=max;
  par[1]=mean;
  par[2]=rms;
  std::cout<< "mean: "<<mean<<" rms: "<<rms<<" max: "<<max<<std::endl;
  return par;
}

 
void computeQuantiles(double* y, TH1F* hist, double band_lim){

  int max_bin = hist->GetMaximumBin();

  y[0]=hist->GetBinCenter(max_bin);
  
  double integral = hist->Integral();
  
  int i_min=max_bin; 
  int i_max=max_bin;
      
  int n_move=0;
  
  while(hist->Integral(i_min,i_max)<band_lim*integral){
    
    if(n_move%2==1)  i_min=std::max(--i_min,1); //i_min--; 
    else  i_max=std::min(++i_max,hist->GetNbinsX()-1);
    n_move++;
  }
  
  if(hist->Integral(i_min,i_max)/integral > band_lim*1.05){
    std::cerr << "[WARNING] " << std::endl;
    
    std::cout << i_min << "\t" << i_max << "\t" << hist->Integral(i_min,i_max)/integral << std::endl;
  }
  y[1]=fabs(hist->GetBinCenter(i_min)-y[0]);
  y[2]=fabs(hist->GetBinCenter(i_max)-y[0]);

}









//-------------------------------- N Wind StuDY ------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//





void  MakeSinglePlotBias( Int_t mass ,  Int_t width,Int_t cat, Int_t winSize, std::string genFcn, std::string fitFcn, double* par) {
  int ngen=0;
  ngen=0;
  TPaveText* label_cms = get_labelCMS(0,  true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TCanvas* ctmp = new TCanvas("ctmp", "ctmp", 1);
  ctmp->cd();
  std::string range;
  //mass=450
  //     |130|200|250|300|
  //     ________________
  // 1000| 0 | 1 | 2 | 3 |
  //     ________________  
  //  800| 4 | 5 | 6 | 7 |
  //     ________________
  //  700| 8 | 9 |10 |11 |
  //     ________________
  //  600|12 |13 |14 |15 |

  if(winSize==0) range= "130_1000";
  if(winSize==1) range= "200_1000";
  if(winSize==2) range= "250_1000";
  if(winSize==3) range= "300_1000";
  if(winSize==4) range= "130_800";
  if(winSize==5) range= "200_800";
  if(winSize==6) range= "250_800";
  if(winSize==7) range= "300_800";
  if(winSize==8) range= "130_700";
  if(winSize==9) range= "200_700";
  if(winSize==10) range= "250_700";
  if(winSize==11) range= "300_700";
  if(winSize==12) range= "130_600";
  if(winSize==13) range= "200_600";
  if(winSize==14) range= "250_600";
  if(winSize==15) range= "300_600";
  //mass=250
  //     |130 |150 |200 |
  //     ________________
  //  550| 16 | 17 | 18 |
  //     ________________  
  //  450| 19 | 20 | 21 |
  //     ________________
  //  350| 22 | 23 | 24 |
  if(winSize==16) range= "130_550";
  if(winSize==17) range= "150_550";
  if(winSize==18) range= "200_500";
  if(winSize==19) range= "130_450";
  if(winSize==20) range= "150_450";
  if(winSize==21) range= "200_450";
  if(winSize==22) range= "130_350";
  if(winSize==23) range= "150_350";
  if(winSize==24) range= "200_350";
  //mass=650
  //     |350 |450 |550 |
  //     ________________
  //  1000| 25 | 26 | 27 |
  //     ________________  
  //  900| 28 | 29 | 30 |
  //     ________________
  //  800| 31 | 32 | 33 |
  if(winSize==25) range= "350_1000";
  if(winSize==26) range= "450_1000";
  if(winSize==27) range= "550_1000";
  if(winSize==28) range= "350_900";
  if(winSize==29) range= "450_900";
  if(winSize==30) range= "550_900";
  if(winSize==31) range= "350_800";
  if(winSize==32) range= "450_800";
  if(winSize==33) range= "550_800";

  //mass=150
  //     |350 |
  //     _____
  //  250| 34 |
  //     _____  
  //  230| 35 | 
  //     _____
  //  200| 36 | 
  if(winSize==34) range= "130_250";
  if(winSize==35) range= "130_230";
  if(winSize==36) range= "130_200";

  //mass=850 || mass==750
  if(winSize==37) range= "450_1000";
  if(winSize==38) range= "550_1000";
  if(winSize==39) range= "650_1000";



  //mass=350
  //     |200 |250 |300 |
  //     ________________
  //  650| 40 | 41 | 42 |
  //     ________________  
  //  550| 43 | 44 | 45 |
  //     ________________
  //  450| 46 | 47 | 48 |
  if(winSize==40) range= "200_650";
  if(winSize==41) range= "250_650";
  if(winSize==42) range= "300_650";
  if(winSize==43) range= "200_550";
  if(winSize==44) range= "250_550";
  if(winSize==45) range= "300_550";
  if(winSize==46) range= "200_450";
  if(winSize==47) range= "250_450";
  if(winSize==48) range= "300_450";

  //mass=550
  //     |300 |400 |500 |
  //     ________________
  //  1000| 49 | 50 | 51 |
  //     ________________  
  //  800| 52 | 53 | 54 |
  //     ________________
  //  650| 55 | 56 | 57 |
  if(winSize==49) range= "300_1000";
  if(winSize==50) range= "400_1000";
  if(winSize==51) range= "500_1000";
  if(winSize==52) range= "300_800";
  if(winSize==53) range= "400_800";
  if(winSize==54) range= "500_800";
  if(winSize==55) range= "300_650";
  if(winSize==56) range= "400_650";
  if(winSize==57) range= "500_650";

  if(winSize==100)range="FIX";


  TFile* fin= new TFile(Form("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/biasCheck-sigGen_%d-g%s_f%s_m%d_w%d_cat%d_w10_%s.root", ngen,genFcn.c_str(), fitFcn.c_str(),mass,width, cat, range.c_str() ));
  std::cout<<Form("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/biasCheck-sigGen_%d-g%s_f%s_m%d_w%d_cat%d_w10_%s.root", ngen,genFcn.c_str(), fitFcn.c_str(),mass,width, cat, range.c_str() )<<std::endl;
  std::cout<<"after file"<<std::endl;
  
  TTree* tree = (TTree*) fin->Get("muTree");
  //histo Nsig
  std::cout<<"after tree"<<std::endl;
  TH1F* h_nsig = new TH1F("h_nsig", "nsig", 320, -300., 300.);
  h_nsig->Sumw2();
  std::cout<<"after create histo"<<std::endl;
  tree->Draw("sigYieldTot-sigYieldTotTruth>>h_nsig", "sigYieldTotErr>0.00001");
  std::cout<<"after fill histo"<<std::endl;
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_nsig->Draw("hist");
  h_nsig->GetXaxis()->SetTitle("N_{fit}-N_{true}");
  ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"NSig_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_%s.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat, range.c_str()));


  //histo errAss
  TH1F* h_sigerr = new TH1F("h_sigerr", "sigerr", 350, -1., 300.);
  h_sigerr->Sumw2();
  tree->Draw("sigYieldTotErr>>h_sigerr", "sigYieldTotErr>0.00001");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_sigerr->Draw("hist");
  h_sigerr->GetXaxis()->SetTitle("#Delta N_{fit}");
   ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"Err_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_%s.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat, range.c_str()));


  //histo errRel
  TH1F* h_sigerrRel = new TH1F("h_sigerrRel", "sigerrRel", 120, -10., 10.);
  h_sigerr->Sumw2();
  tree->Draw("sigYieldTotErr/sigYieldTot>>h_sigerrRel", "sigYieldTotErr>0.00001");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_sigerrRel->Draw("hist");
  h_sigerrRel->GetXaxis()->SetTitle("#Delta N_{fit}/N_{fit}");
  //  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"ErrRel_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_w%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat, winSize));
  //ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"ErrRel_"+TString::Format("g%s_f%s_m%d_w%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat));


  //histo frac
  TH1F* h_frac = new TH1F("h_frac", "frac", 120, 0., 1.);
  h_frac->Sumw2();
  if(fitFcn=="ExpolPL")tree->Draw("fracExpolPL>>h_frac");
  else if (fitFcn=="DiJetPL") tree->Draw("fracDiJetPL>>h_frac");
  else if (fitFcn=="DiJetEXP") tree->Draw("fracDiJetEXP>>h_frac");
  else if (fitFcn=="DiJetEXPOL") tree->Draw("fracDiJetEXPOL>>h_frac");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_frac->Draw("hist");
  h_frac->GetXaxis()->SetTitle("%Comp 1 Bkg");
  // ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"Err_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_w%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat, winSize));
  // ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"frac_"+TString::Format("g%s_f%s_m%d_w%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat));
 


  //histo bias
  TH1F* h_b = new TH1F("h_b", "bias", 60, -6., 6.);
  h_b->Sumw2();
  tree->Draw("(sigYieldTot-sigYieldTotTruth)/sigYieldTotErr>>h_b", "sigYieldTotErr>0.00001 && abs((sigYieldTot-sigYieldTotTruth)/sigYieldTotErr)<5.5&&chi2>0.001&&fitStatus==0 && migradStatus==0&& abs(sigYieldTot)<999" );
  double integral = h_b->ComputeIntegral(); 
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_b->Draw("hist");
  h_b->GetXaxis()->SetTitle("(N_{fit}-N_{true})/#Delta N_{fit}");
  ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias"+TString::Format("g%s_f%s_m%d_w%d_cat%d_%s.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat, range.c_str()));
 

  std::cout<<"compute bias median"<<std::endl;
  double median;
  int nqm=1;
  double xqm[1];
  double yqm[1];
  xqm[0]=0.5;

  // computeQuantiles(yqm, h_b, 0.5);

  h_b->GetQuantiles(nqm,yqm,xqm);
  median =yqm[0];
  std::cout << "**** bias median **** " << median << std::endl;
  par[0]=median;
    
  std::cout<<"compute mean and sigma"<<std::endl;
  double mean;
  double sigmaD;
  double sigmaU;

   TF1 *g = new TF1("g","gaus",-2,2);
  h_b->Fit("g");
  mean=h_b->GetFunction("g")->GetParameter(1);
  sigmaD=h_b->GetFunction("g")->GetParameter(2)/sqrt(h_b->Integral());
  sigmaU=sigmaD;
  std::cout<<" **** mean ****" << mean <<"  + "<< sigmaU<< " - "<<sigmaD << std::endl;
  par[1]=mean;
  par[2]=sigmaD;
  par[3]=sigmaU;

   std::cout<<"68% bands"<<std::endl;
  double mean68;
  double sigmaD68;
  double sigmaU68;

  int nq68=3;
  double xq68[3];
  double yq68[3];
  xq68[0]=0.5-0.68/2.;
  xq68[1]=0.5;
  xq68[2]=0.5+0.68/2;

  h_b->GetQuantiles(nq68,yq68,xq68);

  //computeQuantiles(yq68, h_b, 0.68);
  mean68=yq68[1];
  sigmaD68=fabs(yq68[0]-yq68[1]);
  sigmaU68=fabs(yq68[2]-yq68[1]);
  par[4]=mean68;
  par[5]=sigmaD68;
  par[6]=sigmaU68;


  std::cout<<"95% bands"<<std::endl;
  double mean95;
  double sigmaD95;
  double sigmaU95;

  int nq95=3;
  double xq95[3];
  double yq95[3];
  xq95[0]=0.5-0.95/2.;
  xq95[1]=0.5;
  xq95[2]=0.5+0.95/2;

  h_b->GetQuantiles(nq95,yq95,xq95);
  //  computeQuantiles(yq95, h_b,0.95 );
  mean95=yq95[1];
  sigmaD95=fabs(yq95[0]-yq95[1]);
  sigmaU95=fabs(yq95[2]-yq95[1]);

  par[7]=mean95;
  par[8]=sigmaD95;
  par[9]=sigmaU95;



  std::cout<<"compute nsig mean"<<std::endl;
  double nsig_mean;
  nsig_mean=h_nsig->GetMean();
  par[10]=nsig_mean;

  std::cout<<"compute nsig err mean"<<std::endl;
  double sigerr_mean;
  sigerr_mean=h_sigerr->GetMean();
  par[11]=sigerr_mean;

  std::cout<<"compute nsig errRel mean"<<std::endl;
  double sigerrRel_mean;
  sigerrRel_mean=h_sigerrRel->GetMean();
  par[12]=sigerrRel_mean;

  std::cout<<"compute frac mean"<<std::endl;
  double frac_mean;
  frac_mean=h_frac->GetMean();
  par[13]=frac_mean;


  std::cout<<"mass: "<<mass<<" cat: "<<cat<<std::endl;
  std::cout<<"median: "<<par[0]<<"  mean error: "<<par[11]<<std::endl;
  std::cout<<"95% D: "<<par[8]<<" 68% D: "<<par[5]<<std::endl;
  std::cout<<"68% U: "<<par[6]<<" 95% U: "<<par[9]<<std::endl;


 
}

 
 
void MakeGraphBias(std::string genFcn, std::string fitFcn,Int_t mass, Int_t width,Int_t cat,double* par_w1,double* par_w2,double* par_w3,double* par_w4,bool isVarMax, int var, double bias[4]) {

  //se la massa e' 250 gli ultimi due punti sono uguale poi taglio il range del plot.

  double win[4]={0, 1, 2, 3};
  double winErrU[4]={0., 0., 0., 0.};
  double winErrD[4]={0., 0., 0., 0.};


  //fill vectors
  double median[4] = {par_w1[0], par_w2[0], par_w3[0], par_w4[0]};
  for(int i=0; i<4;i++) bias[i]=median[i];
  double mu[4] = {par_w1[1], par_w2[1], par_w3[1], par_w4[1]};
  double sigmaD[4] = {par_w1[2], par_w2[2], par_w3[2], par_w4[2]};
  double sigmaU[4] = {par_w1[3], par_w2[3], par_w3[3], par_w4[3]};
  double m68[4] = {par_w1[4], par_w2[4], par_w3[4], par_w4[4]};
  double s68D[4] = {par_w1[5], par_w2[5], par_w3[5], par_w4[5]};
  double s68U[4] = {par_w1[6], par_w2[6], par_w3[6], par_w4[6]};
  double m95[4] = {par_w1[7], par_w2[7], par_w3[7], par_w4[7]};
  double s95D[4] = {par_w1[8], par_w2[8], par_w3[8], par_w4[8]};
  double s95U[4] = {par_w1[9], par_w2[9], par_w3[9], par_w4[9]};
  double nsig[4] = {par_w1[10], par_w2[10], par_w3[10], par_w4[10]};
  double sigerr[4] = {par_w1[11], par_w2[11], par_w3[11], par_w4[11]};
 

  int npt;
  if(mass==450)npt=4;
  else npt = 3;
  std::string srange;
  TGraphErrors* meadianGraph = new TGraphErrors(npt,win, median, winErrD, winErrU);
  TH1F* medianHisto= new TH1F("medianHisto", "medianHisto", 4, -0.5, 3.5 );
  if(mass==150){
  if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("130-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("130-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("130-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-200 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-230 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-250 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }else if(mass==450){
  if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("130-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("200-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("250-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(4, TString::Format("300-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(4, TString::Format("%d-600 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-700 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-800 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-1000 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }else if(mass==250){
 if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("130-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("150-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("200-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-350 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-450 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-550 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }else if(mass==650){
 if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("350-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("450-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("550-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-800 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-900 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-1000 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }else if(mass==850|| mass==750){
 if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("450-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("550-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("650-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-1000 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-1000 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-1000 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }else if(mass==350){
 if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("200-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("250-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("300-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-450 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-550 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-650 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }else if(mass==550){
 if(isVarMax){
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("300-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("400-%d GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("500-%d GeV",var));
    srange=TString::Format("MAX_%d",var);
  }else{
    medianHisto->GetXaxis()->SetBinLabel(3, TString::Format("%d-650 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(2, TString::Format("%d-800 GeV",var));
    medianHisto->GetXaxis()->SetBinLabel(1, TString::Format("%d-1000 GeV",var));
    srange=TString::Format("MIN_%d",var);
  }
  }

  TGraphAsymmErrors* m68Graph = new TGraphAsymmErrors(npt,win, m68, winErrD, winErrU, s68D, s68U);
  TGraphAsymmErrors* m95Graph = new TGraphAsymmErrors(npt,win, m95, winErrD, winErrU, s95D, s95U);

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  //meadianGraph->SetMarkerStyle(20);
  meadianGraph->SetMarkerColor(kBlue);
  meadianGraph->SetLineColor(kBlue);
  meadianGraph->SetLineWidth(2);
  meadianGraph->GetYaxis()->SetTitleFont(42);
  meadianGraph->GetXaxis()->SetTitleFont(42);
  meadianGraph->GetYaxis()->SetTitle("median((N_{fit}-N_{true})/#Delta N_{fit})");
      
  m68Graph->SetFillColor(kYellow);
	  
  m95Graph->SetFillColor(kGreen);
  m95Graph->GetYaxis()->SetRangeUser(-4.5,4.5);
  m95Graph->GetXaxis()->SetTitle("FitRange");
  m95Graph->GetYaxis()->SetTitleFont(42);
  m95Graph->GetXaxis()->SetTitleFont(42);
  m95Graph->GetYaxis()->SetTitle("(N_{fit}^{sig}-N_{true})/#Delta N_{fit}");

  TLegend* leg= new TLegend(0.6190805,0.7377622,0.8074713,0.9143357,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);	  
  //leg->SetTextSize(0.075);

  char line[300];
  
  leg->AddEntry(m95Graph,"95%","f");
  leg->AddEntry(m68Graph,"68%","f");
  leg->AddEntry(meadianGraph, "median","pl");
  medianHisto->SetMarkerSize(0);
  medianHisto->GetYaxis()->SetTitle("(N_{fit}^{sig}-N_{true})/#Delta N_{fit}");
  medianHisto->Draw("P");
  if(mass==250)medianHisto->GetXaxis()->SetRangeUser(0, 2);
  if(mass==450)medianHisto->GetXaxis()->SetRangeUser(0, 3);
  medianHisto->GetYaxis()->SetRangeUser(-2.5, 2.5);
 
  m95Graph->Draw("3same");
  
  m68Graph->Draw("3same");
  if(mass==250)m95Graph->GetXaxis()->SetRangeUser(0, 2);
  if(mass==450)m95Graph->GetXaxis()->SetRangeUser(0, 3);
  if(mass==250)m68Graph->GetXaxis()->SetRangeUser(0, 2);
  if(mass==450)m68Graph->GetXaxis()->SetRangeUser(0, 3);

  int Lmax;
  if(mass==450)Lmax=3;
  else Lmax=2;
  
  TLine* tlineU = new TLine(0,0.15,Lmax,0.15);
  tlineU->SetLineStyle(7);
  tlineU->SetLineWidth(2);
  tlineU->SetLineColor(kBlack);
 
  TLine* tlineD = new TLine(0,-0.15,Lmax,-0.15);
  tlineD->SetLineStyle(7);
  tlineD->SetLineWidth(2);
  tlineD->SetLineColor(kBlack);
  
 
 
  meadianGraph->Draw("PEL same");
  if(mass==250)meadianGraph->GetXaxis()->SetRangeUser(0, 2);
  if(mass==450)meadianGraph->GetXaxis()->SetRangeUser(0, 3);
  tlineU->Draw("Lsame");
  tlineD->Draw("Lsame");
  

  leg->Draw();      
  TPaveText *pt = new TPaveText(0.398851,0.7744755,0.5994253,0.9248252,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);	  
  pt->SetTextFont(42);

  sprintf(line,"Gen fun: %s", genFcn.c_str() );
  pt->AddText(line);
  sprintf(line,"Fit fun: %s",  fitFcn.c_str() );
  pt->AddText(line);
  sprintf(line,"Signal Mass: %d", mass );
  pt->AddText(line);
  sprintf(line,"Signal Width: %d", width  );
  pt->AddText(line);
  sprintf(line, "Category: %d", cat);
  pt->AddText(line);
  pt->Draw();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();

  // c->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias_vs_win_"+TString::Format("g%s_f%s_m%d_w%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat));
c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias_vs_win_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE_%s.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat, srange.c_str()));
  

  //make nsig and sigerr trend vs win
  TGraphErrors* nsigGraph = new TGraphErrors(4,win, nsig, winErrD, winErrU);
  TGraphErrors* sigerrGraph = new TGraphErrors(4,win, sigerr, winErrD, winErrU);
 

  nsigGraph->SetMarkerColor(kRed);
  nsigGraph->SetLineColor(kRed);
  nsigGraph->SetLineWidth(2);
  nsigGraph->GetYaxis()->SetTitleFont(42);
  nsigGraph->GetXaxis()->SetTitleFont(42);
  //nsigGraph->GetYaxis()->SetRangeUser(-10.,10.);
  nsigGraph->GetYaxis()->SetTitle("mean(N_{fit}^{sig})");
  nsigGraph->GetXaxis()->SetTitle("FitRange");
  nsigGraph->Draw("APXL");  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  pt->Draw("same");

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"nsig_vs_win_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE_%s.png", genFcn.c_str(), fitFcn.c_str(), mass,width,cat, srange.c_str()));
 

  sigerrGraph->SetMarkerColor(kSpring-6);
  sigerrGraph->SetLineColor(kSpring-6);
  sigerrGraph->SetLineWidth(2);
  sigerrGraph->GetYaxis()->SetTitleFont(42);
  sigerrGraph->GetXaxis()->SetTitleFont(42);
  nsigGraph->GetYaxis()->SetRangeUser(-10.,10.);
  sigerrGraph->GetYaxis()->SetTitle("mean(#Delta N_{fit}^{sig})");
  sigerrGraph->GetXaxis()->SetTitle("FitRange ");
  sigerrGraph->Draw("APXL");  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  pt->Draw("same");
  
  

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"sigerr_vs_win_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE_%s.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat, srange.c_str()));

 


}





void MakePlotBiasOneComb_varMAX(std::string genFcn, std::string fitFcn, int cat, Int_t width, Int_t mass, int var, double bias[4]){

  int r0,r1, r2, r3;
  if(mass==450){
 if(var==1000){
    r0=0;
    r1=1;
    r2=2;
    r3=3;}
  if(var==800){
    r0=4; 
    r1=5; 
    r2=6; 
    r3=7;}
  if(var==700){
    r0=8; 
    r1=9;
    r2=10;
    r3=11;}
  if(var==600){
    r0=12; 
    r1=13; 
    r2=14;
    r3=15;}
  }else if(mass==250){
    if(var==550){
      r0=16;
      r1=17;
      r2=18;
      r3=18;
    }
    if(var==450){
      r0=19;
      r1=20;
      r2=21;
      r3=21;
    }
    if(var==350){
      r0=22;
      r1=23;
      r2=24;
      r3=24;
    }
  }else if(mass==650){
    if(var==1000){
      r0=25;
      r1=26;
      r2=27;
      r3=27;
    }
    if(var==900){
      r0=28;
      r1=29;
      r2=30;
      r3=30;
    }
    if(var==800){
      r0=31;
      r1=32;
      r2=33;
      r3=33;
    }
  }else if(mass==850|| mass==750){
    if(var==1000){
      r0=37;
      r1=38;
      r2=39;
      r3=39;
    }
  }else if(mass==350){
    if(var==650){
      r0=40;
      r1=41;
      r2=42;
      r3=42;
    }
    if(var==550){
      r0=43;
      r1=44;
      r2=45;
      r3=45;
    }
    if(var==450){
      r0=46;
      r1=47;
      r2=48;
      r3=48;
    }
  }else if(mass==550){
    if(var==1000){
      r0=49;
      r1=50;
      r2=51;
      r3=51;
    }
    if(var==800){
      r0=52;
      r1=53;
      r2=54;
      r3=54;
    }
    if(var==650){
      r0=55;
      r1=56;
      r2=57;
      r3=57;
    }
  }
  double par_Range0[14];
  MakeSinglePlotBias(mass, width,cat,r0, genFcn.c_str(), fitFcn.c_str(), par_Range0);
  double par_Range1[14];
  MakeSinglePlotBias(mass, width,cat,r1, genFcn.c_str(), fitFcn.c_str(), par_Range1);
  double par_Range2[14];
  MakeSinglePlotBias(mass, width,cat,r2, genFcn.c_str(), fitFcn.c_str(), par_Range2);
  double par_Range3[14];
  MakeSinglePlotBias(mass, width,cat,r3, genFcn.c_str(), fitFcn.c_str(), par_Range3);
  
  MakeGraphBias(genFcn.c_str(), fitFcn.c_str(),mass, width,cat, par_Range0, par_Range1, par_Range2, par_Range3, true, var , bias);

}



void MakePlotBiasOneComb_varMIN(std::string genFcn, std::string fitFcn, int cat, Int_t width, Int_t mass, int var, double bias[4]){

  int r0,r1, r2, r3;
  if(mass==150){
  if(var==130){
    r0=34;
    r1=35;
    r2=36;
    r3=36;}
 
} else if(mass==450){
  if(var==130){
    r0=0;
    r1=4;
    r2=8;
    r3=12;}
  if(var==200){
    r0=1; 
    r1=5; 
    r2=9; 
    r3=13;}
  if(var==250){
    r0=2; 
    r1=6;
    r2=10;
    r3=14;}
  if(var==300){
    r0=3; 
    r1=7; 
    r2=11;
    r3=15;}
}else if(mass==250){
    if(var==130){
      r0=16;
      r1=19;
      r2=22;
      r3=22;
    }
    if(var==150){
      r0=17;
      r1=20;
      r2=23;
      r3=23;
    }
    if(var==200){
      r0=18;
      r1=21;
      r2=24;
      r3=24;
    }
  }else if(mass==650){
    if(var==350){
      r0=25;
      r1=28;
      r2=31;
      r3=31;
    }
    if(var==450){
      r0=26;
      r1=29;
      r2=32;
      r3=32;
    }
    if(var==550){
      r0=27;
      r1=30;
      r2=33;
      r3=33;
    }
  }else if(mass==350){
    if(var==200){
      r0=40;
      r1=43;
      r2=46;
      r3=46;
    }
    if(var==250){
      r0=41;
      r1=44;
      r2=47;
      r3=47;
    }
    if(var==300){
      r0=42;
      r1=45;
      r2=48;
      r3=48;
    }
  }else if(mass==550){
    if(var==300){
      r0=49;
      r1=52;
      r2=55;
      r3=55;
    }
    if(var==400){
      r0=50;
      r1=53;
      r2=56;
      r3=56;
    }
    if(var==500){
      r0=51;
      r1=54;
      r2=57;
      r3=57;
    }
  }
  double par_Range0[14];
  MakeSinglePlotBias(mass, width,cat,r0, genFcn.c_str(), fitFcn.c_str(), par_Range0);
  double par_Range1[14];
  MakeSinglePlotBias(mass, width,cat,r1, genFcn.c_str(), fitFcn.c_str(), par_Range1);
  double par_Range2[14];
  MakeSinglePlotBias(mass, width,cat,r2, genFcn.c_str(), fitFcn.c_str(), par_Range2);
  double par_Range3[14];
  MakeSinglePlotBias(mass, width,cat,r3, genFcn.c_str(), fitFcn.c_str(), par_Range3);
 
  MakeGraphBias(genFcn.c_str(), fitFcn.c_str(),mass, width,cat, par_Range0, par_Range1, par_Range2, par_Range3, false, var , bias);

}


void MakeBiasGrid450(std::string genFcn, std::string fitFcn, int cat, Int_t width){

  Int_t mass=450;
  //mass==450 
  double bias_MAX1000[4];
  double bias_MAX800[4];
  double bias_MAX700[4];
  double bias_MAX600[4]; 
  //righe
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,1000,bias_MAX1000);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,800,bias_MAX800);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,700,bias_MAX700);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,600,bias_MAX600);


  for(int i=0; i<4; i++)  std::cout<<" MAX 1000 Bias: "<<bias_MAX1000[i]<<std::endl;
  for(int i=0; i<4; i++)  std::cout<<" MAX 800 Bias: "<<bias_MAX800[i]<<std::endl;
  for(int i=0; i<4; i++)  std::cout<<" MAX 700 Bias: "<<bias_MAX700[i]<<std::endl;
  for(int i=0; i<4; i++)  std::cout<<" MAX 600 Bias: "<<bias_MAX600[i]<<std::endl;


  TH2F* g2= new TH2F("g2", "g2", 4, -0.5, 3.5, 4, -0.5, 3.5);
  g2->SetBinContent(1, 1, bias_MAX1000[0]);
  g2->SetBinContent(2, 1, bias_MAX1000[1]);
  g2->SetBinContent(3, 1, bias_MAX1000[2]);
  g2->SetBinContent(4, 1, bias_MAX1000[3]);

  g2->SetBinContent(1, 2, bias_MAX800[0]);
  g2->SetBinContent(2, 2, bias_MAX800[1]);
  g2->SetBinContent(3, 2, bias_MAX800[2]);
  g2->SetBinContent(4, 2, bias_MAX800[3]);

  g2->SetBinContent(1, 3, bias_MAX700[0]);
  g2->SetBinContent(2, 3, bias_MAX700[1]);
  g2->SetBinContent(3, 3, bias_MAX700[2]);
  g2->SetBinContent(4, 3, bias_MAX700[3]);

  g2->SetBinContent(1, 4, bias_MAX600[0]);
  g2->SetBinContent(2, 4, bias_MAX600[1]);
  g2->SetBinContent(3, 4, bias_MAX600[2]);
  g2->SetBinContent(4, 4, bias_MAX600[3]);

  g2->GetXaxis()->SetBinLabel(1,"MIN 130 GeV");
  g2->GetXaxis()->SetBinLabel(2,"MIN 200 GeV");
  g2->GetXaxis()->SetBinLabel(3,"MIN 250 GeV");
  g2->GetXaxis()->SetBinLabel(4,"MIN 300 GeV");
  g2->GetYaxis()->SetBinLabel( 1,"MAX 1000 GeV");
  g2->GetYaxis()->SetBinLabel( 2,"MAX 800 GeV");
  g2->GetYaxis()->SetBinLabel( 3,"MAX 700 GeV");
  g2->GetYaxis()->SetBinLabel( 4,"MAX 600 GeV");

  g2->GetZaxis()->SetRangeUser(-1.4, 1.4);

 TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();

  g2->Draw("COLZTEXT");

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"2DBias_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat));

}





void MakeBiasGrid250(std::string genFcn, std::string fitFcn, int cat, Int_t width){

  Int_t mass=250;
  //mass==450 
  double bias_MAX550[3];
  double bias_MAX450[3];
  double bias_MAX350[3];
 
  //righe
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,550,bias_MAX550);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,450,bias_MAX450);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,350,bias_MAX350);



  for(int i=0; i<3; i++)  std::cout<<" MAX 550 Bias: "<<bias_MAX550[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 450 Bias: "<<bias_MAX450[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 350 Bias: "<<bias_MAX350[i]<<std::endl;



  TH2F* g2= new TH2F("g2", "g2", 3, -0.5, 2.5, 3, -0.5, 2.5);
  g2->SetBinContent(1, 1, bias_MAX550[0]);
  g2->SetBinContent(2, 1, bias_MAX550[1]);
  g2->SetBinContent(3, 1, bias_MAX550[2]);

  g2->SetBinContent(1, 2, bias_MAX450[0]);
  g2->SetBinContent(2, 2, bias_MAX450[1]);
  g2->SetBinContent(3, 2, bias_MAX450[2]);

  g2->SetBinContent(1, 3, bias_MAX350[0]);
  g2->SetBinContent(2, 3, bias_MAX350[1]);
  g2->SetBinContent(3, 3, bias_MAX350[2]);

  g2->GetXaxis()->SetBinLabel(1,"MIN 130 GeV");
  g2->GetXaxis()->SetBinLabel(2,"MIN 150 GeV");
  g2->GetXaxis()->SetBinLabel(3,"MIN 200 GeV");
  g2->GetYaxis()->SetBinLabel( 1,"MAX 550 GeV");
  g2->GetYaxis()->SetBinLabel( 2,"MAX 450 GeV");
  g2->GetYaxis()->SetBinLabel( 3,"MAX 350 GeV");

  g2->GetZaxis()->SetRangeUser(-1.4, 1.4);

 TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();

  g2->Draw("COLZTEXT");

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"2DBias_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat));

}





void MakeBiasGrid650(std::string genFcn, std::string fitFcn, int cat, Int_t width){

  Int_t mass=650;
  //mass==450 
  double bias_MAX1000[3];
  double bias_MAX900[3];
  double bias_MAX800[3];
 
  //righe
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,1000,bias_MAX1000);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,900,bias_MAX900);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,800,bias_MAX800);



  for(int i=0; i<3; i++)  std::cout<<" MAX 1000 Bias: "<<bias_MAX1000[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 900 Bias: "<<bias_MAX900[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 800 Bias: "<<bias_MAX800[i]<<std::endl;



  TH2F* g2= new TH2F("g2", "g2", 3, -0.5, 2.5, 3, -0.5, 2.5);
  g2->SetBinContent(1, 1, bias_MAX1000[0]);
  g2->SetBinContent(2, 1, bias_MAX1000[1]);
  g2->SetBinContent(3, 1, bias_MAX1000[2]);

  g2->SetBinContent(1, 2, bias_MAX900[0]);
  g2->SetBinContent(2, 2, bias_MAX900[1]);
  g2->SetBinContent(3, 2, bias_MAX900[2]);

  g2->SetBinContent(1, 3, bias_MAX800[0]);
  g2->SetBinContent(2, 3, bias_MAX800[1]);
  g2->SetBinContent(3, 3, bias_MAX800[2]);

  g2->GetXaxis()->SetBinLabel(1,"MIN 350 GeV");
  g2->GetXaxis()->SetBinLabel(2,"MIN 450 GeV");
  g2->GetXaxis()->SetBinLabel(3,"MIN 550 GeV");
  g2->GetYaxis()->SetBinLabel( 1,"MAX 1000 GeV");
  g2->GetYaxis()->SetBinLabel( 2,"MAX 900 GeV");
  g2->GetYaxis()->SetBinLabel( 3,"MAX 800 GeV");

  g2->GetZaxis()->SetRangeUser(-1.4, 1.4);

 TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();

  g2->Draw("COLZTEXT");

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"2DBias_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat));

}






void MakeBiasGrid350(std::string genFcn, std::string fitFcn, int cat, Int_t width){

  Int_t mass=350;
  //mass==450 
  double bias_MAX650[3];
  double bias_MAX550[3];
  double bias_MAX450[3];
 
  //righe
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,650,bias_MAX650);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,550,bias_MAX550);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,450,bias_MAX450);



  for(int i=0; i<3; i++)  std::cout<<" MAX 650 Bias: "<<bias_MAX650[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 550 Bias: "<<bias_MAX550[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 450 Bias: "<<bias_MAX450[i]<<std::endl;



  TH2F* g2= new TH2F("g2", "g2", 3, -0.5, 2.5, 3, -0.5, 2.5);
  g2->SetBinContent(1, 1, bias_MAX650[0]);
  g2->SetBinContent(2, 1, bias_MAX650[1]);
  g2->SetBinContent(3, 1, bias_MAX650[2]);

  g2->SetBinContent(1, 2, bias_MAX550[0]);
  g2->SetBinContent(2, 2, bias_MAX550[1]);
  g2->SetBinContent(3, 2, bias_MAX550[2]);

  g2->SetBinContent(1, 3, bias_MAX450[0]);
  g2->SetBinContent(2, 3, bias_MAX450[1]);
  g2->SetBinContent(3, 3, bias_MAX450[2]);

  g2->GetXaxis()->SetBinLabel(1,"MIN 200 GeV");
  g2->GetXaxis()->SetBinLabel(2,"MIN 250 GeV");
  g2->GetXaxis()->SetBinLabel(3,"MIN 300 GeV");
  g2->GetYaxis()->SetBinLabel( 1,"MAX 650 GeV");
  g2->GetYaxis()->SetBinLabel( 2,"MAX 550 GeV");
  g2->GetYaxis()->SetBinLabel( 3,"MAX 450 GeV");

  g2->GetZaxis()->SetRangeUser(-1.4, 1.4);

 TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();

  g2->Draw("COLZTEXT");

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"2DBias_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat));

}






void MakeBiasGrid550(std::string genFcn, std::string fitFcn, int cat, Int_t width){

  Int_t mass=550;
  //mass==450 
  double bias_MAX1000[3];
  double bias_MAX800[3];
  double bias_MAX650[3];
 
  //righe
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,1000,bias_MAX1000);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,800,bias_MAX800);
  MakePlotBiasOneComb_varMAX( genFcn, fitFcn,  cat,  width, mass,650,bias_MAX650);



  for(int i=0; i<3; i++)  std::cout<<" MAX 1000 Bias: "<<bias_MAX1000[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 800 Bias: "<<bias_MAX800[i]<<std::endl;
  for(int i=0; i<3; i++)  std::cout<<" MAX 650 Bias: "<<bias_MAX650[i]<<std::endl;



  TH2F* g2= new TH2F("g2", "g2", 3, -0.5, 2.5, 3, -0.5, 2.5);
  g2->SetBinContent(1, 1, bias_MAX1000[0]);
  g2->SetBinContent(2, 1, bias_MAX1000[1]);
  g2->SetBinContent(3, 1, bias_MAX1000[2]);

  g2->SetBinContent(1, 2, bias_MAX800[0]);
  g2->SetBinContent(2, 2, bias_MAX800[1]);
  g2->SetBinContent(3, 2, bias_MAX800[2]);

  g2->SetBinContent(1, 3, bias_MAX650[0]);
  g2->SetBinContent(2, 3, bias_MAX650[1]);
  g2->SetBinContent(3, 3, bias_MAX650[2]);

  g2->GetXaxis()->SetBinLabel(1,"MIN 300 GeV");
  g2->GetXaxis()->SetBinLabel(2,"MIN 400 GeV");
  g2->GetXaxis()->SetBinLabel(3,"MIN 500 GeV");
  g2->GetYaxis()->SetBinLabel( 1,"MAX 1000 GeV");
  g2->GetYaxis()->SetBinLabel( 2,"MAX 800 GeV");
  g2->GetYaxis()->SetBinLabel( 3,"MAX 650 GeV");

  g2->GetZaxis()->SetRangeUser(-1.4, 1.4);

 TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();

  g2->Draw("COLZTEXT");

  c->SaveAs("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/winSize10/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"2DBias_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_RANGE.png", genFcn.c_str(), fitFcn.c_str(),mass, width, cat));

}


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//





void  MakeSinglePlotBias_wind( Int_t mass ,  Int_t width,Int_t cat, Int_t winSize, std::string genFcn, std::string fitFcn, double* par) {
  int ngen=0;
  if(mass==150)ngen=57.274;
  if(mass==200)ngen=73.980;
  if(mass==250)ngen=86.126;
  if(mass==300)ngen=99.099;
  if(mass==350)ngen=107.403;
  if(mass>=400)ngen=115.707;
  ngen=1;
  TPaveText* label_cms = get_labelCMS(0,  true);
  TPaveText* label_sqrt = get_labelSqrt(0);
  TCanvas* ctmp = new TCanvas("ctmp", "ctmp", 1);
  ctmp->cd();
  TFile* fin= new TFile(Form("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/biasCheck-sigGen_%d-g%s_f%s_m%d_w%d_cat%d_w%d_FIX.root", ngen,genFcn.c_str(), fitFcn.c_str(),mass,width, cat, winSize ));
  std::cout<<Form("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/biasCheck-sigGen_%d-g%s_f%s_m%d_w%d_cat%d_w%d_FIX.root", ngen,genFcn.c_str(), fitFcn.c_str(),mass,width, cat, winSize )<<std::endl;
  std::cout<<"after file"<<std::endl;
  
  TTree* tree = (TTree*) fin->Get("muTree");
  //histo Nsig
  std::cout<<"after tree"<<std::endl;
  TH1F* h_nsig = new TH1F("h_nsig", "nsig", 320, -300., 300.);
  h_nsig->Sumw2();
  std::cout<<"after create histo"<<std::endl;
  tree->Draw("sigYieldTot-sigYieldTotTruth>>h_nsig", "sigYieldTotErr>0.00001");
  std::cout<<"after fill histo"<<std::endl;
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_nsig->Draw("hist");
  h_nsig->GetXaxis()->SetTitle("N_{fit}-N_{true}");
  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"NSig_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_w%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat, winSize));
  ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"NSig_"+TString::Format("sigGen_%d_g%s_f%s_m%d_w%d_cat%d.png", ngen,genFcn.c_str(), fitFcn.c_str(), mass,width, cat));


  //histo errAss
  TH1F* h_sigerr = new TH1F("h_sigerr", "sigerr", 350, -1., 300.);
  h_sigerr->Sumw2();
  tree->Draw("sigYieldTotErr>>h_sigerr", "sigYieldTotErr>0.00001");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_sigerr->Draw("hist");
  h_sigerr->GetXaxis()->SetTitle("#Delta N_{fit}");
  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"Err_"+TString::Format("g%s_f%s_m%d_w%d_cat%d_w%d.png", genFcn.c_str(), fitFcn.c_str(), mass, width,cat, winSize));
  ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"Err_"+TString::Format("sigGen_%d_g%s_f%s_m%d_w%d_cat%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), mass,width, cat));


  //histo errRel
  TH1F* h_sigerrRel = new TH1F("h_sigerrRel", "sigerrRel", 120, -10., 10.);
  h_sigerr->Sumw2();
  tree->Draw("sigYieldTotErr/sigYieldTot>>h_sigerrRel", "sigYieldTotErr>0.00001");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_sigerrRel->Draw("hist");
  h_sigerrRel->GetXaxis()->SetTitle("#Delta N_{fit}/N_{fit}");
  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"ErrRel_"+TString::Format("sigGen_%d_g%s_f%s_m%d_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), mass,width, cat, winSize));
  //ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"ErrRel_"+TString::Format("g%s_f%s_m%d_w%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat));


  //histo frac
  TH1F* h_frac = new TH1F("h_frac", "frac", 120, 0., 1.);
  h_frac->Sumw2();
  if(fitFcn=="ExpolPL")tree->Draw("fracExpolPL>>h_frac");
  else if (fitFcn=="DiJetPL") tree->Draw("fracDiJetPL>>h_frac");
  else if (fitFcn=="DiJetEXP") tree->Draw("fracDiJetEXP>>h_frac");
  else if (fitFcn=="DiJetEXPOL") tree->Draw("fracDiJetEXPOL>>h_frac");
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_frac->Draw("hist");
  h_frac->GetXaxis()->SetTitle("%Comp 1 Bkg");
  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"Err_"+TString::Format("sigGen_%d_g%s_f%s_m%d_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), mass,width, cat, winSize));
  // ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"frac_"+TString::Format("g%s_f%s_m%d_w%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat));
 


  //histo bias
  TH1F* h_b = new TH1F("h_b", "bias", 60, -6., 6.);
  h_b->Sumw2();

  tree->Draw("(sigYieldTot-sigYieldTotTruth)/sigYieldTotErr>>h_b", "sigYieldTotErr>0.00001 && abs((sigYieldTot-sigYieldTotTruth)/sigYieldTotErr)<5.5&&chi2>0.001&&fitStatus==0 && migradStatus==0&& abs(sigYieldTot)<9999" );
  double integral = h_b->ComputeIntegral(); 
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  h_b->Draw("hist");
  h_b->GetXaxis()->SetTitle("(N_{fit}-N_{true})/#Delta N_{fit}");
  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias"+TString::Format("sigGen_%d_g%s_f%s_m%d_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), mass, width,cat, winSize));
  ctmp->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", winSize)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias"+TString::Format("g%s_f%s_m%d_w%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(), mass,width, cat));
 

  std::cout<<"compute bias median"<<std::endl;
  double median;
  int nqm=1;
  double xqm[1];
  double yqm[1];
  xqm[0]=0.5;

  // computeQuantiles(yqm, h_b, 0.5);

  h_b->GetQuantiles(nqm,yqm,xqm);
  median =yqm[0];
  std::cout << "**** bias median **** " << median << std::endl;
  par[0]=median;
    
  std::cout<<"compute mean and sigma"<<std::endl;
  double mean;
  double sigmaD;
  double sigmaU;

   TF1 *g = new TF1("g","gaus",-2,2);
  h_b->Fit("g");
  mean=h_b->GetFunction("g")->GetParameter(1);
  sigmaD=h_b->GetFunction("g")->GetParameter(2)/sqrt(h_b->Integral());
  sigmaU=sigmaD;
  std::cout<<" **** mean ****" << mean <<"  + "<< sigmaU<< " - "<<sigmaD << std::endl;
  par[1]=mean;
  par[2]=sigmaD;
  par[3]=sigmaU;

   std::cout<<"68% bands"<<std::endl;
  double mean68;
  double sigmaD68;
  double sigmaU68;

  int nq68=3;
  double xq68[3];
  double yq68[3];
  xq68[0]=0.5-0.68/2.;
  xq68[1]=0.5;
  xq68[2]=0.5+0.68/2;

  h_b->GetQuantiles(nq68,yq68,xq68);

  //computeQuantiles(yq68, h_b, 0.68);
  mean68=yq68[1];
  sigmaD68=fabs(yq68[0]-yq68[1]);
  sigmaU68=fabs(yq68[2]-yq68[1]);
  par[4]=mean68;
  par[5]=sigmaD68;
  par[6]=sigmaU68;


  std::cout<<"95% bands"<<std::endl;
  double mean95;
  double sigmaD95;
  double sigmaU95;

  int nq95=3;
  double xq95[3];
  double yq95[3];
  xq95[0]=0.5-0.95/2.;
  xq95[1]=0.5;
  xq95[2]=0.5+0.95/2;

  h_b->GetQuantiles(nq95,yq95,xq95);
  //  computeQuantiles(yq95, h_b,0.95 );
  mean95=yq95[1];
  sigmaD95=fabs(yq95[0]-yq95[1]);
  sigmaU95=fabs(yq95[2]-yq95[1]);

  par[7]=mean95;
  par[8]=sigmaD95;
  par[9]=sigmaU95;



  std::cout<<"compute nsig mean"<<std::endl;
  double nsig_mean;
  nsig_mean=h_nsig->GetMean();
  par[10]=nsig_mean;

  std::cout<<"compute nsig err mean"<<std::endl;
  double sigerr_mean;
  sigerr_mean=h_sigerr->GetMean();
  par[11]=sigerr_mean;

  std::cout<<"compute nsig errRel mean"<<std::endl;
  double sigerrRel_mean;
  sigerrRel_mean=h_sigerrRel->GetMean();
  par[12]=sigerrRel_mean;

  std::cout<<"compute frac mean"<<std::endl;
  double frac_mean;
  frac_mean=h_frac->GetMean();
  par[13]=frac_mean;


  std::cout<<"mass: "<<mass<<" cat: "<<cat<<std::endl;
  std::cout<<"median: "<<par[0]<<"  mean error: "<<par[11]<<std::endl;
  std::cout<<"95% D: "<<par[8]<<" 68% D: "<<par[5]<<std::endl;
  std::cout<<"68% U: "<<par[6]<<" 95% U: "<<par[9]<<std::endl;


 
}






//faccio il trend del bias con la larghezza 130-1000 per i miei  punti di massa

double MakeGraphBias_wind(std::string genFcn, std::string fitFcn,Int_t cat,Int_t width,Float_t percent,Int_t wind, double* par_m150,double* par_m250,double* par_m350,double* par_m450,double* par_m550,double* par_m650,double* par_m750,double* par_m850,double* par ) {

  double mass[8]={150,250,350,450,550,650,750,850.};
  double massErrU[8]={0., 0., 0.,0., 0., 0., 0.,0.};
  double massErrD[8]={0., 0., 0.,0., 0., 0., 0.,0.};


  //fill vectors
  double median[8] = {par_m150[0],par_m250[0],par_m350[0], par_m450[0], par_m550[0], par_m650[0],par_m750[0], par_m850[0]};
  double mu[8] = {par_m150[1],par_m250[1],par_m350[1], par_m450[1],  par_m550[1], par_m650[1],par_m750[1], par_m850[1]};
  double sigmaD[8] = {par_m150[2],par_m250[2],par_m350[2],par_m450[2],  par_m550[2], par_m650[2],par_m750[2], par_m850[2]};
  double sigmaU[8] = {par_m150[3],par_m250[3],par_m350[3],par_m450[3], par_m550[3], par_m650[3],par_m750[3], par_m850[3]};
  double m68[8] = {par_m150[4],par_m250[4],par_m350[4],par_m450[4],  par_m550[4], par_m650[4],par_m750[4], par_m850[4]};
  double s68D[8] = {par_m150[5],par_m250[5],par_m350[5],par_m450[5],  par_m550[5], par_m650[5],par_m750[5], par_m850[5]};
  float ntoy = 1000;
  double medianErrD[8] = {par_m150[5]/sqrt(ntoy),par_m250[5]/sqrt(ntoy),par_m350[5]/sqrt(ntoy),par_m450[5]/sqrt(ntoy),  par_m550[5]/sqrt(ntoy), par_m650[5]/sqrt(ntoy),par_m750[5]/sqrt(ntoy), par_m850[5]/sqrt(ntoy)};
  double s68U[8] = {par_m150[6],par_m250[6],par_m350[6], par_m450[6], par_m550[6], par_m650[6],par_m750[6], par_m850[6]};
  double medianErrU[8] = {par_m150[6]/sqrt(ntoy),par_m250[6]/sqrt(ntoy),par_m350[6]/sqrt(ntoy),par_m450[6]/sqrt(ntoy),  par_m550[6]/sqrt(ntoy), par_m650[6]/sqrt(ntoy),par_m750[6]/sqrt(ntoy), par_m850[6]/sqrt(ntoy)};
  double m95[8] = {par_m150[7],par_m250[7],par_m350[7], par_m450[7],  par_m550[7], par_m650[7],par_m750[7], par_m850[7]};
  double s95D[8] ={par_m150[8],par_m250[8],par_m350[8], par_m450[8],  par_m550[8], par_m650[8],par_m750[8], par_m850[8]};
  double s95U[8] ={par_m150[9],par_m250[9],par_m350[9], par_m450[9],  par_m550[9], par_m650[9],par_m750[9], par_m850[9]};
  double nsig[8] ={par_m150[10],par_m250[10],par_m350[10], par_m450[10],  par_m550[10], par_m650[10],par_m750[10], par_m850[10]};
  double sigerr[8] ={par_m150[11],par_m250[11],par_m350[11],par_m450[11],  par_m550[11], par_m650[11],par_m750[11], par_m850[11]};
  double sigerrRel[8] ={par_m150[12],par_m250[12],par_m350[12], par_m450[12],par_m550[12], par_m650[12],par_m750[12], par_m850[12]};
  double frac[8] = {par_m150[13],par_m250[13],par_m350[13],par_m450[13],  par_m550[13], par_m650[13],par_m750[13], par_m850[13]};

  for(int i = 0; i< 8; i++){

    std::cout<< mass[i] << "   "<< median[i]<<"   "<< width<<std::endl;

  }
  

  TGraphAsymmErrors* meadianGraph = new TGraphAsymmErrors(8,mass, median, massErrD, massErrU, medianErrD,medianErrU);
  TGraphAsymmErrors* m68Graph = new TGraphAsymmErrors(8,mass, m68, massErrD, massErrU, s68D, s68U);
  TGraphAsymmErrors* m95Graph = new TGraphAsymmErrors(8,mass, m95, massErrD, massErrU, s95D, s95U);

  TCanvas* c = new TCanvas("c", "c", 1);
  c->cd();

  //meadianGraph->SetMarkerStyle(20);
  meadianGraph->SetMarkerColor(kBlue);
  meadianGraph->SetLineColor(kBlue);
  meadianGraph->SetLineWidth(2);
  meadianGraph->GetYaxis()->SetTitleFont(42);
  meadianGraph->GetXaxis()->SetTitleFont(42);
  meadianGraph->GetYaxis()->SetTitle("median((N_{fit}-N_{true})/#Delta N_{fit})");
      
  m68Graph->SetFillColor(kYellow);
	  
  m95Graph->SetFillColor(kGreen);
  m95Graph->GetYaxis()->SetRangeUser(-2.5,2.5);
  m95Graph->GetXaxis()->SetTitle("m_{H} [GeV]");
  m95Graph->GetYaxis()->SetTitleFont(42);
  m95Graph->GetXaxis()->SetTitleFont(42);
  m95Graph->GetYaxis()->SetTitle("median((N_{fit}-N_{true})/#Delta N_{fit})");

  TLegend* leg= new TLegend(0.6190805,0.7377622,0.8074713,0.9143357,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColor(0);	  
  //leg->SetTextSize(0.075);

  char line[300];
  
  leg->AddEntry(m95Graph,"95%","f");
  leg->AddEntry(m68Graph,"68%","f");
  leg->AddEntry(meadianGraph, "median","pl");
  m95Graph->Draw("a3");
  m68Graph->Draw("3same");

  TLine* tlineU = new TLine(150.,0.15,850,0.15);
  tlineU->SetLineStyle(7);
  tlineU->SetLineWidth(2);
  tlineU->SetLineColor(kBlack);
 
  TLine* tlineD = new TLine(150.,-0.15,850.,-0.15);
  tlineD->SetLineStyle(7);
  tlineD->SetLineWidth(2);
  tlineD->SetLineColor(kBlack);
  
 

  meadianGraph->Draw("PEL same");
  tlineU->Draw("Lsame");
  tlineD->Draw("Lsame");
  

  leg->Draw();      
  TPaveText *pt = new TPaveText(0.398851,0.7744755,0.5994253,0.9248252,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillColor(0);	  
  pt->SetTextFont(42);

  sprintf(line,"Gen fun: %s", genFcn.c_str() );
  pt->AddText(line);
  sprintf(line,"Fit fun: %s",  fitFcn.c_str() );
  pt->AddText(line);
  if(width==0) sprintf(line, "Width: 0.1 GeV", percent);
  else sprintf(line, "Width: %.2f * M", percent);
  pt->AddText(line);
  sprintf(line, "Category: %d", cat);
  pt->AddText(line);
  pt->Draw();

  TPaveText* label_cms = get_labelCMS(0, true);
  TPaveText* label_sqrt = get_labelSqrt(0);  
  label_cms->Draw();
  label_sqrt->Draw();
 
  int ngen = 1;
  c->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias_vs_mass_"+TString::Format("sigGen_%dg%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), width,cat, wind));
    c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", wind)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"bias_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), width,cat, wind));

  //make nsig and sigerr trend vs mass
  TGraphErrors* nsigGraph = new TGraphErrors(8,mass, nsig, massErrD, massErrU);
  TGraphErrors* sigerrGraph = new TGraphErrors(8,mass, sigerr, massErrD, massErrU);
  TGraphErrors* sigerrRelGraph = new TGraphErrors(8,mass, sigerrRel, massErrD, massErrU);
  TGraphErrors* fracGraph = new TGraphErrors(8,mass, frac, massErrD, massErrU);

  nsigGraph->SetMarkerColor(kRed);
  nsigGraph->SetLineColor(kRed);
  nsigGraph->SetLineWidth(2);
  nsigGraph->GetYaxis()->SetTitleFont(42);
  nsigGraph->GetXaxis()->SetTitleFont(42);
  // nsigGraph->GetYaxis()->SetRangeUser(-4.,4.);
  nsigGraph->GetYaxis()->SetTitle("mean(N_{fit}-N_{true})");
  nsigGraph->GetXaxis()->SetTitle("m_{H} [GeV]");
  nsigGraph->Draw("APXL");  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  pt->Draw("same");


  c->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"nsig_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));

  c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", wind)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"nsig_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(), width,cat, wind));


  sigerrGraph->SetMarkerColor(kSpring-6);
  sigerrGraph->SetLineColor(kSpring-6);
  sigerrGraph->SetLineWidth(2);
  sigerrGraph->GetYaxis()->SetTitleFont(42);
  sigerrGraph->GetXaxis()->SetTitleFont(42);
  //  nsigGraph->GetYaxis()->SetRangeUser(-10.,10.);
  sigerrGraph->GetYaxis()->SetTitle("mean(#Delta N_{fit})");
  sigerrGraph->GetXaxis()->SetTitle("m_{H} [GeV]");
  sigerrGraph->Draw("APXL");  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  pt->Draw("same");

  
  c->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"sigerr_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));
  c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", wind)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"sigerr_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));


  sigerrRelGraph->SetMarkerColor(kAzure-2);
  sigerrRelGraph->SetLineColor(kAzure-2);
  sigerrRelGraph->SetLineWidth(2);
  sigerrRelGraph->GetYaxis()->SetTitleFont(42);
  sigerrRelGraph->GetXaxis()->SetTitleFont(42);
  //nsigGraph->GetYaxis()->SetRangeUser(-10.,10.);
  sigerrRelGraph->GetYaxis()->SetTitle("mean(#Delta N_{fit}/N_{fit})");
  sigerrRelGraph->GetXaxis()->SetTitle("m_{H} [GeV]");
  sigerrRelGraph->Draw("APXL");  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  pt->Draw("same");

 
  c->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"sigerrRel_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));
  c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", wind)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"sigerrRel_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));
 

  fracGraph->SetMarkerColor(kViolet-2);
  fracGraph->SetLineColor(kViolet-2);
  fracGraph->SetLineWidth(2);
  fracGraph->GetYaxis()->SetTitleFont(42);
  fracGraph->GetXaxis()->SetTitleFont(42);
  //fracGraph->GetYaxis()->SetRangeUser(0.94,1.01);
 
  fracGraph->GetYaxis()->SetTitle("mean( %)");
  fracGraph->GetXaxis()->SetTitle("m_{H} [GeV]");
  fracGraph->Draw("APXL");  
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  pt->Draw("same");

  
  c->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_6_1_1/src/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes_wind/"+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"frac_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));
  c->SaveAs(TString::Format("/afs/cern.ch/user/s/soffi/www/BiasRes/Winds/sliding9/", wind)+TString::Format("g%s_f%s/", genFcn.c_str(), fitFcn.c_str())+"frac_vs_mass_"+TString::Format("sigGen_%d_g%s_f%s_w%d_cat%d_w%d.png",ngen, genFcn.c_str(), fitFcn.c_str(),width, cat, wind));//sliding9 
 
     //return max  median over mass range

     getPar(median, par);
     double max;
     return max;
  

}




double* MakePlotBiasOneComb_vsMass_wind(std::string genFcn, std::string fitFcn, int cat, Int_t width, Int_t wind, Float_t percent, double* par){


  double par_m150[14];
  MakeSinglePlotBias_wind(150, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m150);
  double par_m250[14];
  MakeSinglePlotBias_wind(250, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m250);
  double par_m350[14];
  MakeSinglePlotBias_wind(350, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m350);
  double par_m450[14];
  MakeSinglePlotBias_wind(450, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m450);
  double par_m550[14];
  MakeSinglePlotBias_wind(550, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m550);
  double par_m650[14];
  MakeSinglePlotBias_wind(650, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m650);
  double par_m750[14];
  MakeSinglePlotBias_wind(750, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m750);
  double par_m850[14];
  MakeSinglePlotBias_wind(850, width,cat,wind, genFcn.c_str(), fitFcn.c_str(), par_m850);

  
  
  double maxbias= MakeGraphBias_wind( genFcn,  fitFcn,cat, width, percent, wind, par_m150,par_m250,par_m350, par_m450,par_m550, par_m650,  par_m750, par_m850,par);

 
}



void MakePlotBiasAllComb_vsMass_wind( Int_t cat, Int_t width, Int_t wind, Float_t percent){
  
  double maxRooKey3DiJetEXP[3];
  double maxRooKey4DiJetEXP[3];
  double maxExpolPLDiJetEXP[3];
  double maxDiJetEXPDiJetEXP[3];
  double maxDiJetEXPOLDiJetEXP[3];
  double maxDiJetPLDiJetEXP[3];
    
  MakePlotBiasOneComb_vsMass_wind("RooKey3","DiJetEXP" , cat, width, wind, percent,maxRooKey3DiJetEXP);
   MakePlotBiasOneComb_vsMass_wind("RooKey4","DiJetEXP" , cat, width, wind, percent,maxRooKey4DiJetEXP);
   MakePlotBiasOneComb_vsMass_wind("ExpolPL","DiJetEXP", cat, width, wind, percent,maxExpolPLDiJetEXP);
   MakePlotBiasOneComb_vsMass_wind("DiJetEXP","DiJetEXP", cat, width, wind, percent,maxDiJetEXPDiJetEXP);
   MakePlotBiasOneComb_vsMass_wind("DiJetEXPOL","DiJetEXP", cat, width, wind, percent,maxDiJetEXPOLDiJetEXP);
   MakePlotBiasOneComb_vsMass_wind("DiJetPL","DiJetEXP", cat, width, wind, percent,maxDiJetPLDiJetEXP);
 
 
  double maxRooKey3DiJetEXPOL[3];
  double maxRooKey4DiJetEXPOL[3];
  double maxExpolPLDiJetEXPOL[3];
  double maxDiJetEXPDiJetEXPOL[3];
  double maxDiJetEXPOLDiJetEXPOL[3];
  double maxDiJetPLDiJetEXPOL[3];
  
   MakePlotBiasOneComb_vsMass_wind("RooKey3","DiJetEXPOL" , cat, width, wind, percent,maxRooKey3DiJetEXPOL);
   MakePlotBiasOneComb_vsMass_wind("RooKey4","DiJetEXPOL" , cat, width, wind, percent,maxRooKey4DiJetEXPOL);
   MakePlotBiasOneComb_vsMass_wind("ExpolPL","DiJetEXPOL", cat, width, wind, percent,maxExpolPLDiJetEXPOL);
   MakePlotBiasOneComb_vsMass_wind("DiJetEXP","DiJetEXPOL", cat, width, wind, percent,maxDiJetEXPDiJetEXPOL);
   MakePlotBiasOneComb_vsMass_wind("DiJetEXPOL","DiJetEXPOL", cat, width, wind, percent,maxDiJetEXPOLDiJetEXPOL);
   MakePlotBiasOneComb_vsMass_wind("DiJetPL","DiJetEXPOL", cat, width, wind, percent,maxDiJetPLDiJetEXPOL);
 
  
  double maxRooKey3DiJetPL[3];
  double maxRooKey4DiJetPL[3];
  double maxExpolPLDiJetPL[3];
  double maxDiJetEXPDiJetPL[3];
  double maxDiJetEXPOLDiJetPL[3];
  double maxDiJetPLDiJetPL[3];
  /*
  MakePlotBiasOneComb_vsMass_wind("RooKey3","DiJetPL" , cat, width, wind, percent,maxRooKey3DiJetPL);
  MakePlotBiasOneComb_vsMass_wind("RooKey4","DiJetPL" , cat, width, wind, percent,maxRooKey4DiJetPL);
  MakePlotBiasOneComb_vsMass_wind("ExpolPL","DiJetPL", cat, width, wind, percent,maxExpolPLDiJetPL);
  MakePlotBiasOneComb_vsMass_wind("DiJetEXP","DiJetPL", cat, width, wind, percent,maxDiJetEXPDiJetPL);
  MakePlotBiasOneComb_vsMass_wind("DiJetEXPOL","DiJetPL", cat, width, wind, percent,maxDiJetEXPOLDiJetPL);
  MakePlotBiasOneComb_vsMass_wind("DiJetPL","DiJetPL", cat, width, wind, percent,maxDiJetPLDiJetPL);
  */
 /*double maxRooKey3ExpolPL[3];
 double maxRooKey4ExpolPL[3];
 double maxExpolPLExpolPL[3];
 double maxDiJetEXPExpolPL[3];
 double maxDiJetEXPOLExpolPL[3];
 double maxDiJetPLExpolPL[3];
 
 MakePlotBiasOneComb_vsMass_wind("RooKey3","ExpolPL" , cat, width, wind, percent,maxRooKey3ExpolPL);
 MakePlotBiasOneComb_vsMass_wind("RooKey4","ExpolPL" , cat, width, wind, percent,maxRooKey4ExpolPL);
 MakePlotBiasOneComb_vsMass_wind("ExpolPL","ExpolPL", cat, width, wind, percent,maxExpolPLExpolPL);
 MakePlotBiasOneComb_vsMass_wind("DiJetEXP","ExpolPL", cat, width, wind, percent,maxDiJetEXPExpolPL);
 MakePlotBiasOneComb_vsMass_wind("DiJetEXPOL","ExpolPL", cat, width, wind, percent,maxDiJetEXPOLExpolPL);
 MakePlotBiasOneComb_vsMass_wind("DiJetPL","ExpolPL", cat, width, wind, percent,maxDiJetPLExpolPL);





 std::cout<<"\\hline"<<std::endl; 
 std::cout<<TString::Format("Cat: %d, Width: %d", cat, width)+" & RooKey3  & RooKey4 & ExpolPL  & DiJetEXP  & DiJetEXPOL & DiJetPL \\\\"<<std::endl;
   std::cout<<"\\hline"<<std::endl; 
      std::cout<<" DiJetEXP & "<< toPrecision(maxRooKey3DiJetEXP[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey3DiJetEXP[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey3DiJetEXP[0],2);
  std::cout<<         " &  "<< toPrecision(maxRooKey4DiJetEXP[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey4DiJetEXP[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey4DiJetEXP[0],2);
  std::cout<<         " &  "<< toPrecision(maxExpolPLDiJetEXP[1],2)<<"$\\backslash$"<< toPrecision(maxExpolPLDiJetEXP[2],2)<<"$\\backslash$"<< toPrecision(maxExpolPLDiJetEXP[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPDiJetEXP[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPDiJetEXP[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPDiJetEXP[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPOLDiJetEXP[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLDiJetEXP[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLDiJetEXP[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetPLDiJetEXP[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLDiJetEXP[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLDiJetEXP[0],2);
  std::cout<< "\\\\ "<<std::endl;
  std::cout<<"\\hline"<<std::endl;

  std::cout<<" DiJetEXPOL & "<< toPrecision(maxRooKey3DiJetEXPOL[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey3DiJetEXPOL[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey3DiJetEXPOL[0],2);
  std::cout<<         " &  "<< toPrecision(maxRooKey4DiJetEXPOL[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey4DiJetEXPOL[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey4DiJetEXPOL[0],2);
  std::cout<<         " &  "<< toPrecision(maxExpolPLDiJetEXPOL[1],2)<<"$\\backslash$"<< toPrecision(maxExpolPLDiJetEXPOL[2],2)<<"$\\backslash$"<< toPrecision(maxExpolPLDiJetEXPOL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPDiJetEXPOL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPDiJetEXPOL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPDiJetEXPOL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPOLDiJetEXPOL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLDiJetEXPOL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLDiJetEXPOL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetPLDiJetEXPOL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLDiJetEXPOL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLDiJetEXPOL[0],2);
  std::cout<< "\\\\ "<<std::endl;
  std::cout<<"\\hline"<<std::endl;
   
  
  std::cout<<" DiJetPL & "<< toPrecision(maxRooKey3DiJetPL[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey3DiJetPL[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey3DiJetPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxRooKey4DiJetPL[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey4DiJetPL[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey4DiJetPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxExpolPLDiJetPL[1],2)<<"$\\backslash$"<< toPrecision(maxExpolPLDiJetPL[2],2)<<"$\\backslash$"<< toPrecision(maxExpolPLDiJetPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPDiJetPL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPDiJetPL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPDiJetPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPOLDiJetPL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLDiJetPL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLDiJetPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetPLDiJetPL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLDiJetPL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLDiJetPL[0],2);
  std::cout<< "\\\\ "<<std::endl;
  std::cout<<"\\hline"<<std::endl;*/
  /*
std::cout<<" ExpolPL & "<< toPrecision(maxRooKey3ExpolPL[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey3ExpolPL[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey3ExpolPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxRooKey4ExpolPL[1],2)<<"$\\backslash$"<< toPrecision(maxRooKey4ExpolPL[2],2)<<"$\\backslash$"<< toPrecision(maxRooKey4ExpolPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxExpolPLExpolPL[1],2)<<"$\\backslash$"<< toPrecision(maxExpolPLExpolPL[2],2)<<"$\\backslash$"<< toPrecision(maxExpolPLExpolPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPExpolPL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPExpolPL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPExpolPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetEXPOLExpolPL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLExpolPL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetEXPOLExpolPL[0],2);
  std::cout<<         " &  "<< toPrecision(maxDiJetPLExpolPL[1],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLExpolPL[2],2)<<"$\\backslash$"<< toPrecision(maxDiJetPLExpolPL[0],2);
  std::cout<< "\\\\ "<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  */

}



void MakeAllBiasPlots(){

  double par[4];
  int  width;
  for(int cat=0; cat<4;cat++ ){
    
    width = 0;
    MakePlotBiasOneComb_vsMass_wind("Expol", "ExpPAR", cat, width, 10, width, par);
    MakePlotBiasOneComb_vsMass_wind("ExpPAR", "ExpPAR", cat, width, 10, width, par);
    MakePlotBiasOneComb_vsMass_wind("DiJet", "ExpPAR", cat, width, 10, width, par);
    width = 10;
    MakePlotBiasOneComb_vsMass_wind("Expol", "ExpPAR", cat, width, 10, width, par);
    MakePlotBiasOneComb_vsMass_wind("ExpPAR", "ExpPAR", cat, width, 10, width, par);
    MakePlotBiasOneComb_vsMass_wind("DiJet", "ExpPAR", cat, width, 10, width, par);
   
  }

}
