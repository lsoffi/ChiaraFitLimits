
 

using namespace RooFit;
using namespace RooStats ;

static const Int_t NCAT = 4;  // chiara

std::string filePOSTfix="";
double signalScaler=1.00;

void  MakeSinglePlotBias( Int_t mass , Int_t cat, std::string genFcn, std::string fitFcn, double* par) ;


TPaveText* get_labelCMS( int legendQuadrant = 0 , std::string year="2012", bool sim=false) {

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
  label_sqrt->SetTextAlign(31); // align right
  label_sqrt->AddText("#sqrt{s} = 8 TeV");
  return label_sqrt;

}








void  MakeSinglePlotBias( Int_t mass , Int_t cat, std::string genFcn, std::string fitFcn, double* par) {

  TFile* fin= new TFile(Form("/afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/TestSWResults/mergeDir/biasCheck_g%s_f%s_m%d_cat%d.root", genFcn.c_str(), fitFcn.c_str(),mass, cat ));
  
  
  TTree* tree = fin->Get("muTree");

  TH1F* h_b = new TH1F("h_b", "bias", 60, -10., 10.);
  h_b->Sumw2();

  tree->Draw("(bkgTrue1fwhm-bkgSig1fwhm)/bkgErrSig1fwhm>>h_b");
 

  TPaveText* label_cms = get_labelCMS(0, "2012", false);
  TPaveText* label_sqrt = get_labelSqrt(0);
  label_cms->Draw("same");
  label_sqrt->Draw("same");
  
  TCanvas* ctmp = new TCanvas("ctmp", "ctmp", 1);
  ctmp->cd();
  h_b->Draw("hist");
  h_b->GetXaxis()->SetTitle("(N_{true}^{FWHM}-N_{fit}^{FWHM})/#Delta N_{fit}^{FWHM}");
  ctmp->SaveAs("/afs/cern.ch/work/s/soffi/CMSSW_5_2_3-globe/src/HiggsAnalysis/HiggsTo2photons/h2gglobe/ChiaraFitLimits/BiasStudy/BiasRes/bias"+TString::Format("g%s_f%s_m%d_cat%d.png", genFcn.c_str(), fitFcn.c_str(), mass, cat));


  //compute median
  double median;
  int nqm=1;
  double xqm[1];
  double yqm[1];
  xqm[0]=0.5;
  h_b->GetQuantiles(nqm,yqm,xqm);
  median =yqm[0];
  std::cout << "**** median **** " << median << std::endl;
  par[cat][0]=median;
    
  //compute mean and sigma
  double mean;
  doulbe sigmaD;
  doulbe sigmaU;

  h_b->Fit("gaus","Q");
  mean=h_b->GetFunction("gaus")->GetParameter(1);
  sigmaD=h_b->GetFunction("gaus")->GetParameter(2)/sqrt(h_b->Integral());
  sigmaU=sigmaD;
  std::cout<<" **** mean ****" << mean <<"  + "<< sigmaU<< " - "<<sigmaD << std::endl;
  par[cat][1]=mean;
  par[cat][2]=sigmaD;
  par[cat][3]=sigmaU;

  //68% bands
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
  mean68=yq68[1];
  sigmaD68=TMath::Abs(yq68[1]-yq68[0]);
  sigmaU68=TMath::Abs(yq68[2]-yq68[1]);
  par[cat][4]=mean68;
  par[cat][5]=sigmaD68;
  par[cat][6]=sigmaU68;


//95% bands
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
  mean95=yq95[1];
  sigmaD95=TMath::Abs(yq95[1]-yq95[0]);
  sigmaU95=TMath::Abs(yq95[2]-yq95[1]);
  par[cat][7]=mean95;
  par[cat][8]=sigmaD95;
  par[cat][9]=sigmaU95;

 
}


void MakeAllPlotBias(){

  double par[4][10];//ncatXnPar(mediana, mu, sigmaU, sigmaD, m68, s68D, s68U, m95, s95D, s95u)

  MakeSinglePlotBias(300, 0, "2exp", "2exp", *par);
  MakeSinglePlotBias(300, 1, "2exp", "2exp", *par);
  MakeSinglePlotBias(300, 2, "2exp", "2exp", *par);
  MakeSinglePlotBias(300, 3, "2exp", "2exp", *par);

}


