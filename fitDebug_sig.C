#include "RooPolynomial.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "TCanvas.h"
#include "RooProdPdf.h"
#include "RooExtendPdf.h"
#include "string"
#include "TFile.h"
#include "TH1F.h"
#include "RooAbsReal.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "TChain.h"
#include "TTree.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooRealConstant.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooFitResult.h"
#include "RooMCStudy.h"
#include "AsymmCBShapeVittorio.h"
#include "RooCrujffPdf.h"
#include "RooCBCrujffPdf.h"
#include "TLegend.h"
#include "RooWorkspace.h"

using namespace RooFit;

void fitDebug_sig(string cut, string filename){
 
  TChain mc("finalTree");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/G_Pt-120to170_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/G_Pt-170to300_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/G_Pt-30to50_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/G_Pt-50to80_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/G_Pt-80to120_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/QCDEM_Pt_170_250_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/QCDEM_Pt_20_30_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/QCDEM_Pt_250_350_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/QCDEM_Pt_30_80_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/QCDEM_Pt_350_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  mc.Add("root://pccmsrm27.cern.ch///cms/local/vtavolar/GammaJets/output_newPreselLooseIso2/QCDEM_Pt_80_170_8TeV_pythia6_hltcut30_hltiso0_mvaWP4.root");
  
  RooRealVar combinedPfIso03Phot("combinedPfIso03Phot", "combinedPfIso03Phot", -7., 15.);
  RooRealVar etaPhot("etaPhot", "etaPhot", -2.5, 2.5);
  RooRealVar mvaIdPhot("mvaIdPhot", "mvaIdPhot", -1.,1.);
  RooRealVar isMatchedPhot("isMatchedPhot","isMatchedPhot", -1., 2.);
  RooRealVar ptPhot("ptPhot", "ptPhot", 0., 1000.);
  RooRealVar weight("weight","weight", 0., 100.);

  RooArgSet argSet("argSet");

  argSet.add(combinedPfIso03Phot);
  argSet.add(etaPhot);
  argSet.add(mvaIdPhot);
  argSet.add(isMatchedPhot);
  argSet.add(ptPhot);
  argSet.add(weight);

  combinedPfIso03Phot.setBins(121);
  etaPhot.setBins(120);
  mvaIdPhot.setBins(180);
  isMatchedPhot.setBins(3);
  ptPhot.setBins(1200);
  weight.setBins(1000);

  std::cout<<"set binning"<<std::endl;


  RooDataSet allSet("allSet", "allSet", argSet, RooFit::WeightVar("weight"), RooFit::Import(mc));

  std::cout<<"created complete dataset"<<endl;

  std::cout<<allSet.GetName()<<std::endl;

  RooDataSet* d_s = (RooDataSet*)allSet.reduce((cut+" && mvaIdPhot>0.711099").c_str());
  std::cout<<"d_s entries: "<<d_s->sumEntries()<<std::endl;
  std::cout<<"created reduced dataset"<<std::endl;

  TH1F* h_set_s = (TH1F*)d_s->createHistogram("h_set_s", combinedPfIso03Phot, RooFit::Binning(combinedPfIso03Phot.getBinning()));
  std::cout<<"integral h_set_s"<<h_set_s->Integral()<<std::endl;


  RooDataHist dh_s("dh_s", "dh_s", combinedPfIso03Phot, Import(*h_set_s));
  std::cout<<dh_s.sum(kTRUE)<<std::endl;
  std::cout<<"created reduced datahist, scut"<<std::endl;
  
  RooRealVar crujffmean("crujffmean", "crujffmean", -1., -10., 0.);
  RooRealVar crujffsigmaL("crujffsigmaL", "crujffsigmaL", 1., 0., 10.);
  RooRealVar crujffsigmaR("crujffsigmaR", "crujffsigmaR", 1., 0., 10.);
  RooRealVar crujffalphaL("crujffalphaL", "crujffalphaL", 0.8, 0., 10.);
  RooRealVar crujffalphaR("crujffalphaR", "crujffalphaR", 0.8, 0., 10.);

  RooCrujffPdf crujff_s("crujff_s", "crujff_s", combinedPfIso03Phot, crujffmean, crujffsigmaL, crujffsigmaR, crujffalphaL, crujffalphaR);


  RooRealVar CBC_mean("CBC_mean", "CBC_mean", -1., -10., 0.);
  RooRealVar CBC_sigma("CBC_sigma", "CBC_sigma", 1., 0., 10.);
  RooRealVar CBC_alphaC("CBC_alphaC", "CBC_alphaC", 0.7, 0., 10.);
  RooRealVar CBC_alphaCB("CBC_alphaCB", "CBC_alphaCB", 0.7, 0., 10.);
  RooRealVar CBC_n("CBC_n", "CBC_n", 9., 0., 300.);

  RooCBCrujffPdf CBCrujff_s("CBCrujff_s", "CBCrujff_s", combinedPfIso03Phot, CBC_mean, CBC_sigma, CBC_alphaC, CBC_alphaCB, CBC_n);

  RooRealVar gaussmean("gaussmean","gaussmean", -1.33, -10., 0.);
  RooRealVar gausssigma("gausssigma", "gausssigma", 1.1, 0., 5.);

  RooGaussian my_gauss_s("my_gauss_s", "my_gauss_s", combinedPfIso03Phot, gaussmean, gausssigma);


  RooRealVar nsig("nsig", "expexted number of ev for sig", 250000., 10., 100000000.);

  RooRealVar nbg_s("nbg_s", "expexted number of ev for bg, scut", 250000., 10., 100000000.);

  RooRealVar cbmean("cbmean", "cbmean", -0.872, -10., 1.);
  RooRealVar cbsigma1("cbsigma1", "cbsigma1", 0.9, 0., 5.);
  RooRealVar cbsigma2("cbsigma2", "cbsigma2", 1.1, 1.1, 1.1);

  RooRealVar cbalpha_s("cbalpha_s", "cbaplha_s", -0.840, -10., 0.);
  RooRealVar cbn_s("cbn_s","cbn_s", 9., 0., 300.);

  RooCBShape my_cb_s("my_cb_s", "my_cb_s",  combinedPfIso03Phot, cbmean, cbsigma1, cbalpha_s, cbn_s);
  AsymmCBShapeVittorio my_asymm_cb("my_cb_s", "my_cb_s",  combinedPfIso03Phot, gaussmean, cbsigma1, gausssigma, cbalpha_s, cbn_s);

  RooRealVar frac_s("frac_s", "frac_s", 0.8, 0., 1.);

  //Gaussian to constrain fraction f in [f*CB + (1-f)*Gauss]
  
  RooGaussian frac_constraint("frac_constraint", "frac_constraint", frac_s, RooConst(0.977),RooConst(0.1));
 
  RooAddPdf my_add_s("my_add_s", "my_add_s", my_asymm_cb, my_gauss_s, frac_s);
  //  RooExtendPdf ext_pdf_s("ext_pdf_s", "ext_pdf_s", my_add_s, nbg_s);

  RooProdPdf my_add_s_constrained("my_add_s_constrained", "my_add_s_constrained", RooArgSet(my_add_s, frac_constraint));

  /*RooDataSet* d = my_add_s.generate(RooArgSet(combinedPfIso03Phot), 100000);
    RooDataHist* dh_s = d->binnedClone();*/
  RooFitResult* result =  CBCrujff_s.fitTo(*d_s/*, Constrain(frac_s)*/, Save(), Range(-5.,15.), /* Extended(kTRUE), */ SumW2Error(kTRUE));

  RooPlot* frame_s = combinedPfIso03Phot.frame(RooFit::Title("Fit to combinedPfIso03Phot, scut region"));

  d_s->plotOn(frame_s, Name("dh_s"));
  CBCrujff_s.plotOn(frame_s,Name("pdf_s"),LineColor(kCyan)/*, Normalization(1.,RooAbsReal::RelativeExpected)*/) ;
  //my_asymm_cb.plotOn(frame_s, Name("my_cb"),RooFit::LineColor(kMagenta), Normalization(frac_s.getVal()));
  //  my_gauss_s.plotOn(frame_s, Name("my_gaussian"), RooFit::LineColor(kMagenta), Normalization((1-frac_s.getVal())));


  frame_s->SetMinimum(0.00001);
  
  TLegend* a = new TLegend(0.63,0.68, 0.88, 0.88);
  a->SetBorderSize(0);
  a->SetFillColor(0);
  a->SetFillStyle(0);
  a->SetTextSize(0.038);
  a->AddEntry(frame_s->findObject("dh_s"), "MC","p");
  a->AddEntry(frame_s->findObject("pdf_s"), "fitting PDF","l");
  //  a->AddEntry(frame_s->findObject("my_cb"), "PDF components","l");
  
  TCanvas* c = new TCanvas();
  c->SetTitle(frame_s->GetTitle());
  frame_s->Draw("");
  a->Draw();
  c->SaveAs((filename+"_s.png").c_str());
  c->SaveAs((filename+"_s.root").c_str());

  Double_t chi2 = frame_s->chiSquare("pdf_s", "dh_s", 8);

  std::cout<<"ChiSquared value, scut: "<<chi2<<std::endl;

  c->SetLogy();
  c->SetTitle(frame_s->GetTitle());
  frame_s->Draw("");
  a->Draw();
  c->SaveAs((filename+"_s_log.png").c_str());
  c->SaveAs((filename+"_s_log.root").c_str());

  RooWorkspace* w_sig = new RooWorkspace("w_sig", "workspace");
  
  w_sig->import(dh_s);
  w_sig->import(CBCrujff_s);

  w_sig->Print();

  w_sig->writeToFile(("workspace_"+filename+".root").c_str());

  TFile* f_fitRes = new TFile(("fitResult_"+filename+".root").c_str(), "RECREATE");
  result->Write();
  //  f_fitRes->Write();
  f_fitRes->Close();
  

}
