#include "HighMass-HggFitter_mgg_Bias.cc"
//

void ProduceWorkspaces_Bias(int mass, Int_t perc, Int_t min, Int_t max){
  gSystem->SetIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.00-cms5/include");
  gSystem->Load("libRooFit");
  gROOT->ProcessLine(".x RooCBCrujffPdf.cxx+");
  cout << "mass = " << mass << endl; 
  // runfits(mass, true);
  runfits(mass, true, perc, min, max);
  
 
}
