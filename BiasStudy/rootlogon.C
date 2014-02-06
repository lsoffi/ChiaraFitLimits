{

  cout << endl << "Welcome to my rootlogon.C" << endl;
  cout << "reading PhysicsTools/Utilities/macros/setTDRStyle.C" << endl;
  cout << "and some personal modifications." << endl << endl;

  gROOT->SetStyle("Plain");
//   setTDRStyle();
  gStyle->SetOptTitle(0);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.15); // per la paletta !
  gStyle->SetTitleXOffset(1.0);
  gStyle->SetTitleYOffset(1.3);
  gStyle->SetNdivisions(505, "X");
  gStyle->SetErrorX(0.5);
  //  gStyle->SetPalette(1,0);

///////// pretty palette ///////////

   const Int_t NRGBs = 5;
  const Int_t NCont = 255;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  gStyle->SetNumberContours(NCont);
  
/////////////////////////////////////

  gROOT->ForceStyle();
  //gROOT->SetStyle("tdrStyle");


  //  std::cout << "Set defaults" << std::endl;
  gSystem->Load("libPhysics");
  gSystem->Load("libRooFit") ;
  gSystem->AddIncludePath("-I/afs/cern.ch/cms/slc5_amd64_gcc462/lcg/roofit/5.32.03-cms13/include");

}
