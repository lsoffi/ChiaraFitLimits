PhotonsMass[0,500];

PhotonsMass_sig_m0[300, 297.0, 303.0];
PhotonsMass_sig_sigma0[2., 0.0, 5.0];
PhotonsMass_sig_m1[300, 297.0, 303.0];
PhotonsMass_sig_sigma1[2.5, 0.0, 5.0];
PhotonsMass_sig_frac[0.5, 0., 1.0];

PhotonsMassGaussSig = Gaussian(PhotonsMass, PhotonsMass_sig_m0, PhotonsMass_sig_sigma0);
PhotonsMassGaussSig_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1, PhotonsMass_sig_sigma1);
PhotonsMassSig      = AddPdf(PhotonsMassGaussSig, PhotonsMassGaussSig_bis,PhotonsMass_sig_frac);




PhotonsMass_sig_m0_cat0[299.75, 299., 300.5];
PhotonsMass_sig_sigma0_cat0[3.2, 2.9, 3.8];
PhotonsMass_sig_m1_cat0[299.75, 299., 300.5];
PhotonsMass_sig_sigma1_cat0[9.7, 9., 11.];
PhotonsMass_sig_frac_cat0[0.89, 0.7, 1.0];

PhotonsMassGaussSig_cat0 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat0, PhotonsMass_sig_sigma0_cat0);
PhotonsMassGaussSig_cat0_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat0, PhotonsMass_sig_sigma1_cat0);
#PhotonsMassSig_cat0 = AddPdf(PhotonsMassGaussSig_cat0,PhotonsMassGaussSig_cat0_bis,PhotonsMass_sig_frac_cat0);


PhotonsMass_sig_m0_cat1[299.54, 299., 300.5];
PhotonsMass_sig_sigma0_cat1[3.79, 3.2, 4.5];
PhotonsMass_sig_m1_cat1[299.54, 299., 300.5];
PhotonsMass_sig_sigma1_cat1[11.97,  9., 13.];
PhotonsMass_sig_frac_cat1[0.7, 0.5, 1.];

PhotonsMassGaussSig_cat1 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat1, PhotonsMass_sig_sigma0_cat1);
PhotonsMassGaussSig_cat1_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat1, PhotonsMass_sig_sigma1_cat1);	
#PhotonsMassSig_cat1      = AddPdf(PhotonsMassGaussSig_cat1,PhotonsMassGaussSig_cat1_bis,PhotonsMass_sig_frac_cat1);





PhotonsMass_sig_m0_cat2[299.72, 299., 300.5];
PhotonsMass_sig_sigma0_cat2[5.5, 5., 6.];
PhotonsMass_sig_m1_cat2[299.496, 299., 300.5];
PhotonsMass_sig_sigma1_cat2[8.87, 8.2, 9.1];
PhotonsMass_sig_frac_cat2[0.59, 0.5,1.];	

PhotonsMassGaussSig_cat2 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat2, PhotonsMass_sig_sigma0_cat2);
PhotonsMassGaussSig_cat2_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat2, PhotonsMass_sig_sigma1_cat2);	
#PhotonsMassSig_cat2      = AddPdf(PhotonsMassGaussSig_cat2,PhotonsMassGaussSig_cat2_bis,PhotonsMass_sig_frac_cat2 );




PhotonsMass_sig_m0_cat3[299.53, 299., 300.5];
PhotonsMass_sig_sigma0_cat3[6.11, 5.5, 6.5];
PhotonsMass_sig_m1_cat3[299.53, 299., 300.5];
PhotonsMass_sig_sigma1_cat3[11.7, 11.5, 12.3];
PhotonsMass_sig_frac_cat3[0.86, 0.7, 1.];

PhotonsMassGaussSig_cat3 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat3, PhotonsMass_sig_sigma0_cat3);
PhotonsMassGaussSig_cat3_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat3, PhotonsMass_sig_sigma1_cat3);
#PhotonsMassSig_cat3      = AddPdf(PhotonsMassGaussSig_cat3,PhotonsMassGaussSig_cat3_bis,PhotonsMass_sig_frac_cat3);







nsig_cat0[50, 0, 100];
PhotonsMassSigExt_cat0 = RooExtendPdf(PhotonsMassSig_cat0, nsig_cat0);
nsig_cat1[50, 0, 100];
PhotonsMassSigExt_cat1 = RooExtendPdf(PhotonsMassSig_cat1, nsig_cat1);
nsig_cat2[50, 0, 100];
PhotonsMassSigExt_cat2 = RooExtendPdf(PhotonsMassSig_cat2, nsig_cat2);
nsig_cat3[50, 0, 100];
PhotonsMassSigExt_cat3 = RooExtendPdf(PhotonsMassSig_cat3, nsig_cat3);



PhotonsMass_sig_alphaCB_cat0[1.54, 0.9, 2.];
PhotonsMass_sig_alphaC_cat0[0.127, 0.08, 0.16];
PhotonsMass_sig_mean_cat0[299.81, 299., 300.5];
PhotonsMass_sig_n_cat0[10.2, 8., 12];
PhotonsMass_sig_sigma_cat0[3.12, 2.5,4.];


PhotonsMass_sig_alphaCB_cat1[1.66, 0.9, 2.];
PhotonsMass_sig_alphaC_cat1[0.131, 0.08, 0.16];
PhotonsMass_sig_mean_cat1[299.81, 299., 300.5];
PhotonsMass_sig_n_cat1[11.39, 8., 13];
PhotonsMass_sig_sigma_cat1[3.63, 2.5,4.5];


PhotonsMass_sig_alphaCB_cat2[1.55, 1.1, 1.7];
PhotonsMass_sig_alphaC_cat2[0.076, 0.05, 0.1];
PhotonsMass_sig_mean_cat2[299.81, 299., 300.5];
PhotonsMass_sig_n_cat2[29.24, 20., 33];
PhotonsMass_sig_sigma_cat2[6.14, 5.3,7.];


PhotonsMass_sig_alphaCB_cat3[1.64, 0.9, 2.];
PhotonsMass_sig_alphaC_cat3[0.09, 0.05, 0.16];
PhotonsMass_sig_mean_cat3[299.81, 299., 300.5];
PhotonsMass_sig_n_cat3[43.24, 25., 55.];
PhotonsMass_sig_sigma_cat3[6.1, 5.,7.];


ReducedMass_sig_alphaCB_cat0[2.1, 1., 3.];
ReducedMass_sig_alphaC_cat0[0.063, 0., 0.2];
ReducedMass_sig_mean_cat0[0., -0.1, 0.1];
ReducedMass_sig_n_cat0[4.59, 0., 10];
ReducedMass_sig_sigma_cat0[0.012 0.,3.];

ReducedMass_sig_alphaCB_cat1[2.23, 1., 3.5];
ReducedMass_sig_alphaC_cat1[0.0675, 0., 0.2];
ReducedMass_sig_mean_cat1[0., -0.1, 0.1];
ReducedMass_sig_n_cat1[4.26, 0., 10];
ReducedMass_sig_sigma_cat1[0.014, 0.,3.];

ReducedMass_sig_alphaCB_cat2[2.38, 1., 3.5];
ReducedMass_sig_alphaC_cat2[0.0611, 0., 0.2];
ReducedMass_sig_mean_cat2[0., -0.1, 0.1];
ReducedMass_sig_n_cat2[3.21, 0., 10];
ReducedMass_sig_sigma_cat2[0.021, 0.,3..];

ReducedMass_sig_alphaCB_cat3[2.35, 1., 3.5];
ReducedMass_sig_alphaC_cat3[0.064, 0., 0.2];
ReducedMass_sig_mean_cat3[0., -0.1, 0.1];
ReducedMass_sig_n_cat3[3.55, 0., 10];
ReducedMass_sig_sigma_cat3[0.022, 0.,3.];







nsigCBC_cat0[50, 0, 100];
nsigCBC_cat1[50, 0, 100];
nsigCBC_cat2[50, 0, 100];
nsigCBC_cat3[50, 0, 100];





wei[1,0,10];

PhotonsMass_bkg_exp[-0.026,-0.05, 0.];
PhotonsMassBkg = Exponential(PhotonsMass, PhotonsMass_bkg_exp);

PhotonsMass_bkg_exp_cat0[-0.018,-0.05, 0.];
PhotonsMassBkg_cat0_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat0);

PhotonsMass_bkg_exp_cat1[-0.016,-0.05, 0.];
PhotonsMassBkg_cat1_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat1);

PhotonsMass_bkg_exp_cat2[-0.015,-0.05, 0.];
PhotonsMassBkg_cat2_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat2);

PhotonsMass_bkg_exp_cat3[-0.013,-0.05, 0.];
PhotonsMassBkg_cat3_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat3);

PhotonsMass_bkg_8TeV_c0_cat0[4.75, 3., 8.];
PhotonsMass_bkg_8TeV_c1_cat0[0.48, -2.,2.];
PhotonsMass_bkg_8TeV_c2_cat0[0.93, 0.,2.];
PhotonsMass_bkg_8TeV_c3_cat0[0.32, 0.,2.];


PhotonsMass_bkg_8TeV_c0_cat1[5., 3., 8.];
PhotonsMass_bkg_8TeV_c1_cat1[0.47, -2.,2.];
PhotonsMass_bkg_8TeV_c2_cat1[0.8, -2.,2.];
PhotonsMass_bkg_8TeV_c3_cat1[0.32, 0.,2.];


PhotonsMass_bkg_8TeV_c0_cat2[4.84, 2., 10.];
PhotonsMass_bkg_8TeV_c1_cat2[0.48, 0.,3.];
PhotonsMass_bkg_8TeV_c2_cat2[0.93, 0.,3.];
PhotonsMass_bkg_8TeV_c3_cat2[0.32, 0.,3.];


PhotonsMass_bkg_8TeV_c0_cat3[4.84, 2., 10.];
PhotonsMass_bkg_8TeV_c1_cat3[0.48, -2.,3.];
PhotonsMass_bkg_8TeV_c2_cat3[0.93, 0.,3.];
PhotonsMass_bkg_8TeV_c3_cat3[0.32, 0.,3.];






PhotonsMass_bkg_8TeV_slope1_cat0[5.7,1., 20.];
PhotonsMass_bkg_8TeV_slope2_cat0[0.000002, 0., 0.00005];
PhotonsMass_bkg_8TeV_slope3_cat0[0.24, 0., 1.];
PhotonsMass_bkg_8TeV_norm_cat0[ 14739. ,1.e+04,1.e+05];


PhotonsMass_bkg_8TeV_slope1_cat1[9.3,1., 20.];
PhotonsMass_bkg_8TeV_slope2_cat1[0.000002, 0., 0.00005];
PhotonsMass_bkg_8TeV_slope3_cat1[0.6, 0., 1.];
PhotonsMass_bkg_8TeV_norm_cat1[ 25418. ,2.e+04,2.e+05];

PhotonsMass_bkg_8TeV_slope1_cat2[8.4,1., 20.];
PhotonsMass_bkg_8TeV_slope2_cat2[0.000002, 0., 0.00005];
PhotonsMass_bkg_8TeV_slope3_cat2[0.6, 0., 1.];
PhotonsMass_bkg_8TeV_norm_cat2[ 15416. ,1.e+04,1.e+05];

PhotonsMass_bkg_8TeV_slope1_cat3[9.44,3., 20.];
PhotonsMass_bkg_8TeV_slope2_cat3[0.000002, 0., 0.00005];
PhotonsMass_bkg_8TeV_slope3_cat3[0.24, 0., 1.];
PhotonsMass_bkg_8TeV_norm_cat3[ 34022 ,3.e+04,3.e+05];


PhotonsMass_bkg_8TeV_expol1_2_cat0[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat0[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_2_cat0[0.24, 0., 10.];  



PhotonsMass_bkg_8TeV_expol1_2_cat1[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat1[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_2_cat1[0.24, 0., 10.];  


PhotonsMass_bkg_8TeV_expol1_2_cat2[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat2[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_2_cat2[0.24, 0., 10.];  


PhotonsMass_bkg_8TeV_expol1_2_cat3[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat3[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_2_cat3[0.24, 0., 10.];  


sqrtS[8000., 8000., 8000.]
