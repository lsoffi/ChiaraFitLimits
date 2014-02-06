PhotonsMass[0,500];
nvtx[0,50];
evweight[0, 1000];


PhotonsMass_sig_m0[400, 397.0, 403.0];
PhotonsMass_sig_sigma0[2., 0.0, 5.0];
PhotonsMass_sig_m1[400, 397.0, 403.0];
PhotonsMass_sig_sigma1[2.5, 0.0, 5.0];
PhotonsMass_sig_frac[0.5, 0.0, 1.0];

PhotonsMassGaussSig = Gaussian(PhotonsMass, PhotonsMass_sig_m0, PhotonsMass_sig_sigma0);
PhotonsMassGaussSig_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1, PhotonsMass_sig_sigma1);
PhotonsMassSig      = AddPdf(PhotonsMassGaussSig, PhotonsMassGaussSig_bis,PhotonsMass_sig_frac);



PhotonsMass_sig_m0_cat0[399.48, 399., 400.5];
PhotonsMass_sig_sigma0_cat0[4.1, 3.5, 4.5];
PhotonsMass_sig_m1_cat0[398.609, 398.5, 401.0];
PhotonsMass_sig_sigma1_cat0[10.2, 9.8, 10.7];
PhotonsMass_sig_frac_cat0[0.835, 0.5, 1.0];
nsig_cat0[50, 0, 100];


PhotonsMassGaussSig_cat0 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat0, PhotonsMass_sig_sigma0_cat0);
PhotonsMassGaussSig_cat0_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat0, PhotonsMass_sig_sigma1_cat0);
#PhotonsMassSig_cat0 = AddPdf(PhotonsMassGaussSig_cat0,PhotonsMassGaussSig_cat0_bis,PhotonsMass_sig_frac_cat0);
PhotonsMassSigExt_cat0 = RooExtendPdf(PhotonsMassSig_cat0, nsig_cat0);


PhotonsMass_sig_m0_cat1[399.31, 399., 400.5];
PhotonsMass_sig_sigma0_cat1[4.62,4.2, 5.0];
PhotonsMass_sig_m1_cat1[398.686, 398.5, 401.0];
PhotonsMass_sig_sigma1_cat1[12.7, 10., 14.0];
PhotonsMass_sig_frac_cat1[0.848, 0.7, 1.];
nsig_cat1[50, 0, 100];


PhotonsMassGaussSig_cat1 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat1, PhotonsMass_sig_sigma0_cat1);
PhotonsMassGaussSig_cat1_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat1, PhotonsMass_sig_sigma1_cat1);	
#PhotonsMassSig_cat1      = AddPdf(PhotonsMassGaussSig_cat1,PhotonsMassGaussSig_cat1_bis,PhotonsMass_sig_frac_cat1);
PhotonsMassSigExt_cat1 = RooExtendPdf(PhotonsMassSig_cat1, nsig_cat1);





PhotonsMass_sig_m0_cat2[399.995, 399., 400.5];
PhotonsMass_sig_sigma0_cat2[7.1, 6.7, 7.5];
PhotonsMass_sig_m1_cat2[399.496, 399., 400.5];
PhotonsMass_sig_sigma1_cat2[11.49, 11., 12.];
PhotonsMass_sig_frac_cat2[0.64, 0.5, 0.8];	
nsig_cat2[50, 0, 100];

PhotonsMassGaussSig_cat2 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat2, PhotonsMass_sig_sigma0_cat2);
PhotonsMassGaussSig_cat2_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat2, PhotonsMass_sig_sigma1_cat2);	
#PhotonsMassSig_cat2      = AddPdf(PhotonsMassGaussSig_cat2,PhotonsMassGaussSig_cat2_bis,PhotonsMass_sig_frac_cat2 );
PhotonsMassSigExt_cat2 = RooExtendPdf(PhotonsMassSig_cat2, nsig_cat2);



PhotonsMass_sig_m0_cat3[399.75, 399., 400.5];
PhotonsMass_sig_sigma0_cat3[10.58, 9.5, 11.7];
PhotonsMass_sig_m1_cat3[399.480, 399., 400.5];
PhotonsMass_sig_sigma1_cat3[6.3, 5.5, 7.5];
PhotonsMass_sig_frac_cat3[0.51, 0.5, 1.];
nsig_cat3[50, 0, 100];

PhotonsMassGaussSig_cat3 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat3, PhotonsMass_sig_sigma0_cat3);
PhotonsMassGaussSig_cat3_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m1_cat3, PhotonsMass_sig_sigma1_cat3);
#PhotonsMassSig_cat3      = AddPdf(PhotonsMassGaussSig_cat3,PhotonsMassGaussSig_cat3_bis,PhotonsMass_sig_frac_cat3);
PhotonsMassSigExt_cat3 = RooExtendPdf(PhotonsMassSig_cat3, nsig_cat3);











PhotonsMass_sig_alphaCB_cat0[1.67, 1., 2.];
PhotonsMass_sig_alphaC_cat0[0.12, 0.05, 0.2];
PhotonsMass_sig_mean_cat0[399.86, 399., 400.5];
PhotonsMass_sig_n_cat0[8.24, 7.2, 9.2];
PhotonsMass_sig_sigma_cat0[4.11, 3.,5.];


PhotonsMass_sig_alphaCB_cat1[1.61, 1., 2.];
PhotonsMass_sig_alphaC_cat1[0.14, 0.05, 0.2];
PhotonsMass_sig_mean_cat1[399.58, 399., 400.5];
PhotonsMass_sig_n_cat1[10.49, 8., 12];
PhotonsMass_sig_sigma_cat1[4.59, 3.5,5.5];


PhotonsMass_sig_alphaCB_cat2[1.47, 1., 2.];
PhotonsMass_sig_alphaC_cat2[0.082, 0.04, 0.15];
PhotonsMass_sig_mean_cat2[399.81, 399., 400.5];
PhotonsMass_sig_n_cat2[15.24, 10., 20];
PhotonsMass_sig_sigma_cat2[7.53, 6.5,8.5];


PhotonsMass_sig_alphaCB_cat3[1.45, 1., 2.];
PhotonsMass_sig_alphaC_cat3[0.09, 0.05, 0.15];
PhotonsMass_sig_mean_cat3[399.81, 399., 400.5];
PhotonsMass_sig_n_cat3[35.24, 10., 45.];
PhotonsMass_sig_sigma_cat3[7.39, 6.5,8.5];


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










#PhotonsMass_bkg_exp[-0.02,-0.05, 0.];
#PhotonsMassBkg = Exponential(PhotonsMass, PhotonsMass_bkg_exp);

#PhotonsMass_bkg_exp_cat0[-0.02,-0.05, 0.];
#PhotonsMassBkg_cat0 = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat0);

#PhotonsMass_bkg_exp_cat1[-0.025,-0.05, 0.];
#PhotonsMassBkg_cat1 = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat1);

#PhotonsMass_bkg_exp_cat2[-0.02,-0.05, 0.];
#PhotonsMassBkg_cat2 = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat2);

#PhotonsMass_bkg_exp_cat3[-0.025,-0.05, 0.];
#PhotonsMassBkg_cat3 = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat3);



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

PhotonsMass_bkg_8TeV_expol1_cat0[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_cat0[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_cat0[0.24, 0., 10.];  



PhotonsMass_bkg_8TeV_expol1_cat1[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_cat1[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_cat1[0.24, 0., 10.];  


PhotonsMass_bkg_8TeV_expol1_cat2[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_cat2[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_cat2[0.24, 0., 10.];  


PhotonsMass_bkg_8TeV_expol1_cat3[5.7,0 ., 40.];	  
PhotonsMass_bkg_8TeV_expol2_cat3[0.000002, 0.,20];
PhotonsMass_bkg_8TeV_expol3_cat3[0.24, 0., 10.];  



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



wei[1,0,10];

sqrtS[8000., 8000., 8000.]
