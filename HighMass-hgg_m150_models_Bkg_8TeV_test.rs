PhotonsMass[0,500];
nvtx[0,50];
evweight[0, 1000];

PhotonsMass_sig_m0[150, 147.0, 153.0];
PhotonsMass_sig_sigma0[2., 0.0, 5.0];
PhotonsMass_sig_m1[150, 147.0, 153.0];
PhotonsMass_sig_sigma1[2.5, 0.0, 5.0];
PhotonsMass_sig_frac[0.5, 0.0, 1.0];

PhotonsMassGaussSig = Gaussian(PhotonsMass, PhotonsMass_sig_m0, PhotonsMass_sig_sigma0);
PhotonsMassGaussSig_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m0, PhotonsMass_sig_sigma1);
PhotonsMassSig      = AddPdf(PhotonsMassGaussSig, PhotonsMassGaussSig_bis,PhotonsMass_sig_frac);




PhotonsMass_sig_m0_cat0[149.93, 149., 150.5];
PhotonsMass_sig_sigma0_cat0[1.63, 1.3, 2.1];
PhotonsMass_sig_m1_cat0[149.93, 148.5, 151.0];
PhotonsMass_sig_sigma1_cat0[4.51, 3.8, 5.2];
PhotonsMass_sig_frac_cat0[0.83, 0., 1.0];


PhotonsMassGaussSig_cat0 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat0, PhotonsMass_sig_sigma0_cat0);
PhotonsMassGaussSig_cat0_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat0, PhotonsMass_sig_sigma1_cat0);
#PhotonsMassSig_cat0 = AddPdf(PhotonsMassGaussSig_cat0,PhotonsMassGaussSig_cat0_bis,PhotonsMass_sig_frac_cat0);


PhotonsMass_sig_m0_cat1[149.887, 149.0, 150.5];
PhotonsMass_sig_sigma0_cat1[1.954, 1.0, 3.0];
PhotonsMass_sig_m1_cat1[148.887, 148.0, 151.0];
PhotonsMass_sig_sigma1_cat1[5.78, 5.0, 6.5];
PhotonsMass_sig_frac_cat1[0.847, 0.5, 1.0];
	

PhotonsMassGaussSig_cat1 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat1, PhotonsMass_sig_sigma0_cat1);
PhotonsMassGaussSig_cat1_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat1, PhotonsMass_sig_sigma1_cat1);	
#PhotonsMassSig_cat1      = AddPdf(PhotonsMassGaussSig_cat1,PhotonsMassGaussSig_cat1_bis,PhotonsMass_sig_frac_cat1);




PhotonsMass_sig_m0_cat2[149.992, 149., 150.5];
PhotonsMass_sig_sigma0_cat2[3.4, 3., 3.8];
PhotonsMass_sig_m1_cat2[149.992, 149.7, 150.3];
PhotonsMass_sig_sigma1_cat2[7.3, 6.8, 7.8];
PhotonsMass_sig_frac_cat2[0.7, 0.5, 1.];
	

PhotonsMassGaussSig_cat2 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat2, PhotonsMass_sig_sigma0_cat2);
PhotonsMassGaussSig_cat2_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat2, PhotonsMass_sig_sigma1_cat2);	
#PhotonsMassSig_cat2      = AddPdf(PhotonsMassGaussSig_cat2,PhotonsMassGaussSig_cat2_bis,PhotonsMass_sig_frac_cat2 );




PhotonsMass_sig_m0_cat3[149.82, 149., 150.5];
PhotonsMass_sig_sigma0_cat3[3.48, 3., 4.];
PhotonsMass_sig_m1_cat3[149.82, 149., 150.5];
PhotonsMass_sig_sigma1_cat3[7.68, 7.2, 8.];
PhotonsMass_sig_frac_cat3[0.7, 0., 1.];


PhotonsMassGaussSig_cat3 = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat3, PhotonsMass_sig_sigma0_cat3);
PhotonsMassGaussSig_cat3_bis = Gaussian(PhotonsMass, PhotonsMass_sig_m0_cat3, PhotonsMass_sig_sigma1_cat3);
#PhotonsMassSig_cat3      = AddPdf(PhotonsMassGaussSig_cat3,PhotonsMassGaussSig_cat3_bis,PhotonsMass_sig_frac_cat3);

nsig_cat0[50, 0, 100];
PhotonsMassSigExt_cat0 = RooExtendPdf(PhotonsMassSig_cat0, nsig_cat0);
nsig_cat1[50, 0, 100];
PhotonsMassSigExt_cat1 = RooExtendPdf(PhotonsMassSig_cat1, nsig_cat1);
nsig_cat2[50, 0, 100];
PhotonsMassSigExt_cat2 = RooExtendPdf(PhotonsMassSig_cat2, nsig_cat2);
nsig_cat3[50, 0, 100];
PhotonsMassSigExt_cat3 = RooExtendPdf(PhotonsMassSig_cat3, nsig_cat3);





PhotonsMass_sig_alphaCB_cat0[1.50, 1.4, 1.6];
PhotonsMass_sig_alphaC_cat0[0.127, 0.11, 0.14];
PhotonsMass_sig_mean_cat0[149.81, 149., 150.5];
PhotonsMass_sig_n_cat0[7.4, 7., 7.8];
PhotonsMass_sig_sigma_cat0[1.62, 1.,2.];



PhotonsMass_sig_alphaCB_cat1[1.55, 1.3, 1.7];
PhotonsMass_sig_alphaC_cat1[0.13, 0.11, 0.14];
PhotonsMass_sig_mean_cat1[149.81, 149., 150.5];
PhotonsMass_sig_n_cat1[11.24, 9.5, 12.5];
PhotonsMass_sig_sigma_cat1[2.12, 1.2,2.8];


PhotonsMass_sig_alphaCB_cat2[1.57, 1.1, 1.9];
PhotonsMass_sig_alphaC_cat2[0.09, 0.05, 0.13];
PhotonsMass_sig_mean_cat2[149.81, 149., 150.5];
PhotonsMass_sig_n_cat2[43.24, 40., 70.];
PhotonsMass_sig_sigma_cat2[3.12, 2.2,4.];


PhotonsMass_sig_alphaCB_cat3[1.78, 1.6, 1.9];
PhotonsMass_sig_alphaC_cat3[0.09, 0.07, 0.11];
PhotonsMass_sig_mean_cat3[149.81, 149., 150.5];
PhotonsMass_sig_n_cat3[58.24, 50., 60.];
PhotonsMass_sig_sigma_cat3[3.32, 2.8,4.];



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


#PhotonsMass_bkg_exp[-0.014,-0.05, 0.];
#PhotonsMassBkg = Exponential(PhotonsMass, PhotonsMass_bkg_exp);
 
#PhotonsMass_bkg_exp_cat0[-0.026,-0.05, 0.];
#PhotonsMassBkg_cat0_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat0);
 
#PhotonsMass_bkg_exp_cat1[-0.028,-0.05, 0.];
#PhotonsMassBkg_cat1_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat1);
 
#PhotonsMass_bkg_exp_cat2[-0.023,-0.05, 0.];
#PhotonsMassBkg_cat2_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat2);
 
#PhotonsMass_bkg_exp_cat3[-0.025,-0.05, 0.];
#PhotonsMassBkg_cat3_exp = Exponential(PhotonsMass, PhotonsMass_bkg_exp_cat3);

PhotonsMass_bkg_8TeV_c0_cat0[4.75, 0., 20.];	
PhotonsMass_bkg_8TeV_c1_cat0[0.48, -20.,20.];	
PhotonsMass_bkg_8TeV_c2_cat0[0.93, 0.,20.];	
PhotonsMass_bkg_8TeV_c3_cat0[0.32, -20.,20.];	
PhotonsMass_bkg_8TeV_c4_cat0[0.32, 0.,20.];    	
PhotonsMass_bkg_8TeV_c5_cat0[0.32, -20.,20.];    		
PhotonsMass_bkg_8TeV_c6_cat0[0.32, 0.,20.];    		


PhotonsMass_bkg_8TeV_c0_cat1[4.75, -20., 20.];	
PhotonsMass_bkg_8TeV_c1_cat1[0.48, -20.,20.];	
PhotonsMass_bkg_8TeV_c2_cat1[0.93, -20.,20.];	
PhotonsMass_bkg_8TeV_c3_cat1[0.32, -20.,20.];	
PhotonsMass_bkg_8TeV_c4_cat1[0.32, -20.,20.]; 
PhotonsMass_bkg_8TeV_c5_cat1[0.32, -20.,20.];    		
PhotonsMass_bkg_8TeV_c6_cat1[0.32, 0.,20.];    		
   

PhotonsMass_bkg_8TeV_c0_cat2[4.75, -20., 20.];	
PhotonsMass_bkg_8TeV_c1_cat2[0.48, -20.,20.];	
PhotonsMass_bkg_8TeV_c2_cat2[0.93, -20.,20.];	
PhotonsMass_bkg_8TeV_c3_cat2[0.32, -20.,20.];	
PhotonsMass_bkg_8TeV_c4_cat2[0.32, -20.,20.];   
PhotonsMass_bkg_8TeV_c5_cat2[0.32, -20.,20.];    		
PhotonsMass_bkg_8TeV_c6_cat2[0.32, 0.,20.];    		


PhotonsMass_bkg_8TeV_c0_cat3[4.75, -20., 20.];	
PhotonsMass_bkg_8TeV_c1_cat3[0.48, -20.,20.];	
PhotonsMass_bkg_8TeV_c2_cat3[0.93, -20.,20.];	
PhotonsMass_bkg_8TeV_c3_cat3[0.32, -20.,20.];	
PhotonsMass_bkg_8TeV_c4_cat3[0.32, -20.,20.]; 
PhotonsMass_bkg_8TeV_c5_cat3[0.32, -20.,20.];    		
PhotonsMass_bkg_8TeV_c6_cat3[0.32, 0.,20.];    	

  






PhotonsMass_bkg_8TeV_slope1_cat0[5.7,1., 20.];
PhotonsMass_bkg_8TeV_slope2_cat0[0.000002, 0., 0.005];
PhotonsMass_bkg_8TeV_slope3_cat0[0.24, 0., 5.];
PhotonsMass_bkg_8TeV_norm_cat0[ 14739. ,1.e+04,1.e+05];


PhotonsMass_bkg_8TeV_slope1_cat1[9.3,1., 20.];
PhotonsMass_bkg_8TeV_slope2_cat1[0.000002, 0., 0.005];
PhotonsMass_bkg_8TeV_slope3_cat1[0.6, 0., 5.];
PhotonsMass_bkg_8TeV_norm_cat1[ 25418. ,2.e+04,2.e+05];

PhotonsMass_bkg_8TeV_slope1_cat2[8.4,1., 20.];
PhotonsMass_bkg_8TeV_slope2_cat2[0.000002, 0., 0.005];
PhotonsMass_bkg_8TeV_slope3_cat2[0.6, 0., 5.];
PhotonsMass_bkg_8TeV_norm_cat2[ 15416. ,1.e+04,1.e+05];

PhotonsMass_bkg_8TeV_slope1_cat3[9.44,3., 20.];
PhotonsMass_bkg_8TeV_slope2_cat3[0.000002, 0., 0.005];
PhotonsMass_bkg_8TeV_slope3_cat3[0.24, 0., 5.];
PhotonsMass_bkg_8TeV_norm_cat3[ 34022 ,3.e+04,3.e+05];



PhotonsMass_bkg_8TeV_expol1_2_cat0[5.7,0 ., 100.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat0[0.000002, 0.,50];
PhotonsMass_bkg_8TeV_expol3_2_cat0[0.24, 0., 30.];  


PhotonsMass_bkg_8TeV_expol1_2_cat1[5.7,0 ., 100.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat1[0.000002, 0.,50];
PhotonsMass_bkg_8TeV_expol3_2_cat1[0.24, 0., 30.];  


PhotonsMass_bkg_8TeV_expol1_2_cat2[5.7,0 ., 100.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat2[0.000002, 0.,50];
PhotonsMass_bkg_8TeV_expol3_2_cat2[0.24, 0., 30.];  


PhotonsMass_bkg_8TeV_expol1_2_cat3[5.7,0 .,100.];	  
PhotonsMass_bkg_8TeV_expol2_2_cat3[0.000002, 0.,50];
PhotonsMass_bkg_8TeV_expol3_2_cat3[0.24, 0., 30.];  



sqrtS[8000., 8000., 8000.]




