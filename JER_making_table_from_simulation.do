#delimit ;

clear all ;
set memory 1000m;
clear mata;
clear matrix;




local p_gamma0=1;
local p_gamma1=1;
local  p_gamma2=-1;
local p_gamma3=-1;
local  p_beta0=1;
local p_beta1=-1;
local p_beta2=1;
local p_alpha=1;


foreach p_rho of numlist -2 -1 0 1 2{;
foreach p_mu of numlist  0 0.5 1  1.25 { ;

use  result_rho`p_rho'mu`p_mu'.dta;

drop  x1 u1 y1 y2 z1_hat x2_hat z1_hathat y1_star z2_hat ln_resid_sq correlation1;
drop if _n>=501;
*from here it is about the RV;

gen error_alpha_rv=alphahat_rv1-`p_alpha';
gen error_beta0_rv=beta0hat_rv1-`p_beta0';
gen error_beta1_rv=beta1hat_rv1-`p_beta1';
gen error_rho_rv=rhohat_rv1-`p_rho';


egen bias_alpha_rv=mean(error_alpha_rv);
egen bias_beta0_rv=mean(error_beta0_rv);
egen bias_beta1_rv=mean(error_beta1_rv);
egen bias_rho_rv=mean(error_rho_rv);


egen mse_alpha_rv=mean(error_alpha_rv^2);
egen mse_beta0_rv=mean(error_beta0_rv^2);
egen mse_beta1_rv=mean(error_beta1_rv^2);
egen mse_rho_rv=mean(error_rho_rv^2);

gen rmse_alpha_rv=sqrt(mse_alpha_rv);
gen rmse_beta0_rv=sqrt(mse_beta0_rv);
gen rmse_beta1_rv=sqrt(mse_beta1_rv);
gen rmse_rho_rv=sqrt(mse_rho_rv);

*from here it is about the gmm;

gen error_alpha_gmm2=alphahat_gmm2-`p_alpha';
gen error_beta0_gmm2=beta0hat_gmm2-`p_beta0';
gen error_beta1_gmm2=beta1hat_gmm2-`p_beta1';
gen error_rho_gmm2=rhohat_gmm2-`p_rho';




egen bias_alpha_gmm2=mean(error_alpha_gmm2);
egen bias_beta0_gmm2=mean(error_beta0_gmm2);
egen bias_beta1_gmm2=mean(error_beta1_gmm2);
egen bias_rho_gmm2=mean(error_rho_gmm2);


egen mse_alpha_gmm2=mean(error_alpha_gmm2^2);
egen mse_beta0_gmm2=mean(error_beta0_gmm2^2);
egen mse_beta1_gmm2=mean(error_beta1_gmm2^2);
egen mse_rho_gmm2=mean(error_rho_gmm2^2);

gen rmse_alpha_gmm2=sqrt(mse_alpha_gmm2);
gen rmse_beta0_gmm2=sqrt(mse_beta0_gmm2);
gen rmse_beta1_gmm2=sqrt(mse_beta1_gmm2);
gen rmse_rho_gmm2=sqrt(mse_rho_gmm2);

generate rho=`p_rho';
generate mu=`p_mu';

keep if _n==1;

keep rho mu bias_alpha_rv rmse_alpha_rv 
bias_beta0_rv rmse_beta0_rv 
bias_beta1_rv rmse_beta1_rv
bias_rho_rv rmse_rho_rv

bias_alpha_gmm2 rmse_alpha_gmm2 
bias_beta0_gmm2 rmse_beta0_gmm2 
bias_beta1_gmm2 rmse_beta1_gmm2
bias_rho_gmm2 rmse_rho_gmm2
;
save temp0.dta, replace;


clear;
use temp0.dta;
keep  rho mu bias_alpha_rv rmse_alpha_rv 
bias_beta0_rv rmse_beta0_rv 
bias_beta1_rv rmse_beta1_rv
bias_rho_rv rmse_rho_rv;

generate method="rv";

order rho mu method bias_alpha_rv rmse_alpha_rv 
bias_beta0_rv rmse_beta0_rv 
bias_beta1_rv rmse_beta1_rv
bias_rho_rv rmse_rho_rv;

rename bias_alpha_rv bias_alpha;
rename bias_beta0_rv bias_beta0;
rename bias_beta1_rv bias_beta1;
rename bias_rho_rv bias_rho;
rename rmse_alpha_rv rmse_alpha;
 rename rmse_beta0_rv rmse_beta0;
rename rmse_beta1_rv rmse_beta1;
rename rmse_rho_rv rmse_rho;

 
save temp1.dta, replace;
clear;
use temp0.dta;


keep  rho mu bias_alpha_gmm2 rmse_alpha_gmm2 
bias_beta0_gmm2 rmse_beta0_gmm2 
bias_beta1_gmm2 rmse_beta1_gmm2
bias_rho_gmm2 rmse_rho_gmm2;

generate method="gmm";

order rho mu method bias_alpha_gmm2 rmse_alpha_gmm2 
bias_beta0_gmm2 rmse_beta0_gmm2 
bias_beta1_gmm2 rmse_beta1_gmm2
bias_rho_gmm2 rmse_rho_gmm2;

rename bias_alpha_gmm2 bias_alpha;
rename bias_beta0_gmm2 bias_beta0;
rename bias_beta1_gmm2 bias_beta1;
rename bias_rho_gmm2 bias_rho;
rename rmse_alpha_gmm2 rmse_alpha;
 rename rmse_beta0_gmm2 rmse_beta0;
rename rmse_beta1_gmm2 rmse_beta1;
rename rmse_rho_gmm2 rmse_rho;



save temp2.dta, replace;
clear;

use temp1.dta;
append using temp2.dta;

save temp_rho`p_rho'mu`p_mu'.dta,replace;
};
};

clear;
use temp_rho-1mu0.dta;

foreach p_rho of numlist -2 -1 0 1 2{;
foreach p_mu of numlist  0 0.5 1  1.5 { ;
append using temp_rho`p_rho'mu`p_mu'.dta;
};
};
sort rho mu method;
;
outsheet using "E:\Data\Economics\stata\montecarlo\gmm_probit\JER_2015\table1.csv", comma replace;


