* This is the simulation for heteroscascity


#delimit ;
clear all ;
set more off;
capture log close;
capture log using "mylog.smcl" ,replace; 
clear mata;
clear matrix;
local clear;
set memory 2000m;
set matsize 10000;
local n_obs 500;
*local p_rho=-2;
local n_replic=500;
local p_gamma0=1.0;
local p_gamma1=1.0;
local  p_gamma2=-1.0;
local p_gamma3=-1.0;
local p_gamma4=1.0;
local  p_beta0=1.0;
local p_beta1=-1.0;
local p_beta2=1.0;
local p_alpha=1.0;
*local lambda=1.25;
local opt_k 14;

timer on 1;
set obs `n_obs';



matrix alphahat_rv=J(`n_replic',1,.);
matrix beta1hat_rv=J(`n_replic',1,.);
matrix beta0hat_rv=J(`n_replic',1,.);
matrix rhohat_rv=J(`n_replic',1,.);
matrix alphahat_gmm1=J(`n_replic',1,.);
matrix beta1hat_gmm1=J(`n_replic',1,.);
matrix beta0hat_gmm1=J(`n_replic',1,.);
matrix rhohat_gmm1=J(`n_replic',1,.);
matrix gamma0hat_gmm1=J(`n_replic',1,.);
matrix gamma1hat_gmm1=J(`n_replic',1,.);
matrix gamma2hat_gmm1=J(`n_replic',1,.);


matrix alphahat_gmm2=J(`n_replic',1,.);
matrix beta1hat_gmm2=J(`n_replic',1,.);
matrix beta0hat_gmm2=J(`n_replic',1,.);
matrix rhohat_gmm2=J(`n_replic',1,.);
matrix alphahat_gmm3=J(`n_replic',1,.);
matrix beta1hat_gmm3=J(`n_replic',1,.);
matrix beta0hat_gmm3=J(`n_replic',1,.);

matrix beta0_stder_eqn=J(`n_replic',1,.);
matrix alpha_stder_eqn=J(`n_replic',1,.);
matrix beta1_stder_eqn=J(`n_replic',1,.);
matrix rho_stder_eqn=J(`n_replic',1,.);



matrix alpha_stder_gmm1=J(`n_replic',1,.);
matrix beta1_stder_gmm1=J(`n_replic',1,.);
matrix beta0_stder_gmm1=J(`n_replic',1,.);
matrix rho_stder_gmm1=J(`n_replic',1,.);
matrix gamma0_stder_gmm1=J(`n_replic',1,.);
matrix gamma1_stder_gmm1=J(`n_replic',1,.);
matrix gamma2_stder_gmm1=J(`n_replic',1,.);


matrix correlation=J(`n_replic',1,.);



matrix omega=J(1,1,.);




gen double u1=1;
gen double y1=1;
gen double y2=1;

gen double y1_star=1;

gen double ln_resid_sq=1;

foreach p_rho of numlist   2 1 -2 -1 0 {;
foreach lambda of numlist  0 0.5  { ;
set seed 10000;

forvalues i=1/`n_replic' {;

display "replication=`i'";



matrix sigma=(1.0, 0.5,0.5, 0, 0\0.5,1.0, 0.5, 0, 0\0.5,0.5,1.0,0,0\0,0,0,1.0,0\0,0,0,0,1.0);
drawnorm x1 z1  z2 v1 eta1,  n(`n_obs') cov(sigma) double;




quietly replace eta1=exp(`lambda'*z1)*eta1;


quietly replace u1=v1+`p_rho'*eta1;



quietly replace y2=`p_gamma0'+`p_gamma1'*x1+`p_gamma2'*z1+`p_gamma3'*z2+eta1;
quietly  replace  y1_star=`p_beta0'+`p_alpha'*y2+`p_beta1'*x1+u1;
quietly replace y1=1 if y1_star>0;
quietly replace y1=0 if y1_star<=0;

gen orig_index=_n;
gen const=1;
sort orig_index;
quietly reg y2  x1 z1 z2 ;
predict double resid1,r;
gen double sigma2_sq=resid1^2;
gen expand_index=_N;

expand expand_index;
sort orig_index;
bysort orig_index: gen new_index=_n;
generate double new_x1=x1[expand_index*new_index];
generate double new_z1=z1[expand_index*new_index];
generate double new_z2=z2[expand_index*new_index];
generate double new_sigma2_sq=sigma2_sq[expand_index*new_index];
display "finish expansion";
sort orig_index new_index;
order orig_index new_index x1 new_x1 z1 new_z1 z2 new_z2;
 gen double dist=sqrt((x1-new_x1)^2+(z1-new_z1)^2    +(z2-new_z2)^2);

quietly bysort orig_index (dist) :gen dist_rank=sum(dist  != dist[_n-1]);


sort orig_index dist_rank;

quietly gen double cv=. ;
quietly gen double resid_sq=. ;
quietly gen k_index=.;



/* local near_k 3*/;

/* from here, it is CV; 

display "finding optimal k";
forvalue near_k=2/40 {;
quietly by orig_index: egen avg_sigma2_sq=mean(new_sigma2_sq) if dist_rank>=2 & dist_rank<=`near_k' ;




quietly replace  resid_sq=(sigma2_sq-avg_sigma2_sq)^2  if dist_rank==2;
/*we calculate the residual for each original point */

quietly egen double total_resid=sum(resid_sq);
/* in the above command, all rows becomes the same value for total resid */

quietly replace cv=total_resid if _n==`near_k';
quietly replace k_index=`near_k' if _n==`near_k';
/* put the CV value at the k-th row */


sort orig_index new_index;
 
drop avg_sigma2 total_resid ;
};

sort cv;

local opt_k=k_index[1];
matrix opt_khat[`i', 1]=`opt_k' ;

*/
drop cv  k_index;



bysort orig_index: egen avg_sigma2_sq=mean(new_sigma2_sq) if dist_rank>=2 & dist_rank<=`opt_k' ;
bysort orig_index (dist_rank): replace avg_sigma2_sq=avg_sigma2_sq[2];
keep if dist_rank==1;

drop new_x1 new_z1  new_z2 new_sigma2_sq new_index;

generate sigma2=sqrt(avg_sigma2_sq);

generate double  x1_hat=x1/avg_sigma2_sq;
generate double z1_hat=z1/avg_sigma2_sq;
generate double z2_hat=z2/avg_sigma2_sq;
 
 










 
 * quietly reg y2 x1 z1 z2;
 * predict double resid1, r;
 * matrix est0=e(b);
 quietly probit y1 y2 x1 resid1;
 matrix  est1=e(b);
generate double index1=est1[1,1]*y2+est1[1,2]*x1+est1[1,3]*resid1+est1[1,4];


/*calculate normal density based on acculate resid2 */
gen double normden_hat=normalden(index1);
replace normden_hat=1e-15 if normden_hat==0;
gen double normcdf_hat=normal(index1);
 replace normcdf_hat=1-1e-15 if normcdf_hat==1;
 replace normcdf_hat=1e-15 if normcdf_hat==0;

 

 
 /*generate weight for optimal instrument for equation 1*/
generate double weight= normden_hat/ (normcdf_hat*(1-normcdf_hat) );
 
 

quietly gen double ins_alpha=-weight*y2;
quietly gen double ins_beta0=-weight*const;
quietly gen double ins_beta1=-weight*x1;
quietly gen double ins_rho=-weight*resid1;
quietly gen double ins_resid1=weight;
quietly gen double ins_resid2=-1/avg_sigma2_sq;
quietly gen double ins_gamma0=const;
quietly gen double ins_gamma1=x1;
quietly gen double ins_gamma2=z1;
quietly gen double ins_gamma3=z2;
 capture noisily  gmm (eq1:y1-normal({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
  (eq2:  ins_resid1* ( y1-normal({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 + ins_resid2*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) )  ,
	 instruments(eq1:ins_beta0 ins_alpha ins_beta1 ins_rho ,noconstant )
	 instruments(eq2:  ins_gamma0 ins_gamma1 ins_gamma2 ins_gamma3, noconstant)
	 derivative(eq1/ beta0= -normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1 -{gamma3}*z2  ) ) )
	 derivative(eq1/alpha=-y2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2   ) ))
	 derivative(eq1/beta1=-x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)   ))
	 derivative(eq1/rho=-(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)*normalden({beta0}+{alpha}*y2
	 +{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)))
	 derivative(eq1/gamma0={rho}*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)) )
	 derivative(eq1/gamma1={rho}*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 derivative(eq1/gamma2={rho}*z1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 derivative(eq1/gamma3={rho}*z2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 derivative(eq2/beta0=-ins_resid1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/alpha=-ins_resid1*y2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/beta1=-ins_resid1*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/rho=-ins_resid1*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)*normalden({beta0}+{alpha}*y2+{beta1}*x1
	 +{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/gamma0=ins_resid1*{rho}*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2)
	 derivative(eq2/gamma1=ins_resid1*{rho}*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2*x1)
	 derivative(eq2/gamma2=ins_resid1*{rho}*z1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2*z1)
	  derivative(eq2/gamma3=ins_resid1*{rho}*z2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2*z2)
	 conv_vtol(1e-14) conv_ptol(1e-14) conv_nrtol(1e-14) conv_maxiter(30) winitial(identity)  onestep technique(gn ) 
 /* from(alpha 1 beta1 1 rho 1 beta0 1 gamma0 1 gamma1 0 gamma2 0 gamma3 -1) */
	 tracelevel(`tolerance')   ;

scalar conv=e(converged);
scalar Q_value=e(Q);

 if conv==1 & Q_value<1e-10 { ;
 matrix est2=e(b);
 matrix var2=e(V);
 matrix beta0hat_gmm1[`i', 1]=est2[1,1] ;
 matrix alphahat_gmm1[`i', 1]=est2[1,2] ;
 matrix beta1hat_gmm1[`i',1]=est2[1,3];
 matrix rhohat_gmm1[`i',1]=est2[1,4];
 matrix gamma0hat_gmm1[`i',1]=est2[1,5];
  matrix gamma1hat_gmm1[`i',1]=est2[1,6];
 matrix gamma2hat_gmm1[`i',1]=est2[1,7];
 
 matrix drop est2;
 
 matrix beta0_stder_gmm1[`i',1]=sqrt(var2[1,1]);
 matrix alpha_stder_gmm1[`i',1]=sqrt(var2[2,2]);
 matrix beta1_stder_gmm1[`i',1]=sqrt(var2[3,3]);
 matrix rho_stder_gmm1[`i',1]=sqrt(var2[4,4]);
 matrix gamma0_stder_gmm1[`i',1]=sqrt(var2[5,5]);
 matrix gamma1_stder_gmm1[`i',1]=sqrt(var2[6,6]);
 matrix gamma1_stder_gmm1[`i',1]=sqrt(var2[7,7]);
 matrix drop var2;
 
quietly probit y1 y2 x1 resid1;
matrix  est3=e(b);

 matrix alphahat_rv[`i', 1]=est3[1,1] ;
 matrix beta1hat_rv[`i',1]=est3[1,2];
 matrix rhohat_rv[`i',1]=est3[1,3];
 matrix beta0hat_rv[`i', 1]=est3[1,4] ;
 matrix drop est3;
 
quietly ivreg y2 (x1 z1 z2 =x1_hat z1_hat z2_hat);
predict double resid2_hat,r;
 quietly probit y1 y2 x1 resid2_hat;
 matrix est4=e(b);
  matrix alphahat_gmm2[`i', 1]=est4[1,1] ;
 matrix beta1hat_gmm2[`i',1]=est4[1,2];
 matrix rhohat_gmm2[`i',1]=est4[1,3];
 matrix beta0hat_gmm2[`i', 1]=est4[1,4] ;
 drop resid2_hat;
 

};


if conv==0 | Q_value>=1e-10 { ;
 	 
	 capture noisily  gmm (eq1:y1-normal({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
  (eq2:  ins_resid1* ( y1-normal({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 + ins_resid2*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) )  ,
	 instruments(eq1:ins_alpha ins_beta1 ins_rho  )
	 instruments(eq2:  ins_gamma1 ins_gamma2 ins_gamma3 )
	 derivative(eq1/ beta0= -normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1 -{gamma3}*z2  ) ) )
	 derivative(eq1/alpha=-y2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2   ) ))
	 derivative(eq1/beta1=-x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)   ))
	 derivative(eq1/rho=-(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)*normalden({beta0}+{alpha}*y2
	 +{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)))
	 derivative(eq1/gamma0={rho}*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)) )
	 derivative(eq1/gamma1={rho}*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 derivative(eq1/gamma2={rho}*z1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 derivative(eq1/gamma3={rho}*z2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2))  )
	 derivative(eq2/beta0=-ins_resid1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/alpha=-ins_resid1*y2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/beta1=-ins_resid1*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/rho=-ins_resid1*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2)*normalden({beta0}+{alpha}*y2+{beta1}*x1
	 +{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) )
	 derivative(eq2/gamma0=ins_resid1*{rho}*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2)
	 derivative(eq2/gamma1=ins_resid1*{rho}*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2*x1)
	 derivative(eq2/gamma2=ins_resid1*{rho}*z1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2*z1)
	  derivative(eq2/gamma3=ins_resid1*{rho}*z2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*z1-{gamma3}*z2) ) 
	 -ins_resid2*z2)
	 
	 
conv_vtol(1e-14) conv_ptol(1e-14) conv_nrtol(1e-14) conv_maxiter(600) winitial(identity) onestep technique(nr 30 dfp 30 bfgs 30 )  
	 tracelevel(`tolerance') 
	/*  from(alpha -1 beta1 -1 rho 1 beta0 1 gamma0 1 gamma1 0 gamma2 0 gamma3 -1)  */;
 
 scalar conv=e(converged);
scalar Q_value=e(Q);


 if conv==1 & Q_value<1e-10 { ;
 matrix est2=e(b);
 matrix var2=e(V);
 matrix beta0hat_gmm1[`i', 1]=est2[1,1] ;
 matrix alphahat_gmm1[`i', 1]=est2[1,2] ;
 matrix beta1hat_gmm1[`i',1]=est2[1,3];
 matrix rhohat_gmm1[`i',1]=est2[1,4];
 matrix gamma0hat_gmm1[`i',1]=est2[1,5];
  matrix gamma1hat_gmm1[`i',1]=est2[1,6];
 matrix gamma2hat_gmm1[`i',1]=est2[1,7];
 
 matrix drop est2;
 quietly probit y1 y2 x1 resid1;
matrix  est3=e(b);

 matrix alphahat_rv[`i', 1]=est3[1,1] ;
 matrix beta1hat_rv[`i',1]=est3[1,2];
 matrix rhohat_rv[`i',1]=est3[1,3];
 matrix beta0hat_rv[`i', 1]=est3[1,4] ;
 matrix drop est3;
 
 ivreg y2 (x1 z1 z2 =x1_hat z1_hat z2_hat);
predict double resid2_hat,r;
 probit y1 y2 x1 resid2_hat;
 matrix est4=e(b);
  matrix alphahat_gmm2[`i', 1]=est4[1,1] ;
 matrix beta1hat_gmm2[`i',1]=est4[1,2];
 matrix rhohat_gmm2[`i',1]=est4[1,3];
 matrix beta0hat_gmm2[`i', 1]=est4[1,4] ;
 
 drop resid2_hat;
 
 matrix beta0_stder_gmm1[`i',1]=sqrt(var2[1,1]);
 matrix alpha_stder_gmm1[`i',1]=sqrt(var2[2,2]);
 matrix beta1_stder_gmm1[`i',1]=sqrt(var2[3,3]);
 matrix rho_stder_gmm1[`i',1]=sqrt(var2[4,4]);
 matrix gamma0_stder_gmm1[`i',1]=sqrt(var2[5,5]);
 matrix gamma1_stder_gmm1[`i',1]=sqrt(var2[6,6]);
 matrix gamma1_stder_gmm1[`i',1]=sqrt(var2[7,7]);
 
 matrix drop var2;
 
 
 
 /* from here, I calclucate the SE by equation */
 
 

/*this is the end of "conv==1 & Q_value<1e-10 { ; */
};
 

 /* this is the end of "if conv==0 | Q_value>=1e-10 { ;  "   */
 };







 
 
 /*
 gen sder_beta0=1/sigma1;
 gen sder_alpha=normden_hat*y2/sigma1;
 gen sder_beta1=normden_hat*x2/sigma1;
 gen sder_rho=sigma1*resid2/sigma1;
  gen sder_gamma1=sigma2*x2/sigma2_sq;
 gen sder_gamma2=sigma2*z1/sigma2_sq;
 gen sder_gamma3=sigma2*z2/sigma2_sq;
 gen sder_gamma0=sigma2/sigma2_sq;
 
 matrix accum A=ins_alpha ins_beta1 ins_rho ins_beta0, noconstant;
 quietly matrix list A;
 matrix accum B=ins_gamma1 ins_gamma2 ins_gamma3 ins_gamma0, noconstant;
 matrix zero1=J(4,4,0);
 matrix zero2=J(4,4,0);
 
 matrix D=(A,zero2\zero1, B);
 
 matrix stder_matrix=invsym(D);
di "so far it is ok";  
 
 */
 
 
 
 
  
 
 
 
 
 
 drop x1 z1 z2 v1 eta1 sigma2_sq resid1 index1 normden_hat  normcdf_hat weight ins_alpha ins_beta0 ins_beta1 
 ins_rho ins_resid1 ins_resid2 ins_gamma0 ins_gamma1 ins_gamma2 ins_gamma3 sigma2 x1_hat z1_hat z2_hat 
 orig_index const expand_index dist dist_rank resid_sq avg_sigma2_sq;
 
 ;






/* this is the end of forvalue loop */
};




timer off 1;


 
 svmat alphahat_gmm1;
svmat beta1hat_gmm1;
 svmat beta0hat_gmm1;
svmat rhohat_gmm1;

svmat alphahat_rv;
svmat beta1hat_rv;
 svmat beta0hat_rv;
svmat rhohat_rv;

svmat alphahat_gmm2;
svmat beta1hat_gmm2;
 svmat beta0hat_gmm2;
svmat rhohat_gmm2;

 svmat alpha_stder_gmm1;
svmat beta1_stder_gmm1;
 svmat beta0_stder_gmm1;
svmat rho_stder_gmm1;
 

 
 
 

summarize alphahat_gmm1;
summarize alpha_stder_gmm1;
sum alphahat_rv;
sum alphahat_gmm2;

save result_knn_rho`p_rho'lambda`lambda'.dta, replace;
drop alphahat_gmm1 beta1hat_gmm1 beta0hat_gmm1 rhohat_gmm1



alphahat_gmm2 beta1hat_gmm2 beta0hat_gmm2 rhohat_gmm2 alphahat_rv beta1hat_rv
beta1hat_rv beta0hat_rv rhohat_rv alpha_stder_gmm1 beta1_stder_gmm1 beta0_stder_gmm1 rho_stder_gmm1;
 



;

 };
 };

;


log close;
