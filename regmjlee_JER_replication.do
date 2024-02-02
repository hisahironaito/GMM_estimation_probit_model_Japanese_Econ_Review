* This do file calculates the opitmal k of the k-nn in the empirical model


#delimit;
clear all;
set memory 6000m;


set more off;
use mjlee.dta;
local n_replic=1;
des;
sum;
matrix alpha_stder_eqn=J(`n_replic',1,.);
matrix beta1_stder_eqn=J(`n_replic',1,.);
matrix beta2_stder_eqn=J(`n_replic',1,.);
local opt_k 44;

generate x0=1;
 rename emp y1;
 rename ohhinc y2;
 rename age x1;
 rename edu x2;
 rename children_0_5 x3;
 rename nonwhite x4;
 rename hus_manager z1;
 rename hus_sales  z2;
 rename hus_farm z3;
 
 
 
gen orig_index=_n;
gen const=1;
sort orig_index;
quietly reg y2  x1  x2 x3 x4 z1 z2 z3 ;
predict double resid1,r;
gen double sigma2_sq=resid1^2;
gen expand_index=_N;



foreach j in x1 x2 x3 x4 z1 z2 z3{;
summarize `j';
scalar var_`j'=r(Var);
};


expand expand_index;
sort orig_index;
bysort orig_index: gen new_index=_n;
generate double new_x1=x1[expand_index*new_index];
generate double new_x2=x2[expand_index*new_index];
generate double new_x3=x3[expand_index*new_index];
generate double new_x4=x4[expand_index*new_index];
generate double new_z1=z1[expand_index*new_index];
generate double new_z2=z2[expand_index*new_index];
generate double new_z3=z3[expand_index*new_index];
generate double new_sigma2_sq=sigma2_sq[expand_index*new_index];
display "finish expansion";
sort orig_index new_index;
order orig_index new_index x1 new_x1    x2 new_x2 x3 new_x3 x4 new_x4 z1 new_z1 z2 new_z2 z3 new_z3 ;

gen double dist=sqrt(((x1-new_x1)^2)/var_x1+  ((x2-new_x2)^2)/var_x2+((x3-new_x3)^2)/var_x3
 +((x4-new_x4)^2)/var_x4+ ((z1-new_z1)^2)/var_z1+((z2-new_z2)^2)/var_z2+((z3-new_z3)^2)/var_z3 );



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
display "finished `near_k' iteration";
};

sort cv;

local opt_k=k_index[1];
generate opt_k=`opt_k';
*matrix opt_khat[`i', 1]=`opt_k' ;

drop cv  k_index;

*/;

quietly by orig_index: egen avg_sigma2_sq_temp=mean(new_sigma2_sq) if dist_rank>=2 & dist_rank<=`opt_k' & dist!=0;
quietly by orig_index: egen avg_sigma2_sq=mean(avg_sigma2_sq_temp);





keep if orig_index==new_index;


drop new_x1 new_x2 new_x3 new_x4 new_z1 new_z2 new_z3    new_sigma2_sq new_index;


/* from here, I try hybrid model */;

gen double x1_hat=x1/avg_sigma2_sq;
gen double x2_hat=x2/avg_sigma2_sq;
gen double x3_hat=x3/avg_sigma2_sq;
gen double x4_hat=x4/avg_sigma2_sq;
gen double z1_hat=z1/avg_sigma2_sq;
gen double z2_hat=z2/avg_sigma2_sq;
gen double z3_hat=z3/avg_sigma2_sq;


ivregress 2sls y2  (x1 x2 x3 x4 z1 z2 z3=x1_hat x2_hat x3_hat x4_hat z1_hat z2_hat z3_hat);
predict resid2,r;

probit y1 y2 x1 x2 x3 x4 resid2;




probit y1 y2 x1 x2 x3 x4 resid1;
matrix  est1=e(b);
generate double index1=est1[1,1]*y2+est1[1,2]*x1+est1[1,3]*x2 + est1[1,4]*x3+est1[1,5]*x4+est1[1,6]*resid1+est1[1,7];


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
quietly gen double ins_beta2=-weight*x2;
quietly gen double ins_beta3=-weight*x3;
quietly gen double ins_beta4=-weight*x4;
quietly gen double ins_rho=-weight*resid1;
quietly gen double ins_resid1=weight;
 gen double ins_resid2=-1/avg_sigma2_sq;
quietly gen double ins_gamma0=const;
quietly gen double ins_gamma1=x1;
quietly gen double ins_gamma2=x2;
quietly gen double ins_gamma3=x3;
quietly gen double ins_gamma4=x4;
quietly gen double ins_gamma5=z1;
quietly gen double ins_gamma6=z2;
quietly gen double ins_gamma7=z3;

 gmm (eq1:y1-normal({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4
+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )
  (eq2:  ins_resid1* ( y1-normal({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4
  +{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )
	 + ins_resid2*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3 ) )  ,
	 instruments(eq1:ins_beta0 ins_alpha ins_beta1  ins_beta2 ins_beta3 ins_beta4 ins_rho ,noconstant )
	 instruments(eq2:  ins_gamma0 ins_gamma1 ins_gamma2 ins_gamma3 ins_gamma4 ins_gamma5 ins_gamma6 ins_gamma7, noconstant)
	 
	 
	 
	 derivative(eq1/ beta0= -normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3  ) ) )

	 derivative(eq1/alpha=-y2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4
	 +{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3   ) ))

	 derivative(eq1/beta1=-x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)   ))

 derivative(eq1/beta2=-x2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)   ))

 derivative(eq1/beta3=-x3*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)   ))

 derivative(eq1/beta4=-x4*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)   ))


	 derivative(eq1/rho=-(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)*normalden({beta0}+{alpha}*y2
	 +{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)))

	 derivative(eq1/gamma0={rho}*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)) )

	 derivative(eq1/gamma1={rho}*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )

	 derivative(eq1/gamma2={rho}*x2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )


derivative(eq1/gamma3={rho}*x3*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )

derivative(eq1/gamma4={rho}*x4*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )


	 derivative(eq1/gamma5={rho}*z1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )

 derivative(eq1/gamma6={rho}*z2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )

 derivative(eq1/gamma7={rho}*z3*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3))  )



	 derivative(eq2/beta0=-ins_resid1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )
	 derivative(eq2/alpha=-ins_resid1*y2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )
	 derivative(eq2/beta1=-ins_resid1*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )
	 derivative(eq2/beta2=-ins_resid1*x2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )
	 derivative(eq2/beta3=-ins_resid1*x3*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )
	 derivative(eq2/beta4=-ins_resid1*x4*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )

derivative(eq2/rho=-ins_resid1*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3)*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4
	 +{rho}*(y2-{gamma0}-{gamma1}*x1--{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) )
	 
	 derivative(eq2/gamma0=ins_resid1*{rho}*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2)
	 
	 derivative(eq2/gamma1=ins_resid1*{rho}*x1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*x1)
	 
	 derivative(eq2/gamma2=ins_resid1*{rho}*x2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*x2)
	 
	 derivative(eq2/gamma3=ins_resid1*{rho}*x3*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*x3)
	 
	 derivative(eq2/gamma4=ins_resid1*{rho}*x4*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*x4)
	 
	 derivative(eq2/gamma5=ins_resid1*{rho}*z1*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*z1)
	 
	 derivative(eq2/gamma6=ins_resid1*{rho}*z2*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*z2)
	 
	 derivative(eq2/gamma7=ins_resid1*{rho}*z3*normalden({beta0}+{alpha}*y2+{beta1}*x1+{beta2}*x2+{beta3}*x3+{beta4}*x4+{rho}*(y2-{gamma0}-{gamma1}*x1-{gamma2}*x2-{gamma3}*x3-{gamma4}*x4
-{gamma5}*z1-{gamma6}*z2-{gamma7}*z3) ) 
	 -ins_resid2*z3)
	 
	 
	 from(alpha 0.002 beta0 0.996 beta1 -0.038  beta2 0.108 beta3 -0.510 beta4 0.209 rho -0.009 
  gamma0 -28.696 gamma1 0.702  gamma2 0 gamma2 2.254 gamma3  1.330 gamma4 -4.550 gamma5 12.085 gamma6 5.082 gamma7 -2.552)
	 conv_vtol(1e-14) conv_ptol(1e-14)  conv_maxiter(300) winitial(identity)  onestep technique(gn)  tracelevel(`tolerance')   ;
  
	

/* alpha is our interest. which is a continous edogenous regressor. The value of alpha and SE should be -.00013  SE= (.0015166) */

