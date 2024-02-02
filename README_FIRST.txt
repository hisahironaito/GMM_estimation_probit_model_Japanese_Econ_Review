To replicate the simulation result, first please run "JER_replication_simulation.do". It will take 4-6 hours depending on CPU power. Then, please run "JER_making_table_from_simulation.do".
 Depending on the number of core and stata version, the result does not become exactly the same. But The result should be almost the same.
For example, please look at the caes where rho =2 and lambda=1. Our paper reported that in this case, the efficiency gain of using GMM probit is much larger in terms of RMSE.
The replication from the simulation should show such a result.   RMSE of RV method should be 0.3270477 whiele RMSE of using GMM probit should be 0.2085812


To replicate the empirical part, please run "regmjlee_JER_replication.do". The key paramter of the interest is alpha, the coefficient of the continous endogenous variable.
 The value of alpha should be -.00013 and SE should be SE= (.0015166) 
