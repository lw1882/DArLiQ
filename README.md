# DArLiQ

Dynamic Autoregressive Liquidity (DArLiQ)

Example: Apple and Bitcoin

func_utils.r: functions e.g. numerical gradient used in other functions.

func_step1_g\_bandwidth.r: bandwidth selection using derived rule of thumb in DArLiq based on the assumption that g(u)=exp(a0+a1\*u). Higher-order polynomials can also be used.

func_step1_g\_trend.r: function to estimate the long-run trend function g via local linear estimator.

func_step2_lambda_update.r: updating equation of the lambda process.

func_step2_lambda_GMM.r: functions related to the estimation of the dynamic parameters of the lambda process via GMM.

func_step2_lambda_MLE.r: functions related to the estimation of the dynamic parameters of the lambda process via MLE.

main_data.r: construct the illiquidity series.

main_estimate.r: example to estimate the g function and the dynamic parameters via GMM/MLE.
