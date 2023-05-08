# endogenous stratification - full sample
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator

# Outputs:
## tau: estimation of treatment effects
## tau_itv: treatment effects of each subgroup
## beta_p: coefficient vector for estimating prognostic scores
## interval: cut-offs on prognostic scores for subgroup

strat_full <- function(X, y, z){
  control = which(z == 0)
  treatment = which(z == 1)
  
  temp = lm(y[control] ~ X[control,])
  beta_hat_p = temp$coefficients
  p_hat = drop(cbind(rep(1, length(y)), X) %*% beta_hat_p)
  
  itv = c(min(p_hat)-1, quantile(p_hat, c(1/3, 2/3)), max(p_hat) + 1)
  
  group = as.numeric(cut(p_hat, itv))
  
  tau_group = rep(0, 3)
  for (i in 1:3){
    tau_group[i] = sum(y[group == i & z == 1]) / sum(as.numeric(group == i & z == 1)) - sum(y[group == i & z == 0]) / sum(as.numeric(group == i & z == 0))
  }
  
  tau = tau_group[group]
  return(list(tau = tau_group[group], tau_itv = tau_group, beta_p = beta_hat_p, interval = itv))
}
