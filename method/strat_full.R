# endogenous stratification - full sample
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator
## K: number of subgroups

# Outputs:
## tau: estimation of treatment effects


strat_full <- function(X, y, z, K){
  control = which(z == 0)
  treatment = which(z == 1)
  
  temp = lm(y[control] ~ X[control,])
  beta_hat_p = temp$coefficients
  p_hat = drop(cbind(rep(1, length(y)), X) %*% beta_hat_p)
  
  itv = c(min(p_hat)-1, quantile(p_hat, (1:(K-1))/K), max(p_hat) + 1)
  
  group = as.numeric(cut(p_hat, itv))
  
  tau_group = rep(0, K)
  for (i in 1:K){
    tau_group[i] = sum(y[group == i & z == 1]) / sum(as.numeric(group == i & z == 1)) - sum(y[group == i & z == 0]) / sum(as.numeric(group == i & z == 0))
  }
  
  tau = tau_group[group]
  return(tau)
}
