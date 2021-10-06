# endogenous stratification - leave-one-out
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator
## K: number of subgroups

# Outputs:
## tau: estimation of treatment effects

strat_loo <- function(X, y, z, K){
  n = length(y)
  control = which(z == 0)
  treatment = which(z == 1)
  
  p_hat = rep(0, n)
  
  temp = lm(y[control] ~ X[control,])
  beta_hat_p = as.vector(temp$coefficients)
  p_hat[treatment] = drop(cbind(rep(1, length(treatment)), X[treatment,]) %*% beta_hat_p)
  
  for (j in 1:length(control)){
    temp = lm(y[control][-j] ~ X[control,][-j,])
    coef = temp$coefficients
    p_hat[control[j]] = drop(c(1, X[control,][j,]) %*% coef)
  }
  
  itv = c(min(p_hat)-1, quantile(p_hat, (1:(K-1))/K), max(p_hat) + 1)
  
  group = as.numeric(cut(p_hat, itv))
  
  tau_group = rep(0, K)
  for (i in 1:K){
    tau_group[i] = sum(y[group == i & z == 1]) / sum(as.numeric(group == i & z == 1)) - sum(y[group == i & z == 0]) / sum(as.numeric(group == i & z == 0))
  }
  
  tau = tau_group[group]
  return(tau)
}
