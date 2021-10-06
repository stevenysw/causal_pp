# prognostic score matching
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator
# Outputs:
## tau: estimation of treatment effects

prognostic_match <- function(X, y, z){
  n = length(y)
  
  control = which(z == 0)
  treatment = which(z == 1)
  
  ## estimated prognostic score
  temp = lm(y[control] ~ X[control,] - 1)
  beta_hat_p = temp$coefficients
  p_hat = drop(X %*% beta_hat_p)
  
  y_match = rep(0, n)
  
  for(i in 1:n){
    if (z[i] == 1){
      temp = get.knnx(p_hat[control], p_hat[i], k = 5)
      match_idx = temp$nn.index
      y_match[i] = mean(y[control][match_idx])
    }
    if (z[i] == 0){
      temp = get.knnx(p_hat[treatment], p_hat[i], k = 5)
      match_idx = temp$nn.index
      y_match[i] = mean(y[treatment][match_idx])
    }
  }
  
  tau = (y - y_match)*(2*z-1)
  return(tau)
}