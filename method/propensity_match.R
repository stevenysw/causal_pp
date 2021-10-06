# propensity score matching (PSM)
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator
# Outputs:
## tau: estimation of treatment effects

propensity_match <- function(X, y, z){
  n = length(y)
  
  control = which(z == 0)
  treatment = which(z == 1)
  
  ### estimated propensity score
  temp = glm(z ~ X - 1, family = "binomial")
  e_hat = temp$fitted.values
  
  y_match = rep(0, n)
  
  for(i in 1:n){
    if (z[i] == 1){
      temp = get.knnx(e_hat[control], e_hat[i], k = 5)
      match_idx = temp$nn.index
      y_match[i] = mean(y[control][match_idx])
    }
    if (z[i] == 0){
      temp = get.knnx(e_hat[treatment], e_hat[i], k = 5)
      match_idx = temp$nn.index
      y_match[i] = mean(y[treatment][match_idx])
    }
  }
  
  tau = (y - y_match)*(2*z-1)
  return(tau)
}