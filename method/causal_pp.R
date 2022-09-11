# score-based estimation of heterogeneous treatment effects
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator
# K: number of nearest neighbors
# Outputs:
## tau: estimation of treatment effects
## beta_e: coefficient vector for estimating propensity scores
## beta_p: coefficient vector for estimating prognostic scores
## model: prediction model

require(MatchIt)

causal_pp <- function(X, y, z, K){
  n = length(y)
  
  control = which(z == 0)
  treatment = which(z == 1)
  
  ### estimated propensity score
  temp = glm(z ~ X, family = "binomial")
  beta_hat_e = temp$coefficients
  e_hat = temp$fitted.values
  
  ## estimated prognostic score
  temp = lm(y[control] ~ X[control,])
  beta_hat_p = temp$coefficients
  p_hat = drop(cbind(rep(1, n), X) %*% beta_hat_p)
  
  ## knn mahalanobis distance match
  y_match = rep(0, n)

  match_temp_1 = matchit(z ~ e_hat + p_hat, method = "nearest", ratio = K, distance = "mahalanobis", replace = T)
  match_temp_2 = matchit(1 - z ~ e_hat + p_hat, method = "nearest", ratio = K, distance = "mahalanobis", replace = T)
  idx_mat = rbind(match_temp_1$match.matrix, match_temp_2$match.matrix)
  idx_mat = idx_mat[order(as.numeric(rownames(idx_mat))),]
  
  for(i in 1:n){
    y_match[i] = mean(y[as.numeric(idx_mat[i,])])
  }
   
  ytilde = (y - y_match)*(2*z-1)

  cart_mod = rpart(ytilde ~ e_hat + p_hat)
  tau_final = predict(cart_mod)
  
  return(list(tau = tau_final, beta_e = beta_hat_e, beta_p = beta_hat_p, model = cart_mod))
}
