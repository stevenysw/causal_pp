# score-based estimation of heterogeneous treatment effects for high-dimensional data
# Inputs:
## X: covariates
## y: outcome variable
## z: treatment indicator
## K: number of nearest neighbors
# Outputs:
## tau: estimation of treatment effects
## beta_e: coefficient vector for estimating propensity scores
## beta_p: coefficient vector for estimating prognostic scores
## model: prediction model

require(glmnet)
require(MatchIt)

causal_pp_high <- function(X, y, z, K){
  n = length(y)
  
  control = which(z == 0)
  treatment = which(z == 1)
  
  ### estimated propensity score
  cv.lasso <- cv.glmnet(X, z, alpha = 1, family = "binomial")
  lam = cv.lasso$lambda.min
  temp = glmnet(X, z, alpha = 1, family = "binomial", lambda = cv.lasso$lambda.min)
  beta_hat_e = coef(temp)
  logit = cbind(rep(1, n), X) %*% beta_hat_e
  e_hat = as.vector(1 / (1 + exp(-logit)))
  
  ## estimated prognostic score
  X_con = X[control,]
  y_con = y[control]
  cv_model <- cv.glmnet(X_con, y_con, alpha = 1)
  best_lambda <- cv_model$lambda.min
  temp = glmnet(X_con, y_con, alpha = 1, lambda = best_lambda)
  beta_hat_p = coef(temp)
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
