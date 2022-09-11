require(grf)

# simulation 1
n = 1000 ## sample size
d = 50 ## ambient dimension
N = 1000 ## number of repetitions

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  beta_e = runif(d, -1, 1)
  beta_p = runif(d, -1, 1)

  X = matrix(runif(n*d), n, d)
  
  ## treatment indicators
  temp = drop(X %*% beta_e)
   
  ## propensity score
  e = exp(temp) / (1 + exp(temp))
  z = rbinom(rep(1,n), rep(1,n), prob = e)
 
  ## prognostic score
  p = drop(X %*% beta_p)

  ## response
  tau0 = as.numeric(e < 0.6 & p < 0)
  
  y = p + tau0 * z + rnorm(n, mean = 0, sd = 1)

  ## causal pp estimation
  tau_cart = causal pp(X,y,z, round(log(n)))
  tau_ca = tau_cart$tau
  mse_ca[j] = mean((tau_ca - tau0)^2)
  
  ## endogenous stratification
  tau_st = strat_loo(X, y, z)
  mse_st[j] = mean((tau_st - tau0)^2)
  
  ## causal forests
  tau_forest = causal_forest(X, y, z, tune.parameters = "all")
  pred = predict(tau_forest)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0)^2)
}

## plot of comparison
boxplot(mse_ca_mh, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 1, d = 50, n = 1000")
