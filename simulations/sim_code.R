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
  tau_cart = causal_pp(X,y,z, round(log(n)))
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
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 1, d = 50, n = 1000")

# simulation 2
n = 5000 ## sample size
d = 10 ## ambient dimension
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
  temp = drop(X %*% beta_e) + 
    0.5 * beta_e[1] * X[,1] * X[,3] +
    0.7 * beta_e[2] * X[,2] * X[,4] +
    0.5 * beta_e[3] * X[,3] * X[,5] +
    0.7 * beta_e[4] * X[,4] * X[,6] +
    0.5 * beta_e[5] * X[,5] * X[,7] +
    0.5 * beta_e[1] * X[,1] * X[,6] +
    0.7 * beta_e[2] * X[,2] * X[,3] +
    0.5 * beta_e[3] * X[,3] * X[,4] +
    0.5 * beta_e[4] * X[,4] * X[,5] +
    0.5 * beta_e[5] * X[,5] * X[,6] +
    beta_e[2] * X[,2] * X[,2] +
    beta_e[4] * X[,4] * X[,4] +
    beta_e[10] * X[,10] * X[,10]
   
  ## propensity score
  e = exp(temp) / (1 + exp(temp))
  z = rbinom(rep(1,n), rep(1,n), prob = e)
 
  ## prognostic score
  p = drop(X %*% beta_p) +
    0.5 * beta_p[1] * X[,1] * X[,3] +
    0.7 * beta_p[2] * X[,2] * X[,4] +
    0.5 * beta_p[3] * X[,3] * X[,8] +
    0.7 * beta_p[4] * X[,4] * X[,9] +
    0.5 * beta_p[5] * X[,8] * X[,10] +
    0.5 * beta_p[1] * X[,1] * X[,9] +
    0.7 * beta_p[2] * X[,2] * X[,3] +
    0.5 * beta_p[3] * X[,3] * X[,4] +
    0.5 * beta_p[4] * X[,4] * X[,8] +
    0.5 * beta_p[5] * X[,8] * X[,9] +
    beta_p[2] * X[,2] * X[,2] +
    beta_p[4] * X[,4] * X[,4] +
    beta_p[10] * X[,10] * X[,10]

  ## response
  tau0 = as.numeric(e < 0.6 & p < 0)
  
  y = p + tau0 * z + rnorm(n, mean = 0, sd = 1)

  ## causal pp estimation
  tau_cart = causal_pp(X,y,z, round(log(n)))
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
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 2, d = 10, n = 5000")

# simulation 3
n = 5000 ## sample size
d = 10 ## ambient dimension
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
  e1 = 0.6
  p1 = 0

  tau0 = rep(0, n)
  for (i in 1:n){
    if (e[i] <= e1 & p[i] <= p1) tau0[i] = 0
    if (e[i] > e1 & p[i] <= p1) tau0[i] = 1
    if (e[i] <= e1 & p[i] > p1) tau0[i] = 1
    if (e[i] > e1 & p[i] > p1) tau0[i] = 2
  }
  
  y = p + tau0 * z + rnorm(n, mean = 0, sd = 1)

  ## causal pp estimation
  tau_cart = causal_pp(X,y,z, round(log(n)))
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
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 3, d = 10, n = 5000")

# simulation 4


# simulation 5


# simulation 6

# simulation 7


