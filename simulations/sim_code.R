require(grf)

# simulation 1
n_train = 1000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 50 ## ambient dimension
N = 1000 ## number of repetitions

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  beta_e = runif(d, -1, 1)
  beta_p = runif(d, -1, 1)
  
  X_train = matrix(runif(n_train*d), n_train, d)
  X_test = matrix(runif(n_test*d), n_test, d)
  
  ## propensity score
  temp_train = drop(X_train %*% beta_e)
  temp_test = drop(X_test %*% beta_e)
  e_train = exp(temp_train) / (1 + exp(temp_train))
  z_train = rbinom(rep(1,n_train), rep(1,n_train), prob = e_train)
  e_test = exp(temp_test) / (1 + exp(temp_test))
  
  ## prognostic score
  p_train = drop(X_train %*% beta_p)
  p_test = drop(X_test %*% beta_p)
  
  ## response
  tau0_train = as.numeric(e_train < 0.6 & p_train < 0)
  tau0_test = as.numeric(e_test < 0.6 & p_test < 0)
  
  y_train = p_train + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 1)
  
  ## causal pp estimation
  tau_cart = causal_pp(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat = exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_loo(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 1, d = 50, n = 1000")

# simulation 2
n_train = 5000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  beta_e = runif(d, -1, 1)
  beta_p = runif(d, -1, 1)
  
  X_train = matrix(runif(n_train*d), n_train, d)
  X_test = matrix(runif(n_test*d), n_test, d)
  
  ## propensity score
  temp_train = drop(X_train %*% beta_e) + 
    0.5 * beta_e[1] * X_train[,1] * X_train[,3] +
    0.7 * beta_e[2] * X_train[,2] * X_train[,4] +
    0.5 * beta_e[3] * X_train[,3] * X_train[,5] +
    0.7 * beta_e[4] * X_train[,4] * X_train[,6] +
    0.5 * beta_e[5] * X_train[,5] * X_train[,7] +
    0.5 * beta_e[1] * X_train[,1] * X_train[,6] +
    0.7 * beta_e[2] * X_train[,2] * X_train[,3] +
    0.5 * beta_e[3] * X_train[,3] * X_train[,4] +
    0.5 * beta_e[4] * X_train[,4] * X_train[,5] +
    0.5 * beta_e[5] * X_train[,5] * X_train[,6] +
    beta_e[2] * X_train[,2] * X_train[,2] +
    beta_e[4] * X_train[,4] * X_train[,4] +
    beta_e[10] * X_train[,10] * X_train[,10]
  
  temp_test = drop(X_test %*% beta_e) + 
    0.5 * beta_e[1] * X_test[,1] * X_test[,3] +
    0.7 * beta_e[2] * X_test[,2] * X_test[,4] +
    0.5 * beta_e[3] * X_test[,3] * X_test[,5] +
    0.7 * beta_e[4] * X_test[,4] * X_test[,6] +
    0.5 * beta_e[5] * X_test[,5] * X_test[,7] +
    0.5 * beta_e[1] * X_test[,1] * X_test[,6] +
    0.7 * beta_e[2] * X_test[,2] * X_test[,3] +
    0.5 * beta_e[3] * X_test[,3] * X_test[,4] +
    0.5 * beta_e[4] * X_test[,4] * X_test[,5] +
    0.5 * beta_e[5] * X_test[,5] * X_test[,6] +
    beta_e[2] * X_test[,2] * X_test[,2] +
    beta_e[4] * X_test[,4] * X_test[,4] +
    beta_e[10] * X_test[,10] * X_test[,10]
  
  e_train = exp(temp_train) / (1 + exp(temp_train))
  z_train = rbinom(rep(1,n_train), rep(1,n_train), prob = e_train)
  e_test = exp(temp_test) / (1 + exp(temp_test))
  
  ## prognostic score
  p_train = drop(X_train %*% beta_p) +
    0.5 * beta_p[1] * X_train[,1] * X_train[,3] +
    0.7 * beta_p[2] * X_train[,2] * X_train[,4] +
    0.5 * beta_p[3] * X_train[,3] * X_train[,8] +
    0.7 * beta_p[4] * X_train[,4] * X_train[,9] +
    0.5 * beta_p[5] * X_train[,8] * X_train[,10] +
    0.5 * beta_p[1] * X_train[,1] * X_train[,9] +
    0.7 * beta_p[2] * X_train[,2] * X_train[,3] +
    0.5 * beta_p[3] * X_train[,3] * X_train[,4] +
    0.5 * beta_p[4] * X_train[,4] * X_train[,8] +
    0.5 * beta_p[5] * X_train[,8] * X_train[,9] +
    beta_p[2] * X_train[,2] * X_train[,2] +
    beta_p[4] * X_train[,4] * X_train[,4] +
    beta_p[10] * X_train[,10] * X_train[,10]
  
  p_test = drop(X_test %*% beta_p) +
    0.5 * beta_p[1] * X_test[,1] * X_test[,3] +
    0.7 * beta_p[2] * X_test[,2] * X_test[,4] +
    0.5 * beta_p[3] * X_test[,3] * X_test[,8] +
    0.7 * beta_p[4] * X_test[,4] * X_test[,9] +
    0.5 * beta_p[5] * X_test[,8] * X_test[,10] +
    0.5 * beta_p[1] * X_test[,1] * X_test[,9] +
    0.7 * beta_p[2] * X_test[,2] * X_test[,3] +
    0.5 * beta_p[3] * X_test[,3] * X_test[,4] +
    0.5 * beta_p[4] * X_test[,4] * X_test[,8] +
    0.5 * beta_p[5] * X_test[,8] * X_test[,9] +
    beta_p[2] * X_test[,2] * X_test[,2] +
    beta_p[4] * X_test[,4] * X_test[,4] +
    beta_p[10] * X_test[,10] * X_test[,10]
  
  ## response
  tau0_train = as.numeric(e_train < 0.6 & p_train < 0)
  tau0_test = as.numeric(e_test < 0.6 & p_test < 0)
  
  y_train = p_train + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 1)
  
  ## causal pp estimation
  tau_cart = causal_pp(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat = exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_loo(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 2, d = 10, n = 5000")

# simulation 3
n_train = 5000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  beta_e = runif(d, -1, 1)
  beta_p = runif(d, -1, 1)
  
  X_train = matrix(runif(n_train*d), n_train, d)
  X_test = matrix(runif(n_test*d), n_test, d)
  
  ## propensity score
  temp_train = drop(X_train %*% beta_e)
  temp_test = drop(X_test %*% beta_e)
  e_train = exp(temp_train) / (1 + exp(temp_train))
  z_train = rbinom(rep(1,n_train), rep(1,n_train), prob = e_train)
  e_test = exp(temp_test) / (1 + exp(temp_test))
  
  ## prognostic score
  p_train = drop(X_train %*% beta_p)
  p_test = drop(X_test %*% beta_p)
  
  ## response
  e1 = 0.6
  p1 = 0
  
  tau0_train = rep(0, n_train)
  for (i in 1:n_train){
    if (e_train[i] <= e1 & p_train[i] <= p1) tau0_train[i] = 0
    if (e_train[i] > e1 & p_train[i] <= p1) tau0_train[i] = 1
    if (e_train[i] <= e1 & p_train[i] > p1) tau0_train[i] = 1
    if (e_train[i] > e1 & p_train[i] > p1) tau0_train[i] = 2
  }
  
  tau0_test = rep(0, n_test)
  for (i in 1:n_test){
    if (e_test[i] <= e1 & p_test[i] <= p1) tau0_test[i] = 0
    if (e_test[i] > e1 & p_test[i] <= p1) tau0_test[i] = 1
    if (e_test[i] <= e1 & p_test[i] > p1) tau0_test[i] = 1
    if (e_test[i] > e1 & p_test[i] > p1) tau0_test[i] = 2
  }
  
  y_train = p_train + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 1)
  
  ## causal pp estimation
  tau_cart = causal_pp(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat = exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_loo(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 3, d = 10, n = 5000")

# simulation 4
n_train = 5000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

beta = rep(1, d)
tau0_train = rep(0, n_train)
tau0_test = rep(0, n_test)

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  X_train = matrix(rnorm(n_train*d), n_train, d)
  X_test = matrix(rnorm(n_test*d), n_test, d)
  z_train = rep(0, n_train)
  z_train[sample(n_train, n_train/2)] = 1
  
  y_train = 1 + drop(X_train %*% beta) + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 100 - d)
  
  ## causal pp estimation
  tau_cart = causal_pp(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat = exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_loo(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 4, d = 10, n = 5000")

# simulation 5
n_train = 5000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

tau0_train = rep(0, n_train)
tau0_test = rep(0, n_test)

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  X_train = matrix(runif(n_train*d), n_train, d)
  X_test = matrix(runif(n_test*d), n_test, d)
  
  p_train = 2 * X_train[,1] - 1
  p_test = 2 * X_test[,1] - 1
  e_train = 0.25 * (1 + dbeta(X_train[,1], 2, 4))
  e_test = 0.25 * (1 + dbeta(X_test[,1], 2, 4))
  z_train = rbinom(rep(1,n_train), rep(1,n_train), prob = e_train)
  
  y_train = p_train + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 1)
  
  ## causal pp estimation
  tau_cart = causal_pp(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat= exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_loo(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 5, d = 10, n = 5000")

# simulation 6
n_train = 5000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  X_train = matrix(runif(n_train*d), n_train, d)
  X_test = matrix(runif(n_test*d), n_test, d)
  
  beta_e = runif(d, -1, 1)
  beta_p = runif(d, -1, 1)
  
  ## propensity score
  temp_train = drop(X_train %*% beta_e)
  temp_test = drop(X_test %*% beta_e)
  e_train = exp(temp_train) / (1 + exp(temp_train))
  z_train = rbinom(rep(1,n_train), rep(1,n_train), prob = e_train)
  e_test = exp(temp_test) / (1 + exp(temp_test))
  
  ## prognostic score
  p_train = drop(X_train %*% beta_p)
  p_test = drop(X_test %*% beta_p)
  
  ## treatment effect
  tau0_train = (1 + 1 / (1 + exp(-20*(e_train-1/3)))) * (1 + 1 / (1 + exp(-20*(p_train-1/3))))
  tau0_test = (1 + 1 / (1 + exp(-20*(e_test-1/3)))) * (1 + 1 / (1 + exp(-20*(p_test-1/3))))
  
  y_train = p_train + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 1)
  
  ## causal pp estimation
  tau_cart = causal_pp(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat= exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_loo(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 6, d = 10, n = 5000")

# simulation 7
n_train = 3000 ## train sample size
n_test = round(n_train / 10) ## test sample size
d = 5000 ## ambient dimension
N = 1000 ## number of repetitions

num.folds <- 10
folds <- sort(seq(n_train) %% num.folds) + 1

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

beta_e = c(0.4, 0.9, -0.4, -0.7, -0.3, 0.6, rep(0, d - 6))
beta_p = c(0.9, -0.9, 0.2, -0.2, 0.9, -0.9, rep(0, d - 6))

for (j in 1:N){
  X_train = matrix(runif(n_train*d), n_train, d)
  X_test = matrix(runif(n_test*d), n_test, d)
  
  ## propensity score
  temp_train = drop(X_train %*% beta_e)
  temp_test = drop(X_test %*% beta_e)
  e_train = exp(temp_train) / (1 + exp(temp_train))
  z_train = rbinom(rep(1,n_train), rep(1,n_train), prob = e_train)
  e_test = exp(temp_test) / (1 + exp(temp_test))
  
  ## prognostic score
  p_train = drop(X_train %*% beta_p)
  p_test = drop(X_test %*% beta_p)
  
  ## response
  tau0_train = as.numeric(e_train < 0.6 & p_train < 0)
  tau0_test = as.numeric(e_test < 0.6 & p_test < 0)
  
  y_train = p_train + tau0_train * z_train + rnorm(n_train, mean = 0, sd = 1)
  
  ## causal pp estimation
  tau_cart = causal_pp_high(X_train, y_train, z_train, round(log(n_train)))
  beta_e_test = as.vector(tau_cart$beta_e)
  temp_test = drop(cbind(1,X_test) %*% beta_e_test)
  e_hat= exp(temp_test) / (1 + exp(temp_test))
  beta_p_test = as.vector(tau_cart$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  data_test = data.frame(e_hat = e_hat, p_hat = p_hat)
  tau_ca = predict(tau_cart$model, data_test)
  mse_ca[j] = mean((tau_ca - tau0_test)^2)
  
  ## endogenous stratification
  tau_end = strat_high(X_train, y_train, z_train, 3)
  beta_p_test = as.vector(tau_end$beta_p)
  p_hat = drop(cbind(1,X_test) %*% beta_p_test)
  group = as.numeric(cut(p_hat, tau_end$interval))
  tau_st = tau_end$tau_itv[group]
  mse_st[j] = mean((tau_st - tau0_test)^2)
  
  ## causal forests
  tau_forest = causal_forest(X_train, y_train, z_train, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
  pred = predict(tau_forest, X_test)
  tau_cf = pred$predictions
  mse_cf[j] = mean((tau_cf - tau0_test)^2)
}

## plot of comparison
boxplot(mse_ca, mse_st, mse_cf, 
        names = c("PP", "ST", "CF"), 
        ylab = "MSE", main = "Scenario 7, d = 5000, n = 3000")

# bootstrap
## Scenario 1
set.seed(1)
n = 5000
d = 10

num.folds <- 10
folds <- sort(seq(n) %% num.folds) + 1

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

B = 1000

pred_pp = matrix(0, nrow = n, ncol = B)
pred_st = matrix(0, nrow = n, ncol = B)
idx_mat = matrix(0, nrow = n, ncol = B)

for (i in 1:B){
  idx = sample(n, n, replace = T)
  idx_mat[,i] = idx
  
  y_sample = y[idx]
  z_sample = z[idx]
  X_sample = X[idx,]
  
  tau_boot = causal_pp(X_sample, y_sample, z_sample, round(log(n)))
  pred_pp[,i] = tau_boot$tau
  pred_st[,i] = strat_full(X_sample, y_sample, z_sample)
}

tau_forest = causal_forest(X, y, z, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
pred = predict(tau_forest, X, estimate.variance = T)

ci_mat = matrix(0, nrow = n, ncol = 7)
ci_mat[,5] = pred$predictions - 1.96 * sqrt(pred$variance.estimates)
ci_mat[,6] = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
ci_mat[,7] = tau0

for (j in 1:n){
  ci_mat[j,1:2] = quantile(pred_pp[idx_mat == j], c(0.025,0.975))
  ci_mat[j,3:4] = quantile(pred_st[idx_mat == j], c(0.025,0.975))
}

mean(tau0 >= ci_mat[,1] & tau0 <= ci_mat[,2])
mean(tau0 >= ci_mat[,3] & tau0 <= ci_mat[,4])
mean(tau0 >= ci_mat[,5] & tau0 <= ci_mat[,6])

## Scenario 2
n = 5000
d = 10

num.folds <- 10
folds <- sort(seq(n) %% num.folds) + 1

beta_e = runif(d, -1, 1)
beta_p = runif(d, -1, 1)

X = matrix(runif(n*d), n, d)

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

e = exp(temp) / (1 + exp(temp))
z = rbinom(rep(1,n), rep(1,n), prob = e)

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

tau0 = as.numeric(e < 0.6 & p < 0)

y = p + tau0 * z + rnorm(n, mean = 0, sd = 1)

B = 1000

pred_pp = matrix(0, nrow = n, ncol = B)
pred_st = matrix(0, nrow = n, ncol = B)
idx_mat = matrix(0, nrow = n, ncol = B)

for (i in 1:B){
  idx = sample(n, n, replace = T)
  idx_mat[,i] = idx
  
  y_sample = y[idx]
  z_sample = z[idx]
  X_sample = X[idx,]
  
  tau_boot = causal_pp(X_sample, y_sample, z_sample, round(log(n)))
  pred_pp[,i] = tau_boot$tau
  pred_st[,i] = strat_full(X_sample, y_sample, z_sample)
}

tau_forest = causal_forest(X, y, z, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
pred = predict(tau_forest, X, estimate.variance = T)
ci_mat = matrix(0, nrow = n, ncol = 7)
ci_mat[,5] = pred$predictions - 1.96 * sqrt(pred$variance.estimates)
ci_mat[,6] = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
ci_mat[,7] = tau0

for (j in 1:n){
  ci_mat[j,1:2] = quantile(pred_pp[idx_mat == j], c(0.025,0.975))
  ci_mat[j,3:4] = quantile(pred_st[idx_mat == j], c(0.025,0.975))
}

mean(tau0 >= ci_mat[,1] & tau0 <= ci_mat[,2])
mean(tau0 >= ci_mat[,3] & tau0 <= ci_mat[,4])
mean(tau0 >= ci_mat[,5] & tau0 <= ci_mat[,6])

## Scenario 3
n = 5000
d = 10

num.folds <- 10
folds <- sort(seq(n) %% num.folds) + 1

beta_e = runif(d, -1, 1)
beta_p = runif(d, -1, 1)

X = matrix(runif(n*d), n, d)

temp = drop(X %*% beta_e)

e = exp(temp) / (1 + exp(temp))
z = rbinom(rep(1,n), rep(1,n), prob = e)

p = drop(X %*% beta_p)

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

B = 1000

pred_pp = matrix(0, nrow = n, ncol = B)
pred_st = matrix(0, nrow = n, ncol = B)
idx_mat = matrix(0, nrow = n, ncol = B)

for (i in 1:B){
  idx = sample(n, n, replace = T)
  idx_mat[,i] = idx
  
  y_sample = y[idx]
  z_sample = z[idx]
  X_sample = X[idx,]
  
  tau_boot = causal_pp(X_sample, y_sample, z_sample, round(log(n)))
  pred_pp[,i] = tau_boot$tau
  pred_st[,i] = strat_full(X_sample, y_sample, z_sample)
}

tau_forest = causal_forest(X, y, z, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
pred = predict(tau_forest, X, estimate.variance = T)

ci_mat = matrix(0, nrow = n, ncol = 7)
ci_mat[,5] = pred$predictions - 1.96 * sqrt(pred$variance.estimates)
ci_mat[,6] = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
ci_mat[,7] = tau0

for (j in 1:n){
  ci_mat[j,1:2] = quantile(pred_pp[idx_mat == j], c(0.025,0.975))
  ci_mat[j,3:4] = quantile(pred_st[idx_mat == j], c(0.025,0.975))
}

mean(tau0 >= ci_mat[,1] & tau0 <= ci_mat[,2])
mean(tau0 >= ci_mat[,3] & tau0 <= ci_mat[,4])
mean(tau0 >= ci_mat[,5] & tau0 <= ci_mat[,6])

## Scenario 4
n = 5000
d = 10

num.folds <- 10
folds <- sort(seq(n) %% num.folds) + 1

beta = rep(1, d)
tau0 = rep(0, n)
X = matrix(rnorm(n*d), n, d)
z = rep(0, n)
z[sample(n, n/2)] = 1

y = 1 + drop(X %*% beta) + tau0 * z + rnorm(n, mean = 0, sd = 100 - d)

B = 1000

pred_pp = matrix(0, nrow = n, ncol = B)
pred_st = matrix(0, nrow = n, ncol = B)
idx_mat = matrix(0, nrow = n, ncol = B)

for (i in 1:B){
  idx = sample(n, n, replace = T)
  idx_mat[,i] = idx
  
  y_sample = y[idx]
  z_sample = z[idx]
  X_sample = X[idx,]
  
  tau_boot = causal_pp(X_sample, y_sample, z_sample, round(log(n)))
  pred_pp[,i] = tau_boot$tau
  pred_st[,i] = strat_full(X_sample, y_sample, z_sample)
}

tau_forest = causal_forest(X, y, z, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
pred = predict(tau_forest, X, estimate.variance = T)

ci_mat = matrix(0, nrow = n, ncol = 7)
ci_mat[,5] = pred$predictions - 1.96 * sqrt(pred$variance.estimates)
ci_mat[,6] = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
ci_mat[,7] = tau0

for (j in 1:n){
  ci_mat[j,1:2] = quantile(pred_pp[idx_mat == j], c(0.025,0.975))
  ci_mat[j,3:4] = quantile(pred_st[idx_mat == j], c(0.025,0.975))
}

mean(tau0 >= ci_mat[,1] & tau0 <= ci_mat[,2])
mean(tau0 >= ci_mat[,3] & tau0 <= ci_mat[,4])
mean(tau0 >= ci_mat[,5] & tau0 <= ci_mat[,6])

## Scenario 5
n = 5000
d = 10

num.folds <- 10
folds <- sort(seq(n) %% num.folds) + 1

X = matrix(runif(n*d), n, d)

tau0 = rep(0, n)
p = 2 * X[,1] - 1
e = 0.25 * (1 + dbeta(X[,1], 2, 4))
z = rbinom(rep(1,n), rep(1,n), prob = e)

y = p + tau0 * z + rnorm(n, mean = 0, sd = 1)

B = 1000

pred_pp = matrix(0, nrow = n, ncol = B)
pred_st = matrix(0, nrow = n, ncol = B)
idx_mat = matrix(0, nrow = n, ncol = B)

for (i in 1:B){
  idx = sample(n, n, replace = T)
  idx_mat[,i] = idx
  
  y_sample = y[idx]
  z_sample = z[idx]
  X_sample = X[idx,]
  
  tau_boot = causal_pp(X_sample, y_sample, z_sample, round(log(n)))
  pred_pp[,i] = tau_boot$tau
  pred_st[,i] = strat_full(X_sample, y_sample, z_sample)
}

tau_forest = causal_forest(X, y, z, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
pred = predict(tau_forest, X, estimate.variance = T)

ci_mat = matrix(0, nrow = n, ncol = 7)
ci_mat[,5] = pred$predictions - 1.96 * sqrt(pred$variance.estimates)
ci_mat[,6] = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
ci_mat[,7] = tau0

for (j in 1:n){
  ci_mat[j,1:2] = quantile(pred_pp[idx_mat == j], c(0.025,0.975))
  ci_mat[j,3:4] = quantile(pred_st[idx_mat == j], c(0.025,0.975))
}

mean(tau0 >= ci_mat[,1] & tau0 <= ci_mat[,2])
mean(tau0 >= ci_mat[,3] & tau0 <= ci_mat[,4])
mean(tau0 >= ci_mat[,5] & tau0 <= ci_mat[,6])

## Scenario 6
n = 1000
d = 10

num.folds <- 10
folds <- sort(seq(n) %% num.folds) + 1

X = matrix(runif(n*d), n, d)

beta_e = runif(d, -1, 1)
beta_p = runif(d, -1, 1)
temp = drop(X %*% beta_e)
e = exp(temp) / (1 + exp(temp))
p = drop(X %*% beta_p)

tau0 = (1 + 1 / (1 + exp(-20*(e-1/3)))) * (1 + 1 / (1 + exp(-20*(p-1/3))))
z = rbinom(rep(1,n), rep(1,n), prob = e)

y = p + tau0 * z + rnorm(n, mean = 0, sd = 1)

B = 1000

pred_pp = matrix(0, nrow = n, ncol = B)
pred_st = matrix(0, nrow = n, ncol = B)
idx_mat = matrix(0, nrow = n, ncol = B)

for (i in 1:B){
  idx = sample(n, n, replace = T)
  idx_mat[,i] = idx
  
  y_sample = y[idx]
  z_sample = z[idx]
  X_sample = X[idx,]
  
  tau_boot = causal_pp(X_sample, y_sample, z_sample, round(log(n)))
  pred_pp[,i] = tau_boot$tau
  pred_st[,i] = strat_full(X_sample, y_sample, z_sample)
}

tau_forest = causal_forest(X, y, z, tune.parameters = "all", clusters = folds, equalize.cluster.weights = TRUE)
pred = predict(tau_forest, X, estimate.variance = T)

ci_mat = matrix(0, nrow = n, ncol = 7)
ci_mat[,5] = pred$predictions - 1.96 * sqrt(pred$variance.estimates)
ci_mat[,6] = pred$predictions + 1.96 * sqrt(pred$variance.estimates)
ci_mat[,7] = tau0

for (j in 1:n){
  ci_mat[j,1:2] = quantile(pred_pp[idx_mat == j], c(0.025,0.975))
  ci_mat[j,3:4] = quantile(pred_st[idx_mat == j], c(0.025,0.975))
}

mean(tau0 >= ci_mat[,1] & tau0 <= ci_mat[,2])
mean(tau0 >= ci_mat[,3] & tau0 <= ci_mat[,4])
mean(tau0 >= ci_mat[,5] & tau0 <= ci_mat[,6])
