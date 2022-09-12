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
n = 5000 ## sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  X = matrix(rnorm(n*d), n, d)
  z = rep(0, n)
  z[sample(n, n/2)] = 1

  y = 1 + drop(X %*% beta) + tau0 * z + rnorm(n, mean = 0, sd = 100 - d)

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
        ylab = "MSE", main = "Scenario 4, d = 10, n = 5000")

# simulation 5
n = 5000 ## sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

tau0 = rep(0, n)

for (j in 1:N){
  X = matrix(runif(n*d), n, d)
  
  p = 2 * X[,1] - 1
  e = 0.25 * (1 + dbeta(X[,1], 2, 4))
  z = rbinom(rep(1,n), rep(1,n), prob = e)
  
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
        ylab = "MSE", main = "Scenario 5, d = 10, n = 5000")

# simulation 6
n = 5000 ## sample size
d = 10 ## ambient dimension
N = 1000 ## number of repetitions

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

for (j in 1:N){
  X = matrix(runif(n*d), n, d)
  
  beta_e = runif(d, -1, 1)
  beta_p = runif(d, -1, 1)
  
  ## propsensity score
  temp = drop(X %*% beta_e)
  e = exp(temp) / (1 + exp(temp))
  z = rbinom(rep(1,n), rep(1,n), prob = e)
  
  ## prognostic score
  p = drop(X %*% beta_p)
  
  ## treatment effect
  tau0 = (1 + 1 / (1 + exp(-20*(e-1/3)))) * (1 + 1 / (1 + exp(-20*(p-1/3))))
  
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
        ylab = "MSE", main = "Scenario 6, d = 10, n = 5000")

# simulation 7
n = 3000 ## sample size
d = 5000 ## ambient dimension
N = 1000 ## number of repetitions

set.seed(123)
mse_ca = rep(0, N)
mse_st = rep(0, N)
mse_cf = rep(0, N)

beta_e = c(0.4, 0.9, -0.4, -0.7, -0.3, 0.6, rep(0, d - 6))
beta_p = c(0.9, -0.9, 0.2, -0.2, 0.9, -0.9, rep(0, d - 6))

for (j in 1:N){
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
        ylab = "MSE", main = "Scenario 7, d = 5000, n = 3000")

# rhc
require(rpart.plot)
Hmisc::getHdata(rhc,where="https://hbiostat.org/data/repo/")
dim(rhc) # number of obs
summary(rhc$swang1) # number of treated vs. untreated
rhc <- rhc[,names(rhc)[order(names(rhc))]] # order column names alphabetically
for (cc in 1:ncol(rhc)) {
  cat(names(rhc)[cc], "----------------------------------------------------\n")
  print(summary(rhc[,cc]))
}
rhc$cat2[is.na(rhc$cat2)] <- "Missing" # remove NAs from cat2
## relabelled based on Hirano and Imbens (2002) Table 1
rhc.dt <- data.frame(
  "treat"=rhc$swang1=="RHC",
  "Y"=rhc$death=="Yes",
  "age"=rhc$age,
  "sex"=rhc$sex=="Female",
  # dummy variables for race
  "raceblack"=rhc$race=="black",
  "raceother"=rhc$race=="other",
  "edu"=rhc$edu,
  # dummy variables for income
  "income1"=rhc$income=="$11-$25k",
  "income2"=rhc$income=="$25-$50k",
  "income3"=rhc$income=="> $50k",
  # dummy variables for Medical insurance
  "ins_care"=rhc$ninsclas=="Medicare",
  "ins_pcare"=rhc$ninsclas=="Private & Medicare",
  "ins_caid"=rhc$ninsclas=="Medicaid",
  "ins_no"=rhc$ninsclas=="No insurance",
  "ins_carecaid"=rhc$ninsclas=="Medicare & Medicaid",
  # dummy variables for Primary disease category
  "cat1_copd"=rhc$cat1=="COPD",
  "cat1_mosfsep"=rhc$cat1=="MOSF w/Sepsis",
  "cat1_mosfmal"=rhc$cat1=="MOSF w/Malignancy",
  "cat1_chf"=rhc$cat1=="CHF",
  "cat1_coma"=rhc$cat1=="Coma",
  "cat1_cirr"=rhc$cat1=="Cirrhosis",
  "cat1_lung"=rhc$cat1=="Lung Cancer",
  "cat1_colon"=rhc$cat1=="Colon Cancer",
  # dummy variables for Secondary disease category
  "cat2_mosfsep"=rhc$cat2=="MOSF w/Sepsis",
  "cat2_coma"=rhc$cat2=="Coma",
  "cat2_mosfmal"=rhc$cat2=="MOSF w/Malignancy",
  "cat2_lung"=rhc$cat2=="Lung Cancer",
  "cat2_cirr"=rhc$cat2=="Cirrhosis",
  "cat2_colon"=rhc$cat2=="Colon Cancer", 
  # Categories of admission diagnosis
  "resp"=rhc$resp=="Yes",
  "card"=rhc$card=="Yes",
  "neuro"=rhc$neuro=="Yes",
  "gastr"=rhc$gastr=="Yes",
  "renal"=rhc$renal=="Yes",
  "meta"=rhc$meta=="Yes",
  "hema"=rhc$hema=="Yes",
  "seps"=rhc$seps=="Yes",
  "trauma"=rhc$trauma=="Yes",
  "ortho"=rhc$ortho=="Yes",
  rhc$das2d3pc,
  "dnr1"=rhc$dnr1=="Yes",
  "ca_yes"=rhc$ca=="Yes",
  "ca_meta"=rhc$ca=="Metastatic",
  rhc$surv2md1,
  rhc$aps1,
  rhc$scoma1,
  rhc$wtkilo1,
  rhc$temp1,
  rhc$meanbp1,
  rhc$resp1,
  rhc$hrt1,
  rhc$pafi1,
  rhc$paco21,
  rhc$ph1,
  rhc$wblc1,
  rhc$hema1,
  rhc$sod1,
  rhc$pot1,
  rhc$crea1,
  rhc$bili1,
  rhc$alb1,
  rhc$urin1, # not included in Hirano and Imbens
  # Categories of comorbidities illness:
  rhc$cardiohx,
  rhc$chfhx,
  rhc$dementhx,
  rhc$psychhx,
  rhc$chrpulhx,
  rhc$renalhx,
  rhc$liverhx,
  rhc$gibledhx,
  rhc$malighx,
  rhc$immunhx,
  rhc$transhx,
  rhc$amihx,
  "wt0"=rhc$wtkilo1==0
)
names(rhc.dt) <- sapply(names(rhc.dt), function(x) {
  ifelse(grepl("rhc.",x),strsplit(x,"rhc.")[[1]][2], x)
})
names(rhc.dt)

l <- rhc.dt[,3:ncol(rhc.dt)]
Lnames <- data.frame(cbind("L.idx"=1:ncol(l),"L.names"=names(l)))
# means of the covariates by treatment status
## Hirano and Imbens (2002) Table 2
l.ttmtdiff <- list()
for (li in 1:ncol(l)) {
  li.a0 <- l[rhc.dt$treat==FALSE,li]
  li.a1 <- l[rhc.dt$treat==TRUE,li]
  l.ttmtdiff[[li]] <- c(mean(li.a0,na.rm=TRUE),mean(li.a1,na.rm=TRUE),
                        t.test(x=li.a1,y=li.a0)[["statistic"]])
}
l.ttmtdiff <- do.call(rbind,l.ttmtdiff)
rownames(l.ttmtdiff) <- names(l)

complete_cases <- rowSums(is.na(l))==0
l <- l[complete_cases,]
## remove covariates with singular value in data
l <- l[,-which(apply(l, 2, function(x) length(unique(x)))==1)]

# convert any logical covariates to integer
l <- apply(l, 2, function(x) as.numeric(x))
summary(l)
l <- as.data.frame(l)
colnames(l) <- NULL
mydata <- data.frame("i"=1:nrow(l),
                     "L"=l,
                     "treat"=as.integer(rhc.dt$treat[complete_cases]),
                     "Y"=as.integer(rhc.dt$Y[complete_cases]))


if (FALSE) {
  #### testing cov.sel
  var.list <- paste0("L.",1:p)
  cov.sel(T = mydata$treat, Y = mydata$Y, X = mydata[,var.list], 
          type = "dr", alg = 1, trace = 0)
  cov.sel(T = mydata$treat, Y = mydata$Y, X = mydata[,var.list], 
          type = "dr", alg = 2, trace = 0)
  cov.sel(T = mydata$treat, Y = mydata$Y, X = mydata[,var.list], 
          type = "np", alg = 1, trace = 0)
  cov.sel(T = mydata$treat, Y = mydata$Y, X = mydata[,var.list], 
          type = "np", alg = 2, trace = 0)
}

X = as.matrix(mydata[,2:73])
z = mydata$treat
y = mydata$Y
           
rpart.plot(tau_cart$model)
           
# nmes
## preprocess the data
smoke <- read.csv("nmesdata.txt", header = T)
smoke <- subset(smoke, TOTALEXP > 0)
smoke <- subset(smoke, select = c(TOTALEXP, packyears, LASTAGE, AGESMOKE,
                    MALE, RACE3, marital, educate,
                     SREGION, POVSTALB, beltuse, yearsince))
smoke <- na.omit(smoke)
smoke$weight <- nmes_data[row.names(smoke),]$HSQACCWT
smoke <- na.omit(smoke)

X = as.matrix(smoke[,-c(1,2)])
y = log(smoke$TOTALEXP)
z = as.numeric(smoke$packyears > 17)

rpart.plot(tau_cart$model)

