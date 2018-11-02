#install.packages("COUNT")
library("COUNT")
library("rjags")
data("badhealth")

head(badhealth)

b0 <- 1.5
b1 <- -.3
b2 <- 1.0
x1 <- .8
x2 <- 1.2

log_lam <- b0 + b1*x1 + b2*x2
(lam <- exp(log_lam))


mod_string = " model {
  for (i in 1:length(numvisit)) {
    numvisit[i] ~ dpois(lam[i])
    log(lam[i]) = int + b_badh*badh[i] + b_age*age[i] + b_intx*age[i]*badh[i]
  }
  
  int ~ dnorm(0.0, 1.0/1e6)
  b_badh ~ dnorm(0.0, 1.0/1e4)
  b_age ~ dnorm(0.0, 1.0/1e4)
  b_intx ~ dnorm(0.0, 1.0/1e4)
  } "

set.seed(102)

data_jags = as.list(badhealth)

params = c("int", "b_badh", "b_age", "b_intx")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

colMeans(mod_csim)

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)
dic


mod_string2 = " model {
  for (i in 1:length(numvisit)) {
    numvisit[i] ~ dpois(lam[i])
    log(lam[i]) = int + b_badh*badh[i] + b_age*age[i]
  }
  
  int ~ dnorm(0.0, 1.0/1e6)
  b_badh ~ dnorm(0.0, 1.0/1e4)
  b_age ~ dnorm(0.0, 1.0/1e4)
  } "

set.seed(102)

params2 = c("int", "b_badh", "b_age")

mod2 = jags.model(textConnection(mod_string2), data=data_jags, n.chains=3)
update(mod2, 1e3)

mod_sim2 = coda.samples(model=mod2,
                       variable.names=params2,
                       n.iter=5e3)
mod_csim2 = as.mcmc(do.call(rbind, mod_sim2))

dic2 = dic.samples(mod2, n.iter=1e3)
dic2

dic


ppois(21, 15*2)

plot(table(rpois(1e3, 30)))


dat = read.csv(file="callers.csv", header=TRUE)

head(dat)

boxplot(calls/days_active ~ isgroup2, data=dat)


mod_string = " model {
  for (i in 1:length(calls)) {
    calls[i] ~ dpois(lam[i])
    log(lam[i]) = b0 + b[1] * age[i] + b[2] * isgroup2[i]
  }
  
  b0 ~ dnorm(0.0, 1.0/1e3)
  for (i in 1:2){
    b[i] ~ dnorm(0.0, 1.0 / 1e3)
  }
} "

set.seed(102)

data_jags = as.list(dat)

params = c("b0", "b")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

summary(mod_sim)

## convergence diagnostics
plot(mod_sim)

gelman.diag(mod_sim)
autocorr.diag(mod_sim)
autocorr.plot(mod_sim)
effectiveSize(mod_sim)

## compute DIC
dic = dic.samples(mod, n.iter=1e3)
dic


## Compute predicted mean number of calls
pmean_coef <- colMeans(mod_csim)
X <- as.matrix(dat[,c(4, 3)])
head(X)
llam_hat = pmean_coef["b0"] + X %*% pmean_coef[c("b[1]", "b[2]")]
lam_hat = exp(llam_hat)
hist(lam_hat)

resid = dat$calls - lam_hat
plot(resid)
plot(lam_hat, dat$calls)

mean(mod_csim[,c("b[2]")] > 0)





