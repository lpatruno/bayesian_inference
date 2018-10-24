library("car")
library("rjags")
data("Anscombe")


head(Anscombe)
pairs(Anscombe)

lin_mod <- lm(education ~ income + young + urban,
              data=Anscombe)


mod_string = " model {
    for (i in 1:length(education)) {
      education[i] ~ dnorm(mu[i], prec)
      mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*urban[i]
    }
    
    b0 ~ dnorm(0.0, 1.0/1.0e6)
    for (i in 1:3) {
      b[i] ~ dnorm(0.0, 1.0/1.0e6)
    }

    prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
      ## Initial guess of variance based on overall
      ## variance of education variable. Uses low prior
      ## effective sample size. Technically, this is not
      ## a true 'prior', but it is not very informative.
    sig2 = 1.0 / prec
    sig = sqrt(sig2)
  } "

data_jags = as.list(Anscombe)

set.seed(72)

params1 = c("b0", "b", "sig")

inits1 = function() {
  inits = list("b0" = rnorm(1, 0.0, 100.0),
               "b" = rnorm(3,0.0,100.0),
               "prec"=rgamma(1,1.0,1.0))
}

mod1 = jags.model(textConnection(mod_string), data=data_jags, inits=inits1, n.chains=3)


update(mod1, 1000) # burn-in

mod1_sim = coda.samples(model=mod1,
                        variable.names=params1,
                        n.iter=5000)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim)) # combine multiple chains

plot(mod1_csim)

gelman.diag(mod1_sim)

autocorr.plot(mod1_csim)

plot(lin_mod)


######## Quiz 7B

dic.samples(mod1, n.iter=1e5)

mod2_string = " model {
  for (i in 1:length(education)) {
    education[i] ~ dnorm(mu[i], prec)
    mu[i] = b0 + b[1]*income[i] + b[2]*young[i]
  }

  b0 ~ dnorm(0.0, 1.0/1.0e6)
  for (i in 1:2) {
    b[i] ~ dnorm(0.0, 1.0/1.0e6)
  }

  prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    ## Initial guess of variance based on overall
    ## variance of education variable. Uses low prior
    ## effective sample size. Technically, this is not
    ## a true 'prior', but it is not very informative.
  sig2 = 1.0 / prec
  sig = sqrt(sig2)
} "

set.seed(72)

data2_jags = list()
data2_jags[["education"]] <- data_jags[["education"]]
data2_jags[["income"]] <- data_jags[["income"]]
data2_jags[["young"]] <- data_jags[["young"]]

params2 = c("b0", "b", "sig")

inits2 = function() {
  inits = list("b0" = rnorm(1, 0.0, 100.0),
               "b" = rnorm(2, 0.0,100.0),
               "prec"=rgamma(1, 1.0,1.0))
}

mod2 = jags.model(textConnection(mod2_string), data=data2_jags,
                  inits=inits2, n.chains=3)


update(mod1, 1000) # burn-in

mod2_sim = coda.samples(model=mod2,
                        variable.names=params2,
                        n.iter=5000)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim)) # combine multiple chains



mod3_string = " model {
  for (i in 1:length(education)) {
    education[i] ~ dnorm(mu[i], prec)
    mu[i] = b0 + b[1]*income[i] + b[2]*young[i] + b[3]*income[i]*young[i]
  }

  b0 ~ dnorm(0.0, 1.0/1.0e6)
  for (i in 1:3) {
    b[i] ~ dnorm(0.0, 1.0/1.0e6)
  }

  prec ~ dgamma(1.0/2.0, 1.0*1500.0/2.0)
    ## Initial guess of variance based on overall
    ## variance of education variable. Uses low prior
    ## effective sample size. Technically, this is not
    ## a true 'prior', but it is not very informative.
  sig2 = 1.0 / prec
  sig = sqrt(sig2)
} "

data3_jags <- data2_jags

set.seed(72)

params3 = c("b0", "b", "sig")

inits3 = function() {
  inits = list("b0" = rnorm(1, 0.0, 100.0),
               "b" = rnorm(3,0.0,100.0),
               "prec"=rgamma(1,1.0,1.0))
}

mod3 = jags.model(textConnection(mod3_string), data=data3_jags,
                  inits=inits3, n.chains=3)


update(mod3, 1000) # burn-in

mod3_sim = coda.samples(model=mod3,
                        variable.names=params3,
                        n.iter=5000)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim)) # combine multiple chains


dic.samples(mod1, n.iter=1e5)
dic.samples(mod2, n.iter=1e5)
dic.samples(mod3, n.iter=1e5)

b1 <- mod1_csim[,c(1)]
mean(b1 > 0)

colMeans(mod1_csim)


