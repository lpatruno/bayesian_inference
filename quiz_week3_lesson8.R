library("rjags")

data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)


mod_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dnorm(mu[grp[i]], prec)
}

for (j in 1:3) {
mu[j] ~ dnorm(0.0, 1.0/1.0e6)
}

prec ~ dgamma(5/2.0, 5*1.0/2.0)
sig = sqrt( 1.0 / prec )
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
              grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(1,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

dic1 <- dic.samples(mod, n.iter=1e3)

pm_params = colMeans(mod_csim)

### HPD interval /mu_3 - /mu_1

HPDinterval(mod_csim[,3] - mod_csim[,1])


mod_string = " model {
  for (i in 1:length(y)) {
    y[i] ~ dnorm(mu[grp[i]], prec[grp[i]])
  }
  
  for (j in 1:3) {
    mu[j] ~ dnorm(0.0, 1.0/1.0e6)
  }
  
  for (j in 1:3) {
    prec[j] ~ dgamma(5/2.0, 5*1.0/2.0)
    sig[j] = sqrt( 1.0 / prec[j] )
  }
} "

set.seed(82)
str(PlantGrowth)
data_jags = list(y=PlantGrowth$weight, 
                 grp=as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
    inits = list("mu"=rnorm(3,0.0,100.0), "prec"=rgamma(3,1.0,1.0))
}

mod = jags.model(textConnection(mod_string), data=data_jags, inits=inits, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains

mod_csim

dic2 <- dic.samples(mod, n.iter=1e3)

pm_params = colMeans(mod_csim)

pm_params

PlantGrowth$group
as.numeric(PlantGrowth$group)

dic1 - dic2


mod_cm = lm(weight ~ -1 + group, data=PlantGrowth)
summary(mod_cm)
