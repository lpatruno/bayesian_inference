library(rjags)
library("MASS")
data("OME")
?OME # background on the data
head(OME)

plot(table(OME$ID))

length(unique(OME$ID))
nrow(OME)

dat = subset(OME, OME != "N/A")
dat$OME = factor(dat$OME) # relabel OME
dat$ID = as.numeric(factor(dat$ID)) # relabel ID so there are no gaps in numbers (they now go from 1 to 63)

## Original reference model and covariate matrix
mod_glm = glm(Correct/Trials ~ Age + OME + Loud + Noise, data=dat, weights=Trials, family="binomial")
X = model.matrix(mod_glm)[,-1]

## Original model (that needs to be extended)
mod_string = " model {
for (i in 1:length(y)) {
y[i] ~ dbin(phi[i], n[i])
logit(phi[i]) = a[ID[i]] + b[1]*Age[i] + b[2]*OMElow[i] + b[3]*Loud[i] + b[4]*Noiseincoherent[i]
}

for (j in 1:max(ID)) {
a[j] ~ dnorm(mu, 1.0/tau)
}

mu ~ dnorm(0, 1/10^2)
tau ~ dgamma(2, 2)
for (j in 1:4) {
b[j] ~ dnorm(0.0, 1.0/4.0^2)
}

} "

data_jags = as.list(as.data.frame(X))
data_jags$y = dat$Correct
data_jags$n = dat$Trials
data_jags$ID = dat$ID
str(data_jags)

params = c("a", "b", "mu", "tau")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod1_sim = coda.samples(model=mod,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

plot(mod1_sim, ask=TRUE)
gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)

dic <- dic.samples(mod, n.iter=5e3)
dic
