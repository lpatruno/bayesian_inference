library("rjags")

dat <- read.csv(file="pctgrowth.csv", header=TRUE)

head(dat)
str(dat)

table(dat$grp)

boxplot(y ~ grp, data=dat)

## Modeling

mod_string = " model {
for (i in 1:length(y)) {
  y[i] ~ dnorm(theta[grp[i]], 1/sig^2)
}

for (j in 1:max(grp)) {
  theta[j] ~ dnorm(mu, 1/tau^2)
}


sig ~ dgamma(1.0, 1.0)
mu ~ dnorm(0, 1/1e6)
tau ~ dgamma(2.0, 2.0/3.0)
} "

set.seed(113)

data_jags = as.list(dat)

params = c("theta", "mu", "sig", "tau")

mod = jags.model(textConnection(mod_string), data=data_jags, n.chains=3)
update(mod, 1e3)

mod_sim = coda.samples(model=mod,
                       variable.names=params,
                       n.iter=5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim))

pos_means <- colMeans(mod_csim)
means_theta <- pos_means[paste0(rep("theta[", 5), 1:5, rep("]", 5))]


# Approx posterior under cell means model
means_anova = tapply(dat$y, INDEX=dat$grp, FUN=mean)

plot(means_anova)
points(means_theta, col="red")
