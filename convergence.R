source("metropolis_hastings.R")

set.seed(61)
post0 = mh(n=n, ybar=ybar, n_iter=10e3, mu_init=0.0, cand_sd=0.9)
coda::traceplot(as.mcmc(post0$mu[-c(1:500)]))

set.seed(61)
post1 = mh(n=n, ybar=ybar, n_iter=1e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post1$mu[-c(1:500)]))

set.seed(61)
post2 = mh(n=n, ybar=ybar, n_iter=100e3, mu_init=0.0, cand_sd=0.04)
coda::traceplot(as.mcmc(post2$mu))

coda::autocorr.plot(as.mcmc(post0$mu))

coda::autocorr.diag(as.mcmc(post0$mu))

coda::autocorr.plot(as.mcmc(post1$mu))

coda::autocorr.diag(as.mcmc(post1$mu))

str(post2)

coda::effectiveSize(as.mcmc(post2$mu)) # effective sample size of ~350

## thin out the samples until autocorrelation is essentially 0. This will leave you with approximately independent samples. The number of samples remaining is similar to the effective sample size.
coda::autocorr.plot(as.mcmc(post2$mu), lag.max=500)

thin_interval = 400 # how far apart the iterations are for autocorrelation to be essentially 0.
thin_indx = seq(from=thin_interval, to=length(post2$mu), by=thin_interval)
head(thin_indx)

post2mu_thin = post2$mu[thin_indx]
traceplot(as.mcmc(post2$mu))

traceplot(as.mcmc(post2mu_thin))

coda::autocorr.plot(as.mcmc(post2mu_thin), lag.max=10)

effectiveSize(as.mcmc(post2mu_thin))

length(post2mu_thin)

str(post0) # contains 10,000 iterations

coda::effectiveSize(as.mcmc(post0$mu)) # effective sample size of ~2,500

?effectiveSize

raftery.diag(as.mcmc(post0$mu))

raftery.diag(as.mcmc(post0$mu), q=0.005, r=0.001, s=0.95)

set.seed(62)
post3 = mh(n=n, ybar=ybar, n_iter=500, mu_init=10.0, cand_sd=0.3)
coda::traceplot(as.mcmc(post3$mu))

set.seed(61)

nsim = 500
post1 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=15.0, cand_sd=0.4)
post1$accpt

post2 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-5.0, cand_sd=0.4)
post2$accpt

post3 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=7.0, cand_sd=0.1)
post3$accpt

post4 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=23.0, cand_sd=0.5)
post4$accpt

post5 = mh(n=n, ybar=ybar, n_iter=nsim, mu_init=-17.0, cand_sd=0.4)
post5$accpt

pmc = mcmc.list(as.mcmc(post1$mu), as.mcmc(post2$mu), 
                as.mcmc(post3$mu), as.mcmc(post4$mu), as.mcmc(post5$mu))
str(pmc)

coda::traceplot(pmc)

coda::gelman.diag(pmc)

coda::gelman.plot(pmc)

nburn = 1000 # remember to discard early iterations
post0$mu_keep = post0$mu[-c(1:1000)]
summary(as.mcmc(post0$mu_keep))

mean(post$mu_keep > 1.0) # posterior probability that mu  > 1.0
