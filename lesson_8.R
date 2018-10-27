library("rjags")

data("PlantGrowth")
?PlantGrowth
head(PlantGrowth)

boxplot(weight ~ group, data=PlantGrowth)

### Modeling

## Reference model
lmod = lm(weight ~ group, data=PlantGrowth)
summary(lmod)

anova(lmod)

plot(lmod) # for graphical residual analysis

model.matrix(lmod)


## Cell means model in jags

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
data_jags = list(y = PlantGrowth$weight, 
                 grp = as.numeric(PlantGrowth$group))

params = c("mu", "sig")

inits = function() {
  inits = list("mu" = rnorm(3, 0.0, 100.0),
               "prec" = rgamma(1, 1.0, 1.0))
}

mod = jags.model(textConnection(mod_string),
                 data = data_jags,
                 inits = inits,
                 n.chains = 3)
update(mod, 1e3)

mod_sim = coda.samples(model = mod,
                       variable.names = params,
                       n.iter = 5e3)
mod_csim = as.mcmc(do.call(rbind, mod_sim)) # combined chains


## Visualize model convergence

plot(mod_sim)
gelman.diag(mod_sim)
autocorr.diag(mod_sim)
effectiveSize(mod_sim)
pm_params = colMeans(mod_csim) # Posterior mean of parameters
pm_params

## Visualize residuals to ensure our choice of model is correct
yhat = pm_params[1:3][data_jags$grp]
resid = data_jags$y - yhat
plot(resid)
plot(yhat, resid)


## Results
summary(mod_sim)

HPDinterval(mod_csim)

# There is a high posterior probability that the mean yield for 
# treatment 2 is greater than the mean yield for the control group.
mean(mod_csim[,3] > mod_csim[,1])


# What is the posterior probability that the increase is at least 10%?
mean(mod_csim[,3] > 1.1 * mod_csim[,1])

