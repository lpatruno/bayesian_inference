library("boot")
library("corrplot")
library("rjags")
data("urine")
?urine
head(urine)

str(urine)

dat = na.omit(urine)

pairs(dat)


Cor = cor(dat)
corrplot(Cor, type="upper", method="ellipse", tl.pos="d")
corrplot(Cor, type="lower", method="number", col="black", 
         add=TRUE, diag=FALSE, tl.pos="n", cl.pos="n")

### Variable Selection

X = scale(dat[,-1], center=TRUE, scale=TRUE)
head(X[,"gravity"])
head(X)

colMeans(X)
apply(X, 2, sd)

### Model

ddexp = function(x, mu, tau) {
  0.5 * tau * exp(-tau * abs(x - mu)) 
}
# double exponential distribution
curve(ddexp(x, mu=0.0, tau=1.0),
      from=-5.0, to=5.0, ylab="density",
      main="Double exponential\ndistribution")
# normal distribution
curve(dnorm(x, mean=0.0, sd=1.0),
      from=-5.0, to=5.0, lty=2, add=TRUE)
legend("topright", legend=c("double exponential", "normal"), lty=c(1,2), bty="n")

# Double Exponential Prior on the \beta coefficients
mod1_string = " model {
    for (i in 1:length(y)) {
        y[i] ~ dbern(p[i])
        logit(p[i]) = int + b[1]*gravity[i] + b[2]*ph[i] + b[3]*osmo[i] + b[4]*cond[i] + b[5]*urea[i] + b[6]*calc[i]
    }
    int ~ dnorm(0.0, 1.0/25.0)
    for (j in 1:6) {
        b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
    }
} "

set.seed(92)
head(X)

data_jags = list(y=dat$r, gravity=X[,"gravity"],
                 ph=X[,"ph"], osmo=X[,"osmo"],
                 cond=X[,"cond"], urea=X[,"urea"],
                 calc=X[,"calc"])

params = c("int", "b")

mod1 = jags.model(textConnection(mod1_string),
                  data=data_jags, n.chains=3)
update(mod1, 1e3)

mod1_sim = coda.samples(model=mod1,
                        variable.names=params,
                        n.iter=5e3)
mod1_csim = as.mcmc(do.call(rbind, mod1_sim))

## convergence diagnostics
plot(mod1_sim, ask=TRUE)

gelman.diag(mod1_sim)
autocorr.diag(mod1_sim)
autocorr.plot(mod1_sim)
effectiveSize(mod1_sim)

## calculate DIC
dic1 = dic.samples(mod1, n.iter=1e3)
dic1

summary(mod1_sim)

par(mfrow=c(3,2))
densplot(mod1_csim[,1:6], xlim=c(-3.0, 3.0))

colnames(X) # variable names


mod2_string = " model {
  for (i in 1:length(y)) {
    y[i] ~ dbern(p[i])
    logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
  }
  int ~ dnorm(0.0, 1.0/25.0)
  for (j in 1:3) {
    b[j] ~ dnorm(0.0, 1.0/25.0) # noninformative for logistic regression
  }
} "

mod2 = jags.model(textConnection(mod2_string), data=data_jags, n.chains=3)

update(mod2, 1e3)

mod2_sim = coda.samples(model=mod2,
                        variable.names=params,
                        n.iter=5e3)
mod2_csim = as.mcmc(do.call(rbind, mod2_sim))

plot(mod2_sim, ask=TRUE)

gelman.diag(mod2_sim)
autocorr.diag(mod2_sim)
autocorr.plot(mod2_sim)
effectiveSize(mod2_sim)

dic2 = dic.samples(mod2, n.iter=1e3)
dic2

mod3_string = " model {
    for (i in 1:length(y)) {
y[i] ~ dbern(p[i])
logit(p[i]) = int + b[1]*gravity[i] + b[2]*cond[i] + b[3]*calc[i]
}
int ~ dnorm(0.0, 1.0/25.0)
for (j in 1:3) {
b[j] ~ ddexp(0.0, sqrt(2.0)) # has variance 1.0
}
} "

mod3 = jags.model(textConnection(mod3_string), data=data_jags, n.chains=3)

update(mod3, 1e3)

mod3_sim = coda.samples(model=mod3,
                        variable.names=params,
                        n.iter=5e3)
mod3_csim = as.mcmc(do.call(rbind, mod3_sim))

plot(mod3_sim, ask=TRUE)

gelman.diag(mod3_sim)
autocorr.diag(mod3_sim)
autocorr.plot(mod3_sim)
effectiveSize(mod3_sim)

dic3 = dic.samples(mod3, n.iter=1e3)
dic3

summary(mod3_sim)

HPDinterval(mod3_csim)

dic1 - dic2
dic1 - dic3

### Prediction

pm_coef = colMeans(mod2_csim)
pm_coef

pm_Xb = pm_coef["int"] + X[,c(1,4,6)] %*% pm_coef[1:3]
phat = 1.0 / (1.0 + exp(-pm_Xb))
head(phat)

par(mfrow=c(1,1))
plot(phat, jitter(dat$r))

(tab0.5 = table(phat > 0.5, data_jags$y))

sum(diag(tab0.5)) / sum(tab0.5)

(tab0.3 = table(phat > 0.3, data_jags$y))

sum(diag(tab0.3)) / sum(tab0.3)
