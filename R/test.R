priors = prior()
priors$X0 <- priors$Theta <- prior.normal()
priors$H <- priors$Sigma_x <- prior.invgamma(2)
class(priors) <- "prior"

priors_tr <- lapply(priors, varchange)
priors_tr

lupost1 = lupost_factory(modelOU, Xsim, tree, priors)
pars_init <- c(0,0,0,0)
out = mcmc::metrop(lupost1, pars_init, nbatch = 1000, scale = 0.3)

post_pars <- out$batch

# transform back
post_pars[,2] = exp(post_pars[,2])
post_pars[,4] = exp(post_pars[,4])


hist(post_pars[,1], breaks = 20)
hist(post_pars[,2], breaks = 20)
hist(post_pars[,3], breaks = 20)
hist(post_pars[,4], breaks = 20)


dinvgamma(X, 2, log = TRUE)

Y = log(X - 0)

X = rinvgamma(10, 2)
Y = log(X)

priors$Sigma_x(X)
exp(priors$Sigma_x(X))


priors_tr$Sigma_x(log(0.5))


xgrid = seq(0.001, 3, length.out = 1000)
ygrid = log(xgrid)


plot(xgrid, exp(sapply(xgrid, priors$Sigma_x)))

#plot(ygrid, exp(sapply(ygrid, priors_tr$Sigma_x)))

priors_tr$Sigma_x

plot(exp(ygrid), exp(exp(sapply(ygrid, priors_tr$Sigma_x))))


xgrid = seq(-4, 4, length.out = 1000)
ygrid = exp(xgrid)

plot(xgrid, sapply(xgrid, dnorm))
#plot(ygrid, sapply(ygrid, dlnorm))

plot(log(ygrid), log(sapply(ygrid, dlnorm)))

Y = rlnorm(1000)
hist(Y, breaks = 50)

hist(log(Y), breaks = 50, freq = FALSE)
lines(xgrid, dnorm(xgrid))

xgrid =  seq(0.001, 5, length.out = 1000)
f1 = mcmc::metrop(function(x) {dlnorm(x)}, initial = exp(0), nbatch = 1000, scale = 0.001)
hist(log(f1$batch), breaks = 50, freq = FALSE)
lines(xgrid, dinvgamma(xgrid, 2))

f2 = mcmc::metrop(f1, nbatch = 1000, scale = 0.01)

hist(log(f2$batch), breaks = 50, freq = FALSE)
lines(xgrid, dinvgamma(xgrid, 2))

hist(exp(f2$batch), breaks = 50, freq = FALSE)
lines(xgrid, dinvgamma(xgrid, 2))




hist(Y, freq = F, breaks = 100, xlim = c(0,3))

ygrid = seq(-4, 4, length.out = 1000)

lpost = function(x){
  if (x <0){
    return(-1e5)
  }
  else{
    dnorm(log(x), log = T) - log(x)
  }
}



{
x = c(exp(0))
nsteps = 1e4
lp0 = dlnorm(x[1], log = TRUE)

for (i in 2:nsteps){
  x1 = rnorm(1, x[i-1], 0.1)
  lp1 = lpost(x1)

  if ((lp1-lp0) > log(runif(1))){
    x[i] = x1
    lp0 = lp1
  }
  else{
    x[i] = x[i-1]
  }
}

hist(x, freq = F, breaks = 50)
lines(xgrid, dlnorm(xgrid))
}


hist(log(x), freq = F, breaks = 50, xlim = c(-4,4))
lines(ygrid, dnorm(ygrid))
























