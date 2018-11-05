require(mgcv)
require(rjags)
require(R2jags)
require(visreg)
require(cplm)
require(ggplot2)
source("jagsresults.R")
#####################
## Generating data ##
#####################

set.seed(42)
nobs <- 500

sdobsx <- 0.1 
truex <- sort(rnorm(nobs,0,2)) # true covariate
errx <- abs(rnorm(nobs,0,sdobsx)) # measurement errors of known variance
obsx <- rnorm(nobs,truex, errx) # observed covariate

sdy <- 0.3 # intrinsic scatter around truey
mfunction <- function(x){
  f <- sin(pi*(x/2))/(1+(2*(x^2)*(sign(x)+1)))
  return(f)}

xb <- mfunction(truex) # non-linear predictor

Y <- rnorm(nobs, xb, sdy) # true y, for now ignoring measurement errors in Y

dat <- data.frame(x=obsx,y=Y)
#### End of generating data


#basisobj = create.bspline.basis(range(obsx),nbasis=10)
#bb <- getbasismatrix(obsx, basisobj, nderiv=0)

bb <- gam(rep(1,nobs)~s(sort(obsx),k=10), fit=FALSE)$X
#bb <- bs(obsx, df=10,intercept=T)

jags.dat <- list( n = nobs,
                  y = Y,
                  X = bb
                  
)




load.module("glm") ## improved samplers for GLMs often worth loading

Model <- "model{

mu <- X %*% b ## expected response

for (i in 1:n) {
y[i] ~ dnorm(mu[i],tau) 
} 

## response 
scale <- 1/tau ## convert tau to standard GLM scale
tau ~ dgamma(.05,.005) ## precision parameter prior 
## Parametric effect priors CHECK tau=1/79^2 is appropriate!
for (i in 1:1) { b[i] ~ dnorm(0,0.00011) } 
## prior for s(x)... 
for (i in 2:9) { b[i] ~ dnorm(0, lambda[1]) }
for (i in 10:10) { b[i] ~ dnorm(0, lambda[2]) }
## smoothing parameter priors CHECK...
for (i in 1:2) {  
lambda[i] ~ dgamma(.05,.005)
#rho[i] <- log(lambda[i])
}
}"


jm  <- jags(data = jags.dat,
            parameters.to.save=c("mu"),
            model.file = textConnection(Model),
            n.chains=3,
            n.iter=1000,
            n.burnin=500)


fit <- jagsresults(x = jm, params = c("mu"),probs = c(0.16, 0.5, 0.84))
fit_data <- data.frame(x=obsx,mean=fit[,"50%"],lwr1=fit[,"16%"],upr1=fit[,"84%"])

ggplot(data=dat,aes(x=x,y=y)) +
  geom_point() +
  geom_ribbon(data=fit_data,aes(x=x,ymin=lwr1, ymax=upr1,y=NULL),fill=c("gray50"),alpha=0.5) +
  geom_line(data=fit_data,aes(x=x,y=mean),color="cyan2") + 
  theme_bw()