library(fda)
require(ggplot2)
require(lessR)
require(BayesianTools)

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
#sdobsy <- 2.5

mfunction <- function(x){
  f <- sin(pi*(x/2))/(1+(2*(x^2)*(sign(x)+1)))
  return(f)}

xb <- mfunction(truex) # non-linear predictor

Y <- rnorm(nobs, xb, sdy) # true y, for now ignoring measurement errors in Y
#### End of generating data

obsdata <- data.frame(x=obsx,y=Y)

b0 <- gam(Y~s(obsx,k=10),method="REML")
visreg(b0,points=list(cex=0.5))


xx <- seq(min(truex),max(truex),length.out = nobs)

X <- gam(rep(1,nobs)~s(truex,k=5), fit = FALSE)$X

likelihood <- function(par){
  ######  Latent variables ###### 
  x1 = par[1:nobs]
  sd_y = par[nobs + 1]
  mu = par[(nobs + 2):(2*nobs + 1)]
  b = par[(2*nobs + 2):(2*nobs + 6)]
 # lambda = par[(2*nobs + 12):(2*nobs + 13)]
  
  ######  GAM likelihood  ######  
#  X <- bs(x1, 30)
  
  mu <- X %*% b
  lly <- sum(dnorm(Y,mean = mu,sd = sd_y,log=T))
  
  llb <- sum(dnorm(b,0,10,log=T))

#  lllambda <- sum(dgamma(lambda[1:2],0.05,0.005,log=T))
  ######  Measurement errors in covariates  ###### 
  llx1  <- sum(dnorm(obsx,mean = x1, sd = errx,log = T))  

  ###### Prediction ###### 
  #  XX <- gam(rep(1,nobs)~s(xx,k=10), fit=FALSE)$X
  # ytrue =  XX %*% b
  ######  End block  ######  
return(lly + llx1 + llb)
}  

low <- c(obsx - errx,0.01,Y - 2*sdy,rep(-5,5))

up <-  c(obsx + errx,10,Y + 2*sdy,rep(5,5))

prior <-  createUniformPrior(low, up, best = NULL)

setup <- createBayesianSetup(likelihood = likelihood,prior = prior,
                             names = c(to("x1", nobs),"sdy",to("ytrue", nobs),
                                       to("b", 5)))

settings <- list(iterations = 5E4,thin = 1,
                 burnin = 1E3, message = T,nrChains = 1)

system.time(
  res <- runMCMC(bayesianSetup = setup,  settings = settings,sampler = "DREAMzs")
)

codaObject = getSample(res, start = 1E3, coda = TRUE)

getmcmc_var <- function(outjags=outjags,vars = vars){
  as.data.frame(do.call(rbind, outjags[,vars]))
}
sDat <- getmcmc_var(codaObject,vars = to("ytrue", nobs))

Qdat <- apply(sDat, 2, quantile)

ylow <-  Qdat["25%",]
ymean <- Qdat["50%",]
yup <-   Qdat["75%",]

comb <- data.frame(x=obsx,yfit=ymean,ylow =ylow,yup =yup )

ggplot(comb, aes(x=x, y=ymean)) +
  geom_point(data=obsdata,aes(x=x,y=y))+
  geom_line() +
 # geom_ribbon(aes(x=obsx,ymin=ylow, ymax=yup,y=NULL),  fill = c("#fdae6b"),show.legend=FALSE) +
  xlab("x")+ylab("y") +
  theme(legend.text = element_text(colour="gray40"), text = element_text(size=20), legend.position=c(0.1,0.75), axis.line = element_line(color = 'black')) +
  theme_bw() +
  scale_alpha(guide="none")



PSplineSetup <- function(x,
                         k = 15, spline.deg = 3, diff.ord = 2) {
  x.min = min(x)
  x.max = max(x)
  # B-spline basis
  dx <- (x.max-x.min)/(k-spline.deg)
  knots <- seq(x.min-spline.deg*dx, x.max+spline.deg*dx, by = dx)
  B <- spline.des(knots = knots, x = x, ord = spline.deg+1, outer.ok = TRUE)$design
  # Difference operator matrix
  D <- diff(diag(k), diff = diff.ord)
  # Re-parameterize B and D into X (fixed effects) and Z (random effects) 
  X <- B%*%outer(knots[1:k], 0:(diff.ord-1), "^")
  Z <- B%*%t(D)%*%solve(tcrossprod(D))
  return(list(X = X, Z = Z, n = nrow(X), q = ncol(X), m = ncol(Z)))
}


