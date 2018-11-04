library(fda)


#####################
## Generating data ##
#####################

set.seed(30)
n <- 50
m <- 2
k <- 30
Xtrue <- sort(rnorm(n,0,1))
GRID <- 200

mfunction <- function(x){
            f <- sin(pi*(x/2))/(1+(2*(x^2)*(sign(x)+1)))
            return(f)}

Y <- mfunction(Xtrue) + rnorm(n,mean=0,sd=.3)

Wmatrix <- sapply(Xtrue,rnorm,n=m,sd=0.8) #each column is for one i ### In this case we have replicates for W
dim(Wmatrix)


#### End of generating data

#############################################################################################
## First Part of ICM algorithm: finding sigma2_U, sigma2_E and naive m function (or g^(0)) ##   
#############################################################################################

### estimate sigma2_U, sigma2_X and mu_X
sigma2_U <- mean(apply(Wmatrix,2,var))  ##pooled variance

## starting values for X and g 
X0 <- colMeans(Wmatrix)

## Grid and knots, Xp is the grid of values over x... keep track of the spline values
Xp <- seq(from=min(X0),to=max(X0),length=GRID)
k=30
Quantiles <- round(seq(from=0,to=1,length=k+2),4)[-c(1,k+2)]
KNOTS <- as.vector(quantile(X0,probs=Quantiles))
rng <- c(-5,5)

## Fitting a naive smoothing splines to obtain g^(0)
#PenaltyMatrixD <- getbasispenalty(create.bspline.basis(rangeval=rng, norder=4, breaks=KNOTS), Lfdobj=int2Lfd(2))
ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=KNOTS)

Lambda <- df2lambda(X0, ybasis, Lfdobj=int2Lfd(2), df=6.112844)
yfdPar <- fdPar(ybasis, Lfdobj=int2Lfd(2), lambda=Lambda)
yfdu <- smooth.basis(X0,Y,yfdPar)$fd

ghat <- eval.fd(X0,yfdu) # evaluating g^(0) at X0
ghat2 <- eval.fd(Xp,yfdu) # evaluating g^(0) at a finer grid

plot(Xtrue,Y,xlim=c(-3.5,3),ylab='m(x)',xlab='x')
lines(X0,Y,pch=20,type='p')
lines(Xp,mfunction(Xp),lwd=2)
lines(Xp,ghat2,lty=2,lwd=2)

### I think I used the code below to obtain the initial smoothing paramater 
## initial smoothing
#ghat <- smooth.spline(X0,Y,all.knots=T,cv=T)
#ghat <- smooth.spline(X0,Y,nknots=k)
#DF <- ghat$df
#ghat$fit$knot
#ghat$lambda
#spar <- ghat$spar
#ghat$cv.crit
#sum(ghat$lev)
#length(mhat$x)


### the DF = 6.112844 can be obtain by calculating A(alpha)
sigma2_E <- sum((ghat-Y)^2)/(n-6.112844)

### initial values for mu_X and sigma2_X
mu_X <- mean(c(Wmatrix[1,],Wmatrix[2,]))
sigma2_X <- var(c(Wmatrix[1,],Wmatrix[2,]))

#################################################################################################
## Second part of ICM: iterating between 1.) find X given g fixed and 2.) find g with X fixed. ##
#################################################################################################

XBig <- NULL
iter <- 20
for(i in 1:iter){
         
### Very simple grid search looking at each coordinate 
Xnew <- vector(mode='numeric',length=n)
  for (j in 1:n){                     
  posterior <- (-1/(2*sigma2_E))*((Y[j]-ghat2)^2) + 
                     ((-1/(2*sigma2_U))*(((Wmatrix[,j][1]-Xp)^2)+((Wmatrix[,j][2]-Xp)^2))) + 
                     ((-1/(2*sigma2_X))*((Xp-mu_X)^2)) ### this is only the part of the posterior that depends on X
  Xnew[j] <- Xp[which(posterior==max(posterior))]
  }
  

XBig <- cbind(XBig,Xnew)

yfdunew <- smooth.basis(Xnew,Y,yfdPar)$fd ### I did not change lambda, same as the initial one 
gnew <- eval.fd(Xnew,yfdunew)

#gnew <- smooth.spline(Xnew,Y,spar=spar,nknots=k)
c1 <- mean(abs(eval.fd(Xp,yfdunew)-ghat2)) ### calculating delta for convergence 
print(c1)

ghat2 <- eval.fd(Xp,yfdunew)

if(i==1){
  ICM1 <- ghat2  
}

}


lines(Xp,ICM1,lty=3,lwd=2,col="red")
lines(Xp,ghat2,col='gray',lwd=2)
legend("topright",c('True Curve','Naive','ICM-1','ICM'),lty=c(1,2,3,1),col=c(1,1,1,'gray'),lwd=c(2,2,2,2))


