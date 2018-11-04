mfunction <- function(x){
            f <- sin(pi*(x/2))/(1+(2*(x^2)*(sign(x)+1)))
            return(f)}
            
mfunction2 <- function(x){
            f <- sin(pi*(x/2))/(1+((2*x)*(sign(x)+1)))
            return(f)}
            
cvFunction <- function(x,y){
#indices <- floor(seq(from=1,to=n,length=k))
Knots <- KNOTS
#rng <- c(-18,18)
ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=Knots)
lambdagrid <- seq(from=0.001, to=1,by=0.005)
cv <- vector(mode='numeric',length=length(lambdagrid))
for(i in 1:length(lambdagrid)){
    yfd <- fdPar(ybasis, Lfdobj=int2Lfd(2), lambdagrid[i])
    cv[i] <- lambda2gcv(log10lambda=log10(lambdagrid[i]), argvals=x, y=y, fdParobj=yfd)
    }
l <- lambdagrid[which(cv==min(cv))]
DF <- lambda2df(x, ybasis,Lfdobj=int2Lfd(2),lambda=l)
out <- list(lambda=l,DF=DF)
out
}

smooth.Camila <- function(x,y,Lambda){
#indices <- floor(seq(from=1,to=n,length=k))
#Knots <- sort(x[indices])
Knots <- KNOTS
#rng <- c(-18,18)
ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=Knots)
yfdPar <- fdPar(ybasis, Lfdobj=int2Lfd(2), lambda=Lambda)
yfdu <- smooth.basis(x,y,yfdPar)$fd
#out <- yfdu
#if (Penalty==TRUE){
#D <-  smooth.basis(x,y,yfdPar)$penmat
#B <- getbasismatrix(x, ybasis, nderiv=0)
#alpha <- Lambda
#Betahat <- solve(t(B)%*%B+alpha*D)%*%t(B)%*%y
#penalty <- t(Betahat)%*%D%*%Betahat
#out <- list(yfd=yfdu,Pen=penalty)
return(yfdu)
#}
#out
}


gfunction <- function(x,y,k,Gam,s2E){
Lambda <- Gam*s2E
indices <- floor(seq(from=1,to=n,length=k))
Knots <- sort(x[indices])
rng <- c(min(x),max(x))
ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=Knots)
#yfdPar <- fdPar(ybasis, Lfdobj=int2Lfd(2), lambda=Lambda)
#yfdu <- smooth.basis(x,y,yfdPar)$fd
#out <- yfdu
#if (Penalty==TRUE){
D <- getbasispenalty(ybasis, Lfdobj=int2Lfd(2))
#D <-  smooth.basis(x,y,yfdPar)$penmat
B <- getbasismatrix(x, ybasis, nderiv=0)
#A <-  B%*%solve(t(B)%*%B+Lambda*D)%*%t(B)
#A <-  B%*%chol2inv(chol(t(B)%*%B+Lambda*D))%*%t(B)
#rmvnorm(1,mean=A%*%y,sigma=s2E*A)
Betasample <- rmvnorm(1,mean=chol2inv(chol(t(B)%*%B+Lambda*D))%*%t(B)%*%y,sigma=s2E*chol2inv(chol(t(B)%*%B+Lambda*D)))
gsample <- B%*%t(Betasample)
#Betahat <- solve(t(B)%*%B+alpha*D)%*%t(B)%*%y
penalty <- (Betasample)%*%D%*%t(Betasample)
out <- list(yfd=gsample,Pen=penalty)
#}
out
}


getBMatrix <- function(x){
#indices <- floor(seq(from=1,to=n,length=k))
#Knots <- sort(x[indices])
#rng <- c(min(-18),max(18))
ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=KNOTS)
B <- getbasismatrix(x, ybasis, nderiv=0)
return(B)}

#getPenMatrix <- function(x){
##indices <- floor(seq(from=1,to=n,length=k))
#Knots <- sort(x[indices])
##rng <- c(min(-18),max(18))
#ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=Knots)
#D <- getbasispenalty(ybasis, Lfdobj=int2Lfd(2))
#return(D)}

gfunction2 <- function(x,y,Gam,s2E){
Lambda <- Gam*s2E
D <- PenaltyMatrixD
B <- getbasismatrix(x, ybasis, nderiv=0)
#Betasample <- solve(((t(B)%*%B)+(Lambda*D)))%*%t(B)%*%y
Betasample <- rmvnorm(1,mean=solve(((t(B)%*%B)+(Lambda*D)))%*%t(B)%*%y,sigma=s2E*solve((t(B)%*%B)+(Lambda*D)))
#gsample <- B%*%t(Betasample)
#penalty <- (Betasample)%*%D%*%t(Betasample)
#out <- list(yfd=gsample,Pen=penalty,betahat=Betasample)
#out <- list(Pen=penalty,betahat=Betasample)
#out
return(Betasample)
}


sigma2EPosterior <- function(g,x,y){ ## x and y must be the whole vectors
        A <- A_E + n/2
        B <- ( (1/B_E) + 0.5*sum((y-eval.fd(x,fit))^2) )^(-1)
        #sigma2E <- 1/rgamma(1,A,B)
        sigma2E <- 1/rgamma(1,shape=A,scale=B)
        return(sigma2E)}
        
sigma2EPosterior2 <- function(g,x,y){ ## x and y must be the whole vectors
        A <- A_E + n/2
        B <- ( (1/B_E) + 0.5*sum((y-g)^2) )^(-1)
        #sigma2E <- 1/rgamma(1,A,B)
        sigma2E <- 1/rgamma(1,shape=A,scale=B)
        return(sigma2E)}
        
sigma2EPosterior3 <- function(BETA,x,y){ ## x and y must be the whole vectors
        gX <- getBMatrix(x)%*%BETA
        A <- A_E + n/2
        B <- ( (1/B_E) + 0.5*sum((y-gX)^2) )^(-1)
        #sigma2E <- 1/rgamma(1,A,B)
        sigma2E <- 1/rgamma(1,shape=A,scale=B)
        #sigma2E <- rinvgamma(1,shape=A,scale=B)
        return(sigma2E)}

#sumOverW <- function(w,x){  #W is 2 by 1 vector and x is only a scalar
#            f <- sum(w-x)^2
#            return(f)}

sigma2UPosterior <- function(W,x){
        A <- A_U + (1/2)*n*m
        B <- ( (1/B_U) + 0.5*sum((W[1,]-x)^2 + (W[2,]-x)^2))^(-1)
        sigma2U <- 1/rgamma(1,shape=A,scale=B)
        #sigma2U <- rinvgamma(1,shape=A,scale=B)
        return(sigma2U)}

gammaPosterior <- function(BETA,M){
              Penalty <- t(BETA)%*%PenaltyMatrixD%*%BETA
              A <- Agamma + M/2
              B <- ( (1/Bgamma) + 0.5*Penalty )^(-1)
              #Gam <- rgamma(1,A,B)
              Gam <- rgamma(1,shape=A,scale=B)
              return(Gam)}

muxPosterior <- function(x){
              Mean <- (n*mean(X)*t2_x + d_x*sigma2_x)/(n*t2_x + sigma2_x)
              Var <- (sigma2_x*t2_x)/(n*t2_x + sigma2_x)
              mu <- rnorm(1,mean=Mean,sd=sqrt(Var))
              return(mu)
              }

sigma2xPosterior <- function(x){
                  A <- A_x+ (n/2)
                  B <- (1/B_x + 1/2*sum((x-mu_x)^2) )^(-1)
                  sigma2x <- 1/rgamma(1,shape=A,scale=B)
                  #sigma2x <- rinvgamma(1,shape=A,scale=B)
                  return(sigma2x)
                  }
                  

X_MetropolisHastings3 <- function(W,Cx,y,s2E,s2U,s2x,mux,BETA,nInt){ ## nInt = number of interations sounds like 5000 is good

sd.prop = sqrt(s2U/2)
X <- Cx
XBig <- matrix(0,nrow=n,ncol=nInt)
#ac <- vector(mode='numeric',length=n)
ac <- vector(mode='numeric',length=n)
set.seed(1)

for(s in 1:nInt) {
gX <- getBMatrix(x=X)%*%BETA
X.p <- rnorm(n, X, sd=sd.prop)
gX.p <- getBMatrix(x=X.p)%*%BETA

#print(s)
#print(X.p)

for (i in 1:n) {
  lhr <- ((-0.5/s2U)*sum(((W[1,i]-X.p[i])^2 + (W[2,i]-X.p[i])^2))) - ((-0.5/s2U)*sum((W[1,i]-X[i])^2 + (W[2,i]-X[i])^2)) +
         ((-0.5/s2E)*sum((y[i]-gX.p[i])^2)) - ((-0.5/s2E)*sum((y[i]-gX[i])^2)) +
         ((-0.5/s2x)*sum((X.p[i]-mux)^2)) -   ((-0.5/s2x)*sum((X[i]-mux)^2))

      if( log(runif(1)) < lhr ) { X[i] <- X.p[i] ; ac[i] <- ac[i] + 1 }
      }
      #print(ac)

XBig[,s] <- X
                    }
posteriorMean <- apply(XBig[,-(1:500)],1,mean)
out <- list(postMean=posteriorMean,AC=(ac/nInt))
out
#return(posteriorMean)
                    }
                    
X_MetropolisHastings4 <- function(W,Cx,y,s2E,s2U,s2x,mux,BETA,nInt){ ## nInt = number of interations sounds like 5000 is good

sd.prop = sqrt(s2U/2)
X <- Cx  # initial values for X
mean.X <- Cx # fixed mean for my proposal distribution, the current value of X
XBig <- matrix(0,nrow=n,ncol=nInt)
#ac <- vector(mode='numeric',length=n)
ac <- vector(mode='numeric',length=n)
set.seed(1)

for(s in 1:nInt) {
gX <- getBMatrix(x=X)%*%BETA
X.p <- rnorm(n, mean.X, sd=sd.prop)
gX.p <- getBMatrix(x=X.p)%*%BETA

#print(s)
#print(X.p)

for (i in 1:n) {
  lhr <- ((-0.5/s2U)*sum(((W[1,i]-X.p[i])^2 + (W[2,i]-X.p[i])^2))) - ((-0.5/s2U)*sum((W[1,i]-X[i])^2 + (W[2,i]-X[i])^2)) +
         ((-0.5/s2E)*sum((y[i]-gX.p[i])^2)) - ((-0.5/s2E)*sum((y[i]-gX[i])^2)) +
         ((-0.5/s2x)*sum((X.p[i]-mux)^2)) -   ((-0.5/s2x)*sum((X[i]-mux)^2))  +
         ((-0.5/(sd.prop^2))*sum((X.p[i]-mean.X)^2)) - ((-0.5/(sd.prop^2))*sum((X[i]-mean.X)^2))

      if( log(runif(1)) < lhr ) { X[i] <- X.p[i] ; ac[i] <- ac[i] + 1 }
      }
      #print(ac)

XBig[,s] <- X                   
                    }
posteriorMean <- apply(XBig[,-(1:500)],1,mean)
out <- list(postMean=posteriorMean,AC=(ac/nInt))
out
#return(posteriorMean)
                    }
                    
X_MetropolisHastings5 <- function(W,Cx,y,s2E,s2U,s2x,mux,BETA,nInt){ ## nInt = number of interations sounds like 5000 is good

sd.prop = 2*sqrt(s2U/2)
X <- Cx
XBig <- matrix(0,nrow=nInt,ncol=1)
#ac <- vector(mode='numeric',length=n)
ac <- 0
set.seed(1)

for(s in 1:nInt) {
  #print(s)
gX <- getBMatrix(x=X)%*%BETA
X.p <- rnorm(1, X, sd=sd.prop)
gX.p <- getBMatrix(x=X.p)%*%BETA

  lhr <- ((-0.5/s2U)*sum(((W[1]-X.p)^2 + (W[2]-X.p)^2))) - ((-0.5/s2U)*sum((W[1]-X)^2 + (W[2]-X)^2)) +
         ((-0.5/s2E)*sum((y-gX.p)^2)) - ((-0.5/s2E)*sum((y-gX)^2)) +
         ((-0.5/s2x)*sum((X.p-mux)^2)) -   ((-0.5/s2x)*sum((X-mux)^2))

      if( log(runif(1)) < lhr ) { X <- X.p }
      #; ac[i] <- ac[i] + 1 }
      #print(ac)

XBig[s] <- X
                    }
posteriorMean <- mean(XBig[-(1:500)])
return(postMean=posteriorMean)
#return(posteriorMean)
                    }