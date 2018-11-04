library(mvtnorm)
library(fda)
#source('functions.r')
source('functions3.r')

############################
## Generating data      ####
############################

set.seed(30)
n <- 50
Xtrue <- sort(rnorm(n,0,1))
GRID <- 200
XICM <- scan("XICM.txt") ### X's obtained from ICM algorithm
k <- 30

Y <- mfunction(Xtrue) + rnorm(n,mean=0,sd=.3)

m <- 2
Wmatrix <- sapply(Xtrue,rnorm,n=m,sd=0.8) #each column is for one i

###########################
## prior hyperparameters ##
###########################

A_E <- 1
B_E <- 1
#plot(density(rgamma(1000, shape=2, scale = 1)),xlim=c(0,200))
#lines(density(1/rgamma(1000,1,scale=1)),col=2)
A_U <- 1
B_U <- 1
Agamma <- 3
Bgamma <- 1000
d_x <-  0
t2_x <- 1
A_x <-  1
B_x <-  1

#####################
## starting values ##
#####################

sigma2_U <- mean(apply(Wmatrix,2,var))  #pooled variance
#sigma2_U <- 0.64

Wbar <- colMeans(Wmatrix)
Quantiles <- round(seq(from=0,to=1,length=k+2),4)[-c(1,k+2)]
KNOTS <- as.vector(quantile(Wbar,probs=Quantiles))
rng <- c(-18,18)
PenaltyMatrixD <- getbasispenalty(create.bspline.basis(rangeval=rng, norder=4, breaks=KNOTS), Lfdobj=int2Lfd(2))
ybasis <- create.bspline.basis(rangeval=rng, norder=4, breaks=KNOTS)

#BWBar <- getBMatrix(x=Wbar)

mu_x <- mean(c(Wmatrix[1,],Wmatrix[2,]))
sigma2_x <- var(c(Wmatrix[1,],Wmatrix[2,]))

#mu_x <- 0
#sigma2_x <- 1

#Xgrid <- seq(from=min(X),to=max(X),length=GRID)
#X <- Wbar 
#CV <-  cvFunction(x=X,y=Y)
#lambda <- CV$lambda
#fit <- smooth.Camila(X,Y,lambda)
#eval.fd(6,fit)
#ghat <- eval.fd(X,fit)
#ghat2 <- eval.fd(Xgrid, fit)
#DF <- CV$DF
#sigma2_E <- sum((Y-ghat)^2)/(n-DF)

sigma2_E <- 0.21
#sigma2_E <- 0.1

lambda <- df2lambda(Wbar, ybasis, Lfdobj=int2Lfd(2), df=4.606075)

Gamma <- lambda/sigma2_E ### lambda here is alpha in the paper

X <- XICM  ## that's my initial value for X the ones I got from ICM

S <- 10

XBIG <- matrix(0,nrow=n,ncol=S)
LAMBDA <- vector(mode='numeric',S)
#FITTED.VALUES <- matrix(0,nrow=n,ncol=S)
BETAHAT <- matrix(0,nrow=32,ncol=S)
SIGMA2E <- vector(mode='numeric',S)

#ptm <- proc.time()

for(i in 1:S){

print(i)

betahat <- t(gfunction2(x=X,y=Y,Gam=Gamma,s2E=sigma2_E))

BETAHAT[,i] <- betahat
LAMBDA[i] <- Gamma*sigma2_E
XBIG[,i] <-  X
SIGMA2E[i] <- sigma2_E

Xnew <- vector(mode='numeric',length=50)
for (j in 1:n){
  print(j)
Xnew[j] <- X_MetropolisHastings5(W=Wmatrix[,j],Cx=X[j],y=Y[j],s2E=sigma2_E,s2U=sigma2_U,s2x=sigma2_x,mux=mu_x,BETA=betahat,nInt=1000)
}

#Xnew <- X_MetropolisHastings4(W=Wmatrix,xicm=XICM,Cx=X,y=Y,s2E=sigma2_E,s2U=sigma2_U,s2x=sigma2_x,mux=mu_x,BETA=betahat,nInt=2000)
#print(Xnew$AC)
#Xnew <- Xnew$postMean


#Xnew <- vector(mode='numeric',length=n)
#for(k in 1:n){
#    print(k)
#    Xnew[k] <- X_MetropolisHastings(W=Wmatrix[,k],y=Y[k],s2E=sigma2_E,s2U=sigma2_U,s2x=sigma2_x,mux=mu_x,g=fit,nInt=2000)
#    }
    
X <- Xnew

sigma2_E <- sigma2EPosterior3(BETA=betahat,x=X,y=Y)
print(sigma2_E)
sigma2_U <- sigma2UPosterior(W=Wmatrix,x=X)
print(sigma2_U)
#Gamma <- gammaPosterior(Penalty=Smooth$Pen,M=30)
Gamma <- gammaPosterior(BETA=betahat,M=30)
print(Gamma)

mu_x <- muxPosterior(x=X)

sigma2_x <- sigma2xPosterior(x=X)
print(sigma2_x)

}

#proc.time() - ptm

#write.table(BETAHAT,'betahat.txt',row.names=FALSE,col.names=FALSE)
#write.table(XBIG,'xbig2.txt',row.names=FALSE,col.names=FALSE)
#write(LAMBDA,'lambda2.txt')
#write(SIGMA2E,'sigma2E.txt')

Burnin <- 20
(mean.sigma2E <- mean(SIGMA2E[-(1:Burnin)]))
(CIsigma2E <- quantile(SIGMA2E[-(1:Burnin)],probs=c(0.025,0.975)))
mean.beta <- apply(BETAHAT[,-(1:Burnin)],1,mean)

#Xgrid <- seq(from=min(Xtrue),to=max(Xtrue),length=GRID)
XGRID <- XBIG[,-(1:Burnin)]
BETA.hat <- BETAHAT[,-(1:Burnin)]
#Bgrid <- getBMatrix(Xgrid)

getfitted <- function(B,b){
            g <- B%*%b
            return(g)}
fittedvalues <- matrix(0,nrow=200,ncol=(S-Burnin))

for(i in 1:(S-Burnin)){
      fittedvalues[,i] <- getBMatrix(Xgrid)%*%BETA.hat[,i]
      }

XGRID2 <- apply(XBIG[,-(1:Burnin)],1,mean)
mean.fitted <- apply(fittedvalues,1,mean)
CI.fitted <- apply(fittedvalues,1,quantile,probs=c(0.025,0.975))
#dim(fittedvalues)
#mean.fitted <-rowMeans(fittedvalues)

#mean.fitted <- apply(FITTED.VALUES[,-(1:Burnin)],1,mean)
(lambda.mean <- mean(LAMBDA[-(1:Burnin)]))
#CI <- apply(FITTED.VALUES[,-(1:Burnin)],1,quantile,probs=c(0.025,0.975))
teste <- smooth.Camila(x=XGRID2,y=Y,Lambda=lambda.mean)
teste.fit <- eval.fd(XGRID2,teste)

pdf('feio.pdf')
plot(Xtrue,Y,xlim=c(-3,3),xlab='x',ylab='m(x)')
#lines(Xgrid,mfunction(Xgrid),col=1,lwd=1)
lines(sort(XGRID2),teste.fit[order(XGRID2)],col=3)
lines(Xgrid,mean.fitted,lty=1,col=2,lwd=1)
lines(Xgrid,CI.fitted[1,],col=4)
#lines(sort(XGRID2),teste2[order(XGRID2)],col=3,lwd=1)
#lines(sort(X),ghat[order(X)],col=4,lwd=1)
#lines(order(XGRID2),mean.fitted[order(XGRID2)],col=2,,lty=1,lwd=1)
dev.off()

lines(X,CI[1,],lty=2,col=2,lwd=1)
lines(X,CI[2,],lty=2,col=2,lwd=1)
#dev.off()
#legend('bottomright',c('True Curve','Classical Spline','Bayesian Spline','95% Credible Intervals'),
#        col=c(1,4,2,2),lty=c(1,1,1,2))



#write.table(FITTED.VALUES,'fitted.values.txt',row.names=FALSE,col.names=FALSE)
#write.table(XBIG,'xbig.txt',row.names=FALSE,col.names=FALSE)
#write(LAMBDA,'lambda.txt')









#R <- eval.penalty(ybasis,Lfdobj=int2Lfd(2),rng= c(min(X),max(X)))
#dim(R) ### that's the matrix D using Berry's notation
#
#B <- getbasismatrix(X, ybasis, nderiv=0)
#D <- R
#alpha <- lambda
#
#ghat <- B%*%solve(t(B)%*%B+alpha*D)%*%t(B)%*%Y
#lines(X,ghat,col=3)
#
#
#matrixA <- B%*%solve(t(B)%*%B+alpha*D)%*%t(B)
#sum(diag(matrixA)) #trace of A(alpha)
#lambda2df(X, ybasis,Lfdobj=int2Lfd(2),lambda=lambdagrid[which(cv==min(cv))]) ## this agrees with trace of A(alpha) :-)
#
#mhat$df ##awesome it also agrees :-)




