library(fda)
library(mgcv)
library(Matrix)


#####################
## Generating data ##
#####################

set.seed(30)
n <- 50
m <- 1 ## no replicates 
k <- 30

GRID <- 200

#=== true X's 
X1true <- sort(rnorm(n,0,1))
X2true <- sort(rnorm(n,3,1))

#=== true functions
m1function <- function(x){
            f <- sin(pi*(x/2))/(1+(2*(x^2)*(sign(x)+1)))
            return(f)}


m2function <- function(x){
  #f <- (1/2)*exp(-(1/2)*((x-3)^2))
  f <- (x-3)^2
  return(f)}

par(mfrow=c(1,2))
plot(X1true,m1function(X1true),type="l")
plot(X2true,m2function(X2true),type="l")

#=== generating Y as an additive model
Y <- m1function(X1true) + m2function(X2true) + rnorm(n,mean=0,sd=.01)

#=== only W1 and W2 are observed

sd_u1 <- round(rnorm(n,mean=0.8,sd=0.2),2)
sd_u2 <- round(rnorm(n,mean=1,sd=0.2),2)

W1 <- sapply(X1true,rnorm,n=m,sd=sd_u1) #no replicates, m=1
W2 <- sapply(X2true,rnorm,n=m,sd=sd_u2) #no replicates, m=1

#### End of generating data

#### Testing gam function versus B-splines 

test <- FALSE
if(test == TRUE){

data_test <- as.data.frame(cbind(Y, X1true, X2true))
b <- gam(Y~ s(X1true,bs="ps")+s(X2true,bs="ps"),data=data_test,method="GCV.Cp")
summary(b)
par(mfrow=c(2,2))
plot(X1true,m1function(X1true),type="l")
plot(X2true,m2function(X2true),type="l")
plot(b,all.terms=TRUE,scheme=c(2,1))

b$sp ### smoothing parameters

coef1 <- b$coefficients[2:10]
M1 <- model.matrix(b)[,2:10]
m1_hat <- M1%*%coef1
par(mfrow=c(2,2))
plot(X1true,m1function(X1true),type="l",ylim=c(-1,2))
lines(X1true, m1_hat,lty=2)

coef2 <- b$coefficients[11:19]
M2 <- model.matrix(b)[,11:19]
m2_hat <- M2%*%coef2
plot(X2true,m2function(X2true),type="l",ylim=c(-1,7))
lines(X2true, m2_hat,lty=2)

# data_test <- as.data.frame(cbind(Y, X1true, X2true))
# b <- gam(Y~s(X1true,bs="cr",k=15)+s(X2true,bs="cr",k=15),data=data_test,method="REML")
# summary(b)
# par(mfrow=c(2,2))
# plot(X1true,m1function(X1true),type="l")
# plot(X2true,m2function(X2true),type="l")
# plot(b,all.terms=TRUE,scheme=c(2,1))


### creating B-splines 

create.bsplines.f <- function(knots,range,order){

  Bbasis_obj <- create.bspline.basis(rangeval=range, norder=order, breaks=knots)   
  return(Bbasis_obj) 
}

create.knots <- function(x,k){
  
  Quantiles <- round(seq(from=0,to=1,length=k+2),4)[-c(1,k+2)]
  KNOTS <- as.vector(quantile(x,probs=Quantiles))
  return(KNOTS)
  
}


Bbasis_obj_1 <- create.bsplines.f(knots=create.knots(x=X1true,k=30),range=c(-5,5),order=4)
Bbasis_obj_2 <- create.bsplines.f(knots=create.knots(x=X2true,k=30),range=c(-0.5,6.6),order=4)
#plot(Bbasis_obj_1)
#plot(Bbasis_obj_2)

#### end of creating basis objects

### evaluating basis functions at x 
eval.Bbasis <- function(x,Bbasis_obj){
  
  B <- getbasismatrix(x, Bbasis_obj, nderiv=0)
  return(B)
  
}

B1 <- eval.Bbasis(X1true,Bbasis_obj_1)
B2 <- eval.Bbasis(X2true,Bbasis_obj_2)

### Penalty matrices
D1 <- getbasispenalty(Bbasis_obj_1, Lfdobj=int2Lfd(2))
D2 <- getbasispenalty(Bbasis_obj_2, Lfdobj=int2Lfd(2))


Penalty_array <- array( c( D1 , D2) , dim = c( dim(D1)[1] , dim(D1)[2] , 2 ) )

DB_matrix <- cbind(B1,B2)

#alpha <- c(0.02,0.08) 
alpha <- b$sp

beta_estimation <- function(y,DB_matrix, penalty_matrices, alpha){
 
  diag_term <- bdiag(alpha[1]*penalty_matrices[,,1],alpha[2]*penalty_matrices[,,2]) 

  beta <- solve(((t(DB_matrix)%*%DB_matrix)+(diag_term)))%*%t(DB_matrix)%*%y
  
}

beta_hat <- beta_estimation(y=Y,DB_matrix=DB_matrix,penalty_matrices=Penalty_array, alpha=alpha)

fit1 <- B1%*%beta_hat[1:32,1]
plot(X1true,m1function(X1true),type="l",ylim=c(-1,2))
lines(X1true,fit1,col=2)

fit2 <- B2%*%beta_hat[33:64,1]
plot(X2true,m2function(X2true),type="l",ylim=c(-1,7))
lines(X2true,fit2)

}

### End of testing 

#############################################################################################
## First Part of ICM algorithm: finding sigma2_U, sigma2_E and naive m function (or g^(0)) ##   
#############################################################################################

### estimate sigma2_U, sigma2_X and mu_X
sigma2_U1 <- sd_u1^2
sigma2_U2 <- sd_u2^2

## starting values for X and g 
X10 <- W1
X20 <- W2

sigma2_X1 <- 1
sigma2_X2 <- 1

mu_X1 <- 0
mu_X2 <- 3

### fitting the initial GAM as if X = W. ### stopped here Nov 5th

data_test <- as.data.frame(cbind(Y, X10, X20))
b0 <- gam(Y~ s(X10,bs="ps")+s(X20,bs="ps"),data=data_test,method="GCV.Cp")
summary(b0)
par(mfrow=c(2,2))
plot(X1true,m1function(X1true),type="l")
plot(X2true,m2function(X2true),type="l")
plot(b0,all.terms=TRUE,scheme=c(2,1))

b0$sp ### smoothing parameters
coef1 <- b0$coefficients[2:10]
M1 <- model.matrix(b0)[,2:10]
m1_hat <- M1%*%coef1
par(mfrow=c(2,2))
plot(X1true,m1function(X1true),type="l",ylim=c(-5,10))
lines(X10, m1_hat,type="p")

coef2 <- b0$coefficients[11:19]
M2 <- model.matrix(b0)[,11:19]
m2_hat <- M2%*%coef2
plot(X2true,m2function(X2true),type="l",ylim=c(-1,7))
lines(X20, m2_hat,type="p")

df = n - b0$df.residual ### effective degrees of freedom 

### 11.89 is the total df 
sigma2_E <- sum((b0$fitted.values-Y)^2)/(n - df)
## or 
#sigma2_E <- b0$sig2

sigm2_E <- 1e-04

#################################################################################################
## Second part of ICM: iterating between 1.) find X given g fixed and 2.) find g with X fixed. ##
#################################################################################################

ghat <- b0$fitted.values ### ghat = m1_hat + m2hat

b <- b0

X1Big <- NULL
X2Big <- NULL
iter <- 20

#X1p <- seq(from=-5,to=5,length=200)
#X2p <- seq(from=0,to=7,length=200)

posterior.f <- function(x, y, sum_hat, W, mu_X, s2_X, s2_U, s2_E){
  g <- sum_hat
  posterior <-  (-1/(2*sigma2_E))*((y-g)^2) + (-1/2)*sum((((W-x)^2)/s2_U)) + (-1/2)*sum((((x-mu_X)^2)/s2_X))
  return(posterior)
    }

### testing posterior function
#k=2
#posterior.f(x=Xinitial, y=Y[k], sum_hat=ghat[k], W=c(W1[k],W2[k]), mu_X=c(mu_X1,mu_X2), s2_X=c(sigma2_X1,sigma2_X2), 
#      s2_U=c(sigma2_U1[k],sigma2_U2[k]), s2_E=sigma2_E)
# end of testing


for(i in 1:iter){
  
### finding new X         
### Very simple grid search looking at each coordinate 
Xnew <- NULL
  for (k in 1:n){ 
    print(k)
    Xinitial <- c(X10[k],X20[k]) ### not sure which initial values to use
    X_opt <- optim(par=Xinitial,y=Y[k], sum_hat = ghat[k], W=c(W1[k],W2[k]), mu_X=c(mu_X1,mu_X2), s2_X=c(sigma2_X1,sigma2_X2), 
                    s2_U=c(sigma2_U1[k],sigma2_U2[k]), s2_E=sigma2_E, fn=posterior.f, method="L-BFGS-B",upper=c(5,7),lower=c(-5,0))$par    
    
    Xnew <- rbind(Xnew,X_opt)
    rownames(Xnew)=NULL
  }
  

X1Big <- cbind(X1Big,Xnew[,1])
X2Big <- cbind(X2Big,Xnew[,2])

data_test <- as.data.frame(cbind(Y[1:5], Xnew[,1], Xnew[,2]))
colnames(data_test) <- c("Y","X1new","X2new")
bnew <- gam(Y~ s(X1new,bs="ps")+s(X2new,bs="ps"),data=data_test,method="GCV.Cp")

#summary(b0)
if(i==20){
par(mfrow=c(2,2))
plot(X1true,m1function(X1true),type="l")
plot(X2true,m2function(X2true),type="l")
plot(bnew,all.terms=TRUE,scheme=c(2,1))
}

c1 <- mean(abs(bnew$fitted.values-ghat)) ### calculating delta for convergence 
print(c1)

ghat <- bnew$fitted.values

b <- bnew

X10 <- Xnew[,1]
X20 <- Xnew[,2]

}




