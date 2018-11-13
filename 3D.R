library(fda)
library(mgcv)
library(Matrix)
require(visreg)
require("plot3D")

#####################
## Generating data ##
#####################

set.seed(30)
n <- 5000
m <- 1 ## no replicates 
k <- 30

GRID <- 200

#=== true X's 
X1true <- sort(runif(n,-5,5))
X2true <- sort(runif(n,-8,8))

#=== true functions
m1function <- function(x){
  f <- sin(pi*(x/2))/(1+(2*(x^2)*(sign(x)+1)))
  return(f)}


m2function <- function(x){
  f <- (1/2)*exp(-(1/2)*((x-3)^2))
  return(f)}

par(mfrow=c(2,1))
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


data_test <- as.data.frame(cbind(Y, X1true, X2true))
b <- gam(Y~s(X1true,bs="cr",k=15)+s(X2true,bs="cr",k=15),data=data_test,method="REML")
summary(b)
par(mfrow=c(2,2))
plot(X1true,m1function(X1true),type="l")
plot(X2true,m2function(X2true),type="l")
visreg(b)


b_rec <- function(obj,N){
  b <- obj$coefficients[1:N]
  Q <-  model.matrix(obj)[,1:N]
  rec <-  rowSums(t(b[1:N]*t(Q[,1:N])))
  return(rec)
  }


plot(X1true,b_rec(b,5))

#plot(b,all.terms=TRUE,scheme=c(2,1))



# predict values on regular xy grid
grid.lines = 50
x1.pred <- seq(1.01*min(X1true), 0.99*max(X1true), length.out = grid.lines)
x2.pred <- seq(1.01*min(X2true), 0.99*max(X2true), length.out = grid.lines)
x1x2 <- expand.grid( X1true = x1.pred, X2true = x2.pred)

y.pred <- matrix(predict(b, newdata = x1x2,type="response"), 
                 nrow = grid.lines, ncol = grid.lines)


y.true <- matrix(m1function(x1x2$X1true) + m2function(x1x2$X2true),
                 nrow = grid.lines, ncol = grid.lines)


# fitted points for droplines to surface
fitpoints <- predict(b,type="response")

par(mfrow=c(1,2))
persp3D(x=x1.pred, y=x2.pred, z = y.true ,   cex = 1, cex.lab=1.5,type="p",pch = 19,alpha=0.5,
        theta = 30, phi = 30, ticktype = "detailed",col="red",bty = "b2",
        xlab="x1",
        ylab="x2",
        zlab="y")

persp3D(x=x1.pred, y=x2.pred, z= y.pred ,   cex = 1, cex.lab=1.5,type="p",pch = 19,alpha=0.5,
          theta = 30, phi = 30, ticktype = "detailed",col="red",bty = "b2",
          xlab="x1",
          ylab="x2",
          zlab="y")


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
plot(Bbasis_obj_1)
plot(Bbasis_obj_2)

#### end of creating basis objects

### evaluating basis functions at x 

eval.Bbasis <- function(x,Bbasis_obj){
  
  B <- getbasismatrix(x, Bbasis_obj, nderiv=0)
  return(B)
  
}

B1 <- eval.Bbasis(X1true,Bbasis_obj_1)
B2 <- eval.Bbasis(X2true,Bbasis_obj_2)

### Penalty matrices, fixed across iterations
D1 <- getbasispenalty(Bbasis_obj_1, Lfdobj=int2Lfd(2))
D2 <- getbasispenalty(Bbasis_obj_2, Lfdobj=int2Lfd(2))

Penalty_array <- array( c( D1 , D2) , dim = c( dim(D1)[1] , dim(D1)[2] , 2 ) )

DB_matrix <- cbind(B1,B2)

alpha <- c(0.2,0.8) 

beta_estimation <- function(y,DB_matrix, penalty_matrices, alpha){
  
  diag_term <- bdiag(alpha[1]*penalty_matrices[,,1],alpha[2]*penalty_matrices[,,2]) 
  
  beta <- solve(((t(DB_matrix)%*%DB_matrix)+(diag_term)))%*%t(DB_matrix)%*%y
 return(beta) 
}

beta <- beta_estimation(Y,DB_matrix, Penalty_array, alpha) 

fit2 <- B2%*%beta[33:64,1]
plot(X2true,m2function(X2true))
lines(X2true,fit2)


fit1 <- B2%*%beta[1:32,1]
plot(X1true,m1function(X1true))
lines(X1true,fit1)