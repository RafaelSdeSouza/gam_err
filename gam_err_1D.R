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