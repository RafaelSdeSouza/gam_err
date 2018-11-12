library("splines")
library("rstan")
X <- seq(from=-5, to=5, by=.1) # generating inputs
B <- t(bs(X, knots=seq(-5,5,1), degree=3, intercept = TRUE)) # creating the B-splines
N <- length(X)
num_basis <- nrow(B)
a0 <- 0.2 # intercept
a <- rnorm(num_basis, 0, 1) # coefficients of B-splines
Y_true <- as.vector(a0*X + a%*%B) # generating the output
Y <- Y_true + rnorm(length(X),0,.2) # adding noise

# Grid for prediction values

XX <- seq(from=-5, to=5, by=.05) # generating inputs
BB <- t(bs(XX, knots=seq(-5,5,1), degree=3, intercept = TRUE)) # creating the B-splines
M <- length(XX)
num_basis2 <- nrow(BB)

gdat <- data.frame(y=Y,x=X)

stan_data  <- list(Y = Y,
                   X = X,
                   B = B,
                   XX = XX,
                   BB = BB,
                   N = N,
                   M = M,
                   num_basis = num_basis,
                   num_basis2 = num_basis2)



stan_model= "data {
  int N;
  int M;
  int num_basis;
  int num_basis2;
  vector[N] Y;
  vector[N] X;
  vector[M] XX;
  matrix[num_basis, N] B;
  matrix[num_basis2, M] BB;
}

parameters {
row_vector[num_basis] a_raw;
real a0;
real<lower=0> sigma;
real<lower=0> tau;

}

transformed parameters {
row_vector[num_basis] a;
vector[N] Y_hat;
a = a_raw*tau; 
Y_hat = a0*X + to_vector(a*B);
}

model {
a_raw ~ normal(0, 1);
tau ~ cauchy(0, 1);
sigma ~ cauchy(0, 1);
Y ~ normal(Y_hat, sigma);
}

generated quantities{
vector[M] mu_pred;
// Posterior parameter distribution of the mean

mu_pred = a0*XX + to_vector(a*BB);

}


"

fit <- stan(model_code = stan_model,
              data = stan_data,
              seed = 42,
              chains = 3,
              iter =1500,
              cores= 3,
              warmup=750)

Y_mean <- extract(fit, "mu_pred")
Y_mean_cred <- apply(Y_mean$mu_pred, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$mu_pred, 2, mean)


fitdat <- data.frame(x = XX, y = Y_mean_mean, lwr1 = Y_mean_cred[1,],upr1 = Y_mean_cred[2,])

ggplot(data=gdat,aes(x=x,y=y)) +
  geom_point() +
  geom_line(data=fitdat,aes(x=x,y=y)) +
  geom_ribbon(data=fitdat,aes(x=x,ymin=lwr1, ymax=upr1,y=NULL),fill=c("gray50"),alpha=0.5) +
  theme_bw()

