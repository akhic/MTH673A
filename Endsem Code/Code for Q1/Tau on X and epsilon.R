#Obtaining the approximate test statistics for N = 100
library(VGAM)
taus <- vector()

for(t in 1:1000){
  N <- 1000
  x <- runif(N,-3,3)
  epsilon <- rnorm(N)
  y <- exp(-x)+epsilon
  
  tau <- kendall.tau(x, epsilon, exact = FALSE, max.n = 3000)
  #tau <- tau/choose(N,4)
  taus <- c(taus,tau)
}
hist(taus)
quantile(taus,probs = c(0.05, 0.95,0.025,0.975))
# The critical values we get from here are -0.008740910 and 0.016835389 respectively
#hist(taus,freq=FALSE)
hist(taus,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,"", ", for N =1000, under H0 ")),xlab=expression(paste(tau,"")))

library(ddalpha)
library(LaplacesDemon)
library(mvtnorm)
library(MASS)
#Power Calculation
mu <- c(0,0)
Sigma <- matrix(c(1, 2/3, 2/3, 1), 2)
taua <- vector()
for(t in 1:1000){
  bivn <- mvrnorm(N, mu = mu, Sigma = Sigma )
 
  x <- bivn[,1]
  epsilon <- bivn[,2]
  y <- exp(-x)+epsilon
  
  tau <- kendall.tau(x, epsilon, exact = TRUE, max.n = 3000)
  
  taua <- c(taua,tau)
  
}

a1 <- quantile(taus,prob=0.025)
b1 <- quantile(taus,prob=0.975)
l1 <- taua >a1
l2 <- taua <b1
l <- l1 & l2
power <- 1-sum(l)/N
power

hist(taua,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,"", ", for N =1000, under HA ")),xlab=expression(paste(tau,"")))

Sigma <- matrix(c(1, 0, 0, 1), 2)
taub <- vector()

for(t in 1:1000){
  data_cauchy <- rmvc(n=N,mu,Sigma)
  x <- data_cauchy[,1]
  epsilon <- data_cauchy[,2]
  y <- exp(-x)+epsilon
  
  
  tau <- kendall.tau(x, epsilon, exact = TRUE, max.n = 3000)
  
  
  
  taub <- c(taub,tau)
  
}
a1 <- quantile(taus,prob=0.025)
b1 <- quantile(taus,prob=0.975)
l1 <- taub >a1
l2 <- taub <b1
l <- l1 & l2
power <- 1-sum(l)/N
power
hist(taub,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,"", ", for N =1000, under HB")),xlab=expression(paste(tau,"")))
