#Obtaining the approximate test statistics for N = 100
a <- function(x1,x2,x3,x4){
  return (sign(abs(x1-x2)+abs(x3-x4)-abs(x1-x3)-abs(x2-x4)))
}


library(VGAM)
taus <- vector()
for(t in 1:100){
  N <- 100
  x <- runif(N,-3,3)
  x <- sort(x)
  epsilon <- rnorm(N)
  y <- exp(-x)+epsilon
  ys <- rep(0,N)
  ys[1] <- y[2]-y[1]
  ys[N] <- ys[N-1] - y[N]
  for(j in 2:(N-1)){
    ys[j] <- y[j-1] - 2*y[j]+y[j+1]
  }
  
  tau <- 0 
  for(i in 1:N){
    if(i>N-3 || i<1){
      break
    }
    for(j in i+1:N){
      if(j>N-2 || j<1){
        break
      }
      for(k in j+1:N){
        if(k>N-1 || k<1){
          break
        }
        for(l in k+1:N){
          if(l>N || l<1){
            break
          }
          n <- n +1
          tau <- tau + a(x[i],x[j],x[k],x[l])*a(ys[i],ys[j],ys[k],ys[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  
  
  taus <- c(taus,tau)
}
quantile(taus,probs = c(0.05, 0.95,0.025,0.975))

hist(taus,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,
                                           "*", ", for N =100, under H0 ")),
     xlab=expression(paste(tau,"*")))

library(ddalpha)
library(LaplacesDemon)
library(mvtnorm)
library(MASS)
#Power Calculation
mu <- c(0,0)
Sigma <- matrix(c(1, 2/3, 2/3, 1), 2)
taua <- vector()
for(t in 1:100){
  bivn <- mvrnorm(N, mu = mu, Sigma = Sigma )
  bivn <- bivn[order(bivn[,1]),]
  x <- bivn[,1]
  epsilon <- bivn[,2]
  y <- exp(-x)+epsilon
  ys <- rep(0,N)
  ys[1] <- y[2]-y[1]
  ys[N] <- ys[N-1] - y[N]
  for(j in 2:(N-1)){
    ys[j] <- y[j-1] - 2*y[j]+y[j+1]
  }
  tau <- 0 
  for(i in 1:N){
    if(i>N-3 || i<1){
      break
    }
    for(j in i+1:N){
      if(j>N-2 || j<1){
        break
      }
      for(k in j+1:N){
        if(k>N-1 || k<1){
          break
        }
        for(l in k+1:N){
          if(l>N || l<1){
            break
          }
          n <- n +1
          tau <- tau + a(x[i],x[j],x[k],x[l])*a(ys[i],ys[j],ys[k],ys[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  taua <- c(taua,tau)
  
}

a1 <- quantile(taus,prob=0.025)
b1 <- quantile(taus,prob=0.975)
l1 <- taua >a1
l2 <- taua <b1
l <- l1 & l2
power <- 1-sum(l)/N
power

hist(taua,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,
                                           "*", ", for N =100, under HA "))
     ,xlab=expression(paste(tau,"*")))

Sigma <- matrix(c(1, 0, 0, 1), 2)
taub <- vector()

for(t in 1:100){
  data_cauchy <- rmvc(n=N,mu,Sigma)
  data_cauchy <- data_cauchy[order(data_cauchy[,1]),]
  x <- data_cauchy[,1]
  epsilon <- data_cauchy[,2]
  y <- exp(-x)+epsilon
  ys <- rep(0,N)
  ys[1] <- y[2]-y[1]
  ys[N] <- ys[N-1] - y[N]
  for(j in 2:(N-1)){
    ys[j] <- y[j-1] - 2*y[j]+y[j+1]
  }
  
  tau <- 0 
  for(i in 1:N){
    if(i>N-3 || i<1){
      break
    }
    for(j in i+1:N){
      if(j>N-2 || j<1){
        break
      }
      for(k in j+1:N){
        if(k>N-1 || k<1){
          break
        }
        for(l in k+1:N){
          if(l>N || l<1){
            break
          }
          tau <- tau + a(x[i],x[j],x[k],x[l])*a(ys[i],ys[j],ys[k],ys[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  
  
  
  taub <- c(taub,tau)
  
}
a1 <- quantile(taus,prob=0.025)
b1 <- quantile(taus,prob=0.975)
taub <- taub[!is.na(taub)]
l1 <- taub >a1
l2 <- taub <b1
l <- l1 & l2
power <- 1-sum(l)/length(taub)
power
hist(taub,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,"*"
                                           , ", for N =100, under HB")),
     xlab=expression(paste(tau,"*")))

