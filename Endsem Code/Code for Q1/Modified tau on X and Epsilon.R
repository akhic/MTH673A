#Obtaining the approximate test statistics for N = 100
a <- function(x1,x2,x3,x4){
  return (sign(abs(x1-x2)+abs(x3-x4)-abs(x1-x3)-abs(x2-x4)))
}
taus <- vector()

for(t in 1:100){
  N <- 100
  x <- runif(N,-3,3)
  epsilon <- rnorm(N)
  y <- exp(-x)+epsilon
  
  
  tau <- 0 
  n <- 0
  
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
          tau <- tau + a(x[i],x[j],x[k],x[l])*
            a(epsilon[i],epsilon[j],epsilon[k],epsilon[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  taus <- c(taus,tau)
}
quantile(taus,probs = c(0.05, 0.95,0.025,0.975))
# The critical values we get from here are 
#-0.008740910 and 0.016835389 respectively

hist(taus,freq=FALSE,main=expression(paste("Approximated distribution of "
                                           , tau,"*", ", for N =100, under H0 "))
     ,xlab=expression(paste(tau,"*")))

library(ddalpha)
library(LaplacesDemon)


#Power Calculation
mu <- c(0,0)
Sigma <- matrix(c(1, 2/3, 2/3, 1), 2)
taua <- vector()
for(t in 1:100){
  bivn <- mvrnorm(N, mu = mu, Sigma = Sigma )
  x <- bivn[,1]
  epsilon <- bivn[,2]
  y <- exp(-x)+epsilon
  
  tau <- 0 
  n <- 0
  
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
          tau <- tau + a(x[i],x[j],x[k],x[l])*
            a(epsilon[i],epsilon[j],epsilon[k],epsilon[l])
          
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

taua
hist(taua,freq=FALSE,main=expression(paste("Approximated distribution of ", 
                                           tau,"*", ", for N =100, under HA ")),
     xlab=expression(paste(tau,"*")))

Sigma <- matrix(c(1, 0, 0, 1), 2)
taub <- vector()

for(t in 1:100){
  data_cauchy <- rmvc(n=N,mu,Sigma)
  x <- data_cauchy[,1]
  epsilon <- data_cauchy[,2]
  y <- exp(-x)+epsilon
  
  tau <- 0 
  n <- 0
  
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
          tau <- tau + a(x[i],x[j],x[k],x[l])*
            a(epsilon[i],epsilon[j],epsilon[k],epsilon[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  taub <- c(taub,tau)
  
}
taub
a1 <- quantile(taus,prob=0.025)
b1 <- quantile(taus,prob=0.975)
l1 <- taub >a1
l2 <- taub <b1
l <- l1 & l2
power <- sum(l)/N
power
sum(l)
a1
hist(taub,freq=FALSE,main=expression
     (paste("Approximated distribution of ", tau,"*", ", for N =100, under HB"))
     ,xlab=expression(paste(tau,"*")))