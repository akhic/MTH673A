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
          tau <- tau + a(x[i],x[j],x[k],x[l])*a(epsilon[i],epsilon[j],epsilon[k],epsilon[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  taus <- c(taus,tau)
}
quantile(taus,probs = c(0.05, 0.95,0.025,0.975))
# The critical values we get from here are -0.008740910 and 0.016835389 respectively
#hist(taus,freq=FALSE)
hist(taus,freq=FALSE,main=expression(paste("Approximated distribution of ", tau,"", ", for N =100 ")),xlab=expression(paste(tau,"")))

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
          tau <- tau + a(x[i],x[j],x[k],x[l])*a(epsilon[i],epsilon[j],epsilon[k],epsilon[l])
          
        }
      }
    }
  }
  tau <- tau/choose(N,4)
  taua <- c(taua,tau)
  
}
  