library(MASS)
library(ddalpha)
library(LaplacesDemon)
library(matrixcalc)
Sigma <- matrix(c(1,0,0,0,1,0,0,0,1),3,3,3)
Sigma
mu = c(0,0,0)
N_sample = 100
data1 <- mvrnorm(n=N_sample,mu,Sigma)
data2 <- mvrnorm(n=N_sample,mu,Sigma)
#depth.space.halfspace(data1, c(10000, 10000,10000), exact = TRUE)
#x=c(0,0,0)
L = 100
x <- seq(-5,5,length.out =L)
y <- seq(-5,5,length.out =L)
z <- seq(-5,5,length.out =L)
pt <- matrix(,nrow=L*L*L,ncol=3)
pts = L*L*L
for (i in 1:L) {
  for(j in 1:L){
    for(k in 1:L){
      pt[(i-1)*L^2+(j-1)*L+k,1] <- x[i]
      pt[(i-1)*L^2+(j-1)*L+k,2] <- y[j]      
      pt[(i-1)*L^2+(j-1)*L+k,3] <- z[k]
    }
  }
}
pt

#I will try 2 methods
# 1) Getting depth for both on the generate points
# 2) Getting depth for both on randomly generated points
dep1 = depth.halfspace(data1,data1)
plot(dep1)
dep2 = depth.halfspace(data1,data2)
plot(dep1,dep2)
dep1 = c(dep1,depth.halfspace(data2,data1))
dep2 = c(dep2,depth.halfspace(data2,data2))
plot(dep1,dep2,main="DD plot for 2 standard Gaussians",xlab='Depth 1',ylab='Depth 2')
lines(abline(0,1,col="red",lwd=3))
data3 <- rmvc(n=N_sample,mu)
data3
dep31 = depth.halfspace(data1,data3)
dep31 = c(dep31,depth.halfspace(data3,data3))
dep13 = depth.halfspace(data1,data1)
dep13 = c(dep13,depth.halfspace(data3,data1))
plot(dep13,dep31,main="DD plot for standard Gaussian with Standard Cauchy",xlab='Depth Normal',ylab='Depth Cauchy')
lines(abline(0,1,col="red",lwd=3))

mu1 = c(1,1,1)
data4 <-mvrnorm(N_sample,mu1,Sigma)
data5 <- rbind(data1,data4)
data5
dep14 = depth.halfspace(data5,data1)
depth41 = depth.halfspace(data5,data4)
plot(dep14,depth41,main="DD plot for standard Gaussian with shifted gaussian",xlab='Depth Normal',ylab='Depth Shifted Normal')
lines(abline(0,1,col="red",lwd=3))

Sigma1 <- matrix(c(3,0,0,0,3,0,0,0,3),3,3,3)
data6<-mvrnorm(N_sample,mu,Sigma1)
data7 <- rbind(data1,data6)
dep15 = depth.halfspace(data7,data1)
dep51 = depth.halfspace(data7,data6)
plot(dep15,dep51,,main="DD plot for standard Gaussian with non standard gaussian",xlab='Depth Normal',ylab='Depth Flattened Normal')
lines(abline(0,1,col="red",lwd=3))

Sigma2 <- matrix(c(1,0,1,0,1,0,1,0,1),3,3,3)
data8 <-mvrnorm(N_sample,mu,Sigma2)
data9 <- rbind(data1,data8)
dep16 = depth.halfspace(data9,data1)
dep61 = depth.halfspace(data9,data8)
plot(dep16,dep61,,main="DD plot for standard Gaussian with non standard gaussian(Correlation)",xlab='Depth Normal',ylab='Depth (Correlated) Normal')
lines(abline(0,1,col="red",lwd=3))
is.positive.semi.definite(Sigma2, tol=1e-8)

Sigma3 <- matrix(c(3,0,3,0,3,0,3,0,3),3,3,3)
data10 <-mvrnorm(N_sample,mu1,Sigma3)
data11 <- rbind(data1,data10)
dep17 = depth.halfspace(data11,data1)
dep71 = depth.halfspace(data11,data10)
plot(dep17,dep71,,main="DD plot for standard Gaussian with non standard gaussian(Correlation,mu)",xlab='Depth Normal',ylab='Depth (Correlated,Shifted) Normal')
lines(abline(0,1,col="red",lwd=3))

mu2 = c(1.2,-0.9,0.65)
Sigma4 <- matrix(c(2.11,0.63,0.67,0.63,3.32,1.23,0.67,1.23,3),3,3,3)
is.positive.semi.definite(Sigma4, tol=1e-8)
data12 <-mvrnorm(N_sample,mu2,Sigma4)
data13 <- rbind(data1,data12)
dep18 = depth.halfspace(data13,data1)
dep81 = depth.halfspace(data13,data12)
plot(dep18,dep81,,main="DD plot for standard Gaussian with non standard gaussian(Correlation,mu)",xlab='Depth Normal',ylab='Depth (Correlated,Shifted) Normal')
lines(abline(0,1,col="red",lwd=3))



