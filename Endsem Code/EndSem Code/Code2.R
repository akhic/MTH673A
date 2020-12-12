library(MLmetrics)
library(quantreg )
library(robustbase)
library(ggplot2)
library(MASS)
library(gamlss)
N <- 10000

Mean_statistics <- vector()
Median_statistics <- vector()
LTS_statistics <- vector()
LMS_statistics <- vector()
Mode_statistics <- vector()
x <- rep(0,N)

for(t in 1:1000){
  y <- rnorm(N)
  Mean_statistics <- c(Mean_statistics,mean(y))
  
  Median_statistics <- c(Median_statistics,median(y))
  
  # Calculating mode using a gaussian KDE  
  density_estimate <- density(y)
  
  mode_value <- density_estimate$x[which.max(density_estimate$y)]
  Mode_statistics <- c(Mode_statistics,mode_value)
  
  lts <- ltsReg(y = y,intercept = TRUE)  
  LTS_statistics <- c(LTS_statistics,lts$coefficients)
  
  
  ysort <- sort(y)
  len <- 1000000
  lms <- 9999999
  currlen <- 0
  avg <- 9999999
  # Calculating median using shortest half(Algorithm provided in class)
  for(i in 1:(N/2)){
    
    currlen <- ysort[N/2+i]-ysort[i]
    avg <- (ysort[i+N/2]+ysort[i])/2
    if(currlen<len){
      len <- currlen
      lms <- avg
    }
  }  
  LMS_statistics <- c(LMS_statistics,lms)
  
}


hist(LMS_statistics)
hist(LTS_statistics)
hist(Mean_statistics)
hist(Median_statistics)
hist(Mode_statistics)

var(LMS_statistics)
var(LTS_statistics)
var(Mean_statistics)
var(Median_statistics)
var(Mode_statistics)



