library(MLmetrics)
library(quantreg )
library(robustbase)


  
tau = seq(0.05,0.5,0.05)

filtered <- function(alpha,x,y)
{
  reg1 <- rq(y ~ x,alpha)
  reg2 <- rq(y ~ x,1-alpha)
  #print(summary(reg1))
  #print(summary(reg2))
  beta1 <- (reg1$coefficients[2])
  c1 <- reg1$coefficients[1]
  beta2 <- reg2$coefficients[2]
  c2 <- reg2$coefficients[1]
  #print(beta1)
  #print(beta2)
  
  boolx1 <- (y-beta1*x-c1)>0
  #print(sum(boolx1))
  boolx2 <- (y-beta2*x-c2)<0
  #print(sum(boolx2))
  boolx <- boolx1 & boolx2
  newx <- vector()
  newy <- vector()
  for(i in 1:length(x)){
    if(boolx[i]){
      newx <- c(newx,x[i])
      newy <- c(newy,y[i])
    }
  }
  
  xy <- cbind(newx,newy)
  return (xy)
  
}
MSE_TM <- vector()
MSE_LTS <- vector()
for(i in 1:100)
{
  epsilon <- rnorm(100,mean=0,sd=3)
  x <- runif(100,0,100)
  y <- 0.5*x+epsilon
  x <- c(x,runif(10,0,100))
  y <- c(y,rnorm(10,0,10))
  minMSE1 <- 999999999999
  minMSE2 <- 999999999999
  for(t in tau){
    if(t!=0.5){
      xy <- filtered(t,x,y)
      x1 <- xy[,1]
      y1 <- xy[,2]
      reg <- lm(y1~x1)
      beta1 <- (reg$coefficients[2])
      c1 <- reg$coefficients[1]
      ypred <- beta1*x+c1
      currmse <- MSE(y_pred = ypred, y_true = y)
      minMSE1 <- min(currmse,minMSE1)
    }
    
    reglts <- ltsReg(x,y,alpha=1-t)
    currmse <- MSE(y_pred = (reglts$fitted.values), y_true = y)
    minMSE2 <- min(currmse,minMSE2)
  }
  MSE_TM <- c(MSE_TM,minMSE1)
  MSE_LTS <- c(MSE_LTS,minMSE2)
}

MSE_TM
MSE_LTS
result <- MSE_LTS < MSE_TM
print(sum(result))


