library(boot)
library(MLmetrics)
library(quantreg )
library(robustbase)
data <- read.csv("DJFranses.csv")
nboot <- 100
#data <-  as.data.frame(data)

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

tau = seq(0.05,0.5,0.05)

boo <- function(data,indices)
{
  dt <- data[indices,]
  #print(dt)
  x <- as.matrix(dt["time"])
  y <- as.matrix(dt["value"])
  #print(y)
  minMSE1 <- 999999999
  minMSE2 <- 999999999
  b1 <- 0
  b2 <- 0
  co1 <- 0
  co2 <- 0
  for(t in tau)
  {
    if(t!=0.5){
      filter <- filtered(t,x,y)
      x1 <- filter[,1]
      y1 <- filter[,2]
      reg = lm(y1~x1)
      beta1 <- (reg$coefficients[2])
      c1 <- reg$coefficients[1]
      ypred <- beta1*x+c1
      currmse <- MSE(y_pred = ypred, y_true = y)
      minMSE1 <- min(currmse,minMSE1)
      if(minMSE1 == currmse){
        b1 <- beta1
        co1 <- c1 
      }
      
    }
    reglts <- ltsReg(x,y,alpha=1-t)
    currmse <- MSE(y_pred = (reglts$fitted.values), y_true = y)
    minMSE2 <- min(currmse,minMSE2)
    if(minMSE2 ==currmse){
      b2 <- reglts$coefficients[2]
      co2 <- reglts$coefficients[1]
    }
    
  }
  c(minMSE1,minMSE2,b1,co1,b2,co2)
  
}
myBootstrap <- boot(data, boo, R=nboot)
myBootstrap$t
MSE_TM <- myBootstrap$t[,1]
MSE_LTS <- myBootstrap$t[,2]
result <- MSE_LTS < MSE_TM
sum(result)

b1 <- myBootstrap$t[,3]
c1 <- myBootstrap$t[,4]
b2 <- myBootstrap$t[,5]
c2 <- myBootstrap$t[,6]

x <- as.matrix(data["time"])
y <- as.matrix(data["value"])

bet1 <- 0
bet2 <- 0
coe1 <- 0
coe2 <- 0
minMSE1 <- 99999999
minMSE2 <- 99999999
for(t in tau)
{
  if(t!=0.5){
    filter <- filtered(t,x,y)
    x1 <- filter[,1]
    y1 <- filter[,2]
    reg = lm(y1~x1)
    beta1 <- (reg$coefficients[2])
    cd1 <- reg$coefficients[1]
    ypred <- beta1*x+cd1
    currmse <- MSE(y_pred = ypred, y_true = y)
    minMSE1 <- min(currmse,minMSE1)
    if(minMSE1 == currmse){
      bet1 <- beta1
      coe1 <- cd1 
    }
    
  }
  reglts <- ltsReg(x,y,alpha=1-t)
  currmse <- MSE(y_pred = (reglts$fitted.values), y_true = y)
  minMSE2 <- min(currmse,minMSE2)
  if(minMSE2 ==currmse){ 
    bet2 <- reglts$coefficients[2]
    coe2 <- reglts$coefficients[1]
  }
  
}

BootstrapMSE1 <- (sum((b1-bet1)^2)+sum((c1-coe1)^2))/nboot
BootstrapMSE2 <- (sum((b2-bet2)^2)+sum((c2-coe2)^2))/nboot
BootstrapMSE1
BootstrapMSE2
print(BootstrapMSE1<BootstrapMSE2)

plot(x,y,main="TM and LTS Regression for DJIA historical price data",xlab="Time",ylab="Price")
lines(abline(reg,col="red",lty=1))
lines(abline(reglts,col="blue",lty=2))
legend(1980,4000,legend=c("Trimmed Mean Regression","LTS Regression"), col=c("red","blue"),lty=c(1,2), ncol=1)
minMSE1
minMSE2
