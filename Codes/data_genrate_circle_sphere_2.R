# variogram estimation on each individual circle on the sphere, using the 3rd covarince function  

library(binhf)
nn=200
nl =1 ; nL = 24; N = nL/2;  m = 0:N;
C1 =1 ; C2= 1 ; a=1 ;p=.5 

phi = 20

# mycov <- function(phi){
#   cm=C1*(exp(-a*abs(lambda)/cos(phi)))
#   #cm=C1*(1-exp(-a*abs(phi)))*(1-p^2)/(1-2*p*cos(lambda)+p^2)
#   ;cm}

phi <- phi*pi/180
lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )

cov <- C1*(1-exp(-a*abs(phi)))*(1-p^2)/(1-2*p*cos(lambda)+p^2)

circulant <- function(x, nrow = length(x)) {
  n <- length(x)
  matrix(x[(1:n - rep(1:nrow, each=n)) %% n + 1L], ncol=n, byrow=TRUE)
}

Cm <- circulant(cov)
r= round(eigen(Cm)$values,4)
v= eigen(Cm)$vectors

sqCm <- round(v%*%diag(sqrt(r))%*%t(v),8)

cov_sim <-array(0,c(nn,nL))

for (hh in 1:nn){
  set.seed(12345+hh)
  x <-rnorm(nL,mean = 0,sd = 1)
  x <- (x-mean(x))/sd(x)
  # mean(x)
  # sd(x)
  X= t(sqCm%*%x)
  
  for(jj in 1:nL){Y<-NULL
                  #jj=2
                  Y<-rbind(X,shift(X,-(jj-1)))
                  cov_sim[hh,jj]<-(1/nL)*sum((Y[1,]-Y[2,])^2)  #-mean(Y[1,])*mean(Y[2,])                
  }
}

cov_ave <- colMeans(cov_sim) # this is the varince


var=2*12*(1-exp(-a*abs(phi)))*(1-cos(lambda))/(5-4*cos(lambda))

data<-data.frame(cbind(lambda, var, cov_ave))
data<-data[order(data$lambda),]

plot(data$lambda, data$var, xaxt = "n", ylim=c(0,max(cov_ave)),type="l", bty="n", lty=1, main=paste("Variogram; Lat=",phi*180/pi,", nL=",nL))
lines(data$lambda, data$cov_ave, col="red",bty="n", xaxt="n", type="l", lty=2)
axis(1, at=data$lambda, labels = data$lambda*180/pi)
legend("bottomright", legend=c("formula","simulation"), col=1:2, lty=c(1,2), bty = "n")

