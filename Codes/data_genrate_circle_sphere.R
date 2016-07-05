# variogram estimation on each individual circle on the sphere 

# library(binhf)
# nn=100
# nl =1 ; nL = 40; N = nL/2;  m = 0:N;
# C1 =1 ; C2= 1 ; a=1 ; 
# 
# phi =70
# 
# mycov <- function(phi){
#   cm <- C1*(exp(-a*abs(lambda)))
#   #cm=C1*(exp(-a*abs(lambda)/cos(phi)))
#   #cm=C1*(1-exp(-a*abs(phi)))*(1-p^2)/(1-2*p*cos(lambda)+p^2)
#   ;cm}

phi <- phi*pi/180
lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )

cov <- mycov(phi)

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

# myvar <- function(lambda){
#   v=2*(1-exp(-a*abs(lambda)))
#   v=2*(1-exp(-a*abs(lambda)/cos(phi)))
#   #v=2*12*(1-exp(-a*abs(phi)))*(1-cos(lambda))/(5-4*cos(lambda))
#   ;v
# }

var <- myvar(lambda)
# var=2*(1-exp(-a*abs(lambda)/cos(phi)))
# var=2*12*(1-exp(-a*abs(phi)))*(1-cos(lambda))/(5-4*cos(lambda))

plot(lambda, var, xaxt = "n", ylim=c(0,max(cov_ave)), bty="n", pch=19, main=paste("Variogram; Lat=",phi*180/pi,", nL=",nL))
lines(lambda, cov_ave, col="red",bty="n", xaxt="n", type="p", pch=17)
axis(1, at=lambda, labels = lambda*180/pi)
legend("bottomright", legend=c("formula","simulation"), col=1:2, pch=c(19,17), bty = "n")

