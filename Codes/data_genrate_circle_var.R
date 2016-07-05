# library(binhf)
# nn=100
# nl =1 ; nL = 24; N = nL/2;  m = 0:N;
# C1=1 ; C2= 1 ; a=1 ; 


lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )

cov <- C1*(exp(-a*abs(lambda)))

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


var=2*(1-exp(-a*abs(lambda)))
plot(lambda, var, xaxt = "n", ylim=c(min(cov_ave),3), bty="n", pch=19, main=paste("Variogram; C1=",C1, "nL=",nL))
lines(lambda, cov_ave, col="red",bty="n", xaxt="n", type="p", pch=17)
axis(1, at=lambda, labels = lambda*180/pi)
legend("topright", legend=c("formula","simulation"), col=1:2, pch=c(19,17), bty = "n")


