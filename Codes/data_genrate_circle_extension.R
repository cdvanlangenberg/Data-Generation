# library(binhf)
# nn=400
# nl =1 ; nL = 16; N = nL/2;  m = 0:N;
# C1=1 ; C2= 1 ; a=1 ; 

lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )

d <- matrix(2*pi*(0:(nL-1))/nL, ncol=1) %*% matrix(1:(N-1), ncol=(N-1))

# a0 <- (1/(a*pi))*(1-exp(-pi*a))
# 
# cov <- C1*( (exp(-a*abs(lambda))) - a0 )

# cov <- C1*(exp(-a*abs(lambda)))


circulant <- function(x, nrow = length(x)) {
  n <- length(x)
  matrix(x[(1:n - rep(1:nrow, each=n)) %% n + 1L], ncol=n, byrow=TRUE)
}

Cm <- circulant(cov, length(cov))
r= round(eigen(Cm)$values,4)
v= eigen(Cm)$vectors

sqCm <- round(v%*%diag(sqrt(r))%*%t(v),8)

cov_sim <-array(0,c(nn,nL))
c0 <- array(0,nn)

for (hh in 1:nn){X<-NULL
                 set.seed(12345+hh)
                 xx <-rnorm(nL,mean = 0,sd = 1)
                 
                 X <- t(sqCm%*%xx)
                 #X <- array(scale(X[1,]))
                 y <- t(rbind(lambda1*180/pi,X))
                 y <- cbind(y, cos(d), sin(d))
                 #y <- y[sample(1:nL,size = round(.9*nL), replace = F),-c(1,2,nL)]
                 reg<-lm(y[,2]~y[,-c(1,2)])
                 c0[hh] <- reg$coefficients[1]
                 
                 for(jj in 1:nL){Y<-NULL
                                 Y<-rbind(X,shift(X,-(jj-1)))
                                 cov_sim[hh,jj]<-(1/nL)*(sum(Y[1,]*Y[2,])) -mean(Y[1,])*mean(Y[2,]) 
                                   }
}

cov_ave <- colMeans(cov_sim) 
var(c0)
a0


#cov_ave

#print(cov_ave)
plot(lambda, cov, xaxt = "n", ylim=c(min(cov_ave),max(cov_ave, cov)), bty="n", pch=19)
lines(lambda, cov_ave, col="red",bty="n", xaxt="n", type="p", pch=17)
#lines(lambda, cov_ave+mean(abs(c0)), col="blue",bty="n", xaxt="n", type="p", pch=17)
axis(1, at=lambda, labels = lambda*180/pi)
legend("topright", legend=c("formula","simulation"), col=1:2, pch=c(19,17), bty = "n")



#reg<-lmer(X2~1|ID, data = y, REML = F)

########################### inverse fourier 
if(FALSE){
fourier<- fft(z = Z, inverse = F)
Z <- .5*cos(2*lambda)+.25*sin(3*lambda)
plot(lambda, Z, xaxt="n", bty="n")

get.trajectory <- function(X.k,ts,acq.freq) {
  
  N   <- length(ts)
  i   <- complex(real = 0, imaginary = 1)
  x.n <- rep(0,N)           
  ks  <- 0:(length(X.k)-1)
  
  for(n in 0:(N-1)) {      
    x.n[n+1] <- sum(X.k * exp(i*2*pi*ks*n/N)) / N
  }
  
  x.n * acq.freq 
}

X.k <- fft(c(X,X))                  

time     <- 4                           
acq.freq <- 100                         
ts  <- seq(-pi, pi, 1/acq.freq)

x.n <- get.trajectory(X.k,ts,acq.freq) 
plot(ts, Re(x.n))

}
