
# get W_m ---------------------------------------------

#W_m <- fft(X[1,,10], inverse = T)/nL

critical <- qchisq(.95, 2*nl*(N))

Hy <- array(NA, nn)
for(kk in 1:nn){
W_m <- Conj(t(apply(X[,,kk], MARGIN = 1, function(x) fft(x, inverse = T)/nL))[,1:length(m)])
W_m <- as.matrix(t(cbind(W_m[,-ncol(W_m)],W_m[,ncol(W_m)]/2)))
V_m <- as.matrix(cbind(Re(W_m), Im(W_m)))

colnames(V_m) = colnames(d)
rownames(V_m) = rownames(d)



Y_m <- array(0, length(m))

#Y_m[1] <- t(Re(W_m[1,]))%*%solve(RCm[,,1])%*%Re(W_m[1,])
Y_m[length(m)] <- t(Re(W_m[length(m),]))%*%solve(RCm[,,length(m)])%*%Re(W_m[length(m),])

for(ii in 2:N){

  #if (m[ii] %in% c(0,N)){Y_m[ii] <- t(Re(W_m[ii,]))%*%solve(RCm[,,ii])%*%Re(W_m[ii,])}
  
    Y_m[ii] <- t(V_m[ii,])%*%solve(Kv[,,ii])%*%(V_m[ii,])
  
}
#print(Y_m)

Hy[kk] <- ifelse(sum(Y_m)>=critical, 1, 0)

}

sum(Hy)
sum(Hy)/nn


# Generate data from R(P,Q) and compare with mom estimator --------------------------------------

# library(foreach)
# library(doParallel)
# library(binhf)

# initialize the parameters  
nn = 100
nl = 3 ; nL = 8 ; N = nL/2;  m = 0:(1*N);
C1=1 ; C2= 1 ; a=1 ; u=1; p=.5; m1 =2:3

phi <- c(30, 40, 60)*pi/180

lat1 <- rep(phi, nl)
lat2 =rep(phi, each=nl)
lambda <- 2*(1:nL-1)*pi/nL


theta <- array(rep(u*(lat1-lat2),nL), c(nl,nl,nL))+rep(lambda, each=nl*nl)
theta <- ifelse(theta<0, theta+2*pi, theta )
Cm <- C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
RR <- array(rep(Cm,nL), c(nl,nl,nL))* (1-p^2)/(1-2*p*cos(theta)+p^2)

x <- matrix((1:nL - rep(1:nL, each=nL)) %% nL + 1L, byrow=T,ncol=nL)

# get the block circulant RPQ matrix

RPQ <-as.data.frame(RR[,,1:nL]) 

for(ii in 2:nL)
  {
  RPQ <- rbind(RPQ, as.data.frame(RR[,,x[ii,]])) 
}

rm(x)

# SVD to get R^(1/2)
r <- eigen(data.matrix(RPQ))$values
v <- eigen(data.matrix(RPQ))$vectors

sqRPQ <- v%*%diag(sqrt(r))%*%t(v)

set.seed(123)
Z <- as.data.frame(replicate(nn, expr = sqRPQ%*%matrix(rnorm(nl*nL), ncol=1)))

# assign 3 cores to process the following

cl <- makeCluster(3)
registerDoParallel(cl)  


MM1 <- foreach( hh = 1:nn, .combine=list ) %dopar% {

  library(binhf)

  mom <- array(NA, c(nL, nl*(nl+1)/2))

for(kk in 1:nL){
  Y <- matrix(data.matrix(Z[,hh]), ncol=nL) 
  Yt<- t(Y[, shift(1:nL,-(kk-1))])
  mo <- (1/nL)*(as.matrix(Y)%*%as.matrix(Yt)) #- rowMeans(Y)%*%t(rowMeans(Y))
  mom[kk,] <- mo[lower.tri(mo, diag = T)]
}

; mom

}


MM <- array(unlist(MM1), c(nL,nl*(nl+1)/2,nn))

mycov <- function(lat1, lat2) {
  cm=C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
;cm
}

# a = phi1 and b = phi2
a = 3
b = 3

# to get the corresponding latitude combination 
x <- matrix(0, ncol=3, nrow = 3)
x[lower.tri(x, diag = T)] <- 1:(nl*(nl+1)/2)

var_phi <- mycov(phi[a],phi[b])*(1-p^2)/(1-2*p*cos(theta)+p^2)
var_sim <- apply(MM, MARGIN = c(1,2), mean)

plot(lambda, var_phi[b,a,], type="l", lty=2, col=3, xaxt = "n", bty="n",ylab="Covariance",
     xlab=expression(lambda), main=paste("nL =", nL, ",latitudes ",phi[a]*180/pi,"&", phi[b]*180/pi), sub=paste("#simulations=",nn))
axis(1, at=lambda, labels=round(lambda*180/pi,0))
lines(lambda, var_sim[,(x[b,a])], type="l", lty=3, col=4, lwd=2)
legend("top", c("Theoretical","Emperical"), lty=c(2,3), col=3:4, bty="n", lwd=1:2)




# generate Cm from generated data ---------------------

# source("data_genrate_modified.R")

library(foreach)
library(doParallel)

cl <- makeCluster(3)
registerDoParallel(cl)  


mom <- array(NA, c(nL, nl*(nl+1)/2, nn))

MM1 <- foreach( hh = 1:nn, .combine=list ) %dopar% {

  library(binhf)
  
  #mom <- array(NA, c(nL, nl*(nl+1)/2))
  
  for(kk in 1:nL){
    Y<- data.frame(X[,,hh])
    Yt<- t(Y[, shift(names(Y),-(kk-1))])
    mo <- (1/nL)*(as.matrix(Y)%*%as.matrix(Yt))
    mom[kk,,hh] <- mo[lower.tri(mo, diag = T)]
  }
  
  #mom[,,hh]
  Cm.full <- (apply(mom[,,hh], MARGIN=2, FUN=function(x) fft(x, inverse = T)/nL)[1:length(m),])
    }

MM <- array(unlist(MM1), c(length(m),nl*(nl+1)/2,nn))

Sim.RCm <- array(0, c(length(m),nl*(nl+1)/2,nn))
Sim.ICm <- array(0, c(length(m),nl*(nl+1)/2,nn))

for(ii in 1:nn){
Sim.RCm[,,ii] <- t(apply(MM[,,ii], MARGIN = 1, Re))
Sim.ICm[,,ii] <- t(apply(MM[,,ii], MARGIN = 1, Im))
}

# average the simulated over replications 

S.RCm <- apply(Sim.RCm, MARGIN = c(1,2), FUN=mean)
S.ICm <- apply(Sim.ICm, MARGIN = c(1,2), FUN=mean)

# rearrange true values of RCm and ICm

T.RCm <- array(0, c(length(m),nl*(nl+1)/2))
T.ICm <- array(0, c(length(m),nl*(nl+1)/2))

for( ii in seq_along(m)) {
  T.RCm[ii,] <- RCm[,,ii][lower.tri(RCm[,,ii], diag = T)]
  T.ICm[ii,] <- ICm[,,ii][lower.tri(ICm[,,ii], diag = T)]
}

# to get the corresponding latitude combination 
x <- matrix(0, nl, nl)
x[lower.tri(x, diag = T)] <- 1:(nl*(nl+1)/2)

# a = phi1 and b = phi2
a = 1 # a is the column number and b is the row number
b = 4

if(a>b|b>nl|a>nl) stop("Incorrect selection \n check: a < b and a,b < latitudes")

# pdf("Cm plots2.pdf", width = 7, height = 7)
# par(mfrow=c(2,1))

plot(m, T.RCm[,x[b,a]], type="l", lty=2, col=3, bty="n",ylab="Covariance", xaxt="n",
     ylim = c(min(T.RCm[,x[b,a]], S.RCm[,x[b,a]]), max(T.RCm[,x[b,a]], S.RCm[,x[b,a]])),
     xlab="m", main=paste("RCm, nL =", nL, ",latitudes ",phi[a]*180/pi,"&", phi[b]*180/pi), sub=paste("#simulations=",nn))
axis(1, at=m, labels = m)
lines(m, S.RCm[,x[b,a]], type="l", lty=3, col=4, lwd=2)
legend("top", c("Theoretical","Simulated"), lty=c(2,3), col=3:4, bty="n", lwd=1:2)

plot(m, T.ICm[,x[b,a]], type="l", lty=2, col=3, bty="n",ylab="Covariance",xaxt="n",
     ylim = c(min(T.ICm[,x[b,a]], S.ICm[,x[b,a]]), max(T.ICm[,x[b,a]], S.ICm[,x[b,a]])),
     xlab="m", main=paste("ICm, nL =", nL, ",latitudes ",phi[a]*180/pi,"&", phi[b]*180/pi), sub=paste("#simulations=",nn))
axis(1, at=m, labels=m)
lines(m, S.ICm[,x[b,a]], type="l", lty=3, col=4, lwd=2)
legend("bottom", c("Theoretical","Simulated"), lty=c(2,3), col=3:4, bty="n", lwd=1:2)

# dev.off()



























