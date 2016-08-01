# switching the estimation on the sphere to estimate variogram using cross covarince 


# library(binhf)
# nn=200 # number of iterations
# 
# nl =2 ; nL = 48; N = nL/2;  m = 0:N;
# C1=1 ; C2= 1 ; a=1 ; u=1; p=.5
# n=2/m
# #n = (p^m)/m
# 
# phi <- c(10,20)*pi/180
# mycov <-function(lat1, lat2)
# {#cm=C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
#  cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
#  ;cm}

lat1 <- rep(phi, nl)
lat2 =rep(phi, each=nl)

lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )
#lambda <- 2*(1:nL-1)*pi/nL

X <- array(0,c(nl, nL, nn)) # the generated data
var_sim <- array(0, c(nn, nL))


for(hh in 1:nn){
  
  RCm <- array(0,c(rep(length(phi),2),length(m)))
  ICm <- array(0,c(rep(length(phi),2),length(m)))
  Kv <-  array(0,c(rep(2*length(phi),2),length(m)))
  d <-  array(0,c(length(m), 2*length(phi)))
  r <-  array(0,c(length(m), 2*length(phi)))
  
  for(kk in 1:length(m)){
    Cm = mycov(lat1, lat2)
    
    RCm[,,kk] <- t(array(round(Cm*cos(m[kk]*u*(lat1-lat2))*(n[kk]),6), dim = rep(length(phi),2)))
    ICm[,,kk] <- t(array(round(Cm*sin(m[kk]*u*(lat1-lat2))*(n[kk]),6), dim = rep(length(phi),2)))
    
    if(m[kk]%in% c(0)){RCm[,,kk]=ICm[,, kk] = 0}
    x1 <- cbind(RCm[,,kk], -ICm[,,kk])
    x2 <- cbind(ICm[,,kk], RCm[,,kk])
    Kv[, , kk] <- round(.5*(rbind(x1, x2)),8)
    colnames(d) <- c(paste0("R",1:nl), paste0("I",1:nl))
    
    if(m[kk]%in% c(0)){Kv[, , kk] = 0
                         Kv[1:nl,1:nl,kk]=RCm[,,kk]}
    
    r[kk,]= round(eigen(Kv[,,kk])$values,4)
    v= eigen(Kv[,,kk])$vectors
    
    sqKv <- round(v%*%diag(sqrt(r[kk,]))%*%t(v),8)
    set.seed(12345+hh)
    x <-rnorm(2*nl)
    
    d[kk,] <- t(sqKv %*% x) #*(cos(lat1))^m[kk]
    rownames(d) <- m
    colnames(d) <- c(paste0("R",1:nl), paste0("I",1:nl))
  }
  
  
  for(ii in 1:nl){
    for(jj in 1:nL){
      X[ii,jj,hh] <- (d[1,ii]+2*sum(d[2:length(m),ii]*round(cos(m[-1]*lambda[jj]),4)
                                   -d[2:length(m),ii+nl]*round(sin(m[-1]*lambda[jj]),4) )) #*cos(phi[ii])
    }
  }
  
  rownames(X)<-round(phi*180/pi,0)
  colnames(X)<- round(lambda*180/pi,0)


#   for(jj in 1:nL){Y<-NULL
#   Y<-rbind(X[1,,hh],shift(X[2,,hh],-(jj-1)))
#   #var_sim[hh,jj]<-(1/nL)*(sum(Y[1,]*Y[2,]))-mean(Y[1,])*mean(Y[2,]) # covarince
#   #var_sim[hh,jj]<-(1/nL)*sum((Y[1,]-Y[2,])^2)
#   }
  
  for(jj in 1:nL){Y1=Y2<-NULL            
                  Y1<-array(shift(X[1,,hh],-(jj-1))-X[1,,hh], nL)
                  Y2<-array(shift(X[2,,hh],-(jj-1))-X[2,,hh], nL)
                  var_sim[hh,jj]<-(1/(2*nL))*sum(Y1*Y2)
                  
  }
}  
var_sim<-colMeans(var_sim)
#var_sim

lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )

# theta <- u*(phi[1]-phi[2])
# theta <- ifelse(theta<0, theta+2*pi, theta )

l1 <- lambda[1:(1+nL/2)]
l2 <- c(0,lambda[nL:(1+nL/2)])

v <- mycov(phi[1],phi[2])
v0 <- 1*v*log(1/(2*(1-cos(u*(phi[1]-phi[2])))))
v1 <- 1*v*log(1/(2*(1-cos(l1+u*(phi[1]-phi[2])))))
v2 <- 1*v*log(1/(2*(1-cos(l2+u*(phi[1]-phi[2])))))
var_phi <- .5*(v0 - .5*(v1+v2))



plot(lambda[lambda>=0], var_sim[lambda>=0], type="o", lty=2, col=3, xaxt = "n", bty="n",ylab="Variance",ylim=c(min(var_sim, var_phi),max(var_sim,var_phi)),pch=19,
     xlab=expression(paste(Delta,lambda)),main=paste("nl=",nl,",nL=",nL, ",phi=",phi[1]*180/pi," & ", phi[2]*180/pi,",p=",p,", adjustment"), sub=paste("#simulations=",nn))
axis(1, at=lambda, labels=round(lambda*180/pi,0))
lines(l1, var_phi, type="o", lty=3, col=4, lwd=2, pch=15)
legend("bottomright", c("simulation","formula"), pch=c(19,15), col=3:4, bty="n", lwd=1:2)



