####### this code will generate the grided data for 2 latitudes X given number of longitudes
####### and compare the cross covarince

library(binhf)
nn=100 # number of iterations

nl =2 ; nL = 36; N = nL/2;  m = 0:(1*N);
C1=1 ; C2= 1 ; a=1 ; u=1; p=.5; m1 =2:3

######## select a model ###################
# n <- p^m        # model1
n <- p^m/m      # model2
# n <- (1/m)^4    # model3
# n <- 1/(2*m)    # model4

######## select two latitudes #############
phi <- c(30,60)*pi/180

mycov <-function(lat1, lat2)
  {cm=C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
   #cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
;cm}

lat1 <- rep(phi, nl)
lat2 =rep(phi, each=nl)
lambda <- 2*(1:nL-1)*pi/nL

Cm = mycov(lat1, lat2)
X <- array(0,c(nl, nL, nn)) # the generated data
var_sim <- array(0, c(nn, nL))
mean_sim <- array(0, c(nn, nl))
cor_sim <- array(0, c(nn, 1))

for(hh in 1:nn){
  
  RCm <- array(0,c(rep(length(phi),2),length(m)))
  ICm <- array(0,c(rep(length(phi),2),length(m)))
  Kv <-  array(0,c(rep(2*length(phi),2),length(m)))
  d <-  array(0,c(length(m), 2*length(phi)))
  r <-  array(0,c(length(m), 2*length(phi)))
  
  for(kk in 1:length(m)){

    RCm[,,kk] <- t(array(Cm*cos(m[kk]*u*(lat1-lat2))*(1*n[kk]), dim = rep(length(phi),2)))
    ICm[,,kk] <- t(array(Cm*sin(m[kk]*u*(lat1-lat2))*(1*n[kk]), dim = rep(length(phi),2)))
    
    if(n==Inf && m[kk]%in% c(0)){RCm[,,kk]=ICm[,, kk] = 0} # make this adjustment when C0 = 0
    
    x1 <- cbind(RCm[,,kk], -ICm[,,kk])
    x2 <- cbind(ICm[,,kk], RCm[,,kk])
    Kv[, , kk] <- .5*(rbind(x1, x2))
    colnames(d) <- c(paste0("R",1:nl), paste0("I",1:nl))
    
    if(m[kk]%in% c(0)){Kv[, , kk] = 0
                       Kv[1:nl,1:nl,kk]=RCm[,,kk]}
    
    r[kk,]= eigen(Kv[,,kk])$values
    r[kk,] <- ifelse(r[kk,]<0, 0, r[kk,])
    v= eigen(Kv[,,kk])$vectors
    
    sqKv <- v%*%diag(sqrt(r[kk,]))%*%t(v)
    set.seed(12345+hh)
    
    x <-rnorm(n=2*nl,mean = 0,sd = 1)
    
    d[kk,] <- t(sqKv %*% x)  #*(cos(lat1))^m[kk]
    rownames(d) <- m
    colnames(d) <- c(paste0("R",1:nl), paste0("I",1:nl))
  }
  
  
  for(ii in 1:nl){
    for(jj in 1:nL){
      X[ii,jj,hh] <- (d[1,ii]+2*sum(d[2:length(m),ii]*round(cos(m[-1]*lambda[jj]),4)
                                    -d[2:length(m),ii+nl]*round(sin(m[-1]*lambda[jj]),4) )) + m1[ii]  #*cos(phi[ii])
    }
  }
  
  rownames(X)<-round(phi*180/pi,0)
  colnames(X)<- round(lambda*180/pi,0)
  
  
  for(jj in 1:nL){Y<-NULL
                  Y<-rbind(X[1,,hh],shift(X[2,,hh],-(jj-1)))
                  var_sim[hh,jj]<-(1/nL)*(sum(Y[1,]*Y[2,])) -mean(Y[1,])*mean(Y[2,])
  }
  
mean_sim[hh,] <-apply(X[,,hh] ,MARGIN = 1, mean )
cor_sim[hh,] <-cor(X[1,,hh], X[2,,hh])

}  

var_sim<-colMeans(var_sim)

theta <- lambda+u*(phi[1]-phi[2])
theta <- ifelse(theta<0, theta+2*pi, theta )

######## select the corresponding theoritical covarince model #############
######## for the above selected model                         #############

# var_phi <- mycov(phi[1],phi[2])*(1-p^2)/(1-2*p*cos(theta)+p^2) #model1
var_phi <- mycov(phi[1],phi[2])*log(1/(1-2*p*cos(lambda+u*(phi[1]-phi[2]))+p^2)) #model2
# var_phi <- 2*mycov(phi[1],phi[2])*(pi^4/90 - (pi*theta)^2/12 + pi*theta^3/12 -(theta)^4/48) #model3
# var_phi <- .5*mycov(phi[1],phi[2])*log(1/(2*(1-cos(theta)))) #model4

########### estimate C0
# c0<-( (mean_sim[,1]-0)*(mean_sim[,2]-0))
# c0<-mean(c0)
# mycov(phi[1],phi[2])
# c0
#######################


plot(lambda, var_phi, type="l", lty=2, col=3, ylim=c(min(var_phi, var_sim), max(var_phi, var_sim)), xaxt = "n", bty="n",ylab="Covariance",
     xlab=expression(lambda), main=paste("nL =", nL, ",latitudes ",phi[1]*180/pi,"&", phi[2]*180/pi), sub=paste("#simulations=",nn))
axis(1, at=lambda, labels=round(lambda*180/pi,0))
lines(lambda, var_sim, type="l", lty=3, col=4, lwd=2)
legend("top", c("Theoretical","Emperical"), lty=c(2,3), col=3:4, bty="n", lwd=1:2)


