nn=100 # number of iterations
# set parameter values
nl = 36 ; nL = 48  ; N = nL/2;  m = 0:N;
C1=1 ; C2= 1 ; a=1 ; u=1; p=.5
n <- p^m

#generate values for phi 0 to pi
lim_phi = 2*pi/3 # boundy for phi

phi <- seq( -lim_phi/2, lim_phi/2, by=lim_phi/(nl+0))[-c(1,nl+(2))]
if(nl%%2==1){
phi <- seq( -lim_phi/2, lim_phi/2, by=lim_phi/(nl+1))[-c(1,nl+(2))]
}

# phi <-phi+pi/2

# phi <- seq( -pi/(nl+0), pi, by=pi/(nl+0))[-c(1,nl+(2))]
#phi <- c(40,20)*pi/180
mt <-length(m)
lat1 <- rep(phi, nl)
lat2 =rep(phi, each=nl)
lambda <- 2*(1:nL-1)*pi/nL
#lambda <-ifelse(lambda>pi, lambda-2*pi, lambda)

X <- array(0,c(nl, nL, nn)) # the generated data
var_sim <- array(0, c(nn, nl))
mean_sim <- array(0, c(nn, nl))
Cm = C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
# Cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
nvar <- function(x){(1/(length(x)))*sum(x^2)} # variance function, without mean  
#nvar <- function(x){(1/(length(x)-0))*sum((x-mean(x))^2)}

for(hh in 1:nn){
  
RCm <- array(0,c(rep(length(phi),2),length(m)))
ICm <- array(0,c(rep(length(phi),2),length(m)))
Kv <-  array(0,c(rep(2*length(phi),2),length(m)))
d <-  array(0,c(length(m), 2*length(phi)))
# r <-  matrix(0,c(length(m), 2*length(phi)))
r <-  matrix(0,nrow = length(m), ncol=2*length(phi))

for(kk in 1:length(m)){
#kk=3
  # RCm[,,kk] <- t(array(round(Cm*cos(m[kk]*u*(lat1-lat2))*(p^m[kk]),6), dim = rep(length(phi),2)))
  # ICm[,,kk] <- t(array(round(Cm*sin(m[kk]*u*(lat1-lat2))*(p^m[kk]),6), dim = rep(length(phi),2)))
  RCm[,,kk] <- t(array(Cm*cos(m[kk]*u*(lat1-lat2))*(1*n[kk]), dim = rep(length(phi),2)))
  ICm[,,kk] <- t(array(Cm*sin(m[kk]*u*(lat1-lat2))*(1*n[kk]), dim = rep(length(phi),2)))

if(m[kk]%in% c(0, N)){ICm[,, kk] = 0}
x1 <- cbind(RCm[,,kk], -ICm[,,kk])
x2 <- cbind(ICm[,,kk], RCm[,,kk])
Kv[, , kk] <- round(.5*(rbind(x1, x2)),8)
colnames(d) <- c(paste0("R",1:nl), paste0("I",1:nl))

if(m[kk]%in% c(0,N)){Kv[, , kk] = 0
Kv[1:nl,1:nl,kk]=RCm[,,kk]}

r[kk,]= round(eigen(Kv[,,kk])$values,8)
r[kk,] <- ifelse(r[kk,]<0, 0, r[kk,])
v= eigen(Kv[,,kk])$vectors

# r[1,1:nl]= round(eigen(RCm[,,1])$values,8)
# r[2,]= round(eigen(Kv[,,2])$values,8)
# if(kk>2)
# r[kk,]= r[(kk-1),]/2
# v= eigen(Kv[,,kk])$vectors

sqKv <- round(v%*%diag(sqrt(r[kk,]))%*%t(v),8)
set.seed(12345+hh)
x <-rnorm(n = 2*nl,mean = 0,sd = 1)

d[kk,] <- t(sqKv %*% x)  #*(cos(lat1))^m[kk] ## d is same as the 
rownames(d) <- m
colnames(d) <- c(paste0("R",1:nl), paste0("I",1:nl))
}

#d
   
# for(ii in 1:nl){
#   for(jj in 1:nL){
#     X[ii,jj,hh] <- (d[1,ii]+2*sum(d[2:(mt-1),ii]*cos(m[-c(1,mt)]*lambda[jj])
#                          -d[2:(mt-1),ii+nl]*sin(m[-c(1,mt)]*lambda[jj]) ) + d[mt,ii])
#   }
# }

for(ii in 1:nl){
  for(jj in 1:nL){
    X[ii,jj,hh] <- (d[1,ii]+2*sum(d[2:length(m),ii]*cos(m[-1]*lambda[jj])
                                  -d[2:length(m),ii+nl]*sin(m[-1]*lambda[jj]) ))
  }
}


rownames(X)<-round(phi*180/pi,0)
colnames(X)<- round(lambda*180/pi,0)

var_sim[hh,] <-apply(X[,,hh] ,MARGIN = 1, nvar )
mean_sim[hh,] <-apply(X[,,hh] ,MARGIN = 1, mean )
}


var_sim<-colMeans(var_sim)
mean_sim<-colMeans(mean_sim)


var_phi <- C1*6*(1-exp(-abs(phi)))

# var_phi <- 6*(1-sqrt(1/(1+phi^2))) # for second Cm function

plot(phi, var_phi, type="l", lty=1, col=3, ylim=c(min(var_sim,var_phi), max(var_sim,var_phi)), xaxt = "n", bty="n",ylab="Variance",
     xlab=expression(phi), main=paste("nl=",nl,"nL=",nL, "C1=",C1,"C2=",C2, "p=",p ), sub=paste("# of simulations=",nn))
axis(1, at=phi, labels=round(phi*180/pi,0))
lines(phi, var_sim, type="l", lty=2, col=4, lwd=2)
legend("top", c("formula","simulation"), lty=c(1,2), col=3:4, bty="n", lwd=1:2)



