
phi
lat1 <- rep(phi, nl)
lat2 =rep(phi, each=nl)

lambda <- 2*(1:nL-1)*pi/nL

theta <- array(rep(u*(lat1-lat2),nL), c(nl,nl,nL))+rep(lambda, each=nl*nl)
theta <- ifelse(theta<0, theta+2*pi, theta )
Cm <- C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
RPQ = array(rep(Cm,nL), c(nl,nl,nL))* (1-p^2)/(1-2*p*cos(theta)+p^2)

for(hh in 1:nn){
  
  cov_sim <- array(NA, c(nl,nl,nL))
  for(kk in 1:nL){
    #Y<- data.frame(X[,,hh])
    Y <- data.frame(DD[,,hh])
    Yt<- t(Y[, shift(names(Y),-(kk-1))])
    cov_sim[,,kk] <- (1/nL)*(as.matrix(Y)%*%as.matrix(Yt))
  
    
    }
  

}


f1 <- function(x) x/8
f2<-function(x) x^2/16
curve(f1, 0,2, xlim=c(0,4), ylim=c(0,1))
curve(f2, 2,4, add = T)


f1 <- function(x) 1/8
f2<-function(x) x/8
curve(f2, 2,4, xlim=c(0,4), ylim=c(0,1))
abline(h=1/8, v=2)
abline(v=c(1,3), lty=2)
curve(f2, 2,4, add = T)












theta <- lambda+u*(phi[1]-phi[2])
theta <- ifelse(theta<0, theta+2*pi, theta )
RPQ <- Cm[2]*(1-p^2)/(1-2*p*cos(theta)+p^2)


# cov_simulate <- array(0, c(nn,nL))

cov_simulate <- array(0, c(nn,nn, nL))

for(hh in 1:nn){

cov_sim <- array(NA, c(nl,nl,nL))
for(kk in 1:nL){
     Y<- data.frame(X[,,hh])
     Yt<- t(Y[, shift(names(Y),-(kk-1))])
      cov_sim[,,kk] <- (1/nL)*(as.matrix(Y)%*%as.matrix(Yt))
}

cov_simulate[hh,] <-cov_sim[1,2,]
}

d1<-colMeans(cov_simulate)
plot(d1)


min.RSS <-function(data, par) {
  # par=c(1,1,1,1,.5)
  Cm <- par[1]*(par[2]-exp(-par[3]*abs(lat1))-exp(-par[3]*abs(lat2))+exp(-par[3]*abs(lat1-lat2)))
  theta <- array(rep(par[4]*(lat1-lat2),nL), c(nl,nl,nL))+rep(lambda, each=nl*nl)
  theta <- ifelse(theta<0, theta+2*pi, theta )
  RPQ = array(rep(Cm,nL), c(nl,nl,nL))* (1-par[5]^2)/(1-2*par[5]*cos(theta)+par[5]^2)[2,1,]
  sum((RPQ[1,2,]-data)^2)
}

OP <- array(0, c(nn, 5))
for(hh in 1:nn){
OPT<- optim(par = c(.9,.9,.9,.9,.4), fn = min.RSS, data=cov_simulate[hh,])
OP[hh,]<-OPT$par
#OPT$value
}


#########################


theta <- lambda+par[4]*(phi[1]-phi[2])
theta <- ifelse(theta<0, theta+2*pi, theta )
PRQ <- Cm[2]*(1-par[5]^2)/(1-2*par[5]*cos(theta)+par[5]^2)


# for(kk in 1:nL){
#   for(ii in 1:nl){
#     for(jj in 1:nl){
#       Y<-rbind(X[ii,,1],shift(X[jj,,1],-(kk-1)))
#       cov_sim[ii,jj,kk] <- (1/nL)*(sum(Y[1,]*Y[2,]))
#     }
#   }
# }
