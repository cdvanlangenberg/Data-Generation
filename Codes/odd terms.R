# nl =2 ; nL = 76; N = nL/2;  m = 0:N;
# C1=1 ; C2= 1 ; a=1 ; u=1; p=.5
# 
# phi <- c(10,60)*pi/180
# mycov <-function(lat1, lat2)
# {#cm=C1*(C2-exp(-a*abs(lat1))-exp(-a*abs(lat2))+exp(-a*abs(lat1-lat2)))
#   cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
#   ;cm}

lat1 <- rep(phi, nl)
lat2 =rep(phi, each=nl)

lambda1 <- 2*(1:nL-1)*pi/nL
lambda <- ifelse(lambda1>pi, lambda1-2*pi, lambda1 )

l1 <- lambda[1:(1+nL/2)]
l2 <- c(0,lambda[nL:(1+nL/2)])


par(mfrow=c(2,2), cex.main = .9 )

################## covarince function 1

v <- mycov(phi[1],phi[2])
v1 <- v*(1-p^2)/(1-2*p*cos(l1+u*(phi[1]-phi[2]))+p^2)
v2 <- v*(1-p^2)/(1-2*p*cos(l2+u*(phi[1]-phi[2]))+p^2)
plot(l1, .5*(v1+v2), pch=19, col=6, xaxt="n", ylab=" ", xlab=expression(paste(Delta,lambda)), main="Covariance function 1", bty="n", ylim=c(min(.5*(v1-v2),.5*(v1+v2)), max(.5*(v1-v2),.5*(v1+v2))))
axis(1, at=l1, labels=round(l1*180/pi,0))
lines(l1, .5*(v1-v2), pch=15, col=7)
legend("topright", cex=.7, legend = c(".5(C(+h)+C(-h))",".5(C(+h)-C(-h))"), pch=c(19,15), col=6:7, bty = "n")
abline(h = 0)


################## covarince function 2

v1 <- 1*v*log(1/(2*(1-cos(l1+u*(phi[1]-phi[2])))))
v2 <- 1*v*log(1/(2*(1-cos(l2+u*(phi[1]-phi[2])))))

plot(l1, .5*(v1+v2), pch=19, col=6, xaxt="n", ylab=" ", xlab=expression(paste(Delta,lambda)), main="Covariance function 2", bty="n", ylim=c(min(.5*(v1-v2),.5*(v1+v2)), max(.5*(v1-v2),.5*(v1+v2))))
axis(1, at=l1, labels=round(l1*180/pi,0))
lines(l1, .5*(v1-v2), pch=15, col=7)
legend("topright", cex=.7, legend = c(".5(C(+h)+C(-h))",".5(C(+h)-C(-h))"), pch=c(19,15), col=6:7, bty = "n")
abline(h = 0)

################## covarince function 3
theta <- u*(phi[1]-phi[2])
theta <- ifelse(theta<0, theta+2*pi, theta )

theta1 <- ifelse( (l1+u*(phi[1]-phi[2]))<0, (l1+u*(phi[1]-phi[2]))+2*pi, (l1+u*(phi[1]-phi[2])) )
theta2 <- ifelse((l2+u*(phi[1]-phi[2]))<0, (l2+u*(phi[1]-phi[2]))+2*pi, (l2+u*(phi[1]-phi[2])) )

v1 <- v*(pi^4/90 - (pi*theta1)^2/12 + pi*theta1^3/12 -(theta1)^4/48)
v2 <- v*(pi^4/90 - (pi*theta2)^2/12 + pi*theta2^3/12 -(theta2)^4/48)

plot(l1, .5*(v1+v2), pch=19, col=6, xaxt="n", ylab=" ", xlab=expression(paste(Delta,lambda)), main="Covariance function 3", bty="n", ylim=c(min(.5*(v1-v2),.5*(v1+v2)), max(.5*(v1-v2),.5*(v1+v2))))
axis(1, at=l1, labels=round(l1*180/pi,0))
lines(l1, .5*(v1-v2), pch=15, col=7)
legend("topright", cex=.6, legend = c(".5(C(+h)+C(-h))",".5(C(+h)-C(-h))"), pch=c(19,15), col=6:7, bty = "n")
abline(h = 0)

################## covarince function 4

v1 <- .5*v*log(1/(1-2*p*cos(theta1)+p^2))
v2 <- .5*v*log(1/(1-2*p*cos(theta2)+p^2))

plot(l1, .5*(v1+v2), pch=19, col=6, xaxt="n", ylab=" ", xlab=expression(paste(Delta,lambda)), main="Covariance function 4", bty="n", ylim=c(min(.5*(v1-v2),.5*(v1+v2)), max(.5*(v1-v2),.5*(v1+v2))))
axis(1, at=l1, labels=round(l1*180/pi,0))
lines(l1, .5*(v1-v2), pch=15, col=7)
legend("topright", cex=.6, legend = c(".5(C(+h)+C(-h))",".5(C(+h)-C(-h))"), pch=c(19,15), col=6:7, bty = "n")
abline(h = 0)

mtext(side=3,line = -1, outer = T, text = paste("Odd terms for Phi =", phi[1]*180/pi,"&",phi[2]*180/pi), cex=1.2 )

