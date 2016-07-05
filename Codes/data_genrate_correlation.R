n <- 5
p <- .6
k <- 400
cov <- cbind(diag(n), diag(p,n))
cov <- rbind(cov,cbind(diag(p,n),diag(n)))


r= round(eigen(cov)$values,4)
v= eigen(cov)$vectors

sqcov <- round(v%*%diag(sqrt(r))%*%t(v),8)

p_hat <- array(0,k)
for(ii in 1:k){x<-NULL
x <- rnorm(2*n)
X <- sqcov%*%x
X1 <- X[1:n,]
X2 <- X[(n+1):(2*n),]
p_hat[ii] <- (1/n)*sum((X1-mean(X1))*(X2-mean(X2)))
}

mean(p_hat)
