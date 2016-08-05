####### this code will generate the grided data for 2 latitudes X given number of longitudes
####### and compare with cross variogram 

if (!require("binhf")) install.packages("binhf")
library(binhf)

nn = 300 # number of iterations

nl = 2
nL = 100
N = nL / 2
m = 0:(4 * N)

C1 = 1
C2 = 1
a = 1
u = 1
p = .5
m1 = 2:3

######## select a model ###################
n <- p ^ m        # model1
# n <- p ^ m / m      # model2
# n <- (1/m) ^ 4    # model3
# n <- 1/(2 * m)    # model4
# n<- 1/m      # modified model 4 (HZ)

######## select two latitudes #############

phi <- c(50, 10) * pi / 180

# choose covarince function
mycov <- function(lat1, lat2)
{
  cm = C1 * (C2 - exp(-a * abs(lat1)) - exp(-a * abs(lat2)) + exp(-a * abs(lat1 -lat2)))
  #cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
  
  return(cm)
}

lat1 <- rep(phi, nl)
lat2 = rep(phi, each = nl)
lambda <- 2 * (1:nL - 1) * pi / nL

var.sim <- array(0, c(nn, nL))
mean_sim <- array(0, c(nn, nl))
cor_sim <- array(0, c(nn, 1))
X <- array(0, c(nl, nL, nn)) # the generated data
l1 <- lambda[1:(1 + nL / 2)]
l2 <- c(0, lambda[nL:(1 + nL / 2)])

#For model 1 ============================================================== 
v <- mycov(phi[1], phi[2])
v0 <- 1 * v * (1-p^2) / (1 - 2*p*cos(u * (phi[1] - phi[2]))+p^2)
v1 <- 1 * v * (1-p^2) / (1 - 2*p*cos(l1 + u * (phi[1] - phi[2]))+p^2)
v2 <- 1 * v * (1-p^2) / (1 - 2*p*cos(l2 + u * (phi[1] - phi[2]))+p^2)
var_phi <- 1*(v0 - .5 * (v1 + v2))
#==========================================================================

#For model 2 ============================================================== 
#v <- mycov(phi[1], phi[2])
#v0 <- 1 * v * log(1 / sqrt(1 - 2*p*cos(u * (phi[1] - phi[2]))+p^2))
#v1 <- 1 * v * log(1 / sqrt(1 - 2*p*cos(l1 + u * (phi[1] - phi[2]))+p^2))
#v2 <- 1 * v * log( 1/ sqrt(1 - 2*p*cos(l2 + u * (phi[1] - phi[2]))+p^2))
#var_phi <- 1*(v0 - .5 * (v1 + v2))
#==========================================================================

circulant <- function(x, nrow = length(x)) {
  n <- length(x)
  matrix(x[(1:n - rep(1:nrow, each=n)) %% n + 1L], ncol=n, byrow=TRUE)
}


# The following code is to generate gridded random observations on two latitudes
# we first specify the covariance matrices. According to Li's dissertation, the 
# covariance matrix is a circulant block matrix, where each block is the corresponding
# longitude correlated with other longitudes.

# we first create block matrices, then we form the circulant block matrix
block.list = list()
block.matrix <- array(0, c(2, 2))
for(j in 1:nL){
  # For model 1 =======================================================================================================
    block.matrix[1, 1] = mycov(phi[1], phi[1]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[1] - phi[1]))+p^2)
    block.matrix[1, 2] = mycov(phi[1], phi[2]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[1] - phi[2]))+p^2)
    block.matrix[2, 1] = mycov(phi[2], phi[1]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[2] - phi[1]))+p^2)
    block.matrix[2, 2] = mycov(phi[2], phi[2]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[2] - phi[2]))+p^2)
  #====================================================================================================================  

    # The following code gives a much better result, but it seems the above formulas are correct   
    # For model 1 =====================================================================================================
    block.matrix[1, 1] = mycov(phi[1], phi[2]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[1] - phi[1]))+p^2)
    block.matrix[1, 2] = mycov(phi[1], phi[2]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[1] - phi[2]))+p^2)
    block.matrix[2, 1] = mycov(phi[1], phi[2]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[2] - phi[1]))+p^2)
    block.matrix[2, 2] = mycov(phi[1], phi[2]) * (1-p^2)/(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[2] - phi[2]))+p^2)
    #==================================================================================================================
    
    
    
    # for model 2 =====================================================================================
    # block.matrix[1, 1] = mycov(phi[1], phi[1]) * log(1/sqrt(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[1] - phi[1]))+p^2))
    # block.matrix[1, 2] = mycov(phi[1], phi[2]) * log(1/sqrt(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[1] - phi[2]))+p^2))
    # block.matrix[2, 1] = mycov(phi[2], phi[1]) * log(1/sqrt(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[2] - phi[1]))+p^2))
    # block.matrix[2, 2] = mycov(phi[2], phi[2]) * log(1/sqrt(1 - 2*p*cos((lambda[1]-lambda[j])+u * (phi[2] - phi[2]))+p^2))  
    #================================================================================================
    block.list[[j]] <- block.matrix
}

#Now we form a circulant block matrix
#block.list <- as.list(block.matrix)
#require(Matrix)

bcirc <- function(list.blocks){
  P <- lapply(seq_along(list.blocks), function(x,y) x ==y, x = circulant(seq_along(list.blocks)))
  Reduce('+',Map(P = P, A=list.blocks, f = function(P,A) kronecker(P,A)))
}

R.matrix = bcirc(block.list)  #it is a real symmetric matrix
R.eigen = eigen(R.matrix)$values
R.eigen <- ifelse(R.eigen < 0, 0, R.eigen)
R.vectors = eigen(R.matrix)$vectors

sq.R.matrix <- R.vectors %*% diag(sqrt(R.eigen)) %*% t(R.vectors)
#real.sq.R.matrix <- Re(sq.R.matrix)
#im.sq.R.matrix <- Im(sq.R.matrix)  # should be a zero matrix
hh=1

for(hh in 1:nn){
  set.seed(12345 + hh*42112)
  x <- rnorm(n = nl * nL,
             mean = 0,
             sd = 1)
  
  x.data <- t(sq.R.matrix %*% x)  #*(cos(lat1))^m[kk]
  for (ii in 1:nl) {
      X[ii, , hh] <- x.data[(((ii-1)+1:(2*nL))%%nl) == 1]
  }
  rownames(X) <- round(phi * 180 / pi, 0)
  colnames(X) <- round(lambda * 180 / pi, 0)
  
  for (jj in 1:nL) {
    Y1 = Y2 <- NULL
    Y1 <- array(shift(X[1, , hh], -(jj - 1)) - X[1, , hh], nL)
    Y2 <- array(shift(X[2, , hh], -(jj - 1)) - X[2, , hh], nL)
    var.sim[hh, jj] <- (1 / (2 * nL)) * sum(Y1 * Y2)
  }
  
  
  mean_sim[hh, ] <- apply(X[, , hh] , MARGIN = 1, mean)
  cor_sim[hh, ] <- cor(X[1, , hh], X[2, , hh])
  
}


var_sim <- colMeans(var.sim)

######## select the corresponding theoritical varigram model #############
######## for the above selected model                         #############


l1 <- lambda[1:(1 + nL / 2)]
l2 <- c(0, lambda[nL:(1 + nL / 2)])


#For model 4 (?)======================================================
v <- mycov(phi[1], phi[2])
v0 <- 1 * v * log(1 / (2 * (1 - cos(u * (phi[1] - phi[2])))))
v1 <- 1 * v * log(1 / (2 * (1 - cos(l1 + u * (phi[1] - phi[2])))))
v2 <- 1 * v * log(1 / (2 * (1 - cos(l2 + u * (phi[1] - phi[2])))))
var_phi <- 1*(v0 - .5 * (v1 + v2))
#=====================================================================

#For model 1 ============================================================== 
v <- mycov(phi[1], phi[2])
v0 <- 1 * v * (1-p^2) / (1 - 2*p*cos(u * (phi[1] - phi[2]))+p^2)
v1 <- 1 * v * (1-p^2) / (1 - 2*p*cos(l1 + u * (phi[1] - phi[2]))+p^2)
v2 <- 1 * v * (1-p^2) / (1 - 2*p*cos(l2 + u * (phi[1] - phi[2]))+p^2)
var_phi <- 1*(v0 - .5 * (v1 + v2))
#==========================================================================

#For model 2 ============================================================== 
#v <- mycov(phi[1], phi[2])
#v0 <- 1 * v * log(1 / sqrt(1 - 2*p*cos(u * (phi[1] - phi[2]))+p^2))
#v1 <- 1 * v * log(1 / sqrt(1 - 2*p*cos(l1 + u * (phi[1] - phi[2]))+p^2))
#v2 <- 1 * v * log( 1/ sqrt(1 - 2*p*cos(l2 + u * (phi[1] - phi[2]))+p^2))
#var_phi <- 1*(v0 - .5 * (v1 + v2))
#==========================================================================


plot(lambda[lambda %in% l1], var_sim[lambda %in% l1], type = "o", lty = 2, col = 3,
  xaxt = "n", bty = "n", ylab = "Variance", pch = 19, xlab = expression(paste(Delta, lambda)),
  ylim = c(min(var_phi, var_sim),max(var_sim,var_phi)), 
  main = paste("nl=", nl, ",nL=", nL, ",phi=", phi[1] * 180 / pi, " & ", phi[2] * 180 / pi, ",p=", p,    ", adjustment"),  
  sub = paste("#simulations=", nn))
axis(1, at = lambda, labels = round(lambda * 180 / pi, 0))
lines(l1, var_phi, type = "o", lty = 3, col = 4, lwd = 2, pch = 15)
legend("bottomright", c("simulation", "formula"), pch = c(19, 15), col = 3:4, bty = "n", lwd = 1:2)




