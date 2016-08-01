list.packages <-
  c("plyr", "foreach", "doParallel", "binhf")

new.packages <-
  list.packages[!(list.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

# Load packages into session
sapply(list.packages, require, character.only = TRUE)

nn = 400
nl = 1
nL = 48
N = nL / 2
m = 0:N

C1 = 1
C2 = 1
a = 1


lambda1 <- 2 * (1:nL - 1) * pi / nL
lambda <- ifelse(lambda1 > pi,
                 lambda1 - 2 * pi, lambda1)

cov <- C1 * (exp(-a * abs(lambda)))

## define a new covariance function by substracting a0
# a0 <- (1 / (a * pi)) * (1 - exp(-pi * a))
# cov <- C1 * ((exp(-a * abs(lambda))) - a0)

#cov <- C1*(exp(-a*abs(lambda))^alpha)
#cov1 <-function(l=lambda, b=a){
#   cov=(1-(3*abs(l)/(2*b))+(abs(l)^3/(2*b^3)))
#   cov[abs(l)>b] =0
#   ;cov
# }
# cov <- cov1()

circulant <- function(x, nrow = length(x)) {
  n <- length(x)
  matrix(x[(1:n - rep(1:nrow, each = n)) %% n + 1L], ncol = n, byrow = TRUE)
}

Cm <- circulant(cov, length(cov))
r = round(eigen(Cm)$values, 8)
v = eigen(Cm)$vectors

sqCm <- round(v %*% diag(sqrt(r)) %*% t(v), 8)

cov_sim <- array(0, c(nn, nL))
dd <- rep(0, nn)

pb <-
  txtProgressBar(min = 0, max = nn, style = 3) # this is not usefull if processing parallel

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

hh <- as.integer()

foreach(hh = 1:nn) %do% {
  # sapply(list.packages, require, character.only = TRUE)

  X <- NULL
  set.seed(12345 + hh)
  x <- rnorm(nL, mean = 0, sd = 1)
  
  X <- t(sqCm %*% x)
  #X <- array(scale(X[1,]))
  
  for (jj in 1:nL) {
    Y <- NULL
    Y <- rbind(X, shift(X, -(jj - 1)))
    cov_sim[hh, jj] <-
      (1 / nL) * (sum(Y[1,] * Y[2,]))  - mean(Y[1, ]) * mean(Y[2, ])
    dd[hh] <- mean(X) ^ 2
    #print(mean(Y[1,]))
  }
}

stopCluster(cl)
close(pb)

cov_ave <- colMeans(cov_sim)
#cov_ave

#print(cov_ave)
plot(
  lambda,
  cov,
  xaxt = "n",
  ylim = c(min(cov_ave), max(cov_ave, cov)),
  bty = "n",
  pch = 19
)
lines(
  lambda,
  cov_ave,
  col = "red",
  bty = "n",
  xaxt = "n",
  type = "p",
  pch = 17
)
axis(1, at = lambda, labels = lambda * 180 / pi)
#axis(2, at=seq(-.2,1,by=.2), labels = seq(-.2,1,by=.2))
legend(
  "topright",
  legend = c("formula", "simulation"),
  col = 1:2,
  pch = c(19, 17),
  bty = "n"
)
