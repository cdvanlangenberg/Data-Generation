####### this code will generate the grided data for 2 latitudes X given number of longitudes
####### and compare with cross variogram

# list.packages <- c("plyr", "foreach", "doParallel", "binhf")
list.packages <- c("foreach", "binhf")

new.packages <-
  list.packages[!(list.packages %in% installed.packages()[, "Package"])]
if (length(new.packages))
  install.packages(new.packages)

# Load packages into session
sapply(list.packages, require, character.only = TRUE)

nn = 400 # number of iterations

nl = 2
nL = 100
N = nL / 2
m = 0:(1 * N)

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

######## select two latitudes #############

phi <- c(20, 80) * pi / 180

# choose covarince function
mycov <- function(lat1, lat2)
{
  cm = C1 * (C2 - exp(-a * abs(lat1)) - exp(-a * abs(lat2)) + exp(-a * abs(lat1 - lat2)))
  #cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
  
  cm
}


lat1 <- rep(phi, nl)
lat2 = rep(phi, each = nl)
lambda <- 2 * (1:nL - 1) * pi / nL

Cm = mycov(lat1, lat2)
X <- array(0, c(nl, nL, nn)) # the generated data
var_sim <- matrix(0, nrow = nn, ncol = nL)
mean_sim <- array(0, c(nn, nl))
cor_sim <- array(0, c(nn, 1))

pb <-
  txtProgressBar(min = 0, max = nn, style = 3) # this is not usefull if processing parallel

cores <- detectCores()
cl <- makeCluster(cores)
registerDoParallel(cl)

hh <- as.integer()

foreach(hh = 1:nn) %do% {
  sapply(list.packages, require, character.only = TRUE)
  
  setTxtProgressBar(pb, hh)
  
  # for (hh in 1:nn) {
  RCm <- array(0, c(rep(length(phi), 2), length(m)))
  ICm <- array(0, c(rep(length(phi), 2), length(m)))
  Kv <-  array(0, c(rep(2 * length(phi), 2), length(m)))
  d <-  array(0, c(length(m), 2 * length(phi)))
  r <-  array(0, c(length(m), 2 * length(phi)))
  
  for (kk in 1:length(m)) {
    RCm[, , kk] <-
      t(array(Cm * cos(m[kk] * u * (lat1 - lat2)) * (1 * n[kk]), dim = rep(length(phi), 2)))
    ICm[, , kk] <-
      t(array(Cm * sin(m[kk] * u * (lat1 - lat2)) * (1 * n[kk]), dim = rep(length(phi), 2)))
    
    if (n == Inf &&
        m[kk] %in% c(0)) {
      RCm[, , kk] = ICm[, , kk] = 0
    } # make this adjustment when C0 = 0
    
    x1 <- cbind(RCm[, , kk], -ICm[, , kk])
    x2 <- cbind(ICm[, , kk], RCm[, , kk])
    Kv[, , kk] <- .5 * (rbind(x1, x2))
    colnames(d) <- c(paste0("R", 1:nl), paste0("I", 1:nl))
    
    if (m[kk] %in% c(0)) {
      Kv[, , kk] = 0
      Kv[1:nl, 1:nl, kk] = RCm[, , kk]
    }
    
    r[kk,] = eigen(Kv[, , kk])$values
    r[kk,] <- ifelse(r[kk,] < 0, 0, r[kk,])
    v = eigen(Kv[, , kk])$vectors
    
    sqKv <- v %*% diag(sqrt(r[kk,])) %*% t(v)
    set.seed(12345 + hh * 53221 + kk * 4346)
    
    x <- rnorm(n = 2 * nl,
               mean = 0,
               sd = 1)
    
    d[kk,] <- t(sqKv %*% x)  #*(cos(lat1))^m[kk]
    rownames(d) <- m
    colnames(d) <- c(paste0("R", 1:nl), paste0("I", 1:nl))
  }
  
  
  for (ii in 1:nl) {
    for (jj in 1:nL) {
      X[ii, jj, hh] <-
        (d[1, ii] + 2 * sum(d[2:length(m), ii] * cos(m[-1] * lambda[jj])
                            - d[2:length(m), ii + nl] * sin(m[-1] * lambda[jj])))  + m1[ii]  #*cos(phi[ii])
    }
  }
  
  rownames(X) <- round(phi * 180 / pi, 0)
  colnames(X) <- round(lambda * 180 / pi, 0)
  
  # the emperical variogram formula given by Hans Wackernagel 20.15
  
    for (jj in 1:nL) {
    Y1 = Y2 <- NULL
    Y1 <- array(shift(X[1, , hh], -(jj - 1)) - X[1, , hh], nL)
    Y2 <- array(shift(X[2, , hh], -(jj - 1)) - X[2, , hh], nL)
    var_sim[hh, jj] <- (1 / (2 * nL)) * sum(Y1 * Y2)
  }
  
  mean_sim[hh,] <- apply(X[, , hh] , MARGIN = 1, mean)
  cor_sim[hh,] <- cor(X[1, , hh], X[2, , hh])

invisible(gc())
   
}

stopCluster(cl)
close(pb)

var_sim <- colMeans(var_sim)

######## select the corresponding theoritical varigram model #############
######## for the above selected model                         #############


l1 <- lambda[1:(1 + nL / 2)]
l2 <- c(0, lambda[nL:(1 + nL / 2)])

#For model 1 ==============================================================
v <- mycov(phi[1], phi[2])
v0 <-
  1 * v * (1 - p ^ 2) / (1 - 2 * p * cos(u * (phi[1] - phi[2])) + p ^ 2)
v1 <-
  1 * v * (1 - p ^ 2) / (1 - 2 * p * cos(l1 + u * (phi[1] - phi[2])) + p ^ 2)
v2 <-
  1 * v * (1 - p ^ 2) / (1 - 2 * p * cos(l2 + u * (phi[1] - phi[2])) + p ^ 2)
var_phi <- 1 * (v0 - .5 * (v1 + v2))
#==========================================================================

#For model 4 (?)======================================================
# v <- mycov(phi[1], phi[2])
# v0 <- 1 * v * log(1 / (2 * (1 - cos(u * (phi[1] - phi[2])))))
# v1 <- 1 * v * log(1 / (2 * (1 - cos(l1 + u * (phi[1] - phi[2])))))
# v2 <- 1 * v * log(1 / (2 * (1 - cos(l2 + u * (phi[1] - phi[2])))))
# var_phi <- 1*(v0 - .5 * (v1 + v2))
#==========================================================================


plot(
  lambda[lambda %in% l1],
  var_sim[lambda %in% l1],
  type = "o",
  lty = 2,
  col = 3,
  xaxt = "n",
  bty = "n",
  ylab = "Variance",
  pch = 19,
  xlab = expression(paste(Delta, lambda)),
  ylim = c(min(var_sim, var_phi), max(var_sim, var_phi)),
  main = paste(
    "nl=",
    nl,
    ",nL=",
    nL,
    ",phi=",
    phi[1] * 180 / pi,
    " & ",
    phi[2] * 180 / pi,
    ",p=",
    p,
    ", adjustment"
  ),
  sub = paste("#simulations=", nn)
)
axis(1, at = lambda, labels = round(lambda * 180 / pi, 0))
lines(
  l1,
  var_phi,
  type = "o",
  lty = 3,
  col = 4,
  lwd = 2,
  pch = 15
)
legend(
  "bottomright",
  c("simulation", "formula"),
  pch = c(19, 15),
  col = 3:4,
  bty = "n",
  lwd = 1:2
)
