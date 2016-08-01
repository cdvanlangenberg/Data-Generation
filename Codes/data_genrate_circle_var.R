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
nL = 100
N = nL / 2
m = 0:N

C1 = 1
C2 = 1
a = 1

lambda1 <- 2 * (1:nL - 1) * pi / nL
lambda <- ifelse(lambda1 > pi, lambda1 - 2 * pi, lambda1)

cov <- C1 * (exp(-a * abs(lambda)))

circulant <- function(x, nrow = length(x)) {
  n <- length(x)
  matrix(x[(1:n - rep(1:nrow, each = n)) %% n + 1L], ncol = n, byrow = TRUE)
}

Cm <- circulant(cov)
r = round(eigen(Cm)$values, 4)
v = eigen(Cm)$vectors

sqCm <- round(v %*% diag(sqrt(r)) %*% t(v), 8)

var_sim <- array(0, c(nn, nL))

pb <-
  txtProgressBar(min = 0, max = nn, style = 3) 

# cores <- detectCores()
# cl <- makeCluster(cores)
# registerDoParallel(cl)

hh <- as.integer()

foreach(hh = 1:nn) %do% {
setTxtProgressBar(pb, hh)
  set.seed(12345 + hh)
  x <- rnorm(nL, mean = 0, sd = 1)
  # x <- (x - mean(x)) / sd(x)
  X = t(sqCm %*% x)
  
  for (jj in 1:nL) {
    Y <- NULL
    Y <- rbind(X, shift(X, -(jj - 1)))
    var_sim[hh, jj] <-
      (1 / nL) * sum((Y[1, ] - Y[2, ]) ^ 2)  #- mean(Y[1,])*mean(Y[2,])
  }
}

# stopCluster(cl)
close(pb)


var_ave <- colMeans(var_sim) # this is the varince


var = 2 * (1 - exp(-a * abs(lambda)))

# pdf("Results/variogram_plot_4000.pdf",height = 5, width = 7)

plot(
  lambda,
  var,
  xaxt = "n",
  ylim = c(min(var_ave), max(var_ave)),
  ylab = "Variance",
  bty = "n",
  pch = 19,
  main = paste("Variogram; C1=", C1, "nL=", nL),
  sub = paste("# of simulations = ", nn)
)
lines(
  lambda,
  var_ave,
  col = "red",
  bty = "n",
  xaxt = "n",
  type = "p",
  pch = 17
)
axis(1, at = lambda, labels = lambda * 180 / pi)
legend(
  "bottomright",
  legend = c("formula", "simulation"),
  col = 1:2,
  pch = c(19, 17),
  bty = "n"
)

# dev.off()