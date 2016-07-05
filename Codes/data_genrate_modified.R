nn = 200 # number of iterations

nl = 4;
nL = 36;
N = nL / 2;
m = 0:(1 * N);
C1 = 1;
C2 = 1;
a = 1;
u = 1;
p = .5;
m1 = 2:3

###### generate values for phi fom -pi/2 to pi/2
# lim_phi = 2*pi/3 # boundy for phi
# 
# phi <- seq( -lim_phi/2, lim_phi/2, by=lim_phi/(nl+0))[-c(1,nl+(2))]
# if(nl%%2==1){
#   phi <- seq( -lim_phi/2, lim_phi/2, by=lim_phi/(nl+1))[-c(1,nl+(2))]
# }

phi <- c(10, 20, 30, 70) * pi / 180

if (length(phi) != nl)
	stop("incorrect number of latitudes \n check: nl and phi")

mycov <- function(lat1, lat2) {
	cm = C1 * (C2 - exp( - a * abs(lat1)) - exp( - a * abs(lat2)) + exp( - a * abs(lat1 - lat2)))
	#cm = C1*(C2-(1/sqrt(a^2+lat1^2))-(1/sqrt(a^2+lat2^2))+(1/sqrt(a^2+(lat1-lat2)^2)))
	;
	cm
}

nvar <- function(x) {
	(1 / (length(x))) * sum(x ^ 2)
}
# variance function, without mean  

lat1 <- rep(phi, nl)
lat2 = rep(phi, each = nl)
lambda <- 2 * (1:nL - 1) * pi / nL

Cm = mycov(lat1, lat2)

############################## select a model ##################### 
n <- p ^ m # model1
var_phi <- C1 * 6 * (1 - exp( - abs(phi))) #model 1
#
# var_phi <- 4*log(2)*(1-exp(-abs(phi))) #model2
# n <- p^m/m      # model2
#
# n <- (1/m)^4    # model3
# var_phi <- (2/45)*(pi^4)*(1-exp(-abs(phi))) #model3
###################################################################

########## Initialize the vectors/matirces 

X <- array(0, c(nl, nL, nn)) # the generated data
var_sim <- array(0, c(nn, nl))
mean_sim <- array(0, c(nn, nl))
cor_sim <- array(0, c(nn, 1))

for (hh in 1:nn) {

	RCm <- array(0, c(rep(length(phi), 2), length(m)))
	ICm <- array(0, c(rep(length(phi), 2), length(m)))
	Kv <- array(0, c(rep(2 * length(phi), 2), length(m)))
	d <- array(0, c(length(m), 2 * length(phi))) # this is w_m real and complex
	r <- array(0, c(length(m), 2 * length(phi)))

	for (kk in 1:length(m)) {

		RCm[,, kk] <- t(array(Cm * cos(m[kk] * u * (lat1 - lat2)) * (1 * n[kk]), dim = rep(length(phi), 2)))
		ICm[,, kk] <- t(array(Cm * sin(m[kk] * u * (lat1 - lat2)) * (1 * n[kk]), dim = rep(length(phi), 2)))

		if (n == Inf && m[kk] %in% c(0)) {
			RCm[,, kk] = ICm[,, kk] = 0
		}
		# make this adjustment when C0 = 0

		x1 <- cbind(RCm[,, kk], - ICm[,, kk])
		x2 <- cbind(ICm[,, kk], RCm[,, kk])
		Kv[,, kk] <- .5 * (rbind(x1, x2))
		colnames(d) <- c(paste0("R", 1:nl), paste0("I", 1:nl))

		if (m[kk] %in% c(0, N)) {
			Kv[,, kk] = 0
			Kv[1:nl, 1:nl, kk] = RCm[,, kk]
		}

		r[kk,] = eigen(Kv[,, kk])$values
		r[kk,] <- ifelse(r[kk,] < 0, 0, r[kk,])
		v = eigen(Kv[,, kk])$vectors

		sqKv <- v %*% diag(sqrt(r[kk,])) %*% t(v)

		set.seed(12345 + hh)
		x <- rnorm(n = 2 * nl, mean = 0, sd = 1)

		d[kk,] <- t(sqKv %*% x) #*(cos(lat1))^m[kk]
		rownames(d) <- m
		colnames(d) <- c(paste0("R", 1:nl), paste0("I", 1:nl))
	}


	for (ii in 1:nl) {
		for (jj in 1:nL) {
			X[ii, jj, hh] <- (d[1, ii] + 2 * sum(d[2:length(m), ii] * round(cos(m[-1] * lambda[jj]), 4)
			- d[2:length(m), ii + nl] * round(sin(m[-1] * lambda[jj]), 4))) #+ m1[ii]  #*cos(phi[ii])
		}
	}

	rownames(X) <- round(phi * 180 / pi, 0)
	colnames(X) <- round(lambda * 180 / pi, 0)

	var_sim[hh,] <- apply(X[,, hh], MARGIN = 1, nvar)
	mean_sim[hh,] <- apply(X[,, hh], MARGIN = 1, mean)
}

var_sim <- colMeans(var_sim)


plot(phi, var_phi, type = "l", lty = 1, col = 3, ylim = c(min(var_sim, var_phi), max(var_sim, var_phi)), xaxt = "n", bty = "n", ylab = "Variance",
	 xlab = expression(phi), main = paste("nl=", nl, "nL=", nL, "C1=", C1, "C2=", C2, "p=", p), sub = paste("# of simulations=", nn))
axis(1, at = phi, labels = round(phi * 180 / pi, 0))
lines(phi, var_sim, type = "l", lty = 2, col = 4, lwd = 2)
legend("top", c("Theoretical", "Emperical"), lty = c(2, 3), col = 3:4, bty = "n", lwd = 1:2)








