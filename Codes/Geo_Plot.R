
library(geoR)
library(RandomFields)
library(akima)
library(scatterplot3d)
library(fields)
library(lattice)

DD<- readRDS("C:/Users/cdvanlan/Dropbox/Research/Spatial/Simulations/DataGeneration/Data/Data_120X180_model3.rds")
ph <- as.numeric(rownames(DD[,,1]))
lam <- as.numeric(colnames(DD[,,1]))

pdf("Data_sample_120_model3.pdf", width = 8, height = 8)
par(mfrow=c(2,2), mar= c(3,3,2,0.4), las=1, mgp=c(2,1,0))

for(ii in c(1,3,4,7)){
#data <- scale(X[,,ii])
data <- DD[,,ii]
data <- data.frame(res=matrix(data, ncol=1))
#data <- data.frame(res=rnorm(nl*nL))
data$lati<-rep(ph,nL)
data$long<-rep(lam,each=nl)
data <-data[order(-data$lati, data$long),c(2,3,1)]

# data1<-read.csv("Data/axial_data_51X125.csv", header=T)
# data <- (apply(X[,,], MARGIN = c(1,2), FUN=mean)) # get averages of the genetrated data
# rownames(data) <-NULL
# colnames(data) <-NULL

data.geo<-as.geodata(data)

long = data.geo$coords[,2]
lati = data.geo$coords[,1]
data.interp<-interp(y=lati,x=long,z=data$res)
plot.surface(data.interp, type="I", main = " ", zlim=c(-10,10), xlab="longitude", ylab=" " ) #, yaxt="n"
#axis(2, at=seq(-60, 60, by=20), labels = seq(30, 150, by=20) )
}

mtext(text ="Axially symmetric data", side = 3, outer = T, line = -1)

par(mfrow=c(1,1))

dev.off()

