####### this code will generate the grid plot

library(geoR)
library(lattice)
library(akima)
library(fields)

#library(RandomFields)
#library(scatterplot3d)


DD<- readRDS("Data/Data_120X180_model1.rds") ## give data loaction

ph <- as.numeric(rownames(DD[,,1]))
lam <- as.numeric(colnames(DD[,,1]))

#pdf("Data_sample_120_model3.pdf", width = 8, height = 8)
par(mfrow=c(2,2), mar= c(3,3,2,0.4), las=1, mgp=c(2,1,0))

for(ii in c(1,3,4,7)){
data <- DD[,,ii]
data <- data.frame(res=matrix(data, ncol=1))
data$lati<-rep(ph,nL)
data$long<-rep(lam,each=nl)
data <-data[order(-data$lati, data$long),c(2,3,1)]

data.geo<-as.geodata(data)

long = data.geo$coords[,2]
lati = data.geo$coords[,1]
data.interp<-interp(y=lati,x=long,z=data$res)
plot.surface(data.interp, type="I", main = " ", zlim=c(-10,10), xlab="longitude", ylab=" " ) #, yaxt="n"
}

mtext(text ="Axially symmetric data", side = 3, outer = T, line = -1)

par(mfrow=c(1,1))

dev.off()

