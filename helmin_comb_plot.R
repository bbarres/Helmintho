###############################################################################
###############################################################################
#Comparison of the different level of resistance by strains
###############################################################################
###############################################################################

#loading the libraries
library(dplyr)
#loading the dataset of the combined results of the 4 SDHI resistance test
EC50_pop<-read.table("data/EC50_byPOP.txt",header=TRUE)


###############################################################################
#plot with individual group by pop
###############################################################################

op<-par(mfrow=c(2,2),mar=c(2,2.5,3,1))
EC50bosc<-EC50_pop[EC50_pop$SA_ID=="boscalid",]
plot(EC50bosc$ED50,col=as.character(EC50bosc$pop_col),main="Boscalid")
abline(0.39,0,col="green3",lwd=2)
abline(3.9,0,col="orange3",lwd=2)

EC50bixa<-EC50_pop[EC50_pop$SA_ID=="bixafen",]
plot(EC50bixa$ED50,col=as.character(EC50bixa$pop_col),main="Bixafen")
abline(0.08,0,col="green3",lwd=2)
abline(0.8,0,col="orange3",lwd=2)

EC50fluo<-EC50_pop[EC50_pop$SA_ID=="fluopyram",]
plot(EC50fluo$ED50,col=as.character(EC50fluo$pop_col),main="Fluopyram")
abline(0.44,0,col="green3",lwd=2)
abline(4.4,0,col="orange3",lwd=2)

EC50flux<-EC50_pop[EC50_pop$SA_ID=="fluxapyroxade",]
plot(EC50flux$ED50,col=as.character(EC50flux$pop_col),main="Fluxapyroxade")
abline(0.21,0,col="green3",lwd=2)
abline(2.1,0,col="orange3",lwd=2)

par(op)


###############################################################################
#plot with individual group by pop with FR
###############################################################################

op<-par(mfrow=c(2,2),mar=c(2,2.5,3,1))

plot(EC50bosc$ED50[order(EC50bosc$ED50)]/0.39,main="Boscalid",
     xlab="Souches ID",ylab="FR",las=1)
abline(0.39/0.39,0,col="green4",lwd=2)
abline(3.9/0.39,0,col="red",lwd=2)

plot(EC50bixa$ED50[order(EC50bixa$ED50)]/0.08,main="Bixafen",
     xlab="Souches ID",ylab="FR",las=1)
abline(0.08/0.08,0,col="green4",lwd=2)
abline(0.8/0.08,0,col="red",lwd=2)

plot(EC50fluo$ED50[order(EC50fluo$ED50)]/0.44,main="Fluopyram",
     xlab="Souches ID",ylab="FR",las=1)
abline(0.44/0.44,0,col="green4",lwd=2)
abline(4.4/0.44,0,col="red",lwd=2)

plot(EC50flux$ED50[order(EC50flux$ED50)]/0.21,main="Fluxapyroxade",
     xlab="Souches ID",ylab="FR",las=1)
abline(0.21/0.21,0,col="green4",lwd=2)
abline(2.1/0.21,0,col="red",lwd=2)

par(op)


###############################################################################
#END
###############################################################################