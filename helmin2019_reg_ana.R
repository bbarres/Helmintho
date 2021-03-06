##############################################################################/
##############################################################################/
#Analyzing Helmintosporiose 2019 bioassay results from Arvalis field trial
##############################################################################/
##############################################################################/

#loading the libraries
library(drc)
library(plotrix)
library(gdata)
library(tidyr)
library(RColorBrewer)
library(ade4)
library(factoextra)

#load the global dataset
datamyc<-read.table("data/helmindata2019.txt",header=T,sep=";",
                    stringsAsFactors=TRUE)


##############################################################################/
#Regression analysis of mycelial growth experiment
##############################################################################/

#first we extract the list of the different SA listed in the file
SAlist<-levels(datamyc$pest_sa_id)
CompRez<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=character(),SERR=character())

#we make a subselection of the data according to the SA
pdf(file="output/plot_ASA2019.pdf",width=7)
for (j in 1:length(SAlist)) {
  data_subSA<-datamyc[datamyc$pest_sa_id==SAlist[j],]
  data_subSA$ech_id<-drop.levels(data_subSA$ech_id)
  
  REZSA<-data.frame(Subs_Act=factor(),sample_ID=factor(),
                    ED50=character(),SERR=character())
  
  for (i in 1:dim(table(data_subSA$ech_id))[1]) {
    tempdat<-data_subSA[data_subSA$ech_id==names(table(data_subSA$ech_id))[i],]
    if(tempdat[tempdat$dose==max(tempdat$dose),"rslt_03"]>30) {
      tempx<-data.frame("Subs_Act"=SAlist[j],"sample_ID"=tempdat$ech_id[1],
                        "ED50"=paste(">",max(tempdat$dose),sep=""),
                        "SERR"="NA")
    } else {
      temp.m1<-drm(rslt_03~dose,
                   data=tempdat,
                   fct=LN.3())
      plot(temp.m1,ylim=c(0,110),xlim=c(0,100),
           main=paste(SAlist[j],as.character(tempdat$ech_id[1])))
      temp<-ED(temp.m1,c(50),type="absolute",interval="delta")
      tempx<-data.frame("Subs_Act"=SAlist[j],
                        "sample_ID"=as.character(tempdat$ech_id[1]),
                        "ED50"=as.character(temp[1]),
                        "SERR"=as.character(temp[2]))
    }
    
    REZSA<-rbind(REZSA,tempx)
  }
  CompRez<-rbind(CompRez,REZSA)
}
dev.off()

#exporting the result as a text file
CompRez<-CompRez[order(CompRez$Subs_Act,CompRez$sample_ID),]
write.table(CompRez,file="output/2019_results.txt",
            sep="\t",quote=FALSE,row.names=FALSE)


##############################################################################/
#barplot to compare the ED50 of the different samples####
##############################################################################/

#just a small graphic to gain insight on the first round of results
#first, we replace the ED50 that were too high to be evaluated with an 
#arbitrary value
cooloor<- brewer.pal(12,"Set3")
CompRez$ED50<-as.character(CompRez$ED50)
CompRez[CompRez$ED50==">30","ED50"]<-35
CompRez$ED50<-as.numeric(as.character(CompRez$ED50))

pdf(file="output/histo_AllInd_ASA2019.pdf",width=60,height=8)
op<-par(mfrow=c(1,1))
par(mar=c(8,3,3,0.5))
barplot(as.numeric(as.character(CompRez$ED50)),
        ylim=c(0,35),
        col=cooloor[as.numeric(as.factor(CompRez$Subs_Act))],
        names.arg=CompRez$sample_ID,las=2,
        main="Comparison of the different samples by SA")
abline(h=35,lty=2)
legend(300,47,levels(CompRez$Subs_Act),fill=cooloor,bty="n")
par(op)
dev.off()

#histograms by samples
samplelist<-as.character(names(table(CompRez$sample_ID)))
pdf(file="output/histo_byInd_ASA2019.pdf",width=9,height=20)
op<-par(mfrow=c(8,5))
for (i in (1:length(samplelist))) {
  temp<-merge(as.data.frame(levels(as.factor(CompRez$Subs_Act))),
              CompRez[CompRez$sample_ID==samplelist[i],],
              by.x=1,by.y=1,
              all.x=TRUE)
  barplot(temp$ED50,col=cooloor,las=1,main=samplelist[i],
          ylim=c(0,52))
}
par(op)
dev.off()
#export to pdf 12 x 16


##############################################################################/
#correlation between ED50 estimated for different active substances####
##############################################################################/

temp<-CompRez[,c(1,2,3)]
temp<-spread(temp,Subs_Act,ED50)

#a function to compute the absolute correlation between pairs of variables
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y,use="pairwise.complete.obs"))
  txt <- format(c(r, 2), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

pairs(temp[,c(2:5)],las=1,main="Correlation between ActSubst",
      lower.panel=panel.smooth, upper.panel=panel.cor)

pairs(log(temp[,c(2:5)]),las=1,main="Correlation between log(ActSubst)",
      lower.panel=panel.smooth, upper.panel=panel.cor)
#export to pdf 11 x 11 inches


##############################################################################/
#Analyzing the multisensitivity profile of the strains####
##############################################################################/

#Clusterization based on scoring of 4 SA
row.names(temp)<-temp$sample_ID

#PCA for the scoring on 4 SA
truc<-dudi.pca(temp[,-c(1)],
               scannf=FALSE,nf=3)
scatter(truc)
#determining the optimal number of clusters
fviz_nbclust(temp[,c(2:5)],kmeans,method="gap_stat")
clust<-kmeans(temp[,c(2:5)],3)
fviz_cluster(clust,data=temp[,c(2:5)])
plot(truc$li[,c(1,2)],col=brewer.pal(5,"Dark2")[clust$cluster],
     pch=19,cex=2)

hclu<-hclust(dist(scale(temp[,c(2:5)]),
                  method="euclidean"),
             method="ward.D2")
plot(hclu)
fviz_dend(hclu,k=3,cex=0.5,rect=TRUE,
          k_colors=brewer.pal(5,"Dark2"))
#export to pdf 9 x 6 inches


##############################################################################/
#plot of the distribution of IC50 for each active substance####
##############################################################################/

#preparing the dataset
temp<-CompRez[,c(1,2,3)]
temp<-spread(temp,Subs_Act,ED50)

#distribution of the IC50 by Active Substance
op<-par(mfrow=c(2,2))
plot(temp[order(c(temp$BIXAFENE)),"BIXAFENE"],
     main="BIXAFENE IC50",bg=cooloor[1],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$BOSCALID)),"BOSCALID"],
     main="BOSCALID IC50",bg=cooloor[2],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$FLUOPYRAM)),"FLUOPYRAM"],
     main="FLUOPYRAM IC50",bg=cooloor[5],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
plot(temp[order(c(temp$FLUXAPYROXAD)),"FLUXAPYROXAD"],
     main="FLUXAPYROXAD IC50",bg=cooloor[3],pch=21,cex=2,las=1,
     ylab="IC50",ylim=c(0,52))
par(op)
#export to pdf 10 x 8 inches


##############################################################################/
#END
##############################################################################/