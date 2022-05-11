library(spatstat)
library(sp)
library(rgeos)
library(aspace)

urm = read.csv("FILE_PATH")


par(mfrow=c(1,2))
plot(urm$xcoord, urm$ycoord, pch=16, col='orange', main = "URM")

#convert URM to ppp object
URMppp <- ppp(urm$xcoord, urm$ycoord,range(urm$xcoord),range(urm$ycoord))
#G function
urmG<-envelope(URMppp,fun=Gest)
par(mfrow=c(1,1))
plot(urmG, main = "G URM")

urmF <- envelope(URMppp, fun=Fest)
plot(urmF, main = "F URM")

urmK <- envelope(URMppp, fun=Kest)
plot(urmK, main = "K URM")

plot(urm$xcoord, urm$ycoord, pch=16, col='orange', main = "URM")

#subset and run clusters
urmXY = urm[,56:57]
urmClust = kmeans(urmXY, 8)

#plot clusters
plot(urmXY[,1],urmXY[,2],col=urmClust$cluster)
points(urmClust$centers[,1],urmClust$centers[,2],col=c(1,2,3,4,5,6,7,8),pch=3,cex=1.5,lwd=2)

#subset clusters
urmc1=which(urmClust$cluster==1)
urmc2=which(urmClust$cluster==2)
urmc3=which(urmClust$cluster==3)
urmc4=which(urmClust$cluster==4)
urmc5=which(urmClust$cluster==5)
urmc6=which(urmClust$cluster==6)
urmc7=which(urmClust$cluster==7)
urmc8=which(urmClust$cluster==8)

#make some ellipses, make sure to install calc_sde2
urmE1=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc1,]),verbose=FALSE)

urmE2=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc2,]),verbose=FALSE)

urmE3=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc3,]),verbose=FALSE)

urmE4=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc4,]),verbose=FALSE)

urmE5=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc5,]),verbose=FALSE)

urmE6=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc6,]),verbose=FALSE)

urmE7=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc7,]),verbose=FALSE)

urmE8=calc_sde2(id=1,filename="elipse_temp.txt",centre.xy=NULL,calccentre=TRUE,weighted=FALSE,
weights=NULL,points=as.matrix(urmXY[urmc8,]),verbose=FALSE)

#draw ellipses
polygon(urmE1[,2],urmE1[,3],border=1)
polygon(urmE2[,2],urmE2[,3],border=2)
polygon(urmE3[,2],urmE3[,3],border=3)
polygon(urmE4[,2],urmE4[,3],border=4)
polygon(urmE5[,2],urmE5[,3],border=5)
polygon(urmE6[,2],urmE6[,3],border=6)
polygon(urmE7[,2],urmE7[,3],border=7)
polygon(urmE8[,2],urmE8[,3],border=8)

#read in URM neighborhood data
urmNabes = read.csv("/Users/Zack/Desktop/Pratt/Classes/SpatialStats/PDX_neighborhood_urm.csv")
#drop county column
urmNabes = urmNabes[-c(5)]
#make correlation matrix
image(cor(urmNabes),axes=FALSE)
axis(1, at=seq(0,1,length=39), labels=names(urmNabes), las = 3)
axis(2, at=seq(0,1,length=39), labels=names(urmNabes), las = 1)

# run some cor tests
cor.test(urmNabes$URM_ratio,urmNabes$DisplacedPop)
cor.test(urmNabes$URM_ratio,urmNabes$CasDayTotal)
cor.test(urmNabes$URM_ratio,urmNabes$BldgLoss)
cor.test(urmNabes$URM_ratio,urmNabes$Debris)

# make a linear model for URM_ratio for all variables
m1 <- lm(urmNabes$URM_ratio ~ urmNabes$DisplacedPop+urmNabes$CasDayTotal+urmNabes$CasNightTotal+
urmNabes$Debris+urmNabes$BldgLoss+urmNabes$PDsExtensive+urmNabes$PDsComplete+
urmNabes$CasNightL2+urmNabes$CasNightL3+urmNabes$CasNightL4)

m2 <- lm(urmNabes$URM_ratio ~ urmNabes$DisplacedPop+urmNabes$Debris+urmNabes$BldgLoss+urmNabes$PDsExtensive+
urmNabes$PDsComplete+urmNabes$normDisPop+urmNabes$CasNightL3+urmNabes$CasNightL4)

m3 <- lm(urmNabes$URM_ratio ~ urmNabes$DisplacedPop+urmNabes$Debris+urmNabes$BldgLoss+urmNabes$PDsExtensive+
urmNabes$PDsComplete+urmNabes$CasNightL3+urmNabes$CasNightL4)

#check for spatial autocorrelation
library(rgdal)
library('ncf')
library('gstat')

nabes = read.csv("/Users/Zack/Desktop/Pratt/Classes/SpatialStats/PDX_neighborhood_units_cent.csv")

#better colors
cgen2<-function(v){
n = (v-min(v))/(max(v)-min(v))
dd=cut(seq(.01,1,by=.01),breaks=100,dig.lab=2,labels=1:100)
dd2=cut(n,breaks=100,dig.lab=2,labels=1:100)
rbPal <- colorRampPalette(c("darkblue", "blue", "cyan"))(100)
ncol=c()
for (i in 1:length(dd2))
{
ncol[i]=rbPal[which(dd == dd2[i])]
}
return(ncol)
}

plot(nabes$xcoord,nabes$ycoord)

cols1=cgen2(nabes$URM_ratio)
cols2=cgen2(nabes$DisplacedP)
#compare variables
par(mfrow=c(1,2))
plot(nabes$xcoord,nabes$ycoord,main="URM_ratio",pch=16,col=cols1)
plot(nabes$xcoord,nabes$ycoord,main="Displaced Pop",pch=16,col=cols2)

#plot morans I autocorrelation
corrDisp <- correlog(x = nabes$xcoord, y = nabes$ycoord, z = nabes$DisplacedP,
increment = 5, resamp = 100, quiet = TRUE)

plot(corrDisp)

corrURM <- correlog(x = nabes$xcoord, y = nabes$ycoord, z = nabes$URM_ratio,
increment = 5, resamp = 100, quiet = TRUE)

plot(corrURM)

par(mfrow=c(1,2))
plot(corrDisp,main="Displaced People")
plot(corrURM,main="URM Ratio")

#make spatial points object for variogram
xyNabes=nabes[,40:41]
names(xyNabes)<-c("x","y")
coordinates(xyNabes)= ~x+y

varioDisp <- variogram(nabes$DisplacedP ~ 1, xyNabes, cloud = TRUE)
plot(varioDisp)

varioURM <- variogram(nabes$URM_ratio ~ 1, xyNabes, cloud = TRUE)
plot(varioURM)

#binned Variograms
varioDisp2 <- variogram(nabes$DisplacedP ~ 1, xyNabes)
plot(varioDisp2, pch = 16, cex = 1.5)

varioUrm2 <- variogram(nabes$URM_ratio ~ 1, xyNabes)
plot(varioUrm2, pch = 16, cex = 1.5, main="URM")

# Directional Variograms
varDispDir<-variogram(nabes$DisplacedP ~ 1, xyNabes,alpha=c(0,45,90,135))
plot(varDispDir, cex = 2.5, pch = 16)

varURMDir<-variogram(nabes$URM_ratio ~ 1, xyNabes,alpha=c(0,45,90,135))
plot(varURMDir, cex = 2.5, pch = 16)

#spatial regression
library(sp)
library(rgeos)
library(rgdal)
library(maptools)
library(spdep)

NabesSHP <-  readOGR("/Users/Zack/Desktop/Pratt/Classes/SpatialStats/PDX_neighborhood_units_comp.shp","PDX_neighborhood_units_comp")
IDs<-row.names(as(NabesSHP,"data.frame"))
NSHPxy<-coordinates(NabesSHP)

#k nearest neighbor
NN1<-knn2nb(knearneigh(NSHPxy,k=1),row.names=IDs)
NN2<-knn2nb(knearneigh(NSHPxy,k=2),row.names=IDs)

par(mfrow=c(1,2))
plot(NabesSHP,main="First Nearest Neighbor")
plot(NN1,NSHPxy,add=TRUE)
plot(NabesSHP,main="First Two Nearest Neighbors")
plot(NN2,NSHPxy,add=TRUE)

#queen Neighbors
nabesQueenKNN<-poly2nb(NabesSHP)
par(mfrow=c(1,1))
plot(NabesSHP,main="Queen")
plot(nabesQueenKNN,NSHPxy,add=TRUE)

#spatial weight matrices
#binary
W1ContBin_nabes<-nb2listw(nabesQueenKNN,style="B",zero.policy=TRUE)
print.listw(W1ContBin_nabes,zero.policy=TRUE)
#row standardized
W2ContRS_nabes<-nb2listw(nabesQueenKNN, zero.policy=TRUE)
print.listw(W2ContRS_nabes,zero.policy=TRUE)

#test spatial dependence
NabesSHP2 = NabesSHP
moran.test(NabesSHP2$DisplacedP,listw=W2ContRS_nabes,zero.policy=TRUE)
moran.plot(NabesSHP2$DisplacedP,listw=W2ContRS_nabes,zero.policy=TRUE)

#OLS
#make sure everything is good to go
NabesSHP2@data$DisplacedP<-as.numeric(as.character(NabesSHP2$DisplacedP))
NabesSHP2@data$Debris<-as.numeric(as.character(NabesSHP2$Debris))
NabesSHP2@data$BldgLoss<-as.numeric(as.character(NabesSHP2$BldgLoss))
NabesSHP2@data$PDsExtensi<-as.numeric(as.character(NabesSHP2$PDsExtensi))
NabesSHP2@data$PDsComplet<-as.numeric(as.character(NabesSHP2$PDsComplet))
NabesSHP2@data$CasNightL3<-as.numeric(as.character(NabesSHP2$CasNightL3))
NabesSHP2@data$CasNightL4<-as.numeric(as.character(NabesSHP2$CasNightL4))
NabesSHP2@data$URM_ratio<-as.numeric(as.character(NabesSHP2$URM_ratio))

#model
MS1<-lm(NabesSHP2$URM_ratio~NabesSHP2$DisplacedP+NabesSHP2$Debris+NabesSHP2$BldgLoss+NabesSHP2$PDsExtensi+
  NabesSHP2$PDsComplet+NabesSHP2$CasNightL3+NabesSHP2$CasNightL4)

#Refine model
MS2<-lm(NabesSHP2$URM_ratio~NabesSHP2$DisplacedP+NabesSHP2$BldgLoss+
  NabesSHP2$PDsComplet+NabesSHP2$CasNightL3+NabesSHP2$CasNightL4)

#morans test
lm.morantest(MS2, W2ContRS_nabes,zero.policy=TRUE)

#testing for alternate models
lm.LMtests(MS2, W2ContRS_nabes, test=c("LMlag", "LMerr"),zero.policy=TRUE)

#trying the lmerr model
Merr<-errorsarlm(NabesSHP2$URM_ratio~NabesSHP2$DisplacedP+NabesSHP2$BldgLoss+
NabesSHP2$PDsComplet+NabesSHP2$CasNightL3+NabesSHP2$CasNightL4,data=NabesSHP2,W2ContRS_nabes,zero.policy=TRUE)

#LMLag model
Mlag<-lagsarlm(NabesSHP2$URM_ratio~NabesSHP2$DisplacedP+NabesSHP2$BldgLoss+
NabesSHP2$PDsComplet+NabesSHP2$CasNightL3+NabesSHP2$CasNightL4,data=NabesSHP2,W2ContRS_nabes,zero.policy=TRUE)

#Breusch-Pagan
bptest.sarlm(Mlag)
bptest.sarlm(Merr)
