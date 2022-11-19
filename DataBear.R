library(ggplot2)
library(gganimate)
library(gifski)
library(geosphere)
library(ggOceanMaps)

BearGPS <- read.csv("PB_GPS.csv")
BearARGOS <-read.csv("PB_Argos.csv")
# creating a data.frame with your lat/lon points

#png("lol.png")
df <- as.data.frame(cbind(Lon = BearGPS$Lon,Lat =BearGPS$Lat))
bbox <-c(min(df$Lon)-15,max(df$Lon)+15,min(df$Lat),max(df$Lat))


#basemap(limits = c(-180,180,50,80)) + 
#  geom_spatial_path(data = df[1:20,], aes(x = Lon, y = Lat),crs = 4326, color = "blue")+
#  geom_spatial_point(data = df[1,], aes(x = Lon, y = Lat),crs=4326, color = "red")
#dev.off() 

BearGPS$MotionAngle <- bearing(cbind(BearGPS$Lon,BearGPS$Lat))*pi/180
BearARGOS$MotionAngle <-bearing(cbind(BearARGOS$Lon,BearARGOS$Lat))*pi/180
BearARGOS$distancekm <- array(0,dim=length(BearARGOS$Lon))
BearGPS$distancekm <- array(0,dim=length(BearGPS$Lon))

for (i in 2:length(BearGPS$Lon)){
  BearGPS$distancekm[i]<-distm(c(BearGPS$Lon[i-1], BearGPS$Lat[i-1]), c(BearGPS$Lon[i], BearGPS$Lat[i]), fun = distHaversine)/1000
}
for (i in 2:length(BearARGOS$Lon)){
  BearARGOS$distancekm[i]<-distm(c(BearARGOS$Lon[i-1], BearARGOS$Lat[i-1]), c(BearARGOS$Lon[i], BearARGOS$Lat[i]), fun = distHaversine)/1000
}
#Turning angle
BearGPS$Tangle<- rbind(0,matrix(BearGPS$MotionAngle[2:357]-BearGPS$MotionAngle[1:356],ncol=1))
BearARGOS$Tangle <-rbind(0,matrix(BearARGOS$MotionAngle[2:length(BearARGOS$MotionAngle)]-BearARGOS$MotionAngle[1:(length(BearARGOS$MotionAngle)-1)],ncol=1))

BearGPS<-BearGPS[1:356,]
BearGPS$TangleDiff<- rbind(0,matrix(BearGPS$Tangle[2:356]-BearGPS$Tangle[1:355],ncol=1))

#We only look at the turning angle and distance km: they are the hidden states.
VMLogLikelihood <- function(k) -sum(log(besselI(k,0)))+k*sum(cos(BearGPS$TangleDiff))
A1 <- matrix(1,ncol=1,nrow=1)
B1<-  matrix(0,ncol=1,nrow=1)
KMLE <- maxLik(VMLogLikelihood,grad=NULL,hess=NULL,c(1),method='BFGS',list(ineqA=A1, ineqB=B1))
KMLE <-KMLE$estimate

GammaMomentshape <- mean(BearGPS$distancekm)^2/var(BearGPS$distancekm)
GammaMomentscale <- var(BearGPS$distancekm)/mean(BearGPS$distancekm)
