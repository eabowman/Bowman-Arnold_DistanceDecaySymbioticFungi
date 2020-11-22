library(raster)
library(rgdal)
library(scales)
library(maptools)
library(tidyverse)

sites <- read.csv('data/SiteCoordinates.csv', as.is = T)

# create column indicating sampling, shape for plotting points
EM.only <- c('M4','M5','M6','M7','M8','M9','H4','H5','H6','A1','A2','A3')
sites$shape <- NA
for(i in 1:nrow(sites)){
  if(sites[i, 'site'] %in% EM.only){sites[i, 'shape'] <- 1}
  else{sites[i,'shape'] <- 16}
}

sites %>%
  dplyr::select(site, lat, long, shape) %>%
  distinct(site, .keep_all = T) -> sites

# remove non focal SCM sites
sites <- sites[!sites$site %in% c('Pullout south of Willow Canyon Rd.','Rose Canyon Rd.',
                                 'Marshall Gulch','Green Mountain/Bug Springs','Bear Wallow',
                                 'Parking pullout South of Rose Canyon','Middle Bear'),]

#get USA map from GADM, level=0 no state border, level=1 state border, level=2 conty border...
#US_0<-getData('GADM', country="USA", level=0)
US_1<-getData('GADM', country="USA", level=1)
#US_2<-getData('GADM', country="USA", level=2)
#US_3<-getData('GADM', country="USA", level=3)
#US_4<-getData('GADM', country="USA", level=4)

#check the data class
class(US_1)
#check data fram names
names(US_1)
#list the all Name of states
unique(US_1$NAME_1)
#get the Arizona state border map from the whole USA map
Arizona<-subset(US_1, NAME_1=="Arizona")
plot(Arizona)

#read the downloaded WorldClm data (tile12) by raster
AnnPrec <- raster('~/../../../Volumes/THALICTRUM/Worldclim_Arizona/wc2.0_30s_bio/wc2.0_bio_30s_12.tif')
AnnMeanTemp <- raster('~/../../../Volumes/THALICTRUM/Worldclim_Arizona/wc2.0_30s_bio/wc2.0_bio_30s_01.tif')
PrecDriest <- raster('~/../../../Volumes/THALICTRUM/Worldclim_Arizona/wc2.0_30s_bio/wc2.0_bio_30s_17.tif')
#plot(AnnPrec)

# or can directly download the whole bio data using the function below:
#var can be bio1 to bio12, bio is the whole climate variables data set
#res = resolutions of the data
#lon, lat = specify a point then it will download the data tile in that area
Bio<-getData("worldclim", var="bio", res=0.5, lon=-110, lat=33)
#crop the cilmate data with specific region, in this case it's Arizona map we get from US_1
Bio.sub<-crop(Bio, extent(Arizona))
#mask function make sure the climate data only within the Arizona border
Bio.sub.mask<-mask(Bio.sub, Arizona)

#<< MAT >> -------------------
pdf(file = paste0(fig.dir, 'MAT_map.pdf'),
    width = 12.5, height = 8.5)

Temp.col<-colorRampPalette(c("blue","deepskyblue","cyan", "yellow","orange", "darkorange", "red"))
plot(Bio.sub.mask$bio1_12/10, main="Annual Mean Temperature", xlab= "Longitude", ylab="Latitude", axes=FALSE, box=FALSE, col= Temp.col(9))
plot(Arizona,add=TRUE)
axis(1)
axis(2)

#add all points
plot(Bio.sub.mask$bio1_12/10, main="Annual Mean Temperature",legend=FALSE,axes=FALSE, box=FALSE, col=Temp.col(9) ,xlab= "Longitude", ylab="Latitude")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=1, tck=-0.02)
axis(2, las=2, cex.axis=1, tck=-0.02)
points(sites$long, sites$lat, pch=sites$shape, col='black', cex=4)
#pointLabel(sites$Longitude, sites$Latitude, label=rownames(sites), cex=0.1)
#text(sites$Longitude, sites$Latitude, label=rownames(sites), cex=0.7) 
#site.drop.overlaps<-sites[-c(3,7,10,11,14,17,18,21,22,23,24,25,26,27,28,31,33,34,37,40,42,45,46),]

par(xpd=TRUE)
#legend(-108.6,37.1, legend=levels(sites$Area.code),bty="n", pch=19,col=color.vec, pt.cex=1.6, cex=1)
#temperature legend
Bio.sub.mask.range<-c(0,20)
plot(Bio.sub.mask$bio1_12/10, legend.only=TRUE,col=Temp.col(9),
     legend.width=1, legend.shrink=0.5, 
     axis.args=list(at=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2], 10),
                    labels=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2],10), 
                    cex.axis=0.8),
     legend.args=list(text=' Temperature (°C)', side=4, font=2, line=2.5, cex=0.8))
dev.off()

#<< MAP >> -------------------
pdf(file = paste0(fig.dir, 'MAP_map.pdf'),
    width = 12.5, height = 8.5)

Ann.prec.col<-colorRampPalette(c("gold","yellow","cyan","deepskyblue","lightblue","cornflowerblue", "blue"))
plot(Bio.sub.mask$bio12_12, main="Annual Precipitation",xlab= "Longitude", ylab="Latitude", axes=FALSE, box=FALSE, col=Ann.prec.col(7))
plot(Arizona,add=TRUE)
axis(1)
axis(2)

#add sample points to the ann. prec. map, group them by colors based on biotic communities
color.vec<-c("green", "gray","blue")
plot(Bio.sub.mask$bio12_12, legend=FALSE,axes=FALSE, box=FALSE, col=Ann.prec.col(7) ,xlab= "Longitude", ylab="Latitude")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=0.8, tck=-0.02)
axis(2, las=2, cex.axis=0.8, tck=-0.02)
points(sites$long, sites$lat, pch=sites$shape, col='black', cex=4)
#text(sites$Longitude, sites$Latitude, label=rownames(sites), cex=0.5) 
#site.pick<-sites[c(45),]
#text(site.pick$Longitude,site.pick$Latitude, label=rownames(site.pick))

#add group legend outside the map region
par(xpd=TRUE)
legend(-108.6,37.2, legend=levels(sites$Biotic.community),bty="n", pch=19,col=color.vec, pt.cex=1.6, cex=1)
#temperature legend
Bio.sub.mask.range<-c(100,1000)
plot(Bio.sub.mask$bio12_12, legend.only=TRUE,col=Ann.prec.col(7),
     legend.width=2, legend.shrink=0.5,
     axis.args=list(at=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2], 200),
                    labels=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2],200), 
                    cex.axis=0.8),
     legend.args=list(text=' Precipitation (mm)', side=4, font=2, line=2.5, cex=0.8))

dev.off()

#<< BIO17*****Change colors >> -------------------
pdf(file = paste0(fig.dir, 'BIO17_map.pdf'),
    width = 12.5, height = 8.5)

Ann.prec.col<-colorRampPalette(c("gold","yellow","cyan","deepskyblue","lightblue","cornflowerblue", "blue"))
plot(Bio.sub.mask$bio17_12, main="Annual Precipitation",xlab= "Longitude", ylab="Latitude", axes=FALSE, box=FALSE, col=Ann.prec.col(7))
plot(Arizona,add=TRUE)
axis(1)
axis(2)

#add sample points to the ann. prec. map, group them by colors based on biotic communities
plot(Bio.sub.mask$bio17_12, legend=FALSE,axes=FALSE, box=FALSE, col=Ann.prec.col(7) ,xlab= "Longitude", ylab="Latitude")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=0.8, tck=-0.02)
axis(2, las=2, cex.axis=0.8, tck=-0.02)
points(sites$long, sites$lat, pch=sites$shape, col='black', cex=2.2)
#text(sites$Longitude, sites$Latitude, label=rownames(sites), cex=0.5) 
#site.pick<-sites[c(45),]
#text(site.pick$Longitude,site.pick$Latitude, label=rownames(site.pick))

#add group legend outside the map region
par(xpd=TRUE)
legend(-108.6,37.2, legend=levels(sites$Biotic.community),bty="n", pch=19,col=color.vec, pt.cex=1.6, cex=1)
#temperature legend
Bio.sub.mask.range<-c(100,1000)
plot(Bio.sub.mask$bio17_12, legend.only=TRUE,col=Ann.prec.col(7),
     legend.width=2, legend.shrink=0.5,
     axis.args=list(at=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2], 200),
                    labels=seq(Bio.sub.mask.range[1], Bio.sub.mask.range[2],200), 
                    cex.axis=0.8),
     legend.args=list(text='Mean precipitation\nof the driest quarter (mm)', side=4, font=2, line=2.5, cex=0.8))

dev.off()

#<< Elevation >> -------------------
pdf(file = paste0(fig.dir, 'Elevation_map.pdf'),
    width = 12.5, height = 8.5)
Alt<-getData("worldclim", var="alt", res=0.5, lon=-110, lat=33)
Alt.sub<-crop(Alt, extent(Arizona))
Alt.sub.mask<-mask(Alt.sub, Arizona)
Alt.col<-colorRampPalette(c("white", "gray90","gray80","gray70","gray60","gray50","gray30","black"))

# point color indicate geographic area
plot(Alt.sub.mask$alt_12,legend=FALSE,axes=FALSE, box=FALSE, col=Alt.col(8) ,xlab= "Longitude", ylab="Latitude")
plot(Arizona,add=TRUE)
axis(1, pos=31, cex.axis=1, tck=-0.02)
axis(2, las=2, cex.axis=1, tck=-0.02)
points(sites$long, sites$lat, pch=sites$shape, col='black', cex=4)
#text(sites$Longitude, sites$Latitude, label=rownames(sites), cex=0.5) 

par(xpd=TRUE)
legend(-108.5,37.2, legend=levels(sites$Area.code),bty="n", pch=19,col=alpha(color.vec,0.5), pt.cex=1.6, cex=1)
Alt.sub.mask.range<-c(500,maxValue(Alt.sub.mask))
plot(Alt.sub.mask$alt_12, legend.only=TRUE,col=Alt.col(8),
     legend.width=2, legend.shrink=0.5,
     axis.args=list(at=seq(Alt.sub.mask.range[1], Alt.sub.mask.range[2], 1000),
                    labels=seq(Alt.sub.mask.range[1], Alt.sub.mask.range[2],1000), 
                    cex.axis=0.9),
     legend.args=list(text='Elevation (m)', side=4, font=2, line=2.5, cex=0.8))
  
dev.off()
