## Script created by Liz Bowman August 9, 2016
## for downloading and isolating climate data from site coordinates

#=========================================================================================
# Climate data
#=========================================================================================

#--First, you need to go to http://www.worldclim.org/bioclim and download climate data 
#--for the quadrant where you sites are. Click on Download tab at top of page and download
#--past, current, or future climate data.
#--Be sure to download by tile (there should be a link at the top of the page for it). 
#--Bioclim data downloads a ton of different variables (http://www.worldclim.org/bioclim)

#-----------------------------------------------------------------------------------------
# download climate data from Worldclim using raster package
#-----------------------------------------------------------------------------------------
#--path to directory where the climate data downloaded from BIOCLIM site is stored
clim.dir <- "~/Dropbox/TechnicalDevelopment/WorldClimData/wc2.1_2.5m_bio/"

#--install and set library
#--uncomment install packages function if you don't have this package installed
#install.packages ("raster")
library(raster)

# << BIOCLIM VARIABLES >> ----------------------------------------------------------------
#--load world climate data: for 19 bio variables
#--no need to change this if you downloaded tile 12
data.1 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_1.tif'),native = F)
data.2 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_2.tif'),native = F)
data.3 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_3.tif'),native = F)
data.4 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_4.tif'),native = F)
data.5 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_5.tif'),native = F)
data.6 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_6.tif'),native = F)
data.7 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_7.tif'),native = F)
data.8 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_8.tif'),native = F)
data.9 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_9.tif'),native = F)
data.10 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_10.tif'),native = F)
data.11 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_11.tif'),native = F)
data.12 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_12.tif'),native = F)
data.13 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_13.tif'),native = F)
data.14 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_14.tif'),native = F)
data.15 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_15.tif'),native = F)
data.16 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_16.tif'),native = F)
data.17 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_17.tif'),native = F)
data.18 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_18.tif'),native = F)
data.19 <- raster(paste0(clim.dir,'wc2.1_2.5m_bio_19.tif'),native = F)

#<< Individual variables >> --------------------------------------------------------------
#--load world climate data: Precipitation, Tavg
#--there is a file for each month
#--no need to change this if you downloaded for tile 12
#--Precipitation
prec.1 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_01.tif'),native = F)
prec.2 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_02.tif'),native = F)
prec.3 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_03.tif'),native = F)
prec.4 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_04.tif'),native = F)
prec.5 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_05.tif'),native = F)
prec.6 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_06.tif'),native = F)
prec.7 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_07.tif'),native = F)
prec.8 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_08.tif'),native = F)
prec.9 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_09.tif'),native = F)
prec.10 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_10.tif'),native = F)
prec.11 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_11.tif'),native = F)
prec.12 <- raster(paste0(clim.dir, 'prec_12/wc2.0_30s_prec_12.tif'),native = F)
#--Average temperature
Tavg.1 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_01.tif'),native = F)
Tavg.2 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_02.tif'),native = F)
Tavg.3 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_03.tif'),native = F)
Tavg.4 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_04.tif'),native = F)
Tavg.5 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_05.tif'),native = F)
Tavg.6 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_06.tif'),native = F)
Tavg.7 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_07.tif'),native = F)
Tavg.8 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_08.tif'),native = F)
Tavg.9 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_09.tif'),native = F)
Tavg.10 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_10.tif'),native = F)
Tavg.11 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_11.tif'),native = F)
Tavg.12 <- raster(paste0(clim.dir, 'tavg_12/wc2.0_30s_tavg_12.tif'),native = F)

#--import csv data frame with site name, lat, long
#--make sure column names are site, lat, and long for next set of commands
sites <- read.csv(paste0(dat.dir,'20190923_SamplingData.csv'), as.is = T, header = T)

#--combine lat and long
xy <- cbind (sites$long, sites$lat)
#--create object of class SpatialPoints from coordinates
sp <- SpatialPoints(xy)
#--extract information from downloaded files
sites_cli <- as.data.frame(cbind(sites,
                                 "BIO1" = extract(data.1, sp, method = 'bilinear'),
                                 "BIO2" = extract(data.2, sp, method = 'bilinear'),
                                 "BIO3" = extract(data.3, sp, method = 'bilinear'),
                                 "BIO4" = extract(data.4, sp, method = 'bilinear'),
                                 "BIO5" = extract(data.5, sp, method = 'bilinear'),
                                 "BIO6" = extract(data.6, sp, method = 'bilinear'),
                                 "BIO7" = extract(data.7, sp, method = 'bilinear'),
                                 "BIO8" = extract(data.8, sp, method = 'bilinear'),
                                 "BIO9" = extract(data.9, sp, method = 'bilinear'),
                                 "BIO10" = extract(data.10, sp, method = 'bilinear'),
                                 "BIO11" = extract(data.11, sp, method = 'bilinear'),
                                 "BIO12" = extract(data.12, sp, method = 'bilinear'),
                                 "BIO13" = extract(data.13, sp, method = 'bilinear'),
                                 "BIO14" = extract(data.14, sp, method = 'bilinear'),
                                 "BIO15" = extract(data.15, sp, method = 'bilinear'),
                                 "BIO16" = extract(data.16, sp, method = 'bilinear'),
                                 "BIO17" = extract(data.17, sp, method = 'bilinear'),
                                 "BIO18" = extract(data.18, sp, method = 'bilinear'),
                                 "BIO19" = extract(data.19, sp, method = 'bilinear')))

#--average temperature and precipitation
# sites_cli2 <- as.data.frame (cbind (sites,
#                                    "Prec1" = extract(prec.1, sp, method = 'bilinear'),
#                                    "Prec2" = extract(prec.2, sp, method = 'bilinear'),
#                                    "Prec3" = extract(prec.3, sp, method = 'bilinear'),
#                                    "Prec4" = extract(prec.4, sp, method = 'bilinear'),
#                                    "Prec5" = extract(prec.5, sp, method = 'bilinear'),
#                                    "Prec6" = extract(prec.6, sp, method = 'bilinear'),
#                                    "Prec7" = extract(prec.7, sp, method = 'bilinear'),
#                                    "Prec8" = extract(prec.8, sp, method = 'bilinear'),
#                                    "Prec9" = extract(prec.9, sp, method = 'bilinear'),
#                                    "Prec10" = extract(prec.10, sp, method = 'bilinear'),
#                                    "Prec11" = extract(prec.11, sp, method = 'bilinear'),
#                                    "Prec12" = extract(prec.12, sp, method = 'bilinear'),
#                                    "Tavg1" = extract(Tavg.1, sp, method = 'bilinear'),
#                                    "Tavg2" = extract(Tavg.2, sp, method = 'bilinear'),
#                                    "Tavg3" = extract(Tavg.3, sp, method = 'bilinear'),
#                                    "Tavg4" = extract(Tavg.4, sp, method = 'bilinear'),
#                                    "Tavg5" = extract(Tavg.5, sp, method = 'bilinear'),
#                                    "Tavg6" = extract(Tavg.6, sp, method = 'bilinear'),
#                                    "Tavg7" = extract(Tavg.7, sp, method = 'bilinear'),
#                                    "Tavg8" = extract(Tavg.8, sp, method = 'bilinear'),
#                                    "Tavg9" = extract(Tavg.9, sp, method = 'bilinear'),
#                                    "Tavg10" = extract(Tavg.10, sp, method = 'bilinear'),
#                                    "Tavg11" = extract(Tavg.11, sp, method = 'bilinear'),
#                                    "Tavg12" = extract(Tavg.12, sp, method = 'bilinear')))

# #--Take average of precipitation and add it to the sites.cli file
# sites_cli2$prec <- (sites_cli2$Prec1 + sites_cli2$Prec2 + sites_cli2$Prec3 +
#                      sites_cli2$Prec4 + sites_cli2$Prec5 + sites_cli2$Prec6 + sites_cli2$Prec7 +
#                      sites_cli2$Prec8 + sites_cli2$Prec9 + sites_cli2$Prec10 + sites_cli2$Prec11 +
#                      sites_cli2$Prec12)/12
# #--Average temperature
# sites_cli2$Tavg <- (sites_cli2$Tavg1 + sites_cli2$Tavg2 + sites_cli2$Tavg3 +
#                      sites_cli2$Tavg4 + sites_cli2$Tavg5 + sites_cli2$Tavg6 + sites_cli2$Tavg7 +
#                      sites_cli2$Tavg8 + sites_cli2$Tavg9 + sites_cli2$Tavg10 + sites_cli2$Tavg11 +
#                      sites_cli2$Tavg12)/12
# 
# sites_cli2 <- sites_cli2[c(1:26, 51, 52)]


#--write file
write.csv(sites_cli,'data/20201201_ClimateData.csv', row.names = F)
