## Script created by Liz Bowman December 26, 2019
## for analyzing climate data

#=========================================================================================
# Analysis of climate data
#=========================================================================================
#--load libraries
library(ggplot2)
#install.packages('cowplot')
library(cowplot)
#intall.packages('stats')
library(stats)

#-----------------------------------------------------------------------------------------
# load file paths and data
#-----------------------------------------------------------------------------------------
#--file paths
dat.dir <- '~/Documents/PhD/3:4_Combined/data/'
fig.dir <- '~/Documents/PhD/3:4_Combined/figures/'
res.dir <- '~/Documents/PhD/3:4_Combined/results/'
#--load data
sites_cli <- read.csv(paste0(dat.dir, '20190923_SamplingData.csv'), as.is = T, header = T)

#-----------------------------------------------------------------------------------------
# Create plot with average climate data per range plotted as a function of latitude
#-----------------------------------------------------------------------------------------

# << create data frame >> ----------------------------------------------------------------
avg.cl <- data.frame(range = c(rep('chiricahua',3),rep('pinaleno',3),rep('mogollon',3),
                               rep('bradshaw',3),rep('santa.catalina',3),
                               rep('huachuca',3), rep('mingus',3), rep('flagstaff',3),
                               rep('hualapai',3)),
                     lat = unique(sites_cli$lat),
                     tmax = NA, tmin = NA, prec = NA, avg.warm.quarter = NA)
for(l in avg.cl$lat) {
  for(r in avg.cl$range) {
    # avg.cl[avg.cl$range == r & avg.cl $ lat == l, 'tmax'] <-
    #   mean(sites_cli[sites_cli$range == r & sites_cli$lat == l, 'Tmax'])
    # avg.cl[avg.cl$range == r & avg.cl $ lat == l, 'tmin'] <-
    #   mean(sites_cli[sites_cli$range == r & sites_cli$lat == l, 'Tmin'])
    # avg.cl[avg.cl$range == r & avg.cl $ lat == l, 'prec'] <-
    #   mean(sites_cli[sites_cli$range == r & sites_cli$lat == l, 'prec'])
    avg.cl[avg.cl$range == r & avg.cl $ lat == l, 'avg.prec'] <-
      mean(sites_cli[sites_cli$range == r & sites_cli$lat == l, 'BIO12'])
    avg.cl[avg.cl$range == r & avg.cl $ lat == l, 'avg.cold.quarter'] <-
      mean(sites_cli[sites_cli$range == r & sites_cli$lat == l, 'BIO11_red'])
    avg.cl[avg.cl$range == r & avg.cl $ lat == l, 'avg.warm.quarter'] <-
      mean(sites_cli[sites_cli$range == r & sites_cli$lat == l, 'BIO10_red'])
  }
}

# << plot climate variables >> -----------------------------------------------------------
# #--max temperature across sites and ranges
# max <- ggplot(data = avg.cl,
#        mapping = aes(x = lat,
#                      y = tmax,
#                      color = range)) +
#   geom_point() +
#   xlab('Latitude') +
#   ylab('Max. temp. (°C)') +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_hue()
# #--min temperature across sites and ranges
# min <- ggplot(data = avg.cl,
#        mapping = aes(x = lat,
#                      y = tmin,
#                      color = range)) +
#   geom_point() +
#   xlab('Latitude') +
#   ylab('Min. temp. (°C)') +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_hue()
# #--precipitation across sites and ranges
# prec <- ggplot(data = avg.cl,
#        mapping = aes(x = lat,
#                      y = prec,
#                      color = range)) +
#   geom_point() +
#   xlab('Latitude') +
#   ylab('Precipitation (mm)') +
#   theme_bw() + 
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_hue()
#--Avg. temperature warm quarter across sites and ranges
warm <- ggplot(data = avg.cl,
       mapping = aes(x = lat,
                     y = avg.warm.quarter,
                     color = range)) +
  geom_point() +
  xlab('Latitude') +
  ylab('Avg. temp. warm quarter (°C)') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_hue()

#--Avg. temperature coldest quarter across sites and ranges
cold <- ggplot(data = avg.cl,
               mapping = aes(x = lat,
                             y = avg.cold.quarter,
                             color = range)) +
  geom_point() +
  xlab('Latitude') +
  ylab('Avg. temp. cold quarter (°C)') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_hue()

#--Avg. temperature warm quarter across sites and ranges
prec <- ggplot(data = avg.cl,
               mapping = aes(x = lat,
                             y = avg.prec,
                             color = range)) +
  geom_point() +
  xlab('Latitude') +
  ylab('Avg. annual precipiation (mm)') +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_hue()

#--arrange plots into one figure
clim <- plot_grid(max + theme(legend.position = "none"),
          min + theme(legend.position = "none"),
          prec + theme(legend.position = "none"),
          warm+ theme(legend.position = "none"),
          labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2, hjust = -2)
#--add legend
legend <- get_legend(max)
#--add legend to plot
plot_grid(clim, legend, rel_widths = c(1, .3))

# << check for sig differences between ranges >> -----------------------------------------
#--max temperature 
max <- aov(avg.cl$tmax ~ avg.cl$lat)
max <- summary(max)
TukeyHSD(aov(avg.cl$tmax ~ avg.cl$range), ordered = T)

#--min temperature 
min <- aov(avg.cl$tmin ~ avg.cl$lat)
min <- summary(min)
TukeyHSD(aov(avg.cl$tmin ~ avg.cl$range), ordered = T)

#--precipitation
#--max temperature 
prec <- aov(avg.cl$prec ~ avg.cl$lat)
prec <- summary(prec)
TukeyHSD(aov(avg.cl$prec ~ avg.cl$range), ordered = T)

#--Avg. temp. warm quarter
#--max temperature 
warm <- aov(avg.cl$avg.warm.quarter ~ avg.cl$lat)
warm <- summary(warm)
TukeyHSD(aov(avg.cl$avg.warm.quarter ~ avg.cl$range), ordered = T)

