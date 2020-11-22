## Script created by Liz Bowman January 24, 2020
## for analysing community data in Sky Island study
## Includes ectomycorrhizal and endophyte data

## For answering 'How does ectomycorrhizal and endophytic fungi vary taxonomically across AZ?' and
## Does this variation correlate to environment?
# Part A focuses on taxonomic variation across arizona
# Part B focuses on how this aligns with environment which could indicate certain groups that
# are more fit in particular environments (warm, dry at lower lat and cool, wet at higher)

#========================================================================================#
# Load libraries and data------------------
#========================================================================================#

#install.packages(c('dplyr','tidyr','vegan','ggplot2','ecodist'))
library(dplyr);library(tidyr);library(vegan);library(ggplot2)
library(nlme);library(ecodist); library(cluster)

#--file paths
dat.dir <- 'data/'
fig.dir <- '~/Documents/PhD/3:4_Combined/figures/'
res.dir <- '~/Documents/PhD/3:4_Combined/results/'

#--Read in taxonomic data for EM
em.tax.data <- read.csv(paste0(dat.dir, '20190919_AllEMSequences.csv'), as.is = T)
# Remove non focal sites in Santa Catalina Mts.
non.focal.scm <- c('Middle Bear','Bear Wallow','Pullout south of Willow Canyon Rd.',
                   'Parking pullout South of Rose Canyon','Green Mountain/Bug Springs',
                   'Marshall Gulch','Rose Canyon Rd.')
em.tax.data <- em.tax.data[!em.tax.data$Site %in% non.focal.scm, ]

#--environmental data
env.data <- read.csv(paste0(dat.dir,'EM_EnvFactors_site.csv'), as.is = T)
rownames(env.data) <- env.data$site

#--spatial distance
spat.data <- read.csv(paste0(dat.dir, 'DistanceClosestPPForest_site.csv'),as.is = T)

#========================================================================================#
# Plot and Chi square analysis-----------------
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Class level----
#----------------------------------------------------------------------------------------#
# Make table chi-square analyses
em.tax.data %>%
  select(SequenceName, Site, Class, Range, MorphologyAbundance) %>%
  filter(!Class == 'Sordariomycetes ', !is.na(Class)) %>%
  group_by(Range, Class) %>%
  summarize(abund.tip = sum(MorphologyAbundance)) %>%
  spread(key = Class, value = abund.tip, fill = 0) -> class.em

# Make table chi-square analyses
em.tax.data %>%
  select(SequenceName, Site, Class, Range, MorphologyAbundance) %>%
  filter(!Class == 'Sordariomycetes ', !is.na(Class)) %>%
  group_by(Site, Range, Class) %>%
  summarize(abund.tip = sum(MorphologyAbundance)) %>%
  spread(key = Class, value = abund.tip, fill = 0) -> mrm.class.em

#<< Chi-square analysis >> --------------------------------------------------------------#
class.em.chisq <- class.em[-1]

#--Chi square
chisq.test(class.em.chisq)

#<< Multiple regression on distance matrices >> -----------------------------------------#
mrm.class.em.mod <- mrm.class.em[-c(1,2)]

#--community distance matrix
mrm.horn <- vegdist(mrm.class.em.mod, method = 'horn')

#--climate distance matrix
mrm.env <- data.frame(site = mrm.class.em $Site,
                       clim.pca = NA,
                       forest.type = NA)
for(i in mrm.env$site){
  clim.i <- env.data[env.data$site == i, 'clim.pca']
  forest.i <- env.data[env.data$site == i, 'forest']
  dist.i <- spat.data[spat.data$site == i, 'dist']
  range.i <- env.data[env.data$site == i, 'range']
  BIO1.i <- env.data[env.data$site == i, 'BIO1']
  BIO12.i <- env.data[env.data$site == i, 'BIO12']
  BIO16.i <- env.data[env.data$site == i, 'BIO16']
  BIO17.i <- env.data[env.data$site == i, 'BIO17']
  mrm.env[mrm.env$site == i, 'clim.pca'] <- clim.i
  mrm.env[mrm.env$site == i, 'forest.type'] <- forest.i
  mrm.env[mrm.env$site == i, 'closest.PP.forest'] <- dist.i
  mrm.env[mrm.env$site == i, 'range'] <- range.i
  mrm.env[mrm.env$site == i, 'BIO1'] <- BIO1.i
  mrm.env[mrm.env$site == i, 'BIO12'] <- BIO12.i
  mrm.env[mrm.env$site == i, 'BIO16'] <- BIO16.i
  mrm.env[mrm.env$site == i, 'BIO17'] <- BIO17.i
}
# Site
bio1.site.dist <- vegdist(mrm.env$BIO1, method = 'euclidean')
bio12.site.dist <- vegdist(mrm.env$BIO12, method = 'euclidean')
forest.gower <- data.frame(forest.type = mrm.env$forest.type)
forest.site.dist <- daisy(forest.gower, metric = 'gower')
clim.site.dist <- vegdist(mrm.env$clim.pca, method = 'euclidean')

#--Geographical distance to closest PP forest
# Site
spat.site.dist <- vegdist(mrm.env$closest.PP.forest, method = 'euclidean')

#--MRM analysis
mrm.sep.class <- MRM(mrm.horn ~ spat.site.dist * bio12.site.dist * forest.site.dist,
                 method = 'linear')
mrm.clim.class <- MRM(mrm.horn ~ spat.site.dist * clim.site.dist * forest.site.dist,
                 method = 'linear')

# << Class level plot >> ----------------------------------------------------------#

em.tax.data %>%
  filter(!is.na(Class.i),!Class.i == 'Sordariomycetes') -> em.tax.class

#--Ranges ordered by MAT (high to low)--------------------------------------------#
em.class.plot.dist <- em.tax.class

# Add forest data
for(i in unique(em.class.plot.dist$Site)){
  em.class.plot.dist[em.class.plot.dist$Site == i, 'forest.type'] <- env.data[env.data$site == i, 'forest']
}

# Ordered from low to long distance
em.clss.plot.dist$Range <- factor(em.class.plot.dist$Range,
                                   levels = c('MogollonRim','Flagstaff','MtLemmon','Mingus',
                                              'Bradshaw','Pinaleno','Huachucas','Hualapai','Chiricahuas'))

em.class.dist <- ggplot(data = em.class.plot.dist,
       aes(x = Range,
           fill = Class.i)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = c('#a6611a','#dfc27d','#f5f5f5', '#80cdc1','#018571')) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x=element_text(angle = -90, hjust = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=16, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_Class_closestPP.pdf', plot = em.class.dist, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--Ranges ordered by MAP (low to high)--------------------------------------------#
em.class.plot.prec <- em.tax.class

# Add forest data
for(i in unique(em.class.plot.prec$Site)){
  em.class.plot.prec[em.class.plot.prec$Site == i, 'forest.type'] <- env.data[env.data$site == i, 'forest']
}

# Ordered from low MAP to high MAP
em.class.plot.prec$Range <- factor(em.class.plot.prec$Range,
                              levels = c('Hualapai','Flagstaff','MogollonRim','Bradshaw',
                                         'Chiricahuas','Huachucas','Mingus','Pinaleno','MtLemmon'))

em.class.prec <- ggplot(data = em.class.plot.prec,
       aes(x = Range,
           fill = Class.i)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = c('#a6611a','#dfc27d','#f5f5f5', '#80cdc1','#018571')) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x=element_text(angle = -90, hjust = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=16, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_Class_prec.pdf', plot = em.class.prec, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#----------------------------------------------------------------------------------------#
# Order level----
#----------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------#
# Make table for analyses
#----------------------------------------------------------------------------------------#
# Make table for chi-square analyses
em.tax.data %>%
  select(SequenceName, Site, Order, Range, MorphologyAbundance) %>%
  filter(!Order %in% c('Corticiales','Eurotiales','Helotiales','Hypocreales','Sebacinales'),
         !is.na(Order)) %>%
  group_by(Range, Order) %>%
  summarize(abund.tip = sum(MorphologyAbundance)) %>%
  spread(key = Order, value = abund.tip, fill = 0) -> order.em

# Make table for mrm analyses
em.tax.data %>%
  select(SequenceName, Site, Order, Range, MorphologyAbundance) %>%
  filter(!Order %in% c('Corticiales','Eurotiales','Helotiales','Hypocreales','Sebacinales'),
         !is.na(Order)) %>%
  group_by(Site, Range, Order) %>%
  summarize(abund.tip = sum(MorphologyAbundance)) %>%
  spread(key = Order, value = abund.tip, fill = 0) -> mrm.order.em

#<< Chi-square analysis >> --------------------------------------------------------------#
em.order.table <- order.em[-1]

#--Chi square
chisq.test(em.order.table)

#<< Multiple regression on distance matrices >> -----------------------------------------#
mrm.order.table <- mrm.order.em[-c(1,2)]

#--community distance matrix
mrm.horn <- vegdist(mrm.order.table, method = 'horn')

#--Environmental distance matrix
mrm.env <- data.frame(site = mrm.order.em$Site,
                      clim.pca = NA,
                      forest.type = NA)
for(i in mrm.order.em$Site){
  clim.i <- env.data[env.data$site == i, 'clim.pca']
  forest.i <- env.data[env.data$site == i, 'forest']
  dist.i <- spat.data[spat.data$site == i, 'dist']
  range.i <- env.data[env.data$site == i, 'range']
  BIO1.i <- env.data[env.data$site == i, 'BIO1']
  BIO12.i <- env.data[env.data$site == i, 'BIO12']
  BIO16.i <- env.data[env.data$site == i, 'BIO16']
  BIO17.i <- env.data[env.data$site == i, 'BIO17']
  mrm.env[mrm.env$site == i, 'clim.pca'] <- clim.i
  mrm.env[mrm.env$site == i, 'forest.type'] <- forest.i
  mrm.env[mrm.env$site == i, 'closest.PP.forest'] <- dist.i
  mrm.env[mrm.env$site == i, 'range'] <- range.i
  mrm.env[mrm.env$site == i, 'BIO1'] <- BIO1.i
  mrm.env[mrm.env$site == i, 'BIO12'] <- BIO12.i
  mrm.env[mrm.env$site == i, 'BIO16'] <- BIO16.i
  mrm.env[mrm.env$site == i, 'BIO17'] <- BIO17.i
}

# Environmental distance matrices
bio1.site.dist <- vegdist(mrm.env$BIO1, method = 'euclidean')
bio12.site.dist <- vegdist(mrm.env$BIO12, method = 'euclidean')
forest.gower <- data.frame(forest.type = mrm.env$forest.type, metric = 'gower')
forest.site.dist <- daisy(forest.gower, metric = 'gower')
clim.site.dist <- vegdist(mrm.env$clim.pca, method = 'euclidean')

#--Geographical distance to closest PP forest
spat.site.dist <- vegdist(mrm.env$closest.PP.forest, method = 'euclidean')

#--MRM analysis
mrm.sep.order <- MRM(mrm.horn ~ spat.site.dist + bio12.site.dist + forest.site.dist,
                 method = 'linear')
mrm.class.order <- MRM(mrm.horn ~ spat.site.dist * clim.site.dist * forest.site.dist,
                     method = 'linear')

# << Order level plot >> ----------------------------------------------------------#
em.tax.data %>%
  filter(!is.na(Order.i),
         !Order.i %in% c('Corticiales','Eurotiales','Helotiales','Hypocreales','Sebacinales')) %>%
  group_by(Site, Range, Order.i) %>%
  summarize(total = sum(MorphologyAbundance)) %>%
  as.data.frame(.) -> em.tax.order

range <- NA
order <- NA
for(i in 1:nrow(em.tax.order)){
  rep.i <- em.tax.order[i,'total']
  range.i <- em.tax.order[i,'Range']
  order.i <- em.tax.order[i,'Order.i']
  range <- append(range, rep(range.i, rep.i))
  order <-  append(order, rep(order.i, rep.i))  
}

em.tax.order <- data.frame(Range = range, Order.i = order)
em.tax.order <- em.tax.order[!is.na(em.tax.order$Range),]

#--Ranges ordered by Distance to the closest PP forest--------------------------------------#
em.order.plot.dist <- em.tax.order

# Ordered from low to long distance
em.order.plot.dist$Range <- factor(em.order.plot.dist$Range,
                                   levels = c('MogollonRim','Flagstaff','MtLemmon','Mingus',
                                              'Bradshaw','Pinaleno','Huachucas','Hualapai','Chiricahuas'))

em.order.dist <- ggplot(data = em.order.plot.dist,
                        aes(x = Range,
                            fill = Order.i)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = c('#8c510a','#bf812d','#dfc27d','#f6e8c3',
                               '#c7eae5','#80cdc1','#35978f','#01665e')) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x=element_text(angle = -90, hjust = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=16, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_Order_closestPP.pdf', plot = em.order.dist, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--Ranges ordered by MAP (low to high)--------------------------------------------#
em.order.plot.prec <- em.tax.order

# Ordered from low MAP to high MAP
em.order.plot.prec$Range <- factor(em.order.plot.prec$Range,
                                   levels = c('Hualapai','Flagstaff','MogollonRim','Bradshaw',
                                              'Chiricahuas','Huachucas','Mingus','Pinaleno','MtLemmon'))

em.order.prec <- ggplot(data = em.order.plot.prec,
       aes(x = Range,
           fill = Order.i)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = c('#8c510a','#bf812d','#dfc27d','#f6e8c3',
                               '#c7eae5','#80cdc1','#35978f','#01665e')) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x=element_text(angle = -90, hjust = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=16, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_Order_prec.pdf', plot = em.order.prec, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)


#----------------------------------------------------------------------------------------#
# Genus level----
#----------------------------------------------------------------------------------------#
em.tax.data %>%
  filter(!is.na(Genus.i),
#         !Genus.i %in% c('Amanita','Amphinema','Boletus','Entoloma','Genea','Helvella',
#                        'Hydnobolites','Hygrophorus','Peziza','Pseudotomentella',
#                        'Sistotrema','Suillus','Tylospora')) %>%
        Order.i %in% c('Thelephorales','Russulales','Agaricales')) %>%
  group_by(Site, Range, Genus.i) %>%
  summarize(total = sum(MorphologyAbundance)) %>%
  as.data.frame(.) -> em.tax.tip.genus

#--Ranges ordered by Distance to the closest PP forest--------------------------------------#
em.genus.plot.dist <- em.tax.tip.genus

# Ordered from low to long distance
em.genus.plot.dist$Range <- factor(em.genus.plot.dist$Range,
                                   levels = c('MogollonRim','Flagstaff','MtLemmon','Mingus',
                                              'Bradshaw','Pinaleno','Huachucas','Hualapai','Chiricahuas'))

em.Genus.dist <- ggplot(data = em.genus.plot.dist,
                        aes(x = Range,
                            fill = Genus.i)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x=element_text(angle = -90, hjust = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=16, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_Genus_closestPP.pdf', plot = em.Genus.dist, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--Ranges ordered by MAP (low to high)--------------------------------------#

em.genus.plot.prec <- em.tax.genus

# Ordered from low MAP to high MAP
em.genus.plot.prec$Range <- factor(em.order.plot.prec$Range,
                                   levels = c('Hualapai','Flagstaff','MogollonRim','Bradshaw',
                                              'Chiricahuas','Huachucas','Mingus','Pinaleno','MtLemmon'))

em.Genus.prec <- ggplot(data = em.genus.plot.dist,
                        aes(x = Range,
                            fill = Genus.i)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per genus") +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        axis.text.x=element_text(angle = -90, hjust = 0),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=16, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_Genus_Prec.pdf', plot = em.Genus.prec, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)
