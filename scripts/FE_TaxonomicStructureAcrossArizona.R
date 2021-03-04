## Script created by Liz Bowman January 28, 2020
## for analysing community data in Sky Island study
## Includes endophyte data

## For answering 'How do endophytic fungi vary taxonomically across AZ?' and
## Does this variation correlate to environment?
# Part A focuses on taxonomic variation across arizona
# Part B focuses on how this aligns with environment which could indicate certain groups that
# are more fit in particular environments (warm, dry at lower lat and cool, wet at higher)

#========================================================================================#
# Load libraries and data------------------
#========================================================================================#

fe.cf.data <- read.csv(paste0(dat.dir,'FE_CF_SeqSimTaxonomy.csv'), as.is = T)
fe.cb.data <- read.csv(paste0(dat.dir,'FE_CB_SeqSimTaxonomy.csv'), as.is = T)

# remove nonfocal sites 
fe.cf.data <- fe.cf.data[!fe.cf.data$Site %in% c('Bear Wallow','Marshall Gulch',
               'Parking pullout South of Rose Canyon','Rose Canyon Rd.','Pullout south of Willow Canyon Rd.',
               'Middle Bear','Green Mountain/Bug Springs'),]
fe.cb.data <- fe.cb.data[!fe.cb.data$site %in% c('Bear Wallow','Marshall Gulch',
               'Parking pullout South of Rose Canyon','Rose Canyon Rd.','Pullout south of Willow Canyon Rd.',
               'Middle Bear','Green Mountain/Bug Springs'),]

#--environmental data
env.cf.data <- read.csv(paste0(dat.dir,'FE_CF_EnvFactors_site.csv'), as.is = T)
rownames(env.cf.data) <- env.cf.data$site

#--environmental data for culture based data
env.data <- read.csv(paste0(dat.dir,'FE_CB_EnvFactors_site.csv'), as.is = T)
rownames(env.data) <- env.data$site

#--spatial distance
spat.data <- read.csv(paste0(dat.dir, 'DistanceClosestPPForest_site.csv'),as.is = T)

#========================================================================================#
# Plot and Chi square analysis-----------------
#========================================================================================#

#----------------------------------------------------------------------------------------#
# Class level: culture-based----
#----------------------------------------------------------------------------------------#
# Make table for chi-square analyses
fe.cb.data %>%
  filter(!is.na(Class.TBAS),
         !Class.TBAS %in% c('0','unclassified','Taphrinomycetes','Saccharomycetes',
                       'Lecanoromycetes','Schizosaccharomycetes','Xylonomycetes',
                       'Agaricomycetes')) %>%
  count(range,Class.TBAS) %>%
  spread(key = Class.TBAS, value = n, fill = 0) -> fe.cb.chisq

# Make table for mrm analyses
fe.cb.data %>%
  filter(!is.na(Class.TBAS),
         !Class.TBAS %in% c('0','unclassified','Taphrinomycetes','Saccharomycetes',
                            'Lecanoromycetes','Schizosaccharomycetes','Xylonomycetes',
                            'Agaricomycetes')) %>%
  count(site, range,Class.TBAS) %>%
  spread(key = Class.TBAS, value = n, fill = 0) -> fe.cb.mrm

#<< Chi-square analysis >> --------------------------------------------------------------#
fe.cb.chisq.table <- fe.cb.chisq[-1]

#--Chi square
chisq.test(fe.cb.chisq.table)

#<< Multiple regression on distance matrices >> -----------------------------------------#
fe.cb.mrm.table <- fe.cb.mrm[-c(1,2)]

#--community distance matrix
mrm.horn <- vegdist(fe.cb.mrm.table, method = 'horn')
mrm.jaccard <- vegdist(fe.cb.mrm.table, method = 'jaccard', binary = F)
mrm.euclidean <- vegdist(fe.cb.mrm.table, method = 'euclidean')

#--environment distance matrix
mrm.env <- data.frame(site = fe.cb.mrm$site,
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
bio17.site.dist <- vegdist(mrm.env$BIO17, method = 'euclidean')
bio12.site.dist <- vegdist(mrm.env$BIO12, method = 'euclidean')
forest.gower <- data.frame(forest.type = as.factor(mrm.env$forest.type))
forest.site.dist <- daisy(forest.gower, metric = 'gower')

#--Geographical distance to closest PP forest
# Site
spat.site.dist <- vegdist(mrm.env$closest.PP.forest, method = 'euclidean')

#--MRM analysis
mrm.cb.horn <- MRM(mrm.horn ~ forest.site.dist * bio17.site.dist * bio12.site.dist,
                 method = 'linear')
mrm.cb.jaccard <- MRM(mrm.jaccard ~ forest.site.dist * bio17.site.dist * bio12.site.dist,
                   method = 'linear')

# << Class level plot >> ----------------------------------------------------------#

# Filter data for plot
fe.cb.data %>%
  filter(!is.na(Class.TBAS),
         !Class.TBAS %in% c('0','unclassified','Taphrinomycetes','Saccharomycetes',
                    'Lecanoromycetes','Schizosaccharomycetes','Xylonomycetes',
                    'Agaricomycetes')) -> fe.cb.plot

# Add plant community data
for(i in unique(fe.cb.plot$site)){
  forest.i <- env.data[env.data$site == i, 'forest']
  fe.cb.plot[fe.cb.plot$site == i, 'forest'] <- forest.i
}
for(i in 1:nrow(fe.cb.plot)){
  if(fe.cb.plot[i,'forest'] == 'pine'){
    fe.cb.plot[i,'forest'] <- 'Pine'} else{fe.cb.plot[i,'forest'] <- 'Pine-Douglas fir'}
}

# Plot by range, arranged based on BIO17 (low to high)
fe.cb.plot$range <- factor(fe.cb.plot$range,
                           levels = c('Huachuca','Chiricahua','SantaCatalina',
                                      'Pinaleno','Bradshaw','MogollonRim','Mingus'))

fe.cb.range.plot <- ggplot(data = fe.cb.plot,
                           aes(x = range,
                               fill = Class.TBAS)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per class") +
  theme_bw() +
  xlab(element_blank()) +
  scale_fill_manual(values = c('#8c510a','#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
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

ggsave('Fig6F_FE_CB_Class_rangeBIO17.pdf', plot = fe.cb.range.plot, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#----------------------------------------------------------------------------------------#
# Class level: culture-free----
#----------------------------------------------------------------------------------------#
# Make table for chi-square analyses
fe.cf.data %>%
  filter(!is.na(Class.TBAS),
         !Class.TBAS %in% c('0','unclassified','Taphrinomycetes','Saccharomycetes',
                       'Lecanoromycetes','Schizosaccharomycetes','Xylonomycetes',
                       'Pneumocystidomycetes','Lichinomycetes','Orbiliomycetes',
                       'Arthoniomycetes','Geoglossomycetes','Xylonomycetes','Coniocybomycetes',
                       'Taphrinomycetes,Saccharomycetes,Pneumocystidomycetes,Schizosaccharomycetes,unclassified'),
         !Site == 'Butterfly trail') %>%
  dplyr::select(range,Class.TBAS,OTU) %>%
  distinct(OTU, Class.TBAS,range) %>%
  count(range,Class.TBAS) %>%
  spread(key = Class.TBAS, value = n, fill = 0) -> fe.cf.chisq

# Make table for mrm analyses
fe.cf.data %>%
  filter(!is.na(Class.TBAS),
         !Class.TBAS %in% c('0','unclassified','Taphrinomycetes','Saccharomycetes',
                       'Lecanoromycetes','Schizosaccharomycetes','Xylonomycetes',
                       'Pneumocystidomycetes','Lichinomycetes','Orbiliomycetes',
                       'Arthoniomycetes','Geoglossomycetes','Xylonomycetes','Coniocybomycetes',
                       'Taphrinomycetes,Saccharomycetes,Pneumocystidomycetes,Schizosaccharomycetes,unclassified'),
         !Site == 'Butterfly trail') %>%
  dplyr::select(Site,range,Class.TBAS,OTU) %>%
  distinct(OTU,Class.TBAS,range,Site) %>%
  count(Site,range,Class.TBAS) %>%
  spread(key = Class.TBAS, value = n, fill = 0) -> fe.cf.mrm

#<< Chi-square analysis >> --------------------------------------------------------------#
fe.cf.chisq.table <- fe.cf.chisq[-1]

#--Chi square
chisq.test(fe.cf.chisq.table)

#<< Multiple regression on distance matrices >> -----------------------------------------#
fe.cf.mrm.table <- fe.cf.mrm[-c(1,2)]

#--community distance matrix
mrm.jaccard <- vegdist(fe.cf.mrm.table, method = 'jaccard', binary = T)

#--environment distance matrix
mrm.env <- data.frame(site = fe.cf.mrm$Site,
                      clim.pca = NA,
                      forest.type = NA)
for(i in mrm.env$site){
  clim.i <- env.cf.data[env.cf.data$site == i, 'clim.pca']
  forest.i <- env.cf.data[env.cf.data$site == i, 'forest']
  dist.i <- spat.data[spat.data$site == i, 'dist']
  range.i <- env.cf.data[env.cf.data$site == i, 'range']
  BIO1.i <- env.cf.data[env.cf.data$site == i, 'BIO1']
  BIO12.i <- env.cf.data[env.cf.data$site == i, 'BIO12']
  BIO16.i <- env.cf.data[env.cf.data$site == i, 'BIO16']
  BIO17.i <- env.cf.data[env.cf.data$site == i, 'BIO17']
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
bio17.site.dist <- vegdist(mrm.env$BIO17, method = 'euclidean')
bio12.site.dist <- vegdist(mrm.env$BIO12, method = 'euclidean')
forest.gower <- data.frame(forest.type = mrm.env$forest.type)
forest.site.dist <- daisy(forest.gower, metric = 'gower')

#--Geographical distance to closest PP forest
# Site
spat.site.dist <- vegdist(mrm.env$closest.PP.forest, method = 'euclidean')

#--MRM analysis
mrm.cf <- MRM(mrm.jaccard ~ forest.site.dist + bio17.site.dist + bio12.site.dist,
              method = 'linear')

# << Class level plot >> ----------------------------------------------------------#
# Filter data for plot
fe.cf.data %>%
  filter(!is.na(Class.TBAS),
         !Class.TBAS %in% c('0','unclassified','Taphrinomycetes','Saccharomycetes',
                       'Lecanoromycetes','Schizosaccharomycetes','Xylonomycetes',
                       'Pneumocystidomycetes','Lichinomycetes','Orbiliomycetes',
                       'Arthoniomycetes','Geoglossomycetes','Xylonomycetes','Coniocybomycetes',
                       'Taphrinomycetes,Saccharomycetes,Pneumocystidomycetes,Schizosaccharomycetes,unclassified'),
         !Site == 'Butterfly trail')  -> fe.cf.plot

# Plot by range, arranged based on BIO17 (low to high)
fe.cf.plot$range <- factor(fe.cf.plot$range,
                           levels = c('Huachuca','Chiricahua','Santa.Catalina',
                                      'Pinaleno','Bradshaw','MogollonRim','Mingus'))

fe.cf.range.plot <- ggplot(data = fe.cf.plot,
                           aes(x = range,
                               fill = Class.TBAS)) +
  geom_bar(position = "fill") +
  ylab("Proportion of \n sequences per class") +
  theme_bw() +
  xlab(element_blank()) +
  scale_fill_manual(values = c('#d8b365','#f6e8c3','#c7eae5','#5ab4ac','#01665e')) +
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

ggsave('Fig6E_FE_CF_Class_rangeBIO17.pdf', plot = fe.cf.range.plot, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)
