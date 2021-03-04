## Script created by Liz Bowman January 15, 2020
## for analysing community data in Sky Island study
## Includes ectomycorrhizal data only

## For answering 'To what extent does distance to the closest P. ponderosa forest
# structure abundance, diversity, and communities of ectomycorrhizal fungi as opposed 
# to environment'?
# Part A on community composition.
# Part B on diversity.
# Part C on diversity (These analyses were not included in the paper, but are included here for transparency).

#========================================================================================#
# Load libraries and data------------------
#========================================================================================#

#--Site x species matrix with tip abundance per OTU
em.site <- read.csv(paste0(dat.dir,'EM_sitexspecies_TipAb_Site.csv'), as.is = T)
#--Environmental data
env.data <- read.csv(paste0(dat.dir,'EM_EnvFactors_site.csv'), as.is = T)

#--Distance to closest Ponderosa pine forest
dist.data <- read.csv(paste0(dat.dir,'DistanceClosestPPForest_site.csv'), as.is = T)

#--Diversity data
em.div <- read.csv(paste0(dat.dir, 'EM_Diversity_site.csv'), as.is = T)

# #--Abundance data
# em.abund <- read.csv(paste0(dat.dir, 'EM_Abundance_site.csv'), as.is = T)
# em.abund.tree <- read.csv(paste0(dat.dir, 'EM_Abundance_tree.csv'), as.is = T)

#--Results table
em.permanova.results <- data.frame(terms = c('dist.closest.forest','BIO12','forest',
                                             'dist.BIO12','dist.forest','BIO12.forest',
                                             'dist.BIO12.forest','Residuals','Total'),
                                   DF = NA,
                                   Sum.sq = NA,
                                   Mean.sq = NA,
                                   F.model = NA,
                                   R2 = NA,
                                   P = NA)

#========================================================================================#
# Part A: Ordination and NMDS with distance to closest Ponderosa pine forest edge as color----
#========================================================================================#
#----------------------------------------------------------------#
# << Jaccard dissimilarity >> ----
#----------------------------------------------------------------#

#--isolate otu data
comm.matrix <- em.site[17:length(em.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- em.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)
jaccard.otu$stress

##--Distance from closest Ponderosa forest----
pdf(file = paste0(fig.dir, 'FigS2_EM_NMDS_DistanceClosestPonderosa_jaccard.pdf'),
    width = 12.5, height = 8.5)

layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
colfunc <- colorRampPalette(c("white","black"))
par(mar=c(5,6,4,4)+.1)

plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
#--color for groups
# color.vec <- data.frame (color = rep(NA,length(rownames(comm.matrix))),
#                          p.group = nmds.ab$elevation)
# color.vec <- sapply(color.vec$p.group, function (x)
#   if (x == 2425) {color.vec = 'black'}
#   else if (x == 2370) {color.vec = 'grey'}
#   else if (x == 2352) {color.vec = 'aquamarine3'}
#   else if(x == 2343) {color.vec = 'darkgreen'}
#   else if(x == 2201) {color.vec = 'darkorange3'}
#   else if(x == 2170) {color.vec = 'darkslategrey'}
#   else if(x == 2119) {color.vec = 'gold1'}
#   else if(x == 2421) {color.vec = 'darkmagenta'}
#   else if (x == 1790) {color.vec = 'brown4'}
#   else{'yellow2'})

#--color gradient for elevation
colfunc <- colorRampPalette(c('darkblue','lightgrey'))
otu.jaccard$color <- colfunc(5)[as.numeric(cut(otu.jaccard$dist.closest.forest,breaks = 5))]

#--Plant community point shapes
shape.df <- data.frame(shape = rep(NA,length(rownames(comm.matrix))),
                    p.group = otu.jaccard$forest.type)
shape.df$shape <- sapply(shape.df$p.group, function (x)
  if (x == 'pine') {shape = 21}
  else if (x == 'pine-oak') {shape = 22}
  else if (x == 'pine-Doug fir') {shape = 1}
  else{18})
# env.fit <- envfit(jaccard.otu ~ BIO12 +dist.closest.forest,
#                   data = otu.jaccard)
# plot(env.fit, col = 'black',
#      labels = c('',''), lty = 2, lwd = 3)
points(jaccard.otu, display = "sites",
       cex = 3, lwd = 4,
       pch = shape.df$shape,
       col = otu.jaccard$color,
       bg = otu.jaccard$color)


# legend('topleft', legend = c('Pine','Pine-oak','Pine-Douglas fir'),
#        pch = c(16,15,1), bty = 'n', cex = 2, bg = 'black', col = 'black')

#--add color gradient legend
xl <- 1; yb <- 1; xr <- 1.5; yt <- 2
par(mar=c(4,0.5,4.1,0.5))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/5),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/5),-1),
  col=colfunc(5)
)

mtext(c(0,15,30,45,90),
      side = 2, at=tail(seq(yb,yt,(yt-yb)/5),-1)-0.05, las=2, cex=1.5)

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$dist.closest.forest)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

# PERMANOVA: just distance to closest PP forest ----------------------------------------
jaccard.dist.adonis <- adonis(comm.dist.jaccard ~ dist.closest.forest,
                              data = otu.jaccard)
jaccard.dist.adonis

# PERMANOVA: dist., BIO12, plant community ----------------------------------------
jaccard.dist.adonis <- adonis(comm.dist.jaccard ~ dist.closest.forest * BIO12 * forest.type,
                              data = otu.jaccard)
jaccard.dist.adonis

# results table
em.permanova.jaccard <- em.permanova.results
em.permanova.jaccard['DF'] <- jaccard.dist.adonis$aov.tab$Df
em.permanova.jaccard['Sum.sq'] <- jaccard.dist.adonis$aov.tab$SumsOfSqs
em.permanova.jaccard['Mean.sq'] <- jaccard.dist.adonis$aov.tab$MeanSqs
em.permanova.jaccard['F.model'] <- jaccard.dist.adonis$aov.tab$F.Model
em.permanova.jaccard['R2'] <- jaccard.dist.adonis$aov.tab$R2
em.permanova.jaccard['P'] <- jaccard.dist.adonis$aov.tab$`Pr(>F)`

# write.csv(em.permanova.jaccard,paste0(res.dir,'EM_PermanovaJaccard_Results.csv'),
#           row.names = F)

#----------------------------------------------------------------#
# << Morisita-horn dissimilarity >> ----
#----------------------------------------------------------------#

#--Remove outlier, comment to keep in
em.site <- em.site[-11,]

#--isolate otu data
comm.matrix <- em.site[17:length(em.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.horn <- em.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.horn <- vegdist(comm.matrix, method = "horn", binary = TRUE)

#--NMDS analysis
horn.otu <- metaMDS(comm.dist.horn, dist = "bray",
                       try = 100, trymax = 100)
horn.otu$stress

##--Distance from closest Ponderosa forest----
pdf(file = paste0(fig.dir, 'Fig4B_EM_NMDS_DistanceClosestPonderosa_MorisitaHorn.pdf'),
    width = 12.5, height = 8.5)

layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
colfunc <- colorRampPalette(c("white","black"))
par(mar=c(5,6,4,4)+.1)

plot(horn.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
#--color for groups
# color.vec <- data.frame (color = rep(NA,length(rownames(comm.matrix))),
#                          p.group = nmds.ab$elevation)
# color.vec <- sapply(color.vec$p.group, function (x)
#   if (x == 2425) {color.vec = 'black'}
#   else if (x == 2370) {color.vec = 'grey'}
#   else if (x == 2352) {color.vec = 'aquamarine3'}
#   else if(x == 2343) {color.vec = 'darkgreen'}
#   else if(x == 2201) {color.vec = 'darkorange3'}
#   else if(x == 2170) {color.vec = 'darkslategrey'}
#   else if(x == 2119) {color.vec = 'gold1'}
#   else if(x == 2421) {color.vec = 'darkmagenta'}
#   else if (x == 1790) {color.vec = 'brown4'}
#   else{'yellow2'})

#--color gradient for elevation
colfunc <- colorRampPalette(c('darkblue','lightgrey'))
otu.horn$color <- colfunc(5)[as.numeric(cut(otu.horn$dist.closest.forest,breaks = 5))]

#--Plant community point shapes
shape.df <- data.frame(shape = rep(NA,length(rownames(comm.matrix))),
                       p.group = otu.horn$forest.type)
shape.df$shape <- sapply(shape.df$p.group, function (x)
  if (x == 'pine') {shape = 21}
  else if (x == 'pine-oak') {shape = 22}
  else if (x == 'pine-Doug fir') {shape = 1}
  else{18})
points(horn.otu, display = "sites",
       cex = 3.5, lwd = 4,
       pch = shape.df$shape,
       col = otu.horn$color,
       bg = otu.horn$color)
env.fit <- envfit(jaccard.otu ~ BIO12 + dist.closest.forest,
                  data = otu.jaccard)
plot(env.fit, col = 'black',
     labels = c('',''),
     cex = 2)

# legend('bottomright', legend = c('Pine','Pine-oak','Pine-Douglas fir'),
#        pch = c(16,15,1), bty = 'n', cex = 2)

#--add color gradient legend
xl <- 1; yb <- 1; xr <- 1.5; yt <- 2
par(mar=c(4,0.5,4.1,0.5))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(
  xl,
  head(seq(yb,yt,(yt-yb)/5),-1),
  xr,
  tail(seq(yb,yt,(yt-yb)/5),-1),
  col=colfunc(5)
)

mtext(c(0,15,30,45,90),
      side = 2, at=tail(seq(yb,yt,(yt-yb)/5),-1)-0.05, las=2, cex=1.5)

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.horn, group = otu.horn$dist.closest.forest)
horn.betadisper <- anova(betadisp)
horn.betadisper

# PERMANOVA ----------------------------------------
horn.dist.adonis <- adonis(comm.dist.horn ~ dist.closest.forest, data = otu.horn)
horn.dist.adonis

# PERMANOVA: environment ----------------------------------------
horn.dist.adonis <- adonis(comm.dist.horn ~ dist.closest.forest * BIO12 * forest.type,
                           data = otu.horn)
horn.dist.adonis

# results table
em.permanova.horn <- em.permanova.results
em.permanova.horn['DF'] <- horn.dist.adonis$aov.tab$Df
em.permanova.horn['Sum.sq'] <- horn.dist.adonis$aov.tab$SumsOfSqs
em.permanova.horn['Mean.sq'] <- horn.dist.adonis$aov.tab$MeanSqs
em.permanova.horn['F.model'] <- horn.dist.adonis$aov.tab$F.Model
em.permanova.horn['R2'] <- horn.dist.adonis$aov.tab$R2
em.permanova.horn['P'] <- horn.dist.adonis$aov.tab$`Pr(>F)`

# write.csv(em.permanova.horn,paste0(res.dir,'EM_PermanovaMorisitaHorn_Results.csv'),
#            row.names = F)

#========================================================================================#
# Part B: Diversity as a function of distance to closest Ponderosa pine forest-----
#========================================================================================#

em.div$forest.type <- factor(em.div$forest.type)

## Overall
# multiple linear regression
# Distance to closest Ponderosa Pine forest
lm.em.dist <- lm(log.fa ~ dist.closest.forest, data = em.div)
summary(lm.em.dist)

# Climate
lm.em.prec <- lm(log.fa ~ BIO1 + BIO12, data = em.div)
summary(lm.em.prec)

# multiple linear regression
lm.em.mult <- lm(log.fa ~ dist.closest.forest + BIO12 + forest.type, data = em.div)
all.lm <- summary(lm.em.mult)
anova.all <- anova(lm.em.mult)

# Assess models
anova(lm.em.dist, lm.em.prec, lm.em.mult)

# results table
div.anova.results <- data.frame(terms = c('dist.closest.forest','BIO12','forest.type','Residuals'))
div.anova.results['DF'] <- anova.all$Df
div.anova.results['Sum.sq'] <- anova.all$`Sum Sq`
div.anova.results['Mean.sq'] <- anova.all$`Mean Sq`
div.anova.results['F.value'] <- anova.all$`F value`
div.anova.results['P'] <- anova.all$`Pr(>F)`

# write.csv(div.anova.results, paste0(res.dir,'EM_DiversityANOVA_Results.csv'),
#           row.names = F)

#--Plot
em.fishers <- ggplot(em.div, aes(x = BIO12,
                   y = log.fa)) +
  geom_point(size = 4) +
  #facet_grid(. ~ forest.type) +
  geom_smooth(se = F, method = 'lm', color = 'darkgrey') +
  xlab('Mean annual precipitation (mm)') +
  ylab("Log Fisher's alpha") +
  theme_classic() +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.spacing = unit(2, "lines"),
        plot.margin=unit(c(1,1.2,1,1),"cm"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig4A_EM_FishersAlpha_Prec.pdf', plot = em.fishers, 
       device = 'pdf', path = fig.dir,
       units = 'in', width = 10, height = 7)

## With Mogollon Rim removed
em.div.nomr <- em.div[!em.div$range == 'Mogollon Rim',]
# multiple linear regression
# Distance to closest Ponderosa Pine forest
lm.em.dist <- lm(log.fa ~ dist.closest.forest, data = em.div.nomr)
summary(lm.em.dist)

# Climate
lm.em.prec <- lm(log.fa ~ BIO1 + BIO12, data = em.div.nomr)
summary(lm.em.prec)

# multiple linear regression
lm.em.mult <- lm(log.fa ~ dist.closest.forest + BIO12 + forest.type, data = em.div.nomr)
all.lm <- summary(lm.em.mult)
anova.all <- anova(lm.em.mult)

# Assess models
anova(lm.em.dist, lm.em.prec, lm.em.mult)

## Only Mogollon Rim 
em.div.mr <- em.div[em.div$range == 'Mogollon Rim',]
# multiple linear regression
# Distance to closest Ponderosa Pine forest
lm.em.dist <- lm(log.fa ~ dist.closest.forest, data = em.div.mr)
summary(lm.em.dist)

# Climate
lm.em.prec <- lm(log.fa ~ BIO1 + BIO12, data = em.div.mr)
summary(lm.em.prec)


#========================================================================================#
# Part C: Abundance as a function of distance and environment -----
#========================================================================================#

## Site data
# Remove outliers and log transform
# em.abund$log.total.tip <- log(em.abund$total.tip)
# em.abund.out <- em.abund[em.abund$Site %in% c('H1','H5','H6','M1'),]
# 
# # multiple linear regression
# lm.em.mult <- lm(log.total.tip ~ dist.closest.forest + BIO12 + forest.type,
#                  data = em.abund)
# all.lm <- summary(lm.em.mult)
# anova.all <- anova(lm.em.mult)
# 
# ## Tree data
# # Overall 
# # Remove outliers and log transform
# em.abund.tree$log.total.tip <- log(em.abund.tree$total.tip)
# em.abund.tree.out <- em.abund.tree[em.abund.tree$total.tip > 11,]
# 
# # multiple linear regression
# lm.em.mult <- lm(log.total.tip ~ dist.closest.forest + BIO12 + forest,
#                  data = em.abund.tree.out)
# all.lm <- summary(lm.em.mult)
# anova.all <- anova(lm.em.mult)
# 
# # Plot
# em.abund <- ggplot(em.abund.tree.out, aes(x = BIO12,
#                                           y = log.total.tip)) +
#   geom_point(size = 4) +
#   #facet_grid(. ~ forest.type) +
#   geom_smooth(se = F, method = 'lm', color = 'darkgrey') +
#   xlab('Mean annual precipitation (mm)') +
#   ylab("Log EM abundance") +
#   theme_classic() +
#   theme(legend.position='right',
#         # axis.title.x = element_text(margin = margin(t = 30)),
#         axis.title.y = element_text(margin = margin(r = 30)),
#         #axis.text.x = element_text(angle = 35),
#         panel.spacing = unit(2, "lines"),
#         plot.margin=unit(c(1,1.2,1,1),"cm"),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.border = element_blank(), axis.line = element_line(),
#         axis.text = element_text(size=20, color = 'black'),
#         axis.title = element_text(size = 28),
#         strip.text.x = element_text(size = 14))
# 
# ggsave('EM_Abundance_Prec.pdf', plot = em.abund, 
#        device = 'pdf', path = fig.dir,
#        units = 'in', width = 10, height = 7)
# 
# ## With Mogollon rim removed
# em.abund.tree.out.nomr <- em.abund.tree.out[!em.abund.tree.out$Range == 'MogollonRim',]
# # multiple linear regression
# lm.em.mult <- lm(log.total.tip ~ dist.closest.forest + BIO12 + forest,
#                  data = em.abund.tree.out.nomr)
# all.lm <- summary(lm.em.mult)
# anova.all <- anova(lm.em.mult)
# 
# ## Only Mogollon rim 
# em.abund.tree.out.mr <- em.abund.tree.out[em.abund.tree.out$Range == 'MogollonRim',]
# # multiple linear regression
# lm.em.mult <- lm(log.total.tip ~  BIO12 + forest,
#                  data = em.abund.tree.out.mr)
# all.lm <- summary(lm.em.mult)
# anova.all <- anova(lm.em.mult)

