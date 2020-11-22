## Script created by Liz Bowman January 6, 2020
## for analysing community data in Sky Island study
## Includes endophyte data

## For answering 'What environmental variables are important in structuring FE communities across Arizona?'
# Part A focuses on a PERMANOVA
# 1. Are there differences between the ranges?
# 2. What contributes to these differences? 
#   Here I will assess various climate variables, habitat type, etc...

#========================================================================================#
# Load libraries and data------------------
#========================================================================================#

#install.packages(c('dplyr','tidyr','vegan','ggplot2'))
library(dplyr);library(tidyr);library(vegan);library(ggplot2);library(nlme)

#--file paths
dat.dir <- 'data/'
fig.dir <- 'figures/'
res.dir <- 'results/'

#--Site x species matrix ITS2 rarefied
fe.site <- read.csv(paste0(dat.dir,'FE_ITS2rarified_sitexspecies.csv'), as.is = T)
#--environmental data
env.data <- read.csv(paste0(dat.dir,'FE_CF_EnvFactors_site.csv'), as.is = T)
#--diversity
fe.cf.div <- read.csv(paste0(dat.dir,'FE_CF_diversity.csv'), as.is = T)

#--Site x species matrix culture-based data
fe.cb.site <- read.csv(paste0(dat.dir,'FE_CB_sitexspecies.csv'), as.is = T)
#--environmental data
env.cb.site <- read.csv(paste0(dat.dir,'FE_CB_EnvFactors_site.csv'), as.is = T)
#--diversity
fe.cb.div <- read.csv(paste0(dat.dir,'FE_CB_diversity.csv'), as.is = T)

#--Results table
fe.permanova.results <- data.frame(terms = c('BIO17','BIO12','forest',
                                             'BIO17:BIO12','BIO17:forest','BIO12:forest',
                                             'BIO17:BIO12:forest','Residuals','Total'),
                                   DF = NA,
                                   Sum.sq = NA,
                                   Mean.sq = NA,
                                   F.model = NA,
                                   R2 = NA,
                                   P = NA)

#========================================================================================#
# Culture-free FE data------------------
#========================================================================================#

#----------------------------------------------------------------#
# Jaccard dissimilarity
#----------------------------------------------------------------#

#--Remove outliers M21
#fe.site <- fe.site[!fe.site$site %in% c('M21', 'B12n'),]

#--isolate otu data
comm.matrix <- fe.site[15:length(fe.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- fe.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)
jaccard.otu$stress

##--NMDS by forest type and precipitation of the driest quarter----
pdf(file = paste0(fig.dir, 'Fig6C_FE_CF_NMDSbyEnv.pdf'),
    width = 12.5, height = 8.5)

layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
colfunc <- colorRampPalette(c("white","black"))
par(mar=c(5,6,4,4)+.1)

plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
# #--color for range
# color.vec <- data.frame(color = rep(NA,length(rownames(comm.matrix))),
#                          b.group = otu.jaccard$range)
# #--color by range
# color.vec$color <- sapply(1:nrow(color.vec), function (x)
#   if(color.vec[x,'b.group'] == 'Bradshaw') {color.vec[x,'color'] = 'black'}
#   else if(color.vec[x,'b.group'] == 'Chiricahua') {color.vec[x,'color'] = 'grey'}
#   else if(color.vec[x,'b.group'] == 'MogollonRim') {color.vec[x,'color'] = 'blue'}
#   else if(color.vec[x,'b.group'] == 'Mingus') {color.vec[x,'color'] = 'green'}
#   else if(color.vec[x,'b.group'] == 'Santa.Catalina') {color.vec[x,'color'] = 'orange'}
#   else if(color.vec[x,'b.group'] == 'Pinaleno') {color.vec[x,'color'] = 'purple'}
#   else {color.vec[x,'color'] = 'brown'})

#--color gradient for elevation
colfunc <- colorRampPalette(c('black','lightgrey'))
otu.jaccard$color <- colfunc(5)[as.numeric(cut(otu.jaccard$BIO17,breaks = 5))]

#--shape for forest type
shape.vec <- data.frame(shape = rep(NA,length(rownames(comm.matrix))),
                        b.group = otu.jaccard$forest)
#--shape by forest type
shape.vec$shape <- sapply(1:nrow(shape.vec), function (x)
  if(shape.vec[x,'b.group'] == 'pine') {shape.vec[x,'shape'] = 16}
  else {shape.vec[x,'shape'] = 1})

points(jaccard.otu, display = "sites", cex = 3, lwd = 4,
       pch = shape.vec$shape,
       col = otu.jaccard$color,
       bg = otu.jaccard$color)
#ordilabel(jaccard.otu, display = 'sites', labels = otu.jaccard$sample)
# fe.fit <- envfit(jaccard.otu ~ forest * BIO17 * BIO12, data = otu.jaccard)
# plot(fe.fit, labels = list(factors = c('',''),
#                         vectors = c('', '')),
#      col = 'black', cex = 2)

# legend('topleft', legend = c('Pine','Pine-Douglas fir'), pch = c(16,1),
#        bty = 'n', cex = 2)

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

mtext(c(42,'','','',64),
      side = 2, at=tail(seq(yb,yt,(yt-yb)/5),-1)-0.05, las=2, cex=1.5)

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$range)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

#--PERMANOVA
jaccard.adonis <- adonis(comm.dist.jaccard ~ BIO17 * BIO12 * forest, data = otu.jaccard)
jaccard.adonis

# results table
fe.permanova.cf.j <- fe.permanova.results
fe.permanova.cf.j['DF'] <- jaccard.adonis$aov.tab$Df
fe.permanova.cf.j['Sum.sq'] <- jaccard.adonis$aov.tab$SumsOfSqs
fe.permanova.cf.j['Mean.sq'] <- jaccard.adonis$aov.tab$MeanSqs
fe.permanova.cf.j['F.model'] <- jaccard.adonis$aov.tab$F.Model
fe.permanova.cf.j['R2'] <- jaccard.adonis$aov.tab$R2
fe.permanova.cf.j['P'] <- jaccard.adonis$aov.tab$`Pr(>F)`

write.csv(fe.permanova.cf.j, paste0(res.dir,'FE_CF_PermanovaJaccard.csv'),
          row.names = F)

#----------------------------------------------------------------#
# Species richness-----
#----------------------------------------------------------------#

#--Additive
fe.cf.lm.a <- lm(spec.richness ~ BIO12 + BIO17 + forest, data = fe.cf.div)
summary(fe.cf.lm.a)
anova(fe.cf.lm.a)

#--Interaction
fe.cf.lm.b <- lm(spec.richness ~ BIO12 * BIO17 * forest, data = fe.cf.div)
summary(fe.cf.lm.b)
anova(fe.cf.lm.b)

cf.specrich <- ggplot(data = fe.cf.div,
       aes(x = BIO17,
           y = spec.richness)) +
  geom_point(size = 3) +
  #facet_grid(. ~ forest) +
  geom_smooth(se = F, method = 'lm', color = 'darkgrey') +
  ylab("Species Richness") +
  theme_bw() +
  xlab('Mean precipitation\nduring the driest quarter (mm)') +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        plot.margin=unit(c(1,1.2,1,1),"cm"),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig6A_FE_CF_speciesrichness.pdf', plot = cf.specrich, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#========================================================================================#
# Culture-based FE data------------------
#========================================================================================#

#----------------------------------------------------------------#
# Jaccard dissimilarity
#----------------------------------------------------------------#

#--Remove outliers M21
#fe.site <- fe.site[!fe.site$site %in% c('M21', 'B12n'),]

#--isolate otu data
comm.matrix <- fe.cb.site[12:length(fe.cb.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- fe.cb.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)
jaccard.otu$stress

#--NMDS by forest type and precipitation of the driest quarter----
pdf(file = paste0(fig.dir, 'Fig6D_FE_CB_Jaccard_NMDSbyEnv.pdf'),
    width = 12.5, height = 8.5)

layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
colfunc <- colorRampPalette(c("white","black"))
par(mar=c(5,6,4,4)+.1)

plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
# #--color for range
# color.vec <- data.frame(color = rep(NA,length(rownames(comm.matrix))),
#                          b.group = otu.jaccard$range)
# #--color by range
# color.vec$color <- sapply(1:nrow(color.vec), function (x)
#   if(color.vec[x,'b.group'] == 'Bradshaw') {color.vec[x,'color'] = 'black'}
#   else if(color.vec[x,'b.group'] == 'Chiricahua') {color.vec[x,'color'] = 'grey'}
#   else if(color.vec[x,'b.group'] == 'MogollonRim') {color.vec[x,'color'] = 'blue'}
#   else if(color.vec[x,'b.group'] == 'Mingus') {color.vec[x,'color'] = 'green'}
#   else if(color.vec[x,'b.group'] == 'Santa.Catalina') {color.vec[x,'color'] = 'orange'}
#   else if(color.vec[x,'b.group'] == 'Pinaleno') {color.vec[x,'color'] = 'purple'}
#   else {color.vec[x,'color'] = 'brown'})

#--color gradient for elevation
colfunc <- colorRampPalette(c('black','lightgrey'))
otu.jaccard$color <- colfunc(5)[as.numeric(cut(otu.jaccard$BIO17,breaks = 5))]

#--shape for forest type
shape.vec <- data.frame(shape = rep(NA,length(rownames(comm.matrix))),
                        b.group = otu.jaccard$forest)
#--shape by forest type
shape.vec$shape <- sapply(1:nrow(shape.vec), function (x)
  if(shape.vec[x,'b.group'] == 'pine') {shape.vec[x,'shape'] = 16}
  else {shape.vec[x,'shape'] = 1})

points(jaccard.otu, display = "sites", cex = 3, lwd = 4,
       pch = shape.vec$shape,
       col = otu.jaccard$color,
       bg = otu.jaccard$color)
#ordilabel(jaccard.otu, display = 'sites', labels = otu.jaccard$sample)
# fe.fit <- envfit(jaccard.otu ~ forest + BIO17 + BIO12, data = otu.jaccard)
# plot(fe.fit, col = 'black', cex = 2,
#      labels = list(factors = c('',''),
#                    vectors = c('', '')))

# legend('topleft', legend = c('Pine','Pine-Douglas fir'), pch = c(16,1),
#        bty = 'n', cex = 2)

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

mtext(c(42,'','','',64),
      side = 2, at=tail(seq(yb,yt,(yt-yb)/5),-1)-0.05, las=2, cex=1.5)

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$range)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

#--PERMANOVA
jaccard.adonis <- adonis(comm.dist.jaccard ~ BIO17 * BIO12 * forest, data = otu.jaccard)
jaccard.adonis

# results table
fe.permanova.cb.j <- fe.permanova.results
fe.permanova.cb.j['DF'] <- jaccard.adonis$aov.tab$Df
fe.permanova.cb.j['Sum.sq'] <- jaccard.adonis$aov.tab$SumsOfSqs
fe.permanova.cb.j['Mean.sq'] <- jaccard.adonis$aov.tab$MeanSqs
fe.permanova.cb.j['F.model'] <- jaccard.adonis$aov.tab$F.Model
fe.permanova.cb.j['R2'] <- jaccard.adonis$aov.tab$R2
fe.permanova.cb.j['P'] <- jaccard.adonis$aov.tab$`Pr(>F)`

write.csv(fe.permanova.cb.j, paste0(res.dir,'FE_CB_PermanovaJaccard.csv'),
          row.names = F)

#----------------------------------------------------------------#
# Culture-based: Morisita-Horn dissimilarity
#----------------------------------------------------------------#

#--Remove outliers M21
#fe.site <- fe.site[!fe.site$site %in% c('M21', 'B12n'),]

#--isolate otu data
comm.matrix <- fe.cb.site[12:length(fe.cb.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.horn <- fe.cb.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.horn<- vegdist(comm.matrix, method = "horn")

#--NMDS analysis
horn.otu <- metaMDS(comm.dist.horn, dist = "bray",
                       try = 100, trymax = 100)
horn.otu$stress

#--NMDS by forest type and precipitation of the driest quarter
pdf(file = paste0(fig.dir, 'FigS3B_FE_CB_MorisitaHorn_NMDSbyEnv.pdf'),
    width = 12.5, height = 8.5)

layout(matrix(1:2,nrow=1),widths=c(0.8,0.2))
colfunc <- colorRampPalette(c("white","black"))
par(mar=c(5,6,4,4)+.1)

plot(horn.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
# #--color for range
# color.vec <- data.frame(color = rep(NA,length(rownames(comm.matrix))),
#                          b.group = otu.jaccard$range)
# #--color by range
# color.vec$color <- sapply(1:nrow(color.vec), function (x)
#   if(color.vec[x,'b.group'] == 'Bradshaw') {color.vec[x,'color'] = 'black'}
#   else if(color.vec[x,'b.group'] == 'Chiricahua') {color.vec[x,'color'] = 'grey'}
#   else if(color.vec[x,'b.group'] == 'MogollonRim') {color.vec[x,'color'] = 'blue'}
#   else if(color.vec[x,'b.group'] == 'Mingus') {color.vec[x,'color'] = 'green'}
#   else if(color.vec[x,'b.group'] == 'Santa.Catalina') {color.vec[x,'color'] = 'orange'}
#   else if(color.vec[x,'b.group'] == 'Pinaleno') {color.vec[x,'color'] = 'purple'}
#   else {color.vec[x,'color'] = 'brown'})

#--color gradient for elevation
colfunc <- colorRampPalette(c('black','lightgrey'))
otu.horn$color <- colfunc(5)[as.numeric(cut(otu.horn$BIO17,breaks = 5))]

#--shape for forest type
shape.vec <- data.frame(shape = rep(NA,length(rownames(comm.matrix))),
                        b.group = otu.horn$forest)
#--shape by forest type
shape.vec$shape <- sapply(1:nrow(shape.vec), function (x)
  if(shape.vec[x,'b.group'] == 'pine') {shape.vec[x,'shape'] = 16}
  else {shape.vec[x,'shape'] = 1})

points(horn.otu, display = "sites", cex = 3, lwd = 4,
       pch = shape.vec$shape,
       col = otu.horn$color,
       bg = otu.horn$color)
#ordilabel(jaccard.otu, display = 'sites', labels = otu.jaccard$sample)
# fe.fit <- envfit(jaccard.otu ~ forest + BIO17 +, data = otu.jaccard)
# plot(fe.fit, col = 'black', cex = 2,
#      labels = list(factors = c('',''),
#                    vectors = c('', '')))
# 
# legend('topright', legend = c('Pine','Pine-Douglas fir'), pch = c(16,1),
#        bty = 'n', cex = 2)

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

mtext(c(43,'','','',64),
      side = 2, at=tail(seq(yb,yt,(yt-yb)/5),-1)-0.05, las=2, cex=1.5)

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.horn, group = otu.horn$range)
horn.betadisper <- anova(betadisp)
horn.betadisper

#--PERMANOVA
horn.adonis <- adonis(comm.dist.horn ~ BIO17 * BIO12 * forest, data = otu.horn)
horn.adonis

# results table
fe.permanova.cb.mh <- fe.permanova.results
fe.permanova.cb.mh['DF'] <- horn.adonis$aov.tab$Df
fe.permanova.cb.mh['Sum.sq'] <- horn.adonis$aov.tab$SumsOfSqs
fe.permanova.cb.mh['Mean.sq'] <- horn.adonis$aov.tab$MeanSqs
fe.permanova.cb.mh['F.model'] <- horn.adonis$aov.tab$F.Model
fe.permanova.cb.mh['R2'] <- horn.adonis$aov.tab$R2
fe.permanova.cb.mh['P'] <- horn.adonis$aov.tab$`Pr(>F)`

write.csv(fe.permanova.cb.mh, paste0(res.dir,'FE_CB_PermanovaMorisitaHorn.csv'),
          row.names = F)

#----------------------------------------------------------------#
# Fisher's alpha---
#----------------------------------------------------------------#

# remove outlier, P2 and M1
fe.cb.div.out <- fe.cb.div[!fe.cb.div$site %in% c('P2','M1'),]

# log transform
fe.cb.div.out$log.fa <- log(fe.cb.div.out$fa)

#--Multiple linear regression: Fisher's alpha
fe.cb.lm.a <- lm(log.fa ~ BIO12 * BIO17 * forest, data = fe.cb.div.out)
summary(fe.cb.lm.a)
anova(fe.cb.lm.a)

#--Multiple linear regression: Fisher's alpha additive
fe.cb.lm.b <- lm(log.fa ~ BIO12 + BIO17 + forest, data = fe.cb.div.out)
summary(fe.cb.lm.b)
anova(fe.cb.lm.b)

anova(fe.cb.lm.a,fe.cb.lm.b)

#--Plot
cb.fisher <- ggplot(data = fe.cb.div.out,
       aes(x = BIO17,
           y = log.fa)) +
  geom_point(size = 3) +
  #geom_smooth(se = F, method = 'lm', color = 'darkgrey') +
  ylab("Fisher's alpha") +
  theme_bw() +
  xlab('Mean precipitation\nduring the driest quarter (mm)') +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        plot.margin=unit(c(1,1.2,1,1),"cm"),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('FigS3A_FE_CB_FishersAlpha.pdf', plot = cb.fisher, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#----------------------------------------------------------------#
# Species richness----
#----------------------------------------------------------------#

#--Multiple linear regression: Fisher's alpha
fe.cb.lm <- lm(spec.richness.all ~ BIO12 + BIO17 + forest, data = fe.cb.div)
summary(fe.cb.lm)
anova(fe.cb.lm)

cb.sr <- ggplot(data = fe.cb.div,
       aes(x = BIO17,
           y = spec.richness.all)) +
  geom_point(size = 3) +
  #geom_smooth(se = F, method = 'lm', color = 'darkgrey') +
  ylab("Species richness") +
  theme_bw() +
  xlab('Mean precipitation of the\ndriest quarter (mm)') +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        plot.margin=unit(c(1,1.2,1,1),"cm"),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig6B_FE_CB_SpecRichness.pdf', plot = cb.sr, 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)
