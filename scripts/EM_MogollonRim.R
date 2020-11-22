## Script created by Liz Bowman September 24, 2019
## for analysing community data in Sky Island study
## Includes both ectomycorrhizal data

#========================================================================================#
# Load libraries and data for abundance------------------
#========================================================================================#

#install.packages(c('dplyr','tidyr','vegan','ggplot2'))
library(dplyr);library(tidyr);library(vegan);library(ggplot2);library(nlme)

#--file paths
dat.dir <- '~/Documents/PhD/3:4_Combined/data/'
fig.dir <- '~/Documents/PhD/3:4_Combined/figures/'
res.dir <- '~/Documents/PhD/3:4_Combined/results/'

#--Site x species matrix with tip abundance per OTU
em.site <- read.csv(paste0(dat.dir,'EM_sitexspecies_TipAb_Site.csv'), as.is = T)
#isolate Mogollon Rim data
em.site <- em.site[em.site$Range == 'MogollonRim',]
#--Site x species matrix ITS2 rarefied
fe.tree <- read.csv(paste0(dat.dir,'FE_ITS2rarified_treexspecies.csv'), as.is = T)
#isolate Mogollon Rim data
fe.tree <- fe.tree[fe.tree$range == 'MogollonRim',]

#--Distance decay data for EM fungi
dist.decay.em <- read.csv(paste0(dat.dir, 'EM_DistanceDecayMatrix.csv'), as.is = T)

#--Make results table for culture-based data
results.em <- data.frame(measure = c('Jaccard.em', 'Horn.em',
                                     'Jaccard.fe'))

#=========================================================================================#
# Jaccard: Ectomycorrhizal fungi 97%----------
#=========================================================================================#

#--Remove outliers
# otu.jaccard <- otu.jaccard[!otu.jaccard$Tree %in% c('N35','H21','M62',
#                                                     'A32','N15','M52','H15'),]
#--isolate otu data
comm.matrix <- em.site[10:length(em.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- em.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)

#--Stress
results.em[results.em$measure == 'Jaccard.em', 'stress'] <- jaccard.otu$stress

#--GGplot plot
# Create dataframe
data.scores <- as.data.frame(scores(jaccard.otu))
data.scores$site <- rep(c('central','south','north'),each = 3)

# outliers row N15, H21
#data.scores <- data.scores[!data.scores$tree %in% c('N15','H21'),]

#--GGplot: NMDS
em.jaccard <- ggplot() + 
  geom_point(data = data.scores,aes(x = NMDS1,
                                    y = NMDS2,
                                    color=site),size=2) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=28,margin = margin(t = 30)), # remove x-axis labels
        axis.title.y = element_text(size=28,margin = margin(r = 30)), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

em.jaccard

#--Base R plot:
#--Range----
#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'NDMS_MogollonRimJaccard_EM.jpeg'),
     width = 1200, height = 1000,
     quality = 100)

plot(jaccard.otu, display = "sites", cex.lab = 1.5,
     cex.axis = 1.5)
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(otu.jaccard)),
                        p.group = otu.jaccard$Site)

color.vec[color.vec$p.group == 'M1', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M2', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M3', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M4', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M5', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M6', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M7', 'color'] <- 'darkgreen'
color.vec[color.vec$p.group == 'M8', 'color'] <- 'darkgreen'
color.vec[color.vec$p.group == 'M9', 'color'] <- 'darkgreen'
color.vec <- color.vec[order(color.vec$p.group),]

# ordipointlabel(jaccard.otu, display = "sites")

points(jaccard.otu, display = "sites", cex = 4,
       pch = 16,
       col = color.vec$color,
       bg = color.vec$color)

legend("topright", legend = c('North','Central','South'), bty = "n",
       col = c('darkgreen','grey','darkblue'),
       pch = 16, cex = 3)

#--Ordihull variations by range, burn, and both burn and range
ordihull(jaccard.otu, groups = color.vec$color,
         col = 'black')
#--Overlays
fit <- envfit(jaccard.otu ~ prec + Tavg, otu.jaccard)
fit
plot(fit, type = 'n')

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
otu.jaccard$site.loc <- rep(c('Central','South','North'), each = 3)
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$site.loc)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

# ANOSIM -----------------------------------------------------------------------------
# range 
jaccard.dist.anosim <- anosim(comm.matrix, grouping = otu.jaccard$site.loc,
                              distance = "jaccard")
jaccard.dist.anosim

results.em[results.em$measure == 'Jaccard.em', 'p-value'] <- jaccard.dist.anosim$signif
results.em[results.em$measure == 'Jaccard.em', 'R'] <- jaccard.dist.anosim$statistic

#=========================================================================================#
# Morisita-horn: Ectomycorrhizal fungi 97%----------
#=========================================================================================#

#--isolate otu data
comm.matrix <- em.site[10:length(em.site)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 2]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.horn<- em.site[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.horn <- vegdist(comm.matrix, method = "horn", binary = F)

#--NMDS analysis
horn.otu <- metaMDS(comm.dist.horn, dist = "bray",
                       try = 100, trymax = 100)

#--Stress
results.em[results.em$measure == 'Horn.em', 'stress'] <- horn.otu$stress

#--GGplot plot
# Create dataframe
data.scores <- as.data.frame(scores(horn.otu))
data.scores$site <- rep(c('central','south','north'),each = 3)

# outliers row N15, H21
#data.scores <- data.scores[!data.scores$tree %in% c('N15','H21'),]

#--GGplot: NMDS
em.horn <- ggplot() + 
  geom_point(data = data.scores,aes(x = NMDS1,
                                    y = NMDS2,
                                    color=site),size=2) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=28,margin = margin(t = 30)), # remove x-axis labels
        axis.title.y = element_text(size=28,margin = margin(r = 30)), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

em.horn

#--Base R plot:
#--Range----
#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'NDMS_MogollonRimMorisitaHorn_EM.jpeg'),
     width = 1200, height = 1000,
     quality = 100)

plot(horn.otu, display = "sites", type = "n", cex.lab = 1.5,
     cex.axis = 1.5)
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(otu.horn)),
                        p.group = otu.horn$Site)

color.vec[color.vec$p.group == 'M1', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M2', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M3', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M4', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M5', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M6', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M7', 'color'] <- 'darkgreen'
color.vec[color.vec$p.group == 'M8', 'color'] <- 'darkgreen'
color.vec[color.vec$p.group == 'M9', 'color'] <- 'darkgreen'
color.vec <- color.vec[order(color.vec$p.group),]

# ordipointlabel(jaccard.otu, display = "sites")

points(horn.otu, display = "sites", cex = 4,
       pch = 16,
       col = color.vec$color,
       bg = color.vec$color)

legend("topright", legend = c('North','Central','South'), bty = "n",
       col = c('darkgreen','grey','darkblue'),
       pch = 16, cex = 3)

#--Ordihull variations by range, burn, and both burn and range
ordihull(horn.otu, groups = color.vec$color,
         col = 'black')
#--Overlays
fit <- envfit(horn.otu ~ prec + Tavg, otu.horn)
fit
plot(fit, type = 'n')

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
otu.horn$site.loc <- rep(c('Central','South','North'), each = 3)
betadisp <- betadisper(comm.dist.horn, group = otu.horn$site.loc)
horn.betadisper <- anova(betadisp)
horn.betadisper

# ANOSIM -----------------------------------------------------------------------------
# range 
horn.dist.anosim <- anosim(comm.matrix, grouping = otu.horn$site.loc,
                              distance = "horn")
horn.dist.anosim

results.em[results.em$measure == 'Horn.em', 'p-value'] <- horn.dist.anosim$signif
results.em[results.em$measure == 'Horn.em', 'R'] <- horn.dist.anosim$statistic


#=========================================================================================#
# Distance decay: Ectomycorrhizal fungi 97%----------
#=========================================================================================#
#--isolate Mogollon Rim data
mogollon.rim <- c('M1','M2','M3','M4','M5','M6','M7','M8','M9')
dist.decay.mr <- dist.decay.em[dist.decay.em$site1 %in% mogollon.rim &
                              dist.decay.em$site2 %in% mogollon.rim,]
#-Remove self comparisons
dist.decay.mr <-  dist.decay.mr[!dist.decay.mr$jaccard.dissimilarity == 0, ]
#-Add column for similarity
dist.decay.mr$jaccard.similarity <- 1-dist.decay.mr$jaccard.dissimilarity
dist.decay.mr$horn.similarity <- 1-dist.decay.mr$horn.dissimilarity

# Jaccard
mr.jaccard <- lm(jaccard.dissimilarity ~ spatial.site * clim, data = dist.decay.mr)
summary(mr.jaccard)
ggplot(dist.decay.mr, aes(x = clim,
                          y = jaccard.dissimilarity))+
  geom_point()

# Morisita-Horn
mr.horn <- lm(horn.dissimilarity ~ spatial.site * clim, data = dist.decay.mr)
summary(mr.horn)
ggplot(dist.decay.mr, aes(x = spatial.site,
                          y = horn.dissimilarity))+
  geom_point()


#=========================================================================================#
# PERMANOVA: Ectomycorrhizal----------
#=========================================================================================#
# Jaccard
adonis(comm.dist.jaccard ~ dist.closest.forest * prec * forest.type, data = otu.jaccard)

# Morisita horn
adonis(comm.dist.horn ~ dist.closest.forest * prec * forest.type,
       data = otu.horn)

#=========================================================================================#
# Jaccard: Foliar endophytes 95%----------
#=========================================================================================#

#--isolate otu data
comm.matrix <- fe.tree[11:length(fe.tree)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- fe.tree[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)

#--Stress
results.em[results.em$measure == 'Jaccard.fe', 'stress'] <- jaccard.otu$stress

#--GGplot plot
# Create dataframe
data.scores <- as.data.frame(scores(jaccard.otu))
data.scores$tree <- otu.jaccard$tree

# outliers row N15, H21
#data.scores <- data.scores[!data.scores$tree %in% c('N15','H21'),]

#--GGplot: NMDS
em.jaccard <- ggplot() + 
  geom_point(data = data.scores,aes(x = NMDS1,
                                    y = NMDS2,
                                    color = tree),size=2) + # add the point markers
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_blank(),  # remove x-axis text
        axis.text.y = element_blank(), # remove y-axis text
        axis.ticks = element_blank(),  # remove axis ticks
        axis.title.x = element_text(size=28,margin = margin(t = 30)), # remove x-axis labels
        axis.title.y = element_text(size=28,margin = margin(r = 30)), # remove y-axis labels
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  #remove major-grid labels
        panel.grid.minor = element_blank(),  #remove minor-grid labels
        plot.background = element_blank())

em.jaccard

#--Base R plot:
#--Range----
#--format and output NMDS plots to figure folder
jpeg(filename = paste0(fig.dir, 'NDMS_MogollonRimJaccard_FE.jpeg'),
     width = 1200, height = 1000,
     quality = 100)

plot(jaccard.otu, display = "sites", cex.lab = 1.5,
     cex.axis = 1.5)
# colors for points
color.vec <- data.frame(color = rep(NA, nrow(otu.jaccard)),
                        p.group = otu.jaccard$tree)

color.vec[color.vec$p.group == 'M11', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M12', 'color'] <- 'darkblue'
color.vec[color.vec$p.group == 'M21', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M22', 'color'] <- 'grey'
color.vec[color.vec$p.group == 'M31', 'color'] <- 'darkgreen'

color.vec <- color.vec[order(color.vec$p.group),]

# ordipointlabel(jaccard.otu, display = "sites")

points(jaccard.otu, display = "sites", cex = 4,
       pch = 16,
       col = color.vec$color,
       bg = color.vec$color)

legend("topright", legend = c('Site1','Site2','Site3'), bty = "n",
       col = c('darkblue','grey','darkgreen'),
       pch = 16, cex = 3)

#--Ordihull variations by range, burn, and both burn and range
ordihull(jaccard.otu, groups = color.vec$color,
         col = 'black')
#--Overlays
fit <- envfit(jaccard.otu ~ prec + Tavg, otu.jaccard)
fit
#plot(fit, type = 'n')

dev.off()

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$Site)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

# PERMANOVA -----------------------------------------------------------------------------
# range 
jaccard.adonis <- adonis(comm.dist.jaccard ~ Site, data = otu.jaccard,
                         permutations = 2000)
jaccard.adonis

results.em[results.em$measure == 'Jaccard.fe', 'p-value'] <- jaccard.adonis$aov.tab$`Pr(>F)`[1]
results.em[results.em$measure == 'Jaccard.fe', 'R2'] <- jaccard.adonis$aov.tab$R2[1]
results.em[results.em$measure == 'Jaccard.fe', 'F'] <- jaccard.adonis$aov.tab$F.Model[1]
results.em[results.em$measure == 'Jaccard.fe', 'DF1'] <- jaccard.adonis$aov.tab$Df[1]
results.em[results.em$measure == 'Jaccard.fe', 'DF2'] <- jaccard.adonis$aov.tab$Df[3]

write.csv(results.em, paste0(res.dir, 'MogollonRimOrdination.csv'),row.names = F)
