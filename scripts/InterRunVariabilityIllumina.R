## Script created by Liz Bowman January 21, 2020
## for analysing community data in Sky Island study
## Includes culture-free endophyte data

## For assess inter-run variability by assessing via ordination whether samples
## cluster together based on geographic orgin or illumina run

#========================================================================================#
# Load libraries and data------------------
#========================================================================================#

library(vegan)
fig.dir <- 'figures/'

#--otu data grouped by tree
otu.data <- read.csv('data/FE_ITS2rarified_treexspecies.csv', as.is = T)

#--otu data grouped by site
otu.site <- read.csv('data/FE_ITS2rarified_sitexspecies.csv', as.is = T)

#--remove Mingus, Pinaleno, Santa Catalina, and Huachuca
#otu.data <- otu.data[otu.data$range %in% c('MogollonRim','Chiricahua','Bradshaw'),]

#----------------------------------------------------------------#
# Jaccard dissimilarity: grouped by tree
#----------------------------------------------------------------#

#--Remove outliers M21
otu.data <- otu.data[!otu.data$sample %in% c('M21', 'B12n'),]

#--isolate otu data
comm.matrix <- otu.data[11:length(otu.data)]

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- otu.data[row.names(comm.matrix),] # keeps same rows in otu.99 metadata as comm.matrix

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = TRUE)

#--NMDS analysis
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)
jaccard.otu$stress

##--Distance from closest Ponderosa forest----
pdf(file = paste0(fig.dir, 'FE_InterRunVariability_jaccard.pdf'),
    width = 12.5, height = 8.5)

plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
#--color for groups
color.vec <- data.frame (color = rep(NA,length(rownames(comm.matrix))),
                         p.group = otu.jaccard$sample,
                         b.group = otu.jaccard$range)
#--color by range
# color.vec <- sapply(1:nrow(color.vec), function (x)
#   if(color.vec[x,'p.group'] == 'C23n') {color.vec[x,'color'] = 'black'}
#   else if(color.vec[x,'p.group'] == 'M21n') {color.vec[x,'color'] = 'darkblue'}
#   else if(color.vec[x,'b.group'] == 'Chiricahua') {color.vec[x,'color'] = 'grey'}
#   else if(color.vec[x,'b.group'] == 'MogollonRim') {color.vec[x,'color'] = 'blue'}
#   else if(color.vec[x,'b.group'] == 'Mingus') {color.vec[x,'color'] = 'yellow'}
#   else if(color.vec[x,'b.group'] == 'Santa.Catalina') {color.vec[x,'color'] = 'orange'}
#   else if(color.vec[x,'b.group'] == 'Pinaleno') {color.vec[x,'color'] = 'purple'}
#   else if(color.vec[x,'b.group'] == 'Huachuca') {color.vec[x,'color'] = 'brown'}
#   else {color.vec[x,'color'] = 'green'})
#--color by run1 or run2
color.vec <- sapply(1:nrow(color.vec), function (x)
  if(color.vec[x,'p.group'] == 'C23n') {color.vec[x,'color'] = 'black'}
  else if(color.vec[x,'p.group'] == 'M21n') {color.vec[x,'color'] = 'black'}
  else if(color.vec[x,'b.group'] == 'Mingus') {color.vec[x,'color'] = 'black'}
  else if(color.vec[x,'b.group'] == 'Santa.Catalina') {color.vec[x,'color'] = 'black'}
  else {color.vec[x,'color'] = 'grey'})

points(jaccard.otu, display = "sites", cex = 3,
       pch = 19,
       col = color.vec,
       bg = color.vec)
#ordilabel(jaccard.otu, display = 'sites', labels = otu.jaccard$sample)

dev.off()

#--add column for illumina run
otu.jaccard$illumina.run <- NA
for(i in 1:nrow(otu.jaccard)){
  if(otu.jaccard[i, 'range'] %in% c('Santa.Catalina','Mingus')){otu.jaccard[i, 'illumina.run'] <- 'run2'}
  else if(otu.jaccard[i, 'sample'] %in% c('M21n','C23n')){otu.jaccard[i, 'illumina.run'] <- 'run2'}
  else{otu.jaccard[i,'illumina.run'] <- 'run1'}
}

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$illumina.run)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

#--ANOSIM ----------------------------------------
jaccard.anosim <- anosim(comm.matrix, grouping = otu.jaccard$illumina.run, distance = 'jaccard')
jaccard.anosim

#--PERMANOVA
jaccard.adonis <- adonis(comm.dist.jaccard ~ illumina.run * range, data = otu.jaccard)
jaccard.adonis

#----------------------------------------------------------------#
# Jaccard dissimilarity: grouped by site
#----------------------------------------------------------------#

#--isolte community data
comm.matrix <- otu.site[15:length(otu.site)]

#--remove OTUs with less than 8 occurrences 
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]

#--distance matrix based on Jaccard dissimilarity
comm.dist.jaccard <- vegdist(comm.matrix, method = 'jaccard')

#--run nmds using distance matrix generated above
jaccard.otu <- metaMDS(comm.dist.jaccard, dist = "bray",
                       try = 100, trymax = 100)

##--Distance from closest Ponderosa forest----
pdf(file = paste0(fig.dir, 'FE_InterRunVariability_BySite_jaccard.pdf'),
    width = 12.5, height = 8.5)

plot(jaccard.otu, display = "sites", type = "n", cex.lab = 2.5,
     cex.axis = 2.5, xlab = 'Axis 1', ylab = 'Axis 2')
#--color for groups
color.vec <- data.frame (color = rep(NA,length(rownames(comm.matrix))),
                         b.group = otu.jaccard$range)
#--color by run1 or run2
color.vec$color <- sapply(1:nrow(color.vec), function (x)
  if(color.vec[x,'b.group'] == 'Mingus') {color.vec[x,'color'] = 'black'}
  else if(color.vec[x,'b.group'] == 'Santa.Catalina') {color.vec[x,'color'] = 'black'}
  else {color.vec[x,'color'] = 'grey'})

points(jaccard.otu, display = "sites", cex = 3,
       pch = 19,
       col = color.vec$color,
       bg = color.vec$color)
#ordilabel(jaccard.otu, display = 'sites', labels = otu.jaccard$sample)

dev.off()

#--add column for illumina run
otu.jaccard$illumina.run <- NA
for(i in 1:nrow(otu.jaccard)){
  if(otu.jaccard[i, 'range'] %in% c('Santa.Catalina','Mingus')){otu.jaccard[i, 'illumina.run'] <- 'run2'}
  else{otu.jaccard[i,'illumina.run'] <- 'run1'}
}

#--BetaDisper
#--a multivariate analogue of Levene's test for homogeneity of variance
betadisp <- betadisper(comm.dist.jaccard, group = otu.jaccard$illumina.run)
jaccard.betadisper <- anova(betadisp)
jaccard.betadisper

#--ANOSIM ----------------------------------------
jaccard.anosim <- anosim(comm.matrix, grouping = otu.jaccard$illumina.run, distance = 'jaccard')
jaccard.anosim

#--PERMANOVA
jaccard.adonis <- adonis(comm.dist.jaccard ~ illumina.run * range, data = otu.jaccard)
jaccard.adonis
