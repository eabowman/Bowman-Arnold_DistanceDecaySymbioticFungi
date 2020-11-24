## Script created by Liz Bowman Sept.26, 2019
## for analysing abundance data in Sky Island study

#========================================================================================#
# Load libraries and data for abundance------------------
#========================================================================================#

#install.packages(c('dplyr','tidyr','vegan','ggplot2'))
library(dplyr);library(tidyr);library(vegan);library(ggplot2);library(nlme)

#--file paths
dat.dir <- 'data/'
dat.out <- 'data_output/'
fig.dir <- 'figures/'
res.dir <- 'results/'

#--Site x species matrix with tip abundance per OTU
otu.97.site <- read.csv('../3:4_Combined/data/EM_sitexspecies_TipAb_Site.csv', as.is = T)
otu.97.site[otu.97.site$Site == 'Butterfly trail','Site'] <- 'Butterfly.trail'
otu.97.site[otu.97.site$Site == 'Sunset Trail','Site'] <- 'Sunset.Trail'
otu.97.site[otu.97.site$Site == 'Upper Solders Camp Rd.','Site'] <- 'Upper.Solders.Camp.Rd.'

#--Site x species matrix Culture-free
otu.95.site.cf <- read.csv('../3:4_Combined/data/FE_ITS2rarified_sitexspecies.csv', as.is = T)
otu.95.site.cf[otu.95.site.cf$Site == 'Sunset Trail','Site'] <- 'Sunset.Trail'
otu.95.site.cf[otu.95.site.cf$Site == 'Upper Solders Camp Rd.','Site'] <- 'Upper.Solders.Camp.Rd.'

#--Site x species matrix Culture-based
otu.95.site.cb <- read.csv('../3:4_Combined/data/FE_CB_sitexspecies.csv', as.is = T)
otu.95.site.cb[otu.95.site.cb$site == 'Sunset Trail','site'] <- 'Sunset.Trail'
otu.95.site.cb[otu.95.site.cb$site == 'Upper Solders Camp Rd.','site'] <- 'Upper.Solders.Camp.Rd.'
otu.95.site.cb[otu.95.site.cb$site == 'Butterfly trail','site'] <- 'Butterfly.trail'

em.pairwise <- read.csv('data/EM_PairwiseDissimilarityMatrix_site.csv', as.is = T)
em.pairwise[em.pairwise$site1 == 'Sunset Trail','site1'] <- 'Sunset.Trail'
em.pairwise[em.pairwise$site2 == 'Sunset Trail','site2'] <- 'Sunset.Trail'
em.pairwise[em.pairwise$site1 == 'Upper Solders Camp Rd.','site1'] <- 'Upper.Solders.Camp.Rd.'
em.pairwise[em.pairwise$site2 == 'Upper Solders Camp Rd.','site2'] <- 'Upper.Solders.Camp.Rd.'
em.pairwise[em.pairwise$site1 == 'Butterfly trail','site1'] <- 'Butterfly.trail'
em.pairwise[em.pairwise$site2 == 'Butterfly trail','site2'] <- 'Butterfly.trail'

fe.cb.pairwise <- read.csv('data/FE_CB_PairwiseDissimilarityMatrix_site.csv',
                           as.is = T)
fe.cb.pairwise[fe.cb.pairwise$site1 == 'Sunset Trail','site1'] <- 'Sunset.Trail'
fe.cb.pairwise[fe.cb.pairwise$site2 == 'Sunset Trail','site2'] <- 'Sunset.Trail'
fe.cb.pairwise[fe.cb.pairwise$site1 == 'Upper Solders Camp Rd.','site1'] <- 'Upper.Solders.Camp.Rd.'
fe.cb.pairwise[fe.cb.pairwise$site2 == 'Upper Solders Camp Rd.','site2'] <- 'Upper.Solders.Camp.Rd.'
fe.cb.pairwise[fe.cb.pairwise$site1 == 'Butterfly trail','site1'] <- 'Butterfly.trail'
fe.cb.pairwise[fe.cb.pairwise$site2 == 'Butterfly trail','site2'] <- 'Butterfly.trail'

fe.cf.pairwise <- read.csv('data/FE_CF_PairwiseDissimilarityMatrix_site_ITS2rarefied.csv',
                           as.is = T)
fe.cf.pairwise[fe.cf.pairwise$site1 == 'Sunset Trail','site1'] <- 'Sunset.Trail'
fe.cf.pairwise[fe.cf.pairwise$site2 == 'Sunset Trail','site2'] <- 'Sunset.Trail'
fe.cf.pairwise[fe.cf.pairwise$site1 == 'Upper Solders Camp Rd.','site1'] <- 'Upper.Solders.Camp.Rd.'
fe.cf.pairwise[fe.cf.pairwise$site2 == 'Upper Solders Camp Rd.','site2'] <- 'Upper.Solders.Camp.Rd.'
fe.cf.pairwise[fe.cf.pairwise$site1 == 'Butterfly trail','site1'] <- 'Butterfly.trail'
fe.cf.pairwise[fe.cf.pairwise$site2 == 'Butterfly trail','site2'] <- 'Butterfly.trail'

#========================================================================================#
# Generate dissimilarity matrices for environmental data----
#========================================================================================#
em.env <- read.csv('data/EM_EnvFactors_site.csv', as.is = T)
fe.cb.env <- read.csv('data/FE_CB_EnvFactors_site.csv', as.is = T)
fe.cf.env <- read.csv('data/FE_CF_EnvFactors_site.csv', as.is = T)

# << EM >> ----
em.clim.pca <- dist(em.env$clim.pca, )

#========================================================================================#
# EM Community dissimilarity as a function of distance to contiguous forests: Grouped by site -------
#========================================================================================#
#-----------------------------------------------------------------#
#<< Create similarity index >> ---------------------------
#-----------------------------------------------------------------#

# With all OTUs --------
#--isolate otu data
comm.matrix <- otu.97.site[17:length(otu.97.site)]
rownames(comm.matrix) <- otu.97.site$Site

comm.matrix <- comm.matrix[colSums(comm.matrix) > 0]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F)
comm.dist.matrix <- as.data.frame(as.matrix(comm.dist.jaccard))
comm.dist.horn <- vegdist(comm.matrix, method = 'horn', binary = F)
comm.dist.horn.matrix <- as.data.frame(as.matrix(comm.dist.horn))

#--Make data frame
comm.dist.matrix$X <- rownames(comm.dist.matrix)
sim.jaccard <- gather(comm.dist.matrix, 'site2','jaccard.dissimilarity', -X)
names(sim.jaccard)[1] <- 'site1'

for(i in unique(sim.jaccard$site1)){
  for(t in unique(sim.jaccard$site2)){
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'env'] <-
      env.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'clim'] <-
      clim.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'spatial.site'] <-
      spatial.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'horn.dissimilarity'] <-
      comm.dist.horn.matrix[i, t]
    sim.jaccard[sim.jaccard$site1 == i, 'range'] <- otu.97.site[otu.97.site$Site == i, 'Range']
  }
}

# Remove samples with exactly the same similarity (0)
sim.jaccard.out <- sim.jaccard[sim.jaccard$dissimilarity > 0,]
sim.jaccard.out <- sim.jaccard.out[sim.jaccard.out$dissimilarity < 1,]
# Check normality
source('scripts/functions.R')
normtest(sim.jaccard.out, 'dissimilarity')
normtest(sim.jaccard.out, 'spatial')
normtest(sim.jaccard.out, 'clim')

#-log transform
sim.jaccard.out$log.clim <- log(sim.jaccard.out$clim)
sim.jaccard.out$log.spatial <- log(sim.jaccard.out$spatial)
sim.jaccard.out$log.env <- log(sim.jaccard.out$env)
sim.jaccard.out$log.dissimilarity <- log(sim.jaccard.out$dissimilarity)

normtest(sim.jaccard.out, 'log.spatial')
normtest(sim.jaccard.out, 'log.clim')
normtest(sim.jaccard.out, 'log.env')
normtest(sim.jaccard.out, 'log.dissimilarity')

# inverse of dissimilarity
sim.jaccard.out$dissimilarity.inv <- 1-sim.jaccard.out$dissimilarity
# outliers
sim.jaccard.out <- sim.jaccard.out[!sim.jaccard.out$dissimilarity == 0.4,]

#--------------------------------------------------------------------------------------#
# Bootstrapping multiple regression------------
#--------------------------------------------------------------------------------------#
#install.packages('boot')
library(boot)

# bootstrap function modified from Chpt. 7 in Discovering Statistics Using R
# formula is a regression formula
bootReg <- function(formula, data, indices){
  d <- data[indices,] # subset of dataframe, i refers to a particular bootstrap sample
  fit <- lm(formula, data = d)
  return(coef(fit)) # coef extractions the coefficients from a regression object
  # this will return the intercept and any slope coefficients for predictors in the model
  #return(summary(fit)$r.square)
}

# use boot() function to obtain bootstrap samples
bootResults <- boot(data = sim.jaccard.out,
                    statistic = bootReg, # function created above
                    formula = dissimilarity.inv ~ spatial * clim, 
                    R = 2000) # replications

# to look at CI using the return(coef(fit)) portion of bootReg function
boot.ci(bootResults, type = 'bca', index = 1) # bootstrap CI for intercept, index = 1
boot.ci(bootResults, type = 'bca', index = 2) # coefficient for spatial
boot.ci(bootResults, type = 'bca', index = 3) # coefficient for clim
boot.ci(bootResults, type = 'bca', index = 4) # coefficient for their interaction

# for use with the return(summary(fit)$r.squared) part of bootReg function
bootResults
plot(bootResults)
# get 95% confidence interval
boot.ci(bootResults, type="bca")

#<< Plot: Jaccard >> --------------------------
ggplot(sim.jaccard.out,
       aes(x = spatial,
           y = dissimilarity.inv)) +
  geom_point(size = 4) +
  geom_smooth(method="lm", color = 'darkgrey', se = F)  +
  ylab('Jaccard similarity') +
  xlab('Pairwise distance (km)') +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24),
        title = element_text(size = 22))

ggsave('EM_DistanceDecay_site.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7) #Variation partitioning
response <- otu.97.site[c('dist.cont.central', 'prec','forest.type')]
varpart(comm.dist.jaccard,
        ~ response$dist.cont.central,
        ~ response$prec,
        ~ response$forest.type)

#========================================================================================#
# FE culture-free cmmunity dissimilarity as a function of distance to contiguous forests: Grouped by site -------
#========================================================================================#
#-----------------------------------------------------------------#
#<< Create similarity index >> ---------------------------
#-----------------------------------------------------------------#

# << Greater than or equal to 8 occurences >> ------------
#--isolate otu data
comm.matrix <- otu.95.site.cf[15:length(otu.95.site.cf)]
rownames(comm.matrix) <- otu.95.site.cf$Site

comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F)
comm.dist.matrix <- as.data.frame(as.matrix(comm.dist.jaccard))

#--Make data frame
comm.dist.matrix$X <- rownames(comm.dist.matrix)
sim.jaccard <- gather(comm.dist.matrix, 'site2','jaccard.dissimilarity', -X)
names(sim.jaccard)[1] <- 'site1'
# isolate mogollon rim comparisons
#sim.jaccard <- sim.jaccard[sim.jaccard$site2 %in% c('M1','M2','M3','M4','M5','M6','M7','M8','M9'),]

for(i in unique(sim.jaccard$site1)){
  for(t in unique(sim.jaccard$site2)){
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'env'] <-
      env.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'clim'] <-
      clim.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'spatial.site'] <-
      spatial.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i, 'range'] <- otu.95.site[otu.95.site$Site == i, 'range']
  }
}

# Remove samples with exactly the same similarity (0)
sim.jaccard.out <- sim.jaccard[sim.jaccard$dissimilarity > 0,]
sim.jaccard.out <- sim.jaccard.out[sim.jaccard.out$dissimilarity < 1,]
# Check normality
source('scripts/functions.R')
normtest(sim.jaccard.out, 'dissimilarity')
normtest(sim.jaccard.out, 'spatial')
normtest(sim.jaccard.out, 'clim')

#-log transform
sim.jaccard.out$log.clim <- log(sim.jaccard.out$clim)
sim.jaccard.out$log.spatial <- log(sim.jaccard.out$spatial)
sim.jaccard.out$log.env <- log(sim.jaccard.out$env)
sim.jaccard.out$log.dissimilarity <- log(sim.jaccard.out$dissimilarity)

normtest(sim.jaccard.out, 'log.spatial')
normtest(sim.jaccard.out, 'log.clim')
normtest(sim.jaccard.out, 'log.env')
normtest(sim.jaccard.out, 'log.dissimilarity')

# inverse of dissimilarity
sim.jaccard.out$dissimilarity.inv <- 1-sim.jaccard.out$dissimilarity
# outliers
sim.jaccard.out <- sim.jaccard.out[!sim.jaccard.out$dissimilarity == 0.4,]

#--------------------------------------------------------------------------------------#
# Bootstrapping multiple regression------------
#--------------------------------------------------------------------------------------#
#install.packages('boot')
library(boot)

# bootstrap function modified from Chpt. 7 in Discovering Statistics Using R
# formula is a regression formula
bootReg <- function(formula, data, indices){
  d <- data[indices,] # subset of dataframe, i refers to a particular bootstrap sample
  fit <- lm(formula, data = d)
  return(coef(fit)) # coef extractions the coefficients from a regression object
  # this will return the intercept and any slope coefficients for predictors in the model
  #return(summary(fit)$r.square)
}

# use boot() function to obtain bootstrap samples
bootResults <- boot(data = sim.jaccard.out,
                    statistic = bootReg, # function created above
                    formula = dissimilarity.inv ~ spatial * clim, 
                    R = 2000) # replications

# to look at CI using the return(coef(fit)) portion of bootReg function
boot.ci(bootResults, type = 'bca', index = 1) # bootstrap CI for intercept, index = 1
boot.ci(bootResults, type = 'bca', index = 2) # coefficient for spatial
boot.ci(bootResults, type = 'bca', index = 3) # coefficient for clim
boot.ci(bootResults, type = 'bca', index = 4) # coefficient for their interaction

# for use with the return(summary(fit)$r.squared) part of bootReg function
bootResults
plot(bootResults)
# get 95% confidence interval
boot.ci(bootResults, type="bca")

#<< Plot: Jaccard >> --------------------------
ggplot(sim.jaccard.out,
       aes(x = spatial,
           y = dissimilarity.inv)) +
  geom_point(size = 4) +
  geom_smooth(method="lm", color = 'darkgrey', se = F)  +
  ylab('Jaccard similarity') +
  xlab('Pairwise distance (km)') +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24),
        title = element_text(size = 22))

ggsave('EM_DistanceDecay_site.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)#Variation partitioning
response <- otu.97.site[c('dist.cont.central', 'prec','forest.type')]
varpart(comm.dist.jaccard,
        ~ response$dist.cont.central,
        ~ response$prec,
        ~ response$forest.type)

#========================================================================================#
# FE culture based Community dissimilarity as a function of distance to contiguous forests: Grouped by site -------
#========================================================================================#
#-----------------------------------------------------------------#
#<< Create similarity index >> ---------------------------
#-----------------------------------------------------------------#

# << all data >>----
#--isolate otu data
comm.matrix <- otu.95.site.cb[12:length(otu.95.site.cb)]
rownames(comm.matrix) <- otu.95.site.cb$site

#--comment to include singletons
comm.matrix <- comm.matrix[colSums(comm.matrix) > 0]
comm.matrix <- comm.matrix[rowSums(comm.matrix) > 0, ] # remove rows with sums of 0
otu.jaccard <- otu.95.site.cb[rownames(otu.95.site.cb) %in% rownames(comm.matrix),]
otu.horn <- otu.95.site.cb[rownames(otu.95.site.cb) %in% rownames(comm.matrix),]

#--distance matrix using jaccard index
comm.dist.jaccard <- vegdist(comm.matrix, method = "jaccard", binary = F)
comm.dist.matrix <- as.data.frame(as.matrix(comm.dist.jaccard))
#--distance matrix using Morisita horn index
comm.dist.horn <- vegdist(comm.matrix, method = 'horn')
comm.dist.horn <- as.data.frame(as.matrix(comm.dist.horn))

#--Make data frame
comm.dist.matrix$X <- rownames(comm.dist.matrix)
sim.jaccard <- gather(comm.dist.matrix, 'site2','jaccard.dissimilarity', -X)
names(sim.jaccard)[1] <- 'site1'
# isolate mogollon rim comparisons
#sim.jaccard <- sim.jaccard[sim.jaccard$site2 %in% c('M1','M2','M3','M4','M5','M6','M7','M8','M9'),]

for(i in unique(sim.jaccard$site1)){
  for(t in unique(sim.jaccard$site2)){
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'spatial.site'] <-
      spatial.dist.site[i, t]
    sim.jaccard[sim.jaccard$site1 == i, 'range.1'] <- otu.95.site[otu.95.site$site == i, 'range']
    sim.jaccard[sim.jaccard$site2 == t, 'range.2'] <- otu.95.site[otu.95.site$site == t, 'range']
    sim.jaccard[sim.jaccard$site1 == i & sim.jaccard$site2 == t, 'horn.dissimilarity'] <-
      comm.dist.horn[i, t]
  }
}


# range comparison
for(i in 1:nrow(sim.jaccard)){
  if(sim.jaccard[i,'range.1'] == sim.jaccard[i,'range.2']){sim.jaccard[i,'range.comp'] <- 'Same'
  } else{sim.jaccard[i,'range.comp'] <- 'Different'}
}

# Remove samples with exactly the same similarity (0)
sim.jaccard.out <- sim.jaccard[sim.jaccard$dissimilarity > 0,]
sim.jaccard.out <- sim.jaccard.out[sim.jaccard.out$dissimilarity < 1,]
# Check normality
source('scripts/functions.R')
normtest(sim.jaccard.out, 'dissimilarity')
normtest(sim.jaccard.out, 'spatial')
normtest(sim.jaccard.out, 'clim')

#-log transform
sim.jaccard.out$log.clim <- log(sim.jaccard.out$clim)
sim.jaccard.out$log.spatial <- log(sim.jaccard.out$spatial)
sim.jaccard.out$log.env <- log(sim.jaccard.out$env)
sim.jaccard.out$log.dissimilarity <- log(sim.jaccard.out$dissimilarity)

normtest(sim.jaccard.out, 'log.spatial')
normtest(sim.jaccard.out, 'log.clim')
normtest(sim.jaccard.out, 'log.env')
normtest(sim.jaccard.out, 'log.dissimilarity')

# inverse of dissimilarity
sim.jaccard.out$dissimilarity.inv <- 1-sim.jaccard.out$dissimilarity
# outliers
sim.jaccard.out <- sim.jaccard.out[!sim.jaccard.out$dissimilarity == 0.4,]

#--------------------------------------------------------------------------------------#
# Bootstrapping multiple regression------------
#--------------------------------------------------------------------------------------#
#install.packages('boot')
library(boot)

# bootstrap function modified from Chpt. 7 in Discovering Statistics Using R
# formula is a regression formula
bootReg <- function(formula, data, indices){
  d <- data[indices,] # subset of dataframe, i refers to a particular bootstrap sample
  fit <- lm(formula, data = d)
  return(coef(fit)) # coef extractions the coefficients from a regression object
  # this will return the intercept and any slope coefficients for predictors in the model
  #return(summary(fit)$r.square)
}

# use boot() function to obtain bootstrap samples
bootResults <- boot(data = sim.jaccard.out,
                    statistic = bootReg, # function created above
                    formula = dissimilarity.inv ~ spatial * clim, 
                    R = 2000) # replications

# to look at CI using the return(coef(fit)) portion of bootReg function
boot.ci(bootResults, type = 'bca', index = 1) # bootstrap CI for intercept, index = 1
boot.ci(bootResults, type = 'bca', index = 2) # coefficient for spatial
boot.ci(bootResults, type = 'bca', index = 3) # coefficient for clim
boot.ci(bootResults, type = 'bca', index = 4) # coefficient for their interaction

# for use with the return(summary(fit)$r.squared) part of bootReg function
bootResults
plot(bootResults)
# get 95% confidence interval
boot.ci(bootResults, type="bca")

#<< Plot: Jaccard >> --------------------------
ggplot(sim.jaccard.out,
       aes(x = spatial,
           y = dissimilarity.inv)) +
  geom_point(size = 4) +
  geom_smooth(method="lm", color = 'darkgrey', se = F)  +
  ylab('Jaccard similarity') +
  xlab('Pairwise distance (km)') +
  theme_classic() +
  theme(axis.text=element_text(size=20),
        axis.title=element_text(size=24),
        title = element_text(size = 22))

ggsave('EM_DistanceDecay_site.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)#Variation partitioning
response <- otu.97.site[c('dist.cont.central', 'prec','forest.type')]
varpart(comm.dist.jaccard,
        ~ response$dist.cont.central,
        ~ response$prec,
        ~ response$forest.type)

