## Script created for analyzing community data 
## Includes ectomycorrhizal data and endophyte data
## Presents an alternative method to the partial mantel test
## Script written by Dr. Elizabeth Bowman November 22, 2020
## All code and methodology was based on Chapter 7: Spatial analysis of 
## ecological data in Numerical Ecology with R (2nd edition)
## by D Borcard, F Gillet, P Legendre 2018.

## The output from these analyses can be seen in Supplementary Table S6.
 
#========================================================================================#
# Load data------------------
#========================================================================================#

#--Site x species matrix with tip abundance per OTU
em.site <- read.csv(paste0(dat.dir,'EM_sitexspecies_TipAb_Site.csv'), as.is = T)

#--Site x species matrix ITS2 rarefied
fe.cf.site <- read.csv(paste0(dat.dir,'FE_ITS2rarified_sitexspecies.csv'), as.is = T)

#--Site x species matrix culture-based FE
fe.cb.site <- read.csv(paste0(dat.dir, 'FE_CB_sitexspecies.csv'), as.is = T)

#--Environmental data
em.env <- read.csv(paste0(dat.dir, 'EM_EnvFactors_site.csv'), as.is = T)
fe.cf.env <- read.csv(paste0(dat.dir, 'FE_CF_EnvFactors_site.csv'), as.is = T)
fe.cf.env <- fe.cf.env[!fe.cf.env$site == 'P3',]
fe.cb.env <- read.csv(paste0(dat.dir, 'FE_CB_EnvFactors_site.csv'), as.is = T)

# << Read in functions >> ----
# varpart2.MEM.R was written by Legendre et al. 2012 'Variation partitioning
# involving orthogonal spatial eigenfunction submodels'

source('scripts/varpart2.MEM.R')

#========================================================================================#
# distance based Moran's Eigenvector Map: Variation partitioning----
#========================================================================================#

# << Ectomycorrhizal fungi: overall >> -----
# add longitudinal data to em.site
for(i in em.site$Site){
  em.site[em.site$Site == i, 'long'] <- em.env[em.env$Site == i, 'long']
}
em.site <- em.site[c(1:16,298, 17:297)] # reorganize column order

# isolate community data
comm.matrix <- em.site[18:length(em.site)]

# Remove rare OTU (i.e. OTU with < 4 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]

# isolate location data
em.xy <- em.site[c('lat','long')]

## construct the matrix of dbMEM variables
em.dbmem.temp <- dbmem(em.xy, silent = F)
em.dbmem <- as.data.frame(em.dbmem.temp)
(thr <- give.thresh(dist(em.xy)))
# display and count the eigenvalues
attributes(em.dbmem.temp)$values

# Hellinger transformation of hte species data 
em.hel <- decostand(comm.matrix, 'hel')
em.horn <- vegdist(em.hel, method = 'horn')
em.jac <- vegdist(em.hel, method = 'jaccard')

## 1. Is there linear trend in EM data?
anova(rda(em.hel, em.xy)) # result: significant trend
# computate linearly detrended EM data
em.hel.det <- resid(lm(as.matrix(em.hel) ~ ., data = em.xy))

## 2. Test and forward selection of the environmental variables
# Recode forest data to dummy variables
forest.type <- model.matrix( ~ em.env[,11])[,-1]
em.env2 <- cbind(em.env[,c(6,8,9,10)], forest.type)
colnames(em.env2) <- c('BIO1','BIO12','BIO16','BIO17',
                       'forest_PineDoug','forest_PineOak')

# Forward selection of the environmental variables
em.env.rda <- rda(em.hel ~ ., data = em.env2)
(em.env.R2a <- RsquareAdj(em.env.rda)$adj.r.squared)
em.env.fwd <- forward.sel(em.hel, em.env2,
                          adjR2thresh = em.env.R2a,
                          nperm = 9999)
sort(em.env.fwd$order) # BIO16
env.red <- em.env2[,'BIO16']

## 3. test and foward selection of the dbMEM variables
# Run the global dbMEM analysis on the detrended mite data
em.det.dbmem.rda <- rda(em.hel.det ~ ., em.dbmem)
anova(em.det.dbmem.rda) # not significant

## 4. variation partitioning
# Jaccard
em.jac.MEM.varpart2 <- varpart2.MEM(em.jac, em.dbmem, em.env2$BIO16,
                                is.MEM = c(1,2), method = 'h')
em.jac.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(em.jac.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# Morisita-horn
em.horn.MEM.varpart2 <- varpart2.MEM(em.horn, em.dbmem, em.env2$BIO16,
                                is.MEM = c(1,2), method = 'h')
em.horn.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(em.horn.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# << Ectomycorrhizal fungi: regional, non-Mogollon Rim sites >> -----
# Exclude Mogollon Rim data
em.site.nonmr <- filter(em.site, !Range == 'MogollonRim')
em.env.nonmr <- filter(em.env, !range == 'mogollon')

# isolate community data
comm.matrix <- em.site.nonmr[18:length(em.site.nonmr)]

# Remove rare OTU (i.e. OTU with < 4 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]

# isolate location data
em.xy <- em.site.nonmr[c('lat','long')]

## construct the matrix of dbMEM variables
em.dbmem.temp <- dbmem(em.xy, silent = F)
em.dbmem <- as.data.frame(em.dbmem.temp)
(thr <- give.thresh(dist(em.xy)))
# display and count the eigenvalues
attributes(em.dbmem.temp)$values

# Hellinger transformation of hte species data 
em.hel <- decostand(comm.matrix, 'hel')
em.horn <- vegdist(em.hel, method = 'horn')
em.jac <- vegdist(em.hel, method = 'jaccard')

## 1. Is there linear trend in EM data?
anova(rda(em.hel, em.xy)) # result: significant trend
# computate linearly detrended EM data
em.hel.det <- resid(lm(as.matrix(em.hel) ~ ., data = em.xy))

## 2. Test and forward selection of the environmental variables
# Recode forest data to dummy variables
forest.type <- model.matrix( ~ em.env.nonmr[,11])[,-1]
em.env2 <- cbind(em.env.nonmr[,c(6,8,9,10)], forest.type)
colnames(em.env2) <- c('BIO1','BIO12','BIO16','BIO17',
                       'forest_PineDoug','forest_PineOak')

# Forward selection of the environmental variables
em.env.rda <- rda(em.hel ~ ., data = em.env2)
(em.env.R2a <- RsquareAdj(em.env.rda)$adj.r.squared)
em.env.fwd <- forward.sel(em.hel, em.env2,
                          adjR2thresh = em.env.R2a,
                          nperm = 9999)
sort(em.env.fwd$order) # BIO16
env.red <- em.env2[,'BIO16']

## 3. test and foward selection of the dbMEM variables
# Run the global dbMEM analysis on the detrended mite data
em.det.dbmem.rda <- rda(em.hel.det ~ ., em.dbmem)
anova(em.det.dbmem.rda) # not significant

## 4. variation partitioning
# Jaccard
em.jac.MEM.varpart2 <- varpart2.MEM(em.jac, em.dbmem, em.env2$BIO16,
                                is.MEM = c(1,2), method = 'h')
em.jac.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(em.jac.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# Morisita-horn
em.horn.MEM.varpart2 <- varpart2.MEM(em.horn, em.dbmem, em.env2$BIO16,
                                is.MEM = c(1,2), method = 'h')
em.horn.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(em.horn.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# << Ectomycorrhizal fungi: regional, Mogollon Rim sites included>> -----
# Include only Mogollon Rim data
em.site.mr <- filter(em.site, Range == 'MogollonRim')
em.env.mr <- filter(em.env, range == 'mogollon')

# isolate community data
comm.matrix <- em.site.mr[18:length(em.site.mr)]

# Remove rare OTU (i.e. OTU with < 4 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]

# isolate location data
em.xy <- em.site.mr[c('lat','long')]

## construct the matrix of dbMEM variables
em.dbmem.temp <- dbmem(em.xy, silent = F)
em.dbmem <- as.data.frame(em.dbmem.temp)
(thr <- give.thresh(dist(em.xy)))
# display and count the eigenvalues
attributes(em.dbmem.temp)$values

# Hellinger transformation of hte species data 
em.hel <- decostand(comm.matrix, 'hel')
em.horn <- vegdist(em.hel, method = 'horn')
em.jac <- vegdist(em.hel, method = 'jaccard')

## 1. Is there linear trend in EM data?
anova(rda(em.hel, em.xy)) # result: significant trend
# computate linearly detrended EM data
em.hel.det <- resid(lm(as.matrix(em.hel) ~ ., data = em.xy))

## 2. Test and forward selection of the environmental variables
# Recode forest data to dummy variables
forest.type <- model.matrix( ~ em.env.mr[,11])[,-1]
em.env2 <- cbind(em.env.mr[,c(6,8,9,10)], forest.type)
colnames(em.env2) <- c('BIO1','BIO12','BIO16','BIO17',
                       'forest_PineDoug','forest_PineOak')

# Forward selection of the environmental variables
em.env.rda <- rda(em.hel ~ ., data = em.env2)
(em.env.R2a <- RsquareAdj(em.env.rda)$adj.r.squared)
em.env.fwd <- forward.sel(em.hel, em.env2,
                          adjR2thresh = em.env.R2a,
                          nperm = 9999)
sort(em.env.fwd$order) # BIO16
env.red <- em.env2[,'BIO16']

## 3. test and foward selection of the dbMEM variables
# Run the global dbMEM analysis on the detrended mite data
em.det.dbmem.rda <- rda(em.hel.det ~ ., em.dbmem)
anova(em.det.dbmem.rda) # not significant

## 4. variation partitioning
# Jaccard
em.jac.MEM.varpart2 <- varpart2.MEM(em.jac, em.dbmem, em.env2$BIO16,
                                is.MEM = c(1,2), method = 'h')
em.jac.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(em.jac.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# Morisita-horn
em.horn.MEM.varpart2 <- varpart2.MEM(em.horn, em.dbmem, em.env2$BIO16,
                                is.MEM = c(1,2), method = 'h')
em.horn.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(em.horn.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# << Endophytic fungi: culture-free data, overall >> -----
# add longitudinal data to fe.cf.site
for(i in fe.cf.site$Site){
  fe.cf.site[fe.cf.site$Site == i, 'long'] <- fe.cf.env[fe.cf.env$site == i, 'long']
}
fe.cf.site <- fe.cf.site[c(1:3,1197,4:1196)] # reorganize column order

# isolate community data
comm.matrix <- fe.cf.site[16:length(fe.cf.site)]

# Remove rare OTU (i.e. OTU with < 8 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]

# isolate location data
fe.cf.xy <- fe.cf.site[c('lat','long')]

## construct the matrix of dbMEM variables
fe.cf.dbmem.temp <- dbmem(fe.cf.xy, silent = F)
fe.cf.dbmem <- as.data.frame(fe.cf.dbmem.temp)
(thr <- give.thresh(dist(fe.cf.xy)))
# display and count the eigenvalues
attributes(fe.cf.dbmem.temp)$values

# Hellinger transformation of hte species data 
fe.cf.hel <- decostand(comm.matrix, 'hel')
fe.cf.horn <- vegdist(fe.cf.hel, method = 'horn')
fe.cf.jac <- vegdist(fe.cf.hel, method = 'jaccard')

## 1. Is there linear trend in EM data?
anova(rda(fe.cf.hel, fe.cf.xy)) # result: significant trend
# computate linearly detrended EM data
fe.cf.hel.det <- resid(lm(as.matrix(fe.cf.hel) ~ ., data = fe.cf.xy))

## 2. Test and forward selection of the environmental variables
# Recode forest data to dummy variables
forest.type <- model.matrix( ~ fe.cf.env[,11])[,-1]
fe.cf.env2 <- cbind(fe.cf.env[,c(6,8,9,10)], forest.type)
colnames(em.env2) <- c('BIO1','BIO12','BIO16','BIO17',
                       'forest_PineDoug')

# Forward selection of the environmental variables
fe.cf.env.rda <- rda(fe.cf.hel ~ ., data = fe.cf.env2)
(fe.cf.env.R2a <- RsquareAdj(fe.cf.env.rda)$adj.r.squared)
fe.cf.env.fwd <- forward.sel(fe.cf.hel, fe.cf.env2,
                          adjR2thresh = fe.cf.env.R2a,
                          nperm = 9999)
temp <- sort(fe.cf.env.fwd$order) # 
fe.env2.for <- fe.cf.env2[temp]

## 3. test and foward selection of the dbMEM variables
# Run the global dbMEM analysis on the detrended mite data
fe.cf.det.dbmem.rda <- rda(fe.cf.hel.det ~ ., fe.cf.dbmem)
anova(fe.cf.det.dbmem.rda) # not significant

## 4. variation partitioning
# Jaccard
fe.cf.jac.MEM.varpart2 <- varpart2.MEM(fe.cf.jac, fe.cf.dbmem, fe.env2.for,
                                is.MEM = c(1,2), method = 'h')
fe.cf.jac.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(fe.cf.jac.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# << Endophytic fungi: culture-based data, overall >> -----
# add longitudinal data to fe.cb.site
for(i in fe.cb.site$site){
  fe.cb.site[fe.cb.site$site == i, 'long'] <- fe.cb.env[fe.cb.env$site == i, 'long']
}
fe.cb.site <- fe.cb.site[c(1:3,110,4:109)] # reorganize column order

# isolate community data
comm.matrix <- fe.cb.site[13:length(fe.cb.site)]

# Remove rare OTU (i.e. OTU with < 4 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]

# isolate location data
fe.cb.xy <- fe.cb.site[c('lat','long')]

## construct the matrix of dbMEM variables
fe.cb.dbmem.temp <- dbmem(fe.cb.xy, silent = F)
fe.cb.dbmem <- as.data.frame(fe.cb.dbmem.temp)
(thr <- give.thresh(dist(fe.cb.xy)))
# display and count the eigenvalues
attributes(fe.cb.dbmem.temp)$values

# Hellinger transformation of hte species data 
fe.cb.hel <- decostand(comm.matrix, 'hel')
fe.cb.horn <- vegdist(fe.cb.hel, method = 'horn')
fe.cb.jac <- vegdist(fe.cb.hel, method = 'jaccard')

## 1. Is there linear trend in EM data?
anova(rda(fe.cb.hel, fe.cb.xy)) # result: not significant

## 2. Test and forward selection of the environmental variables
# Recode forest data to dummy variables
forest.type <- model.matrix( ~ fe.cb.env[,11])[,-1]
fe.cb.env2 <- cbind(fe.cb.env[,c(6,8,9,10)], forest.type)
colnames(fe.cb.env2) <- c('BIO1','BIO12','BIO16','BIO17',
                       'forest_PineDoug')

# Forward selection of the environmental variables
fe.cb.env.rda <- rda(fe.cb.hel ~ ., data = fe.cb.env2)
(fe.cb.env.R2a <- RsquareAdj(fe.cb.env.rda)$adj.r.squared)
fe.cb.env.fwd <- forward.sel(fe.cb.hel, fe.cb.env2,
                          adjR2thresh = fe.cb.env.R2a,
                          nperm = 9999)
sort(fe.cb.env.fwd$order) # 


## 3. test and foward selection of the dbMEM variables
# Run the global dbMEM analysis on the detrended mite data
fe.cb.det.dbmem.rda <- rda(fe.cb.hel ~ ., fe.cb.dbmem)
anova(fe.cb.det.dbmem.rda) # not significant

## 4. variation partitioning
# Jaccard
fe.cb.jac.MEM.varpart2 <- varpart2.MEM(fe.cb.jac, fe.cb.dbmem, fe.cb.env2,
                                is.MEM = c(1,2), method = 'h')
fe.cb.jac.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(fe.cb.jac.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))

# Morisita-horn
fe.cb.horn.MEM.varpart2 <- varpart2.MEM(fe.cb.horn, fe.cb.dbmem, fe.cb.env2,
                                is.MEM = c(1,2), method = 'h')
fe.cb.horn.MEM.varpart2
showvarparts(2, bg=c(1,3))
plot(fe.cb.horn.MEM.varpart2, Xnames = c('spatial', 'env.'),
     bg=c(1,3))
