## Script created for analyzing community data 
## Includes ectomycorrhizal data and endophyte data
## Presents an alternative method to the partial mantel test
## Script created by Dr. Elizabeth Bowman November 22, 2020

## We followed methods written in "Moran's Eigenvector Maps and related methods
## for the spatial multiscale analysis of ecological data" by Stéphane Dray
## https://cran.r-project.org/web/packages/adespatial/vignettes/tutorial.html#multiscale-analysis-with-mem

#========================================================================================#
# Load data------------------
#========================================================================================#

# install.packages('adespatial', dep = T)
library(adespatial)
# install.packages('ade4', dep = T)
library(ade4)
# install.packages(c('adegraphics','spdep','maptools'))
library(adegraphics); library(spdep); library(maptools)
# install.packages(c('vegan','packfor','PCNM'))
library(vegan); library(PCNM)

#--file paths
dat.dir <- './data/'
fig.dir <- './figures/'
res.dir <- './results/'

#--Site x species matrix with tip abundance per OTU
em.site <- read.csv(paste0(dat.dir,'EM_sitexspecies_TipAb_Site.csv'), as.is = T)

#--Site x species matrix ITS2 rarefied
fe.cf.site <- read.csv(paste0(dat.dir,'FE_ITS2rarified_sitexspecies.csv'), as.is = T)

#--Site x species matrix culture-based FE
fe.cb.site <- read.csv(paste0(dat.dir, 'FE_CB_sitexspecies.csv'), as.is = T)

#--Environmental data
em.env <- read.csv(paste0(dat.dir, 'EM_EnvFactors_site.csv'), as.is = T)
fe.cf.env <- read.csv(paste0(dat.dir, 'FE_CF_EnvFactors_site.csv'), as.is = T)
fe.cb.env <- read.csv(paste0(dat.dir, 'FE_CB_EnvFactors_site.csv'), as.is = T)


#========================================================================================#
# distance based Moran's Eigenvector Map: Variation partitioning----
#========================================================================================#
# << Ectomycorrhizal fungi >> -----
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

# Compute MEM eigenfunctions
em.dbMEM <- PCNM(dist(em.xy), dbMEM = T)
em.MEM <- em.dbMEM$vectors

# Hellinger transformation of hte species data 
em.hel <- decostand(comm.matrix, 'hel')

# isolate environmental data
em.env <- em.site[c('BIO1','BIO12','BIO16','BIO17','forest.type')]

# RDA
temp <- rda(em.hel, em.MEM)
RsquareAdj(temp)

# Hierarchical partiitoning of the shared fractions
em.MEM.varpart2 <- varpart2.MEM(em.hel, em.MEM, em.env,
                                is.MEM = c(1,2), method = 'h')
em.MEM.varpart2
plot(em.MEM.varpart2, Xnames = c('spatial', 'environment'))

# << Endophytic fungi: culture-free data >> -----
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

# Compute MEM eigenfunctions
fe.cf.dbMEM <- PCNM(dist(fe.cf.xy), dbMEM = T)
fe.cf.MEM <- fe.cf.dbMEM$vectors

# Hellinger transformation of hte species data 
fe.cf.hel <- decostand(comm.matrix, 'hel')

# isolate environmental data
fe.cf.env <- fe.cf.site[c('BIO1','BIO12','BIO16','BIO17','forest')]

# RDA
temp <- rda(fe.cf.hel, fe.cf.MEM)
RsquareAdj(temp)

# Hierarchical partiitoning of the shared fractions
fe.cf.MEM.varpart2 <- varpart2.MEM(fe.cf.hel, fe.cf.MEM, fe.cf.env,
                                is.MEM = c(1,2), method = 'h')
fe.cf.MEM.varpart2
plot(fe.cf.MEM.varpart2, Xnames = c('spatial', 'environment'))

# << Endophytic fungi: culture-based data >> -----
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

# Compute MEM eigenfunctions
fe.cb.dbMEM <- PCNM(dist(fe.cb.xy), dbMEM = T)
fe.cb.MEM <- fe.cb.dbMEM$vectors

# Hellinger transformation of hte species data 
fe.cb.hel <- decostand(comm.matrix, 'hel')

# isolate environmental data
fe.cb.env <- fe.cb.site[c('BIO1','BIO12','BIO16','BIO17','forest')]

# RDA
temp <- rda(fe.cb.hel, fe.cb.MEM)
RsquareAdj(temp)

# Hierarchical partiitoning of the shared fractions
fe.cb.MEM.varpart2 <- varpart2.MEM(fe.cb.hel, fe.cb.MEM, fe.cb.env,
                                is.MEM = c(1,2), method = 'h')
fe.cb.MEM.varpart2
plot(fe.cb.MEM.varpart2, Xnames = c('spatial', 'environment'))



#========================================================================================#
# distance based Moran's Eigenvector Map----
#========================================================================================#

# << Ectomycorrhizal fungi >> -----
# add longitudinal data to em.site
for(i in em.site$Site){
  em.site[em.site$Site == i, 'long'] <- em.env[em.env$Site == i, 'long']
}
em.site <- em.site[c(1:16,298, 17:297)] # reorganize column order

# isolate community data
comm.matrix <- em.site[18:length(em.site)]

# Remove rare OTU (i.e. OTU with < 4 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]

# Hellinger transform data
comm.matrix.hell <- decostand(comm.matrix, method = 'hellinger')

# isolate geographic data
em.xy <- em.site[c('lat','long')]
em.xy.matrix <- as.matrix(em.xy)

# Final code for creation of the spatial weighting matrix
library(adespatial);library(sp);library(spdep)
nb <- chooseCN(coordinates(em.xy.matrix), type = 1, plot.nb = FALSE)
distnb <- nbdists(nb, em.xy.matrix)
fdist <- lapply(distnb, function(x) 1 - x/max(dist(em.xy.matrix)))
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)
 
# isolate environmental data
em.env.mem <- em.site[c('BIO1','BIO12','BIO16','BIO17')]

#--Multispati analysis
# when multivariate data are considered, it isi possible to search for spatial
# structures by computing the univariate statistics (e.g. Moran's Coefficient)
# on each variable separately. Another alternative is to summarize data by 
# multivariate methods and then detect spatial structures using the output
# of the analysis. For instance, we applied a centered PCA onthe abundance data.
pca.hell <- dudi.pca(comm.matrix.hell, scale = F, scannf = F, nf = 2)
moran.randtest(pca.hell$li, listw = lw)

ms.hell <- multispati(pca.hell, listw = lw, scannf = F)
summary(ms.hell)

#--Selection of spatial weighting matrix (SWM) and Moran's Eigenvector Maps (MEM)
mem.gab.sel <- mem.select(pca.hell$tab, listw = lw)
mem.gab.sel$global.test
mem.gab.sel$summary

cand.lw <- listw.candidates(em.xy.matrix, nb = c("del", "gab"),
                            weights = c("bin", "flin"))
sel.lw <- listw.select(pca.hell$tab,
                       candidates = cand.lw,
                       nperm = 1000)
sel.lw$candidates
lw.best <- cand.lw[[sel.lw$best.id]]

#--Canonical analysis
rda.hell <- pcaiv(pca.hell, sel.lw$best$MEM.select, scannf = F)
test.rda <- randtest(rda.hell)
test.rda ## Variation explained by the spatial predictors (R^2) is highly sig,
plot(test.rda)

#--Variation partitioning
# isolate environmental data
em.env.mem <- em.site[c('BIO1','BIO12','BIO16','BIO17','forest.type')]
em.clim.mem <- em.site[c('BIO1','BIO12','BIO16','BIO17')]
em.for.mem <- em.site[c('forest.type')]

vp1.em <- varpart(Y = pca.hell$tab,
                  em.clim.mem,
                  em.for.mem,
                  sel.lw$best$MEM.select)
vp1.em
plot(vp1.em, bg = c(3,5), Xnames = c('climate','forest.type','spatial'))

# << Endophytic fungi: culture-free data >> -----
# add longitudinal data to fe.cf.site
for(i in fe.cf.site$Site){
  fe.cf.site[fe.cf.site$Site == i, 'long'] <- fe.cf.env[fe.cf.env$site == i, 'long']
}
fe.cf.site <- fe.cf.site[c(1:3,1197,4:1196)] # reorganize column order

# isolate community data
comm.matrix <- fe.cf.site[16:length(fe.cf.site)]

# Remove rare OTU (i.e. OTU with < 8 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 8]

# Hellinger transform data
comm.matrix.hell <- decostand(comm.matrix, method = 'hellinger')

# isolate geographic data
fe.cf.xy <- fe.cf.site[c('lat','long')]
fe.cf.xy.matrix <- as.matrix(fe.cf.xy)

# Final code for creation of the spatial weighting matrix
library(adespatial);library(sp);library(spdep)
nb <- chooseCN(coordinates(fe.cf.xy.matrix), type = 1, plot.nb = FALSE)
distnb <- nbdists(nb, fe.cf.xy.matrix)
fdist <- lapply(distnb, function(x) 1 - x/max(dist(fe.cf.xy.matrix)))
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)
 
# isolate environmental data
fe.cf.env.mem <- fe.cf.site[c('BIO1','BIO12','BIO16','BIO17')]

#--Multispati analysis
# when multivariate data are considered, it isi possible to search for spatial
# structures by computing the univariate statistics (e.g. Moran's Coefficient)
# on each variable separately. Another alternative is to summarize data by 
# multivariate methods and then detect spatial structures using the output
# of the analysis. For instance, we applied a centered PCA onthe abundance data.

# Check autocorrelation of community and environment
pca.hell <- dudi.pca(comm.matrix.hell, scale = F, scannf = F, nf = 2)
pca.env <- dudi.pca(fe.cf.env.mem, scale = F, scannf = F, nf = 2)
moran.randtest(pca.hell$li, listw = lw)
moran.randtest(fe.cf.env.mem, lw)

ms.hell <- multispati(pca.hell, listw = lw, scannf = F)
summary(ms.hell)
plot(ms.hell)

env.hell <- multispati(pca.env, listw = lw, scannf = F)
summary(env.hell)
plot(env.hell)

#--Selection of spatial weighting matrix (SWM) and Moran's Eigenvector Maps (MEM)
mem.gab.sel <- mem.select(pca.hell$tab, listw = lw)
mem.gab.sel$global.test
mem.gab.sel$summary

cand.lw <- listw.candidates(fe.cf.xy.matrix, nb = c("del", "gab"),
                            weights = c("bin", "flin"))
sel.lw <- listw.select(pca.hell$tab,
                       candidates = cand.lw,
                       nperm = 1000)
sel.lw$candidates
lw.best <- cand.lw[[sel.lw$best.id]]

#--Canonical analysis
rda.hell <- pcaiv(pca.hell, sel.lw$best$MEM.select, scannf = F)
test.rda <- randtest(rda.hell)
test.rda ## Variation explained by the spatial predictors (R^2) is highly sig,
plot(test.rda)

#--Variation partitioning
# isolate environmental data
fe.cf.env.mem <- fe.cf.site[c('BIO1','BIO12','BIO16','BIO17','forest')]
fe.cf.clim.mem <- fe.cf.site[c('BIO1','BIO12','BIO16','BIO17')]
fe.cf.for.mem <- fe.cf.site[c('forest')]

vp1.fe.cf <- varpart(Y = pca.hell$tab,
                     fe.cf.clim.mem,
                     fe.cf.for.mem,
                     sel.lw$best$MEM.select)
vp1.fe.cf
plot(vp1.fe.cf, bg = c(3,5), Xnames = c('climate','forest','spatial'))

# << Endophytic fungi: culture-based data >> -----
# add longitudinal data to fe.cb.site
for(i in fe.cb.site$site){
  fe.cb.site[fe.cb.site$site == i, 'long'] <- fe.cb.env[fe.cb.env$site == i, 'long']
}
fe.cb.site <- fe.cb.site[c(1:3,110,4:109)] # reorganize column order

# isolate community data
comm.matrix <- fe.cb.site[13:length(fe.cb.site)]

# Remove rare OTU (i.e. OTU with < 4 occurrences)
comm.matrix <- comm.matrix[colSums(comm.matrix) >= 4]

# Hellinger transform data
comm.matrix.hell <- decostand(comm.matrix, method = 'hellinger')

# isolate geographic data
fe.cb.xy <- fe.cb.site[c('lat','long')]
fe.cb.xy.matrix <- as.matrix(fe.cb.xy)

# Final code for creation of the spatial weighting matrix
library(adespatial);library(sp);library(spdep)
nb <- chooseCN(coordinates(fe.cb.xy.matrix), type = 1, plot.nb = FALSE)
distnb <- nbdists(nb, fe.cb.xy.matrix)
fdist <- lapply(distnb, function(x) 1 - x/max(dist(fe.cb.xy.matrix)))
lw <- nb2listw(nb, style = 'W', glist = fdist, zero.policy = TRUE)
 
# isolate environmental data
fe.cb.env.mem <- fe.cb.site[c('BIO1','BIO12','BIO16','BIO17')]

#--Multispati analysis
# when multivariate data are considered, it isi possible to search for spatial
# structures by computing the univariate statistics (e.g. Moran's Coefficient)
# on each variable separately. Another alternative is to summarize data by 
# multivariate methods and then detect spatial structures using the output
# of the analysis. For instance, we applied a centered PCA on the abundance data.
pca.hell <- dudi.pca(comm.matrix.hell, scale = F, scannf = F, nf = 2)
moran.randtest(pca.hell$li, listw = lw)

ms.hell <- multispati(pca.hell, listw = lw, scannf = F)
summary(ms.hell)

#--Selection of spatial weighting matrix (SWM) and Moran's Eigenvector Maps (MEM)
mem.gab.sel <- mem.select(pca.hell$tab, listw = lw)
mem.gab.sel$global.test
mem.gab.sel$summary

cand.lw <- listw.candidates(fe.cb.xy.matrix, nb = c("del", "gab"),
                            weights = c("bin", "flin"))
sel.lw <- listw.select(pca.hell$tab,
                       candidates = cand.lw,
                       nperm = 1000)
sel.lw$candidates
lw.best <- cand.lw[[sel.lw$best.id]]

#--Canonical analysis
rda.hell <- pcaiv(pca.hell, sel.lw$best$MEM.select, scannf = F)
test.rda <- randtest(rda.hell)
test.rda ## Variation explained by the spatial predictors (R^2) is highly sig,
plot(test.rda)

#--Variation partitioning
# isolate environmental data
fe.cb.env.mem <- fe.cb.site[c('BIO1','BIO12','BIO16','BIO17','forest')]
fe.cb.clim.mem <- fe.cb.site[c('BIO1','BIO12','BIO16','BIO17')]
fe.cb.for.mem <- fe.cb.site[c('forest')]

vp1.fe.cb <- varpart(Y = pca.hell$tab, 
                     fe.cb.clim.mem,
                     fe.cb.for.mem,
                     sel.lw$best$MEM.select)
vp1.fe.cb
plot(vp1.fe.cb, bg = c(3,5), Xnames = c('climate','forest type','spatial'))
