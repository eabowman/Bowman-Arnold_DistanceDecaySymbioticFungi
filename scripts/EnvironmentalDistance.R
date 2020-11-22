## Script created by Liz Bowman Oct. 9, 2019
## for analysing creating pair-wise distance matrices

#========================================================================================#
# Load libraries and data ------------------
#========================================================================================#

#install.packages(c('dplyr','tidyr','vegan','ggplot2'))
library(dplyr);library(tidyr);library(vegan);library(ggplot2);library(nlme)

#--file paths
dat.dir <- '~/Documents/PhD/3:4_Combined/data/'
fig.dir <- '~/Documents/PhD/3:4_Combined/figures/'
res.dir <- '~/Documents/PhD/3:4_Combined/results/'

#--Load sampling site coordinates and climate data
clim.data <- read.csv(paste0(dat.dir, '20200107_ClimateData.csv'), as.is = T)

# Site level environmental data
# BIO1 = annual mean temperature
# BIO7 = Temperature annual range (max temp. of warmest month - max temp of coldest month)
# BIO12 = annual precipitation
# BIO16 = precipitation of wettest quarter
# BIO17 = precipitiation of driest quarter
clim.data %>%
  select(site, range, lat, long, masl, BIO1, BIO7, BIO12, BIO16, BIO17, forest) %>%
  distinct(.) -> clim.data.site

#========================================================================================#
# Environmental matrix: single factor and multifactor ------------------
#========================================================================================#

#<< Multifactor >> ----
# Site level ----
# filter FE sites
# fe.sites <- c('B1','B2','B3','C1','C2','C3','H1','H2','H3',
#               'M1','M2','M3','N1','N2','N3','P1','P2','P3',
#               'Sunset Trail','Upper Solders Camp Rd.','Butterfly Trail')
# clim.data.site <- clim.data.site[clim.data.site$site %in% fe.sites,]

#Remove those sites that are not in the EM data
scm.remove <- c('Bear Wallow', 'Middle Bear', 'Parking pullout South of Rose Canyon',
                'Green Mountain/Bug Springs','Marshall Gulch','Rose Canyon Rd.',
                'Pullout south of Willow Canyon Rd.')
clim.data.site <- clim.data.site[!clim.data.site$site %in% scm.remove,]

pca.env <- clim.data.site[c('BIO1','BIO12','BIO16','BIO17','forest')]
pca.env$forest <- as.numeric(factor(pca.env$forest, labels = c(1:3))) # make numeric
rownames(pca.env) <- clim.data.site$site

pca.env.analysis <- prcomp(pca.env,
                           center = TRUE,
                           scale. = TRUE)

summary(pca.env.analysis)
env.eigenvector <- scores(pca.env.analysis, choices = c(1:4))
env.eigenvector <- data.frame(env.eigenvector)
names(env.eigenvector) <- c('PCA1','PCA2')
env.eigenvector$site <- rownames(pca.env)

clim.data.site['env.pca'] <- NA
for(i in unique(clim.data.site$site)){
  clim.data.site[clim.data.site$site == i, 'env.pca'] <- 
    env.eigenvector[env.eigenvector$site == i, 'PCA1']
}

write.csv(clim.data.site, 'data/EM_EnvFactors_site.csv')

# distance matrix
multi.env.site <- dist(clim.data.site$env.pca, method = 'euclidean')
multi.env.site <- as.matrix(multi.env.site)
rownames(multi.env.site) <- clim.data.site$site
colnames(multi.env.site) <- clim.data.site$site

#write.csv(multi.env.site, 'data/DistanceMatrix/20200117_FE_EnvDistMatrix_site.csv',row.names = T)
write.csv(multi.env.site, 'data/DistanceMatrix/20200117_EM_EnvDistMatrix_site_pca.csv',row.names = T)

# Tree level ----
# Perform PCA of environmental factors
# filter tree level
# fe.trees <- c('B11','B12','B12','B21','B22','B31','C11','C12','C13','C14','C21',
#               'C22','C23','C23','C32','C33','C34','H11','H12','H21','H22','H31',
#               'LB031','LB033','LB048','M11','M12','M21','M21','M22','M31','N11',
#               'N24','N25','N32','P11','P12','P21','P22')
# clim.data <- clim.data[clim.data$tree %in% fe.trees,]

pca.env <- clim.data[c('BIO1','BIO12','BIO16','BIO17','forest')]
pca.env$forest <- as.numeric(factor(pca.env$forest, labels = c(1:2))) # make numeric
rownames(pca.env) <- clim.data$tree

pca.env.analysis <- prcomp(pca.env,
                           center = TRUE,
                           scale. = TRUE)

summary(pca.env.analysis)
env.eigenvector <- scores(pca.env.analysis, choices = c(1:4))
env.eigenvector <- data.frame(env.eigenvector)
names(env.eigenvector) <- c('PCA1','PCA2')
env.eigenvector$tree <- rownames(pca.env)

clim.data['env.pca'] <- NA
for(i in unique(clim.data$tree)){
  clim.data[clim.data$tree == i, 'env.pca'] <- 
    env.eigenvector[env.eigenvector$tree == i, 'PCA1']
}

write.csv(clim.data, 'data/FE_EnvFactors_tree.csv')

# distance matrix
multi.env.tree <- dist(clim.data$env.pca, method = 'euclidean')
multi.env.tree <- as.matrix(multi.env.tree)
rownames(multi.env.tree) <- clim.data$tree
colnames(multi.env.tree) <- clim.data$tree

write.csv(multi.env.tree, 'data/DistanceMatrix/20200107_FE_EnvDistMatrix_tree.csv',row.names = T)

#<< Climate >> ----
# Site level ----
# distance matrix
clim.site <- dist(clim.data.site[c('BIO1','BIO12','BIO16','BIO17')], method = 'euclidean')
clim.site <- as.matrix(clim.site)
rownames(clim.site) <- clim.data.site$site
colnames(clim.site) <- clim.data.site$site

#write.csv(clim.site, 'data/DistanceMatrix/20200107_FE_ClimateDistMatrix_site_nopca.csv',row.names = T)
write.csv(clim.site, 'data/DistanceMatrix/20200107_EM_ClimateDistMatrix_site_nopca.csv',row.names = T)

em.site.env <- read.csv('data/EM_EnvFactors_site.csv', as.is = T)
clim.site.em <- as.matrix(dist(em.site.env$clim.pca, method = 'euclidean'))
rownames(clim.site.em) <- em.site.env$site
colnames(clim.site.em) <- em.site.env$site

write.csv(clim.site.em, 'data/DistanceMatrix/20200112_EM_ClimateDistMatrix_site_pca.csv')

fe.site.env <- read.csv('data/FE_EnvFactors_site.csv', as.is = T)
clim.site.fe <- as.matrix(dist(fe.site.env$clim.pca, method = 'euclidean'))
rownames(clim.site.fe) <- fe.site.env$site
colnames(clim.site.fe) <- fe.site.env$site

write.csv(clim.site.fe, 'data/DistanceMatrix/20200112_FE_ClimateDistMatrix_site_pca.csv')

# Tree level ----
# distance matrix
clim.tree <- dist(clim.data[c('BIO1','BIO12','BIO16','BIO17')], method = 'euclidean')
clim.tree <- as.matrix(clim.tree)
rownames(clim.tree) <- clim.data$tree
colnames(clim.tree) <- clim.data$tree

write.csv(clim.tree, 'data/DistanceMatrix/20200107_ClimateDisteMatrix_tree.csv',row.names = T)
