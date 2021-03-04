## Script created by Liz Bowman January 6, 2020
## for analyzing community data in Sky Island study
## Includes ectomycorrhizal data and endophyte data

## For answering 'To what extent are EM and FE dispersal limited?
# Part A focuses at Distance decay
# 1. within isolated montane forests (mogollon rim removed)
# 2. between isolated montane forests (mogollon rim removed)
# 3. within a contiguous forest (mogollon rim only, EM only for this analysis)
# The analysis of distance decay was done in Jmp using Spearman's rho.

# Part B focuses on parsing out whether dispersal limitation is due to distance alone or
#   due to environmental decay
# 1. isolated montane forests (mogollon rim removed), unable to separate local and regional**
# 2. within a contiguous forest (mogollon rim only, em only)

# dbMEM.R script presents an alternative statistical method to partial mantel 

#========================================================================================#
# Load data------------------
#========================================================================================#

#--Site x species matrix with tip abundance per OTU
em.site <- read.csv(paste0(dat.dir,'EM_sitexspecies_TipAb_Site.csv'), as.is = T)

#--Site x species matrix ITS2 rarefied
fe.cf.site <- read.csv(paste0(dat.dir,'FE_ITS2rarified_sitexspecies.csv'), as.is = T)

#--Site x species matrix culture-based FE
fe.cb.site <- read.csv(paste0(dat.dir, 'FE_CB_sitexspecies.csv'), as.is = T)

#--Pairwise Dissimilarity for EM fungi
dist.decay.em <- read.csv(paste0(dat.dir, 'EM_PairwiseDissimilarityMatrix_site.csv'), as.is = T)
# for(i in unique(dist.decay.em$site2)){
#   dist.decay.em[dist.decay.em$site2 == i, 'range.2'] <- em.site[em.site$Site == i, 'Range']
# }
# #--assess if within same or different range origin
# for(i in 1:nrow(dist.decay.em)){
#   if(dist.decay.em[i, 'range.1'] == dist.decay.em[i, 'range.2']){
#     dist.decay.em[i, 'comp.range'] <- 'Same'
#   } else{dist.decay.em[i, 'comp.range'] <- 'Different'}
# }

#--Pairwise Dissimilarity for culture-free FE
dist.decay.cf <- read.csv(paste0(dat.dir, 'FE_CF_PairwiseDissimilarityMatrix_site_ITS2rarefied.csv'),
                          as.is = T)
# for(i in unique(dist.decay.fe$site2)){
#   dist.decay.fe[dist.decay.fe$site2 == i, 'range.2'] <- fe.site[fe.site$Site == i, 'range']
# }
# #--assess if within same or different range origin
# for(i in 1:nrow(dist.decay.fe)){
#   if(dist.decay.fe[i, 'range.1'] == dist.decay.fe[i, 'range.2']){
#     dist.decay.fe[i, 'comp.range'] <- 'Same'
#   } else{dist.decay.fe[i, 'comp.range'] <- 'Different'}
# }

#--Pairwise Dissimilarity for culture-based FE
dist.decay.cb <- read.csv(paste0(dat.dir, 'FE_CB_PairwiseDissimilarityMatrix_site.csv'),
                          as.is = T)

#========================================================================================#
# Part A: Distance decay ----
# analyses completed in JMP using Spearman's rho
# Below is the code to create Fig. 2
#========================================================================================#

#--------------------------------------------------------------------#
# Ectomycorrhizal fungi ----
#--------------------------------------------------------------------#

# << Distance decay: overall, Mogollon Rim removed  >> ------
# Remove any distances that are equal to 0
dist.decay.em.overall <- dist.decay.em[dist.decay.em$spatial.site > 0,]
# Remove Mogollon Rim
dist.decay.em.overall <- dist.decay.em.overall[!dist.decay.em.overall$range.1 == 'MogollonRim',]
dist.decay.em.overall <- dist.decay.em.overall[!dist.decay.em.overall$range.2 == 'MogollonRim',]
# Regional
dist.decay.em.regional <- dist.decay.em.overall[dist.decay.em.overall$comp.range == 'Different',]
# Local
dist.decay.em.local <- dist.decay.em.overall[dist.decay.em.overall$comp.range == 'Same',]
# Remove outliers
dist.decay.em.overall <- dist.decay.em.overall[dist.decay.em.overall$jaccard.dissimilarity > 0.7,]

#--plot Jaccard local
dist.decay.em.local %>% 
  filter(jaccard.dissimilarity >= 0.85) -> dist.decay.em.local.out
ggplot(data = dist.decay.em.local.out,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey', size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.85,1.0) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_EM_distdecayJaccard_local.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Jaccard regional
dist.decay.em.regional %>% 
  filter(jaccard.dissimilarity >= 0.92) -> dist.decay.em.regional.out
ggplot(data = dist.decay.em.regional.out,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey', size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.85,1.0) +
  xlim(0, 500) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_EM_distdecayJaccard_regional.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Jaccard overall
dist.decay.em.overall %>%
  filter(jaccard.dissimilarity > 0.925) -> dist.decay.em.overall.out
ggplot(data = dist.decay.em.overall.out,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey', size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  xlim(0, 500) +
   ylim(0.85,1.0) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_EM_distdecayJaccard_overall.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita-Horn local
dist.decay.em.local %>%
  filter(horn.dissimilarity > 0.825) -> dist.decay.em.local.out
ggplot(data = dist.decay.em.local,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey') +
  ylab("Morisita-Horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.85,1.0) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_distdecayHorn_local.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita-horn regional
dist.decay.em.regional %>%
  filter(horn.dissimilarity > 0.875) -> dist.decay.em.regional.out
ggplot(data = dist.decay.em.regional.out,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey') +
  ylab("Morisita-horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  ylim(0.85,1.0) +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_distdecayHorn_regional.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita horn overall
dist.decay.em.overall %>%
  filter(horn.dissimilarity > 0.85) -> dist.decay.em.overall.out
ggplot(data = dist.decay.em.overall,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey') +
  ylab("Morisita-horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  ylim(0.85,1.0) +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('EM_distdecayHorn_overall.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

# << Distance decay: Mogollon Rim (Fig. 3)  >> ------
# Remove any distances that are equal to 0
dist.decay.em.mogollon <- dist.decay.em[dist.decay.em$spatial.site > 0,]
# Remove Mogollon Rim
dist.decay.em.mogollon <- dist.decay.em.mogollon[dist.decay.em.mogollon$range.1 == 'MogollonRim',]
dist.decay.em.mogollon <- dist.decay.em.mogollon[dist.decay.em.mogollon$range.2 == 'MogollonRim',]
# Remove outliers
dist.decay.em.mogollon <- dist.decay.em.mogollon[dist.decay.em.mogollon$horn.dissimilarity > 0.45,]

#--plot Jaccard
ggplot(data = dist.decay.em.mogollon,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  ylim(0.5,1.0) +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig3A_EM_distdecayJaccard_Mogollon.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita horn
ggplot(data = dist.decay.em.mogollon,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Morisita-horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  ylim(0.5,1.0) +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig3B_EM_distdecayHorn_Mogollon.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--------------------------------------------------------------------#
# Endophytic fungi ----
#--------------------------------------------------------------------#
# << Distance decay: Culture-free overall  >> ------
# Remove any distances that are equal to 0
dist.decay.cf.overall <- dist.decay.cf[dist.decay.cf$spatial.site > 0,]
# Remove Mogollon Rim
dist.decay.cf.overall <- dist.decay.cf.overall[!dist.decay.cf.overall$range.1 == 'MogollonRim',]
dist.decay.cf.overall <- dist.decay.cf.overall[!dist.decay.cf.overall$range.2 == 'MogollonRim',]
# Regional
dist.decay.cf.regional <- dist.decay.cf.overall[dist.decay.cf.overall$comp.range == 'Different',]
# Remove outliers
#dist.decay.cf.overall <- dist.decay.cf.overall[dist.decay.cf.overall$jaccard.dissimilarity > 0.6,]
# Local
dist.decay.cf.local <- dist.decay.cf.overall[dist.decay.cf.overall$comp.range == 'Same',]

#--plot Jaccard regional
dist.decay.cf.regional %>%
  filter(jaccard.dissimilarity >= 0.80) -> dist.decay.cf.regional.out
ggplot(data = dist.decay.cf.regional.out,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey', size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.4,1.0) +
  xlim(0,400) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_FE_distdecayJaccard_CFoverall_regional.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Jaccard local
ggplot(data = dist.decay.cf.local,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey', size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.4,1.0) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_FE_distdecayJaccard_CFoverall_local.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Jaccard overall
dist.decay.cf.overall %>%
  filter(jaccard.dissimilarity >= 0.80) -> dist.decay.cf.overall.out
ggplot(data = dist.decay.cf.overall,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey', size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  ylim(0.4,1.0) +
  xlim(0,400) +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_FE_distdecayJaccard_CFoverall.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

# << Distance decay: Culture-based overall  >> ------
# Remove any distances that are equal to 0
dist.decay.cb.overall <- dist.decay.cb[dist.decay.cb$spatial.site > 0,]
# Remove Mogollon Rim
dist.decay.cb.overall <- dist.decay.cb.overall[!dist.decay.cb.overall$range.1 == 'MogollonRim',]
dist.decay.cb.overall <- dist.decay.cb.overall[!dist.decay.cb.overall$range.2 == 'MogollonRim',]
# Regional
dist.decay.cb.regional <- dist.decay.cb.overall[dist.decay.cb.overall$range.comp == 'Different',]
# Remove outliers
#dist.decay.cf.overall <- dist.decay.cf.overall[dist.decay.cf.overall$jaccard.dissimilarity > 0.6,]
# Local
dist.decay.cb.local <- dist.decay.cb.overall[dist.decay.cb.overall$range.comp == 'Same',]

#--plot Jaccard
dist.decay.cb.regional %>%
  filter(jaccard.dissimilarity > 0.60) -> dist.decay.cb.regional.out
ggplot(data = dist.decay.cb.regional.out,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.4,1.0) +
  xlim(0,400) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_FE_distdecayJaccard_CBoverall_regional.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Jaccard local
ggplot(data = dist.decay.cb.local,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.4,1.0) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_FE_distdecayJaccard_CBoverall_local.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Jaccard overall
dist.decay.cb.overall %>%
  filter(jaccard.dissimilarity > 0.5) -> dist.decay.cb.overall.out
ggplot(data = dist.decay.cb.overall.out,
       aes(x = spatial.site,
           y = jaccard.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Jaccard dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  ylim(0.4,1.0) +
  xlim(0,400) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=24, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('Fig2_FE_distdecayJaccard_CBoverall.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita-horn
dist.decay.cb.regional %>%
  filter(horn.dissimilarity >= 0.1) -> dist.decay.cb.regional.out
ggplot(data = dist.decay.cb.regional.out,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  geom_smooth(method = 'lm', se = F, color = 'darkgrey') +
  ylab("Morisita-Horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  ylim(0.0,1.0) +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('FigS1_FE_distdecayHorn_CBoverall_regional.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita-horn local
ggplot(data = dist.decay.cb.local,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Morisita-Horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('FigS1_FE_distdecayHorn_CBoverall_local.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#--plot Morisita-horn overall
ggplot(data = dist.decay.cb.overall,
       aes(x = spatial.site,
           y = horn.dissimilarity)) +
  geom_point(size = 3) +
  ylab("Morisita-Horn dissimilarity") +
  xlab('Pairwise distance between sites (km)') +
  theme_bw() +
  xlab(element_blank()) +
  #ggtitle("Proportion of Classes by Topography") +
  guides (fill=guide_legend(title=NULL)) +
  theme(legend.position='right',
        # axis.title.x = element_text(margin = margin(t = 30)),
        axis.title.y = element_text(margin = margin(r = 30)),
        #axis.text.x = element_text(angle = 35),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(), axis.line = element_line(),
        axis.text = element_text(size=20, color = 'black'),
        axis.title = element_text(size = 28),
        strip.text.x = element_text(size = 14))

ggsave('FigS1_FE_distdecayHorn_CBoverall.pdf', plot = last_plot(), 
       device = 'pdf', path = './figures/',
       units = 'in', width = 10, height = 7)

#========================================================================================#
# Part B: Partial Mantel ----
#========================================================================================#

# The output of this data can be found in Table 3 and Supplementary Table S5.

#--------------------------------------------------------------------#
# Ectomycorrhizal fungi ----
#--------------------------------------------------------------------#

#--Mogollon Rim sites
mogollon.rim <- c('M1','M2','M3','M4','M5','M6','M7','M8','M9')

#--climate distance matrix
em.clim <- read.csv(paste0(dat.dir, 'EM_EnvFactors_site.csv'), as.is = T)
# non mogollon rim
em.clim.non <- em.clim[!em.clim$Site %in% mogollon.rim,]
em.clim.dist.non <- vegdist(em.clim.non$clim.pca, method = 'euclidean')
# mogollon rim only
em.clim.mr <- em.clim[em.clim$Site %in% mogollon.rim,]
em.clim.dist.mr <- vegdist(em.clim.mr$clim.pca, method = 'euclidean')

#--Environmental matrix
# non mogollon rim
em.env.non <- data.frame(clim.pca = em.clim.non$clim.pca,
                     forest.type = as.factor(em.clim.non$forest))
em.env.dist.non <- daisy(em.env.non, metric = 'gower')
# mogollon rim only
em.env.mr <- data.frame(clim.pca = em.clim.mr$clim.pca,
                        forest.type = as.factor(em.clim.mr$forest))
em.env.dist.mr <- daisy(em.env.mr, metric = 'gower')

#--spatial matrix
em.spatial <- read.csv(paste0(dat.dir, 'DistanceMatrix/EM_GeoDistanceMatrix_site.csv'),
                       row.names = 1)
# remove mogollon rim
em.spatial.non <- em.spatial[!rownames(em.spatial) %in% mogollon.rim,]
em.spatial.non <- em.spatial.non[!colnames(em.spatial.non) %in% mogollon.rim]
# isolate mogollon rim
em.spatial.mr <- em.spatial[rownames(em.spatial) %in% mogollon.rim,]
em.spatial.mr <- em.spatial.mr[colnames(em.spatial.mr) %in% mogollon.rim]

#--EM site x species matrix, remove mogollon rim
em.site.non <- em.site[!em.site$Site %in% mogollon.rim,]
# isolate community data, with mogollon rim removed
em.cm <- em.site.non[17:length(em.site.non)]
# remove species with less than 2 occurrences
em.cm <- em.cm[colSums(em.cm) > 2]
# Jaccard
jaccard.em.dist <- vegdist(em.cm,method = 'jaccard')
# Morisita horn
horn.em.dist <- vegdist(em.cm, method = 'horn')

#--EM site x species matrix, isolate mogollon rim
em.site.mr <- em.site[em.site$Site %in% mogollon.rim,]
# isolate community data, with mogollon rim isolated
em.cm.mr <- em.site.mr[17:length(em.site.mr)]
# remove species with less than 2 occurrences
em.cm.mr <- em.cm.mr[colSums(em.cm.mr) > 2]
# Jaccard
jaccard.em.dist.mr <- vegdist(em.cm.mr,method = 'jaccard')
# Morisita horn
horn.em.dist.mr <- vegdist(em.cm.mr, method = 'horn')

#--<< Isolated communities only, with Mogollon Rim removed >> ------------
partial.mantel.noMR <- data.frame(test = c('mp.em.spatial.envrem.jaccard','mo.em.spatial.envrem.horn',
                                           'mp.em.spatial.climrem.jaccard','mo.em.spatial.climrem.horn',
                                           'mp.em.clim.spatialrem.jaccard','mo.em.clim.spatialrem.horn',
                                           'mp.em.env.spatialrem.jaccard','mo.em.env.spatialrem.horn'),
                                  mantel.stat.r = NA,
                                  p.value = NA)

#--spatial with effect of environment removed
mp.em.spatial.envrem.jaccard <- mantel.partial(jaccard.em.dist, em.spatial.non, em.env.dist.non)
partial.mantel.noMR[1, 'mantel.stat.r'] <- mp.em.spatial.envrem.jaccard$statistic
partial.mantel.noMR[1, 'p.value'] <- mp.em.spatial.envrem.jaccard$signif

mo.em.spatial.envrem.horn <- mantel.partial(horn.em.dist, em.spatial.non, em.env.dist.non)
partial.mantel.noMR[2, 'mantel.stat.r'] <- mo.em.spatial.envrem.horn$statistic
partial.mantel.noMR[2, 'p.value'] <- mo.em.spatial.envrem.horn$signif

#--spatial with effect of climate removed
mp.em.spatial.climrem.jaccard <- mantel.partial(jaccard.em.dist, em.spatial.non, em.clim.dist.non)
partial.mantel.noMR[3, 'mantel.stat.r'] <- mp.em.spatial.climrem.jaccard$statistic
partial.mantel.noMR[3, 'p.value'] <- mp.em.spatial.climrem.jaccard$signif

mo.em.spatial.climrem.horn <- mantel.partial(horn.em.dist, em.spatial.non, em.clim.dist.non)
partial.mantel.noMR[4, 'mantel.stat.r'] <- mo.em.spatial.climrem.horn$statistic
partial.mantel.noMR[4, 'p.value'] <- mo.em.spatial.climrem.horn$signif

#--climate with effect of spatial removed
mp.em.clim.spatialrem.jaccard <- mantel.partial(jaccard.em.dist, em.clim.dist.non, em.spatial.non)
partial.mantel.noMR[5, 'mantel.stat.r'] <- mp.em.clim.spatialrem.jaccard$statistic
partial.mantel.noMR[5, 'p.value'] <- mp.em.clim.spatialrem.jaccard$signif

mo.em.clim.spatialrem.horn <- mantel.partial(horn.em.dist, em.clim.dist.non, em.spatial.non)
partial.mantel.noMR[6, 'mantel.stat.r'] <- mo.em.clim.spatialrem.horn$statistic
partial.mantel.noMR[6, 'p.value'] <- mo.em.clim.spatialrem.horn$signif

#--environment with effect of spatial removed
mp.em.env.spatialrem.jaccard <- mantel.partial(jaccard.em.dist, em.env.dist.non, em.spatial.non)
partial.mantel.noMR[7, 'mantel.stat.r'] <- mp.em.env.spatialrem.jaccard$statistic
partial.mantel.noMR[7, 'p.value'] <- mp.em.env.spatialrem.jaccard$signif

mo.em.env.spatialrem.horn <- mantel.partial(horn.em.dist, em.env.dist.non, em.spatial.non)
partial.mantel.noMR[8, 'mantel.stat.r'] <- mo.em.env.spatialrem.horn$statistic
partial.mantel.noMR[8, 'p.value'] <- mo.em.env.spatialrem.horn$signif

write.csv(partial.mantel.noMR, paste0(res.dir,'EM_PartialMantel_WinBtwn_NoMR.csv'), row.names = F)

#--<< Connected communities only, with Mogollon Rim only >> ------------
partial.mantel.MR <- data.frame(test = c('mp.em.spatial.envrem.jaccard.mr','mo.em.spatial.envrem.horn.mr',
                                           'mp.em.spatial.climrem.jaccard.mr','mo.em.spatial.climrem.horn.mr',
                                           'mp.em.clim.spatialrem.jaccard.mr','mo.em.clim.spatialrem.horn.mr',
                                           'mp.em.env.spatialrem.jaccard.mr','mo.em.env.spatialrem.horn.mr'),
                                  mantel.stat.r = NA,
                                  p.value = NA)

#--spatial with effect of environment removed
mp.em.spatial.envrem.jaccard.mr <- mantel.partial(jaccard.em.dist.mr, em.spatial.mr, em.env.dist.mr)
partial.mantel.MR[1, 'mantel.stat.r'] <- mp.em.spatial.envrem.jaccard.mr$statistic
partial.mantel.MR[1, 'p.value'] <- mp.em.spatial.envrem.jaccard.mr$signif

mo.em.spatial.envrem.horn.mr <- mantel.partial(horn.em.dist.mr, em.spatial.mr, em.env.dist.mr)
partial.mantel.MR[2, 'mantel.stat.r'] <- mo.em.spatial.envrem.horn.mr$statistic
partial.mantel.MR[2, 'p.value'] <- mo.em.spatial.envrem.horn.mr$signif

#--spatial with effect of climate removed
mp.em.spatial.climrem.jaccard.mr <- mantel.partial(jaccard.em.dist.mr, em.spatial.mr, em.clim.dist.mr)
partial.mantel.MR[3, 'mantel.stat.r'] <- mp.em.spatial.climrem.jaccard.mr$statistic
partial.mantel.MR[3, 'p.value'] <- mp.em.spatial.climrem.jaccard.mr$signif

mo.em.spatial.climrem.horn.mr <- mantel.partial(horn.em.dist.mr, em.spatial.mr, em.clim.dist.mr)
partial.mantel.MR[4, 'mantel.stat.r'] <- mo.em.spatial.climrem.horn.mr$statistic
partial.mantel.MR[4, 'p.value'] <- mo.em.spatial.climrem.horn.mr$signif

#--climate with effect of spatial removed
mp.em.clim.spatialrem.jaccard.mr <- mantel.partial(jaccard.em.dist.mr, em.clim.dist.mr, em.spatial.mr)
partial.mantel.MR[5, 'mantel.stat.r'] <- mp.em.clim.spatialrem.jaccard.mr$statistic
partial.mantel.MR[5, 'p.value'] <- mp.em.clim.spatialrem.jaccard.mr$signif

mo.em.clim.spatialrem.horn.mr <- mantel.partial(horn.em.dist.mr, em.clim.dist.mr, em.spatial.mr)
partial.mantel.MR[6, 'mantel.stat.r'] <- mo.em.clim.spatialrem.horn.mr$statistic
partial.mantel.MR[6, 'p.value'] <- mo.em.clim.spatialrem.horn.mr$signif

#--environment with effect of spatial removed
mp.em.env.spatialrem.jaccard.mr <- mantel.partial(jaccard.em.dist.mr, em.env.dist.mr, em.spatial.mr)
partial.mantel.MR[7, 'mantel.stat.r'] <- mp.em.env.spatialrem.jaccard.mr$statistic
partial.mantel.MR[7, 'p.value'] <- mp.em.env.spatialrem.jaccard.mr$signif

mo.em.env.spatialrem.horn.mr <- mantel.partial(horn.em.dist.mr, em.env.dist.mr, em.spatial.mr)
partial.mantel.MR[8, 'mantel.stat.r'] <- mo.em.env.spatialrem.horn.mr$statistic
partial.mantel.MR[8, 'p.value'] <- mo.em.env.spatialrem.horn.mr$signif

write.csv(partial.mantel.MR, paste0(res.dir,'EM_PartialMantel_MRonly.csv'), row.names = F)

#--------------------------------------------------------------------#
# Endophytic fungi: Culture-free ----
#--------------------------------------------------------------------#

#--Mogollon Rim sites
mogollon.rim <- c('M1','M2','M3','P3')

#--climate distance matrix
fe.cf.env <- read.csv(paste0(data.dir,'FE_CF_EnvFactors_site.csv'))
# mogollon rim removed
fe.cf.env.non <- fe.cf.env[!fe.cf.env$site %in% mogollon.rim,]
fe.clim.dist.non <- vegdist(fe.cf.env.non$clim.pca, method = 'euclidean')å

#--Environmental matrix 
# mogollon rim removed
fe.env.non <- data.frame(clim.pca = fe.cf.env.non$clim.pca,
                         forest.type = as.factor(fe.cf.env.non$forest))
fe.env.dist.non <- daisy(fe.env.non, metric = 'gower')

#--spatial matrix
fe.cf.spatial <- read.csv(paste0(dat.dir,'DistanceMatrix/FE_CF_GeoDistanceMatrix_site.csv'),
                          row.names = 1)
# Mogollon Rim removed
fe.spatial.non <- fe.cf.spatial[!colnames(fe.cf.spatial) %in% mogollon.rim]
fe.spatial.non <- fe.spatial.non[!rownames(fe.spatial.non) %in% mogollon.rim, ]

#--EM site x species matrix, remove mogollon rim
fe.cf.site.non <- fe.cf.site[!fe.cf.site$Site %in% mogollon.rim,]
# isolate community data, with mogollon rim removed
fe.cm <- fe.cf.site.non[15:length(fe.cf.site.non)]
# remove species with less than 8 occurrences
fe.cm <- fe.cm[colSums(fe.cm) > 8]
# Jaccard
jaccard.fe.dist <- vegdist(fe.cm,method = 'jaccard')

#--<< Isolated communities only, with Mogollon Rim removed >> ------------
partial.mantel.fe.noMR <- data.frame(test = c('mp.fe.spatial.envrem.jaccard',
                                           'mp.fe.spatial.climrem.jaccard',
                                           'mp.fe.clim.spatialrem.jaccard',
                                           'mp.fe.env.spatialrem.jaccard'),
                                  mantel.stat.r = NA,
                                  p.value = NA)

#--spatial with effect of environment removed
mp.fe.spatial.envrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.spatial.non, fe.env.dist.non)
partial.mantel.fe.noMR[1, 'mantel.stat.r'] <- mp.fe.spatial.envrem.jaccard$statistic
partial.mantel.fe.noMR[1, 'p.value'] <- mp.fe.spatial.envrem.jaccard$signif

#--spatial with effect of climate removed
mp.fe.spatial.climrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.spatial.non, fe.clim.dist.non)
partial.mantel.fe.noMR[2, 'mantel.stat.r'] <- mp.fe.spatial.climrem.jaccard$statistic
partial.mantel.fe.noMR[2, 'p.value'] <- mp.fe.spatial.climrem.jaccard$signif

#--climate with effect of spatial removed
mp.fe.clim.spatialrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.clim.dist.non, fe.spatial.non)
partial.mantel.fe.noMR[3, 'mantel.stat.r'] <- mp.fe.clim.spatialrem.jaccard$statistic
partial.mantel.fe.noMR[3, 'p.value'] <- mp.fe.clim.spatialrem.jaccard$signif

#--environment with effect of spatial removed
mp.fe.env.spatialrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.env.dist.non, fe.spatial.non)
partial.mantel.fe.noMR[4, 'mantel.stat.r'] <- mp.fe.env.spatialrem.jaccard$statistic
partial.mantel.fe.noMR[4, 'p.value'] <- mp.fe.env.spatialrem.jaccard$signif

write.csv(partial.mantel.fe.noMR, paste0(res.dir,'FE_CF_PartialMantel_WinBtwn_NoMR.csv'),
          row.names = F)

#--------------------------------------------------------------------#
# Endophytic fungi: Culture-based ----
#--------------------------------------------------------------------#

#--Mogollon Rim sites
mogollon.rim <- c('M1','M2','M3')

#--climate distance matrix
env.data <- read.csv(paste0(dat.dir,'FE_CB_EnvFactors_site.csv'))
# remove mogollon rim 
env.data.non <- env.data[!env.data$range == 'mogollon',]
fe.clim.dist.non <- vegdist(env.data.non['clim.pca'], method = 'euclidean')

#--Environmental distance matrix 
fe.env.non <- data.frame(clim.pca = env.data.non$clim.pca,
                         forest.type = as.factor(env.data.non$forest))
fe.env.dist.non <- daisy(fe.env.non, metric = 'gower')

#--spatial matrix
fe.spatial <- read.csv(paste0(dat.dir,'DistanceMatrix/FE_CB_GeoDistanceMatrix_site.csv'),
                       row.names = 1)
# Mogollon Rim removed
fe.spatial.non <- fe.spatial[!colnames(fe.spatial) %in% mogollon.rim]
fe.spatial.non <- fe.spatial.non[!rownames(fe.spatial.non) %in% mogollon.rim,]

#--FE site x species matrix, remove mogollon rim
fe.cb.site.non <- fe.cb.site[!fe.cb.site$site %in% mogollon.rim,]
rownames(fe.cb.site.non) <- fe.cb.site.non$site
# isolate community data, with mogollon rim removed
fe.cm <- fe.cb.site.non[12:length(fe.cb.site.non)]
# remove species with less than 8 occurrences
fe.cm <- fe.cm[colSums(fe.cm) >= 4]
# Jaccard
jaccard.fe.dist <- vegdist(fe.cm,method = 'jaccard')
# Morisita-Horn
horn.fe.dist <- vegdist(fe.cm, method = 'horn')

#--<< Isolated communities only, with Mogollon Rim removed >> ------------
partial.mantel.fe.noMR <- data.frame(test = c('mp.fe.cb.spatial.envrem.jaccard','mo.fe.cb.spatial.envrem.horn',
                                           'mp.fe.cb.spatial.climrem.jaccard','mo.fe.cb.spatial.climrem.horn',
                                           'mp.fe.cb.clim.spatialrem.jaccard','mo.fe.cb.clim.spatialrem.horn',
                                           'mp.fe.cb.env.spatialrem.jaccard','mo.fe.cb.env.spatialrem.horn'),
                                  mantel.stat.r = NA,
                                  p.value = NA)

#--spatial with effect of environment removed
mp.fe.spatial.envrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.spatial.non, fe.env.dist.non)
partial.mantel.fe.noMR[1, 'mantel.stat.r'] <- mp.fe.spatial.envrem.jaccard$statistic
partial.mantel.fe.noMR[1, 'p.value'] <- mp.fe.spatial.envrem.jaccard$signif

mp.fe.spatial.envrem.horn <- mantel.partial(horn.fe.dist, fe.spatial.non, fe.env.dist.non)
partial.mantel.fe.noMR[2, 'mantel.stat.r'] <- mp.fe.spatial.envrem.horn$statistic
partial.mantel.fe.noMR[2, 'p.value'] <- mp.fe.spatial.envrem.horn$signif

#--spatial with effect of climate removed
mp.fe.spatial.climrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.spatial.non, fe.clim.dist.non)
partial.mantel.fe.noMR[3, 'mantel.stat.r'] <- mp.fe.spatial.climrem.jaccard$statistic
partial.mantel.fe.noMR[3, 'p.value'] <- mp.fe.spatial.climrem.jaccard$signif

mp.fe.spatial.climrem.horn<- mantel.partial(horn.fe.dist, fe.spatial.non, fe.clim.dist.non)
partial.mantel.fe.noMR[4, 'mantel.stat.r'] <- mp.fe.spatial.climrem.horn$statistic
partial.mantel.fe.noMR[4, 'p.value'] <- mp.fe.spatial.climrem.horn$signif

#--climate with effect of spatial removed
mp.fe.clim.spatialrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.clim.dist.non, fe.spatial.non)
partial.mantel.fe.noMR[5, 'mantel.stat.r'] <- mp.fe.clim.spatialrem.jaccard$statistic
partial.mantel.fe.noMR[5, 'p.value'] <- mp.fe.clim.spatialrem.jaccard$signif

mp.fe.clim.spatialrem.horn <- mantel.partial(horn.fe.dist, fe.clim.dist.non, fe.spatial.non)
partial.mantel.fe.noMR[6, 'mantel.stat.r'] <- mp.fe.clim.spatialrem.horn$statistic
partial.mantel.fe.noMR[6, 'p.value'] <- mp.fe.clim.spatialrem.horn$signif

#--environment with effect of spatial removed
mp.fe.env.spatialrem.jaccard <- mantel.partial(jaccard.fe.dist, fe.env.dist.non, fe.spatial.non)
partial.mantel.fe.noMR[7, 'mantel.stat.r'] <- mp.fe.env.spatialrem.jaccard$statistic
partial.mantel.fe.noMR[7, 'p.value'] <- mp.fe.env.spatialrem.jaccard$signif

mp.fe.env.spatialrem.horn <- mantel.partial(horn.fe.dist, fe.env.dist.non, fe.spatial.non)
partial.mantel.fe.noMR[8, 'mantel.stat.r'] <- mp.fe.env.spatialrem.horn$statistic
partial.mantel.fe.noMR[8, 'p.value'] <- mp.fe.env.spatialrem.horn$signif

write.csv(partial.mantel.fe.noMR, paste0(fig.dir,'FE_CB_PartialMantel_WinBtwn_NoMR.csv'),
          row.names = F)
