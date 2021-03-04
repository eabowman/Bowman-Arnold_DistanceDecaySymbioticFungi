## Script created by Liz Bowman July 25, 2019
## for analysing mock community data of second Illumina run

#========================================================================================#
# ITS2 95% maxEE=1.0 : Load libraries and data for mock community analyses---------------
#========================================================================================#
#--Overall
mock.real <- read.csv('./data/MockCommReal.csv', 
                      as.is = T)

#--Data frame output from Usearch
mock.ITS2 <- read.csv('./data/MockCommunity_ITS2.csv',
                      as.is = T)

mock.ITS2 <- distinct(mock.ITS2, mock.member, .keep_all = T)

#--------------------------------------------------------------#
# Taxonomy by Phylum-----
#--------------------------------------------------------------#
# Expected mock community
mock.real$Phylum <- factor(mock.real$Phylum)
mock.input <- ggplot(mock.real, aes(x = Phylum)) +
  geom_bar() +
  ylab('Count') +
  xlab('') +
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
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave('FigS4A_mock.pdf', plot = mock.input, 
       device = 'pdf', path = fig.dir,
       units = 'in', width = 10, height = 7)

# Sequenced mock community
mock.ITS2$Phylum <- factor(mock.ITS2$Phylum,
                           levels = c('Ascomycota','Basidiomycota',
                                      'Chytridiomycota','Mucoromycota',
                                      'Zoopagomycota'))
mock.seq<- ggplot(mock.ITS2, aes(x = Phylum)) +
  geom_bar() +
  ylab('Count') +
  xlab('') +
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
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(drop=F)

ggsave('FigS4B_mock.pdf', plot = mock.seq, 
       device = 'pdf', path = fig.dir,
       units = 'in', width = 10, height = 7)

#--------------------------------------------------------------#
# Taxonomy by Class-----
#--------------------------------------------------------------#

mock.seq<- ggplot(mock.ITS2, aes(x = Class)) +
  geom_bar() +
  ylab('Count') +
  xlab('') +
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
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave('FigS4B_mock.pdf', plot = mock.seq, 
#        device = 'pdf', path = fig.dir,
#        units = 'in', width = 10, height = 7)

mock.input <- ggplot(mock.real, aes(x = Class)) +
  geom_bar() +
  ylab('Count') +
  xlab('') +
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
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(angle = 45, hjust = 1))

# ggsave('FigS4B_mock.pdf', plot = mock.input, 
#        device = 'pdf', path = fig.dir,
#        units = 'in', width = 10, height = 7)
