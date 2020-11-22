## Script created by Liz Bowman January 18, 2018
## for analysing community data in Sky Island study
## Comparison of culture-free and culture-based data to asses overlap between the two

#========================================================================================#
# Load data for comparison of culture-free and culture-based communities------------------
#========================================================================================#
otu.all.data <- read.csv(paste0(dat.dir,'20191028_ITS2_CBandCFfromfiltered_95%.csv'),
                         as.is = T, row.names = 1)

#--TBAS taxonomy of ITS2 CB and CF data
tax.data <- read.csv(paste0(dat.dir,'CultureBasedCultureFree_Assignments.csv'), as.is = T)

#--Make results table for culture-based data
results.cf.cb <- data.frame(measure = c('Jaccard', 'Jaccard.chiricahua'))
#-----------------------------------------------------------------------------------------#
# Remove rows with contaminants and mock community------------------
#-----------------------------------------------------------------------------------------#
# 
# #--remove NTCs, extraction blanks, and contamination row
# otu.all.data <- otu.all.data[!otu.all.data$range %in% c('contaminant','other','mock'),]
# 
# #--For testing within range variation between culture-based and culture-free
# Chiricahua.samples <- c("C13c","C11c","C15c","C22c","C23c","C25c","C31c","C32c","C33c",
#                         "C34c","C35c","C11","C12","C13","C14","C21","C22","C23","C32",
#                         "C33","C34")
# Chir.matrix <- otu.all.data[Chiricahua.samples,]

# #--remove culture-based samples with 1 OTU
# otu.all.data.clean <- otu.all.data[3:length(otu.all.data)]
# otu.all.data.clean <- otu.all.data.clean[rowSums(otu.all.data.clean) > 1,]
# #--Remove OTU with less than or equal to 2 occurences
# otu.all.data.clean <- otu.all.data.clean[colSums(
#   otu.all.data.clean) > 3]


#=========================================================================================#
# Overlap-------------
#=========================================================================================#

cult.based <- otu.all.data[otu.all.data$origin == 'CultureBased',]
cult.based.comm <- cult.based[5:length(cult.based)]
cult.based.comm <- cult.based.comm[,which(colSums(cult.based.comm,na.rm = T) > 1)]

cult.free <- otu.all.data[otu.all.data$origin == 'CultureFree',]
cult.free.comm <- cult.free[5:length(cult.free)]
cult.free.comm <- cult.free.comm[,which(colSums(cult.free.comm[],na.rm = T) > 8)]

shared.colnames <- colnames(cult.based.comm[colnames(cult.based.comm) %in% colnames(cult.free.comm)])
shared.otu <- otu.all.data[colnames(otu.all.data) %in% shared.colnames]

# write.csv(shared.otu, paste0(res.dir, 'ITS2_sharedCBCF.csv'),
#           row.names = T)

#=========================================================================================#
# Taxonomy-------------
#=========================================================================================#

shared.class <- tax.data[tax.data$Otu.95 %in% shared.colnames,
                         c('Otu.95','Most.common.class.level.assignment','Class.level.assignment')]
shared.class <- distinct(shared.class)

shared.class %>% 
  count(Class.level.assignment) %>%
  mutate(rel.ab = n/23) -> relab.shared

ggplot(shared.class, aes(x = Class.level.assignment)) +
  geom_bar()2
