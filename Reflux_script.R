## Reflux study analysis

# title:
# author:
# output:

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(grafify)

ref <- read.csv('Select reflux 01-08-2023.csv')
prot <- read.csv('Progenesis_CF_HV_unique_proteins.csv')
mic1 <- read.csv('GORCF_MAR_2024_Raw_Count.csv')
mic2 <- read.csv('GORCF_MAR_2024_Normalised_Count.csv')
clin <- read.csv('GOR_ClinicalEndpoint.csv')
clin1 <- read.csv('GOR-CF database 07-03-2024.csv')
clin2 <- read.csv('GORCF data 24-05.csv') %>% drop_na()
mic0 <- read.csv('GORCF sputum TEST.csv')

## EXPLORATORY ANALYSIS

##################################

# number of patients with complete set of observations
mic3 <- mic0 %>% filter(Sample == 'Sputum') %>% drop_na()
ref1 <- ref[-1,] 

r <- data.frame(PatID = ref1$PatID, dataR = 'Yes')
p <- data.frame(PatID = colnames(assay), dataP = 'Yes')
m <- data.frame(PatID = mic3$PatID, dataM = 'Yes')

View(r %>% full_join(p) %>% full_join(m) %>% drop_na())
mr <- r %>% full_join(m)
mr[is.na(mr) == TRUE] <- 'No'

# cohort demographics overview
clin2 %>% ggplot(aes(x = Sex, y = Age, fill = Sex)) + geom_boxplot() + theme_bw(base_size = 10) +
  ggtitle('Patient age distribution (18-48 years)') +
  theme(plot.title = element_text(hjust = 0.7, size = 10),
        legend.position = 'none') +
  scale_x_discrete(labels = c('Female \n n = 12', 'Male \n n = 39')) +
  scale_fill_manual(values = c('#F8766D', "#999999"))
  

ggsave('AgevsSex.tiff', dpi = 300, width = 65, height = 55, units = 'mm')

range(clin2$Age)

# Clinical data overview
library(forcats)
library(stringr)

colnames(clin)[1] <- 'PatID'

clin <- clin %>% mutate(FEV_change = X1y_priorFEV1. - BaselineFEV1pp)

clinl1 <- clin %>% dplyr::select(!FEV_change) %>% pivot_longer(cols = !PatID, names_to = 'variable', values_to = 'value') 

fev <- clinl1 %>% filter(variable %ni% c('IV.Days.12months', 'IVCourses12months')) %>% 
  ggplot(aes(x = variable, y = value)) + 
  geom_boxplot() + geom_jitter() + theme_bw(base_size = 12) + 
  ylab('% Predicted FEV') +
  theme(axis.title.x = element_blank())

iv <- clinl1 %>% filter(variable == 'IV.Days.12months') %>% 
  ggplot(aes(x = variable, y = value)) + 
  geom_boxplot() + geom_jitter() + 
  theme_bw(base_size = 12) + ylab('Days') +
  theme(axis.title.x = element_blank())
iv

ivc <- clinl1 %>% filter(variable == 'IVCourses12months') %>% 
  ggplot(aes(x = variable, y = value)) + 
  geom_boxplot() + geom_point(position = position_dodge2(width = 0.5)) + 
  theme_bw(base_size = 12) + ylab('Number') +
  theme(axis.title.x = element_blank())
ivc

g1 <- arrangeGrob(fev, iv, ivc, nrow = 1, widths = c(0.6, 0.2, 0.2),
                  top = grid::textGrob('Clinical variables', gp=grid::gpar(fontsize=10)))
plot(g1)
ggsave('ClinicalBoxplots.tiff', g1, width = 200, height = 70, unit = 'mm', dpi = 300)



fevc <- clinl1 %>% filter(variable %ni% c('IV.Days.12months', 'IVCourses12months')) %>% 
  ggplot(aes(x = variable, y = value, colour = PatID)) + 
  geom_point() + geom_line(aes(group = PatID)) + theme_bw(base_size = 12) + 
  ylab('% Predicted FEV') +
  theme(axis.title.x = element_blank(),
        legend.position = 'none')

fevc

# map microbiome data batches to clinical outcomes
library(stringr)

clin1 <- clin %>% left_join(mr %>% dplyr::select(PatID, dataM)) %>%
  left_join(mics1 %>% mutate(PatID = rownames(.),
                             dataM = 'Yes') %>%
              dplyr::select(PatID, dataM), by = 'PatID') %>%
  mutate(dataM = ifelse(dataM.x == 'Yes' & dataM.y == 'Yes', 'Batch1',
                        ifelse(dataM.x == 'No' & dataM.y == 'Yes', 'Batch2', 'None')))

clin1$dataM[is.na(clin1$dataM) == TRUE] <- 'None'

clinl <- clin1 %>% pivot_longer(cols = c(BaselineFEV1pp, IV.Days.12months), names_to = 'measure', values_to = 'value')

labs <- c('Baseline FEV1%', 'IV days 12m')

cp <- clinl %>%
  ggplot(aes(x = fct_inorder(measure), y = value, colour = dataM)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(dodge.width = 0.75) +
  theme_bw(base_size = 12) + xlab('Lung function measure') +
  scale_x_discrete(labels = str_wrap(labs, width = 15)) +
  scale_colour_manual(values = c("indianred1","grey34", 'black'), 
                    name = 'Batch')

cp

ggsave('MicrobiomevsClin1.tiff', dpi = 300, width = 100, height = 80, unit = 'mm')

# map microbiome data batches to reflux data
library(ggbeeswarm)


rl <- ref1 %>% 
  left_join(mr %>% dplyr::select(PatID, dataM)) %>%
  left_join(mics1 %>% mutate(dataM = 'Yes') %>%
              dplyr::select(PatID, dataM), by = 'PatID') %>%
  mutate(dataM = ifelse(dataM.x == 'Yes' & dataM.y == 'Yes', 'Batch1',
                        ifelse(dataM.x == 'No' & dataM.y == 'Yes', 'Batch2', 'None'))) %>%
  pivot_longer(cols = Proximal.total.events, names_to = 'measure', values_to = 'value')

rl$dataM[is.na(rl$dataM) == TRUE] <- 'None'

rl <- rl %>% mutate(data = ifelse(dataM == 'None', 'Prot', 'Prot & Mic'))

rp <- rl %>%  ggplot(aes(x = measure, y = value, colour = data)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(dodge.width = 0.75) +
  theme_bw(base_size = 12) + xlab('Reflux measure') +
  scale_x_discrete(labels = 'Proximal total events') +
  scale_colour_manual(values = c("indianred1",'black'), 
                    name = 'Batch') #+
  theme(legend.position = 'none')

rp

ggsave('MicrobiomevsRef1.tiff', dpi = 300, width = 100, height = 80, unit = 'mm')



# boxplots of clinical data 
ref2l %>% ggplot(aes(x = measure, y = value)) + geom_boxplot(outlier.shape = NA, aes(fill = category)) + 
  geom_jitter(size = 1) +
  facet_wrap(~factor(range, c('loRange', 'midRange', 'hiRange')), scales = 'free') +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = wrap_format(10))+
  ggtitle('Reflux measures') +
  scale_fill_brewer(palette = 'Set2')

# PCA of clinical data
library(factoextra)
library(ggfortify)
library(ggrepel)

rownames(clin) <- clin$PatID
clin1 <- clin %>% 
  dplyr::select(BaselineFEV1pp, IVCourses12months, IV.Days.12months) %>%
  drop_na()

pcac <- prcomp(clin1, scale = TRUE)

eigc <- get_eigenvalue(pcac)
scoresc <- pcac$x

autoplot(pcac, loadings = TRUE, loadings.label = TRUE, data = annClin, colour = 'Demeester.score') +
  geom_text_repel(aes(label = rownames(clin1)), size = 3)

annClin <- clin1 %>% mutate(PatID = rownames(clin1)) %>%
  left_join(ref3 %>% mutate(PatID = rownames(ref3)))
rownames(annClin) <- annClin$PatID

plot(annClin$Proximal.totals, annClin$IV.Days.12months)

# correlation of clinical variables
library(ggcorrplot)

corc <- cor(clin %>% dplyr::select(!PatID), method = 'pearson',
            use = 'complete.obs')

col <- colorRampPalette(c("red", "white", "blue"))(20)

ggcorrplot(corc, tl.cex = 10, hc.order = TRUE, legend.title = 'r', lab = TRUE,
           lab_size = 2.5) +
  ggtitle('Pearsons correlation of clinical variables') +
  theme(plot.title = element_text(size = 10))

ggsave('CorrClin.tiff', dpi = 300, width = 105, height = 105, unit = 'mm')

# map data available in microbiome dataset to key reflux variables
ref1 %>% full_join(mr %>% dplyr::select(PatID, dataM)) %>%
  dplyr::select(PatID, dataM, Total.events, Proximal.total.events) %>%
  pivot_longer(cols = c(Total.events, Proximal.total.events), names_to = 'variable', values_to = 'value') %>%  
  ggplot(aes(x = variable, y = value, colour = dataM)) + geom_violin() + geom_point(position = position_jitterdodge()) +
  theme_bw(base_size = 10) + scale_colour_discrete(name = 'Microbiome data', labels = c('Yes', 'No')) +
  ylab('count')

ggsave('RefMicrobiome_violin.tiff', width = 135, height = 83, unit = 'mm', dpi = 300)

###############################  

library(stringr)

# Exploratory analysis of reflux data
ref2 <- ref1 %>% drop_na()

# visualise distribution of reflux measures

ref2l <- ref2 %>% pivot_longer(cols = !PatID, names_to = 'measure', values_to = 'value')
ref2sum <-  ref2l %>% group_by(measure) %>% summarise(max = max(value),
                                  median = median(value)) %>%
  mutate(range = ifelse(max < 50, 'loRange', ifelse(
    max > 100, 'hiRange', 'midRange'
  )))

ref2l <- ref2l %>% left_join(ref2sum)
ref2l$range <- as.factor(ref2l$range)
ref2l <- ref2l %>% mutate(category = ifelse(str_detect(measure, 'events') == TRUE, 'event', ifelse(
  str_detect(measure, 'exposure') == TRUE, 'exposure', ifelse(
    str_detect(measure, 'episode') == TRUE, 'episode', 'exposure'
  )
)))

ref2l$measure <- gsub('\\.', ' ', ref2l$measure)
                             
library(scales)  

ref2l %>% ggplot(aes(x = measure, y = value)) + geom_boxplot(outlier.shape = NA, aes(fill = category)) + 
  geom_jitter(size = 1) +
  facet_wrap(~factor(range, c('loRange', 'midRange', 'hiRange')), scales = 'free') +
  theme_bw(base_size = 9) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = 'none',
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = wrap_format(10))+
  ggtitle('Reflux measures') +
  scale_fill_brewer(palette = 'Set2')

ggsave('RefluxBoxplots.tiff', dpi = 300, width = 220, height = 65, unit = 'mm')

# Principal Component Analysis of reflux measures
library(factoextra)

rownames(ref2) <- ref2$PatID
ref3 <- ref2[,-1]

colnames(ref3) <- gsub('\\.episode', '', 
                       gsub('\\.exposure', '', 
                       gsub('\\.events', '', colnames(ref3))))

pca <- prcomp(as.matrix(ref3), scale = TRUE)

scores <- as.data.frame(pca$x)                        
loadings <- as.data.frame(pca$rotation)

eigen <- get_eigenvalue(pca)

barplot(eigen$variance.percent[1:10])

scores %>% ggplot(aes(y = PC1, x = PC2)) + geom_point() + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash')

library(ggfortify)
library(ggrepel)

p <- autoplot(pca, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 4, 
         loadings.label.repel = TRUE) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  ggtitle('PCA biplot of reflux measures') +
  theme(plot.title = element_text(hjust =0.5)) +
  geom_text_repel(aes(label = rownames(ref3)), size = 3)

  
# keep only selected loadings    
p$layers
  
remove <- c('Gass', 'Upright.Longest', 'Supines')

p$layers[[2]]$data <- p$layers[[2]]$data[!(rownames(p$layers[[2]]$data) %in% remove),] 
p$layers[[3]]$data <- p$layers[[3]]$data[!(rownames(p$layers[[3]]$data) %in% remove),]
  
p

# change loadings colors  
ref2l$measure <- gsub('\\.episode', '', 
                      gsub('\\.exposure', '', 
                           gsub('\\.events', '', ref2l$measure)))

library(RColorBrewer)

set2 <- brewer.pal(8, 'Set2')

loads <- data.frame(measure = rownames(loadings)) %>% 
  left_join(ref2l %>% dplyr::select(measure, category) %>% unique()) %>%
  mutate(hue = ifelse(category == 'event', set2[2], ifelse(
    category == 'episode', set2[1], ifelse(
      category == 'exposure', set2[3], set2[3]
    )
  )))

p <- autoplot(pca, loadings = TRUE, loadings.label = TRUE, loadings.label.size = 4, 
              loadings.label.repel = TRUE, loadings.colour = loads$hue,
              loadings.label.colour = loads$hue) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  ggtitle('PCA biplot of reflux measures') +
  theme(plot.title = element_text(hjust =0.5))
#+
geom_text_repel(aes(label = rownames(ref3)), size = 3)

p

ggsave('PCAReflux.tiff', dpi = 300, width = 130, height = 90, unit = 'mm')  
  
# loadings plots
library(gridExtra)
library(forcats)

p1 <- loadings %>% arrange(PC1) %>% ggplot(aes(x = fct_inorder(rownames(.)), y = PC1)) + 
  geom_bar(stat = 'identity') + coord_flip() + theme_bw(base_size = 12) + 
  theme(axis.title.y = element_blank())

p2 <- loadings %>% arrange(PC2) %>% ggplot(aes(x = fct_inorder(rownames(.)), y = PC2)) + 
  geom_bar(stat = 'identity') + coord_flip() + theme_bw(base_size = 12) +
  theme(axis.title.y = element_blank()) + 
  xlab('PC2 loadings')

g <- arrangeGrob(p1, p2, nrow = 1, ncol =2, top = 'Principal Component Loadings')
plot(g)

ggsave('LoadingsRefPCA.tiff', g, width = 140, height = 90, unit = 'mm', dpi = 300)

# correlation matrix of reflux measures
library(ggcorrplot)

cormat <- cor(ref3, method = 'pearson')
col <- colorRampPalette(c("red", "white", "blue"))(20)
ggcorrplot(cormat, tl.cex = 6, hc.order = TRUE, legend.title = 'r', lab = TRUE,
           lab_size = 2)

# create patient groups based on reflux measure separation
acLength <- c('CF036', 'CF029', 'CF023', 'CF050', 'CF022', 'CF011', 'CF008', 'CF007', 'CF012')
proxTot <- c('CF001', 'CF009', 'CF015', 'CF003', 'CF035', 'CF020', 'CF010')

refGroups <- data.frame(PatID = ref1$PatID) %>%
  mutate(refGroup = ifelse(PatID %in% acLength, 'acLength', ifelse(
    PatID %in% proxTot, 'proxTot', ifelse(
      PatID %in% c('CF021', 'CF027'), 'NA', 'other'))))

# plot proximal total events across patients
ref1 %>% ggplot(aes(x = PatID, y = sort(Total.events))) + 
  geom_bar(stat = 'identity') + theme(axis.text.x = element_text(angle = -90))



################################



# Exploratory analysis of proteomics data
library(stringr)

colnames(prot) <- prot[2,]
prot <- prot[-c(1:2),]

assay <- prot %>% dplyr::select(contains('CF'))
rownames(assay) <- prot$Accession  

assay <- assay %>% mutate_if(is.character, as.numeric)

# intensity value distribution and variance
assayl <- assay %>% mutate(prot = rownames(assay)) %>% 
  pivot_longer(cols = !prot, names_to = 'PatID', values_to = 'intensity')

assayl %>% ggplot(aes(x = intensity, colour = PatID)) + geom_density()

assayl %>% group_by(prot) %>% summarise(mean = mean(intensity), sd = sd(intensity)) %>%
  ggplot(aes(y = sd, x = rank(mean))) + geom_point()

# variance stabilisation
library(vsn)
library(limma)
library(gridExtra)
library(geomtextpath)

vsn <- vsn::vsnMatrix(as.matrix(assay))
vsnout <- vsn@hx
vsnout <- as.data.frame(vsnout) %>% mutate(prot = rownames(vsnout)) %>% relocate(prot)


vsnl <- vsnout %>% pivot_longer(cols = !prot, names_to = 'PatID', values_to = 'intensity')

vsnl %>% ggplot(aes(x = intensity, colour = PatID)) + geom_textdensity(aes(label = PatID)) +
  theme_bw(base_size = 10) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.75, size = 10)) +
  ggtitle('Density plot of transformed protein data')

ggsave('DensityProteomicsVSN.tiff', width = 70, height = 60, unit = 'mm', dpi =300)

 
vsnl %>% group_by(prot) %>% summarise(mean = mean(intensity), sd = sd(intensity)) %>%
  ggplot(aes(y = sd, x = rank(mean))) + geom_point() 

# Principal Component Analysis
library(factoextra)
library(ggfortify)

assayT <- t(vsnout[,-1]) 

write.csv(assayT, 'ProteinVSNTransformed.csv')

assayT <- read.csv('ProteinVSNTransformed.csv')
colnames(assayT)[1] <- 'PatID'

pca1 <- prcomp(assayT, scale = TRUE)

scores1 <- as.data.frame(pca1$x)                        
loadings1 <- as.data.frame(pca1$rotation)

eigen1 <- get_eigenvalue(pca1)

barplot(eigen1$variance.percent[1:10])

an <- data.frame(PatID = rownames(assayT)) %>% 
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))

rownames(an) <- an$PatID

scores1 <- scores1 %>% mutate(PatID = rownames(.)) %>% left_join(an)

pc23 <- scores1 %>% ggplot(aes(y = PC3, x = PC2)) + 
  geom_point(aes(fill = ptile_var), size = 1, pch = 21) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  scale_fill_gradient2(low = '#0000FFFF',
                       mid = '#009E61FF',
                       high = '#FFFF00FF',
                       midpoint = 0.5,
                       labels = round(quants,1)) +
  theme_bw(base_size = 9) +
  ylab('PC3 (10.5%)') + labs(fill = 'Proximal \n reflux \n events') +
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.key.height = unit(5, 'mm'),
        legend.key.width = unit(3, 'mm'),
        axis.title.x = element_blank())

pc23
  
pc12 <- scores1 %>% ggplot(aes(y = PC1, x = PC2)) + 
  geom_point(aes(fill = ptile_var), size = 1, pch = 21) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  scale_fill_gradient2(low = '#0000FFFF',
                       mid = '#009E61FF',
                       high = '#FFFF00FF',
                       midpoint = 0.5,
                       labels = round(quants,1)) +
  theme_bw(base_size = 9) +
  ylab('PC1 (21.2%)') + labs(fill = 'Proximal \n reflux \n events') +
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.key.height = unit(5, 'mm'),
        legend.key.width = unit(3, 'mm'),
        axis.title.x = element_blank())


pc12  
library(ggpubr)

protPC <- ggarrange(pc12, pc23, common.legend = TRUE, legend = 'right', nrow = 2)
protPC
annotate_figure(protPC, top = text_grob("PCA plots of sputum proteome data", 
                                      color = "black", face = "bold", size = 10),
                bottom = text_grob('PC2 (10.9%)', size = 9))

ggsave('Proteome_PCA.tiff', dpi = 300, width = 75, height = 90, unit = 'mm')

# patients without proteomic profiles
setdiff(refGroups$PatID, rownames(assayT))

# create quantile colour scale 
testcol <- c('#0000FFFF', '#009E61FF','#FFFF00FF')
colvalquant <- base::seq(from = 0, to = 1, length.out = 5)
quants <- quantile(an$Proximal.total.events, colvalquant)
an$ptile_var <- ecdf(an$Proximal.total.events)(an$Proximal.total.events)

a1 <- autoplot(pca1, data = an, x = 1, y = 2, colour = 'ptile_var', size = 2) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() #+
  #scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis')
  #+
  #geom_text(aes(x = 0.2, y = 0.73), label = 'CF055', size = 4)

a1

a2 <- autoplot(pca1, data = an,
         colour = 'Proximal.total.events', x = 3, y = 2, size = 2) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis') +
  #scale_color_brewer(palette = 'Set2') +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) +
  ggtitle('PCA biplot of sputum protein data') #+
  #geom_text(aes(x = 0.40, y = 0.73), label = 'CF055', size = 4)
  #geom_text_repel(label = rownames(pca1$x))

a2


a <- arrangeGrob(a1, a2, ncol = 2, nrow = 1, top = 'PCA biplots of sputum protein data')

ggsave('PCAProteins.tiff', a, dpi = 300, width = 220, height = 90, unit = 'mm')

# create PatID groups based on proteomics PCA separation
p3a <- c('CF007', 'CF009', 'CF001', 'CF008', 'CF035', 'CF017', 'CF003', 'CF015', 'CF012', 'CF018', 'CF029', 'CF023')
p3out <- c('CF052', 'CF045', 'CF055')

annot <- data.frame(PatID = rownames(scores1)) %>%
  mutate(pc3 = ifelse(PatID %in% p3a, 'p3a', ifelse(
    PatID %in% p3out, 'p3out', 'p3b')))

p1a <- c('CF010', 'CF022', 'CF006', 'CF027', 'CF034', 'CF056', 'CF005', 'CF053', 'CF036', 'CF033')
p1b <- c('CF015', 'CF024', 'CF032', 'CF014', 'CF021')
p1c <- c('CF052', 'CF045', 'CF055')

refGroups <- refGroups %>% mutate(
  pc3 = ifelse(PatID %in% p3a, 'p3a', ifelse(
  PatID %in% p3out, 'p3out', 'p3b')))#,
  pc1 = ifelse(PatID %in% p1a, 'p1a', ifelse(
    PatID %in% p1b, 'p1b', ifelse(
      PatID %in% p1c, 'p1c', 'p1d'
    )
  ))

rownames(refGroups) <- refGroups$PatID
refGroups1 <- refGroups %>% filter(PatID %ni% c("CF011", "CF040", "CF041")) %>%
  dplyr::select(!PatID)

# map proteomic groups to reflux measure PCA
dev.new()

autoplot(pca, #loadings = TRUE, loadings.label = TRUE, loadings.label.size = 3, loadings.label.repel = TRUE,
         data = refGroups %>% filter(PatID %in% rownames(ref3)), colour = 'pc3') + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  geom_text_repel(aes(label = rownames(ref3)), size = 3)

# patients separated along proteomics PC3 at difference risk of experiencing extreme reflux events? 

# plot total events across patients
refGroups %>% left_join(ref1 %>% select(PatID, Total.events)) %>% 
  ggplot(aes(x = PatID, y = sort(Total.events), fill = pc3)) + 
  geom_bar(stat = 'identity') + theme(axis.text.x = element_text(angle = -90))


# proteins contributing most to PC3
p3l <- data.frame(pc3 = loadings1$PC3, 
                  Accession = rownames(loadings1)) %>%
  arrange(pc3)

p3lm <- p3l[c(1:10, 430:440), 1:2]
p3lm <- p3lm %>% left_join(prot %>% select(Accession, Description))

# viusualise specific loadings on a biplot
p0 <- autoplot(pca1, data = assayT,
         x = 2, y = 3, loadings = TRUE, loadings.label = TRUE) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  geom_text_repel(aes(label = rownames(assayT)), size = 3, vjust = -0.5)

p0$layers

p0$layers[[2]]$data <- p0$layers[[2]]$data["P00450",]
p0$layers[[3]]$data<-p0$layers[[3]]$data["P00450",]

p0


# difference according to reflux groups
as.data.frame(assayT) %>% mutate(PatID = rownames(assayT)) %>%
  left_join(ref4 %>% dplyr::select(PatID, RefGroup)) %>%
  ggplot(aes(x = RefGroup, y = A0A0C4DH38)) + geom_boxplot()



#############################



# correlation of reflux and protein PCs 
scoresCor <- scores %>% mutate(PatID = rownames(scores))
scores1Cor <- scores1 %>% mutate(PatID = rownames(scores1))

scoresCor2 <- scoresCor %>% dplyr::select(PatID, PC1, PC2) %>%
  full_join(scores1Cor %>% dplyr::select(PatID, PC1, PC2, PC3), by = 'PatID')

setdiff(scoresCor$PatID, scores1Cor$PatID)
setdiff(scores1Cor$PatID, scoresCor$PatID)
scoresCor2 <- scoresCor2 %>% filter(PatID %ni% c("CF011", "CF040", "CF041", 'CF021', 'CF027'))
colnames(scoresCor2) <- c('PatID', 'PC1.R', 'PC2.R', 'PC1.P', 'PC2.P', 'PC3.P')

library(ggcorrplot)

cormatPC <- cor(scoresCor2 %>% dplyr::select(!PatID), method = 'pearson')
col <- colorRampPalette(c("red", "white", "blue"))(20)
ggcorrplot(cormatPC, tl.cex = 10, hc.order = TRUE, legend.title = 'r', lab = TRUE,
           lab_size = 3) +
  ggtitle('Pearsons correlation of principal component scores') +
  theme(plot.title = element_text(hjust = 0.5, size = 10))

ggsave('CorrPCA.tiff', dpi = 300, width = 90, height = 90, unit = 'mm')

# plot PCA groupings between datasets
annot1 <- data.frame(PatID = rownames(scores)) %>%
  mutate(pc3 = ifelse(PatID %in% p3a, 'p3a', ifelse(
    PatID %in% p3out, 'p3out', 'p3b')))

p <- autoplot(pca, data = annot1, colour = 'pc3') + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  ggtitle('PCA biplot of reflux measures') +
  theme(plot.title = element_text(hjust =0.5),
        legend.position = 'none') +
  #geom_text_repel(label =  rownames(pca$x)) +
  scale_colour_brewer(palette = 'Set2')

p

ap <- arrangeGrob(a2, p, ncol = 2, nrow = 1)

ggsave('PCAProtRef.tiff', ap, dpi = 300, width = 220, height = 90, unit = 'mm')



##################################

# Microbiome data

# Sputum relative abundance data
colnames(mic0)[1] <- 'PatID'
mic <- drop_na(mic0) %>% filter(Cohort == 'CF')

# diversity measures and descriptions
div <- mic[,c(1:3, 229:240)]

# relative abundances of bacterial geni
mic <- mic[,c(1:228)]

# remove taxa with 0 observations
table(colSums(mic1[4:228]) == 0)

micf <- as.data.frame(t(mic[4:228])) %>% filter(rowSums(.) > 0) %>%
  t() %>% as.data.frame() 

rownames(micf) <- mic$PatID 

# relative abundances of bacterial phyla (higher tax level?)
mic2 <- mic[,c(1:3, 241:ncol(mic))]

#
#
#

# Sputum and mouth swirl count data
# filter for patients with reflux data
ref1 <- ref[-1,]

mic1 <- mic1 %>% filter(PatID %in% ref1$PatID)
mic2 <- mic2 %>% filter(PatID %in% ref1$PatID)

# samples without microbiome data according to reflux measures
sampleInfo %>% full_join(mic %>% dplyr::select(PatID, Sample)) %>%
  arrange(Total.events) %>%
  ggplot(aes(x = fct_inorder(PatID), y = Total.events, fill = cca)) + geom_bar(stat = 'identity')


# number of patients and samples
n_distinct(mic2$PatID)
table(mic1$Sample_Type)
View(table(mic1 %>% group_by(PatID, Sample_Type) %>% summarise(n = n())))

# filter sputum samples only
mics1 <- mic1 %>% filter(Sample_Type == 'Sputum')
mics2 <- mic2 %>% filter(Sample_Type == 'Sputum')

# remove taxa with 0 observations
mics1 <- mics1 %>% dplyr::select(!c('OTU_ID', 'Sample_Type_A', 'Sample_Type_B', 'Cohort', 'SampleCohort'))
mics2 <- mics2 %>% dplyr::select(!c('OTU_ID', 'Sample_Type_A', 'Sample_Type_B', 'Cohort', 'SampleCohort'))

table(colSums(mics1[3:ncol(mics1)]) == 0)

mics1 <- as.data.frame(t(mics1[3:ncol(mics1)])) %>% filter(rowSums(.) > 0) %>%
  t() %>% as.data.frame() %>% cbind(mics1$PatID) %>% 
  relocate('mics1$PatID')

colnames(mics1)[1] <- 'PatID'

mics2 <- as.data.frame(t(mics2[3:ncol(mics2)])) %>% filter(rowSums(.) > 0) %>%
  t() %>% as.data.frame() %>% cbind(mics2$PatID) %>% 
  relocate('mics2$PatID')

colnames(mics2)[1] <- 'PatID'

ncol(mics1)
rownames(mics1) <- mics1$PatID
mics1 <- mics1[,-1] # unfiltered non-normalised

rownames(mics2) <- mics2$PatID
mics2 <- mics2[,-1] 

# OPTIONAL - filtering OTUs
# filter out OTUs with count exceeding 0.01% of count sum across all libraries
low.count.removal <- function(
    data, # OTU count df of size n (sample) x p (OTU)
    percent = 0.01 # cutoff chosen
) 
{
  keep.otu = which(colSums(data)*100/(sum(colSums(data))) > percent)
  data.filter = data[,keep.otu]
  return(list(data.filter = data.filter, keep.otu = keep.otu))
}

mic1f <- low.count.removal(mics1, percent = 0.01)
micf <- mic1f$data.filter # filtered % all counts


mic2f <- low.count.removal(mics2, percent = 0.01)
micf2 <- mic2f$data.filter

# filter out OTUs with count exceeding 0.1% of count sum in at least one sample
low.count.id <- lapply(1:nrow(mics2), function(x){
row <- mics2[x,]
row <- (row/rowSums(row))*100
row <- as.data.frame(row > 0.1)})

low.count.df <- bind_rows(low.count.id)
low.count.df <- low.count.df %>% mutate_if(is.logical, as.numeric)

keep.otu <- which(colSums(low.count.df) > 0)
micf2 <- mics2[, keep.otu] # filtered % sample counts


# compare library sizes across samples ordered by increasing Proximal reflux events
# compare library size by batch
library(forcats)


dev.new()
mics1 %>% mutate(libSize = rowSums(.),
                 PatID = rownames(.)) %>% 
  left_join(ref1 %>% 
              dplyr::select(PatID, Proximal.total.events)) %>%
  arrange(Proximal.total.events) %>%
  left_join(mic %>% dplyr::select(PatID, Cohort)) %>%
  ggplot(aes(x = fct_inorder(as.factor(PatID)), y = libSize, fill = Cohort)) +
  geom_col() +
  scale_fill_manual(values = c("indianred1","gray"), 
                    name = 'Batch', labels = c('Batch1', 'Batch2')) +
  theme_bw(base_size = 12) +
  xlab('Patient ID') +
  theme(axis.text.x = element_text(angle = 45, size = 6),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle('Library sizes in 16S dataset')

ggsave('LibSize.tiff', width = 150, height = 90, units = 'mm', dpi = 300)  

table(rowSums(mics1)) # difference between smallest and largest library x14


# evaluate distribution of 0s within OTUs following different filtering
table(micf[micf == 0]) # 38% 0s
table(micf1[micf1 == 0]) # 45% 0s
table(mics1[mics1 == 0]) # 71% 0s

table(micf2[micf2 == 0]) # 
table(micf12[micf12 == 0]) #  59% 0s
table(mics2[mics2 == 0]) #  76% 0s

par(mfrow = c(4, 4))
for (i in 1: 58) {
  hist(micf[,i], main = colnames(micf)[i], cex.main = 0.7, breaks = 10)}

#

# Cumulative Sum Scaling normalisation of library size
library(metagenomeSeq)

data.metagenomeSeq = newMRexperiment(t(micf), # using filtered data
                                     featureData = NULL, libSize = NULL, normFactors = NULL) 
p = cumNormStat(data.metagenomeSeq) # percentile selection - default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p = p)

micf.CSS = t(MRcounts(data.cumnorm, norm = TRUE, log = TRUE))


# ensure the data are in a correct format: number of samples in rows
dim(micf.CSS)  

# effect of log transformation on variance
CSScountl <- as.data.frame(micf1a.CSS) %>% mutate(PatID = rownames(.)) %>%
  pivot_longer(cols = !PatID, names_to = 'otu', values_to = 'count')

CSSsum <- CSScountl %>% group_by(otu) %>% summarise(mean_count = mean(count))

CSScountl %>% left_join(CSSsum) %>% arrange(mean_count) %>%
  ggplot(aes(x = fct_inorder(as.factor(mean_count)), y = count)) + geom_point()

#

CSScountLogl <- as.data.frame(micf1.CSS) %>% mutate(PatID = rownames(.)) %>%
  pivot_longer(cols = !PatID, names_to = 'otu', values_to = 'count')

CSSLogsum <- CSScountLogl %>% group_by(otu) %>% summarise(mean_count = mean(count))

CSScountLogl %>% left_join(CSSLogsum) %>% arrange(mean_count) %>%
  ggplot(aes(x = fct_inorder(as.factor(mean_count)), y = count)) + geom_point()


#

# Centered Log Ratio Transformation with offset 1
library(mixOmics)

micst <- logratio.transfo(mics1, logratio = 'CLR', offset = 1)
micstf <- logratio.transfo(micf, logratio = 'CLR', offset = 1)
micstf1 <- logratio.transfo(micf1, logratio = 'CLR', offset = 1)

class(micst) <- 'matrix'
class(micstf) <- 'matrix'
class(micstf1) <- 'matrix'

micstf2 <- logratio.transfo(micf2, logratio = 'CLR', offset = 1) # rarefied dataset
micstf12 <- logratio.transfo(micf12, logratio = 'CLR', offset = 1)

class(micstf2) <- 'matrix'
class(micstf12) <- 'matrix'

# check CLR
mics1t <- as.data.frame(t(mics1 + 1))
gm <- exp(mean(log(mics1t$CF003))) #geom mean C003
log(mics1t$CF003/gm)

log(mics1t$CF003)

#

# assess variance of M (mean of log-transformed count data within observation)
mics1o <- mics1 + 1
mics1ol <- log(mics1o)

micfo <- micf + 1
micfol <- log(micfo)

micf1o <- micf1 + 1
micf1ol <- log(micf1o)

micf2o <- micf2 + 1
micf2ol <- log(micfo)

micf12o <- micf12 + 1
micf12ol <- log(micf12o)


library(matrixStats)
rowMeans(mics1ol)
rowVars(as.matrix(mics1ol))

dataSum <- function(dataLog){
  dataLong <- dataLog %>% mutate(PatID = rownames(dataLog)) %>%
    pivot_longer(cols = !PatID, names_to = 'OTU', values_to = 'value')
  
  dataM <- dataLong %>% group_by(PatID) %>% 
    summarise(M = mean(value),
              sd = sd(value)) %>%
    left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))
}

micfsum <- dataSum(micf1ol)

var(micfsum$M)

m2 <- micfsum %>% arrange(Proximal.total.events) %>%
  ggplot(aes(x = fct_inorder(as.factor(PatID)), y = M, colour = Proximal.total.events)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = M-sd, ymax = M+sd)) + 
  theme_bw() +
  xlab('PatID') +
  ggtitle('Total count filtered')

m2

pdf('MvsFilteringMicrobiome.pdf')
grid.arrange(m1, m2, m3, nrow = 3, ncol = 1)
dev.off()

#

# alternative robust CLR 
# does not require pseudocount addition and keeps 0s as 0s
library(vegan)
micstr <- decostand(mics1, method = 'rclr')
micstrf <- decostand(micf, method = 'rclr')
micstrf1 <- decostand(micf1, method = 'rclr')


#

# PCA of sputum microbiome data
library(factoextra)
library(ggfortify)

micst # unfiltered clr-transformed
micstf # filtered as % of total count, clr-transformed
micstf1 # filtered as % of sample counts, clr-transformed

micstr # unfiltered rlclr-transformed
micstrf # filtered as % of total count, clr-transformed
micstrf1 # filtered as % of sample counts, clr-transformed

pcam <- pca(micstf1)
plotIndiv(pcam)

annot <- data.frame(PatID = rownames(micf)) %>% left_join(ref1 %>% 
    dplyr::select(PatID, Proximal.total.events, Total.events)) %>%
  left_join(clin %>% dplyr::select(PatID, BaselineFEV1pp, IV.Days.12months)) %>%
  left_join(mr %>% dplyr::select(PatID, dataM)) %>%
  mutate(dataM = ifelse(dataM == 'Yes', 'Batch1', 'Batch2')) %>%
  left_join(micf %>% mutate(PatID = rownames(.),
                             libSize = rowSums(.)) %>%
              dplyr::select(PatID, libSize))



rownames(annot) <- annot$PatID


pcam1 <- prcomp(micstf, center = TRUE, scale = TRUE)

get_eigenvalue(pcam1)$variance.percent[1:10]
barplot(get_eigenvalue(pcam1)$variance.percent[1:10])

p2 <- autoplot(pcam1, data = annot, colour = 'Proximal.total.events', size = 2, x = 1, y = 2) + 
  theme_bw() +
  ggtitle('PCA plot of CLR transformed data') +
  #geom_text_repel(aes(label = rownames(micf1)), size = ) +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  #scale_colour_manual(values = c("indianred1","grey34"), 
                  #  name = 'Batch', labels = c('Batch1', 'Batch2')) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis')


p2
ggsave('PCA_Proximal_TFCLR.tiff', dpi = 300, width = 150, height = 100, unit = 'mm')


micpca <- arrangeGrob(p1, p2, p3, p4, nrow = 2, ncol = 2, top = 'PCA non scaled')
dev.new()
plot(micpca)
ggsave('PCAMicrobiomeCSSScaled.tiff', micpca, width = 400, height = 200, unit = 'mm', dpi = 300)


# batch effect & library size effect (marginal mean of log count data per observation)
annot <- annot %>% full_join(mic %>% dplyr::select(PatID, Cohort)) %>%
  left_join(micsum %>% dplyr::select(PatID, M)) %>%
  left_join(micfsum %>% dplyr::select(PatID, M), by = 'PatID') %>%
  left_join(micf1sum %>% dplyr::select(PatID, M), by = 'PatID') %>%
  rename(M = M.x,
         Mf = M.y,
         Mf1 = M) %>%
  mutate(sum = rowSums(mics1),
         sumf = rowSums(micf),
         sumf1 = rowSums(micf1))

colnames(annot)[6] <- 'libSize'

p2 <- autoplot(pcam1, data = annot, colour = 'sumf', size = 3, x = 1, y = 2) + theme_bw() +
  ggtitle('CSS total count filtered') +
  geom_text_repel(aes(label = rownames(micst)), size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed')


p2

tiff('MvsMicrobiomePCA.tiff', unit = 'mm', res = 300, width = 300, height = 180)
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)
dev.off()

annot %>% ggplot(aes(x = Cohort, y = Proximal.total.events, colour = Cohort)) + geom_violin() + geom_jitter()


# correspondence analyis on untransformed compositional data
library(vegan)

ca <- cca(mics1)
plot(ca, scaling = 'symmetric', display = 'sites', type = 'text', hill = TRUE)

# Principal Coordinate Analysis on count data
library(ape)

micfRA <- decostand(micf, method = 'total')
distRA <- vegdist(micfRA, method = 'bray')

pcoa <- cmdscale(distRA, eig = TRUE)

pcoa <- pcoa(distRA)
biplot(pcoa)

vectors <- as.data.frame(pcoa[['vectors']])
vectors %>% ggplot(aes(x = Axis.1, y = Axis.2, colour = annot$Proximal.total.events)) + 
  geom_point(size = 2) +
  #geom_text_repel(aes(label = rownames(vectors)), size = 2) +
  theme_bw() +
  ggtitle('Principal Coordinate Analysis of relative abundance data') +
  theme(plot.title = element_text(hjust = 0.5)) +
  #scale_colour_manual(values = c("indianred1","grey34"), 
                      #name = 'Batch')
  scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis')
  
ggsave('PCoA_RA_rfelux.tiff',  dpi = 300, width = 150, height = 100, unit = 'mm')

#

# pCoA with vegan
library(plotrix)

colnames(micfRA) <- otuNames2$otu
vpcoa <- capscale(micfRA ~ 1, distance = 'bray')
dev.new()
ordiplot(vpcoa, type = 'text')

summary(vpcoa)

annot <- data.frame(PatID = rownames(micf)) %>% left_join(ref1 %>% dplyr::select(Proximal.total.events, PatID))

colval <- color.scale(annot$Proximal.total.events, c(0,0,1), c(0,1,1), c(1,0,0))
annot <- cbind(annot, colval)

testcol <- c('#0000FFFF', '#002BD4FF', '#0052ADFF', '#009E61FF','#FFFF00FF')
labels = c('0', '5.2', '10.2', '19.7', '66.3')

dev.new()

tiff('S_microbiome_PCoA.tiff', res = 300, width = 130, height = 100, unit = 'mm')
par(mar = c(4, 4, 4, 4))
plot(vpcoa, type = 'n', main = 'PCoA of sputum microbiome data', cex.main = 1, cex.lab = 1,
    yaxt="n", xaxt = 'n', xlab = 'MDS1 (42.7%)', ylab = 'MDS2 (14.1%)')
points(vpcoa, display = "sites", cex = 0.7, pch = 21, 
       bg = colval)
text(vpcoa, display = "sp", cex = 0.4, col = "darkred", select = colSums(micfRA > 0.4) > 0)
color.legend(-3.2,-2.5,-2.9,-1.5, legend = labels, rect.col = testcol, cex = 0.4, gradient = 'y',
             align = 'lt')
axis(1 ,cex.axis=1)
axis(2 ,cex.axis=1)
text(x = -2.9, y = -1.4, 'Proximal_Total_Event_No', cex = 0.4, col = 'black')

dev.off()


#

#

NMS <- metaMDS(distRaw,
          distance = "bray",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

goodness(NMS)
stressplot(NMS)

dev.new()
plot(NMS, 'sites')
orditorp(NMS, 'sites')

NMSscores <- as.data.frame(scores(NMS))

NMSscores %>% mutate(PatID = rownames(.)) %>%
  left_join(annot, by = 'PatID') %>%
  ggplot(aes(x = NMDS1, y = NMDS2, colour = dataM)) + 
  geom_point() +
  geom_text_repel(aes(label = rownames(NMSscores))) +
  ggtitle('CSS') #+
  scale_colour_continuous(type = 'viridis') 


# diagnostics of library size effect
# center the clr transformed data matrix
matCen <- function(dataT){
  s <- lapply(1:ncol(data), function(x){
  data[,x] <- data[,x]-mean(data[,x])})   
  sdf <- bind_cols(s)}


data <- micf1
dataT <- micstf1
dataL <- micf1ol

micstfC <- matCen(dataT)

colnames(micstfC) <- colnames(dataT)
rownames(micstfC) <- rownames(dataT)

# calculate correlation between centered matrix columns and row means of log transformed data
scor <- lapply(colnames(micstfC), function(x){
  col <- micstfC %>% dplyr::select(x)
  col <- cor(col, rowMeans(dataL))
  })

scordf <- bind_cols(scor)
colnames(scordf) <- colnames(micstfC)
scordft <- as.data.frame(t(scordf))
colnames(scordft)[1] <- 'cor'
scordft1 <- cbind(scordft, log(colMeans(data)))
colnames(scordft1)[2] <- 'OTU_logMean'

scordft1 %>% ggplot(aes(x = OTU_logMean, y = cor)) + geom_point() +
  ggtitle('Sample filtered rarefied CLR')


#
#
#

# plotting microbiome composition across patients at the phylum level
ra <- function(data) {
  dataRA <- lapply(1:ncol(data), function(x){
  data[,x] <- (data[,x]/rowSums(data))*100})
  dataRA <- bind_cols(dataRA)}


micf1RA <- ra(micf1)
micf1CSSRA <- ra(micf1a.CSS)

colnames(micf1CSSRA) <- colnames(micf1.CSS)
rownames(micf1CSSRA) <- rownames(micf1.CSS)

library(vegan)

micfRA <- decostand(micf, method = 'total')


micfRAl <- micfRA %>% mutate(PatID = rownames(.)) %>%
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events)) %>%
  pivot_longer(cols = !c(PatID, Proximal.total.events), names_to = 'OTU', values_to = 'RA')

micfRAl <- micfRAl %>% separate(OTU, c('k', 'T1', 'Phylum', 'Class', 'Order', 'Family', 'Genus'), sep = '__') %>%
  mutate(Family = gsub("\\..*","", Family),
          Order = gsub("\\..*","", Order),
          Class = gsub("\\..*","", Class),
         Phylum = gsub("\\..*","", Phylum))

micfRAl <- micfRAl %>% mutate(
  Family = ifelse(Family == '', NA, Family),
  Genus = ifelse(Genus == '', NA, Genus)
)


# order x axis according to reflux variable
Xorder <- micfRAl %>% arrange(Proximal.total.events) %>% dplyr::select(PatID)
Xorder <- unique(Xorder)

# according to PLS Comp1
XorderS <- RMscoresX %>% mutate(PatID = rownames(.)) %>% arrange(comp1)

# clinical variable
Xorder0 <- micf1RAl %>% left_join(clin %>% dplyr::select(PatID, BaselineFEV1pp, IV.Days.12months)) %>%
  arrange(desc(BaselineFEV1pp)) %>% dplyr::select(PatID) %>% unique()

Xorder01 <- micf1RAl %>% left_join(clin %>% dplyr::select(PatID, BaselineFEV1pp, IV.Days.12months)) %>%
  arrange(IV.Days.12months) %>% dplyr::select(PatID) %>% unique()

# order x axis according to sPLS scores
Xorder1 <- RMSscores %>% arrange(comp1) %>% mutate(PatID = rownames(.)) %>%
  dplyr::select(PatID)

# order x axis according to the library size
Xorder2 <- mics1 %>% mutate(PatID = rownames(mics1),
                            libSize = rowSums(mics1)) %>%
  arrange(libSize) %>%
  left_join(annot %>% dplyr::select(PatID, dataM)) %>%
  dplyr::select(PatID, dataM)

B1IDs <- Xorder2 %>% filter(dataM == 'Batch1')
B2IDs <- Xorder2 %>% filter(dataM == 'Batch2')

# create distinct color palette
library(Polychrome)
library(RColorBrewer)

p53 <- createPalette(53,  c("#ff0000", "#00ff00", "#0000ff"))
p53 <- as.vector(t(matrix(p53, ncol=3)))
names(p53) <- NULL

set2 <- brewer.pal(8, 'Set2')

b2 <- micf1RAl %>% 
  filter(PatID %in% B2IDs$PatID) %>%
  mutate(Phylum = ifelse(Phylum %in% 
                           c('Firmicutes', 'Proteobacteria'), 
                         Phylum, 'Other')) %>%
  group_by(PatID, Phylum) %>% summarise(RA = sum(RA)) %>%
  ggplot(aes(x = PatID, y = RA, 
             fill = factor(Phylum, 
                           levels = c('Firmicutes', 'Other', 'Proteobacteria')))) + 
  geom_col(position = 'stack') +
  ylab('Relative abundance (%)') + 
  ggtitle('Batch 2') +
  scale_x_discrete(limits = B2IDs$PatID) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45),
        plot.title = element_text(hjust = 0.5),
        legend.position  = 'none') +
  scale_fill_manual(name = 'Phylum', values = c(set2[2], 'gray80', set2[3]))

ggsave('RAPhylum.tiff', dpi = 300, width = 140, height = 80, unit = 'mm')


b <- arrangeGrob(b2, b1, nrow = 1, widths = c(0.38, 0.62))
plot(b)

ggsave('LibSize_vs_RAPhylum.tiff', b, dpi = 300, width = 200, height = 100, unit = 'mm')

#

micf1RAl %>% 
  mutate(Class = ifelse(Class %in% 
                           c('Bacilli', 'Gammaproteobacteria', 'Betaproteobacteria'), 
                         Class, 'Other')) %>%
  group_by(PatID, Class) %>% summarise(RA = sum(RA)) %>%
  ggplot(aes(x = PatID, y = RA, 
             fill = factor(Class, 
                           levels = c('Bacilli', 'Other', 'Gammaproteobacteria', 'Betaproteobacteria')))) + 
  geom_col(position = 'stack') +
  ylab('Relative abundance (%)') + 
  ggtitle('Microbial composition  of sputum samples (Class level)') +
  scale_x_discrete(limits = Xorder$PatID) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = 'Order', values = c(set2[2], 'gray80', set2[3], set2[4]))

ggsave('RAClass.tiff', dpi = 300, width = 150, height = 80, unit = 'mm')

#
library(scales)
show_col(hue_pal()(18))

micfRAl %>% 
  mutate(Order = ifelse(Order %in% 
                          c('Bacillales', 'Lactobacillales',
                            'Burkholderiales', 'Xanthomonadales', 'Pasteurellales'), 
                        Order, 'Other')) %>%
  group_by(PatID, Order) %>% summarise(RA = sum(RA)) %>%
  ggplot(aes(x = PatID, y = RA, 
             fill = factor(Order, 
                           levels = c('Other', 'Bacillales', 'Lactobacillales', 
                                      'Pasteurellales', 'Xanthomonadales', 'Burkholderiales')))) + 
  geom_col(position = 'stack') +
  ylab('Relative abundance (%)') + 
  ggtitle('Microbial composition  of sputum samples (Order level)') +
  scale_x_discrete(limits = Xorder$PatID) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = 'Order', values =   c('gray80', '#E88526', '#D39200', 
                                                 '#00BFC4', '#C77CFF', '#FF61CC'))
                      

ggsave('RAOrder.tiff', dpi = 300, width = 145, height = 80, unit = 'mm')

#

micfRAl %>% 
  mutate(Family = ifelse(Family %in% 
                          c('Streptococcaceae', 'Carnobacteriaceae', 'Micrococcaceae',
                            'Pasteurellaceae', 'Xanthomonadaceae', 'Oxalobacteraceae'), 
                        Family, 'Other')) %>%
  group_by(PatID, Family) %>% summarise(RA = sum(RA)) %>%
  ggplot(aes(x = PatID, y = RA, 
             fill = factor(Family, 
                           levels = c('Streptococcaceae', 'Carnobacteriaceae', 'Micrococcaceae', 'Other',
                                      'Pasteurellaceae', 'Xanthomonadaceae', 'Oxalobacteraceae')))) + 
  geom_col(position = 'stack') +
  ylab('Relative abundance (%)') + 
  ggtitle('Microbial composition  of sputum samples (Family level)') +
  scale_x_discrete(limits = Xorder$PatID) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45),
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(name = 'Family', values = c('#E88526', '#D39200','#B79F00', 'gray80', 
                                               '#00BFC4', '#00ADFA', '#619cFF'))


ggsave('RAFamily.tiff', dpi = 300, width = 150, height = 90, unit = 'mm')


# show specific taxon representatives only
micf1RAl %>% 
  mutate(newGenus = ifelse(
    Genus %in% c('Streptococcus', 'Rothia', 'Granulicatella', 'Ralstonia', 'Lautropia', 'Mycobacterium'), Genus, 'Other')) %>%
    group_by(PatID, newGenus) %>% summarise(RA = sum(RA)) %>%
  ggplot(aes(x = PatID, y = RA, fill = newGenus)) + geom_col(position = 'stack') +
  ylab('Relative abundance (%)') + ggtitle('Order') +
  theme(legend.position = 'bottom',
        axis.text.x = element_text(angle = 45)) + 
  scale_x_discrete(limits = Xorder$PatID) 

# plot only phyla that are selected by the sPLS (stability > 0.6)
c1s1 <- c1s %>% separate(OTU, c('Order', 'Family', 'Genus'), sep = '_') %>%
  mutate(Family = gsub("\\..*","", Family))

c1ord <- unique(c1s1$Order)

c1mic <- micf1RAl %>% filter(Order %in% c1ord) %>%
  group_by(PatID, Order) %>% summarise(RA = sum(RA))

p <- c1mic %>% group_by(Order) %>%
  do(plots = ggplot(data = ., aes(x = PatID, y = RA)) + geom_bar(stat = 'identity') +
       ggtitle(.$Order) + scale_x_discrete(limits = Xorder1$PatID))

marrangeGrob(p$plots, nrow = 3, ncol = 3)

dev.new()
micfRAl %>%
  filter(Genus == 'Aggregatibacter') %>% arrange(Proximal.total.events) %>%
  ggplot(aes(x = fct_inorder(as.factor(PatID)), y = RA)) + geom_bar(stat = 'identity')


# heatmap of filtered data set at OTU level (relative abundance)
micfRA # total count filtered, proportion (RA)

library(pheatmap)

rownames(micfRA) <- rownames(micf)
colnames(micfRA) <- otuNames2$otu

rownames(annot) <- annot$PatID

# add info about pseudomonas infection
annot1 <- data.frame(PatID = rownames(micfRA)) %>% left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))
rownames(annot1) <- annot1$PatID

annot <- annot %>% left_join(clinu2 %>% dplyr::select(PatID, Pseudomonas.colonisation))
annot$Pseudomonas.colonisation[annot$Pseudomonas.colonisation == 'TBC'] <- 'NA'

ann_colors = list(
  dataM = c(Batch1 = "indianred1", Batch2 = "grey"))

tiff('OTU_RA_heatmap.tiff', res = 300, width = 450, height = 280, unit = 'mm')

heatmap <- pheatmap(as.matrix(t(micfRA)), annotation_col = annot %>% dplyr::select(dataM), 
                    show_rownames = TRUE, cluster_cols =  T, clustering_distance_cols = distRA,
                    show_colnames = TRUE, cluster_rows = F,
                    fontsize = 10, annotation_colors = ann_colors,
                    main = 'Relative abundance of OTUs in sputum samples (Genus level) ') 

dev.new()
heatmap
dev.off()

# according to reflux
micfRA <- micfRA %>% mutate(PatID = rownames(.)) %>% left_join(annot1) %>%
  arrange(Proximal.total.events)

rownames(micfRA) <- micfRA$PatID

micfRA <- micfRA %>% dplyr::select(-c(PatID, Proximal.total.events))

# top loadings from multi-block analysis
Mloadings <- Mloadings %>% arrange(desc(abs(comp1)))
MloadingsT <- Mloadings[1:10,]

micfRAT <- micfRA[,MloadingsT$otu]

# create quantile colour gradient
annot <- annot1

annot <- annot %>% arrange(Proximal.total.events)
annot$ecdf <- 1:length(annot$Proximal.total.events)/length(annot$Proximal.total.events)
annot$ecdf[1] <- 0
colval <- color.scale(annot$ecdf, c(0,0,1), c(0,1,1), c(1,0,0))
annot <- cbind(annot, colval)

cols <- annot$colval
names(cols) <- annot$Proximal.total.events
colsl <- list(Proximal.total.events = cols)

tiff('OTU_RA_TL_heatmap.tiff', res = 300, width = 200, height = 200, unit = 'mm')
heatmap <- pheatmap(as.matrix(t(micfRA)), annotation_col = annot1 %>% dplyr::select(Proximal.total.events), 
                    show_rownames = TRUE, cluster_cols =  T, clustering_distance_cols = distRA,
                    show_colnames = FALSE, cluster_rows = F,
                    fontsize = 7,
                    main = 'Relative abundance of OTUs in sputum samples',
                    annotation_colors = colsl,
                    cellheight = 7, cellwidth = 7) 

heatmap
dev.off()
#
#
#

## Batch effect

batch <- data.frame(PatID = rownames(mics1),
                    dataM1 = 'Yes') %>%
  left_join(mr %>% dplyr::select(PatID, dataM)) %>%
  mutate(Batch = ifelse(dataM == 'Yes' & dataM1 == 'Yes',
                        'Batch1', 'Batch2')) %>%
  dplyr::select(-c(dataM, dataM1))

table(batch$Batch)

# distribution of transformed data

par(mfrow = c(4, 4))
for (i in 1: 58) {
  hist(micstf[,i], main = colnames(micstf)[i], cex.main = 0.7)}

# fit linear regression model of batch effect to each OTU
library(lme4)

micstfl <- as.data.frame(micstf) %>% mutate(PatID = rownames(.)) %>%
  left_join(batch) %>%
  pivot_longer(cols = !c(PatID, Batch), names_to = 'otu', values_to = 'value')

regRes <- micstfl %>% group_by(otu) %>%
  do(as.data.frame(summary(lm(value ~ Batch, data = .))$coefficients))

# remove coefficients related to intercept
regRes <- regRes[seq(2, nrow(regRes), 2),]

# plot distributions of regression coefficients related to batch
dev.new()
hist(regRes$Estimate)

regRes %>% ggplot(aes(x = Estimate)) + geom_density()


# compare descriptive statistics
micstfl %>% group_by(Batch) %>% summarise(median = median(value),
                                          mean  = mean(value),
                                          sd = sd(value))

sumtf <- micstfl %>% group_by(otu) %>% summarise(median = median(value))

# density plots (variable-wise batch effect)
d <- micstfl %>% ggplot(aes(x = value, colour = Batch)) + geom_density(size = 1.5) + theme_bw() +
  scale_colour_manual(values = c("indianred1","gray")) +
  theme(plot.title = element_text(hjust = 0.5))

# box plots
b <- micstfl %>% ggplot(aes(x = Batch, y = value, fill = Batch)) + geom_boxplot() + 
  theme_bw() +  scale_fill_manual(values = c("indianred1","gray"))


db <- arrangeGrob(d, b, top = 'Distribution of CLR transformed data', nrow = 2)
plot(db)
ggsave('BatchDistribution.tiff', db, dpi = 300, width = 110, height = 150, unit = 'mm')

# relative log expression plots
dev.new()
micstfl %>% left_join(sumtf) %>% mutate(valueCen = value - median) %>%
  ggplot(aes(x = PatID, y = valueCen, fill = Batch)) + geom_boxplot()  +
  geom_hline(yintercept = 0, linetype = 'dashed')


# sparsity of data between batches
bmicf <- micf %>% mutate(PatID = rownames(.)) %>% left_join(batch)

b1micf <- bmicf %>% filter(Batch == 'Batch1') %>% dplyr::select(-c(PatID, Batch))
b2micf <- bmicf %>% filter(Batch == 'Batch2') %>% dplyr::select(-c(PatID, Batch))

(table(b1micf[b1micf == 0])/(nrow(b1micf)*ncol(b1micf)))*100 # 35.2 % 0s
(table(b2micf[b2micf == 0])/(nrow(b2micf)*ncol(b2micf)))*100 # 41.8 % 0s


#############################
#########################################################



## INTEGRATIVE ANALYSIS

# regression PLS between proximal events and microbiome data (filtered 0.1% min one library)
micstf1 <- as.data.frame(micstf1) # clr, sample filtering
micstf <- as.data.frame(micstf) # clr, total filtering

micf1.CSS # CSS, sample filtering
micf.CSS # CSS, total filtering

# filter datasets to include the same samples
setdiff(ref1$PatID, rownames(micstf))
vdif <- c("CF001", "CF014", "CF015", "CF021", "CF029", "CF032", "CF035", 
          "CF036", "CF045", "CF055", "CF056")

rownames(ref1) <- ref1$PatID

'%ni%' <- Negate('%in%')

r1 <- ref1 %>% filter(PatID %ni% vdif) %>%
  dplyr::select(Proximal.total.events) %>%
  mutate(Proximal.total.events = ifelse(Proximal.total.events == 0, 0.01, Proximal.total.events))



# format the microbiome variable names for visualisations
otuNames <- as.data.frame(colnames(micf)) 
colnames(otuNames) <- 'otu'
otuNames1 <- otuNames %>% separate(otu, c('k', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6'), sep = '__')
otuNames2 <- otuNames1 %>% mutate(T4 = ifelse(is.na(T4) == TRUE, T3, T4)) %>%
  transmute(T4 = gsub('\\..*', '', T4),
            T5 = T5,
            T6 = T6)

otuNames2[is.na(otuNames2) == TRUE] <- ''
otuNames2 <- otuNames2 %>% transmute(otu = paste(T4, T5, T6, sep = '_'))

# colnames(micstf) <- otuNames2$otu
micstf <- as.data.frame(micstf)

# regression PLS
# tuning number of components (dimensions)
tune.pls <- pls(micstf, r1$Proximal.total.events, ncomp = 4 , mode = 'regression')
q.tune.pls <- perf(tune.pls, validation = 'Mfold', folds = 6, nrepeat = 50)
plot(q.tune.pls, criterion = 'Q2.total')

# tuning number of variables for a sparse model
# tuning number of variables
list.keepX <- c(seq(5, 50, 5))

v.tune.pls <- tune.spls(micstf, r1, test.keepX = list.keepX,
                        ncomp = 2,
                        validation = 'Mfold', folds = 5, nrepeat = 50,
                        mode = 'regression', measure = 'MSE')

dev.new()
plot(v.tune.pls) # 1 component, 25 OTUs
ggsave('TuneVarPLS_50rep.tiff', dpi = 300, width = 200, height = 120, unit = 'mm')

write.csv(v.tune.pls$measure.pred, 'SplsTuningResults.csv')


RMpls <- pls(micstf, r1, ncomp = 1 , mode = 'regression')
RMspls <- spls(micstf, r1, ncomp = 1 , mode = 'regression', keepX = 25)

plotIndiv(RMpls, rep.space = 'X-variate')
plotIndiv(RMpls, rep.space = 'Y-variate')
plotIndiv(RMpls)

RMscoresX <- as.data.frame(RMspls$variates$X)

dev.new()
xPlot <- RMscoresX %>% ggplot(aes(y = comp2, x = comp1, colour = r1$Proximal.total.events)) + 
  geom_point(size = 2) +
  geom_text_repel(label = rownames(RMscoresX)) + theme_bw() +
  theme_bw() +
  ggtitle('PLS scores plot of CLR transformed data') +
  geom_hline(yintercept = 0, linetype = 'dashed') + 
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis') +
  xlab('Comp1 (11.16%)') + ylab('Comp2 (13.52%)')

xPlot

ggsave('RMpls_Microbiome_TFCLR.tiff', dpi = 300, width = 150, height = 100, unit = 'mm')


RMscoresY <- as.data.frame(RMspls$variates$Y)

yPlot <- RMscoresY %>% ggplot(aes(y = comp2, x = comp1, colour = r1$Proximal.total.events)) + geom_point() +
  geom_text_repel(label = rownames(RMscoresX)) + theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') + geom_vline(xintercept = 0, linetype = 'dashed')

grid.arrange(xPlot, yPlot, nrow = 1)


cor(RMscoresX$comp1, RMscoresY$comp1)
tiff(filename = 'RMspls_refluxVSXscores.tiff', res = 300, width = 100, height = 80, unit = 'mm')
plot(RMscoresX$comp1, r1$Proximal.total.events,
     ylab = 'Proximal total events',
     xlab = 'PLS scores for microbiome dataset (Comp1)')
dev.off()

# extract 10 top loadings for each components
dev.new()
plotLoadings(RMspls, block = 1, ndisplay = 20,
             name.var = otuNames2$otu)

RMloadingsX <- as.data.frame(RMspls$loadings$X)
RMloadingsY <- as.data.frame(RMspls$loadings$Y)

library(forcats)

RMloadingsX %>%
  mutate(otu = otuNames2$otu,
        Phylum = gsub('\\..*', '', otuNames1$T2)) %>%
  arrange(desc(abs(comp1))) %>%
  filter(row_number() %in% c(1:12)) %>%
  arrange(comp1) %>% 
  ggplot(aes(y = fct_inorder(otu), x = comp1)) + #fill = Phylum)) + 
  geom_bar(stat = 'identity') + 
  theme_bw(base_size = 10) + ylab('OTU') +
  scale_fill_brewer(palette = 'Greys') +
  theme(legend.position = 'bottom',
        legend.justification = "right",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle('Top PLS loadings')

ggsave('RMSloadingsC1.tiff', dpi = 300, width = 150, height = 110, unit = 'mm')

c1 <- RMloadingsX %>% arrange(desc(abs(comp1))) %>%
  filter(row_number() %in% c(1:20))
c1v <- rownames(c1)

c2 <- RMloadingsX %>% arrange(desc(abs(comp2))) %>%
  filter(row_number() %in% c(1:10))
c2v <- rownames(c2)

plotLoadings(RMspls, block = 'X', ndisplay = 15)

#
dev.new()
plotVar(RMpls, cutoff = 0.6, block = 'X')

# calculate relative abundance
micf1RA <- lapply(1:ncol(micf1), function(x){
  micf1[,x] <- micf1[,x]/rowSums(micf1)*100
})

micf1RA <- bind_cols(micf1RA)
colnames(micf1RA) <- otuNames1$otu

RM <- r1 %>% cbind(micstf1) %>% mutate(PatID = rownames(.))
RMl <- RM %>% pivot_longer(cols = !c(PatID, Proximal.total.events), names_to = 'otu', values_to = 'value')

RMplot <- function(i){
  p <- RMl %>% filter(otu == i) %>%
    ggplot(aes(y = Proximal.total.events,
               x = value)) + geom_point() +
    geom_text_repel(aes(label = PatID), size = 2) +
    ggtitle(i) + theme_bw()
  print(p)
  }

RMplotsC1 <- lapply(c1sv, RMplot) 
marrangeGrob(RMplotsC1, nrow = 4, ncol = 3)


# stability of variables selected by the sparse model
perf.RMspls <- perf(RMspls, validation = 'Mfold', folds = 5, nrepeat = 50)

stab.RMspls <- as.data.frame(perf.RMspls$features$stability.X$comp1)
colnames(stab.RMspls) <- 'comp1'


stab.RMspls <- stab.RMspls %>% mutate(OTU = rownames(.))

RMSloadings <- as.data.frame(RMspls$loadings$X) %>% mutate(OTU = rownames(.)) %>%
  full_join(stab.RMspls, by = 'OTU') %>% arrange(desc(abs(comp1.x))) %>%
  filter(comp1.x != 0)

View(RMSloadings %>% arrange(desc(abs(comp1.x))))
plot(abs(RMSloadings$comp1.x), RMSloadings$comp1.y)

c1s <- RMSloadings %>% arrange(desc(abs(comp1.x))) %>%
  filter(row_number() %in% c(1:18))

c1sv <- c1s$OTU

dev.new()
plotLoadings(RMspls, block = 'X', ndisplay = 23)

RMSscores <- as.data.frame(RMspls$variates$X)
cor(RMSscores$comp1, r1$Proximal.total.events)

# dispersion of the top loadings
topLoads <- RMSloadings$OTU[1:12]
dev.new()
micstf[,topLoads] %>% pivot_longer(cols = topLoads, names_to = 'OTU', values_to = 'abundance') %>%
  ggplot(aes(x = OTU, y = abundance)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45))

library(ggmatplot)
dev.new()
ggmatplot(micfRA, micf) + theme(legend.position = 'none') + xlim(0, 0.02) + ylim(0,500)

# plot top loadings CLR abundance values vs reflux variable
micstfTL <- micstf[,topLoads] %>% mutate(PatID = rownames(.)) %>%
  left_join(r1 %>% mutate(PatID = rownames(r1)))

pdf('RMspls_TopLoadings_CLRvsreflux.pdf')
par(mfrow = c(3, 3))
for (i in 1:22) {
  plot(micstfTL[,i], micstfTL[,'Proximal.total.events'],
       main = colnames(micstfTL)[i],
       ylab = 'Proximal.tot.events',
       xlab = 'CLR',
       cex.main = 0.5)}
dev.off()

micfRATL <- micfRA[,topLoads] %>% mutate(PatID = rownames(.)) %>%
  left_join(r1 %>% mutate(PatID = rownames(r1)))

# create barplots with RA change with reflux for the top 4 loadings?
library(forcats)
library(gridExtra)

p4 <- micfRATL %>% pivot_longer(cols =! c('PatID', 'Proximal.total.events'), names_to = 'OTU', values_to = 'RA') %>%
  filter(OTU == 'k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Oxalobacteraceae.g__') %>%
  arrange(Proximal.total.events) %>%
  ggplot(aes(x = fct_inorder(as.factor(PatID)), y = RA)) + geom_col(fill = '#FF61CC') + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45)) +
  xlab('PatID') + ylab('Relative abundance') +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle('Burkholderiales__Oxalobacteraceae')
  
g <- arrangeGrob(p1, p2, p3, p4, nrow = 2)
plot(g)

ggsave('RMspls_topLoads_RAvsReflux.tiff', g, dpi = 300, width = 300, height = 200, unit = 'mm')


#
#
#

## Testing effect of reflux on microbiome with other sparse models

# selbal 
devtools::install_github(repo = "malucalle/selbal")
library(selbal)

micf.CSS # count data normalised for library size
micf.CSSo <- micf.CSS + 1# offset by 1 to remove 0 values
micf.CSSo <- as.data.frame(micf.CSSo)



RMselbal <- selbal(x = micf.CSSo, y = r1$Proximal.total.events)  

dev.new()
grid.draw(RMselbal$global.plot2)

results_selbal <- selbal_wrapper(result = RMselbal, X = micf.CSSo)

# optimise number of OTUs to select
RMselbal1 <- selbal.cv(x = micf.CSSo, y = r1$Proximal.total.events, 
                       n.fold = 5, n.iter = 10, covar = NULL)

dev.off()
grid.draw(RMselbal1$global.plot)
View(RMselbal1$cv.tab)

RMselbal1$global.balance
RMselbal1$glm
RMselbal1$accuracy.nvar

#

# Lasso-clr
library(glmnet)
library(ipflasso)
library(ncvreg)

RMclrlasso <- glmnet(x = micstf, y = r1$Proximal.total.events, nlambda = 20)

RMclrlasso <- glmnet(x = micstf, y = r1$Proximal.total.events, 
                     lambda = seq(1.4930, 2.425, 0.02))

# plot coefficients for each variable for different values of lambda
plot(RMclrlasso, xvar = 'lambda', label = T)

# numerical output
RMclrlasso

# run model with lambda resulting in a desired number of variables
RMclrlasso1 <- glmnet(x = micstf, y = r1$Proximal.total.events, lambda = 3.399404)

# no of selected variables
RMclrlasso1[['df']]

# names of selected variables
View(as.data.frame(colnames(micstf[,RMclrlasso1[['beta']]@i])))

# coefficients
x <- predict(RMclrlasso1, type = 'coef')

y <- data.frame(name = row.names(x),
                val = x[,1]) %>% 
  filter(val != 0) %>%
  arrange(desc(abs(val)))

row.names(y) <- NULL


# CV optimisation of model parameters
cvfit <- cv.glmnet(x = as.matrix(micstf), y = r1$Proximal.total.events,
                   nfolds = 5)
plot(cvfit)

cvfit1 <- cvr.glmnet(X = as.matrix(micstf), Y = r1$Proximal.total.events,
                     nfolds = 5, ncv = 50, family = 'gaussian', type.measure = 'mse')

cvfit1_res <- data.frame(lambda = cvfit1[['lambda']], MSE = cvfit1[['cvm']])

min <- cvfit1_res[25,]

p1 <- cvfit1_res %>% ggplot(aes(x = log(lambda), y = MSE)) + geom_point() +
  theme_bw() + geom_point(data = min, colour = 'red') + 
  geom_vline(xintercept = log(min$lambda), linetype = 'dashed') +
  theme(plot.title = element_text(hjust = 0.5))

coefs <- cvfit1[['coeff']][-1,]
nvar <- coefs
nvar[nvar != 0] <- 1
nvar1 <- colSums(nvar)

cvfit1_res <- cvfit1_res %>% mutate(nvar = nvar1)

p2 <- cvfit1_res %>% ggplot(aes(x = log(lambda), y = nvar)) + geom_col() +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ylab('Variables') + 
  geom_vline(xintercept = log(min$lambda), linetype = 'dashed')

g <- arrangeGrob(p2, p1, ncol = 1, heights = c(0.3, 0.7), top = 'Cross-validation curve for lasso-clr')
plot(g)

ggsave('CVcurve_lassoCLR.tiff', g, width = 100, height = 80, unit = 'mm', dpi = 300)


# lambda for min MSE
cvfit$lambda.min
coef(cvfit, s = 'lambda.min')

mse <- as.data.frame(cvfit1[['cvm']])
colnames(mse) <- 'MSE'
which(mse$MSE == min(mse))
cvfit1[['lambda']][25]

# coefficients for min MSE
coefs7 <- as.data.frame(cvfit1[['coeff']][,c('s25')])
colnames(coefs7) <- 'coef'
coefs7 <- coefs7 %>% mutate(name = rownames(.)) %>% relocate()
row.names(coef7) <- NULL
View(coefs7 %>% filter(coef != 0) %>% arrange(desc(abs(coef))))

# Proportion of variance explained
rsq = 1 - (cvfit$cvm[25]/var(r1$Proximal.total.events))

cvfit2 <- cv.ncvreg(X = micstf, y = r1$Proximal.total.events, penalty="lasso")
plot(cvfit2, type="rsq")

#
#
#

# Integration of proteome and microbiome data
micstf # microbiome - total sample filtered CLR transformed
assayT <- as.data.frame(assayT) # proteome - vsn transformed

setdiff(rownames(assayT), rownames(micstf))

assayT2 <- assayT %>% filter(rownames(assayT) %ni% c("CF015", "CF029", "CF032", "CF035", 
                            "CF001", "CF014", "CF021", "CF036", "CF045", "CF055", "CF056"))

setdiff(rownames(micstf), rownames(assayT))

micstf2 <- micstf %>% filter(rownames(micstf) %ni% c("CF011", "CF040", "CF041"))
colnames(micstf2) <- otuNames2$otu

# tuning number of components (dimensions)
tune.PMpls <- pls(micstf2, assayT2, ncomp = 4 , mode = 'canonical')
q.tune.PMpls <- perf(tune.PMpls, validation = 'Mfold', folds = 5, nrepeat = 50)
plot(q.tune.PMpls, criterion = 'Q2.total')

PMpls <- pls(micstf2, assayT2, ncomp = 1, mode = 'canonical')

plotIndiv(PMpls)
plotVar(PMpls, cutoff = 0.7)

dev.new()
plotLoadings(PMpls, ndisplay = 20, comp = 1)

PMloadingsY <- as.data.frame(PMpls$loadings$Y) %>%
  mutate(Accession = rownames(.)) %>%
  left_join(prot %>% dplyr::select(Accession, Description))

PMloadingsX <- as.data.frame(PMpls$loadings$X) %>%
  mutate(Accession = rownames(.)) 


PMscoresX <- as.data.frame(PMpls$variates$X)
PMscoresY <- as.data.frame(PMpls$variates$Y)

cor(PMscoresX$comp1, PMscoresY$comp1)

PMpls$prop_expl_var

annot2 <- annot %>% filter(PatID %ni% c("CF015", "CF029", "CF032", "CF035", 
                                        "CF001", "CF014", "CF021", "CF036", 
                                        "CF045", "CF055", "CF056", "CF011", 
                                        "CF040", "CF041"))

PMscoresX %>% ggplot(aes(x = comp1, y = comp2, colour = annot2$Proximal.total.events)) +
  geom_point(size = 3) +
  scale_colour_continuous(name = 'reflux', type = 'viridis') +
  ggtitle('Microbiome')

PMscoresY %>% ggplot(aes(x = comp1, y = comp2, colour = annot2$Proximal.total.events)) +
  geom_point(size = 3) +
  scale_colour_continuous(name = 'reflux', type = 'viridis') +
  ggtitle('Proteome')


###############################

#
#
#

# Reflux and protein data
# Sparse canonical correlation analysis
library(PMA)
library(ade4)

assay <- assay %>% mutate_if(is.character, as.numeric)

assayt <- t(assay)

setdiff(rownames(ref3), rownames(assayT))

'%ni%' <- Negate('%in%')

ref4 <- ref3 %>% filter(rownames(ref3) %ni% c("CF011", "CF040", "CF041"))
ref4 <- ref4 %>% mutate_if(is.character, as.numeric)

setdiff(rownames(assayT), rownames(ref4))

assayT1 <- as.data.frame(assayT) %>% mutate(PatID = rownames(assayT)) %>%
  arrange(PatID) %>%
  filter(PatID %ni% c("CF021", "CF027"))

assayT1 <- assayT1 %>% dplyr::select(!PatID)

pcaR <- dudi.pca(ref4, scale = TRUE, scannf = FALSE)
pcaP <- dudi.pca(assayT1, scale = TRUE, scannf = FALSE)

rv1 <- RV.rtest(pcaR$tab, pcaP$tab, 999)
rv1
plot(rv1)

tref4 <- t(ref4)
tassayT1 <- t(assayT1)

# run standard CCA
ccaRes <- CCA(t(tref4), t(tassayT1), typex = 'standard', typez = 'standard', K = 2)
ccaRes

# choosing the best tuning parameters for sCCA
CCA.permute(ref4, assayT1, nperms = 100, typex = 'standard', typez = 'standard')


perm.out <- CCA.permute(ref4, assayT1, nperms = 100,
            penaltyxs = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4),
            penaltyzs = c(0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4))

perm.out

plot(perm.out)


# perform sCCA with chosen parameters
ccaRes <- CCA(t(tref4), t(tassayT1), penaltyx = 0.15, penaltyz = 0.15) # optimising parameter choice?

ccaRes

# multitable PCA
combined <- cbind(t(tref4[ccaRes$u != 0, ]),
                 t(tassayT1[ccaRes$v != 0, ]))

pcaRes = dudi.pca(combined, scale = TRUE, scannf = FALSE, nf = 3)

# compare cca variables with proteomics PC3 loadings
p3l <- p3l %>% mutate(rank = as.numeric(rownames(p3l)))
View(data.frame(Accession = colnames(combined)) %>% left_join(p3l) %>%
  left_join(prot %>% select(Accession, Description)))

p1l <- data.frame(pc1 = loadings1$PC1, 
                  Accession = rownames(loadings1)) %>%
  arrange(pc1)

p1l <- p1l %>% mutate(rank = rownames(p1l))
View(data.frame(Accession = colnames(combined)) %>% left_join(p1l) %>%
       left_join(prot %>% select(Accession, Description)))

# annotation
sampleInfo  <- data.frame(pcaRes$li) %>% mutate(PatID = rownames(.)) #%>%
  left_join(refGroups)
featureInfo <- data.frame(pcaRes$c1, 
                          feature = ifelse(rownames(pcaRes$c1) == 'V1', 'Ref', 'Prot'))

# multitable interpret pca
library(ggrepel)

ggplot() + geom_point(data = sampleInfo,
                       aes(x = Axis1, y = Axis2), size = 3) +
  geom_label_repel(data = featureInfo,
                   aes(x = 5.5 * CS1, y = 5.5 * CS2, label = feature, fill = feature),
                   size = 2, segment.size = 0.3,
                   label.padding = unit(0.1, "lines"), label.size = 0) +
  geom_point(data = featureInfo,
             aes(x = 5.5 * CS1, y = 5.5 * CS2, fill = feature),
             size = 1, shape = 23, col = "#383838") +
  scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#a6d854", "#e78ac3")) +
  guides(fill = guide_legend(override.aes = list(shape = 32, size = 0))) +
  coord_fixed()+
  labs(x = sprintf("Axis1 [%s%% Variance]",
                   100 * round(pcaRes$eig[1] / sum(pcaRes$eig), 2)),
       y = sprintf("Axis2 [%s%% Variance]",
                   100 * round(pcaRes$eig[2] / sum(pcaRes$eig), 2)),
       fill = "Feature Type", col = "Sample Type")

# with prcomp
library(ggfortify)

pcaRes1 <- prcomp(combined, scale = TRUE)
loadingsC <- pcaRes1$rotation
scoresC <- pcaRes1$x

eigen1 <-  get_eigenvalue(pcaRes1)
barplot(eigen1$variance.percent[1:10])

dev.new()
pcca1 <- autoplot(pcaRes1, data = annotC, colour = 'IVCourses12months', 
                 #loadings = TRUE, loadings.label = TRUE, loadings.label.repel = TRUE,
         loadings.label.size = 3) + #data = sampleInfo %>% select(PatID, pc3), colour = 'pc3') +
  #geom_text(aes(label = rownames(combined)), size = 3) +
    theme_bw(base_size = 12) +
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  ggtitle('IV antibiotic days in 12 months') +
  theme(plot.title = element_text(hjust = 0.5)) +
  #scale_colour_continuous(name = 'Number', type = 'viridis')
  scale_colour_binned(n.breaks = 4, type = 'viridis')

  
pcca$layers[[2]]$data <- pcca$layers[[2]]$data['Totals',] 
pcca$layers[[3]]$data <- pcca$layers[[3]]$data['Totals',]

pcca1

ggsave('PCAsCCA.tiff', dpi = 300, width = 130, height = 90, unit = 'mm')  

p <- grid.arrange(pcca1, pcca2, pcca3, nrow = 1, ncol = 3)
ggsave('PCAsCCAClin.tiff', p, dpi = 300, width = 300, height = 70, unit = 'mm')


annotC <- data.frame(PatID = rownames(scoresC)) %>%
  left_join(clin)

annotC$IVCourses12months <- as.factor(annotC$IVCourses12months)

# annotate variables selected by CCA
featureInfoA <- featureInfo %>% mutate(Accession = rownames(.)) %>% 
  left_join(prot %>% select(Accession, Description))

# create new grouping based on multitable PCA results
sampleInfo <- sampleInfo %>% mutate(cca = ifelse(Axis1 > 0, 'cca1', 'cca2')) #%>%
  left_join(refGroups1)

autoplot(pcaRes1, #loadings = TRUE, loadings.label = TRUE, loadings.label.repel = TRUE,
         loadings.label.size = 3, data = sampleInfo, colour = 'refGroup') +
  geom_text_repel(aes(label = rownames(combined)), size = 3) #+
  stat_ellipse(aes(colour = cca), level = 0.95)

# map groups to total events 
library(forcats)

refGroups1 <- refGroups %>% left_join(sampleInfo %>% dplyr::select(PatID, cca)) %>%
  left_join(ref1 %>% dplyr::select(PatID, Total.events, Proximal.total.events)) %>% arrange(Total.events)

refGroups1 %>% ggplot(aes(x = fct_inorder(PatID), y = Total.events, fill = cca)) + 
  geom_bar(stat = 'identity') + theme(axis.text.x = element_text(angle = -90)) +
  geom_hline(yintercept = 75)

# plot protein abundances between CCA/prot PC3 groups
ccaProts <- as.data.frame(t(assayT1)) %>% filter(rownames(.) %in% rownames(featureInfo)) %>%
  t() %>% as.data.frame() %>% mutate(PatID = rownames(.)) %>% left_join(refGroups1 %>% select(PatID, cca)) %>%
  left_join(ref1 %>% select(PatID, Total.events, Proximal.total.events))

ccaBox <- lapply(colnames(ccaProts %>% select(!c(PatID, cca))), function(i){
  ccaProts %>% ggplot(aes(x = cca, y = get(i))) +
    geom_boxplot(aes(fill = cca)) + ggtitle(i) + scale_fill_brewer(palette = 'Dark2') +
    theme(legend.position = 'none')
})

marrangeGrob(ccaBox, nrow = 4, ncol =4)

regRes <- lapply(colnames(ccaProts %>% select(!c(PatID, cca))), function(i){
  subset <- ccaProts %>% select(i, cca)
  colnames(subset)[1] <- 'var'
  m <- lm(var ~ cca,
          data = subset)
  as.data.frame(summary(m)$coefficients)})

regResD <- bind_rows(regRes, .id = 'id') %>%
  rename(p.value = 'Pr(>|t|)') %>%
  filter(grepl('cca', rownames(.))) %>%
  mutate(adj.p = p.adjust(p.value, method = 'BH'))

View(regResD %>% filter(adj.p < 0.05))

# map cca groups to clinical variables
colnames(clin)[1] <- 'PatID'
clin1 <- clin %>% left_join(sampleInfo %>% dplyr::select(PatID, cca))

ccaClin <- lapply(colnames(clin1 %>% dplyr::select(!c(PatID, cca))), function(i){
  clin1 %>% ggplot(aes(x = cca, y = get(i))) +
    geom_boxplot(aes(fill = cca)) + ggtitle(i) + scale_fill_brewer(palette = 'Dark2') +
    theme(legend.position = 'none')
})

marrangeGrob(ccaClin, nrow = 3, ncol = 2)

clin1 %>% arrange(IVCourses12months) %>% ggplot(aes(x = fct_inorder(PatID), y = IVCourses12months, fill = cca)) +
  geom_bar(stat = 'identity')

clin1 %>% arrange(IV.Days.12months) %>% ggplot(aes(x = fct_inorder(PatID), y = IV.Days.12months, fill = cca)) +
  geom_bar(stat = 'identity')

# correlation between clin variables and total events
library(ggcorrplot)

clin1 <- clin1 %>% left_join(ref1 %>% dplyr::select(PatID, Total.events))
cormat <- cor(clin1 %>% dplyr::select(!c(PatID, cca)), method = 'pearson',
              use = 'pairwise.complete.obs')
col <- colorRampPalette(c("red", "white", "blue"))(20)
ggcorrplot(cormat, tl.cex = 10, hc.order = TRUE, legend.title = 'r', lab = TRUE)

# enrichment analysis of protein loadings according to PCA1

loadingsC1 <- as.data.frame(loadingsC) %>% dplyr::select(PC1) %>%
  mutate(Accession = rownames(loadingsC)) %>%
  filter(Accession %ni% c('Totals', 'Mixeds', 'Uprights')) %>%
  left_join(prot[,1:2])

colnames(loadingsC1)[3] <- 'LeadingProt'

library(clusterProfiler)
library("org.Hs.eg.db")


## convert gene symbols into a character vector
prots <- as.vector(loadingsC1$LeadingProt)

# convert accessions to entrez ID
hs <- org.Hs.eg.db

entrezID <- mapIds(hs, 
                   keys = prots,
                   column = "ENTREZID",
                   keytype = "UNIPROT")

# map ENTREZID back to protein accessions and loadings
loadingsC1$ENTREZID <- entrezID

# divide into loadings > 0 and < 0
neg <- loadingsC1 %>% filter(PC1 <0)
pos <- loadingsC1 %>% filter(PC1 >0)

## remove proteins with ENTREXID NAs
neg <- neg %>% drop_na(ENTREZID)
pos <- pos %>% drop_na(ENTREZID)

neg1 <- neg %>% left_join(prot %>% dplyr::select(Accession, Description))

## create background database - all proteins across 4 time points based on 2up converted to gene symbols with UniProt
bg <- prot[,2]

bgID <- mapIds(hs, 
               keys = bg,
               column = "ENTREZID",
               keytype = "UNIPROT")

bgID1 <- data.frame(Reduce(rbind, bgID))
bgID1 <- bgID1 %>% rename(ENTREZID =  Reduce.rbind..bgID.)
bgID2 <- bgID1 %>% drop_na(ENTREZID)
bgIDv <- as.vector(bgID2$ENTREZID)

## convert entrez IDs into vector
entrezIDp <- as.vector(pos$ENTREZID)
entrezIDn <- as.vector(neg$ENTREZID)

## GO enrichment analysis

ggoP <- groupGO(gene = entrezIDp,
               OrgDb = org.Hs.eg.db, ## database - equine not available
               ont = 'BP', ## subontology (cellular component, mol function, bio process)
               level = 3, ## specificity of the term
               readable = TRUE) 

ggoN <- groupGO(gene = entrezIDn,
                OrgDb = org.Hs.eg.db, ## database - equine not available
                ont = 'BP', ## subontology (cellular component, mol function, bio process)
                level = 3, ## specificity of the term
                readable = TRUE) 

# pathway analysis
## reactome pathway analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

pathwaysN <- enrichPathway(gene = entrezIDn,
                            organism = 'human',
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'fdr',
                            qvalueCutoff = 0.9,
                            readable = TRUE,
                           universe = bgID)
head(summary(pathwaysN))

pathwaysP <- enrichPathway(gene = entrezIDp,
                            organism = 'human',
                            pvalueCutoff = 0.05,
                            pAdjustMethod = 'fdr',
                            qvalueCutoff = 0.9,
                            readable = TRUE,
                           universe = bgID)


# merge significant pathways from both datasets and plot together
pathwaysPs <- pathwaysP@result %>% filter(pvalue < 0.05)
pathwaysNs <- pathwaysN@result %>% filter(pvalue < 0.05)

pathways <- rbind(pathwaysPs, pathwaysNs)

pathways$regulation <- c('Downregulated', 'Downregulated', 'Downregulated', 'Upregulated')

plot <- pathways %>% ggplot(aes(x = GeneRatio, y = Description, 
                                colour = pvalue, size = Count)) +
  geom_point(aes(colour = pvalue)) + theme_bw(base_size = 12) +
  facet_wrap(~regulation, scales = 'free', ncol = 1) +
  ylab('Reactome pathway') +
  scale_colour_continuous(limits = c(0.00, 0.05)) +
  scale_size_continuous(limits = c(1, 8))

plot

ggsave('sCCAPathways.tiff', dpi = 300, unit = 'mm', width = 120, height = 80)



################################



# CCA with mixOmics
BiocManager::install("mixOmics")
library(mixOmics)

# scale and center reflux measures (different units)
ref4s <- scale(ref4)
assayT1s <- scale(assayT1)

# shrinkage method
ccaS <- rcc(ref4s, assayT1s, method = 'shrinkage', ncomp = 3)
plot(ccaS, type= 'barplot')

plotIndiv(ccaS)
plotVar(ccaS, cutoff = 0.5)

# ridge method
# tuning parameters

grid1 <- seq(0, 0.8, length = 5)
grid2 <- seq(0.001, 0.05, length = 5)

cv.tune.ccaR <- tune.rcc(ref4s, assayT1s, grid1 = grid1,
                         grid2 = grid2,
                         validation = 'loo')


cv.tune.ccaR

# lambda 1 - 0.2
# lambda 2 - 0.0005

ccaR <- rcc(ref4s, assayT1s, method = 'ridge', lambda1 = 0.075, lambda2 = 0.025)

plot(ccaR, type = 'barplot')

plotIndiv(ccaR)
plotVar(ccaR)



#################################



# canonical PLS with mixOmics - Reflux and protein datasets
library(mixOmics)

# scales and centers variables by default
pls <- pls(assayT1, ref4, ncomp = 2, mode = 'canonical')
pls

plotVar(pls, cutoff = 0.7)

View(as.data.frame(pls$loadings$X) %>% mutate(Accession = rownames(.)) %>%
       left_join(prot %>% dplyr::select(Accession, Description)))

# sparse PLS
# tuning number of components (dimensions)
tune.pls <- pls(assayT1, ref4, ncomp = 4 , mode = 'canonical')
q.tune.pls <- perf(tune.pls, validation = 'Mfold', folds = 6, nrepeat = 5)
plot(q.tune.pls, criterion = 'Q2.total')
ggsave('QcriterionPLS.tiff', dpi = 300, width = 80, height = 80, unit = 'mm')

# tuning number of variables
list.keepX <- c(seq(5, 50, 5))
list.keepY <- c(1:5)

v.tune.pls <- tune.spls(assayT1, ref4, test.keepX = list.keepX,
                        test.keepY = list.keepY, ncomp = 2,
                        validation = 'Mfold', folds = 6, nrepeat = 100,
                        mode = 'canonical', measure = 'cor')

dev.new()
plot(v.tune.pls)
ggsave('TuneVarPLS_100rep.tiff', dpi = 300, width = 200, height = 120, unit = 'mm')

write.csv(v.tune.pls$measure.pred, 'SplsTuningResults.csv')

# keepX = 35, 45
# keepY = 5, 3

# running optimised model
Spls <- spls(assayT1, ref4, ncomp = 2, keepX = c(35,45), 
             keepY = c(5, 3), mode = 'canonical')


# save the outputs of canonical spls
loadingsProt <- as.data.frame(Spls$loadings$X) %>%
  mutate(Accession = rownames(.)) %>% left_join(prot %>% dplyr::select(Accession, Description)) 

loadingsRef <- as.data.frame(Spls$loadings$Y) %>%
  mutate(Var = rownames(.))

loadingsProtS <- as.data.frame(Spls$loadings.star[[1]]) %>%
  mutate(Var = rownames(.))

loadingsRefS <- as.data.frame(Spls$loadings.star[[2]]) %>%
  mutate(Var = rownames(.))

scoresProt <- as.data.frame(Spls$variates$X) %>%
  mutate(PatID = rownames(.))

scoresRef <- as.data.frame(Spls$variates$Y) %>%
  mutate(PatID = rownames(.))

write.csv(loadingsProt, 'SplsProtLoadings.csv')
write.csv(loadingsRef, 'SplsRefLoadings.csv')

write.csv(loadingsProtS, 'SplsProtLoadings_star.csv')
write.csv(loadingsRefS, 'SplsRefLoadings_star.csv')

write.csv(scoresProt, 'SplsProtScores.csv')
write.csv(scoresRef, 'SplsRefScores.csv')


# sample plots
library(RColorBrewer)
pal <- brewer.pal(11, 'Spectral')

annotC <- data.frame(PatID = rownames(ref4)) %>%
  left_join(clin %>% 
              dplyr::select(PatID, IVCourses12months, IV.Days.12months, BaselineFEV1pp)) %>%
  mutate(Condition = 'CF',
         FEVbin = ifelse(BaselineFEV1pp < 30, '20-30',
                         ifelse(BaselineFEV1pp %in% (30:40), '30-40',
                                ifelse(BaselineFEV1pp %in% (40:50), '40-50',
                                       ifelse(BaselineFEV1pp %in% (50:60), '50-60',
                                              ifelse(BaselineFEV1pp %in% (60:70), '60-70', 
                                                     ifelse(BaselineFEV1pp %in% (70:80), '70-80', '80-90')))))))

joint <- plotIndiv(Spls,
          rep.space = 'XY-variate',
          group = annotC$IVCourses12months,
          title = 'IV antibiotic courses in 12 months',
          legend = TRUE,
          col = pal[1:6],
          legend.title = 'Number',
          pch = as.factor(annotC$Condition))


joint <- plotIndiv(Spls,
                   rep.space = 'XY-variate',
                   group = annotC$FEVbin,
                   title = 'Baseline Forced Expiratory Volume',
                   legend = TRUE,
                   legend.title = 'FEV (%)',
                   col = pal[1:8],
                   pch = as.factor(annotC$Condition))


pdf('sPLSIndPlotJoint1.pdf', width = 8, height = 10)
joint
dev.off()

scoresProt %>% ggplot(aes(x = comp1, y = comp2)) + 
  geom_point(aes(colour = ref4$Proximal.totals), size = 3) + theme_bw()


# not annotated
annotC <- annotC %>% left_join(scoresRef %>% dplyr::select(PatID, comp1)) %>%
  mutate(comp1 = comp1 + 4)

ref <- plotIndiv(Spls,
                 rep.space = 'Y-variate',
                 title = 'Reflux data',
                 group = annotC$comp1,
                 pch = as.factor(annotC$Condition),
                 legend = FALSE)

ref


prot <- plotIndiv(Spls,
                  rep.space = 'X-variate',
                  group = annotC$comp1,
                  title = 'Protein data',
                  pch = as.factor(annotC$Condition),
                  legend = FALSE)

prot



# manual biplot
plotIndiv(Spls, rep.space = 'Y-variate', style = 'graphics')
arrows(0, 0, Spls$loadings.star[[2]][,1]*3, Spls$loadings.star[[2]][,2]*3, length = 0.1)

plotIndiv(Spls, rep.space = 'X-variate', style = 'graphics')
arrows(0, 0, Spls$loadings.star[[1]][,1]*10, Spls$loadings.star[[1]][,2]*10, length = 0.1)

View(Spls$loadings.star[[1]])

# agreement between sample projections into C1/C2 space


plotArrow(Spls,
          legend = FALSE,
          group = annotC$comp1) +
  ggtitle('Agreement between protein and reflux sample plots') 


cor(scoresProt$comp1, scoresRef$comp1)
cor(scoresProt$comp2, scoresRef$comp2)

library(IMIFA)
library(vegan)

proc <- Procrustes(as.matrix(scoresProt[,-3]), as.matrix(scoresRef[,-3]), sumsq = TRUE,
                   translate = TRUE)

proc1 <- vegan::procrustes(Y = as.matrix(scoresRef[,-3]), X = as.matrix(scoresProt[,-3]))
proc1

# procrustes ss = 316.4
ss <- sum((proc1$X - proc1$Yrot)^2)
sm <- sum((proc1$X - colMeans(proc1$X))^2)

d <- ss/sm

# calculate procrustes dissimilarity score d
# calculate distances between X and column means of X

p1  <- plot(proc1, kind = 1, to.target = FALSE, type = 'text', cex = 0.6)
p1

proc2 <- vegan::procrustes(Y = as.matrix(scoresRef[,-3]), X = as.matrix(scoresProt[,-3]),
                           symmetric = TRUE)
proc2

# Mantel correlation of matrices

tiff(filename = 'SplsProcrustes.tiff', width = 180, height = 160, unit = 'mm', res = 300)
p2  <- plot(proc2, kind = 1, type = 'text', to.target = FALSE, cex = 0.6)
p2
dev.off()

perm <- protest(as.matrix(scoresProt[,-3]), as.matrix(scoresRef[,-3]), permutations = 1000000)
perm

# variable correlation plot
plotVar(Spls, cutoff = 0.5)
plotVar(Spls) 

Spls$prop_expl_var

dev.new()
cim(Spls, comp = 1:2)

# loadings vectors for each component
vip.spls <- vip(Spls)

# variables selected for C1
c1Vars <- as.data.frame(vip.spls[selectVar(Spls, comp = 1)$X$name,1])

perf.Spls <- perf(Spls, validation = 'Mfold', folds = 6, nrepeat = 50)

# proteins
stab.SplsC1 <- perf.Spls$features$stability.X$comp1
stab.SplsC2 <- perf.Spls$features$stability.X$comp2

stab1 <- as.data.frame(stab.SplsC1[1:35]) 
stab2 <- as.data.frame(stab.SplsC2[1:45]) 

colnames(stab1) <- 'FreqOccur'
colnames(stab2) <- 'FreqOccur'

write.csv(stab1, 'SplsProtC1Stability.csv')
write.csv(stab2, 'SplsProtC2Stability.csv')

stab1f <- stab1 %>% filter(FreqOccur > 0.7)
stab2f <- stab2 %>% filter(FreqOccur > 0.7)

# reflux
stab.SplsC1R <- perf.Spls$features$stability.Y$comp1
stab.SplsC2R <- perf.Spls$features$stability.Y$comp2

stab1R <- as.data.frame(stab.SplsC1R[1:5]) 
stab2R <- as.data.frame(stab.SplsC2R[1:3]) 

colnames(stab1R) <- 'FreqOccur'
colnames(stab2R) <- 'FreqOccur'

write.csv(stab1R, 'SplsRefC1Stability.csv')
write.csv(stab2R, 'SplsRefC2Stability.csv')

stab1Rf <- stab1R %>% filter(FreqOccur > 0.7)
stab2Rf <- stab2R %>% filter(FreqOccur > 0.7)


# extract stability measure for variables selected by Spls
extr.stab.Spls.C1 <- stab.Spls[selectVar(Spls, comp = 1)$X$name]

# create vector of high stability variables for circle correlation plot
Xnames <- as.data.frame(colnames(Spls$X))
colnames(Xnames) <- 'prots'
Ynames <- as.data.frame(colnames(Spls$Y))
colnames(Ynames) <- 'refs'

stabP <- rbind(stab1f, stab2f) %>%
  mutate(prots = rownames(.),
         VarX = rownames(.))

stabY <- rbind(stab1Rf, stab2Rf) %>%
  mutate(refs = rownames(.),
         VarY = rownames(.))

Xnames <- Xnames %>% full_join(stabP)
Ynames <- Ynames %>% full_join(stabY)

listName <- list(Xnames$VarX, Ynames$VarY) 

# replace accession numbers with gene symbols
Xnames <- Xnames %>% left_join(prot %>% dplyr::select(Accession, Description), by = c('prots' = 'Accession'))
Xnames <- Xnames %>% separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(prots, Gene)
Xnames <- Xnames %>% full_join(stabP)
Xnames <- Xnames %>% mutate(Gene = ifelse(is.na(FreqOccur) == TRUE, NA, Gene))


listName <- list(Xnames$Gene, Ynames$VarY) 

# variable correlation plot
plotVar(Spls, 
        cutoff = 0.5,
        var.names = listName,
        title = 'Correlation plot of the selected variables')

plotVar(Spls,
        var.names = listName,
        title = 'Correlation plot of the most stable selected variables')

View(Spls$loadings$X)

# extract loadings for stable variables
Xloadings <- as.data.frame(Spls$loadings$X) %>% 
  mutate(prots = rownames(.))
Yloadings <- as.data.frame(Spls$loadings$Y) %>%
  mutate(refs = rownames(.))

Xnames <- Xnames %>% left_join(Xloadings) %>% drop_na()
Ynames <- Ynames %>% left_join(Yloadings) %>% drop_na()

XnamesN <- Xnames %>% filter(comp1 < 0) %>% 
  left_join(prot %>% dplyr::select(Accession, Description), by = c('prots' = 'Accession'))
XnamesP <- Xnames %>% filter(comp1 > 0) %>% 
  left_join(prot %>% dplyr::select(Accession, Description), by = c('prots' = 'Accession'))

XnamesC2 <- Xnames %>% filter(comp2 != 0) %>% 
  left_join(prot %>% dplyr::select(Accession, Description), by = c('prots' = 'Accession'))

# compare with CCA output
loadingsC2 <- loadingsC1 %>%
  left_join(prot %>% dplyr::select(Accession, Description))



############################

# predicting one reflux measure based on sputum proteome
assayT <- as.data.frame(assayT) %>% arrange(rownames(.))

y <- as.data.frame(ref1$Proximal.total.events)
colnames(y) <- 'Proximal'
rownames(y) <- ref1$PatID

'%ni%' <- Negate('%in%')

setdiff(rownames(y), rownames(assayT))
y <- y %>% filter(rownames(.) %ni% c("CF011", "CF040", "CF041"))
y1 <- y %>% mutate(Proximal = ifelse(Proximal == 0, 0.01, Proximal))

# sparse PLS
# tuning number of components (dimensions)
tune.pls <- pls(assayT, y1, ncomp = 4 , mode = 'regression')
q.tune.pls <- perf(tune.pls, validation = 'Mfold', folds = 6, nrepeat = 50)
plot(q.tune.pls, criterion = 'Q2.total')
ggsave('QcriterionPLS.tiff', dpi = 300, width = 80, height = 80, unit = 'mm')

# tuning number of variables for sparse PLS
list.keepX <- c(seq(5, 50, 5))

v.tune.pls <- tune.spls(assayT, y1, test.keepX = list.keepX,
                        ncomp = 2,
                        validation = 'Mfold', folds = 6, nrepeat = 50,
                        mode = 'regression', measure = 'MSE')


plot(v.tune.pls)


# running optimised model
spls1 <- spls(assayT, y1, ncomp = 2, mode = 'regression', keepX = 50)
pls1 <- pls(assayT, y1, ncomp = 2, mode = 'regression')

# variation explained
spls1$prop_expl_var$X
tune.pls$prop_expl_var$X
pls1$prop_expl_var$X

# scores plots
plotIndiv(spls1,
          rep.space = 'X-variate',
          group = y$Proximal)
          

plotIndiv(pls1,
          rep.space = 'X-variate',
          group = y1$Proximal,
          legend = TRUE)

scoresX <- as.data.frame(pls1$variates$X) 
scoresY <- as.data.frame(pls1$variates$Y)

y2 <- y1 %>% mutate(PatID = rownames(.))
scoresX <- scoresX %>% mutate(PatID = rownames(.)) %>% left_join(y2)
scoresY <- scoresY %>% mutate(PatID = rownames(.)) %>% left_join(y2)

scoresX %>% ggplot(aes(y = comp2, x = comp1, colour = Proximal)) + geom_point(size = 3) + 
  theme_bw() +  scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis') +
  xlab('Comp1 (10.4%)') + ylab('Comp2 (11.7%)') +
  ggtitle('PLS scores plot of protein data') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('PLSscoresProt.tiff', width = 200, height = 140, unit = 'mm', dpi = 300)


scoresY %>% ggplot(aes(y = comp2, x = comp1, colour = Proximal)) + geom_point(size = 3) + 
  theme_bw() +  scale_colour_continuous(name = 'Proximal.total.events', type = 'viridis') 

loadingsXs <- as.data.frame(spls1$loadings$X) %>% mutate(Accession = rownames(.)) %>% left_join(prot %>% dplyr::select(Accession, Description)) 
loadingsX <- as.data.frame(pls1$loadings$X) %>% mutate(Accession = rownames(.)) %>% 
  left_join(prot %>% dplyr::select(Accession, Description)) 

plotVar(pls1, cutoff = 0.5)

rownames(Xnames) <- Xnames$prot
plotLoadings(pls1, block = 'X', comp = 2, ndisplay = 20, name.var = Xnames$Gene,
             title = 'Loadings on component 2 \n protein data', size.title = 1.5)



# relationship between components
plot(pls1$variates$X, pls1$variates$Y)
plot(spls1$variates$X, spls1$variates$Y)

cor(pls1$variates$X, pls1$variates$Y)
cor(spls1$variates$X, spls1$variates$Y)

# validation of top C1 loadings associated with reflux
dev.new()
cim(pls1, mapping = 'XY')
cim.res <- cim(pls1)

loadingsX <- loadingsProt

loadingsX <- loadingsX %>% arrange(desc(abs(comp1)))
top_prots <- loadingsX[1:10,]
top_prots <- top_prots$Accession
top_prots1 <- assayT[,top_prots]

join <- ref1 %>% dplyr::select(PatID, Total.events, Proximal.total.events) %>% left_join(top_prots1 %>% mutate(PatID = rownames(.))) 
joinl <- join %>% pivot_longer(cols = !c(PatID, Proximal.total.events, Total.events), names_to = 'prot', values_to = 'intensity') %>%
  left_join(Xnames, by = c('prot' = 'prots'))

plots <- joinl %>% group_by(prot) %>%
  do(plots = ggplot(data = ., aes(x = intensity, y = Total.events)) + 
       geom_point() + ggtitle(.$Gene) + theme_bw())

marrangeGrob(plots$plots, nrow = 4, ncol = 3)


# try to run regression reflux vs selected proteins
summary(lm(P99999 ~ Proximal, data = join))


###############################


# Integration of clinical and reflux data (PLS canonical)
library(mixOmics)

clin3 <- clin %>% 
  dplyr::select(BaselineFEV1pp, IV.Days.12months)
rownames(clin3) <- clin$PatID

rownames(ref) <- ref$PatID
ref3 <- ref[-1,-1]

setdiff(rownames(clin3), rownames(ref3))
clin3 <- clin3 %>% filter(rownames(.) %ni% c('CF021', 'CF027'))

pls.cr <- pls(X = ref3, Y = clin3, ncomp = 2, scale = TRUE, mode = 'canonical')

plotVar(pls.cr, cutoff = 0.5)

View(pls.cr$loadings.star[[1]])
vip.pls.cr <- vip(pls.cr)

plot(x = ref3$Total.events, y = clin3$IVCourses12months)
text(x = ref3$Total.events, y = clin3$IVCourses12months, labels = rownames(ref3), 
     cex = 0.7, pos = 1)

plot(x = ref3$Total.events, y = clin3$IV.Days.12months)#, xlim = c(0, 200))
text(x = ref3$Total.events, y = clin3$IV.Days.12months, labels = rownames(ref3), 
     cex = 0.7, pos = 1)

plot(x = ref3$Demeester.score, y = clin3$BaselineFEV1pp)


plotIndiv(pls.cr)
biplot(pls.cr, block = 'X')

cor(clin3$IV.Days.12months[-c(37, 40)], ref3$Total.events[-c(37, 40)])

plotArrow(pls.cr)

join <- cbind(clin3, ref3)

hist(join$IV.Days.12months)
hist(join$Total.events)

join <- join %>% mutate(IVDays1 = ifelse(IV.Days.12months == 0, 0.01, IV.Days.12months))

join$IV.Days.12months <- as.numeric(join$IV.Days.12months)

join <- join %>% mutate(PatID = rownames(join))

clin2 <- clin2 %>% mutate(PatID = ifelse(
  Patient.ID %in% (1:9), paste('CF00', Patient.ID, sep =''),
  paste('CF0', Patient.ID, sep ='')))
clin2$CFTR.Rx[38] <- 'No'

join <- join %>% full_join(clin2 %>% dplyr::select(PatID, Age, Sex, CFTR.Rx, Gene1, Gene2, BMI))
join$CFTR.Rx[38] <- 'No'

n_distinct(join$Gene1)

# remove reflux outlier
r1 <- glm(IVDays1 ~ Total.events,
          data = join,
          family = Gamma(link = 'identity'))

summary(r1)

r <- glm(IVDays1 ~ Total.events + Age + Sex + BMI + CFTR.Rx + Gene1,
         data = join,
         family = Gamma(link = 'identity'))

summary(r)

data

r2 <- glm(IVDays1 ~ Total.events + Age + BMI + CFTR.Rx,
          data = join %>% filter(IVDays1 < 100),
          family = Gamma(link = 'inverse'))

summary(r2)



#################################
#######################################################



# MULTI-BLOCK ANALYSIS

# multiblock PLS/sPLS - Reflux, Proteomic and Microbiome variables
r1 <- ref1 %>%
  dplyr::select(Proximal.total.events) %>%
  mutate(Proximal.total.events = ifelse(Proximal.total.events == 0, 0.01, Proximal.total.events))
rownames(r1) <- ref1$PatID

assayT <- assayT %>% as.data.frame() %>% 
  mutate(PatID = rownames(.)) %>%
  arrange(PatID) 

colnames(micstf) <- otuNames2$otu
micstf <- as.data.frame(micstf)

micstf # clr transformed total filtered microbiome data
r1 # reflux proximal total events
assayT # sputum proteome

# number of samples in each dataset
library(gplots)

list <- list(Reflux = rownames(r1), 
             Protein = rownames(assayT), 
             Microbiome = rownames(micstf))

tiff(filename = 'VennRMP.tiff', width = 100, height = 110, unit = 'mm', res = 300)
v <- venn(list)
v
dev.off()   

print(v)

# vector of complete observations across 3 datasets
tr <- c("CF003", "CF005", "CF006", "CF007", "CF008", "CF009", "CF010", "CF012", "CF017", "CF018", 
        "CF020", "CF022", "CF023", "CF024", "CF025", "CF027", "CF033", "CF034", "CF038", "CF039",
        "CF042", "CF048", "CF050", "CF051", "CF052", "CF053", "CF057")

# create training datasets (n = 27)
rtr <- r1 %>% mutate(PatID = rownames(.)) %>%
  filter(PatID %in% tr) %>%
  arrange(PatID)

rownames(rtr) <- rtr$PatID

rtr <- rtr %>% subset(select = -(PatID))

#

protr <- assayT %>% filter(PatID %in% tr) %>%
  arrange(PatID)

rownames(protr) <- protr$PatID

protr <- protr %>% subset(select = -c(PatID))

#

mictr <- as.data.frame(micstf) %>% mutate(PatID = rownames(.)) %>%
  filter(PatID %in% tr) %>%
  arrange(PatID)

rownames(mictr) <- mictr$PatID
  
mictr <- mictr %>% subset(select = -(PatID))

# code datasets for multiblock analysis
X <- list(prot = protr, mic = mictr) # predictor datasets
Y <- as.matrix(rtr) # outcome dataset (variable)

# create design matrix
design <- matrix(0.1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# run multiblock PLS with full design
library(mixOmics)

RMPpls <- block.pls(X, Y, design = design, ncomp = 2, mode = 'regression')

RMPpls$design

RMPpls$prop_expl_var

# sample ordination
dev.new()
plotIndiv(RMPpls)
plotIndiv(RMPpls, block = 'prot',
          legend = FALSE,
          subtitle = 'Proteome')

plotIndiv(RMPpls, block = 'mic',
          legend = FALSE,
          subtitle = 'Microbiome')

Pscores <- RMPpls$variates$prot %>% as.data.frame() %>% mutate(block = 'Proteome', PatID = rownames(.))
Mscores <- RMPpls$variates$mic %>% as.data.frame() %>% mutate(block = 'Microbiome', PatID = rownames(.))
Scores <- rbind(Pscores, Mscores) %>% left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))

Scores %>% ggplot(aes(x = comp1, y = comp2, colour = Proximal.total.events)) + 
  geom_point(size = 2) + scale_color_continuous(type = 'viridis') + theme_bw() +
  geom_text_repel(label = Scores$PatID, size = 2) +
  facet_wrap(~ block, scales = 'free') +
  labs(title = 'Sample ordination in PLS blocks') +
  theme(plot.title = element_text(hjust = 0.5))

ggsave('BPls_ordinationPlots.tiff', dpi = 300, width = 300, height = 130, unit = 'mm')

#
#
#

# agreement between sample ordination across datasets
dev.new()

plotArrow(RMPpls,
          legend = FALSE) +
  ggtitle('Agreement between protein and reflux sample plots') 


cor(Pscores$comp1, Mscores$comp1)
cor(Pscores$comp2, Mscores$comp2)

cor(Pscores$comp2, Y)
cor(Mscores$comp2, Y)

library(IMIFA)
library(vegan)

proc <- Procrustes(as.matrix(Pscores[,1:2]), as.matrix(Mscores[,1:2]), sumsq = TRUE,
                   translate = TRUE)

proc1 <- vegan::procrustes(Y = as.matrix(Pscores[,1:2]), X = as.matrix(Mscores[,1:2]),
                           symmetric = TRUE)


proc1

# procrustes ss = 171.3
ss <- sum((proc1$X - proc1$Yrot)^2)
sm <- sum((proc1$X - colMeans(proc1$X))^2)

d <- ss/sm

# calculate procrustes dissimilarity score d
# calculate distances between X and column means of X

dev.new()
tiff(filename = 'BplsProcrustes1.tiff', width = 180, height = 160, unit = 'mm', res = 300)
p2  <- plot(proc1, kind = 1, type = 'p', to.target = FALSE, cex = 0.6)
p2
dev.off()

perm <- protest(as.matrix(Mscores[,1:2]), as.matrix(Pscores[,1:2]), permutations = 1000000)
perm

#
#
#


annClin <- data.frame(PatID = rownames(rtr)) %>%
  left_join(clin %>% dplyr::select(PatID, BaselineFEV1pp, IV.Days.12months))


print(r1$Proximal.total.events)
  
rownames(annClin) <- annClin$PatID
cor(mictr$Actinomycetales_Actinomycetaceae.g_Actinomyces, RMPpls$variates$mic)

# correlation between variables and components
plotVar(RMPpls, legend = TRUE, cutoff = 0.7)
corC <- plotVar(RMPpls)
corC1 <-  corC %>% mutate(x = round(x, digits = 2),
                          y = round(y, digits = 2))

write.csv(corC1, 'BPLS_CorrelationCoefs_C1Y.csv')

corCf <- corC1 %>% filter(abs(x) >= 0.7) %>%
  left_join(prot %>% dplyr::select(Accession, Description), by = c('names' = 'Accession'))

corCf <- corCf %>% separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(-c(PE, SV)) %>%
  mutate(Gene = ifelse(is.na(Gene) == TRUE, names, Gene),
         CorRef = ifelse(x > 0, 'pos', 'neg')) 

# heatmap of z-score between variables and the outcome (reflux)
conc <- cbind(protr, mictr)
concf <- conc[,c(corCf$names)]
concfs <- as.data.frame(scale(concf, center = TRUE, scale = TRUE))
colnames(concfs) <- corCf$Gene

concfs <- concfs[match(rownames(as.data.frame(Y) %>% arrange(Proximal.total.events)), rownames(concfs)),]

rows <- data.frame(names = colnames(concf)) %>% 
  left_join(corCf %>% dplyr::select(names, Block, Gene, CorRef))
rownames(rows) <- rows$Gene
rows <- rows %>% dplyr::select(CorRef)

annot1 <- an %>% filter(PatID %in% rownames(Y)) %>% dplyr::select(Proximal.total.events)

dev.new()

annot <- data.frame(PatID = rownames(protr)) %>% left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))
rownames(annot) <- annot$PatID

annot <- annot %>% arrange(Proximal.total.events)

annot$ecdf <- 1:length(annot$Proximal.total.events)/length(annot$Proximal.total.events)
annot$ecdf[1] <- 0
colval <- color.scale(annot$ecdf, c(0,0,1), c(0,1,1), c(1,0,0))
annot <- cbind(annot, colval)

cols <- annot$colval
names(cols) <- annot$ecdf
colsl <- list(ecdf = cols)

tiff('BPls_CorComp07_heatmap_quant_ordi.tiff', res = 300, width = 400, height = 250, unit = 'mm')
heatmap1 <- pheatmap(t(concfs), annotation_col = annot %>% dplyr::select(ecdf), 
                     cluster_cols = T, 
                     annotation_row = rows, cluster_rows = T,
                     show_colnames = F,
                     annotation_colors = colsl,
                     treeheight_col = 70)

heatmap1
dev.off()

plot(concfs$Burkholderiales_Oxalobacteraceae.g_Ralstonia, concfs$`IGKV3-20`)

map <- heatmap1
map[[4]] # plot

ids <- data.frame(PatID = rownames(concfs)) %>%
                    mutate(row = row_number())


idOrd <- ids[match(heatmap1$tree_col$order, ids$row),] %>%
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))

bar <- idOrd %>% ggplot(aes(x = fct_inorder(PatID), y = Proximal.total.events)) + 
  geom_col() + theme(axis.title.x = element_blank(),
                     axis.text.x = element_blank())           
bar

dev.new()
plot(arrangeGrob(bar, map[[4]], nrow = 2))

#
#
#

# enrichment analysis of proteins most correlated with outcome variable
corCf1 <- corCf %>% filter(Block == 'prot')
corN <- corCf1 %>% filter(x < 0)
corP <- corCf1 %>% filter(x > 0)

# repeat hierarchical clustering
distC <- vegdist(concfs, method = 'euclidean')
clustC <- hclust(distC, method = 'complete')

dev.new()
plot(clustC)

# cut the tree and assign cluster membership
clustm <- cutree(clustC, h = 12)
clustM <- data.frame(PatID = clustC[['labels']], cluster = clustm)

# boxplots of reflux events according to clusters (w/o outlier)
# remove the outlier

clustM1 <- clustM %>% left_join(ref1) 

clustM1 %>% ggplot(aes(x = as.factor(cluster), y = Proximal.total.events, colour = as.factor(cluster))) +
  geom_boxplot() +
  geom_jitter() +
  theme_bw()

clustM1l <- clustM1 %>% pivot_longer(cols =! c(PatID, cluster), names_to = 'parameter', values_to = 'value')

p <- clustM1l %>% group_by(parameter) %>%
  do(plots = ggplot(data = ., aes(x = cluster, y = value)) + geom_boxplot() + ggtitle(.$parameter))

marrangeGrob(p$plots, nrow = 4, ncol = 4)

# plot top correlated proteins
corCf2 <- corCf %>% filter(abs(x) > 0.8) 

topProts <- assayT[,corCf2$names] %>%
  mutate(PatID = rownames(.)) %>%
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))

p4 <- topProts %>% ggplot(aes(x = P22894, y = Proximal.total.events)) + geom_point() + theme_bw() +
  ggtitle('GSN') + theme(plot.title = element_text(hjust = 0.5))


p <- grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
ggsave('TopCorProts.tiff', p, dpi = 300, width = 100, height = 100, unit = 'mm')

plot(assayT$P22894, assayT$A0A0C4DH38)

#
#
#


# similarity score between variables from different blocks           
Pnames <- data.frame(names = colnames(protr)) %>%
  left_join(prot %>% dplyr::select(Accession, Description), by = c('names' = 'Accession')) %>% 
  separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(-c(PE, SV))

rownames(Pnames) <- Pnames$names


tiff('BPls_CorVars_heatmap.tiff', res = 300, width = 300, height = 200, unit = 'mm')
dev.new()
cim(RMPpls, cutoff = 0.6, margins = c(7, 25), transpose = TRUE, row.names = Pnames$Gene,
    row.cex = 1.2, col.cex = 1)
dev.off()

# extract similarity scores
corMat <- cim(RMPpls)
corMat1 <- corMat$mat %>% as.data.frame()

# filter based on arbitrary thershold
corF <- which(abs(corMat1) > 0.55, arr.ind = TRUE) %>% as.data.frame() %>% distinct(row) %>%
  mutate(prot = rownames(.))

corF1 <- which(abs(t(corMat1)) > 0.55, arr.ind = TRUE) %>% as.data.frame() %>% distinct(row) %>%
  mutate(otu = rownames(.))

corF$prot <- sub('\\.',';', corF$prot)  

corMatF <- corMat1[corF$prot,corF1$otu]

# repeat hierarchical clustering
distF <- vegdist(t(corMatF), method = 'euclidean')
clustF <- hclust(distF, method = 'complete')

distFp <- vegdist(corMatF, method = 'euclidean')
clustFp <- hclust(distFp, method = 'complete')

dev.new()
plot(clustFp)

# cut the tree and assign cluster membership
clustm <- cutree(clustF, h = 4)
clustM <- data.frame(otu = clustF[['labels']], cluster = clustm)

clustmp <- cutree(clustFp, h = 2)
clustMp <- data.frame(prot = clustFp[['labels']], cluster = clustmp)

clustMp %>% filter(cluster == 2) %>% rownames(.)

#

# loadings
dev.new()
plotLoadings(RMPpls, ndisplay = 10)

Ploadings <- RMPpls$loadings$prot %>% as.data.frame() %>%
  mutate(Accession = rownames(.)) %>% 
  left_join(prot %>% dplyr::select(Accession, Description)) %>%
  separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(-c(PE, SV))

Mloadings <- RMPpls$loadings$mic %>% as.data.frame() %>%
  mutate(otu = rownames(.))

pc1 <- Ploadings %>% arrange(desc(abs(comp1))) %>%
  slice(1:10) %>% arrange(desc(comp1)) %>% 
  ggplot(aes(x = fct_inorder(Gene), y = comp1)) + 
  geom_col() + coord_flip() + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.title.x = element_blank()) +
  ggtitle('Comp1') +
  theme(plot.title = element_text(hjust = 0.5))

pc2 <- Ploadings %>% arrange(desc(abs(comp2))) %>%
  slice(1:10) %>% arrange(desc(comp2)) %>% 
  ggplot(aes(x = fct_inorder(Gene), y = comp2)) + 
  geom_col() + coord_flip() + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.title.x = element_blank()) +
  ggtitle('Comp2') +
  theme(plot.title = element_text(hjust = 0.5))

pc2

pc <- arrangeGrob(pc1, pc2, nrow = 2, top = 'Proteome PLS loadings')
plot(pc)

ggsave('Bpls_prot_loadings.tiff', pc, dpi = 300, width = 80, height = 150, unit = 'mm')

#

mc1 <- Mloadings %>% arrange(desc(abs(comp1))) %>%
  slice(1:10) %>% arrange(desc(comp1)) %>% 
  ggplot(aes(x = fct_inorder(otu), y = comp1)) + 
  geom_col() + coord_flip() + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.title.x = element_blank()) +
  ggtitle('Comp1') +
  theme(plot.title = element_text(hjust = 0.5))

mc1

mc2 <- Mloadings %>% arrange(desc(abs(comp2))) %>%
  slice(1:10) %>% arrange(desc(comp2)) %>% 
  ggplot(aes(x = fct_inorder(otu), y = comp2)) + 
  geom_col() + coord_flip() + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.title.x = element_blank()) +
  ggtitle('Comp2') +
  theme(plot.title = element_text(hjust = 0.5))

mc2

mc <- arrangeGrob(mc1, mc2, nrow = 2, top = 'Microbiome PLS loadings')
plot(mc)

ggsave('Bpls_mic_loadings.tiff', mc, dpi = 300, width = 150, height = 150, unit = 'mm')



# correlation between components
cor(RMPpls$variates$prot, RMPpls$variates$mic) # 0.71, 0.75
cor(RMPpls$variates$prot, RMPpls$variates$Y) # 0.66
cor(RMPpls$variates$mic, RMPpls$variates$Y) # 0.72

dev.new()
plotArrow(RMPpls)

library(MatrixCorrelation)
rv2 <- RV2(mictr, protr, center = TRUE, impute = FALSE)
rv2

rvadj <- RVadj(mictr, protr, version = 'Maye', center = TRUE)
rvadj

# prediction & model performance

pred <- predict(RMPpls, newdata = test)
out <- Y %>% filter(rownames(.) %in% rownames(prtst))

dev.new()
plot(pred$AveragedPredict[,,1], out$Proximal.total.events)
abline(1,1)
