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
mic <- read.csv('GORCF sputum TEST.csv')
clin <- read.csv('GOR_ClinicalEndpoint.csv')
clin2 <- read.csv('GORCF data 24-05.csv') %>% drop_na()

# number of patients with complete set of observations
r <- data.frame(PatID = ref$PatID, dataR = 'Yes')
p <- data.frame(PatID = colnames(assay), dataP = 'Yes')
m <- data.frame(PatID = mic[1:17,1], dataM = 'Yes')

View(r %>% full_join(p) %>% full_join(m) %>% drop_na())
mr <- r %>% full_join(m) %>% drop_na()

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

# map data available in microbiome dataset to clinical
clin <- clin %>% full_join(mr %>% dplyr::select(PatID, dataM))

clinl <- clin %>% pivot_longer(cols = !c(PatID,dataM), names_to = 'measure', values_to = 'value')

labs <- c('Visit FEV', 'Baseline FEV', 'One y prior FEV', 'Two y prior FEV', 'IV courses 12m', 
          'IV days 12m')

clinl %>%
  ggplot(aes(x = fct_inorder(measure), y = value, colour = dataM)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(colour = dataM), position = position_dodge2(width = 0.6)) +
  theme_bw(base_size = 12) + xlab('Lung function measure') +
  ggtitle('Microbiome data matching clinical data') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = str_wrap(labs, width = 10)) 
ggsave('MicrobiomevsClin.tiff', dpi = 300, width = 140, height = 80, unit = 'mm')


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

rownames(clin) <- clin$PatID
clin1 <- clin %>% 
  select(BaselineFEV1pp, IVCourses12months, IV.Days.12months) %>%
  drop_na()

pcac <- prcomp(clin1, scale = TRUE)

eigc <- get_eigenvalue(pcac)
scoresc <- pcac$x

autoplot(pcac, loadings = TRUE, loadings.label = TRUE) +
  geom_text_repel(aes(label = rownames(clin1)), size = 3)

###############################  

# Exploratory analysis of reflux data

ref1 <- ref[-1,] 

# exclude patients with missing reflux measures
'%ni%' <- Negate('%in%')
ref2 <- ref1 %>% filter(PatID %ni% c('CF021', 'CF027'))

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
                       gsub('\\.event', '', colnames(ref3))))

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
                           gsub('\\.event', '', ref2l$measure)))

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

ggsave('PCAReflux.tiff', dpi = 300, width = 130, height = 90, unit = 'mm')  
  
# correlation matrix of reflux measures
library(ggcorrplot)

cormat <- cor(ref3, method = 'pearson')
col <- colorRampPalette(c("red", "white", "blue"))(20)
ggcorrplot(cormat, tl.cex = 10, hc.order = TRUE, legend.title = 'r', lab = TRUE)

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
assayT <- t(vsnout[,-1]) 

pca1 <- prcomp(assayT, scale = TRUE)

scores1 <- as.data.frame(pca1$x)                        
loadings1 <- as.data.frame(pca1$rotation)

eigen1 <- get_eigenvalue(pca1)

barplot(eigen1$variance.percent[1:10])

scores1 %>% ggplot(aes(y = PC1, x = PC2)) + geom_point() + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash')

scores1 %>% ggplot(aes(y = PC2, x = PC3)) + geom_point() + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash')

# patients without proteomic profiles
setdiff(refGroups$PatID, rownames(assayT))


# plot PCA results and map to reflux variables/groups
a1 <- autoplot(pca1, #data = refGroups %>% filter(PatID %ni% c("CF011", "CF040", "CF041")), colour = 'refGroup',
         x = 1, y = 2) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  geom_text(aes(x = 0.2, y = 0.73), label = 'CF055', size = 4)

a2 <- autoplot(pca1, data = annot %>% filter(PatID %in% rownames(scores1)),
         colour = 'pc3', x = 3, y = 2) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  scale_color_brewer(palette = 'Set2')+
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
  )))

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

#############################


# Integration of clinical and reflux data (PLS canonical)

clin3 <- clin %>% 
  select(BaselineFEV1pp, IVCourses12months, IV.Days.12months)

rownames(ref) <- ref$PatID
ref3 <- ref[-1,-1]

pls.cr <- pls(X = ref3, Y = clin3, ncomp = 2, scale = TRUE, mode = 'canonical')
pls.cr2 <- pls(X = ref3, Y = clin3, ncomp = 2, scale = TRUE, mode = 'regression')

plotVar(pls.cr, cutoff = 0.7)

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



r2 <- glm(IVDays1 ~ Total.events + Age + Sex + BMI + CFTR.Rx,
         data = join,
         family = Gamma(link = 'identity'))

summary(r2)

# Sparse canonical correlation analysis
library(PMA)
library(ade4)

assay <- assay %>% mutate_if(is.character, as.numeric)

assayt <- t(assay)

setdiff(rownames(ref3), rownames(assayT))
ref4 <- ref3 %>% filter(rownames(ref3) %ni% c("CF011", "CF040", "CF041"))
ref4 <- ref4 %>% mutate_if(is.character, as.numeric)

setdiff(rownames(assayT), rownames(ref3))
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

# canonical PLS
# scales and centers variables by default
pls <- pls(ref4, assayT1, ncomp = 2, mode = 'canonical')
pls

# sparse PLS
# tuning number of components (dimensions)
tune.pls <- pls(assayT1, ref4, ncomp = 4 , mode = 'canonical')
q.tune.pls <- perf(tune.pls, validation = 'Mfold', folds = 10, nrepeat = 5)
plot(q.tune.pls, criterion = 'Q2.total')
ggsave('QcriterionPLS.tiff', dpi = 300, width = 80, height = 80, unit = 'mm')

# tuning number of variables
list.keepX <- c(seq(5, 50, 5))
list.keepY <- c(1:5)

v.tune.pls <- tune.spls(assayT1, ref4, test.keepX = list.keepX,
                        test.keepY = list.keepY, ncomp = 2,
                        nrepeat = 1, folds = 10,
                        mode = 'regression', measure = 'cor')

plot(v.tune.pls)
ggsave('TuneVarPLS.tiff', dpi = 300, width = 200, height = 120, unit = 'mm')


# running optimised model
dev.new()
Spls <- spls(assayT1, ref4, ncomp = 2, keepX = v.tune.pls$choice.keepX, 
             keepY = v.tune.pls$choice.keepY, mode = 'canonical')


# sample plots
pal <- brewer.pal(11, 'Spectral')

joint <- plotIndiv(Spls,
          rep.space = 'XY-variate',
          group = annotC$IVCourses12months,
          title = 'IV antibiotic courses in 12 months',
          legend = TRUE,
          col = pal[11:6],
          legend.title = 'Number')

joint <- plotIndiv(Spls,
                   rep.space = 'XY-variate',
                   group = annotC$BaselineFEV1pp,
                   title = 'Baseline Forced Expiratory Volume',
                   legend = TRUE,
                   legend.title = 'FEV (%)')

joint

pdf('sPLSIndPlotJoint.pdf', width = 8, height = 10)
joint
dev.off()

ref <- plotIndiv(Spls,
                 rep.space = 'Y-variate',
                 group = rownames(ref4),
                 title = 'Reflux data')

ref


prot <- plotIndiv(Spls,
                  rep.space = 'X-variate',
                  group = rownames(ref4),
                  title = 'Protein data')

prot

# variable correlation plot
plotVar(Spls, cutoff = 0.5)
plotVar(Spls) 

Spls$prop_expl_var

# loadings vectors for each component
vip.spls <- vip(Spls)

# variables selected for C1
c1Vars <- as.data.frame(vip.spls[selectVar(Spls, comp = 1)$X$name,1])

perf.Spls <- perf(Spls, validation = 'Mfold', folds = 10, nrepeat = 5)

# proteins
stab.SplsC1 <- perf.Spls$features$stability.X$comp1
stab.SplsC2 <- perf.Spls$features$stability.X$comp2

stab1 <- as.data.frame(stab.SplsC1[1:v.tune.pls$choice.keepX[1]]) 
stab2 <- as.data.frame(stab.SplsC2[1:v.tune.pls$choice.keepX[1]]) 

colnames(stab1) <- 'FreqOccur'
colnames(stab2) <- 'FreqOccur'

stab1 <- stab1 %>% filter(FreqOccurC1 > 0.8)
stab2 <- stab2 %>% filter(FreqOccurC2 > 0.8)

# reflux
stab.SplsC1R <- perf.Spls$features$stability.Y$comp1
stab.SplsC2R <- perf.Spls$features$stability.Y$comp2

stab1R <- as.data.frame(stab.SplsC1R[1:v.tune.pls$choice.keepY[1]]) 
stab2R <- as.data.frame(stab.SplsC2R[1:v.tune.pls$choice.keepY[1]]) 

colnames(stab1R) <- 'FreqOccur'
colnames(stab2R) <- 'FreqOccur'

stab1 <- stab1 %>% filter(FreqOccurC1 > 0.8)
stab2 <- stab2 %>% filter(FreqOccurC2 > 0.8)


#extract stability measure for variables selected by Spls
extr.stab.Spls.C1 <- stab.Spls[selectVar(Spls, comp =1)$X$name]

# create vector of high stability variables for circle correlation plot
Xnames <- as.data.frame(colnames(Spls$X))
colnames(Xnames) <- 'prots'
Ynames <- as.data.frame(colnames(Spls$Y))
colnames(Ynames) <- 'refs'

stabP <- rbind(stab1, stab2) %>%
  mutate(prots = rownames(.),
         VarX = rownames(.))

stabY <- rbind(stab1R, stab2R) %>%
  mutate(refs = rownames(.),
         VarY = rownames(.))

Xnames <- Xnames %>% full_join(stabP)
Ynames <- Ynames %>% full_join(stabY)

listName <- list(Xnames$VarX, Ynames$VarY) 

# variable correlation plot
plotVar(Spls, 
        cutoff = 0.5,
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

# compare with CCA output
loadingsC2 <- loadingsC1 %>%
  left_join(prot %>% dplyr::select(Accession, Description))

###############################


# Sputum microbiome data
colnames(mic)[1] <- 'PatID'
mic <- drop_na(mic) %>% filter(Cohort == 'CF')

# diversity measures and descriptions
div <- mic[,c(1:3, 229:240)]

# relative abundances of bacterial geni
mic1 <- mic[,c(1:228)]

# remove taxa with 0 observations
table(colSums(mic1[4:228]) == 0)

mic1f <- as.data.frame(t(mic1[4:228])) %>% filter(rowSums(.) > 0) %>%
  t() %>% as.data.frame() 

rownames(mic1f) <- mic1$PatID 

# relative abundances of bacterial phyla (higher tax level?)
mic2 <- mic[,c(1:3, 241:ncol(mic))]

# centered log ratio transformation of compositional data
mic1ft <- logratio.transfo(mic1f, logratio = 'CLR', offset = 1)

# principal component analysis
pca <- rda(mic1[,c(4:ncol(mic1))], scale = FALSE)
biplot(pca, scaling = 'symmetric')

micpca <- prcomp(mic1[,c(4:ncol(mic1))], scale = FALSE) 

'%ni%' <- Negate('%in%')
refGroups1 <- refGroups %>% filter(PatID %in% mic[,1]) 
sampleInfo1 <- sampleInfo %>% right_join(mic1[,c(1:2)])

autoplot(micpca, data = sampleInfo1, colour = 'cca')

# samples without microbiome data according to reflux measures
sampleInfo %>% full_join(mic %>% dplyr::select(PatID, Sample)) %>%
  arrange(Total.events) %>%
  ggplot(aes(x = fct_inorder(PatID), y = Total.events, fill = cca)) + geom_bar(stat = 'identity')


