# Reflux project
# Upper respiratory microbiome analysis

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(mixOmics)
library(vegan)
library(ape)
library(factoextra)
library(ggfortify)
library(pheatmap)
library(MatrixCorrelation)
library(plotrix)

#
#
#

ref <- read.csv('Select reflux 01-08-2023.csv')
mic1 <- read.csv('GORCF_MAR_2024_Raw_Count.csv')

ref1 <- ref[-1,]

mic2 <- mic1 %>% filter(PatID %in% ref1$PatID)
mics2 <- mic2 %>% filter(Sample_Type == 'MouthSwirl')

mics2 <- mics2 %>% dplyr::select(!c('OTU_ID', 'Sample_Type_A', 'Sample_Type_B', 'Cohort', 'SampleCohort'))
table(colSums(mics2[3:ncol(mics2)]) == 0)

mics2 <- as.data.frame(t(mics2[3:ncol(mics2)])) %>% filter(rowSums(.) > 0) %>%
  t() %>% as.data.frame() %>% cbind(mics2$PatID) %>% 
  relocate('mics2$PatID')

colnames(mics2)[1] <- 'PatID'

rownames(mics2) <- mics2$PatID
mics2 <- mics2[,-1]

table(mics2[mics2 == 0])/(nrow(mics2)*ncol(mics2)) # 69% observations = 0

# Filtering low-abundant OTUs
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

mic2f <- low.count.removal(mics2, percent = 0.01)
micf2 <- mic2f$data.filter # filtered % all counts

table(micf2[micf2 == 0])/(nrow(micf2)*ncol(micf2)) # 30% observations = 0

# library sizes
mics2 %>% mutate(libSize = rowSums(.),
                 PatID = rownames(.)) %>% 
  left_join(ref1 %>% 
              dplyr::select(PatID, Proximal.total.events)) %>%
  arrange(Proximal.total.events) %>%
  #left_join(mic %>% dplyr::select(PatID, Cohort)) %>%
  ggplot(aes(x = fct_inorder(as.factor(PatID)), y = libSize)) + #, fill = Cohort)) +
  geom_col() +
  #scale_fill_manual(values = c("indianred1","gray"), 
  #                 name = 'Batch', labels = c('Batch1', 'Batch2')) +
  theme_bw(base_size = 12) +
  xlab('Patient ID') +
  theme(axis.text.x = element_text(angle = 45, size = 6),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle('Library sizes in 16S mouth swirl dataset')

ggsave('LibSize_MS.tiff', width = 110, height = 80, dpi = 300, unit = 'mm')

max(rowSums(mics2))/min(rowSums(mics2))

#
#
#

#
# format the microbiome variable names for visualisations
otuNames <- as.data.frame(colnames(micf2)) 
colnames(otuNames) <- 'otu'
otuNames1 <- otuNames %>% separate(otu, c('k', 'T1', 'T2', 'T3', 'T4', 'T5', 'T6'), sep = '__')
otuNames2 <- otuNames1 %>% mutate(T4 = ifelse(is.na(T4) == TRUE, T3, T4)) %>%
  transmute(T4 = gsub('\\..*', '', T4),
            T5 = T5,
            T6 = T6)

otuNames2[is.na(otuNames2) == TRUE] <- ''
otuNames3 <- otuNames2 %>% transmute(otu = paste(T4, T5, T6, sep = '_'))
otuNames4 <- otuNames2 %>% transmute(otu = paste(T5, T6, sep = '_'))


# principal coordinate analysis of compositional data
colnames(micfRA2) <- otuNames3$otu

micfRA2 <- decostand(micf2, method = 'total')

distRA2 <- vegdist(micfRA2, method = 'bray')

pcoa <- pcoa(distRA2)
dev.new()
biplot(pcoa)
plot(pcoa$values$Eigenvalues)

annot <- data.frame(PatID = rownames(micf2)) %>% 
  left_join(ref1 %>% 
              dplyr::select(PatID, Proximal.total.events, Total.events)) %>%
  left_join(micf2 %>% mutate(libSize = rowSums(micf2),
                            PatID = rownames(micf2)) %>%
              dplyr::select(PatID, libSize))

rownames(annot) <- annot$PatID

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

# pCoA with vegan
vpcoa <- capscale(micfRA2 ~ 1, distance = 'bray')
dev.new()
ordiplot(vpcoa, type = 'text')

annot <- cbind(annot, colval)
testcol <- c('#0000FFFF', '#002DD2FF', '#0051AEFF', '#00946BFF','#FFFF00FF')
labels = c('0', '5.6', '10.1', '18.4', '66.3')
colval <- color.scale(annot$Proximal.total.events, c(0,0,1), c(0,1,1), c(1,0,0))

tiff('MS_microbiome_PCoA.tiff', res = 300, width = 146, height = 110, unit = 'mm')
plot(vpcoa, type = 'n', main = 'PCoA of mouth swirl microbiome data')
points(vpcoa, display = "sites", cex = 0.6, pch = 21,
       bg = color.scale(annot$Proximal.total.events, c(0,0,1), c(0,1,1), c(1,0,0)))
text(vpcoa, display = "sp", cex = 0.4, col = "darkred", select = colSums(micfRA2 > 0.5) > 0)
color.legend(-2.75,-1.65,-2.45,-0.75, legend = labels, rect.col = testcol, cex = 0.4, gradient = 'y',
             align = 'lt')
text(x = -2.5, y = -0.65, 'Proximal_Total_Event_No', cex = 0.4, col = 'black')

dev.off()





#
#
#

# transformation of compositional data into Euclidean space
mictf2 <- logratio.transfo(micf2, logratio = 'CLR', offset = 1)
class(mictf2) <- 'matrix'

# PCA of transformed data
pcam <- prcomp(mictf2, center = TRUE, scale = TRUE)

get_eigenvalue(pcam)$variance.percent[1:10]
barplot(get_eigenvalue(pcam)$variance.percent[1:10])

p2 <- autoplot(pcam, data = annot, colour = 'Proximal.total.events', size = 2, x = 1, y = 2) + 
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

ggsave('MS_PCA.tiff', dpi = 300, width = 150, height = 100, unit  ='mm')

# loadings
load <- pcam$rotation

#
#
#


# heatmap of filtered data set at OTU level (relative abundance)
micfRA2 # total count filtered, proportion (RA)

rownames(micfRA2) <- rownames(micf2)
colnames(micfRA2) <- otuNames3$otu

ann_colors = list(
  dataM = c(Batch1 = "indianred1", Batch2 = "grey"))

rownames(annot) <- annot$PatID

micfRA2a <- micfRA2 %>% mutate(PatID = rownames(.)) %>% 
  left_join(annot %>% dplyr::select(PatID, Proximal.total.events)) %>%
  arrange(Proximal.total.events) 
rownames(micfRA2a) <- micfRA2a$PatID
micfRA2a <- micfRA2a %>% dplyr::select(-c(PatID, Proximal.total.events))

tiff('MS_OTU_RA_heatmap.tiff', res = 300, width = 450, height = 280, unit = 'mm')

heatmap <- pheatmap(as.matrix(t(micfRA2a)), annotation_col = annot %>% dplyr::select(Proximal.total.events), 
                    show_rownames = TRUE, cluster_cols =  F, clustering_distance_cols = distRA2,
                    show_colnames = TRUE, cluster_rows = F,
                    fontsize = 10, #annotation_colors = ann_colors,
                    main = 'Relative abundance of OTUs in mouth swirl samples (Genus level)') 

dev.new()
heatmap
dev.off()


#

# similarity between upper and lower respiratory microbiome
mictf2 <- as.data.frame(mictf2) # clr-transformed filtered oral microbiome
micstf <- as.data.frame(micstf) # clr-transformed filtered sputum microbiome

# RV coefficient
mictf2a <- mictf2 %>% filter(rownames(mictf2) %in% rownames(micstf))

rv2 <- RV2(micstf, mictf2a, center = TRUE, impute = FALSE)
rv2 # 0.26

rvadj <- RVadj(micstf, mictf2a, version = 'Maye', center = TRUE)
rvadj # 0.27

# canonical PLS
MMpls <- pls(X = micstf, Y = mictf2a, mode = 'canonical', ncomp = 2) 

plotIndiv(MMpls)
MMpls$prop_expl_var
plotArrow(MMpls)

plotVar(MMpls, cutoff = 0.7)

dev.new()
plotLoadings(MMpls, ndisplay = 10, comp = 1)

cor(MMpls$variates$X, MMpls$variates$Y)

dev.new()
cim(MMpls, cutoff = 0.45, margins = c(30, 30))

#
#
#

# relationship between oral microbiome and reflux
rownames(ref1) <- ref1$PatID
Y = ref1 %>% filter(rownames(.) %in% rownames(mictf2)) %>% dplyr::select(Proximal.total.events)

tune.pls <- pls(mictf2, Y, ncomp = 4 , mode = 'regression')
q.tune.pls <- perf(tune.pls, validation = 'Mfold', folds = 6, nrepeat = 50)
plot(q.tune.pls, criterion = 'Q2.total')

RMOpls <- pls(X = mictf2, Y = Y, mode = 'regression', ncomp = 1)

cor(RMOpls$variates$X, RMOpls$variates$Y)

RMOpls$prop_expl_var

dev.new()
plotLoadings(RMOpls, block = 'X', ndisplay = 10)

# sparse PLS
list.keepX <- c(seq(2, 20, 1))

Y <- Y %>% mutate(Proximal.total.events = ifelse(Proximal.total.events == 0, 0.01, Proximal.total.events))

v.tune.pls <- tune.spls(mictf2, Y, test.keepX = list.keepX,
                        ncomp = 2,
                        validation = 'Mfold', folds = 6, nrepeat = 50,
                        mode = 'regression', measure = 'MSE')

plot(v.tune.pls)

RMOpls1 <- spls(X = mictf2, Y = Y, mode = 'regression', ncomp = 1, keepX = 4)

dev.new()
plotLoadings(RMOpls1, block = 'X')

# variable stability
perf.RMOpls1 <- perf(RMOpls1, validation = 'Mfold', folds = 5, nrepeat = 50)

stab.RMspls <- as.data.frame(perf.RMOpls1$features$stability.X$comp1)

# only Streptoccocus and Ralstonia have good stability

#
#
#

# visualising change in RA of OTUs related with reflux
micfRA2l <- micfRA2 %>% mutate(PatID = rownames(.)) %>%
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events)) %>%
  pivot_longer(cols = !c(PatID, Proximal.total.events), names_to = 'OTU', values_to = 'RA')

micfRA2l <- micfRA2l %>% separate(OTU, c('Order', 'Family', 'Genus'), sep = '_') %>%
  mutate(Family = gsub("\\..*","", Family),
         Order = gsub("\\..*","", Order))#,
         #Class = gsub("\\..*","", Class),
         #Phylum = gsub("\\..*","", Phylum))

micfRA2l <- micfRA2l %>% mutate(
  Family = ifelse(Family == '', NA, Family),
  Genus = ifelse(Genus == '', NA, Genus)
)

# order x axis according to reflux variable
Xorder <- micfRA2l %>% arrange(Proximal.total.events) %>% dplyr::select(PatID)
Xorder <- unique(Xorder)

# plot RA of individual geni
dev.new()
micfRA2l %>% arrange(Proximal.total.events) %>%
  filter(Genus == 'Streptococcus') %>%
  ggplot(aes(x = fct_inorder(as.factor(PatID)), y = RA)) + geom_bar(stat = 'identity')

#
#
#

# lasso clr
library(glmnet)
library(ipflasso)
library(ncvreg)

RMclrlasso <- glmnet(x = mictf2, y = Y$Proximal.total.events, nlambda = 20)


# plot coefficients for each variable for different values of lambda
dev.new()
plot(RMclrlasso, xvar = 'lambda', label = T)

# numerical output
RMclrlasso

# run model with lambda resulting in a desired number of variables
RMclrlasso1 <- glmnet(x = mictf2, y = Y$Proximal.total.events, lambda = 2.898064)

# no of selected variables
RMclrlasso1[['df']]

# names of selected variables
View(as.data.frame(colnames(mictf2[,RMclrlasso1[['beta']]@i])))

# coefficients
x <- predict(RMclrlasso1, type = 'coef')

y <- data.frame(name = row.names(x),
                val = x[,1]) %>% 
  filter(val != 0) %>%
  arrange(desc(abs(val)))

row.names(y) <- NULL


# CV optimisation of model parameters
cvfit <- cv.glmnet(x = as.matrix(mictf2), y = Y$Proximal.total.events,
                   nfolds = 5)
plot(cvfit)

cvfit1 <- cvr.glmnet(X = as.matrix(mictf2), Y = Y$Proximal.total.events,
                     nfolds = 5, ncv = 50, family = 'gaussian', type.measure = 'mse')

cvfit1_res <- data.frame(lambda = cvfit1[['lambda']], MSE = cvfit1[['cvm']])

min <- cvfit1_res[25,]

p1 <- cvfit1_res %>% ggplot(aes(x = log(lambda), y = MSE)) + geom_point() +
  theme_bw() + geom_point(data = min, colour = 'red') + 
  geom_vline(xintercept = log(min$lambda), linetype = 'dashed') +
  theme(plot.title = element_text(hjust = 0.5))

p1

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
cvfit1[['lambda']][19]

# coefficients for min MSE
coefs7 <- as.data.frame(cvfit1[['coeff']][,c('s19')])
colnames(coefs7) <- 'coef'
coefs7 <- coefs7 %>% mutate(name = rownames(.)) %>% relocate()
row.names(coef7) <- NULL
View(coefs7 %>% filter(coef != 0) %>% arrange(desc(abs(coef))))

# Proportion of variance explained
rsq = 1 - (cvfit$cvm[19]/var(Y$Proximal.total.events))

cvfit2 <- cv.ncvreg(X = mictf2, y = Y$Proximal.total.events, penalty = "lasso")

dev.new()
plot(cvfit2, type = "rsq")

#
#
#

# selbal

# CSS normalisation
library(metagenomeSeq)

data.metagenomeSeq = newMRexperiment(t(micf2), # using filtered data
                                     featureData = NULL, libSize = NULL, normFactors = NULL) 
p = cumNormStat(data.metagenomeSeq) # percentile selection - default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p = p)

micf2.CSS = t(MRcounts(data.cumnorm, norm = TRUE, log = TRUE))

library(selbal)

micf2.CSS <- as.data.frame(micf2.CSS) # count data normalised for library size
micf2.CSSo <- micf2.CSS + 1 # offset by 1 to remove 0 values

RMselbal <- selbal(x = micf2.CSSo, y = Y$Proximal.total.events)  

dev.new()
grid.draw(RMselbal$global.plot2)

results_selbal <- selbal_wrapper(result = RMselbal, X = micf.CSSo)

# optimise number of OTUs to select
RMselbal1 <- selbal.cv(x = micf2.CSSo, y = Y$Proximal.total.events, 
                       n.fold = 5, n.iter = 10, covar = NULL)

dev.off()
grid.draw(RMselbal1$global.plot)
View(RMselbal1$cv.tab)

RMselbal1$global.balance
RMselbal1$glm
RMselbal1$accuracy.nvar

#

