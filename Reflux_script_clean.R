## Analysis of sputum proteome and microbiome in cystic fibrosis patients

# author: Aggie Turlo
# project: Cystic fibrosis reflux study
# date: 05/05/2026

#####################

library(tidyverse)
library(factoextra)
library(ggfortify)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(ggcorrplot)
library(cowplot)
library(ggpubr)
library(vsn)
library(vegan)
library(mixOmics)
library(plotrix)
library(pheatmap)
library(MatrixCorrelation)
library(dendextend)

setwd("~/Cystic_fibrosis_R_Lord")

# load data
ref <- read.csv('Select reflux 01-08-2023.csv')
prot <- read.csv('Progenesis_CF_HV_unique_proteins.csv')
mic1 <- read.csv('GORCF_MAR_2024_Raw_Count.csv')
clin <- read.csv('GOR_ClinicalEndpoint.csv')
clin1 <- read.csv('GOR-CF database 07-03-2024.csv')

#

## Principal Component Analysis of the reflux measures
ref1 <- ref[-1,] # remove control values row
ref2 <- ref1 %>% drop_na() # remove patients with majority of observations missing
rownames(ref2) <- ref2$PatID
colnames(ref2) <- gsub('\\.episode', '', 
                       gsub('\\.exposure', '', 
                            gsub('\\.events', '', colnames(ref2))))

ref2l <- ref2 %>% pivot_longer(cols = !PatID, names_to = 'measure', values_to = 'value')
ref2l <- ref2l %>% mutate(category = 
                            ifelse(str_detect(
                              measure, 'events') == TRUE, 'event', 
                              ifelse(str_detect(
                                measure, 'exposure') == TRUE, 'exposure', 
                                ifelse(str_detect(measure, 'episode') == TRUE, 
                                       'episode', 'exposure')
                                )))
ref3 <- ref2[,-1]

# run PCA
pca <- prcomp(as.matrix(ref3), scale = TRUE)
# export scores and loadings
scores <- as.data.frame(pca$x)                        
loadings <- as.data.frame(pca$rotation)

# PCA score plot
# show only 3 most important loadings
l <- c(rep(NA,4),'blue',rep(NA,2),'blue',rep(NA,3),'blue',rep(NA,6))
p <- autoplot(pca, 
              loadings = TRUE, 
              loadings.label = TRUE, 
              loadings.label.size = 3, 
              loadings.label.repel = TRUE, 
              loadings.colour = l,
              loadings.label.colour = l,
              loadings.cutoff = 0.2) + 
  geom_vline(xintercept = 0, linetype = 'longdash') + 
  geom_hline(yintercept = 0, linetype = 'longdash') +
  theme_bw() +
  ggtitle('PCA biplot of reflux measures') +
  theme(plot.title = element_text(hjust =0.5))

p

# PCA loadings plots
p1 <- loadings %>% arrange(PC1) %>% ggplot(aes(x = fct_inorder(rownames(.)), y = PC1)) + 
  geom_bar(stat = 'identity') + coord_flip() + theme_bw(base_size = 8) + 
  theme(axis.title.y = element_blank())

p2 <- loadings %>% arrange(PC2) %>% ggplot(aes(x = fct_inorder(rownames(.)), y = PC2)) + 
  geom_bar(stat = 'identity') + coord_flip() + theme_bw(base_size = 8) +
  theme(axis.title.y = element_blank()) + 
  xlab('PC2 loadings')

g <- arrangeGrob(p1, p2, nrow = 1, ncol = 2, 
                 top = grid::textGrob('Principal Component Loadings', 
                                      gp=grid::gpar(fontsize=10)))
plot(g)

# correlation matrix of reflux measures
cormat <- cor(ref3, method = 'pearson')
col <- colorRampPalette(c("red", "white", "blue"))(20)
corplot <- ggcorrplot(cormat, tl.cex = 7, hc.order = TRUE, legend.title = 'r', lab = TRUE,
                      lab_size = 1, type = 'lower') + 
  ggtitle('Pearsons correlation between reflux parameters') +
  theme(plot.title = element_text(hjust = 0.5, size = 10),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, 'line'))
corplot

# Assemble Supplementary Figure 1
pdf(file = 'Figure_S1.pdf', height = 5)
ggdraw() +
  draw_plot(p, x = 0, y = 0.5, width = 0.4, height = 0.45) +
  draw_plot(g, x = 0, y = 0, width = 0.5, height = 0.47) +
  draw_plot(corplot, x = 0.45, y = 0.2, width = 0.55, height = 0.8) +
  draw_plot_label(label = c('A', 'B', 'C'), size = 10, 
                  x = c(0, 0, 0.45), y = c(1, 0.5, 0.95))
dev.off()

#
#
#

# Modelling the effect of reflux on clinical outcomes
# create consistent patient id
colnames(clin)[1] <- 'PatID'
clin1$Assignment.Number <- gsub('GOR', 'CF', clin1$Assignment.Number)
colnames(clin1)[1] <- 'PatID'

# create joint clinical-reflux database with relevant variables
data <- clin %>% 
  dplyr::select(PatID, BaselineFEV1pp, IV.Days.12months, 
                X1y_priorFEV1., X2y_prior_FEV1.) %>%
  left_join(ref1 %>% 
              dplyr::select(PatID, Total.events, Proximal.total.events, 
                            Demeester.score, Total.acid.exposure)) %>%
  left_join(clin1 %>% 
              dplyr::select(PatID, Age, BMI, Pseudomonas.colonisation,
                            Pancreatic.insufficency, Diabetes, CFTR.Modulator))

data$BMI <- as.numeric(data$BMI)

# shorten column names
colnames(data) <- str_split(colnames(data), "\\.", simplify=T)[,1]
colnames(data) <- as.character(colnames(data))
colnames(data)[9] <- 'Total.acid'
data <- data %>% mutate(Pseudomonas = 
                          ifelse(Pseudomonas == 'TBC', NA, Pseudomonas))

# modelling FEV1 % predicted
fev1 <- lm(BaselineFEV1pp ~ Proximal + Age + BMI + Diabetes + Pseudomonas,
           data = data)

summary(fev1)

# modelling number of days with IV antibiotic treatment
iv <- glm.nb(IV ~ Proximal + Pseudomonas + Diabetes + BaselineFEV1pp,
             data = data)

summary(iv)

#
#
#

## Principal Component Analysis of protein data
colnames(prot) <- prot[2,]
prot <- prot[-c(1:2),]
assay <- prot %>% dplyr::select(contains('CF'))
rownames(assay) <- prot$Accession  
assay <- assay %>% mutate_if(is.character, as.numeric)

# variance stabilisation
vsn <- vsn::vsnMatrix(as.matrix(assay))
vsnout <- vsn@hx
vsnout <- as.data.frame(vsnout) %>% mutate(prot = rownames(vsnout)) %>% 
  relocate(prot)
vsnl <- vsnout %>% pivot_longer(cols = !prot, 
                                names_to = 'PatID', 
                                values_to = 'intensity')
assayT <- t(vsnout[,-1])

# run PCA
pca1 <- prcomp(assayT, scale = TRUE)
# export scores
scores1 <- as.data.frame(pca1$x)                        
scores2 <- scores1 %>% mutate(PatID = rownames(.)) %>%
  relocate(PatID)

write.csv(scores2, 'Proteome_PCA_scores.csv')
# to be included as Supplemental File 1

# export loadings
loadings1 <- as.data.frame(pca1$rotation)
loadings2 <- loadings1 %>% mutate(Accession = rownames(.)) %>%
  left_join(prot %>% dplyr::select(Accession, Description)) %>%
  separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  mutate(Description = gsub('OS.*', '', Description)) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(!c(PE, SV)) %>%
  relocate(Accession, Description, Gene)

write.csv(loadings2, 'Proteome_PCA_loadings.csv')
# to be included as Supplemental File 1

# export variance explaines
eigen1 <- get_eigenvalue(pca1)

write.csv(eigen1, 'Proteome_PCA_eigenvalues.csv')
# to be included as Supplemental File 1

# visualise PCA scores
an <- data.frame(PatID = rownames(assayT)) %>% 
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))
rownames(an) <- an$PatID

# create quantile colour scale 
testcol <- c('#0000FFFF', '#009E61FF','#FFFF00FF')
colvalquant <- base::seq(from = 0, to = 1, length.out = 5)
quants <- quantile(an$Proximal.total.events, colvalquant)
an$ptile_var <- ecdf(an$Proximal.total.events)(an$Proximal.total.events)

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
  theme(legend.text = element_text(size = 6), 
        legend.title = element_text(size = 8), 
        legend.key.height = unit(5, 'mm'),
        legend.key.width = unit(3, 'mm'),
        axis.title.x = element_blank())

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
  theme(legend.text = element_text(size = 6), 
        legend.title = element_text(size = 8), 
        legend.key.height = unit(5, 'mm'),
        legend.key.width = unit(3, 'mm'),
        axis.title.x = element_blank())

protPC <- ggarrange(pc12, pc23, common.legend = TRUE, legend = 'right', nrow = 2)
protPC <- annotate_figure(protPC, top = text_grob("PCA plots of sputum proteome data", 
                                                  color = "black", face = "bold", size = 8.5),
                          bottom = text_grob('PC2 (10.9%)', size = 9))

protPC

#
#
#

## Exploratory analysis of microbiome data
# microbiome data from patients with reflux data available
mic1 <- mic1 %>% filter(PatID %in% ref1$PatID)
# keep sputum samples only
mics1 <- mic1 %>% filter(Sample_Type == 'Sputum')
# remove taxa with 0 observations
mics1 <- mics1 %>% 
  dplyr::select(!c('OTU_ID', 'Sample_Type_A', 'Sample_Type_B', 
                   'Cohort', 'SampleCohort'))

table(colSums(mics1[3:ncol(mics1)]) == 0)

mics1 <- as.data.frame(t(mics1[3:ncol(mics1)])) %>% 
  filter(rowSums(.) > 0) %>%
  t() %>% as.data.frame() %>% 
  cbind(mics1 %>% dplyr::select(PatID)) %>%
  relocate(PatID)

rownames(mics1) <- mics1$PatID
mics1 <- mics1[,-1] # unfiltered non-normalised OTU table

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

mics1f <- low.count.removal(mics1, percent = 0.01)
micf <- mics1f$data.filter # filtered based on % of all counts

table(micf[micf == 0]) # 38% of data = 0
table(mics1[mics1 == 0]) # 71% of data = 0

#

# format OTU names for visualisations
otuNames <- as.data.frame(colnames(micf)) 
colnames(otuNames) <- 'otu'
otuNames[35,1] <- 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Mogibacteriaceae.g__'
otuNames[36,1] <- 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Mogibacteriaceae.g__Mogibacterium'
otuNames[10,1] <- 'k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Paraprevotellaceae.g__Prevotella'
otuNames[37,1] <- 'k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Tissierellaceae.g__Parvimonas'

otuNames1 <- otuNames %>% separate(otu, c('k', 'T1', 'T2', 'T3', 'T4', 'T5'), sep = '\\.')

otuNames2 <- otuNames1 %>% mutate(T3 = gsub('__', '.', (T3 = ifelse(T3 == 'Other', T2, T3))),
                                  T4 = gsub('__', '.', (T4 = ifelse(T4 == 'Other', T3, T4))),
                                  T5 = gsub('__', '.', (T5 = ifelse(T5 %in% c('g__',  'Other'), 
                                                                    paste(T4, '.uncl.', sep = ''), T5)))) 

otuNames3 <- otuNames2 %>% transmute(otu = paste(T4, T5, sep = '_'))
otuNames4 <- otuNames2 %>% transmute(otu = T5)

#

# relative abundance transformation
micfRA <- decostand(micf, method = 'total')

#

# Centered Log Ratio Transformation with offset of 1
micst <- logratio.transfo(mics1, logratio = 'CLR', offset = 1)
micstf <- logratio.transfo(micf, logratio = 'CLR', offset = 1)

class(micst) <- 'matrix'
class(micstf) <- 'matrix'

#

# Principal Coordinate analysis
colnames(micfRA) <- otuNames4$otu
vpcoa <- capscale(micfRA ~ 1, distance = 'bray')

annot <- data.frame(PatID = rownames(micf)) %>% 
  left_join(ref1 %>% dplyr::select(Proximal.total.events, PatID))
annot <- annot %>% arrange(Proximal.total.events)
rownames(annot) <- annot$PatID
annot$ecdf <- 1:length(annot$Proximal.total.events)/length(annot$Proximal.total.events)
annot$ecdf[1] <- 0
annot <- annot[order(match(annot$PatID, rownames(micfRA))),]

colval <- color.scale(annot$ecdf, c(0,0,1), c(0,1,1), c(1,0,0))
names(colval) <- annot$PatID
annot <- cbind(annot, colval)

testcol <- c('#0000FFFF', '#008877FF', '#00FF00FF', '#88FF00FF','#FFFF00FF')
labels = c('0', '5.5', '10.2', '19.7', '66.3')

dev.new()
par(mar = c(3.9, 3.6, 1, 1))
plot(vpcoa, type = 'n', main = 'PCoA plot of sputum microbiome data', cex.main = 0.75, 
     cex.lab = 0.6, cex.axis = 0.6, xlab = 'MDS1 (42.7%)', ylab = 'MDS2 (14.1%)')
points(vpcoa, display = "sites", cex = 0.6, pch = 21, 
       bg = colval)
text(vpcoa, display = "sp", cex = 0.5, col = "darkred", select = colSums(micfRA > 0.5) > 0)
color.legend(-2.5, 1.0, -2.2, 2.0, legend = labels, rect.col = testcol, cex = 0.5, gradient = 'y',
             align = 'lt')
text(x = -2, y = 2.2, 'Proximal reflux events', cex = 0.6, col = 'black')

res <- recordPlot()
res

#

# hierarchical clustering
micfRA <- micfRA*100
distRA <- vegdist(micfRA, method = 'bray')

annot <- data.frame(PatID = rownames(micf)) %>% 
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))
rownames(annot) <- annot$PatID

annot <- annot %>% arrange(Proximal.total.events)
annot$ecdf <- 1:length(annot$Proximal.total.events)/length(annot$Proximal.total.events)
annot$ecdf[1] <- 0
colval <- color.scale(annot$ecdf, c(0,0,1), c(0,1,1), c(1,0,0))
annot <- cbind(annot, colval)

cols <- annot$colval
names(cols) <- annot$Proximal.total.events
colsl <- list(ecdf = cols)

testcol <- c('#0000FFFF', '#008877FF', '#00FF00FF', '#88FF00FF','#FFFF00FF')
labels = c('0', '5.5', '10.2', '19.7', '66.3')
leg <- data.frame(labels, testcol, x = 1, y = 1:5)
leg$labels <- as.factor(leg$labels)

heatmap <- pheatmap(as.matrix(t(micfRA)), 
                    annotation_col = annot %>% dplyr::select(ecdf), 
                    show_rownames = TRUE, 
                    cluster_cols =  T, 
                    clustering_distance_cols = distRA,
                    show_colnames = FALSE, 
                    cluster_rows = F,
                    fontsize = 5.5,
                    annotation_colors = colsl,
                    cellheight = 5, cellwidth = 5,
                    annotation_legend = FALSE,
                    annotation_names_col = FALSE) 

heatmap
dev.off()

# create figure legend
legGrob <- leg %>% ggplot(aes(x = fct_inorder(as.factor(labels)), y = x, fill = labels)) + geom_tile() +
  scale_fill_manual(values = testcol, breaks = labels, name = 'Proximal \n total events') +
  theme_classic(base_size = 8) +
  ggtitle('Proximal reflux events') +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 8))

legGrob

hm <- as_ggplot(heatmap$gtable)
hm


#

# Assemble Figure 2
protPC # 2 PCA scores plots proteome
res # PCoA scores plot microbiome
hm # RA heatmap microbiome

pdf(file = 'Figure_2_v03.pdf', width = 7, height = 5.8)
ggdraw() +  
  draw_plot(protPC, x = 0.025, y = 0.5, width = 0.4, height = 0.48) +
  draw_plot(res, x = 0, y = 0, width = 0.43, height = 0.4) +
  draw_plot(hm, x = 0.43, y = 0, width = 0.55, height = 1) +
  draw_plot(legGrob, x = 0.76, y = 0.84, width = 0.2, height = 0.08) +
  draw_plot_label(label = c('A','B','C'), size = 10, 
                  x = c(0, 0, 0.42), y = c(1, 0.43, 1)) +
  draw_plot_label(label = 'Heatmap of bacterial relative abundance in sputum',
                  size = 9, x = 0.25, y = 0.98) +
  draw_plot_label(label = 'RA', x = 0.905, y = 0.83, size = 8)
dev.off()

#
#
#

## Multi-block analysis
# multiblock PLS - reflux, proteomic and microbiome data
# create reflux outcome variable
r1 <- ref1 %>%
  dplyr::select(Proximal.total.events) %>%
  mutate(Proximal.total.events = ifelse(Proximal.total.events == 0, 
                                        0.01, Proximal.total.events))
rownames(r1) <- ref1$PatID

# format protein dataset
assayT <- assayT %>% as.data.frame() %>% 
  mutate(PatID = rownames(.)) %>%
  arrange(PatID) 

# fromat microbiome dataset
colnames(micstf) <- otuNames$otu
micstf <- as.data.frame(micstf)

micstf # clr transformed total filtered microbiome data
r1 # reflux proximal total events
assayT # sputum proteome

# identify complete observations across 3 datasets
tr <- intersect(
  intersect(rownames(r1), rownames(assayT)),
  rownames(micstf))

# create training dataset (n = 27)
# reflux
rtr <- r1 %>% mutate(PatID = rownames(.)) %>%
  filter(PatID %in% tr) %>%
  arrange(PatID)

rownames(rtr) <- rtr$PatID
rtr <- rtr %>% 
  subset(select = -(PatID))

# protein
protr <- assayT %>% filter(PatID %in% tr) %>%
  arrange(PatID)
rownames(protr) <- protr$PatID
protr <- protr %>% 
  subset(select = -c(PatID))

# microbiome
mictr <- as.data.frame(micstf) %>% mutate(PatID = rownames(.)) %>%
  filter(PatID %in% tr) %>%
  arrange(PatID)
rownames(mictr) <- mictr$PatID
mictr <- mictr %>% subset(select = -(PatID))

#

# calculate RV coefficient and adjusted RV coefficient
rv2 <- RV2(mictr, protr, center = TRUE, impute = FALSE)
rv2

rvadj <- RVadj(mictr, protr, version = 'Maye', center = TRUE)
rvadj

# code datasets for multiblock analysis
X <- list(prot = protr, mic = mictr) # predictor datasets
Y <- as.matrix(rtr) # outcome dataset (variable)

# create design matrix
design <- matrix(0.1, ncol = length(X), nrow = length(X),
                 dimnames = list(names(X), names(X)))
diag(design) <- 0

# run multiblock PLS with
RMPpls <- block.pls(X, Y, design = design, ncomp = 2, mode = 'regression')
RMPpls$design
RMPpls$prop_expl_var

# sample scores
Pscores <- RMPpls$variates$prot %>% as.data.frame() %>% 
  mutate(block = 'Proteome', PatID = rownames(.))
Mscores <- RMPpls$variates$mic %>% as.data.frame() %>% 
  mutate(block = 'Microbiome', PatID = rownames(.))
Scores <- rbind(Pscores, Mscores) %>% 
  left_join(ref1 %>% dplyr::select(PatID, Proximal.total.events))

Scores1 <- Scores %>% dplyr::select(!c(Proximal.total.events))
write.csv(Scores1, 'Bpls_scores.csv') # save all the sample scores
# included in Supplemental File 2

# visualise score plots
colvalquant <- base::seq(from = 0, to = 1, length.out = 5)
quants <- quantile(Scores$Proximal.total.events, colvalquant)
Scores$ptile_var <- ecdf(Scores$Proximal.total.events)(Scores$Proximal.total.events)

pls_plot <- Scores %>% ggplot(aes(x = comp1, y = comp2)) + 
  geom_point(size = 1.5, pch = 21, aes(fill = ptile_var)) + 
  theme_bw(base_size = 9) +
  facet_wrap(~ block, scales = 'free', ncol = 1) +
  labs(title = 'Block PLS scores plots') +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient2(low = '#0000FFFF',
                       mid = '#009E61FF',
                       high = '#FFFF00FF',
                       midpoint = 0.5,
                       labels = round(quants,1)) + 
  labs(fill = 'Proximal reflux events') +
  theme(legend.text = element_text(size = 8), 
        legend.title = element_text(size = 8), 
        legend.key.height = unit(3, 'mm'),
        legend.key.width = unit(5, 'mm'),
        legend.position = 'bottom',
        plot.title = element_text(size = 10))


pls_plot

# extract model loadings
Ploadings <- RMPpls$loadings$prot %>% as.data.frame() %>%
  mutate(Accession = rownames(.)) %>% 
  left_join(prot %>% dplyr::select(Accession, Description)) %>%
  separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  #dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(-c(PE, SV)) %>%
  mutate(Description = gsub('OS.*', '', Description),
         Accession = gsub(';.*', '', Accession)) 

write.csv(Ploadings, 'Bpls_protein_loadings.csv')
# included in Supplemental File 2

Mloadings <- RMPpls$loadings$mic %>% as.data.frame() %>%
  mutate(otu = otuNames4$otu)

colnames(Mloadings)[3] <- 'OTU'

Loadings_C1 <- rbind(Mloadings %>% 
                       mutate(block = 'Microbiome',
                              var = OTU) %>%
                       arrange(desc(abs(comp1))) %>% dplyr::slice(1:10) %>%
                       arrange(desc(comp1)) %>%
                       dplyr::select(comp1, var, block), 
                     Ploadings %>% mutate(block = 'Proteome',
                                          var = Description) %>%
                       arrange(desc(abs(comp1))) %>% dplyr::slice(1:10) %>%
                       arrange(desc(comp1)) %>%
                       dplyr::select(comp1, var, block))


pc1_plot <- Loadings_C1 %>%
  ggplot(aes(x = fct_inorder(var), y = comp1)) + 
  facet_wrap(~ block, scales = 'free', ncol = 1) +
  theme_bw(base_size = 9) +
  geom_col(width = 0.8) + coord_flip() + 
  theme_bw() + theme(axis.title.y = element_blank(),
                     axis.title.x = element_blank(),
                     axis.text.y = element_text(size = 8),
                     axis.text.x = element_text(size = 8),
                     plot.title = element_text(size = 10)) +
  ggtitle('PLS component 1 loadings') +
  theme(plot.title = element_text(hjust = 0))

pc1_plot

#

# Assemble Figure 3
pdf('Figure_3_v02.pdf', width = 7, height = 4.5)
ggdraw() +
  draw_plot(pls_plot, x = 0.03, y = 0, width = 0.35, height = 0.93) +
  draw_plot(pc1_plot, x = 0.43, y = 0.14, width = 0.52, height = 0.8) +
  draw_plot_label(label = c('A','B'), x = c(0, 0.6), y = c(1, 1))
dev.off()

#

# agreement across datasets
cor(Pscores$comp1, Mscores$comp1)
cor(Pscores$comp2, Mscores$comp2)

cor(Pscores$comp2, Y)
cor(Mscores$comp2, Y)

proc1 <- vegan::procrustes(X = as.matrix(Pscores[,1:2]), Y = as.matrix(Mscores[,1:2]),
                           symmetric = TRUE)


proc1

#

# correlation between individual variables and BPLS components
corC <- plotVar(RMPpls)
corC1 <-  corC %>% mutate(x = round(x, digits = 2),
                          y = round(y, digits = 2))

write.csv(corC1, 'BPLS_CorrelationCoefs.csv')
# included in Supplemental File 2

# filter variables with correlation coefficient > 0.7
corCf <- corC1 %>% filter(abs(x) >= 0.7) %>%
  left_join(prot %>% dplyr::select(Accession, Description), 
            by = c('names' = 'Accession'))

corCf <- corCf %>% 
  separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(-c(PE, SV)) %>%
  mutate(Gene = ifelse(is.na(Gene) == TRUE, names, Gene),
         CorRef = ifelse(x > 0, 'pos', 'neg')) 

corCf$names[24] <- 'k__Bacteria.p__Firmicutes.c__Bacilli.o__Lactobacillales.f__Streptococcaceae.g__'

# heatmap of z-score between variables and the outcome (reflux)
conc <- cbind(protr, mictr)
concf <- conc[,c(corCf$names)]
concfs <- as.data.frame(scale(concf, center = TRUE, scale = TRUE))
colnames(concfs) <- corCf$Gene
concfs <- concfs[match(rownames(as.data.frame(Y) %>% arrange(Proximal.total.events)), rownames(concfs)),]
colnames(concfs)[24] <- 'f.Streptococcaceae'

rows <- data.frame(names = colnames(concf)) %>% 
  left_join(corCf %>% dplyr::select(names, Block, Gene, CorRef))
rownames(rows) <- rows$Gene
rows <- rows %>% dplyr::select(CorRef)
rownames(rows)[24] <- 'f.Streptococcaceae'

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

testcol <- c('#0000FFFF', '#008976FF', '#00FF00FF', '#76FF00FF','#FFFF00FF')
labels = c('0', '5.5', '10.2', '19.7', '66.3')
leg <- data.frame(labels, testcol, x = 1, y = 1:5)
leg$labels <- as.factor(leg$labels)

concfsP <- concfs[,-24]  
concfsM <- as.data.frame(concfs[,24])
rownames(concfsM) <- rownames(concfs)
colnames(concfsM) <- 'f.Streptococcaceae'

breaksList = seq(-3.7, 3.1, by = 0.1)

Bpls_hm <- pheatmap(t(concfs), annotation_col = annot %>% dplyr::select(ecdf),
                    clustering_method = 'ward.D2',
                    cluster_cols = T, 
                    annotation_row = rows, cluster_rows = T,
                    show_colnames = T,
                    annotation_colors = colsl,
                    treeheight_col = 50,
                    treeheight_row = 10,
                    annotation_legend = FALSE,
                    annotation_names_col = FALSE,
                    annotation_names_row = FALSE,
                    fontsize = 5,
                    fontsize_row = 6,
                    color = colorRampPalette(
                      rev(
                        brewer.pal(
                          n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
                    breaks = breaksList) # Sets the breaks of the color scale as in breaksList)

Bpls_hm1 <- as_grob(Bpls_hm$gtable)

#

# create figure legend
legGrob <- leg %>% ggplot(aes(x = fct_inorder(as.factor(labels)), y = x, fill = labels)) + 
  geom_tile() +
  scale_fill_manual(values = testcol, breaks = labels) +
  theme_classic(base_size = 8) +
  ggtitle('Proximal reflux events') +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 6, vjust = -3)) +
  scale_y_discrete(position = 'right')

legGrob

legGrob2 <- data.frame(group = 1, Corr = c('Pos', 'Neg')) %>% 
  ggplot(aes(x = Corr, y = group, fill = Corr)) + 
  geom_tile() +
  scale_fill_manual(values = c('#00BFC4','#F8766D')) +
  theme_classic(base_size = 8) +
  ggtitle('Correlation with C1') +
  theme(axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none',
        plot.title = element_text(size = 6, vjust = -3, hjust = 1.1)) +
  scale_y_discrete(position = 'right')


# Assemble Figure 4
pdf(file = 'Figure_4_v02.pdf', width = 7, height = 5)
ggdraw() +  
  draw_plot(Bpls_hm1, x = 0.15, y = 0.04, width = 0.8, height = 0.8) +
  draw_plot(legGrob, x = 0.76, y = 0.82, width = 0.18, height = 0.11) +
  draw_plot(legGrob2, x = 0.86, y = 0.72, width = 0.08, height = 0.11) +
  draw_plot_label(label = 'z-score', x = 0.86, y = 0.715, size = 7) +
  draw_plot_label(label = 'Variables most strongly correlated with PLS component 1', 
                  x = 0, y = 0.97, size = 10) 
dev.off()

#

# Enrichment analysis of proteins most correlated with outcome variable
# subset into positively and negatively correlated
corCf1 <- corCf %>% filter(Block == 'prot')
corN <- corCf1 %>% filter(x < 0)
corP <- corCf1 %>% filter(x > 0)

posuniprot <- corP$names
posuniprot <- gsub(';.*','', posuniprot) # keep the first UniProt accession

neguniprot <- corN$names
neguniprot <- gsub(';.*','', neguniprot) 

# export for use with online enrichment tool https://tools.dice-database.org/GOnet/
write.csv(posuniprot, 'POSuniprot.csv')
write.csv(neguniprot, 'NEGuniprot.csv')

## create background database of gene symbols
bguniprot <- prot$Accession
bguniprot <- gsub(';.*','', bguniprot)
write.csv(bguniprot, 'BGUniprot.csv') # export for use in GOnet

#
#
#

# repeat hierarchical clustering of variable most associated with reflux
distC <- vegdist(concfs, method = 'euclidean')
clustC <- hclust(distC, method = 'ward.D2')

dend <- as.dendrogram(hclust(distC, method = 'ward.D2'))
plot(dend)
get_leaves_nodePar(dend)[[1]]
par(mfrow = c(1,1), mar = c(4,4,4,4))
dend <- set(dend, "labels_cex", 0.5)
get_leaves_nodePar(dend)[[1]]
dend <- color_branches(dend, k = 2, col = brewer.pal(n = 2, name = 'Set2'))
plot(dend, main = 'Hierarchical clustering of sputum samples', cex.main = 0.8)
dend1 <- recordPlot()

# cut the tree and assign cluster membership
clustm <- cutree(clustC, h = 15)
clustM <- data.frame(PatID = clustC[['labels']], cluster = clustm)
clustM1 <- clustM %>% left_join(ref1) %>% left_join(clin)

# boxplots of reflux events according to clusters (w/o outlier)
clustM1l <- clustM1 %>% 
  dplyr::select(PatID, cluster, Proximal.total.events, 
                BaselineFEV1pp, IV.Days.12months) %>%
  pivot_longer(cols =! c(PatID, cluster), 
               names_to = 'parameter', 
               values_to = 'value') %>%
  mutate(parameter = ifelse(
    parameter == 'Proximal.total.events', 'Proximal reflux events',
    ifelse(parameter == 'BaselineFEV1pp', 'FEV1 (% predicted)',
           'IV days (number)')))

clust_boxplots <- clustM1l %>% 
  ggplot(aes(x = as.factor(cluster), y = value, fill = as.factor(cluster))) +
  facet_wrap(~ factor(parameter, 
                      levels = c('Proximal reflux events',
                                 'FEV1 (% predicted)',
                                 'IV days (number)')), 
             nrow = 3, scale = 'free') +
  geom_boxplot(outliers = FALSE) + 
  geom_jitter(size = 0.8) +
  scale_fill_brewer(palette = 'Set2') +
  theme_bw(base_size = 9) +
  labs(fill = 'cluster') + xlab('cluster') +
  theme(legend.position = 'none',
        axis.title.y = element_blank())

clust_boxplots

# Assemble Supplementary Figure 2
pdf('Figure_S2_v02.pdf', width = 7, height = 5)
ggdraw() + 
  draw_plot(dend1, x = 0, y = 0.2, width = 0.6, height = 0.7) +
  draw_plot(clust_boxplots, x = 0.53, y = 0.1, width = 0.23, height = 0.8) +
  draw_plot_label(label = c('A', 'B'), x = c(0, 0.5), y = c(0.95, 0.95)) +
  draw_plot_label(label = '*', x = 0.608, y = 0.8)
dev.off()  


#
#
#

