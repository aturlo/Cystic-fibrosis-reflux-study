# univariate modelling of proteins selected by sPLS on reflux and sputum proteome data
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(lme4)
library(lmerTest)
library(stringr)
library(gridExtra)
library(effects)

loadingsPLS <- read.csv('C:/Users/q29634at/OneDrive - The University of Manchester/Documents/Cystic_fibrosis_R_Lord/sPLS_canonical/SplsProtLoadings.csv')
stabPLS1 <- read.csv('C:/Users/q29634at/OneDrive - The University of Manchester/Documents/Cystic_fibrosis_R_Lord/sPLS_canonical/SplsProtC1Stability.csv')
stabPLS2 <- read.csv('C:/Users/q29634at/OneDrive - The University of Manchester/Documents/Cystic_fibrosis_R_Lord/sPLS_canonical/SplsProtC2Stability.csv')
assayT <- read.csv('C:/Users/q29634at/OneDrive - The University of Manchester/Documents/Cystic_fibrosis_R_Lord/ProteinVSNTransformed.csv')
ref <- read.csv('C:/Users/q29634at/OneDrive - The University of Manchester/Documents/Cystic_fibrosis_R_Lord/Select reflux 01-08-2023.csv')
prot <- read.csv('C:/Users/q29634at/OneDrive - The University of Manchester/Documents/Cystic_fibrosis_R_Lord/Progenesis_CF_HV_unique_proteins.csv')


ref <- ref[-1,]

colnames(prot) <- prot[2,]
prot <- prot[-c(1:2),]
prot$Accession <- gsub(';', '.', prot$Accession)

colnames(stabPLS1)[1] <- 'Accession'
colnames(stabPLS2)[1] <- 'Accession'
colnames(assayT)[1] <- 'PatID'

loadingsPLS <- loadingsPLS %>% left_join(stabPLS1) # component 1
loadingsPLS <- loadingsPLS %>% left_join(stabPLS2, by = 'Accession') # component 2

loadingsPLS <- loadingsPLS %>% arrange(desc(abs(comp1))) # component 1
loadingsPLS <- loadingsPLS %>% arrange(desc(abs(comp2))) # component 2

loadingsPLS <- loadingsPLS %>% filter(FreqOccur > 0.6)

top_prots <- loadingsPLS[1:10,]
top_prots <- top_prots$Accession
top_prots <- gsub(';', '.', top_prots)

# for pls with proximal reflux measure only
top_genes <- c('IGHV5-51', 'BPIFB1', 'IGKV4-1', 'ALDOA', 'EEF1D', 
               'ANPEP', 'GSN', 'IGLV3-9', 'LTF', 'CD63')
top_prots <- c('A0A0C4DH38', 'Q8TDL5', 'P06312', 'P04075.P05062', 'P29692', 
               'P15144', 'P06396', 'A0A075B6K5.A0A075B6K2', 'P02788', 'P08962')

top_prots1 <- as.data.frame(assayT[,top_prots])
rownames(top_prots1) <- assayT$PatID


join <- ref %>% dplyr::select(PatID, Total.events, Proximal.total.events, Total.acid.exposure) %>% left_join(top_prots1 %>% mutate(PatID = rownames(.)))
joinl <- join %>% pivot_longer(cols = !c(PatID, Total.events, Proximal.total.events, Total.acid.exposure), names_to = 'prot', values_to = 'intensity')

joinl <- joinl %>% left_join(prot %>% dplyr::select(Accession, Description), by = c('prot' = 'Accession'))
joinl <- joinl %>% separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ')

plots <- joinl %>% group_by(Gene) %>%
  do(plots = ggplot(data = ., aes(x = intensity, y = Total.acid.exposure)) + 
       geom_point() + ggtitle(.$Gene) + theme_bw())

marrangeGrob(plots$plots, nrow = 4, ncol = 3)

hist(join$Total.events)
hist(join$Proximal.total.events)
hist(join$Total.acid.exposure)

joinl$Proximal.total.events[joinl$Proximal.total.events == 0] <- 0.01

gamreg <- function(Prot){
subset <- joinl %>% filter(prot == Prot) %>% drop_na(intensity)

mod <- glm(Total.acid.exposure ~ intensity,
           data = subset,
           family = Gamma(link = 'log'))

as.data.frame(coef(summary(mod)))

}

gamreg('Q8TDL5')

regRes <- lapply(top_prots, gamreg)
regResdf <- bind_rows(regRes) %>%
  mutate(coef = rownames(.),
         prot = rep(top_prots, each = 2)) %>%
  inner_join(joinl %>% dplyr::select(prot, Gene) %>% distinct())

colnames(regResdf)[4] <- 'pvalue'

write.csv(regResdf, 'RegResTopC1VSTotal.csv')

regResdf <- regResdf %>% filter(str_detect(coef, 'intensity')) %>%
  mutate(adjp = p.adjust(pvalue, method = 'BH' ))

#
#
#

# visualise predicted values

gamplot <- function(gene){
  subset <- joinl %>% filter(Gene == gene)
  
  mod <- glm(Total.events ~ intensity,
             data = subset,
             family = Gamma(link = 'log'))
  
  pred <- subset %>% mutate(pred = exp(predict(mod, newdata = subset)))
  
  plot <- pred %>% ggplot(aes(x = intensity, y = Total.events)) + geom_point() +
    geom_line(aes(x = intensity, y = pred))
  
  print(plot)
}

gamplot('MUC16')

#
#
#

# Enrichment analysis of selected proteins

library(clusterProfiler)
library("org.Hs.eg.db")

# add leading protein names
loadingsPLS$Accession <- gsub(';', '.', loadingsPLS$Accession)
loadingsPLS <- loadingsPLS %>% left_join(prot[,1:2])
colnames(loadingsPLS)[7] <- 'LeadingProt'

## convert gene symbols into a character vector
prots <- as.vector(loadingsPLS$LeadingProt)


# convert accessions to entrez ID
hs <- org.Hs.eg.db

entrezIDp <- mapIds(hs, 
                   keys = posgenes$Gene,
                   column = "ENTREZID",
                   keytype = "SYMBOL")


# map ENTREZID back to protein accessions and loadings
loadingsPLS$ENTREZID <- entrezID

# fill missing ENTREZID by hand
loadingsPLS <- loadingsPLS %>% mutate(ENTREZID = ifelse(
  LeadingProt =='A0A0C4', '28388', ifelse(
    LeadingProt == 'P06312', '28908', ifelse(
      LeadingProt == 'A0A075', '28804', ifelse(
        LeadingProt == 'A0A075', '28803', ifelse(
          LeadingProt == 'A0A0C4', '28942', ENTREZID
        )
      )
    )
  )))

# component 2
loadingsPLS <- loadingsPLS %>% mutate(ENTREZID = ifelse(
  LeadingProt =='P01782', '28451', ifelse(
    LeadingProt == 'P01714', '28797', ifelse(
      LeadingProt == 'P01876', ' 3493', ENTREZID))))


# divide into loadings > 0 and < 0
neg <- loadingsPLS %>% filter(comp1 <0)
pos <- loadingsPLS %>% filter(comp1 >0)

## remove proteins with ENTREZID NAs
neg <- neg %>% drop_na(ENTREZID)
pos <- pos %>% drop_na(ENTREZID)

neg1 <- neg %>% left_join(prot %>% dplyr::select(Accession, Description))

## create background database - all proteins across 4 time points based on 2up converted to gene symbols with UniProt
bg <- prot[,2]

bgID <- mapIds(hs, 
               keys = bgene$Gene,
               column = "ENTREZID",
               keytype = "SYMBOL")

bgID1 <- data.frame(Reduce(rbind, bgID))
colnames(bgID1) <- 'ENTREZID'
bgID2 <- bgID1 %>% drop_na(ENTREZID)
bgIDv <- as.vector(bgID2$ENTREZID)

## convert entrez IDs into vector
entrezIDp <- as.vector(pos$ENTREZID)
entrezIDn <- as.vector(neg$ENTREZID)

## GO term annotation analysis

ggoP <- groupGO(gene = entrezIDp,
                OrgDb = org.Hs.eg.db, ## database 
                ont = 'BP', ## subontology (cellular component, mol function, bio process)
                level = 3, ## specificity of the term
                readable = TRUE) 

ggoN <- groupGO(gene = entrezIDn,
                OrgDb = org.Hs.eg.db, ## database - equine not available
                ont = 'BP', ## subontology (cellular component, mol function, bio process)
                level = 3, ## specificity of the term
                readable = TRUE) 

egoN <- enrichGO(gene = entrezIDn,
                 OrgDb = org.Hs.eg.db,
                 ont = 'ALL',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'fdr',
                 readable = TRUE,
                 universe = bgID)

write.csv(egoN@result, 'C1NegGoEnrich.csv')

egoP <- enrichGO(gene = entrezIDp,
                 OrgDb = org.Hs.eg.db,
                 ont = 'BP',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'BH',
                 readable = TRUE,
                 universe = bgID)
                 

# pathway analysis
## reactome pathway analysis
BiocManager::install("ReactomePA")
library(ReactomePA)

pathwaysN <- enrichPathway(gene = entrezIDn,
                           organism = 'human',
                           pvalueCutoff = 0.05,
                           pAdjustMethod = 'fdr',
                           qvalueCutoff = 0.9,
                           readable = TRUE)

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

#
#
#

# gene names

# background
bgene <- as.data.frame(prot$Description)
colnames(bgene) <- 'Description'
bgene <- bgene %>% separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(Gene)
write.csv(bgene, 'BGgenes.csv')

# c1 neg 
negenes <- neg %>% separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(Gene)
write.csv(negenes, 'C1NegGenes.csv')

posgenes <- pos %>% separate(Description, c('Description', 'Gene'), sep ='GN=') %>%
  dplyr::select(!Description) %>%
  separate(Gene, c('Gene', 'PE', 'SV'), sep = ' ') %>%
  dplyr::select(Gene)
write.csv(posgenes, 'C1PosGenes.csv')
