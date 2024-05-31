# Reflux project data integration using MOFA2

BiocManager::install("MOFA2")
library(MOFA2)
library(dplyr)
library(tidyr)

prot <- read.csv('ProteinVSNTransformed.csv')
ref <- read.csv('Select reflux 01-08-2023.csv')
clin <- read.csv('GOR_ClinicalEndpoint.csv')
clin2 <- read.csv('GORCF data 24-05.csv') %>% drop_na()

ref <- ref[-1,]
colnames(prot)[1] <- 'sample'

# format datasets to create MOFA object from long dataframe
ref <- ref %>% mutate(view = 'reflux')
refl <- ref %>% pivot_longer(cols = !c(PatID, view), names_to = 'feature', values_to = 'value')
colnames(refl)[1] <- 'sample'

prot <- prot %>% mutate(view = 'protein')
protl <- prot %>% pivot_longer(cols = !c(sample, view), names_to = 'feature', values_to = 'value')

data <- rbind(refl, protl)

mofa <- create_mofa(data)
plot_data_overview(mofa)
