# Reflux study analysis

# title:
# author:
# output:

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(lme4)
library(lmerTest)
library(stringr)
library(gridExtra)
library(effects)

ref <- read.csv('Select reflux 01-08-2023.csv')
clinu <- read.csv('GOR_ClinicalEndpoint.csv')
clinu2 <- read.csv('GOR-CF database 07-03-2024.csv')

ref <- ref[-1,]

# create consistent patient id
clinu2$Assignment.Number <- gsub('GOR', 'CF', clinu2$Assignment.Number)
colnames(clinu2)[1] <- 'PatID'

# create joint clinical-reflux database with relevant variables
data <- clinu %>% dplyr::select(PatID, BaselineFEV1pp, IV.Days.12months, X1y_priorFEV1., X2y_prior_FEV1.) %>%
  left_join(ref %>% 
              dplyr::select(PatID, Total.events, Proximal.total.events, Demeester.score, Total.acid.exposure)) %>%
  left_join(clinu2 %>% dplyr::select(PatID, Age, BMI, Pseudomonas.colonisation, 
                                    Pancreatic.insufficency, Diabetes, CFTR.Modulator))

data$BMI <- as.numeric(data$BMI)


# shorten column names
colnames(data) <- str_split(colnames(data), "\\.", simplify=T)[,1]
colnames(data) <- as.character(colnames(data))

colnames(data)[9] <- 'Total.acid'
data <- data %>% mutate(Pseudomonas = ifelse(Pseudomonas == 'TBC', NA, Pseudomonas))


# plot clinical and reflux continuous variables against each other
tiff('RegCovPlots.tiff', res = 300, width = 130, height = 110, unit = 'mm')
pairs(data[,c(2:3, 5:7)])
dev.off()

cor(data[,2:7], use = 'complete.obs')

# plot relationship between outcomes and categorical variables
fevplot <- data %>% pivot_longer(cols = c(Pseudomonas, Pancreatic, Diabetes), names_to = 'variable', values_to = 'value') %>%
  drop_na(value) %>%
  ggplot(aes(x = variable, y = BaselineFEV1pp, colour = value)) + geom_boxplot() + geom_point(position = position_jitterdodge()) +
  theme_bw(base_size = 8) + scale_color_brewer(palette = 'Set2') + ylab('FEV') + theme(axis.title.x = element_blank())
  
fevplot

ivplot <- data %>% pivot_longer(cols = c(Pseudomonas, Pancreatic, Diabetes), names_to = 'variable', values_to = 'value') %>%
  drop_na(value) %>%
  ggplot(aes(x = variable, y = IV, colour = value)) + geom_boxplot() + geom_point(position = position_jitterdodge()) +
  theme_bw(base_size = 8) + scale_color_brewer(palette = 'Set2') + ylab('IV antibiotic days')

g1 <- arrangeGrob(fevplot, ivplot, ncol = 1 , nrow = 2)
plot(g1)

ggsave('RegVocBox.tiff', g1, width = 120, height = 70, unit = 'mm', dpi = 300)

# plot FEV change in time for each patient
datal <- data %>% pivot_longer(cols = c(BaselineFEV1pp, X1y_priorFEV1, X2y_prior_FEV1), 
                               names_to = 'TimePt', values_to = 'FEV') %>%
  mutate(TimePt = ifelse(TimePt == 'X1y_priorFEV1', -1, ifelse(
    TimePt == 'X2y_prior_FEV1', -2, 0
  )))

fevplot2 <- datal %>% ggplot(aes(x = TimePt, y = FEV, colour = PatID)) + 
  geom_point() + 
  geom_line(aes(group = PatID)) +
  theme(legend.position = 'none') +
  scale_x_discrete(limits=c(-2, -1, 0))

fevplot2a <- datal %>% ggplot(aes(x = TimePt, y = FEV)) + 
  geom_jitter() +
  scale_x_discrete(limits = c(-2, -1, 0)) +
  geom_smooth(method = 'lm')

fevplot2a

datal %>% ggplot(aes(x = as.factor(TimePt), y = FEV, fill = Pseudomonas)) + 
  geom_boxplot()

datal %>% drop_na(Pseudomonas) %>% ggplot(aes(x = TimePt, y = FEV, colour = Pseudomonas)) +
  geom_smooth(method = 'lm') + geom_jitter() +
  scale_x_discrete(limits=as.factor(c(-2, -1, 0))) + xlim(-2, 0)

#

# histograms of outcome variables
tiff(filename = 'HistFEV.tiff', res = 300, width = 90, height = 70, unit = 'mm')
hist(data$BaselineFEV1pp, breaks = 15, main = 'Forced Expiratory Volume 1',
     xlab = '% predicted FEV1', cex.main = 0.9) 
dev.off()

tiff(filename = 'HistIVDays.tiff', res = 300, width = 90, height = 70, unit = 'mm')
hist(data$IV, breaks = 15, main = 'IV antibiotic treatment in 12 months',
     xlab = 'Days', cex.main = 0.9)
dev.off()

# predict FEV
library(caret)
library(leaps)
library(olsrr)

fev1 <- lm(BaselineFEV1pp ~ Proximal + Age + BMI + Pseudomonas +  Diabetes,
          data = data)

fev1 <- lm(BaselineFEV1pp ~ Age + BMI + Proximal + Pseudomonas,
           data = data)

fev1 <- lm(BaselineFEV1pp ~ Proximal  + Age + BMI + Pseudomonas + Diabetes,
           data = data)

summary(fev1)

# all possible predictor combinations
allfev <- as.data.frame(ols_step_all_possible(fev1)[1])

# best subset selection
bestfev <- ols_step_best_subset(fev1)[['metrics']]

cp <- bestfev %>% ggplot(aes(x = n, y = cp)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('Cp')

bic <- bestfev %>% ggplot(aes(x = n, y = sbc)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('BIC')

r2 <- bestfev %>% ggplot(aes(x = n, y = adjr)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('adj R2')

bsfev <- arrangeGrob(cp, bic, r2, nrow = 1, ncol = 3)
plot(bsfev)

ggsave('FEVbestModel.tiff', bsfev, width = 160, height = 60, unit = 'mm', dpi = 300)

modelsfev <- regsubsets(BaselineFEV1pp ~ ., data = data %>% dplyr::select(!c(PatID, IV, Total, CFTR)), nvmax = 6) 

summary(modelsfev)

sumfev <- summary(modelsfev)

# choosing the best model based on test error estimates
cbind( 
  Cp     = sumfev$cp,
  r2     = sumfev$rsq,
  Adj_r2 = sumfev$adjr2,
  BIC    = sumfev$bic
)

cpfev <- data.frame(predictors = c(1:6), cp = sumfev$cp)
bicfev <- data.frame(predictors = c(1:6), bic = sumfev$bic)

cpfev %>% ggplot(aes(x = predictors, y = cp)) + geom_point()
bicfev %>% ggplot(aes(x = predictors, y = bic)) + geom_point()

# final model

fev2 <- lm(BaselineFEV1pp ~ Proximal + Age + BMI,
           data = data)
AIC(fev2)
BIC(fev2)

summary(fev2)

write.csv(summary(fev1)[['coefficients']], 'FEVRegCoefs.csv')

# model diagnostics
par(mfrow = c(2,2))
dev.new()
plot(fev1)


# filter only combinations including proximal variable
allfevprox <- allfev %>% filter(grepl('Proximal', result.predictors))

# choose best model for each variable number
allfevprox %>% group_by(result.n) %>% summarise(best = max(result.aic))

#
#
#

# plot marginal effect of proximal reflux events on FEV

# predict FEV values fixing all the other predictor variables
newdata <- data %>% mutate(Age = mean(na.omit(Age)),
                           BMI = mean(na.omit(BMI)),
                           Pseudomonas = 'No',
                           Diabetes = 'No')


pred <- data %>% mutate(FEV1pred = predict(fev1, newdata = newdata))

pred1 <- predict(fev1,
                 newdata = newdata,
                 interval = 'prediction',
                 level = 0.95) %>%
  cbind(newdata)



dev.new()

rev <- pred1 %>% ggplot(aes(x = Proximal, y = BaselineFEV1pp)) +
  geom_point() + geom_line(aes(x = Proximal, y = fit), colour = 'blue') +
  geom_line(aes(x = Proximal, y = upr), colour = 'blue', linetype = 'dashed') +
  geom_line(aes(x = Proximal, y = lwr), colour = 'blue', linetype = 'dashed') +
  #ggtitle('Marginal effect of reflux on FEV1 % predicted') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  xlab('Proximal reflux events (number)') + ylab('FEV1 % predicted') 

rev

ggsave('RefluxVSFEV.tiff', dpi = 300, width = 100, height = 80, unit = 'mm')

#
#
#

# predicting FEV change in time
# random intercept and slope model 

hist(datal$FEV)

fevC <- lmer(FEV ~ TimePt*(Proximal + Age + BMI + Pseudomonas + Diabetes) + (1 | PatID) + (TimePt-1|PatID), # random slope and intercept uncorrelated 
           data = datal)

fevC1 <- lmer(FEV ~ TimePt*(Proximal + Age + BMI + Pseudomonas + Diabetes) + (1 + TimePt|PatID), # random slope and intercept correlated 
             data = datal)

summary(fevC)
summary(fevC1)

fevC <- lmer(FEV ~ TimePt*Proximal + TimePt*Pseudomonas + Age + BMI + (1 | PatID) + (TimePt-1|PatID), # random slope and intercept uncorrelated 
             data = datal)

summary(fevC)

# predicting IV days
# two-part model

# create two variables describing outcome (presence-absence and magnitude when present)
data <- data %>% mutate(IV1 = ifelse(IV == 0, 0, 1),
                        IV2 = ifelse(IV == 0, NA, IV),
                        IV3 = ifelse(IV == 0, 0.01, IV))


# negative binomial regression
library(MASS)
library(topmodels)

mean(na.omit(data$IV))
var(na.omit(data$IV))
sd(na.omit(data$IV))

iv <- glm.nb(IV ~ Proximal + Pseudomonas + Diabetes + BaselineFEV1pp,
             data = data)

summary(iv)

AIC(iv)

write.csv(summary(iv)[['coefficients']], 'IVRegCoefs.csv')

best <- cbind(Estimate = coef(iv), confint(iv))

# incidence ratios
exp(est)

# residuals
iv$resid

library(DHARMa)

resnb <- simulateResiduals(iv, refit = T)
plot(resnb)

# dispersion
testDispersion(resnb, plot=F)

summary(iv)
rootogram(iv)

ivp <- glm(IV ~ Proximal + Pseudomonas + Diabetes + BaselineFEV1pp,
            data = data,
            family = 'poisson')

summary(ivp)

# plot marginal effect of reflux

newdata1 <- data %>% mutate(Pseudomonas = 'Yes',
                            Diabetes = 'Yes',
                            BaselineFEV1pp = mean(na.omit(BaselineFEV1pp)))


pred2 <- cbind(newdata1, predict(iv, newdata = newdata1, se.fit = TRUE))



est <- cbind(Estimate = coef(iv), confint(iv))
est <- exp(est)

rev <- pred1 %>% ggplot(aes(x = Proximal, y = BaselineFEV1pp)) +
  geom_point() + geom_line(aes(x = Proximal, y = fit), colour = 'blue') +
  geom_line(aes(x = Proximal, y = upr), colour = 'blue', linetype = 'dashed') +
  geom_line(aes(x = Proximal, y = lwr), colour = 'blue', linetype = 'dashed') +

#

riv <- pred1 %>% ggplot(aes(x = Proximal, y = IV)) +
  geom_point() + geom_line(aes(x = Proximal, y = IVpred), colour = 'blue') +
  ggtitle('Marginal effect of reflux on IV day number') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  xlab('Proximal reflux events (number)') + ylab('IV days') 

ggsave('RefluxVSiv.tiff', dpi = 300, width = 100, height = 80, unit = 'mm')

#
#

newdata2 <- data %>% mutate(Pseudomonas = 'Yes',
                            Diabetes = 'Yes',
                            Proximal = mean(na.omit(Proximal)))


pred2 <- data %>% mutate(IVpred = predict(iv, newdata = newdata2),
                         IVpred1 = exp(predict(iv, newdata = newdata2)))

#

pred2 %>% ggplot(aes(x = BaselineFEV1pp, y = IV)) +
  geom_point() + geom_line(aes(x = BaselineFEV1pp, y = IVpred), colour = 'blue') +
  ggtitle('Marginal effect of reflux on IV day number') +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  xlab('Proximal reflux events (number)') + ylab('IV days') 








# TWO-STEP MODEL
# part 1 - logistic regression
data1 <- data %>% dplyr::select(IV1, Proximal, Pseudomonas, Diabetes, BMI, Age) %>% drop_na()
data1$IV1 <- as.factor(data1$IV1)

iv1 <- glm(IV1 ~ Proximal + Age + BMI + Pseudomonas + Diabetes,
           data = data1,
           family = binomial(link = 'logit'))

iv1 <- glm(IV1 ~ Proximal + Pseudomonas + Diabetes,
           data1,
           family = binomial(link = 'logit'))

summary(iv1)

# visualise marginal effects
data1 <- data %>% mutate(Diabetes = 'Yes', Proximal = median(Proximal))
data1 <- data1 %>% mutate(pred = plogis(predict(iv1, newdata = data1)))

effects <- allEffects(iv1)
ci <- data.frame(lower = plogis(effects$Pseudomonas$lower)) %>%
  mutate(Pseudomonas = c('No', 'Yes')) %>%
  mutate(upper = plogis(effects$Pseudomonas$upper))

data1 %>% left_join(ci) %>% drop_na() %>%
  ggplot(aes(y = pred, x = Pseudomonas)) + geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .1) +
  theme_bw()

data1 %>% drop_na() %>% ggplot(aes(y = pred, x = Pseudomonas)) + geom_point() +
  ylim(0, 1)

# all predictor combinations
alliv1 <- as.data.frame(ols_step_all_possible(iv1)[1])

# filter only combinations includign proximal variable
alliv1prox <- alliv1 %>% filter(grepl('Proximal', result.predictors))

# choose best model for each variable number
alliv1prox %>% group_by(result.n) %>% summarise(best = max(result.aic))


# best subset selection
ols_step_best_subset(iv)
bestiv1 <- ols_step_best_subset(iv)[['metrics']]

aic <- bestiv1 %>% ggplot(aes(x = n, y = aic)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('AIC')

bic <- bestiv1 %>% ggplot(aes(x = n, y = sbc)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('BIC')

r2 <- bestiv1 %>% ggplot(aes(x = n, y = adjr)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('adj R2')

bsiv <- arrangeGrob(aic, bic, r2, nrow = 1, ncol = 3)
plot(bsiv)

ggsave('IV1bestModel.tiff', bsiv, width = 160, height = 60, unit = 'mm', dpi = 300)

### final model - step 1
iv1a <- glm(IV1 ~ Pseudomonas + BMI + Pancreatic + Proximal,
           data = data,
           family = binomial(link = 'logit'))

summary(iv1a)

write.csv(summary(iv1a)[['coefficients']], 'IV1RegCoefs.csv')




# part 2 - gamma regression
iv2 <- glm(IV2 ~ Proximal + BMI + Age + Pseudomonas + Diabetes,
           data = data %>% drop_na(IV2),
           family = Gamma(link = 'identity'))

iv2 <- glm(IV2 ~ Proximal + Pseudomonas + Diabetes,
           data = data %>% drop_na(IV2),
           family = Gamma(link = 'identity'))

summary(iv2)

plot(allEffects(iv2))

# all predictor combinations
alliv2 <- as.data.frame(ols_step_all_possible(iv2)[1])

# filter only combinations includign proximal variable
alliv2prox <- alliv2 %>% filter(grepl('Proximal', result.predictors))

# choose best model for each variable number
alliv2prox %>% group_by(result.n) %>% summarise(best = max(result.aic))

# best subset selection
ols_step_best_subset(iv2)
bestiv2 <- ols_step_best_subset(iv2)[['metrics']]

cp <- bestiv2 %>% ggplot(aes(x = n, y = aic)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('AIC')

bic <- bestiv2 %>% ggplot(aes(x = n, y = sbc)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('BIC')

r2 <- bestiv2 %>% ggplot(aes(x = n, y = adjr)) + geom_point(colour = 'royalblue4') + geom_line(colour = 'hotpink4') + theme_bw() +
  xlab('Number of Predictors') + ylab('adj R2')

bsiv2 <- arrangeGrob(cp, bic, r2, nrow = 1, ncol = 3)
plot(bsiv2)

ggsave('IV2bestModel.tiff', bsiv2, width = 160, height = 60, unit = 'mm', dpi = 300)

# final model - step 1
iv2a <- glm(IV2 ~  BMI + Pseudomonas + Diabetes + Age + Proximal,
            data = data %>% drop_na(IV2),
            family = Gamma(link = 'inverse'))

summary(iv2a)
write.csv(summary(iv2a)[['coefficients']], 'IV2RegCoefs.csv')



# plot outcome against reflux values

plot(data$Proximal, data$IV)
ext <- c('CF020', 'CF035', 'CF007', 'CF009', 'CF008', 'CF003', 'CF001')

data <- data %>% mutate(ext = ifelse(PatID %in% ext, 'Ext', 'NExt'))
data %>% ggplot(aes(x = ext, y = IV)) + geom_boxplot() + geom_jitter() + geom_text(aes(label = PatID))


