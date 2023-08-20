

# Load necessary packages
library(tidyverse)
library(magrittr)
library(lattice)
library(afex)
library(emmeans)
library(lme4)
library(car)
library(data.table)
library(buildmer)
library(tidyr)
library(janitor)
library(broom.mixed)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

# Read data
source('./Analysis/Prepare_data_4.R')

#==========================================================================================#
#================================= 1. Explore individual data ==============================
#==========================================================================================#


# Plot RT by FP by participant and model using pooled data
FPfitAll=lm(meanRT ~ foreperiod,
            data=summaryData)

fit.params=tidy(FPfitAll)

summary(FPfitAll)


ggplot(data=summaryData,
       aes(x=foreperiod,
           y=meanRT)) +
  stat_summary(fun="mean", geom="point", size=1.5)+
  geom_abline(intercept=fit.params$estimate[1],
              slope=fit.params$estimate[2],
              color="blue")+
  facet_wrap(~ ID, ncol=6)


# Plot RT by FP by participant and model using individual data
dataGroupedByRT <- summaryData %>% 
  group_by(ID,foreperiod) %>% 
  summarise(meanRT=mean(meanRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)))

data.no_pooling <- dataGroupedByRT %>%
  select(-foreperiod) %>%
  group_by(ID) %>%
  nest(data = c(numForeperiod, meanRT)) %>%
  mutate(fit = map(data, ~ lm(meanRT ~ numForeperiod, data = .)),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, term, estimate) %>%
  complete(ID, term, fill = list(estimate = 0)) %>%
  pivot_wider(names_from = term,
              values_from = estimate) %>% 
  clean_names()


data.no_pooling <- data.no_pooling %>%
  rename(ID=id,
         numForeperiod=num_foreperiod)


ggplot(data = dataGroupedByRT,
       aes(x = numForeperiod, y = meanRT)) + 
  geom_abline(data = data.no_pooling,
              aes(intercept = intercept,
                  slope = numForeperiod),
              color = "blue") +
  geom_point() +
  facet_wrap(~ID, ncol=6) + 
  scale_x_continuous(breaks = 0:4 * 2) +
  theme(strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 12))


#==========================================================================================#
#====================================== 2. Prepare model ===================================
#==========================================================================================#

# Set contrasts for variables used in the models
contrasts(goData$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData$condition) <- c(-1/2, 1/2)
contrasts(goData$oneBacktrialType) <- c(-1/2, 1/2)
contrasts(goData$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(goData$oneBackFPDiff) <- contrastCodes

contrasts(goData2$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData2$condition) <- c(-1/2, 1/2)
contrasts(goData2$oneBacktrialType) <- c(-1/2, 1/2)
contrasts(goData2$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData2$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(goData2$oneBackFPDiff) <- contrastCodes

contrasts(dataAll$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(dataAll$condition) <- c(-1/2, 1/2)
contrasts(dataAll$oneBacktrialType) <- c(-1/2, 1/2)
contrasts(dataAll$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(dataAll$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(dataAll$oneBackFPDiff) <- contrastCodes

#=========================== 3.1. Foreperiod, condition and sequential effects =============================
# We use the strategy of keeping it maximal to find a model that converges and progressively
# remove terms, one of the strategies recommended to avoid overfitting:
# https://rdrr.io/cran/lme4/man/isSingular.html

#  Trying to fit the model using foreperiod and FPn-1 as factors results
# in R hanging during execution; for this reason, we use them as
# numerical variables

# We start by fitting a model using mixed, from afex.

fplmm <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod*condition*numOneBackFP|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                return = 'merMod',
                REML=TRUE)

summary(fplmm)

# Although this solves the fitting problem, it does not solve the singular fit issue, 
# which indicates that the model is too complex (i.e., it is overfitted:
# https://stats.stackexchange.com/questions/378939/dealing-with-singular-fit-in-mixed-models)

goData$scaledNumForeperiod <- scale(goData$numForeperiod)[,1]
goData$scaledNumOneBackFP <- scale(goData$numOneBackFP)[,1]

scalefplmm <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                  (1+scaledNumForeperiod*condition*scaledNumOneBackFP|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)

# This did not help, so we begin removing components

# 3.1.1.2. Remove correlations of mixed part
fplmm1v2 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod*condition*numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)
summary(fplmm1v2)
anova(fplmm1v2)

# Singular fit persists

# 3.1.1.3. Remove interactions of mixed part
fplmm1v3 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)
summary(fplmm1v3)

# We got rid of the singular fit. Additionally, this model incorporates reasonable assumptions
# about the random effects structure, since we have no a priori reasons to assume
# interactions or correlations between random effects in these data

# Choose dependent variable


# Fit models with RT and inverse RT without trimming
fplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

# If we ran this with the maximal structure, it does not converge and returns a
# singular fit. This is the case for all models below we ran with the maximal structure
# (not shown here).

# Let's check that the current structure does not provide a singular fit:
isSingular(fplmm1)

# Now we run the same model with inverse RT as outcome
invfplmm1 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+condition+numOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(invfplmm1)


anova(fplmm1)
anova(invfplmm1)

# Amount of variance accounted for by the model
cor(fitted(fplmm1), goData$RT)^2
cor(fitted(invfplmm1), goData$invRT)^2

# The first model explains a larger amount of the variance 

# Check normality of residuals
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")

# Both models show considerable departures from normality

# Plot residuals
plot(fplmm1, resid(.) ~ fitted(.),
     main="Residuals fplmm1")
plot(invfplmm1, resid(.) ~ fitted(.),
     main="Residuals invfplmm1")

# It appears that residuals correlate somewhat with fitted values; there are also outliers

# Residual histograms
qplot(resid(fplmm1),
      main="Residuals fplmm1")
qplot(resid(invfplmm1),
      main="Residuals invfplmm1")

# Both appear to be relatively normally distributed, although the first has 
# a larger positive skew

# Fit models with RT and inverse RT after outlier trimming
trimfplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

triminvfplmm1 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+condition+numOneBackFP||ID),
                   data = goData2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(trimfplmm1)
isSingular(triminvfplmm1)

anova(trimfplmm1)
anova(triminvfplmm1)
summary(triminvfplmm1)

# Amount of variance accounted for by the model
cor(fitted(trimfplmm1), goData2$RT)^2
cor(fitted(triminvfplmm1), goData2$invRT)^2

# Again, the first plot accounts for a larger amount of the variance

# check normality
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm1")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm1")

# The second model is much closer to normality

# Plot residuals
plot(trimfplmm1, resid(.) ~ fitted(.),
     main="Residuals trimfplmm1")
plot(triminvfplmm1, resid(.) ~ fitted(.),
     main="Residuals triminvfplmm1")

# Outliers are gone, but residuals still appear to correlate with fitted values in
# the first model

qplot(resid(trimfplmm1),
      main="Residuals trimfplmm1")
qplot(resid(triminvfplmm1),
      main="Residuals triminvfplmm1")

# Both still appear relatively normally distributed, with the second model
# performing better

# Use z-scores for centering with no outlier trimming
fplmm3 <- mixed(formula = RTzscore ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(fplmm3)

anova(fplmm3)

# Amount of variance accounted for by the model
cor(fitted(fplmm3), goData$RT)^2

# Very low variance explained

# check normality
qqnorm(resid(fplmm3),
       main="Normal q-qplot fplmm3")

# Far from normality

# Plot residuals
plot(fplmm3, resid(.) ~ fitted(.),
     main="Residuals fplmm3")

# Correlation between residuals and fitted values seems to dissappear

qplot(resid(fplmm3),
      main="Residuals fplmm3")

# Residuals are quite assymnetric

# Use z-scores for centering with trimming
fplmm4 <- mixed(formula = RTzscore ~ numForeperiod*condition*numOneBackFP + 
                  (1+numForeperiod+condition+numOneBackFP||ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(fplmm4)

# Singular fit

# Amount of variance accounted for by the model
cor(fitted(fplmm4), goData2$RT)^2

# check normality
qqnorm(resid(fplmm4),
       main="Normal q-qplot fplmm4")

# Plot residuals
plot(fplmm4, resid(.) ~ fitted(.),
     main="Residuals fplmm4")

qplot(resid(fplmm4),
      main="Residuals fplmm4")

# Only models with RT or inverse RT perform well

# Of those:
# Variance explained is higher with trimming
# Q-q plots are better with trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with inverse RT performs better than the one with RT


# Find random effects structure
triminvfplmm1 <- buildmer(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod*condition*numOneBackFP|ID),
                       data = goData2,
                       buildmerControl = buildmerControl(ddf = "Satterthwaite",
                                                         calc.anova = TRUE))
                      
isSingular(triminvfplmm1)


#==========================================================================================#
#==================================== 3. Model asssessment =================================
#==========================================================================================#

#=============================== 3.1. FP and FP n-1 as numerical ==================================
triminvfplmm1 <- mixed(formula = invRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                         condition + numForeperiod:condition + numOneBackFP:condition + 
                         numForeperiod:numOneBackFP:condition + (1 + condition | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

summary(triminvfplmm1)
anova(triminvfplmm1)

# Model comparisons

# Without FP n-1
triminvfplmm2 <- mixed(formula = invRT ~ 1 + numForeperiod + condition + numForeperiod:condition + 
                         (1 + condition | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

anova(triminvfplmm2, triminvfplmm1)

# Without condition
triminvfplmm3 <- mixed(formula = invRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                          (1 | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

anova(triminvfplmm3, triminvfplmm1)

# Without foreperiod
triminvfplmm4 <- mixed(formula = invRT ~ 1 + numOneBackFP + condition + numOneBackFP:condition + 
                         (1 + condition | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

anova(triminvfplmm4, triminvfplmm1)

# Without three way interaction
triminvfplmm5 <- mixed(formula = invRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                         condition + numForeperiod:condition + numOneBackFP:condition + 
                         (1 + condition | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

anova(triminvfplmm5, triminvfplmm1)

# Amount of variance accounted for by the models
cor(fitted(triminvfplmm1), goData2$RT)^2
cor(fitted(triminvfplmm2), goData2$invRT)^2
cor(fitted(triminvfplmm3), goData2$invRT)^2
cor(fitted(triminvfplmm4), goData2$invRT)^2
cor(fitted(triminvfplmm5), goData2$invRT)^2

# Arrange by BIC
BIC(triminvfplmm1, triminvfplmm2, triminvfplmm3, triminvfplmm4, triminvfplmm5) %>%
  arrange(BIC)

# The fullest model explains most of the variance and BICs are very close between
# this and the model without a random-effect from Fpn-1. In line with the 
# "keep it maximal" strategy, we keep this full model

#====================== 3.2. Using FP n-1 as categorical for emm comparisons ==========================

triminvfplmm6 <- buildmer(invRT ~ numForeperiod * condition * oneBackFP + 
                           (1+numForeperiod*condition*oneBackFP|ID), 
                         data=goData2,
                         buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(triminvfplmm6)
buildform <- formula(triminvfplmm6)


triminvfplmm6 <- mixed(formula = invRT ~  1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
                         condition + numForeperiod:condition + oneBackFP:condition + 
                         numForeperiod:oneBackFP:condition + (1 + condition + numForeperiod | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

# emm
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)
Fp_by_Previous=emtrends(triminvfplmm5, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(triminvfplmm5, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

#========================= 3.3. Both FP and FP n-1 as categorical ===================================

# 3.3.1. Using invRT
# Find optimal structure using buildmer
triminvfplmm7 <- buildmer(invRT ~ foreperiod * condition * oneBackFP + 
                           (1+foreperiod*condition*oneBackFP|ID), 
                         data=goData2,
                         buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmm7)


triminvfplmm7 <- mixed(formula = invRT ~  1 + foreperiod + oneBackFP + foreperiod:oneBackFP + 
                         condition + foreperiod:condition + oneBackFP:condition + 
                         foreperiod:oneBackFP:condition + (1 + condition + foreperiod + oneBackFP || ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

anova(triminvfplmm7)

# Pairwise comparisons by FP (estimates consecutive differences)
emmip(triminvfplmm7, condition ~ oneBackFP|foreperiod, CIs = TRUE)

triminvfplmm7emm <- emmeans(triminvfplmm7, ~ oneBackFP * condition|foreperiod)

triminvfplmm7emm <- emmeans(triminvfplmm7, pairwise ~ oneBackFP * condition|foreperiod)
contrast(triminvfplmm7emm[[1]],
         interaction = c("consec", "consec"),
         #by = "foreperiod",
         adjust = "holm")

triminvfplmm3emm <- emmeans(triminvfplmm3, pairwise ~ condition*oneBackFP|foreperiod)
contrast(triminvfplmm3emm[[1]], interaction = c("consec"), by = c("foreperiod", "oneBackFP"), adjust = "holm")
contrast(triminvfplmm3emm[[1]], interaction = c("consec"), by = c("foreperiod", "condition"), adjust = "holm")

# 3.3.2. Using RT
# Find optimal structure using buildmer
trimfplmm7 <- buildmer(RT ~ foreperiod * condition * oneBackFP + 
                            (1+foreperiod*condition*oneBackFP|ID), 
                          data=goData2,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm7)
formula(trimfplmm7)


trimfplmm7 <- mixed(formula = RT ~  1 + foreperiod + oneBackFP + foreperiod:oneBackFP + 
                         condition + foreperiod:condition + oneBackFP:condition + 
                         foreperiod:oneBackFP:condition + (1 + condition | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

anova(trimfplmm7)

# Pairwise comparisons by FP (estimates consecutive differences)
emmip(trimfplmm7, condition ~ oneBackFP|foreperiod, CIs = TRUE,
      xlab = "FP n-1",
      facelab = "label_both") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_color_manual(values = c("orange", "blue"))
ggsave("./Analysis/Plots/mixed_models_pairwise.png",
       width = 8.5,
       height = 5.7)

trimfplmm7emm <- emmeans(trimfplmm7, ~ oneBackFP * condition|foreperiod)

trimfplmm7emm <- emmeans(trimfplmm7, pairwise ~ oneBackFP * condition|foreperiod)
contrast(trimfplmm7emm[[1]], interaction = c("consec", "consec"), by = "foreperiod", adjust = "holm")

contrast(trimfplmm7emm[[1]], interaction = c("consec"), by = c("foreperiod", "oneBackFP"), adjust = "holm")
contrast(trimfplmm7emm[[1]], interaction = c("consec"), by = c("foreperiod", "condition"), adjust = "holm")

#===================================================================================#
# Sanity check: model comparisons without trimming
#===================================================================================#
invfplmm2 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+condition||ID),
                       data = goData,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(invfplmm2, invfplmm1, refit=FALSE)

invfplmm3 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+numForeperiod+numOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(invfplmm3, invfplmm1, refit=FALSE)

invfplmm4 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                     (1+condition+numOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(invfplmm4, invfplmm1, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(invfplmm1), goData$RT)^2
cor(fitted(invfplmm2), goData$invRT)^2
cor(fitted(invfplmm3), goData$invRT)^2
cor(fitted(invfplmm4), goData$invRT)^2

#=================== 3.3. Use RT and logRT to get priors for next experiment =====================

# logRT and FP and FP n-1 as numerical
trimlogfplmm1 <- buildmer(formula = logRT ~ numForeperiod*condition*numOneBackFP + 
                            (1+numForeperiod*condition*numOneBackFP|ID),
                          data = goData2,
                          buildmerControl = list(#direction='backward',
                            crit='LRT',#ddf = "Satterthwaite",
                            family=gaussian(link = 'identity'),
                            calc.anova = TRUE))

formula(trimlogfplmm1)

isSingular(trimlogfplmm1)

trimlogfplmm1 <- mixed(formula = logRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                         condition + numForeperiod:condition + condition:numOneBackFP + 
                         numForeperiod:condition:numOneBackFP + (1 + condition | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

summary(trimlogfplmm1)

# logRT and FP and FP n-1 as categorical
trimlogfplmm2 <- buildmer(formula = logRT ~ foreperiod*condition*oneBackFP + 
                            (1+foreperiod*condition*oneBackFP|ID),
                          data = goData2,
                          buildmerControl = list(#direction='backward',
                            crit='LRT',#ddf = "Satterthwaite",
                            family=gaussian(link = 'identity'),
                            calc.anova = TRUE))

formula(trimlogfplmm2)

# Fit model manually
trimlogfplmm2 <- mixed(formula = logRT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
                         foreperiod:condition + oneBackFP:condition + foreperiod:oneBackFP:condition + 
                         (1 + condition | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

isSingular(trimlogfplmm2)

summary(trimlogfplmm2)

# RT as linear and and FP and FP n-1 as numerical
trimfplmm1 <- buildmer(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod*condition*numOneBackFP|ID),
                       data = goData2,
                       buildmerControl = list(#direction='backward',
                         crit='LRT',#ddf = "Satterthwaite",
                         family=gaussian(link = 'identity'),
                         calc.anova = TRUE))

formula(trimfplmm1)

isSingular(trimfplmm1)

trimfplmm1 <- mixed(RT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                      condition + numForeperiod:condition + numOneBackFP:condition + 
                      numForeperiod:numOneBackFP:condition + (1 + condition | ID),
                    data=goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod",
                    check_contrasts = FALSE)

options(scipen = 999)
summary(trimfplmm1)


# RT as linear and FP and FP n-1 as categorical
trimfplmm2 <- buildmer(RT ~ foreperiod * oneBackFP * condition + 
                            (1+foreperiod * oneBackFP * condition|ID), 
                          data=goData2,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm2)
formula(trimfplmm2)


trimfplmm2 <- mixed(RT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
                      foreperiod:condition + oneBackFP:condition + foreperiod:oneBackFP:condition + 
                      (1 + condition | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

summary(trimfplmm2)

anova(trimfplmm2)

#============================================================================================================#
#===================================== 4. Sequential effects with n-1 trial type =============================
#============================================================================================================#

#======= Continuous variables ========

# RT
trimfplmmtrialtype1 <- buildmer(RT ~ numForeperiod * condition * oneBackFPGo + 
                                     (1+numForeperiod * condition * oneBackFPGo|ID), 
                                   data=goData2,
                                   buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimfplmmtrialtype1)

isSingular(trimfplmmtrialtype1)

trimfplmmtrialtype1 <- mixed(RT ~ 1 + numForeperiod + oneBackFPGo + numForeperiod:oneBackFPGo + 
                               condition + numForeperiod:condition + oneBackFPGo:condition + 
                               numForeperiod:oneBackFPGo:condition + (1 + condition + oneBackFPGo | 
                                                                        ID),
                                data=goData2,
                                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                progress = TRUE,
                                expand_re = TRUE,
                                method =  'KR',
                                REML=TRUE,
                                return = "merMod",
                                check_contrasts = FALSE)

summary(trimfplmmtrialtype1)

anova(trimfplmmtrialtype1)

# 1/RT
triminvfplmmtrialtype1 <- buildmer(invRT ~ numForeperiod * condition * oneBackFPGo + 
                            (1+numForeperiod * condition * oneBackFPGo|ID), 
                          data=goData2,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmmtrialtype1)

isSingular(triminvfplmmtrialtype1)

triminvfplmmtrialtype1 <- mixed(invRT ~ 1 + numForeperiod + oneBackFPGo + numForeperiod:oneBackFPGo + 
                                  condition + numForeperiod:condition + oneBackFPGo:condition + 
                                  numForeperiod:oneBackFPGo:condition + (1 + condition | ID),
                                data=goData2,
                                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                progress = TRUE,
                                expand_re = TRUE,
                                method =  'KR',
                                REML=TRUE,
                                return = "merMod",
                                check_contrasts = FALSE)

summary(triminvfplmmtrialtype1)

anova(triminvfplmmtrialtype1)

#=========== Categorical variables ============

# RT
trimfplmmtrialtype2 <- buildmer(RT ~ foreperiod * condition * oneBackFPGo + 
                                     (1+foreperiod * condition * oneBackFPGo|ID), 
                                   data=goData2,
                                   buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimfplmmtrialtype2)

isSingular(trimfplmmtrialtype2)

trimfplmmtrialtype2 <- mixed(RT ~ 1 + foreperiod + oneBackFPGo + foreperiod:oneBackFPGo + 
                               condition + foreperiod:condition + oneBackFPGo:condition + 
                               foreperiod:oneBackFPGo:condition + (1 + condition | ID),
                                data=goData2,
                                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                progress = TRUE,
                                expand_re = TRUE,
                                method =  'KR',
                                REML=TRUE,
                                return = "merMod",
                                check_contrasts = FALSE)


summary(trimfplmmtrialtype2)

anova(trimfplmmtrialtype2)

# 1/RT
triminvfplmmtrialtype2 <- buildmer(invRT ~ foreperiod * condition * oneBackFPGo + 
                                     (1+foreperiod * condition * oneBackFPGo|ID), 
                                   data=goData2,
                                   buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmmtrialtype2)

isSingular(triminvfplmmtrialtype2)

triminvfplmmtrialtype2 <- mixed(invRT ~ 1 + foreperiod + oneBackFPGo + foreperiod:oneBackFPGo + 
                                  condition + foreperiod:condition + oneBackFPGo:condition + 
                                  foreperiod:oneBackFPGo:condition + (1 + condition | ID),
                                data=goData2,
                                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                progress = TRUE,
                                expand_re = TRUE,
                                method =  'KR',
                                REML=TRUE,
                                return = "merMod",
                                check_contrasts = FALSE)

summary(triminvfplmmtrialtype2)

anova(triminvfplmmtrialtype2)


# FP n-1 and trial type as separate variables

#=========== Continuous variables ============

# RT
trimfplmmtrialtype3 <- buildmer(RT ~ numForeperiod * condition * numOneBackFP * oneBacktrialType + 
                                  (1 + numForeperiod * condition * numOneBackFPGo * oneBacktrialType |ID), 
                                data=goData2,
                                buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimfplmmtrialtype3)

isSingular(trimfplmmtrialtype3)

trimfplmmtrialtype3 <- mixed(RT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                               condition + oneBacktrialType + numForeperiod:condition + 
                               condition:oneBacktrialType + numForeperiod:oneBacktrialType + 
                               numOneBackFP:oneBacktrialType + numOneBackFP:condition + 
                               numForeperiod:numOneBackFP:condition + numOneBackFP:condition:oneBacktrialType + 
                               (1 + condition + oneBacktrialType | ID),
                             data=goData2,
                             control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                             progress = TRUE,
                             expand_re = TRUE,
                             method =  'KR',
                             REML=TRUE,
                             return = "merMod",
                             check_contrasts = FALSE)

summary(trimfplmmtrialtype3)

anova(trimfplmmtrialtype3)

# 1/RT
triminvfplmmtrialtype3 <- buildmer(invRT ~ numForeperiod * condition * numOneBackFP * oneBacktrialType + 
                                     (1 + numForeperiod * condition * numOneBackFPGo * oneBacktrialType |ID), 
                                   data=goData2,
                                   buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmmtrialtype3)

isSingular(triminvfplmmtrialtype3)

triminvfplmmtrialtype3 <- mixed(invRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
                                  oneBacktrialType + condition + numForeperiod:condition + 
                                  oneBacktrialType:condition + numForeperiod:oneBacktrialType + 
                                  numOneBackFP:oneBacktrialType + numForeperiod:numOneBackFP:oneBacktrialType + 
                                  numOneBackFP:condition + numForeperiod:numOneBackFP:condition + 
                                  numOneBackFP:oneBacktrialType:condition + (1 + condition + 
                                                                               oneBacktrialType | ID),
                                data=goData2,
                                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                progress = TRUE,
                                expand_re = TRUE,
                                method =  'KR',
                                REML=TRUE,
                                return = "merMod",
                                check_contrasts = FALSE)

summary(triminvfplmmtrialtype3)

anova(triminvfplmmtrialtype3)

#=========== Categorical variables ============

# RT
trimfplmmtrialtype4 <- buildmer(RT ~ foreperiod * condition * oneBackFP * oneBacktrialType + 
                                  (1+foreperiod * condition * oneBackFP * oneBacktrialType|ID), 
                                data=goData2,
                                buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimfplmmtrialtype4)

isSingular(trimfplmmtrialtype4)

trimfplmmtrialtype4 <- mixed(RT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
                               oneBacktrialType + foreperiod:condition + condition:oneBacktrialType + 
                               foreperiod:oneBacktrialType + oneBackFP:condition + foreperiod:oneBackFP:condition + 
                               foreperiod:oneBackFP:oneBacktrialType + condition:oneBackFP:oneBacktrialType + 
                               foreperiod:condition:oneBacktrialType +
                               foreperiod:oneBackFP:condition:oneBacktrialType +
                               (1 + condition + oneBacktrialType + condition:oneBacktrialType | 
                                  ID),
                             data=goData2,
                             control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                             progress = TRUE,
                             expand_re = TRUE,
                             method =  'KR',
                             REML=TRUE,
                             return = "merMod")#,
                             #check_contrasts = FALSE)

trimfplmmtrialtype4 <- mixed(RT ~ foreperiod * condition * oneBackFP * oneBacktrialType +
                               (1 + condition + oneBacktrialType + condition:oneBacktrialType | 
                                  ID),
                             data=goData2,
                             control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                             progress = TRUE,
                             expand_re = TRUE,
                             method =  'KR',
                             REML=TRUE,
                             return = "merMod")

summary(trimfplmmtrialtype4)

anova(trimfplmmtrialtype4)

#=========== By condition ===============
tt4action <- buildmer(RT ~ foreperiod * oneBackFP * oneBacktrialType + 
                  (1+foreperiod * oneBackFP * oneBacktrialType|ID), 
                  data=filter(goData2, condition == "action"),
                  buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))


isSingular(tt4action)

formula(tt4action)

tt4action <- mixed(RT ~ 1 + foreperiod + oneBacktrialType + oneBackFP + foreperiod:oneBackFP + 
                     foreperiod:oneBacktrialType + foreperiod:oneBackFP:oneBacktrialType + (1 + oneBacktrialType | ID),
                   data=filter(goData2, condition == "action"),
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod",
                   check_contrasts = FALSE)

summary(tt4action)

anova(tt4action)

# 1/RT
ttinv4action <- buildmer(invRT ~ foreperiod * oneBackFP * oneBacktrialType + 
                                     (1+foreperiod * oneBackFP * oneBacktrialType|ID), 
                                   data=filter(goData2, condition == "action"),
                                   buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(ttinv4action)

formula(ttinv4action)


ttinv4action <- mixed(invRT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + oneBacktrialType + 
                        foreperiod:oneBacktrialType + oneBackFP:oneBacktrialType +
                        foreperiod:oneBackFP:oneBacktrialType + 
                        (1 + oneBacktrialType | ID),
                      data=filter(goData2, condition == "action"),
                      control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method =  'KR',
                      REML=TRUE,
                      return = "merMod",
                      check_contrasts = FALSE)

summary(ttinv4action)

anova(ttinv4action)


# External
tt4external <- buildmer(RT ~ foreperiod * oneBackFP * oneBacktrialType + 
                        (1+foreperiod * oneBackFP * oneBacktrialType|ID), 
                      data=filter(goData2, condition == "external"),
                      buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))


isSingular(tt4external)

formula(tt4external)

tt4external <- mixed(RT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + oneBacktrialType + 
                     oneBackFP:oneBacktrialType + foreperiod:oneBacktrialType + 
                       foreperiod:oneBackFP:oneBacktrialType + 
                       (1 + oneBacktrialType | ID),
                   data=filter(goData2, condition == "external"),
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'KR',
                   REML=TRUE,
                   return = "merMod",
                   check_contrasts = FALSE)

summary(tt4external)

anova(tt4external)


# 1/RT
ttinv4external <- buildmer(invRT ~ foreperiod * oneBackFP * oneBacktrialType + 
                           (1+foreperiod * oneBackFP * oneBacktrialType|ID), 
                         data=filter(goData2, condition == "external"),
                         buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(ttinv4external)

formula(ttinv4external)


ttinv4external <- mixed(invRT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + oneBacktrialType + 
                          oneBackFP:oneBacktrialType + foreperiod:oneBacktrialType + 
                          foreperiod:oneBackFP:oneBacktrialType + 
                          (1 + foreperiod + oneBacktrialType | ID),
                        data=filter(goData2, condition == "external"),
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method =  'KR',
                        REML=TRUE,
                        return = "merMod",
                        check_contrasts = FALSE)

summary(ttinv4external)

anova(ttinv4external)

#============ by trial type ==================

# Go
tt4go <- buildmer(RT ~ foreperiod * oneBackFP * condition + 
                        (1+foreperiod * oneBackFP * condition |ID), 
                      data=filter(goData2, oneBacktrialType == "go"),
                      buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))


isSingular(tt4go)

formula(tt4go)

tt4go <- mixed(RT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
                 foreperiod:condition + oneBackFP:condition + foreperiod:oneBackFP:condition + 
                 (1 + condition | ID),
               data=filter(goData2, oneBacktrialType == "go"),
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'KR',
               REML=TRUE,
               return = "merMod",
               check_contrasts = FALSE)

summary(tt4go)

anova(tt4go)

# 1/RT
ttinv4go <- buildmer(invRT ~ foreperiod * oneBackFP * condition + 
                    (1+foreperiod * oneBackFP * condition |ID), 
                  data=filter(goData2, oneBacktrialType == "go"),
                  buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))


isSingular(ttinv4go)

formula(ttinv4go)

ttinv4go <- mixed(invRT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
                 foreperiod:condition + oneBackFP:condition + foreperiod:oneBackFP:condition + 
                 (1 + condition | ID),
               data=filter(goData2, oneBacktrialType == "go"),
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'KR',
               REML=TRUE,
               return = "merMod",
               check_contrasts = FALSE)

summary(ttinv4go)

anova(ttinv4go)


# No-go
tt4nogo <- buildmer(RT ~ foreperiod * oneBackFP * condition + 
                    (1+foreperiod * oneBackFP * condition |ID), 
                  data=filter(goData2, oneBacktrialType == "no-go"),
                  buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))


isSingular(tt4nogo)

formula(tt4nogo)

tt4nogo <- mixed(RT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + foreperiod:oneBackFP +
                 condition:oneBackFP + foreperiod:oneBackFP:condition +
                 (1 + condition | ID),
               data=filter(goData2, oneBacktrialType == "no-go"),
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'KR',
               REML=TRUE,
               return = "merMod",
               check_contrasts = FALSE)

summary(tt4nogo)

anova(tt4nogo)

# 1/RT
ttinv4nogo <- buildmer(invRT ~ foreperiod * oneBackFP * condition + 
                      (1+foreperiod * oneBackFP * condition |ID), 
                    data=filter(goData2, oneBacktrialType == "no-go"),
                    buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))


isSingular(ttinv4nogo)

formula(ttinv4nogo)

ttinv4nogo <- mixed(RT ~ 1 + foreperiod + condition + oneBackFP + foreperiod:condition + foreperiod:oneBackFP +
                      condition:oneBackFP + foreperiod:oneBackFP:condition +
                      (1 + condition | ID),
                    data=filter(goData2, oneBacktrialType == "no-go"),
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod",
                    check_contrasts = FALSE)

summary(ttinv4nogo)

anova(ttinv4nogo)


#============================================================================================================#
#=========================== 5. Difference between the durations of FP n and FP n-1 as predictor ==============
#============================================================================================================#

# Fullest model using difference as numerical predictor
fpDifflmm1 <- mixed(formula = invRT ~ numForeperiod * condition * numOneBackFPDiff + 
                      (1+numForeperiod*condition*numOneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm2 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

anova(fpDifflmm2)

fpDifflmm3 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm4 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

fpDifflmm5 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod+condition+oneBackFPDiff||ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    return = "merMod")

cor(fitted(fpDifflmm2), goData2$invRT)^2
cor(fitted(fpDifflmm4), goData2$invRT)^2
cor(fitted(fpDifflmm5), goData2$invRT)^2

# Fullest model using foreperiod and one back fp as categorical variables
fpDifflmm1 <- mixed(formula = invRT ~ condition * oneBackFPDiff + 
                      (1+condition*oneBackFPDiff|ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    check_contrasts=FALSE,
                    return = "merMod")

# Remove correlations
fpDifflmm2 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    check_contrasts=FALSE,
                    return = "merMod")

# Remove interactions
fpDifflmm3 <- mixed(formula = invRT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod+condition+oneBackFPDiff||ID),
                    data = goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE,
                    check_contrasts=FALSE,
                    return = "merMod")

# Find optimal structure using buildmer
fpDifflmm <- buildmer(invRT ~ foreperiod * condition * oneBackFPDiff + 
                        (1+foreperiod*condition*oneBackFPDiff|ID), 
                      data=goData2,
                      buildmerControl = buildmerControl(include = ~ foreperiod * condition * oneBackFPDiff, calc.anova = TRUE, ddf = "Satterthwaite"))


#==============================================================================================#
#============================== 6. Model with scaled numerical predictors ======================
#==============================================================================================#

# Models with scaled variables

goData$scaledNumForeperiod <- scale(goData$numForeperiod)[,1]
goData$scaledNumOneBackFP <- scale(goData$numOneBackFP)[,1]
goData2$scaledNumForeperiod <- scale(goData2$numForeperiod)[,1]
goData2$scaledNumOneBackFP <- scale(goData2$numOneBackFP)[,1]

#============================ 6.1. Model with inverse Rt as DV ==================
# Find random factor structure with buildmer
scaledtriminvfplmm1 <- buildmer(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                                  (1+scaledNumForeperiod*condition*scaledNumOneBackFP|ID),
                                data = goData2,
                                buildmerControl = list(#direction='backward',
                                  crit='LRT',#ddf = "Satterthwaite",
                                  family=gaussian(link = 'identity'),
                                  calc.anova = TRUE))

formula(scaledtriminvfplmm1)

scaledtriminvfplmm1 <- mixed(formula = invRT ~ 1 + scaledNumForeperiod + scaledNumOneBackFP + scaledNumForeperiod:scaledNumOneBackFP + 
                               condition + scaledNumForeperiod:condition + scaledNumOneBackFP:condition + 
                               scaledNumForeperiod:scaledNumOneBackFP:condition + (1 + condition | 
                                                                                     ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)

summary(scaledtriminvfplmm1)

#============================ 6.2. Model with logRt as DV ==================

# Find random factor structure with buildmer
scaledtrimlogfplmm1 <- buildmer(formula = logRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                                  (1+scaledNumForeperiod*condition*scaledNumOneBackFP|ID),
                      data = goData2,
                      buildmerControl = list(#direction='backward',
                                           crit='LRT',#ddf = "Satterthwaite",
                                           family=gaussian(link = 'identity'),
                                           calc.anova = TRUE))

formula(scaledtrimlogfplmm1)

# Fit model manually
scaledtrimlogfplmm1 <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + scaledNumOneBackFP + scaledNumForeperiod:scaledNumOneBackFP + 
                               condition + scaledNumForeperiod:condition + scaledNumOneBackFP:condition + 
                               scaledNumForeperiod:scaledNumOneBackFP:condition + (1 + condition + 
                                                                                     scaledNumOneBackFP + scaledNumForeperiod | ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod",
                       check_contrasts = FALSE)


summary(scaledtrimlogfplmm1)


#===============================================================================================#
#===================================== 7. Accuracy ==============================================
#===============================================================================================#

fpacc1stlevel <- glmer(acc_result ~ 1 + scaledNumForeperiod * condition * scaledNumOneBackFP +
                         (1 + condition | ID),
                       data = dataAll,
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer = c("bobyqa"),
                                              optCtrl = list(maxfun = 2e5)))


isSingular(fpacc1stlevel)

summary(fpacc1stlevel)


fpacc2ndlevel <- glmer(acc_result ~ 1 + (scaledNumForeperiod + squaredScaledNumForeperiod) * condition * scaledNumOneBackFP +
                         (1 + condition | ID),
                       data = dataAll,
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer = c("bobyqa"),
                                              optCtrl = list(maxfun = 2e5)))

isSingular(fpacc2ndlevel)

summary(fpacc2ndlevel)

anova(fpacc2ndlevel)

anova(fpacc2ndlevel, fpacc1stlevel)
