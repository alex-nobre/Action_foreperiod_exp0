

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
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

source('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Prepare_data_4.R')

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

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(goData$oneBackFPDiff) <- contrastCodes

contrasts(goData2$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(goData2$condition) <- c(-1/2, 1/2)
contrasts(goData2$oneBacktrialType) <- c(-1/2, 1/2)
contrasts(goData2$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)

contrastCodes <- cbind(c(3/4, -1/4, -1/4, -1/4),
                       c(-1/4, 3/4, -1/4, -1/4),
                       c(-1/4, -1/4, 3/4, -1/4)) %>%
  set_colnames(c('nogoVsShorter', 'longerVsShorter', 'equalVsShorter'))

contrasts(goData2$oneBackFPDiff) <- contrastCodes

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

#==========================================================================================#
#====================================== 3. Build model ===================================
#==========================================================================================#

# Model comparisons
triminvfplmm2 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+condition||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(triminvfplmm2, triminvfplmm1, refit=FALSE)

triminvfplmm3 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+numForeperiod+numOneBackFP||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(triminvfplmm3, triminvfplmm1, refit=FALSE)

triminvfplmm4 <- mixed(formula = invRT ~ numForeperiod*condition*numOneBackFP + 
                         (1+condition+numOneBackFP||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(triminvfplmm4, triminvfplmm1, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(triminvfplmm1), goData2$RT)^2
cor(fitted(triminvfplmm2), goData2$invRT)^2
cor(fitted(triminvfplmm3), goData2$invRT)^2
cor(fitted(triminvfplmm4), goData2$invRT)^2

# Arrange by BIC
BIC(triminvfplmm, triminvfplmm2, triminvfplmm3, triminvfplmm4) %>%
  arrange(BIC)

# The fullest model explains most of the variance and BICs are very close between
# this and the model without a random-effect from Fpn-1. In line with the 
# "keep it maximal" strategy, we keep this full model

# Find optimal structure using buildmer
triminvfplmm <- buildmer(invRT ~ foreperiod * condition * oneBackFP + 
                        (1+foreperiod*condition*oneBackFP|ID), 
                      data=goData2,
                      buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmm)

summary(triminvfplmm)

triminvfplmm <- buildmer(invRT ~ numForeperiod * condition * numOneBackFP + 
                           (1+numForeperiod*condition*numOneBackFP|ID), 
                         data=goData2,
                         buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmm)

summary(triminvfplmm)

# Using FP n-1 as categorical for emm comparisons
triminvfplmm5 <- buildmer(invRT ~ numForeperiod * condition * oneBackFP + 
                           (1+numForeperiod*condition*oneBackFP|ID), 
                         data=goData2,
                         buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(triminvfplmm5)
formula(triminvfplmm5)


triminvfplmm5 <- mixed(invRT ~  1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
                         condition + numForeperiod:condition + oneBackFP:condition + 
                         numForeperiod:oneBackFP:condition + (1 + condition + numForeperiod | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

# emm
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)
Fp_by_Previous=emtrends(triminvfplmm5, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(triminvfplmm5, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")


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

#=================== Use logRT to get priors for next experiment =====================
trimlogfplmm1 <- buildmer(formula = logRT ~ foreperiod*condition*oneBackFP + 
                                  (1+foreperiod*condition*oneBackFP|ID),
                                data = goData2,
                                buildmerControl = list(#direction='backward',
                                  crit='LRT',#ddf = "Satterthwaite",
                                  family=gaussian(link = 'identity'),
                                  calc.anova = TRUE))

formula(trimlogfplmm1)

# Fit model manually
trimlogfplmm1 <- mixed(formula = logRT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
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

isSingular(trimlogfplmm1)

summary(trimlogfplmm1)

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
#=========================== 4. Difference between the durations of FP n and FP n-1 as predictor ==============
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
#============================== 5. Model with scaled numerical predictors ======================
#==============================================================================================#

# Models with scaled variables

goData$scaledNumForeperiod <- scale(goData$numForeperiod)[,1]
goData$scaledNumOneBackFP <- scale(goData$numOneBackFP)[,1]
goData2$scaledNumForeperiod <- scale(goData2$numForeperiod)[,1]
goData2$scaledNumOneBackFP <- scale(goData2$numOneBackFP)[,1]

#============================ 5.1. Model with inverse Rt as DV ==================
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

#============================ 5.2. Model with logRt as DV ==================

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

library(sjPlot)
library(prediction)



