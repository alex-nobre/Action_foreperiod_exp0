


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

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)

# Read data
source('./Analysis/Prepare_data_4.R')

#=========================== 1. Choose dependent variable =============================
# We use the strategy of keeping it maximal to find a model that converges and progressively
# remove terms, one of the strategies recommended to avoid overfitting:
# https://rdrr.io/cran/lme4/man/isSingular.html

#  Trying to fit the model using foreperiod and FPn-1 as factors results
# in R hanging during execution; for this reason, we use them as
# numerical variables


scaledfplmm <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                      (1+scaledNumForeperiod*condition*scaledNumOneBackFP|ID),
                    data = goData,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE)

# This results in a singular fit, so we begin removing components

# 3.1.1.2. Remove correlations of mixed part
scalefplmm1v2 <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                    (1+scaledNumForeperiod*condition*scaledNumOneBackFP||ID),
                  data = goData,
                  control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                  progress = TRUE,
                  expand_re = TRUE,
                  method =  'S',
                  return = 'merMod',
                  REML=TRUE)
summary(fplmm1v2)
anova(fplmm1v2)

# Singular fit persists

# 3.1.1.3. Remove interactions of mixed part
scaledfplmm1v3 <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                    (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
                  data = goData,
                  control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                  progress = TRUE,
                  expand_re = TRUE,
                  method =  'S',
                  REML=TRUE)
summary(scaledfplmm1v3)
anova(scaledfplmm1v3)

# We got rid of the singular fit. Additionally, this model incorporates reasonable assumptions
# about the random effects structure, since we have no a priori reasons to assume
# interactions or correlations between random effects in these data

# Choose dependent variable


# Fit models with RT and inverse RT without trimming
scaledfplmm1 <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                  (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
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
isSingular(scaledfplmm1)

# Now we run the same model with inverse RT as outcome
scaledinvfplmm1 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                           (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

isSingular(scaledinvfplmm1)


anova(scaledfplmm1)
anova(scaledinvfplmm1)

# Amount of variance accounted for by the model
cor(fitted(scaledfplmm1), goData$RT)^2
cor(fitted(scaledinvfplmm1), goData$invRT)^2

# The first model explains a larger amount of the variance 

# Check normality of residuals
qqnorm(resid(scaledfplmm1),
       main="Normal q-qplot scaledfplmm1")
qqnorm(resid(scaledinvfplmm1),
       main="Normal q-qplot scaledinvfplmm1")

# Both models show considerable departures from normality

# Plot residuals
plot(scaledfplmm1, resid(.) ~ fitted(.),
     main="Residuals scaledfplmm1")
plot(scaledinvfplmm1, resid(.) ~ fitted(.),
     main="Residuals scaledinvfplmm1")

# It appears that residuals correlate somewhat with fitted values; there are also outliers

# Residual histograms
qplot(resid(scaledfplmm1),
      main="Residuals scaledfplmm1")
qplot(resid(scaledinvfplmm1),
      main="Residuals scaledinvfplmm1")

# Both appear to be relatively normally distributed, although the first has 
# a larger positive skew

# Fit models with RT and inverse RT after outlier trimming
scaledtrimfplmm <- mixed(formula = RT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                           (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
                   data = goData2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

scaledtriminvfplmm <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                        (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
                      data = goData2,
                      control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method =  'S',
                      REML=TRUE,
                      return = "merMod")

isSingular(scaledtrimfplmm)
isSingular(scaledtriminvfplmm)

anova(scaledtrimfplmm)
anova(scaledtriminvfplmm)

# Amount of variance accounted for by the model
cor(fitted(scaledtrimfplmm), goData2$RT)^2
cor(fitted(scaledtriminvfplmm), goData2$invRT)^2

# Again, the first plot accounts for a larget amount of the variance

# check normality
qqnorm(resid(scaledtrimfplmm),
       main="Normal q-qplot scaledtrimfplmm")
qqnorm(resid(scaledtriminvfplmm),
       main="Normal q-qplot scaledtriminvfplmm")

# The second model is much closer to normality

# Plot residuals
plot(scaledtrimfplmm, resid(.) ~ fitted(.),
     main="Residuals scaledtrimfplmm")
plot(scaledtriminvfplmm, resid(.) ~ fitted(.),
     main="Residuals scaledtriminvfplmm")

# Outliers are gone, but residuals still appear to correlate with fitted values in
# the first model

qplot(resid(scaledtrimfplmm),
      main="Residuals scaledtrimfplmm")
qplot(resid(scaledtriminvfplmm),
      main="Residuals scaledtriminvfplmm")

# Both still appear relatively normally distributed, with the second model
# performing better

# Use z-scores for centering with no outlier trimming
scaledfplmm3 <- mixed(formula = RTzscore ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                  (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(scaledfplmm3)

anova(scaledfplmm3)

# Amount of variance accounted for by the model
cor(fitted(scaledfplmm3), goData$RT)^2

# Very low variance explained

# check normality
qqnorm(resid(scaledfplmm3),
       main="Normal q-qplot scaledfplmm3")

# Far from normality

# Plot residuals
plot(scaledfplmm3, resid(.) ~ fitted(.),
     main="Residuals scaledfplmm3")

# Correlation between residuals and fitted values seems to dissappear

qplot(resid(scaledfplmm3),
      main="Residuals scaledfplmm3")

# Residuals are quite assymnetric

# Use z-scores for centering with trimming
scaledfplmm4 <- mixed(formula = RTzscore ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                  (1+scaledNumForeperiod+condition+scaledNumOneBackFP||ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

isSingular(scaledfplmm4)

# Singular fit

# Amount of variance accounted for by the model
cor(fitted(scaledfplmm4), goData2$RT)^2

# check normality
qqnorm(resid(scaledfplmm4),
       main="Normal q-qplot scaledfplmm4")

# Plot residuals
plot(scaledfplmm4, resid(.) ~ fitted(.),
     main="Residuals scaledfplmm4")

qplot(resid(scaledfplmm4),
      main="Residuals scaledfplmm4")

# Only models with RT or inverse RT perform well

# Of those:
# Variance explained is higher with trimming
# Q-q plots are better with trimming
# Residuals are less correlated with fitted values with trimming
# Residuals are more normally distributed with trimming

# The model with inverse RT performs better than the one with RT

# Model comparisons
scaledtriminvfplmm2 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                         (1+scaledNumForeperiod+condition||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(scaledtriminvfplmm2, scaledtriminvfplmm, refit=FALSE)

scaledtriminvfplmm3 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                         (1+scaledNumForeperiod+scaledNumOneBackFP||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(scaledtriminvfplmm3, scaledtriminvfplmm, refit=FALSE)

scaledtriminvfplmm4 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                         (1+condition+scaledNumOneBackFP||ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")

anova(scaledtriminvfplmm4, scaledtriminvfplmm, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(scaledtriminvfplmm), goData2$RT)^2
cor(fitted(scaledtriminvfplmm2), goData2$invRT)^2
cor(fitted(scaledtriminvfplmm3), goData2$invRT)^2
cor(fitted(scaledtriminvfplmm4), goData2$invRT)^2

# Arrange by BIC
BIC(scaledtriminvfplmm, scaledtriminvfplmm2, scaledtriminvfplmm3, scaledtriminvfplmm4) %>%
  arrange(BIC)

# The fullest model explains most of the variance and BICs are very close between
# this and the model without a random-effect from Fpn-1. In line with the 
# "keep it maximal" strategy, we keep this full model

#==========================================================================================#
#==================================== 2. Model assessment =================================
#==========================================================================================#

#=============================== 2.1. FP and FP n-1 as numerical ==================================


#====================== 2.2. Using FP n-1 as categorical for emm comparisons ==========================

#============= 2.2.1. Using RT ================
sclTrimFpLmm6 <- buildmer(RT ~ scaledNumForeperiod * condition * oneBackFP + 
                            (1+scaledNumForeperiod*condition*oneBackFP|ID), 
                          data=goData2,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(sclTrimFpLmm6)
formula(sclTrimFpLmm6)


sclTrimFpLmm6 <- mixed(formula = RT ~  1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                         condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                         scaledNumForeperiod:oneBackFP:condition + (1 + condition | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

anova(sclTrimFpLmm6)

# Visualize random effects
dotplot(ranef(sclTrimFpLmm6, condVar = TRUE))

# emm
Fp_by_Previous=emtrends(sclTrimFpLmm6, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(sclTrimFpLmm6, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

# Same but without n-1 no-go trials
sclTrimFpLmm6 <- buildmer(RT ~ scaledNumForeperiod * condition * oneBackFP + 
                            (1+scaledNumForeperiod*condition*oneBackFP|ID), 
                          data=goData3,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(sclTrimFpLmm6)
formula(sclTrimFpLmm6)


sclTrimFpLmm6 <- mixed(formula = RT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                         condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                         scaledNumForeperiod:oneBackFP:condition + (1 + condition + scaledNumForeperiod | ID),
                       data=goData3,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

summary(sclTrimFpLmm6)
anova(sclTrimFpLmm6)

# emm
Fp_by_Previous=emtrends(sclTrimFpLmm6, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(sclTrimFpLmm6, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")


# Now without n-1 go trials
trimfplmm6 <- buildmer(RT ~ numForeperiod * condition * oneBackFP + 
                         (1+numForeperiod*condition*oneBackFP|ID), 
                       data=filter(goData2, oneBacktrialType == "no-go"),
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm6)
formula(trimfplmm6)

trimfplmm6 <- mixed(formula = RT ~ 1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
                      condition + numForeperiod:condition + oneBackFP:condition + 
                      numForeperiod:oneBackFP:condition + (1 + condition | ID),
                    data=filter(goData2, oneBacktrialType == "no-go"),
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

summary(trimfplmm6)
anova(trimfplmm6)

# emm
Fp_by_Previous=emtrends(trimfplmm6, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm6, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

#=============== 2.2.2. Using invRT =================
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
Fp_by_Previous=emtrends(triminvfplmm6, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(triminvfplmm6, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

#============= 2.2.3. Using logRT ================
sclTrimLogFpLmm6 <- buildmer(logRT ~ scaledNumForeperiod * condition * oneBackFP + 
                            (1+scaledNumForeperiod*condition*oneBackFP|ID), 
                          data=goData2,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(sclTrimLogFpLmm6)
formula(sclTrimLogFpLmm6)


sclTrimLogFpLmm6 <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                         condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                         scaledNumForeperiod:oneBackFP:condition + (1 + condition | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

anova(sclTrimLogFpLmm6)

# Visualize random effects
dotplot(ranef(sclTrimFpLmm6, condVar = TRUE))

# emm
Fp_by_Previous=emtrends(sclTrimFpLmm6, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(sclTrimFpLmm6, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

# Same but without n-1 no-go trials
sclTrimFpLmm6 <- buildmer(RT ~ scaledNumForeperiod * condition * oneBackFP + 
                            (1+scaledNumForeperiod*condition*oneBackFP|ID), 
                          data=goData3,
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(sclTrimFpLmm6)
formula(sclTrimFpLmm6)


sclTrimFpLmm6 <- mixed(formula = RT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                         condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                         scaledNumForeperiod:oneBackFP:condition + (1 + condition + scaledNumForeperiod | ID),
                       data=goData3,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

summary(sclTrimFpLmm6)
anova(sclTrimFpLmm6)

# emm
Fp_by_Previous=emtrends(sclTrimFpLmm6, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(sclTrimFpLmm6, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")


# Now without n-1 go trials
trimfplmm6 <- buildmer(RT ~ numForeperiod * condition * oneBackFP + 
                         (1+numForeperiod*condition*oneBackFP|ID), 
                       data=filter(goData2, oneBacktrialType == "no-go"),
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm6)
formula(trimfplmm6)

trimfplmm6 <- mixed(formula = RT ~ 1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
                      condition + numForeperiod:condition + oneBackFP:condition + 
                      numForeperiod:oneBackFP:condition + (1 + condition | ID),
                    data=filter(goData2, oneBacktrialType == "no-go"),
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

summary(trimfplmm6)
anova(trimfplmm6)

# emm
Fp_by_Previous=emtrends(trimfplmm6, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm6, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")


#========================= 2.3. Both FP and FP n-1 as categorical ===================================

# 2.3.1. Using invRT
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

# 2.3.2. Using RT
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
scaledinvfplmm2 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                     (1+scaledNumForeperiod+condition||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(scaledinvfplmm2, scaledinvfplmm1, refit=FALSE)

scaledinvfplmm3 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                     (1+scaledNumForeperiod+scaledNumOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(scaledinvfplmm3, scaledinvfplmm1, refit=FALSE)

scaledinvfplmm4 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*scaledNumOneBackFP + 
                     (1+condition+scaledNumOneBackFP||ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

anova(scaledinvfplmm4, scaledinvfplmm1, refit=FALSE)

# Amount of variance accounted for by the models
cor(fitted(scaledinvfplmm1), goData$RT)^2
cor(fitted(scaledinvfplmm2), goData$invRT)^2
cor(fitted(scaledinvfplmm3), goData$invRT)^2
cor(fitted(scaledinvfplmm4), goData$invRT)^2


