
#==========================================================================================#
# Fit the same mixed models as in Mixed_models.R, but using RT in s instead of 
# ms to improve convergeability

# For this we use centered and scaled numerical predictors
#==========================================================================================#

# Load necessary packages

# Data processing and plotting
library(magrittr)
library(tidyverse)
library(lattice)
library(gridExtra)
library(data.table)
library(knitr)

# Linear modeling
library(car)
library(janitor)

# Mixed effects modeling
library(afex)
library(emmeans)
library(lme4)
library(MuMIn)
library(buildmer)
library(broom.mixed)

# Post-hocs
library(marginaleffects)

# Save defaults
graphical_defaults <- par()
options_defaults <- options()

# emm options
emm_options(lmer.df = "satterthwaite", lmerTest.limit = 12000)

# Read data
source('./Analysis/Prepare_data_4.R')

#======================================= 0. Functions ======================================
hist_resid <- function(M,ptitle='Residuals') {
  d <- data.frame(resid=residuals(M)) 
  d  %>% ggplot(aes(x=resid)) + 
    geom_histogram(aes(y=after_stat(density)), bins=75, color='black', fill='grey') + 
    geom_density(color='darkred') + 
    ggtitle(ptitle) -> pl
  return(pl)
}

fitstats = function(M,mname='M') {
  QQ<-qqnorm(residuals(M), plot.it=FALSE)
  R2qq <- cor(QQ$x,QQ$y)^2
  dfqq = data.frame(stat='R2qq', V1=R2qq)
  r2tab <- r.squaredGLMM(M)  %>% 
    t  %>% as.data.frame  %>% rownames_to_column(var='stat')  %>% 
    rbind(.,dfqq)
  r2tab$stat = c("$R^2_m$","$R^2_c$",'$R^2_{qq}$' )
  colnames(r2tab) <- c('stat',mname)
  return(r2tab)
}

#==========================================================================================#
#====================================== 1. Prepare model ===================================
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

#=========================== 1.1. Find maximal converging structure =============================
# We use buildmer for this

# To choose the data transformation that leads to the optimal random effects structure, we fit models including only 
# random intercepts and compare R2 and residuals

fplmm <- buildmer(formula = RT ~ scaledNumForeperiod*condition*oneBackFP + 
                    (1|ID),
                  data = goData2,
                  buildmerControl = list(crit = "LRT",
                                 family = gaussian(link = "identity"),
                                 calc.anova = TRUE))
formula(fplmm)

fplmm <- mixed(formula = RT ~ scaledNumForeperiod*condition*oneBackFP + 
                  (1|ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                return = 'merMod',
                REML=TRUE)

summary(fplmm)

#============================== 1.2. Choose dependent variable ================================

# Fit models with RT, inverse RT, and logRT without trimming
fplmm1 <- mixed(formula = RT ~ scaledNumForeperiod*condition*oneBackFP + 
                  (1|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")


# Now we run the same model with inverse RT and logRT as outcomes
invfplmm1 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*oneBackFP + 
                     (1|ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

logfplmm1 <- mixed(formula = logRT ~ scaledNumForeperiod*condition*oneBackFP + 
                     (1|ID),
                   data = goData,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")


# Amount of variance accounted for by the model
cor(fitted(fplmm1), goData$RT)^2
cor(fitted(invfplmm1), goData$invRT)^2
cor(fitted(logfplmm1), goData$logRT)^2

# The last model explains a larger amount of the variance 

# Check normality of residuals
qqnorm(resid(fplmm1),
       main="Normal q-qplot fplmm1")
qqnorm(resid(invfplmm1),
       main="Normal q-qplot invfplmm1")
qqnorm(resid(logfplmm1),
       main="Normal q-qplot logfplmm1")

# Both models show considerable departures from normality

# Plot residuals
plot(fplmm1, resid(.) ~ fitted(.),
     main="Residuals fplmm1")
plot(invfplmm1, resid(.) ~ fitted(.),
     main="Residuals invfplmm1")
plot(logfplmm1, resid(.) ~ fitted(.),
     main="Residuals logfplmm1")

# It appears that residuals correlate somewhat with fitted values; there are also outliers


# All have a considerable departure from normality, although less so for the log model

# Fit models after outlier trimming
trimfplmm1 <- mixed(formula = RT ~ scaledNumForeperiod*condition*oneBackFP + 
                      (1|ID),
                data = goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE,
                return = "merMod")

triminvfplmm1 <- mixed(formula = invRT ~ scaledNumForeperiod*condition*oneBackFP + 
                         (1|ID),
                   data = goData2,
                   control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                   progress = TRUE,
                   expand_re = TRUE,
                   method =  'S',
                   REML=TRUE,
                   return = "merMod")

trimlogfplmm1 <- mixed(formula = logRT ~ scaledNumForeperiod*condition*oneBackFP + 
                         (1|ID),
                       data = goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'S',
                       REML=TRUE,
                       return = "merMod")


isSingular(trimfplmm1)
isSingular(triminvfplmm1)
isSingular(trimlogfplmm1)

# Amount of variance accounted for by the model
cor(fitted(trimfplmm1), goData2$RT)^2
cor(fitted(triminvfplmm1), goData2$invRT)^2
cor(fitted(trimlogfplmm1), goData2$logRT)^2

# Again, the first plot accounts for a larger amount of the variance

# check normality
qqnorm(resid(trimfplmm1),
       main="Normal q-qplot trimfplmm1")
qqnorm(resid(triminvfplmm1),
       main="Normal q-qplot triminvfplmm1")
qqnorm(resid(trimlogfplmm1),
       main="Normal q-qplot trimlogfplmm1")

# The second model is much closer to normality

# Plot residuals
plot(trimfplmm1, resid(.) ~ fitted(.),
     main="Residuals trimfplmm1")
plot(triminvfplmm1, resid(.) ~ fitted(.),
     main="Residuals triminvfplmm1")
plot(trimlogfplmm1, resid(.) ~ fitted(.),
     main="Residuals trimlogfplmm1")

# Outliers are gone, but residuals still appear to correlate with fitted values in
# the first model

# LIttle difference after trimming, although outliers are rarer; logRT still appears to perform better
grid.arrange(hist_resid(fplmm1, 'RT'),
             hist_resid(invfplmm1, '1/RT'),
             hist_resid(logfplmm1, 'logRT'),
             hist_resid(trimfplmm1, 'trimmed RT'),
             hist_resid(triminvfplmm1, 'trimmed 1/RT'),
             hist_resid(trimlogfplmm1, 'trimmed logRT'),
             ncol=2,
             as.table=FALSE)

# The same pattern is apparent as in the model with no trimming

R2table <- fitstats(fplmm1, 'RT') %>%
  plyr::join(., fitstats(invfplmm1, '1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(logfplmm1, 'log(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimfplmm1, 'trim RT'), by='stat') %>%
  plyr::join(., fitstats(triminvfplmm1, 'trim 1/(RT)'), by='stat') %>%
  plyr::join(., fitstats(trimlogfplmm1, 'trim log(RT)'), by='stat') %>%
  kable(digits=4)


#==========================================================================================#
#==================================== 3. Model assessment =================================
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

#========== 3.2.1. Using RT ===============
trimfplmm <- buildmer(RT ~ scaledNumForeperiod * condition * oneBackFP + 
                         (1+scaledNumForeperiod*condition*oneBackFP|ID), 
                       data=goData2,
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm)
formula(trimfplmm)


trimfplmm <- mixed(formula = RT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                     condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                     scaledNumForeperiod:oneBackFP:condition + (1 + condition + scaledNumForeperiod | ID),
                    data=goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

anova(trimfplmm)

#Visualize random effects
dotplot(ranef(trimfplmm, condVar = TRUE))

# emm
Fp_by_Previous=emtrends(trimfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

# Same but without n-1 no-go trials
trimfplmm <- buildmer(RT ~ scaledNumForeperiod * condition * oneBackFP + 
                         (1 + scaledNumForeperiod * condition * oneBackFP|ID), 
                       data=goData3,
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm)
formula(trimfplmm)


trimfplmm <- mixed(formula = RT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                     condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                     scaledNumForeperiod:oneBackFP:condition + 
                     (1 + condition + scaledNumForeperiod | ID),
                    data=goData3,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

summary(trimfplmm)
anova(trimfplmm)

# emm
Fp_by_Previous=emtrends(trimfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")


# Now without n-1 go trials
trimfplmm <- buildmer(RT ~ numForeperiod * condition * oneBackFP + 
                         (1+numForeperiod*condition*oneBackFP|ID), 
                       data=filter(goData2, oneBacktrialType == "no-go"),
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimfplmm)
formula(trimfplmm)

trimfplmm <- mixed(formula = RT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                     condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                     scaledNumForeperiod:oneBackFP:condition + 
                     (1 + condition + scaledNumForeperiod | ID),
                    data=filter(goData2, oneBacktrialType == "no-go"),
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

summary(trimfplmm)
anova(trimfplmm)

# emm
Fp_by_Previous=emtrends(trimfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

#========== 3.2.2. Using invRT ===============
triminvfplmm <- buildmer(invRT ~ numForeperiod * condition * oneBackFP + 
                           (1+numForeperiod*condition*oneBackFP|ID), 
                         data=goData2,
                         buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(triminvfplmm)
buildform <- formula(triminvfplmm)


triminvfplmm <- mixed(formula = invRT ~  1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
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
Fp_by_Previous=emtrends(triminvfplmm, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(triminvfplmm, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

#========== 3.2.3. Using logRT ===============
trimlogfplmm <- buildmer(logRT ~ scaledNumForeperiod * condition * oneBackFP + 
                         (1+scaledNumForeperiod*condition*oneBackFP|ID), 
                       data=goData2,
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimlogfplmm)
formula(trimlogfplmm)


trimlogfplmm <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                     condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                     scaledNumForeperiod:oneBackFP:condition + 
                     (1 + condition + scaledNumForeperiod | ID),
                    data=goData2,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = FALSE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

anova(trimlogfplmm)


#Visualize random effects
dotplot(ranef(trimlogfplmm, condVar = TRUE))

# Single slopes tests
fp_by_condition <- slopes(trimlogfplmm, by = "condition", variables = "scaledNumForeperiod",
                          p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ condition, var="scaledNumForeperiod")) # equivalent to slopes

fp_by_oneback <- slopes(trimlogfplmm, by = "oneBackFP", variables = "scaledNumForeperiod",
                        p_adjust = "holm")

test(emtrends(trimlogfplmm, ~ oneBackFP, var="scaledNumForeperiod")) # equivalent to slopes

threeway_int <- slopes(trimlogfplmm, by = c("oneBackFP", "condition"), variables = "scaledNumForeperiod",
                        p_adjust = "holm")




# Pairwise comparisons
fp_by_condition_comp <- emtrends(trimlogfplmm, "condition", var = "scaledNumForeperiod")
fp_by_condition_comp
update(pairs(fp_by_condition_comp), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimlogfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

threeway_int_comp = emtrends(trimlogfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
threeway_int_comp
update(pairs(threeway_int_comp), by = NULL, adjust = "holm")

# Same but without n-1 no-go trials
trimlogfplmmGo <- buildmer(logRT ~ scaledNumForeperiod * condition * oneBackFP + 
                         (1 + scaledNumForeperiod * condition * oneBackFP|ID), 
                       data=goData3,
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimlogfplmmGo)
formula(trimlogfplmmGo)


trimlogfplmmGo <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                        condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                        scaledNumForeperiod:oneBackFP:condition +
                        (1 + condition + scaledNumForeperiod | ID),
                    data=goData3,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = FALSE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

anova(trimlogfplmmGo)

# emm
Fp_by_Previous=emtrends(trimlogfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")


threeway_int_comp_Go = emtrends(trimlogfplmmGo, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
threeway_int_comp_Go
update(pairs(threeway_int_comp_Go), by = NULL, adjust = "holm")

# Now without n-1 go trials
trimlogfplmmNoGo <- buildmer(logRT ~ scaledNumForeperiod * condition * oneBackFP + 
                         (1 + scaledNumForeperiod * condition * oneBackFP|ID), 
                       data=filter(goData2, oneBacktrialType == "no-go"),
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

isSingular(trimlogfplmmNoGo)
formula(trimlogfplmmNoGo)

trimlogfplmmNoGo <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + condition + oneBackFP + scaledNumForeperiod:condition + 
                     scaledNumForeperiod:oneBackFP + scaledNumForeperiod:condition:oneBackFP + 
                     (1 + condition | ID),
                    data=filter(goData2, oneBacktrialType == "no-go"),
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'KR',
                    REML=TRUE,
                    return = "merMod")

summary(trimlogfplmmNoGo)
anova(trimlogfplmmNoGo)

# emm
Fp_by_Previous=emtrends(trimlogfplmm, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimlogfplmm, c("condition", "oneBackFP"), var = "scaledNumForeperiod")
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


#=============== Sanity check: model comparisons without trimming =================
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

#=================== 3.3. Use logRT to get priors for next experiment =====================

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



#============================================================================================================#
#===================================== 4. Sequential effects with n-1 trial type =============================
#============================================================================================================#

#============================ 4.1. FP n-1 and trial type in combined variable ============================

#======= 4.1.1. Continuous variables ========

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

#=========== 4.1.2. Categorical variables ============

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


#================================ 4.2. FP n-1 and trial type as separate variables =======================

#=========== 4.2.1. Continuous variables ============

# RT
trimfplmmtrialtype3 <- buildmer(RT ~ numForeperiod * condition * numOneBackFP * oneBacktrialType + 
                                  (1 + numForeperiod * condition * numOneBackFP * oneBacktrialType |ID), 
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
                                     (1 + numForeperiod * condition * numOneBackFP * oneBacktrialType |ID), 
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

#=========== 4.2.2. FP n-1 as categorical and FP n as continuous ============

# 4.2.2.1. RT
trimfplmmtrialtype4 <- buildmer(RT ~ numForeperiod * condition * oneBackFP * oneBacktrialType + 
                                  (1 + numForeperiod * condition * oneBackFP * oneBacktrialType |ID), 
                                data=goData2,
                                buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimfplmmtrialtype4)

isSingular(trimfplmmtrialtype4)

trimfplmmtrialtype4 <- mixed(RT ~  1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
                               condition + oneBacktrialType + numForeperiod:condition + 
                               condition:oneBacktrialType + numForeperiod:oneBacktrialType + 
                               oneBackFP:oneBacktrialType + oneBackFP:condition + numForeperiod:oneBackFP:condition +
                               numForeperiod:oneBackFP:condition:oneBacktrialType +
                               (1 + condition + oneBacktrialType | ID),
                             data=goData2,
                             control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                             progress = TRUE,
                             expand_re = TRUE,
                             method =  'KR',
                             REML=TRUE,
                             return = "merMod",
                             check_contrasts = FALSE)

summary(trimfplmmtrialtype4)

anova(trimfplmmtrialtype4)

# emm
Fp_by_Previous=emtrends(trimfplmmtrialtype4, "oneBackFP", var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimfplmmtrialtype4, c("condition", "oneBackFP"), var = "numForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

# 4.2.2.2. 1/RT
triminvfplmmtrialtype4 <- buildmer(invRT ~ numForeperiod * condition * numOneBackFP * oneBacktrialType + 
                                     (1 + numForeperiod * condition * numOneBackFPGo * oneBacktrialType |ID), 
                                   data=goData2,
                                   buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(triminvfplmmtrialtype4)

isSingular(triminvfplmmtrialtype4)

triminvfplmmtrialtype4 <- mixed(invRT ~ 1 + numForeperiod + numOneBackFP + numForeperiod:numOneBackFP + 
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

summary(triminvfplmmtrialtype4)

anova(triminvfplmmtrialtype4)

# 4.2.2.3. logRT
trimlogfplmmtrialtype <- buildmer(logRT ~ scaledNumForeperiod * condition * oneBacktrialType + 
                                  (1 + scaledNumForeperiod * condition * oneBacktrialType |ID), 
                                data=goData2,
                                buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimlogfplmmtrialtype)

isSingular(trimlogfplmmtrialtype)

trimlogfplmmtrialtype <- mixed(logRT ~ 1 + scaledNumForeperiod + oneBacktrialType + condition + 
                                 scaledNumForeperiod:condition + oneBacktrialType:condition + 
                                 scaledNumForeperiod:oneBacktrialType + scaledNumForeperiod:oneBacktrialType:condition + 
                                 (1 + condition + oneBacktrialType + scaledNumForeperiod + 
                                    condition:oneBacktrialType | ID),
                             data=goData2,
                             control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                             progress = TRUE,
                             expand_re = FALSE,
                             method =  'KR',
                             REML=TRUE,
                             return = "merMod")

summary(trimlogfplmmtrialtype)

anova(trimlogfplmmtrialtype)

# Pairwise comparisons
cond_trialtype <- emmeans(trimlogfplmmtrialtype, ~ oneBacktrialType | condition)
pairs(cond_trialtype)

fp_trialtype <- emtrends(trimlogfplmmtrialtype, "oneBacktrialType", var = "scaledNumForeperiod")
pairs(fp_trialtype)

Fp_by_Previous=emtrends(trimlogfplmmtrialtype, "oneBackFP", var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "holm")

Fp_by_Previous=emtrends(trimlogfplmmtrialtype, c("condition", "oneBacktrialType"), var = "scaledNumForeperiod")
Fp_by_Previous
update(pairs(Fp_by_Previous), by = NULL, adjust = "none")

three_way_trialtype = emtrends(trimlogfplmmtrialtype, ~ condition | oneBacktrialType, var = "scaledNumForeperiod")
con <- pairs(three_way_trialtype, simple = "condition")
contrast(con, "pairwise", by = NULL)


#=========== 4.2.3. Categorical variables ============

# 4.2.3.1. RT
trimfplmmtrialtype5 <- buildmer(RT ~ foreperiod * condition * oneBackFP * oneBacktrialType + 
                                  (1+foreperiod * condition * oneBackFP * oneBacktrialType|ID), 
                                data=goData2,
                                buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite"))

formula(trimfplmmtrialtype5)

isSingular(trimfplmmtrialtype5)

trimfplmmtrialtype5 <- mixed(RT ~ 1 + foreperiod + oneBackFP + foreperiod:oneBackFP + condition + 
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

trimfplmmtrialtype5 <- mixed(RT ~ foreperiod * condition * oneBackFP * oneBacktrialType +
                               (1 + condition + oneBacktrialType + condition:oneBacktrialType | 
                                  ID),
                             data=goData2,
                             control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                             progress = TRUE,
                             expand_re = TRUE,
                             method =  'KR',
                             REML=TRUE,
                             return = "merMod")

summary(trimfplmmtrialtype5)

anova(trimfplmmtrialtype5)

#=========== 4.3. By condition ===============
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

#============ 4.4. by trial type ==================

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


#================= 8.4. Model using prevFPLonger =======================
trimprevfplonglmm1 <- buildmer(RT ~ foreperiod * condition * prevFPLonger + 
                                 (1 + foreperiod * condition * prevFPLonger|ID), 
                               data = goData2,
                               buildmerControl = list(#direction='backward',
                                 crit='LRT',#ddf = "Satterthwaite",
                                 family=gaussian(link = 'identity'),
                                 calc.anova = TRUE))

formula(trimprevfplonglmm1)

# Model using prevFPLonger
trimprevfplonglmm1 <- mixed(formula = RT ~ 1 + foreperiod + prevFPLonger + condition + foreperiod:prevFPLonger + 
                              prevFPLonger:condition + foreperiod:condition + foreperiod:prevFPLonger:condition + 
                              (1 + condition | ID),
                            data = goData3,
                            control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                            progress = TRUE,
                            expand_re = TRUE,
                            method =  'S',
                            REML=TRUE,
                            return = "merMod")

isSingular(trimprevfplonglmm1)

summary(trimprevfplonglmm1)
anova(trimprevfplonglmm1)
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
dataAcc <- dataAcc %>%
  filter(ID != "002")

dataAccGo <- dataAccGo %>%
  filter(ID != "002")

dataAccNoGo <- dataAccNoGo %>%
  filter(ID != "002")

# 7.1. all trials
fpaccglmer <- buildmer(error_result ~ scaledNumForeperiod * condition + 
                        (1+scaledNumForeperiod*condition|ID), 
                      data=dataAcc,
                      family = binomial(link = "logit"),
                      buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite",
                                                        include = ~ scaledNumForeperiod:condition))

isSingular(fpaccglmer)
formula(fpaccglmer)

fpaccglmer <- mixed(formula = error_result ~ 1 + condition + scaledNumForeperiod + 
                         scaledNumForeperiod:condition + (1 + scaledNumForeperiod | ID),
                       data=dataAcc,
                       family = binomial(link = "logit"),
                       control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                       progress = TRUE,
                       expand_re = FALSE,
                       method = "LRT")


anova(fpaccglmer)

test(emtrends(fpaccglmer, "condition", var = "scaledNumForeperiod"))

slopes(fpaccglmer, by = "condition", variables = "scaledNumForeperiod")

marginalmeans(fpaccglmer, type = "response", variables = c("condition", "scaledNumForeperiod"))


# 7.2. only go trials
goaccglmer <- buildmer(error_result ~ scaledNumForeperiod * condition + 
                            (1+scaledNumForeperiod*condition|ID), 
                          data=dataAccGo,
                          family = binomial(link = "logit"),
                          buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite",
                                                            include = ~ scaledNumForeperiod:condition))


isSingular(goaccglmer)
formula(goaccglmer)

goaccglmer <- mixed(formula = error_result ~ 1 + scaledNumForeperiod + condition + scaledNumForeperiod:condition + 
                      (1 + condition | ID),
                    data=dataAccGo,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                    progress = TRUE,
                    expand_re = FALSE,
                    method = "LRT")

anova(goaccglmer)

accgo_fp_condition <- emtrends(goaccglmer, "condition", var = "scaledNumForeperiod")


# 7.2. only no-go trials
nogoaccglmer <- buildmer(error_result ~ scaledNumForeperiod * condition + 
                         (1+scaledNumForeperiod*condition|ID), 
                       data=dataAccNoGo,
                       family = binomial(link = "logit"),
                       buildmerControl = buildmerControl(calc.anova = TRUE, ddf = "Satterthwaite",
                                                         include = ~ scaledNumForeperiod:condition))

isSingular(nogoaccglmer)
formula(nogoaccglmer)

nogoaccglmer <- mixed(formula = error_result ~ 1 + scaledNumForeperiod + condition + scaledNumForeperiod:condition + 
                        (1 + scaledNumForeperiod | ID),
                    data=dataAccNoGo,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                    progress = TRUE,
                    expand_re = FALSE,
                    method = "LRT")

anova(nogoaccglmer)

accnogo_fp_condition <- emtrends(nogoaccglmer, "condition", var = "scaledNumForeperiod")
