
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




#==========================================================================================#
#=========================== 2. Prepare model ==========================
#==========================================================================================#


fplmm <- mixed(formula = logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                 foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                 (1 + condition | ID),
               data=goData2,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = FALSE,
               method =  'KR',
               REML=TRUE,
               return = "merMod")


#==========================================================================================#
#==================================== 3. Model assessment =================================
#==========================================================================================#

#=========== 3.1. Omnibus anova ===============
anova(fplmm)


#Visualize random effects
dotplot(ranef(fplmm, condVar = TRUE))

# Pairwise comparisons by FP (estimates consecutive differences)

#========== 3.2. Two-way interactions ==========
# Compare difference between conditions for each level of FPn
fp_cond_emm <- emmeans(fplmm, ~ condition|foreperiod)
contrast(fp_cond_emm, interaction = c("pairwise"), adjust = "holm")

# Compare difference between consecutive levels of FPn-1 for each level of FPn
fp_onebackfp_emm <- emmeans(fplmm, ~ oneBackFP|foreperiod)
contrast(fp_onebackfp_emm, interaction = c("consec"), adjust = "holm")


#========== 3.3 Three-way interaction =========
# Multiple plots to visualize interactions
emmip(fplmm, condition ~ foreperiod|oneBackFP, CIs = TRUE,
      xlab = "FP",
      facelab = "label_both") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 16, hjust = 0.5)) +
  scale_color_manual(values = c("orange", "blue"))
ggsave("./Analysis/Plots/mixed_models_pairwise.png",
       width = 8.5,
       height = 5.7)

emmip(fplmm, oneBackFP ~ foreperiod|condition, CIs = TRUE,
      xlab = "FP",
      facelab = "label_both") +
  labs(title = "RT pairwise comparisons") +
  theme(plot.title = element_text(size = 16, hjust = 0.5))


# First analysis: compare two-way interactions (FPn x FPn-1) between conditions
# Compare FP linear fits between consecutive levels of FPn-1, separately for each condition
fpemm <- emmeans(fplmm, ~ foreperiod*oneBackFP|condition)
contrast(fpemm, interaction = c("poly", "consec", "pairwise"), adjust = "holm")

# Now compare magnitude of those interactions between conditions
fpemm <- emmeans(fplmm, ~ foreperiod*oneBackFP*condition)
contrast(fpemm, interaction = c("poly", "consec", "pairwise"), adjust = "holm")



# Second analysis: compare FP effect between conditions, separately for each FPn-1
# Test significance of FP linear fits separately for each level of FPn-1 and each condition
fpemm <- emmeans(fplmm, ~ foreperiod|oneBackFP*condition)
contrast(fpemm, interaction = c("poly"), adjust = "holm")


# Now test difference of FP fits between conditions, separately for each level of FPn-1
fpemm <- emmeans(fplmm, ~ foreperiod*condition|oneBackFP)
contrast(fpemm, interaction = c("poly", "pairwise"), adjust = "holm")




# Same but without no-go trials
fplmm_nogo <- buildmer(logRT ~ foreperiod * condition * oneBackFP + 
                             (1 + foreperiod * condition * oneBackFP|ID), 
                           data=goData3,
                           buildmerControl = buildmerControl(crit = "LRT",
                                                              family = gaussian(link = "identity"),
                                                              calc.anova = TRUE,
                                                              include = ~ foreperiod * condition * oneBackFP))

isSingular(fplmm_nogo)
formula(fplmm_nogo)


fplmm_nogo <- mixed(formula = logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                          foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
                          (1 + condition | ID),
                        data=goData3,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = FALSE,
                        method =  'KR',
                        REML=TRUE,
                        return = "merMod")

anova(fplmm_nogo)

# emm
# Test significance of FP linear fits separately for each level of FPn-1 and each condition
fpemm_nogo <- emmeans(fplmm_nogo, ~ foreperiod|oneBackFP*condition)
contrast(fpemm_nogo, interaction = c("poly"), adjust = "holm")


# Now test difference of FP fits between conditions, separately for each level of FPn-1
fpemm_nogo <- emmeans(fplmm_nogo, ~ foreperiod*condition|oneBackFP)
contrast(fpemm_nogo, interaction = c("poly", "pairwise"), adjust = "holm")




# Now without n-1 go trials
fplmm_nogo <- buildmer(logRT ~ foreperiod * condition * oneBackFP + 
                               (1 + foreperiod * condition * oneBackFP|ID), 
                             data=filter(goData2, oneBacktrialType == "no-go"),
                             buildmerControl = buildmerControl(crit = "LRT",
                                                               family = gaussian(link = "identity"),
                                                               calc.anova = TRUE,
                                                               include = ~ foreperiod * condition * oneBackFP))

isSingular(fplmm_nogo)
formula(fplmm_nogo)

fplmm_nogo <- mixed(formula = logRT ~ 1 + foreperiod + condition + foreperiod:condition + oneBackFP + 
                      foreperiod:oneBackFP + condition:oneBackFP + foreperiod:condition:oneBackFP + 
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


#================= 5.2. Model using prevFPLonger =======================
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

# 7.1. Only go trials (misses)

goaccglmer <- mixed(formula = error_result ~ 1 + foreperiod:condition + foreperiod + condition + 
                      (1 | ID),
                    data=dataAccGo,
                    family = binomial(link = "logit"),
                    control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                    progress = TRUE,
                    expand_re = FALSE,
                    method = "LRT")

anova(goaccglmer)


# 7.2. only no-go trials (false alarms)
nogoaccglmer <- mixed(formula = error_result ~ 1 + foreperiod:condition + foreperiod + condition + 
                        (1 | ID),
                      data=dataAccNoGo,
                      family = binomial(link = "logit"),
                      control = glmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5)),
                      progress = TRUE,
                      expand_re = FALSE,
                      method = "LRT")

anova(nogoaccglmer)

accnogo_fp_condition <- emtrends(nogoaccglmer, "condition", var = "scaledNumForeperiod")

#===============================================================================================#
#============================= 8. Analyses using exp(-foreperiod) ==============================
#===============================================================================================#
trimexpfplmm1 <- buildmer(RT ~ expFP * condition * oneBackFP + 
                            (1+expFP*condition*oneBackFP|ID), 
                          data=goData2,
                          buildmerControl = list(direction='backward',
                                                 crit='LRT',#ddf = "Satterthwaite",
                                                 family=gaussian(link = 'identity'),
                                                 calc.anova = TRUE))



#===============================================================================================#
#=================================== 9. Analyses with smaller N =================================
#===============================================================================================#

nPartsExp3 <- 28
nReps <- 1000
pValues <- vector(mode = "numeric", length = nReps)

# Sample 28 participants nReps times and fit mixed model for each sample, extracting p-value of comparison of interest
for(iRep in 1:nReps) {
  
  # Sample participants
  sampParts <- sample(levels(goData2$ID), nPartsExp3)
  sampData <- goData2 %>%
    filter(ID %in% sampParts)
  
  # Recompute scaled variables (which were de-centered by trimming)
  sampData$squaredNumForeperiod <-  (sampData$numForeperiod)^2
  sampData$squaredNumOneBackFP <- (sampData$numOneBackFP)^2
  sampData$scaledNumForeperiod <-  scale(sampData$numForeperiod, scale = FALSE)[,1]
  sampData$squaredScaledNumForeperiod <- (sampData$scaledNumForeperiod)^2
  sampData$scaledNumOneBackFP <- scale(sampData$numOneBackFP, scale = FALSE)[,1]
  
  
  # Fit model using same structure as before
  sampLmm <- mixed(formula = logRT ~ 1 + scaledNumForeperiod + oneBackFP + scaledNumForeperiod:oneBackFP + 
                          condition + scaledNumForeperiod:condition + oneBackFP:condition + 
                          scaledNumForeperiod:oneBackFP:condition + 
                          (1 + condition + scaledNumForeperiod | ID),
                        data=sampData,
                        control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = FALSE,
                        method =  'KR',
                        REML=TRUE,
                        return = "merMod")
  
  
  # Run anova and extract p-value
  sampAnova <- anova(sampLmm)
  #pValues[iRep] <- sampAnova['scaledNumForeperiod:oneBackFP:condition',]$`Pr(>F)`
  pValues[iRep] <- sampAnova['scaledNumForeperiod:condition',]$`Pr(>F)`
  
}

# Get p-value of comparison with full sample
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
stdAnova <- anova(trimlogfplmm)


# Plot sampled values against original values
jpeg("./Analysis/Plots/Smaller_sample_analyses/28_parts_1000_sims.jpg", width = 20, height = 12, units = "cm", res = 300)
hist(pValues, breaks = 100, col = "lightblue",
     main = "Exp. 1: p-Values for sample size = 28 (1000 simulations)")
#abline(v = stdAnova['scaledNumForeperiod:oneBackFP:condition',]$`Pr(>F)`, lwd = 2)
abline(v = stdAnova['scaledNumForeperiod:condition',]$`Pr(>F)`, lwd = 2)
abline(v = 0.001, lwd = 2, lty = 2)
text(x = 0.000055, y = 900, labels = "p-value for full sample")
text(x = 0.001, y = 900, labels = "p = .001")
dev.off()

