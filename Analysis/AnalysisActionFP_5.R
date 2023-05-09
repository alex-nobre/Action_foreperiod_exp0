

#================================================================================#
# Changes
# Manually set contrasts for anova
#================================================================================#

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
library(codingMatrices)
library(performance)


# Functions for raincloud plots
source('./Analysis/Prepare_data_4.R')


#==========================================================================================#
#======================================= 0. Data quality ===================================
#==========================================================================================#
testnormality = function(dfr) return(shapiro.test(dfr$invRT)$p.value)
p = as.vector(by(goData, goData$ID, testnormality))
names(p) = levels(goData$ID)
names(p[p < 0.05])


qqmath(~RT|ID, data=goData)
qqmath(~invRT|ID, data=goData)


#==========================================================================================#
#======================================= 1. Descriptives ==================================
#==========================================================================================#
#========== 1.1. Learning/fatigue effects ==============

# 1.1.1. RTs across blocks (no condition)
ggplot(data=summaryData,
       aes(x=block,
           y=meanRT))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=1))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

ggplot(data=summaryData,
                   aes(x=block,
                       y=meanRT,
                       color=condition))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=condition))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())+
  scale_color_manual(values=c('orange','blue'))

# 1.1.2. RTs across blocks by counterbalancing order
ggplot(data=summaryData2,
       aes(x=block,
           y=meanRT,
           color=counterbalance))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line', linewidth=1, aes(group=counterbalance))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('blue','orange')) +
  labs(title='RT by block split by counterbalancing order')

# Check for influence of external fixation duration
ggplot(data=filter(goData,condition=='external'),
       aes(x=extFixationDuration,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)

# Check for influence of latency of action key press on RT
ggplot(data=filter(goData,condition=='action'),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  facet_wrap(~foreperiod)

# Accuracy
summaryAccData <- data %>%
  group_by(ID,foreperiod,condition) %>%
  summarise(meanAcc=mean(Acc)) %>%
  ungroup()

ggplot(data=summaryAccData,
       aes(x=foreperiod,
           y=meanAcc,
           color=condition)) +
  geom_boxplot()


#==========================================================================================#
#==================================== 2. Basic models ======================================
#==========================================================================================#

# Set constrasts for variables used in ANOVAs
contrasts(summaryData$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryData$condition) <- c(-1/2, 1/2)
contrasts(summaryData$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)

contrasts(summaryData2$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryData2$condition) <- c(-1/2, 1/2)
contrasts(summaryData2$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)

#==================== 2.1. FP x RT by condition ======================
lines_by_condition <- ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") + 
  labs(title = "RT by condition",
       x = "Foreperiod",
       y = "Mean RT") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"))

ggplot2::ggsave("./Analysis/Plots/RT_by_condition.png",
                lines_by_condition)

fpAnova <- aov_ez(id = "ID",
       dv = "meanRT",
       data = summaryData2,
       within = c("foreperiod", "condition"))

check_sphericity(fpAnova)

is_norm <- check_normality(fpAnova)

plot(is_norm)

plot(is_norm, type = "qq")

nice(fpAnova,
     correction='none')
summary(fpAnova)


fpregression <- lm(meanRT ~ condition * foreperiod, data = summaryData)
summary(fpregression)
anova(fpregression)

logfpregression <- lm(meanRT ~ condition * logFP, data = summaryData)
anova(logfpregression)

# fpEmmeans <- emmeans(fpAnova,
#                      pairwise ~ condition|foreperiod,
#                      adjust = 'none')


fpEmmeans <- emmeans(fpAnova,
                     pairwise ~ foreperiod|condition,
                     adjust = 'bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('poly'),
                               adjust='bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('consec'),
                               adjust='bonferroni')

#======================= 2.2. Sequential effects ============================
# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                  dv = "meanRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition", "oneBackFP"))

seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanInvRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBackFP"))

nice(seqEffAnova,
     correction='none')
                  
seqEFffregression <- lm(meanRT ~ foreperiod * oneBackFP * condition, 
                   data = summaryData)
summary(seqEFffregression)
anova(seqEFffregression)

logseqEffregression <- lm(meanRT~logFP*logoneBackFP*condition,
                          data=summaryData) 
anova(logseqEffregression)

# 2.2.2. Lm with difference between the durations of FPn and FPn-1 as regressor
fpDiffRegression <- lm(meanRT ~ foreperiod * condition * oneBackFPDiff,
                       data = summaryData)
summary(fpDiffRegression)
anova(fpDiffRegression)

# 2.2.3. Anova for FP n-2
ntworegression <- lm(meanRT ~ foreperiod * oneBackFP * twoBackFP, 
                     data = summaryData)
summary(ntworegression)
anova(ntworegression)
Anova(ntworegression, type = "II")

#================= 2.3. N-1 trial type, condition and sequential effects =================
# 2.3.1. Anova with foreperiod, condition and n-1 trial type
trialTypeAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData,
                      within = c("foreperiod", "condition", "oneBacktrialType"))

nice(trialTypeAnova,
     correction='none')


# 2.3.2. Anova with foreperiod, n-1 trial type and sequential effects
seqEfftrialTypeAnova <- aov_ez(id = "ID",
                         dv = "meanRT",
                         data = summaryData,
                         within = c("foreperiod", "oneBackFP", "oneBacktrialType"))

nice(seqEfftrialTypeAnova,
     correction='none')

# 2.3.3. Anova with foreperiod, condition, n-1 trial type and sequential effects
fullAnova <- aov_ez(id = "ID",
                    dv = "meanRT",
                    data = summaryData,
                    within = c("foreperiod", "condition",
                               "oneBackFP", "oneBacktrialType"))

nice(fullAnova,
     correction='none')

# Only go-trials
summaryDataGo <- summaryData %>%
  filter(oneBacktrialType=="go")

goAnova <- aov_ez(id = "ID",
                         dv = "meanRT",
                         data = summaryDataGo,
                         #within = c("foreperiod", "condition", "oneBackFP"))
                  within = c("foreperiod", "condition", "oneBackFP"))

nice(goAnova,
     correction='none')

# Sequential effects (separated by condition)
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = oneBackFP)) +
  #stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange', 'green'))


#==================== 2.4. Foreperiod, condition and block ======================
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)

#==========================================================================================#
#====================================== 3. Mixed models ====================================
#==========================================================================================#

#=========================== 3.1. Foreperiod, condition and sequential effects =============================

# ============ 3.1.1. n-1 sequential effect =============
# 3.1.1.1 Fullest model
fplmm1 <- mixed(formula = RT ~ numForeperiod*condition*numOneBackFP + 
                 (1+numForeperiod*condition*numOneBackFP|ID),
               data = goData,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'S',
               REML=TRUE)

fplmm1v2 <- lmer(RT ~ numForeperiod*condition*numOneBackFP + 
                   (1+numForeperiod*condition*numOneBackFP||ID),
                 data=goData,
                 contrasts = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 REML=TRUE)

fplmm1v3 <- lmer(RT ~ numForeperiod*condition*numOneBackFP + 
                   (1+numForeperiod+condition+numOneBackFP||ID),
                 data=goData,
                 contrasts = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 REML=TRUE)


# ============ 3.1.2. difference between the durations of FPn and FPn-1 as regressor ===========
# 3.1.1.1 Fullest model
fpDifflmm1 <- mixed(formula = RT ~ foreperiod * condition * oneBackFPDiff + 
                  (1+foreperiod*condition*oneBackFPDiff|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)

fpDifflmm2 <- mixed(formula = RT ~ foreperiod * condition * oneBackFPDiff + 
                      (1+foreperiod*condition*oneBackFPDiff||ID),
                    data = goData,
                    control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                    progress = TRUE,
                    expand_re = TRUE,
                    method =  'S',
                    REML=TRUE)

#======================= 3.2. Foreperiod, condition and n-1 trial type =====================
# 3.2.1.1. Fullest model
oneBacktrialTypelmm1 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType + 
                                (1+foreperiod*condition*oneBacktrialType|ID),
                              data = goData,
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'S',
                              REML=TRUE)

oneBacktrialTypelmm2 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType + 
                                (1+foreperiod*condition*oneBacktrialType||ID),
                              data = goData,
                              control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                              progress = TRUE,
                              expand_re = TRUE,
                              method =  'S',
                              REML=TRUE)

# 3.2.2
# 3.2.2.1. Remove correlations of mixed part
oneBacktrialTypeSeqlmm2 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType * oneBackFP + 
                                   (1+foreperiod*condition*oneBacktrialType*oneBackFP||ID),
                                 data = goData,
                                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                 progress = TRUE,
                                 expand_re = TRUE,
                                 method =  'S',
                                 REML=TRUE)

# 3.2.2.3. Remove interactions of mixed part
oneBacktrialTypeSeqlmm3 <- mixed(formula = RT ~ foreperiod * condition * oneBacktrialType * oneBackFP + 
                                   (1+foreperiod+condition+oneBacktrialType+oneBackFP||ID),
                                 data = goData,
                                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                                 progress = TRUE,
                                 expand_re = TRUE,
                                 method =  'S',
                                 REML=TRUE)
#====================== 3.4. Foreperiod, counterbalance and block ========================
# 3.4.1. Fullest model
blockFplmm1 <- mixed(formula = RT ~ foreperiod * counterbalance * block + 
                  (1+foreperiod*counterbalance*block|ID),
                data = goData,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'S',
                REML=TRUE)

# 3.4.2. Remove correlations of mixed part
blockFplmm2 <- mixed(formula = RT ~ foreperiod * counterbalance * block + 
                       (1+foreperiod*counterbalance*block||ID),
                     data = goData,
                     control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                     progress = TRUE,
                     expand_re = TRUE,
                     method =  'S',
                     REML=TRUE)

# 3.4.3. Remove interactions of mixed part
blockFplmm3 <- mixed(formula = RT ~ foreperiod * counterbalance * block + 
                       (1+foreperiod+counterbalance+block||ID),
                     data = goData,
                     control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                     progress = TRUE,
                     expand_re = TRUE,
                     method =  'S',
                     REML=TRUE)
summary(blockFplmm3)$varcor
summary(blockFplmm3)
anova(blockFplmm3)


# 3.4.4. Plot
ggplot(data = goData,
       aes(x = block,
           y = RT,
           color = counterbalance)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group = counterbalance)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("blue","orange")) +
  facet_wrap(~foreperiod)

# 3.5. Learning 
#==========================================================================================#
#============================ 4. Two-stage regression ======================================
#==========================================================================================#

#=============== 4.1. Foreperiod, condition and sequential effects ====================
# regSummaryData <- summaryData %>%
#   mutate(foreperiod = as.numeric(foreperiod),
#          oneBackFP = as.numeric(oneBackFP))

# Fit global model to extract coefficients
globallm <- lm(meanRT ~ foreperiod * condition * oneBackFP,
               data = summaryData)

# Create matrix for coefficients (one column by coefficient, one line by participant)
idList <- unique(goData$ID)
coefmatrix <- data.frame(matrix(nrow = length(idList), 
                                ncol = length(names(globallm$coefficients))))
colnames(coefmatrix) <- names(globallm$coefficients)


# Fit lms and store coefficients in data frame
for(sub in 1:length(idList)) {
  subID <- idList[sub]
  subData <- goData %>%
    filter(ID == subID)
  sublm <- lm(RT ~ foreperiod * condition * oneBackFP,
              data = subData)
  subcoefs <- sublm$coefficients
  coefmatrix[sub,] <- subcoefs
}

# FPn 1.2
fp12ttest <- t.test(coefmatrix$foreperiod1.2,
                    alternative = "less",
                    mu = 0)

# FPn 1.8
fp18ttest <- t.test(coefmatrix$foreperiod1.8,
                    alternative = "less",
                    mu = 0)

# Condition (action)
actionttest <- t.test(coefmatrix$conditionaction,
                    alternative = "two.sided",
                    mu = 0)

# FPn-1 1.2
oneBackFP12ttest <- t.test(coefmatrix$oneBackFP1.2,
                          alternative='two.sided',
                          mu=0)

# FPn-1 1.8
oneBackFP18ttest <- t.test(coefmatrix$oneBackFP1.8,
                          alternative='two.sided',
                          mu=0)

# Two-way interaction FPn 1.2 x condition (action)
fp12actionttest <- t.test(coefmatrix$"foreperiod1.2:conditionaction",
                          alternative='two.sided',
                          mu=0)

# Two-way interaction FPn 1.8 x condition (action)
fp18actionttest <- t.test(coefmatrix$"foreperiod1.8:conditionaction",
                          alternative='two.sided',
                          mu=0)

# Two-way interaction FPn 1.2 x FPn-1 1.2
fp12oneBackFP12ttest <- t.test(coefmatrix$"foreperiod1.2:oneBackFP1.2",
                          alternative='two.sided',
                          mu=0)

# Two-way interaction FPn 1.8 x FPn-1 1.2
fp18oneBackFP12ttest <- t.test(coefmatrix$"foreperiod1.8:oneBackFP1.2",
                               alternative='two.sided',
                               mu=0)

# Two-way interaction FPn 1.2 x FPn-1 1.8
fp12oneBackFP18ttest <- t.test(coefmatrix$"foreperiod1.2:oneBackFP1.8",
                               alternative='two.sided',
                               mu=0)

# Two-way interaction FPn 1.8 x FPn-1 1.8
fp18oneBackFP18ttest <- t.test(coefmatrix$"foreperiod1.8:oneBackFP1.8",
                               alternative='two.sided',
                               mu=0)

# Two-way interaction condition (action) x FPn-1 1.2
actiononeBackFP12ttest <- t.test(coefmatrix$"conditionaction:oneBackFP1.2",
                               alternative='two.sided',
                               mu=0)

# Two-way interaction condition (action) x FPn-1 1.8
actiononeBackFP18ttest <- t.test(coefmatrix$"conditionaction:oneBackFP1.8",
                                 alternative='two.sided',
                                 mu=0)

# Three-way interaction FPn 1.2 x condition (action) x FPn-1 1.2
fp12actiononeBackFP12ttest <- t.test(coefmatrix$"foreperiod1.2:conditionaction:oneBackFP1.2",
                                 alternative='two.sided',
                                 mu=0)

# Three-way interaction FPn 1.8 x condition (action) x FPn-1 1.2
fp18actiononeBackFP12ttest <- t.test(coefmatrix$"foreperiod1.8:conditionaction:oneBackFP1.2",
                                     alternative='two.sided',
                                     mu=0)

# Three-way interaction FPn 1.2 x condition (action) x FPn-1 1.8
fp12actiononeBackFP18ttest <- t.test(coefmatrix$"foreperiod1.2:conditionaction:oneBackFP1.8",
                                     alternative='two.sided',
                                     mu=0)

# Three-way interaction FPn 1.8 x condition (action) x FPn-1 1.8
fp18actiononeBackFP18ttest <- t.test(coefmatrix$"foreperiod1.8:conditionaction:oneBackFP1.8",
                                     alternative='two.sided',
                                     mu=0)

# 4.2. Separate tests by foreperiod
summaryDatafp06 <- summaryData %>%
  filter(foreperiod=='0.6')

# 2.2.1. Anova for FP n-1
seqEffAnovafp06 <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryDatafp06,
                      within = c("condition", "oneBackFP"))

nice(seqEffAnovafp06,
     correction='none')

seqEffAnovafp06 <- aov_ez(id='ID',
                          dv='meanSeqEff',
                          data=summaryDatafp06,
                          within=c('condition','oneBackFP'),
                          na.rm=TRUE)

nice(seqEffAnovafp06,
     correction='none')

seqEFffregressionfp06 <- lm(meanRT ~ oneBackFP * condition, 
                        data = summaryData)

Anova(seqEFffregression)
anova(seqEFffregression)

logseqEffregression <- lm(meanRT~logFP*logoneBackFP*condition,
                          data=summaryData) 
anova(logseqEffregression)
#===============================================================#
#============================= Plots ============================
#===============================================================#

# Only by foreperiod
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           group = 1)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line") +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))


# Sequential effects (aggregated across conditions)
ggplot(data = summaryData2,
       aes(x = oneBackFP,
           y = meanRT)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", size = 0.8) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))


# Sequential effects (separated by condition)
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color = oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue','orange','green'))

# Sequential effects (separated by oneBackFP)
ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", size = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  facet_wrap(~oneBackFP) +
  scale_color_manual(values = c('orange','blue'))



# Only action data
ggplot(data = filter(summaryData, condition == "action"),
       aes(x = foreperiod,
           y = meanRT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line" ) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))

# Accuracy
ggplot(data = goData,
       aes(x = foreperiod,
           y = Acc)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))

# Fit by participant
ggplot(data = goData,
       aes(x = foreperiod,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.0))) +
  facet_wrap(~ID) +
  scale_color_manual

# Fit by participant
ggplot(data = goData,
       aes(x = foreperiod,
           y = RT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group=1)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.0)),
        axis.title = element_text(size = rel(1.0))) +
  facet_wrap(~ID) +
  scale_color_manual(values = c("orange","blue"))
  #ylim(0.25,0.55)

# Boxplots
boxplots <- ggplot(data=summaryData,
                   aes(x=foreperiod,
                       y=meanRT,
                       color=condition))+
  geom_boxplot()+
  scale_color_manual(values=c('orange','blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

boxplots


logBoxplots <- ggplot(data=summaryData,
                   aes(x=foreperiod,
                       y=meanLogRT,
                       color=condition))+
  geom_boxplot()+
  scale_color_manual(values=c('orange','blue'))+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

logBoxplots


# Histograms
histograms <- ggplot(data=summaryData,
                     aes(x=meanRT,
                         color=foreperiod))+
  geom_histogram()+
  facet_grid(foreperiod~condition)+
  theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank())

histograms <- ggplot(data=goData,
                     aes(x=RT))+
  geom_histogram()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  facet_wrap(~ID)


ggplot(data = summaryData,
       aes(x = foreperiod,
           y = meanLogRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", size = 1, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", width = 0.2, geom = "errorbar") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5))) +
  scale_color_manual(values = c("orange", "blue"))

