

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
library(ggsignif)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

# Read data
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
  theme(axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8))) +
  labs(x = "External fixation duration",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/extfixduration.tiff",
       width = 13.4,
       height = 10)

# Check for influence of latency of action key press on RT
ggplot(data=filter(goData,condition=='action'),
       aes(x=action_trigger.rt,
           y=RT,
           color=foreperiod))+
  geom_point() +
  theme(axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8))) +
  labs(x = "Action trigger delay",
       y = "RT") +
  facet_wrap(~foreperiod)
ggsave("./Analysis/Plots/actiontrigpress.tiff",
         width = 13.4,
         height = 10)
  
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
contrasts(summaryData$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

contrasts(summaryData2$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryData2$condition) <- c(-1/2, 1/2)
contrasts(summaryData2$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryData2$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

#==================== 2.1. FP x RT by condition ======================
lines_by_condition <- ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.2, geom = "errorbar") + 
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

# Anova with RT
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


# log FP
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanLogRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.2, geom = "errorbar") + 
  labs(title = "LogRT by condition",
       x = "Foreperiod",
       y = "Mean LogRT") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"))

ggplot(data = summaryData2,
       aes(x = condition,
           y = meanLogRT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group=1)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.2, geom = "errorbar") + 
  labs(title = "LogRT by condition",
       x = "condition",
       y = "Mean LogRT") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)))

logFPAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"))

logFPAnova

#======================= 2.2. Sequential effects ============================
#========== 2.2.1. RT ============
# Sequential effects (separated by condition)
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.1, geom = "errorbar") +
  labs(title = "Sequential effects by condition",
       x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-1") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange', 'green'))
ggsave("./Analysis/Plots/SeqEff.png",
       width = 6.7,
       height = 5)

# Only FP n-1 (no FP)
ggplot(data = summaryData2,
       aes(x = oneBackFP,
           y = meanRT,
           color=condition)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.1, geom = "errorbar") +
  labs(title = "FP n-1 x condition",
       x = "FP n-1",
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"))
ggsave("./Analysis/Plots/onebackxCondition.png",
       width = 6.7,
       height = 5)


# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                  dv = "meanRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition", "oneBackFP"))

seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanInvRT",
                      data = summaryData2,
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

# Sequential effects of FP n-2
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=twoBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = twoBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.1, geom = "errorbar") +
  labs(title = "Sequential effects of FP n-2 by condition",
       x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-2") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  facet_wrap(~condition) +
  scale_color_manual(values = c('blue', 'orange', 'green'))
ggsave("./Analysis/Plots/SeqEffTwoBack.png",
       width = 6.7,
       height = 5)

# Anova for FP n-2
seqEfffp2Anova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "twoBackFP"))

seqEffffp2Anova <- aov_ez(id = "ID",
                      dv = "meanInvRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBackFP"))


nice(seqEfffp2Anova,
     correction='none')

#=================== 2.2.2. Accuracy ==========================
# Set constrasts for variables used in ANOVAs
contrasts(summaryDataAll$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryDataAll$condition) <- c(-1/2, 1/2)
contrasts(summaryDataAll$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryDataAll$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)

ggplot(data = summaryDataAll,
       aes(x = foreperiod,
           y = meanAcc,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.05, linewidth = 0.7) +
  scale_color_manual(values = c("orange", "blue")) +
  labs(x = "Foreperiod",
       y = "Mean Acc", 
       color = "Condition") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.8)),
        axis.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8)))
ggsave("./Analysis/Plots/acc_by_condition.png",
       width = 15,
       height = 10)

ggplot(data = summaryDataAll,
       aes(x = foreperiod,
           y = meanAcc,
           color = oneBackFP)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = oneBackFP)) +
  facet_wrap(~ condition)


ggplot(data = summaryDataAll,
       aes(x = foreperiod,
           y = meanAcc,
           color = oneBackFPGo)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = oneBackFPGo)) +
  facet_wrap(~ condition)

ggplot(data = summaryDataAll,
       aes(x = foreperiod,
           y = meanAcc,
           color = oneBackFP)) +
  geom_jitter()

ggplot(data = summaryDataAll) +
  geom_histogram(aes(x = meanAcc)) +
  facet_grid(foreperiod ~ condition)

AccAnova <- aov_ez(id = "ID",
                   dv = "meanAcc",
                   data = summaryDataAll,
                   within = c("foreperiod", "condition", "oneBackFPGo"))


#================= 2.2.3. N-1 trial type, condition and sequential effects =================
ggplot(data = summaryData2,
       aes(x = oneBacktrialType,
           y = meanRT,
           color=condition)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.1, geom = "errorbar") +
  labs(x = "FP n-1 trial type",
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  facet_wrap(~foreperiod) +
  scale_color_manual(values = c('orange', 'blue'))
ggsave("./Analysis/Plots/oneBackTrialType.png",
       width = 7.5,
       height = 5)

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", size = 0.8, width = 0.1, geom = "errorbar") +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-1") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.8)),
        strip.text = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8))) +
  facet_grid(oneBacktrialType ~ condition) +
  scale_color_manual(values = c("blue", "orange", "green"))
ggsave("./Analysis/Plots/RT_condition_trialtype.png",
       width = 13.4,
       height = 10)

ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(oneBacktrialType ~ condition) +
  scale_color_manual(values = c("blue", "orange", "green")) +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-1") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)))

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
                    data = summaryData2,
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


#================= 2.4. Condition and sequential effects including no-go FP n-1 =================
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = oneBackFPGo)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = oneBackFPGo)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
  facet_wrap(~ condition) +
  scale_color_manual(values = c("blue", "orange", "green", "magenta")) +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-1") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.text = element_text(size = rel(1.2)),
        legend.title = element_text(size = rel(1.2)))
ggsave("./Analysis/Plots/seqEffNoGO.png")


# Anova
seqEffGoAnova <- aov_ez(id = "ID",
                        dv = "meanRT",
                        data = summaryData2,
                        within = c("foreperiod", "oneBackFPGo", "condition"))

nice(seqEffGoAnova,
     correction='none')

seqEffGoemm <- emmeans(seqEffGoAnova, pairwise ~ oneBackFPGo * condition| foreperiod)
contrast(seqEffGoemm[[1]], interaction = c("consec", "pairwise"), by = "foreperiod")


#==================== 2.4. Foreperiod, condition and block ======================
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)



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


