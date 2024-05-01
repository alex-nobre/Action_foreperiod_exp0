

# Load necessary packages

# Read and process data
library(tidyverse)
library(magrittr)
library(data.table)

# Plotting
library(ggplot2)
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggpubr)

# Descriptives
library(Hmisc)

# Linear models
library(afex)
library(car)
library(codingMatrices)
library(broom)
library(modelr)

# Mixed models
library(emmeans)
library(lme4)
library(performance)
library(ggsignif)
library(rtdists)

# Bayesian analyses
library(bayestestR)
library(BayesFactor)
library(brms)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

# Read data
source('./Analysis/Prepare_data_4.R')
source("./Analysis/R_rainclouds.R")

#==========================================================================================#
#======================================= 0. Data quality ===================================
#==========================================================================================#
testnormality = function(dfr) return(shapiro.test(dfr$invRT)$p.value)
p = as.vector(by(goData, goData$ID, testnormality))
names(p) = levels(goData$ID)
names(p[p < 0.05])


qqmath(~RT|ID, data=goData)
qqmath(~invRT|ID, data=goData)


# % of extreme RTs removed
(ntrials_before_extrem - ntrials_after_extrem)/ntrials_before_extrem * 100

# % of RTs removed due to trimming by participant
n_notrim <- goData %>%
  group_by(ID) %>%
  summarise(n_trials_no_trim = n())

n_trim <- goData2 %>%
  group_by(ID) %>%
  summarise(n_trials_trim = n())

n_trials_comp <- left_join(n_notrim, n_trim) %>%
  mutate(trimmed = n_trials_no_trim - n_trials_trim,
         percent_trimmed = (trimmed/n_trials_no_trim) * 100) %>%
  summarise(meanT = mean(percent_trimmed), sdT = sd(percent_trimmed))

#==================== Practice effects ====================
# RTs across blocks (no condition)
ggplot(data=summaryData2,
       aes(x=block,
           y=meanRT))+
  stat_summary(fun='mean',geom='point')+
  stat_summary(fun='mean',geom='line',size=1,aes(group=1))+
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

# RTs across blocks (by condition)
ggplot(data=summaryData2,
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

#================= Counterbalancing order ==================
# RTs across blocks by counterbalancing order
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
  

#================ Check for influence of external fixation duration ===================
cor(goData2 %>% filter(condition == "external") %>% pull(RT),
    goData2 %>% filter(condition == "external") %>% pull(extFixationDuration))

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

ggsave("./Analysis/Plots/extfixduration.png",
       width = 13.4,
       height = 10)

#============== Check for influence of latency of action key press on RT ================
meanAT <- mean(goData2 %>% filter(condition == "action") %>% pull(action_trigger.rt))
sdAT <- sd(goData2 %>% filter(condition == "action") %>% pull(action_trigger.rt))

cor(goData2 %>% filter(condition == "action", action_trigger.rt < meanAT + 3 * sdAT) %>% pull(RT),
    goData2 %>% filter(condition == "action", action_trigger.rt < meanAT + 3 * sdAT) %>% pull(action_trigger.rt))


ggplot(data=filter(goData,condition=='action',
                   action_trigger.rt < meanAT + 3 * sdAT),
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

ggsave("./Analysis/Plots/actiontrigpress.png",
       width = 13.4,
       height = 10)

# Correlations between RT and accuracy
ggplot(data=summaryDataAccGo,
       aes(x=meanRT,
           y=meanAcc)) +
  geom_point()

#================ Influence of delays between actions and WS onset ============
# Histograms of delays between action and WS
ggplot(data = goData2) +
  geom_histogram(aes(x = round(delay,3)))

# Frequency of values in delay vector
table(round(goData$delay, 3))

# Plot RT against delay excluding zero-delay trials
ggplot(data = filter(goData2, delay > 0, delay < 0.034), aes(x = delay, y = RT)) +
  geom_point() +
  geom_abline()

# Plot RT against delay with best-fitting regression line
ggplot(data = goData2, aes(x = delay, y = RT)) +
  geom_point() +
  geom_abline() +
  ylim(c(0,1.0))

# Correlation between RT and delay excluding zero-delay trials
cor(filter(goData2, delay > 0)$delay, filter(goData2, delay > 0)$RT)

#================================= 1.2. Stopping-rule ======================================

#================ 1.2.1. Prepare data ================

# Variables used as predictors: numForeperiod and numOneBackFP
# Dependent variable: logRT
# Variables nested by condition and ID

buildmodel <- function(data) {
  lm(logRT ~ numForeperiod*numOneBackFP,
     data = data)
}

nested_data <- goData2 %>%
  select(ID, condition, numForeperiod, numOneBackFP, logRT) %>%
  group_by(ID, condition) %>%
  nest()

fitted_data <- nested_data %>%
  mutate(fit = map(data, buildmodel),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, condition, term, estimate) %>%
  pivot_wider(names_from = term,
              values_from = estimate)


#============================== 1.2.2. Sequential bayes factors ===========================
external_fits <- fitted_data[fitted_data$condition=='external',]
action_fits <- fitted_data[fitted_data$condition=='action',]

srange <- 10:nrow(external_fits)

fp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$numForeperiod[1:range],
                    y = action_fits$numForeperiod[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, fp_bfs)
lines(srange, fp_bfs)

onebackfp_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$numOneBackFP[1:range],
                    y = action_fits$numOneBackFP[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, onebackfp_bfs)
lines(srange, onebackfp_bfs)

interact_bfs <- sapply(srange, function(range) {
  extractBF(ttestBF(x = external_fits$`numForeperiod:numOneBackFP`[1:range],
                    y = action_fits$`numForeperiod:numOneBackFP`[1:range],
                    paired=TRUE),
            onlybf = TRUE)
})

plot(srange, interact_bfs)
lines(srange, interact_bfs)

#============================ 1.2.3. Mixed models BF comparison ============================ 


# Condition x fp n
with_condition <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          (1 + condition + numForeperiod | ID),
                        data = goData2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = logRT ~ 1 + numForeperiod +
                        (1 + numForeperiod | ID),
                      data = goData2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_condition)

bic_to_bf(c(BIC(no_condition),
            BIC(with_condition)),
          denominator = c(BIC(no_condition)))

# Condition x fp n x fp n-1

with_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          oneBackFP + numForeperiod:oneBackFP + condition:oneBackFP + 
                          (1 + condition + numForeperiod | ID),
                        data = data2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_onebackfp)

no_onebackfp <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod +
                        (1 + condition + numForeperiod | ID),
                      data = data2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')


isSingular(no_onebackfp)

BIC(no_onebackfp, with_onebackfp)

bic_to_bf(c(BIC(no_onebackfp),
            BIC(with_onebackfp)),
          denominator = c(BIC(no_onebackfp)))


# Condition x fp n x fp n-1 (with and without condition)
with_condition <- mixed(formula = logRT ~ 1 + condition + numForeperiod + condition:numForeperiod + 
                          oneBackFP + numForeperiod:oneBackFP + condition:oneBackFP + condition:numForeperiod:oneBackFP +
                          (1 + condition + numForeperiod | ID),
                        data = goData2,
                        control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                        progress = TRUE,
                        expand_re = TRUE,
                        method = 'S',
                        REML = TRUE,
                        return = 'merMod')

isSingular(with_condition)

no_condition <- mixed(formula = logRT ~ 1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP +
                        (1 + numForeperiod | ID),
                      data = goData2,
                      control = lmerControl(optimizer = c('bobyqa'), optCtrl = list(maxfun=2e5), calc.derivs = FALSE),
                      progress = TRUE,
                      expand_re = TRUE,
                      method = 'S',
                      REML = TRUE,
                      return = 'merMod')

isSingular(no_condition)

bic_to_bf(c(BIC(no_condition),
            BIC(with_condition)),
          denominator = c(BIC(no_condition)))

#==========================================================================================#
#======================================= 1. Descriptives ==================================
#==========================================================================================#
#========== 1.1. Plots ==============
  
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

# Only by condition
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


# Histograms of slopes
buildmodel <- function(data) {
  lm(logRT ~ numForeperiod*numOneBackFP,
     data = data)
}

nested_data <- goData2 %>%
  select(ID, condition, numForeperiod, numOneBackFP, logRT) %>%
  group_by(ID, condition) %>%
  nest()

fitted_data <- nested_data %>%
  mutate(fit = map(data, buildmodel),
         params = map(fit, tidy)) %>%
  ungroup() %>%
  unnest(c(params)) %>%
  select(ID, condition, term, estimate) %>%
  pivot_wider(names_from = term,
              values_from = estimate)


external_fits <- fitted_data[fitted_data$condition=='external',]
action_fits <- fitted_data[fitted_data$condition=='action',]

ggplot(data = fitted_data, aes(x = numForeperiod)) +
  geom_histogram() +
  facet_wrap(~condition)

# Plot grand means
RT_by_condition <- ggplot(data = summaryData2 %>% 
                               group_by(ID, foreperiod, condition) %>% 
                               summarise(meanRT = mean(meanRT)),
                             aes(x = foreperiod,
                                 y = meanRT,
                                 color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 1.0, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.0, width = 0.1, geom = "errorbar") + 
  labs(title = "RT",
       x = "FP (s)",
       y = "Mean RT (s)",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.key = element_blank()) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.jpg",
                RT_by_condition,
                width = 15,
                height = 10)

# Plot means by participant
rt_by_fp_cond_part <- ggplot(data = summaryData2,
                             aes(x = foreperiod,
                                 y = meanRT,
                                 color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.2, geom = "errorbar") + 
  labs(title = "RT by condition",
       x = "Foreperiod",
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action")) +
  facet_wrap(~ ID)
ggplot2::ggsave("./Analysis/Plots/RT_by_cond_fp_part.png",
                rt_by_fp_cond_part,
                width = 6.7,
                height = 5)




# Boxplot of RT by FP and condition
boxplot_fp_condition <- ggplot(data = summaryData2,
                               aes(x = foreperiod,
                                   y = meanRT,
                                   fill = condition)) +
  geom_boxplot() + 
  labs(x = "Foreperiod",
       y = "Mean RT",
       fill = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_fill_manual(values = c("orange", "blue"), labels = c("External", "Action"))
ggplot2::ggsave("./Analysis/Plots/boxplot_RT_by_condition.png",
                boxplot_fp_condition,
                width = 6.7,
                height = 5)

# log RT by FP and condition
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

# Using exp(-foreperiod)
rt_by_expfp_condition <- ggplot(data = summaryData2 %>% 
                                  group_by(ID, expFP, condition) %>% 
                                  summarise(meanRT = mean(meanRT)),
                                aes(x = expFP,
                                    y = meanRT,
                                    color = condition)) +
  geom_jitter(height = 0, width = 0.005, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group = condition)) +
  stat_summary(fun.data = "mean_se", linewidth = 1.2, width = 0.02, geom = "errorbar") +
  labs(title = "RT",
       x = "exp(-FP)", 
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange","blue"))


# Sequential effects (separated by condition)
ggplot(data = summaryData2 %>%
         group_by(ID, foreperiod, condition, oneBackFP) %>%
         summarise(meanRT = mean(meanRT)),
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  geom_jitter(height = 0, width = 0.30, alpha = 0.5) +
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

ggsave("./Analysis/Plots/SeqEff_by_condition.tiff",
       width = 20,
       height = 11.11,
       unit = "cm")

# Sequential effects separated by FP n-1
seqEff_by_oneback <- ggplot(data = summaryData2 %>%
                              group_by(ID, foreperiod, condition, oneBackFP) %>%
                              summarise(meanRT = mean(meanRT)),
                            aes(x = foreperiod,
                                y = meanRT,
                                color=condition)) +
  geom_jitter(height = 0, width = 0.30, size = 1.0, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.0, width = 0.1, geom = "errorbar") + 
  labs(x = "FP (s)",
       y = "Mean RT (s)",
       color = "Condition") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = rel(0.9)),
        axis.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(0.9)),
        legend.text = element_text(size = rel(0.8)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt")) +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`0.6` = "FP[n-1] == 0.6",
                                      `1.2` = "FP[n-1] == 1.2",
                                      `1.8` = "FP[n-1] == 1.8"),
                                    default = label_parsed)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

ggsave("./Analysis/Plots/SeqEff.jpg",
       seqEff_by_oneback,
       width = 20,
       height = 8.33,
       units = "cm")


# Boxplot
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           fill = oneBackFP)) +
  geom_boxplot() +
  labs(x = "Foreperiod",
       y = "Mean RT",
       fill = "FP n-1") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        strip.background = element_blank(),
        strip.text = element_text(size = 14)) +
  facet_wrap(~condition,
             labeller = as_labeller(c(external = "External",
                                      action = "Action"))) +
  scale_fill_manual(values = c('blue', 'orange', 'green'))
ggsave("./Analysis/Plots/SeqEff_boxplot.pdf",
       width = 6.7,
       height = 5)

# Raincloud
ggplot(data = summaryData2 %>%
         group_by(ID, foreperiod, condition, oneBackFP) %>%
         summarise(meanRT = mean(meanRT)),
       aes(x = foreperiod,
           y = meanRT,
           fill=foreperiod)) +
  stat_halfeye(adjust = 0.8, .width = 0, point_color = NA, position = "dodge") +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = 1)) +
  #stat_dots(side = "left") +
  facet_grid(condition ~ oneBackFP) +
  labs(x = "Foreperiod",
       y = "Mean RT",
       color = "FP n-1") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_fill_manual(values = c('blue', 'orange', 'green'))

ggsave("./Analysis/Plots/SeqEff_raincloud.jpeg",
#ggsave("./Analysis/Plots/SeqEff_raincloud.pdf",
       width = 6.7,
       height = 5)


# Only FP n-1 (no FP n)
ggplot(data = summaryData2,
       aes(x = oneBackFP,
           y = meanRT,
           color=condition)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.1, geom = "errorbar") +
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

# Accuracy
error_plot <- ggplot(data = summaryDataAcc %>% 
                     group_by(ID, foreperiod, condition) %>%
                     summarise(errorRate = mean(errorRate)),
                   aes(x = foreperiod,
                       y = errorRate,
                       color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 1.0, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.0, width = 0.1, geom = "errorbar") +
  scale_color_manual(values = c("orange", "blue"),
                     labels = c("External", "Action")) +
  labs(title = "All trials",
       x = "FP (s)",
       y = "Mean Error Rate", 
       color = "Condition") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.8)),
        axis.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8)),
        plot.title = element_text(hjust = 0.5, size = rel(2.0)),
        legend.key = element_blank())
ggsave("./Analysis/Plots/acc_by_condition.jpg",
       acc_plot,
       width = 15,
       height = 10)

# Accuracy go
error_go_plot <- ggplot(data = summaryDataAccGo %>% 
         group_by(ID, foreperiod, condition) %>%
         summarise(errorRate = mean(errorRate)),
       aes(x = foreperiod,
           y = errorRate,
           color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 1.0, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.0, width = 0.1, geom = "errorbar") + 
  scale_color_manual(values = c("orange", "blue"),
                     labels = c("External", "Action")) +
  labs(title = "Misses",
       x = "FP (s)",
       y = "Mean Error Rate", 
       color = "Condition") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.8)),
        axis.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8)),
        plot.title = element_text(hjust = 0.5, size = rel(2.0)),
        legend.key = element_blank())
ggsave("./Analysis/Plots/error_by_condition_go.jpg",
       error_go_plot,
       width = 15,
       height = 10)

# Accuracy no go
error_nogo_plot <- ggplot(data = summaryDataAccNoGo %>% 
         group_by(ID, foreperiod, condition) %>%
         summarise(errorRate = mean(errorRate)),
       aes(x = foreperiod,
           y = errorRate,
           color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 1.0, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.0, width = 0.1, geom = "errorbar") + 
  scale_color_manual(values = c("orange", "blue"),
                     labels = c("External", "Action")) +
  labs(title = "False alarms",
       x = "FP (s)",
       y = "Mean Error Rate", 
       color = "Condition") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.8)),
        axis.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.8)),
        plot.title = element_text(hjust = 0.5, size = rel(2.0)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        plot.margin = unit(c(5.5, 5.5, 5.5, 1), "pt"))
ggsave("./Analysis/Plots/error_by_condition_nogo.jpg",
       error_nogo_plot,
       width = 15,
       height = 10)

# Error plots in single panel
# Save legend as grob
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(error_nogo_plot + 
                                                          theme(legend.title = element_text(size = rel(1.1)),
                                                                legend.text = element_text(size = rel(0.9))))), "guide-box")

# Visualize
grid.arrange(error_plot + theme(legend.position = "none",
                              axis.text = element_text(size = rel(1.2)),
                              axis.title = element_text(size = rel(1.4)),
                              plot.title = element_text(size = rel(1.5))),
             error_go_plot + theme(legend.position = "none",
                                 axis.text = element_text(size = rel(1.2)),
                                 axis.title = element_text(size = rel(1.4)),
                                 axis.title.y = element_blank(),
                                 plot.title = element_text(size = rel(1.5))),
             error_nogo_plot + theme(legend.position = "none",
                                   axis.text = element_text(size = rel(1.2)),
                                   axis.title = element_text(size = rel(1.4)),
                                   axis.title.y = element_blank(),
                                   plot.title = element_text(size = rel(1.5))),
             cond_legend,
             nrow = 1,
             widths = c(3/10, 3/10, 3/10, 1/10))

# Save
error_all_plots <- arrangeGrob(error_plot + theme(legend.position = "none",
                                              axis.text = element_text(size = rel(1.2)),
                                              axis.title = element_text(size = rel(1.4)),
                                              plot.title = element_text(size = rel(1.5))),
                             error_go_plot + theme(legend.position = "none",
                                                 axis.text = element_text(size = rel(1.2)),
                                                 axis.title = element_text(size = rel(1.4)),
                                                 axis.title.y = element_blank(),
                                                 plot.title = element_text(size = rel(1.5))),
                             error_nogo_plot + theme(legend.position = "none",
                                                   axis.text = element_text(size = rel(1.2)),
                                                   axis.title = element_text(size = rel(1.4)),
                                                   axis.title.y = element_blank(),
                                                   plot.title = element_text(size = rel(1.5))),
                             cond_legend,
                             nrow = 1,
                             widths = c(3/10, 3/10, 3/10, 1/10))

ggsave("./Analysis/Plots/error_all_plots.jpg",
       error_all_plots,
       width = 12,
       height = 6.67)

# RT and accuracy in single panel
# Save legend as grob
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(error_nogo_plot + 
                                                          theme(legend.title = element_text(size = rel(1.1)),
                                                                legend.text = element_text(size = rel(0.9))))), "guide-box")

# Visualize
grid.arrange(RT_by_condition + labs(tag = "A", x = "") +
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(1.1)),
                     axis.title = element_text(size = rel(1.3)),
                     plot.title = element_text(size = rel(1.4)),
                     plot.tag = element_text(size = rel(1.5))),
             error_go_plot + labs(tag = "B") +
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(1.1)),
                     axis.title = element_text(size = rel(1.3)),
                     plot.title = element_text(size = rel(1.4)),
                     plot.tag = element_text(size = rel(1.5))),
             error_nogo_plot + labs(tag = "C", x = "") + 
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(1.1)),
                     axis.title = element_text(size = rel(1.3)),
                     plot.title = element_text(size = rel(1.4)),
                     plot.tag = element_text(size = rel(1.5))) +
               labs(y = ""),
             cond_legend,
             nrow = 1,
             widths = c(5/17, 5/17, 5/17, 2/17))


rt_error_plots <- arrangeGrob(RT_by_condition + labs(tag = "A", x = "") +
                              theme(legend.position = "none",
                                    axis.text = element_text(size = rel(1.1)),
                                    axis.title = element_text(size = rel(1.3)),
                                    plot.title = element_text(size = rel(1.4)),
                                    plot.tag = element_text(size = rel(1.5))),
                            error_go_plot + labs(tag = "B") +
                              theme(legend.position = "none",
                                    axis.text = element_text(size = rel(1.1)),
                                    axis.title = element_text(size = rel(1.3)),
                                    plot.title = element_text(size = rel(1.4)),
                                    plot.tag = element_text(size = rel(1.5))),
                            error_nogo_plot + labs(tag = "C", x = "") + 
                              theme(legend.position = "none",
                                    axis.text = element_text(size = rel(1.1)),
                                    axis.title = element_text(size = rel(1.3)),
                                    plot.title = element_text(size = rel(1.4)),
                                    plot.tag = element_text(size = rel(1.5))) +
                              labs(y = ""),
                            cond_legend,
                            nrow = 1,
                            widths = c(5/17, 5/17, 5/17, 2/17))

ggsave("./Analysis/Plots/rt_error_plots.jpg",
       rt_error_plots,
       width = 20,#16,
       height = 8.33,#11.94,
       units = "cm")


#  n-1 trial type
ggplot(data = summaryData2,
       aes(x = oneBacktrialType,
           y = meanRT,
           color=condition)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.1, geom = "errorbar") +
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


# Sequential effects separated by FP n-1 and n-1 trial type
seqEff_by_oneback_trialtype <- ggplot(data = summaryData2 %>%
         group_by(ID, foreperiod, oneBackFP, oneBacktrialType, condition) %>%
         summarise(meanRT = mean(meanRT)),
       aes(x = foreperiod,
           y = meanRT,
           color=condition)) +
  geom_jitter(height = 0, width = 0.30, size = 1.0, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.2, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.0, width = 0.1, geom = "errorbar") + 
  labs(x = "FP (s)",
       y = "Mean RT (s)",
       color = "Condition") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.text = element_text(size = rel(0.9)),
        axis.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        legend.title = element_text(size = rel(0.9)),
        legend.text = element_text(size = rel(0.8)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt")) +
  facet_grid(oneBacktrialType~oneBackFP,
             labeller = labeller(oneBackFP = as_labeller(c(`0.6` = "FP[n-1] == 0.6",
                                               `1.2` = "FP[n-1] == 1.2",
                                               `1.8` = "FP[n-1] == 1.8"),
                                             default = label_parsed),
                                 oneBacktrialType = c("go" = "Go", "no-go" = "No-go"))) +
  scale_color_manual(values = c('orange', 'blue'), labels = c("External", "Action"))

ggsave("./Analysis/Plots/SeqEff_oneBackTrialType.jpg",
       seqEff_by_oneback_trialtype,
       width = 20,
       height = 14.33,
       units = "cm")

ggsave("./Analysis/Plots/SeqEff_oneBackTrialType.tiff",
       seqEff_by_oneback_trialtype,
       width = 20,
       height = 17.11,
       unit = "cm")

# Sequential effects separated by condition and n-1 trial type
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color=oneBackFP)) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = oneBackFP)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 0.8, width = 0.1, geom = "errorbar") +
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
ggsave("./Analysis/Plots/RT_condition_trialtype.pdf",
       width = 6.7,
       height = 5)


# SD
sd_by_condition <- ggplot(data = goData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(sdRT = sd(RT)),
                          aes(x = foreperiod,
                              y = sdRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.2, width = 0.1, geom = "errorbar") + 
  labs(title = "RT",
       x = "FP (s)",
       y = "SD RT (s)",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2)),
        legend.key = element_blank()) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.pdf",
                sd_by_condition,
                width = 15,
                height = 10)

#==================================== 1.2. Tables ================================
# RT By foreperiod and condition
meandata <- summaryData2 %>%
  group_by(condition, foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandata  

 # RT By foreperiod
meandataFP <- summaryData2 %>%
  group_by(foreperiod) %>%
  summarise(condRT = mean(meanRT),
            varRT = var(meanRT),
            sdRT = sd(meanRT))
meandataFP

# marginal means for acc
margAcc <- summaryDataAcc %>%
  group_by(foreperiod, condition) %>%
  summarise(condAcc = mean(meanAcc),
            varAcc = var(meanAcc),
            sdAcc = sd(meanAcc))
margAcc  

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
# Anova with RT
fpAnova <- aov_ez(id = "ID",
       dv = "meanRT",
       data = summaryData2,
       within = c("foreperiod", "condition"),
       anova_table = list(es = "pes"))


# Check anova requirements
check_sphericity(fpAnova)

is_norm <- check_normality(fpAnova)

plot(is_norm)

plot(is_norm, type = "qq")

summary(fpAnova)

# Pairwise comparisons
fpEmmeans <- emmeans(fpAnova,
                     pairwise ~ foreperiod|condition,
                     adjust = 'bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('poly'),
                               adjust='bonferroni')

fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
                               interaction=c('consec'),
                               adjust='bonferroni')



# Anova with logRT
logFPAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"))

logFPAnova

check_sphericity(logFPAnova)


# Anova with expFP
expFPAnova <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2,
                     within = c("expFP", "condition"))

expFPAnova


# Check anova requirements
check_sphericity(expFPAnova)

check_normality(expFPAnova) |>
  plot()

check_normality(expFPAnova) |>
  plot(type = "qq")



#======================= 2.2. Sequential effects ============================

#========== 2.2.1. RT ============
# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                  dv = "meanRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition", "oneBackFP"),
                  anova_table = list(es = "pes"))


check_sphericity(seqEffAnova)

is_norm <- check_normality(seqEffAnova)

nice(seqEffAnova)

# Pairwise comparisons for conditions wihtin levels of foreperiod
conditionEmmeans <- emmeans(seqEffAnova,
                     pairwise ~ condition|foreperiod,
                     adjust = "holm")


conditionEmmeansContrasts <- contrast(conditionEmmeans[[1]],
                               interaction = c("consec"),
                               adjust = "holm")
                  

# Pairwise comparisons for levels of oneBackFP within foreperiod
oneBackEmmeans <- emmeans(seqEffAnova,
                          pairwise ~ oneBackFP|foreperiod,
                          adjust = "holm")

oneBackfpEmmeansConstrasts <- contrast(oneBackEmmeans[[1]],
                                        interaction = c("consec"),
                                        adjust = "holm")

#========== 2.2.2. 1/RT ============
seqEffInvAnova <- aov_ez(id = "ID",
                         dv = "meanInvRT",
                         data = summaryData2,
                         within = c("foreperiod", "condition", "oneBackFP"))


nice(seqEffInvAnova,
     correction='none')


#========== 2.2.3. logRT ============
seqEffLogAnova <- aov_ez(id = "ID",
                      dv = "meanLogRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBackFP"),
                      anova_table = list(es = "pes"))

check_sphericity(seqEffLogAnova)

is_norm <- check_normality(seqEffLogAnova)

nice(seqEffLogAnova)

# Pairwise comparisons for conditions wihtin levels of foreperiod
conditionEmmeans <- emmeans(seqEffLogAnova,
                            pairwise ~ condition|foreperiod,
                            adjust = "holm")


conditionEmmeansContrasts <- contrast(conditionEmmeans[[1]],
                                      interaction = c("consec"),
                                      adjust = "holm")


# Pairwise comparisons for levels of oneBackFP within foreperiod
oneBackEmmeans <- emmeans(seqEffLogAnova,
                          pairwise ~ oneBackFP|foreperiod,
                          adjust = "holm")

oneBackfpEmmeansConstrasts <- contrast(oneBackEmmeans[[1]],
                                       interaction = c("consec"),
                                       adjust = "holm")



#========== 2.2.4. expFP ============
seqEffExpAnova <- aov_ez(id = "ID",
                         dv = "meanRT",
                         data = summaryData2,
                         within = c("expFP", "condition", "oneBackFP"),
                         anova_table = list(es = "pes"))

check_sphericity(seqEffExpAnova)

check_normality(seqEffExpAnova) |>
  plot()

summary(seqEffExpAnova)


# Without no-go n-1 trials

# RT
seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData3,
                      within = c("foreperiod", "condition", "oneBackFP"),
                      anova_table = list(es = "pes"))
# Anova for FP n-2
seqEff2Anova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "twoBackFP"))



#=================== 2.2.3. Accuracy ==========================
# Set constrasts for variables used in ANOVAs
contrasts(summaryDataAll$foreperiod) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryDataAll$condition) <- c(-1/2, 1/2)
contrasts(summaryDataAll$oneBackFP) <- contr.treatment(3)-matrix(rep(1/3,6),ncol=2)
contrasts(summaryDataAll$oneBackFPGo) <- contr.treatment(4)-matrix(rep(1/4,12),ncol=3)


ggplot(data = summaryDataAcc,
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


#================= 2.2.4. N-1 trial type, condition and sequential effects =================

# 2.2.4.1. Anova with foreperiod, condition and n-1 trial type
trialTypeAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData2,
                      within = c("foreperiod", "condition", "oneBacktrialType"),
                      anova_table = list(es = "pes"))

nice(trialTypeAnova)

ttypEmmeans <- emmeans(trialTypeAnova,
                       pairwise ~ oneBacktrialType|condition*foreperiod,
                       adjust = "holm")

ttypEmmeansContrast <- contrast(ttypEmmeans[[1]],
                                interaction = c("pairwise"),
                                adjust = "holm")

# 2.2.4.2. Anova with foreperiod, n-1 trial type and sequential effects
seqEfftrialTypeAnova <- aov_ez(id = "ID",
                         dv = "meanRT",
                         data = summaryData,
                         within = c("foreperiod", "oneBackFP", "oneBacktrialType"))

nice(seqEfftrialTypeAnova,
     correction='none')

# 2.2.4.3. Anova with foreperiod, condition, n-1 trial type and sequential effects
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


#================= 2.5. Condition and sequential effects including no-go FP n-1 =================

# Anova
seqEffGoAnova <- aov_ez(id = "ID",
                        dv = "meanRT",
                        data = summaryData2,
                        within = c("foreperiod", "oneBackFPGo", "condition"),
                        anova_table = list(es = "pes"))

nice(seqEffGoAnova)

seqEffGoemm <- emmeans(seqEffGoAnova, 
                       pairwise ~ oneBackFPGo| condition * foreperiod,
                       adjust = "holm")

contrast(seqEffGoemm[[1]],
         interaction = c("pairwise"),
         adjust = "holm")


#==================== 2.6. Foreperiod, condition and block ======================
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)


# Separate tests by foreperiod
summaryDatafp06 <- summaryData %>%
  filter(foreperiod=='0.6')

#  Anova for FP n-1
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

#==========================================================================================#
#===================== 3. Plots for between-experiments comparisons ========================
#==========================================================================================#
extrFPData <- summaryData2 %>%
  filter(foreperiod %in% c("0.6", "1.8")) %>%
  mutate(foreperiod = fct_relevel(foreperiod, c("0.6", "1.8")))

rt_extr_fps_part <- ggplot(data = extrFPData,
                              aes(x = foreperiod,
                                  y = meanRT,
                                  color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.8, aes(group = condition)) +
  stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", linewidth = 0.8, width = 0.2) +
  labs(title = "RT by condition",
       x = "Foreperiod",
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text = element_text(size = rel(1.2)),
        axis.title = element_text(size = rel(1.2))) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action")) +
  facet_wrap(~ ID)

ggsave("./Analysis/Plots/rt_extr_fps_part.png",
       rt_extr_fps_part,
       width = 6.7,
       height = 5)

rt_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(meanCorRT = mean(meanCorRT)),
                          aes(x = foreperiod,
                              y = meanCorRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group = condition)) +
  stat_summary(fun.data = "mean_se", linewidth = 1.2, width = 0.1, geom = "errorbar") +
  labs(title = "RT",
       x = "FP", 
       y = "Mean RT",
       color = "Condition") +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  scale_color_manual(values = c("orange","blue"))
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.png",
                rt_by_condition,
                width = 8.5,
                height = 5.7)


ggplot(data = summaryData2 %>%
         group_by(ID, foreperiod, condition, oneBackFP) %>%
         summarise(meanCorRT = mean(meanCorRT)),
       aes(x = foreperiod,
           y = meanCorRT,
           color=oneBackFP)) +
  geom_jitter(height = 0, width = 0.30, alpha = 0.5) +
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

#=================================== 4. Learning =============================================

# RTs by block (no condition)
ggplot(data=summaryData2,
       aes(x=foreperiod,
           y=meanRT))+
  stat_summary(fun='mean',geom='point') +
  stat_summary(fun='mean',geom='line',size=1,aes(group=1)) +
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar') +
  facet_wrap(~block,
             labeller = labeller(block = c(`0` = "Block 1",
                                              `1` = "Block 2",
                                              `2` = "Block 3",
                                              `3` = "Block 4"))) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

# RTs by block (only external)
external_by_block <- ggplot(data=summaryData2 %>%
                              filter(condition == "external"),
                            aes(x=foreperiod,
                                y=meanRT,
                                color = "orange"))+
  stat_summary(fun='mean',geom='point') +
  stat_summary(fun='mean',geom='line',size=1,aes(group=1)) +
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar') +
  facet_wrap(~block,
             labeller = labeller(block = c(`0` = "Block 1",
                                           `1` = "Block 2",
                                           `2` = "Block 3",
                                           `3` = "Block 4"))) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank())

ggplot2::ggsave("./Analysis/Plots/external_by_block.pdf",
                external_by_block,
                width = 15,
                height = 10)

#============================ 4.1. Effects across blocks ==============================
# Foreperiod, condition and block
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData2)

anova(blocklm)


ggplot(data = goData2,
       aes(x = trial,
           y = RT,
           color = condition)) +
  #geom_point() +
  geom_line(aes(group = condition)) +
  geom_smooth(method = 'lm') +
  facet_wrap(~ID) +
  scale_color_manual(values = c('orange', 'blue'))

#========================== 5.2. Split anovas by counterbalancing order ===================

fpAnova_ae <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='action-external',],
                     within = c("foreperiod", "condition"))

fpAnova_ea <- aov_ez(id = "ID",
                     dv = "meanRT",
                     data = summaryData2[summaryData2$counterbalance=='external-action',],
                     within = c("foreperiod", "condition"))

# Examine first trials for possible stabilization of RT
firstBlockData <- goData2 %>%
  filter(block == '0', trial_bl %in% 1:80)


ggplot(data = firstBlockData,
       aes(x = trial_bl,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  geom_smooth() +
  scale_color_manual(values = c('orange', 'blue'))

# Examine learning across all trials
ggplot(data = goData2,
       aes(x = trial,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  geom_smooth() +
  scale_color_manual(values = c('orange', 'blue')) +
  facet_wrap(~counterbalance)

# Examine learning across blocks
ggplot(data = goData2,
       aes(x = trial_bl,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  geom_smooth() +
  scale_color_manual(values = c('orange', 'blue')) +
  facet_grid(counterbalance ~ block)
ggsave("./Analysis/Plots/learning_across_blocks.jpg",
       width = 13.4,
       height = 10)

# Bin trials across all blocks
data2 <- data2 %>%
  group_by(ID) %>%
  mutate(binTrial = ntile(trial, n =4)) %>%
  ungroup()

ggplot(data = data2,
       aes(x = binTrial,
           y = RT)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = 1)) +
  scale_color_manual(values = c('orange', 'blue'))


# Binned trials by condition
data2 <- data2 %>%
  group_by(ID, condition) %>%
  mutate(binTrial = ntile(trial, n =2)) %>%
  mutate(binTrial = as.factor(binTrial)) %>%
  ungroup()


ggplot(data = data2 %>%
         group_by(ID, condition, counterbalance, binTrial) %>%
         summarise(meanRT = mean(RT)),
       aes(x = binTrial,
           y = meanRT,
           color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  scale_color_manual(values = c('orange', 'blue')) +
  facet_wrap(~ counterbalance)

firstBlockData <- firstBlockData %>%
  group_by(ID) %>%
  mutate(binTrial = ntile(trial, n = 8)) %>%
  ungroup()

ggplot(data = firstBlockData,
       aes(x = binTrial,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 0.5, aes(group = condition)) +
  scale_color_manual(values = c('orange', 'blue'))


# Regression by block
blocklm <- lm(meanRT ~ foreperiod * condition * block,
              data = summaryData2)

anova(blocklm)
Anova(blocklm)

ggplot(data = data2,
       aes(x = trial,
           y = RT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1, aes(group = condition)) +
  #stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.2) +
  scale_color_manual(values = c('orange', 'blue')) +
  facet_wrap(~ foreperiod, nrow = 2, ncol = 1)

# 2.4.2. Split participants' trias in 3 bins and compare average RTs
dataBin <- data2 %>%
  group_by(ID, condition) %>%
  mutate(trialBin)


# 2.4.3. Examine fp length order by condition
# By FP order
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = condition)) +
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('orange','blue')) +
  facet_grid(counterbalance ~ fpOrder)

# By FP order and block
ggplot(data = summaryData2,
       aes(x = foreperiod,
           y = meanRT,
           color = condition)) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", aes(group = condition)) +
  stat_summary(fun.data='mean_cl_boot',width=0.2,geom='errorbar')+
  theme(plot.title=element_text(size = rel(2), hjust = 0.5),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.text = element_text(size = rel(1.5)),
        axis.title = element_text(size = rel(1.5)))+
  scale_color_manual(values=c('orange','blue')) +
  facet_grid(counterbalance ~ block)
ggsave("./Analysis/Plots/learning_rt_by_fp_cond.jpg",
       width = 13.4,
       height = 10)
