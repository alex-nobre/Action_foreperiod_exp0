

#================================================================================#
#================================================================================#

# Load necessary packages

# Read and process data
library(readr)
library(magrittr)
library(dplyr)
library(data.table)

# Plotting
library(ggplot2)
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggpubr)

# Linear models
library(car)
library(codingMatrices)

# Mixed models
library(afex)
library(emmeans)
library(lme4)
library(performance)
library(ggsignif)


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

# Check for influence of latency of action key press on RT
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

# Plot grand means
RT_by_condition <- ggplot(data = summaryData2 %>% 
                               group_by(ID, foreperiod, condition) %>% 
                               summarise(meanRT = mean(meanRT)),
                             aes(x = foreperiod,
                                 y = meanRT,
                                 color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.2, width = 0.1, geom = "errorbar") + 
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
ggplot2::ggsave("./Analysis/Plots/RT_by_condition.pdf",
                RT_by_condition,
                width = 15,
                height = 10)


# Histograms
rt_fp_cond_hists <- ggplot(goData2, aes(x = RT)) +
  geom_histogram() +
  facet_grid(foreperiod ~ condition)


# Boxplot
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
ggsave("./Analysis/Plots/SeqEff.jpeg",
#ggsave("./Analysis/Plots/SeqEff.pdf",
       width = 6.7,
       height = 5)

# Sequential effects separated by FP n-1
seqEff_by_oneback <- ggplot(data = summaryData2 %>%
                              group_by(ID, foreperiod, condition, oneBackFP) %>%
                              summarise(meanRT = mean(meanRT)),
                            aes(x = foreperiod,
                                y = meanRT,
                                color=condition)) +
  geom_jitter(height = 0, width = 0.30, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.2, width = 0.1, geom = "errorbar") + 
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
  facet_wrap(~oneBackFP) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

ggsave("./Analysis/Plots/SeqEff.pdf",
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
acc_plot <- ggplot(data = summaryDataAcc %>% 
                     group_by(ID, foreperiod, condition) %>%
                     summarise(errorRate = mean(errorRate)),
                   aes(x = foreperiod,
                       y = errorRate,
                       color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.2, width = 0.1, geom = "errorbar") +
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
ggsave("./Analysis/Plots/acc_by_condition.png",
       acc_plot,
       width = 15,
       height = 10)

# Accuracy go
acc_go_plot <- ggplot(data = summaryDataAccGo %>% 
         group_by(ID, foreperiod, condition) %>%
         summarise(errorRate = mean(errorRate)),
       aes(x = foreperiod,
           y = errorRate,
           color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.2, width = 0.1, geom = "errorbar") + 
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
ggsave("./Analysis/Plots/acc_by_condition_go.png",
       acc_go_plot,
       width = 15,
       height = 10)

# Accuracy no go
acc_nogo_plot <- ggplot(data = summaryDataAccNoGo %>% 
         group_by(ID, foreperiod, condition) %>%
         summarise(errorRate = mean(errorRate)),
       aes(x = foreperiod,
           y = errorRate,
           color = condition)) +
  geom_jitter(height = 0, width = 0.15, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 1.4, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 1.2, width = 0.1, geom = "errorbar") + 
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
ggsave("./Analysis/Plots/acc_by_condition_nogo.png",
       acc_nogo_plot,
       width = 15,
       height = 10)

# Acuracy plots in single panel
# Save legend as grob
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(acc_nogo_plot + 
                                                          theme(legend.title = element_text(size = rel(1.1)),
                                                                legend.text = element_text(size = rel(0.9))))), "guide-box")

# Visualize
grid.arrange(acc_plot + theme(legend.position = "none",
                              axis.text = element_text(size = rel(1.2)),
                              axis.title = element_text(size = rel(1.4)),
                              plot.title = element_text(size = rel(1.5))),
             acc_go_plot + theme(legend.position = "none",
                                 axis.text = element_text(size = rel(1.2)),
                                 axis.title = element_text(size = rel(1.4)),
                                 axis.title.y = element_blank(),
                                 plot.title = element_text(size = rel(1.5))),
             acc_nogo_plot + theme(legend.position = "none",
                                   axis.text = element_text(size = rel(1.2)),
                                   axis.title = element_text(size = rel(1.4)),
                                   axis.title.y = element_blank(),
                                   plot.title = element_text(size = rel(1.5))),
             cond_legend,
             nrow = 1,
             widths = c(3/10, 3/10, 3/10, 1/10))

# Save
acc_all_plots <- arrangeGrob(acc_plot + theme(legend.position = "none",
                                              axis.text = element_text(size = rel(1.2)),
                                              axis.title = element_text(size = rel(1.4)),
                                              plot.title = element_text(size = rel(1.5))),
                             acc_go_plot + theme(legend.position = "none",
                                                 axis.text = element_text(size = rel(1.2)),
                                                 axis.title = element_text(size = rel(1.4)),
                                                 axis.title.y = element_blank(),
                                                 plot.title = element_text(size = rel(1.5))),
                             acc_nogo_plot + theme(legend.position = "none",
                                                   axis.text = element_text(size = rel(1.2)),
                                                   axis.title = element_text(size = rel(1.4)),
                                                   axis.title.y = element_blank(),
                                                   plot.title = element_text(size = rel(1.5))),
                             cond_legend,
                             nrow = 1,
                             widths = c(3/10, 3/10, 3/10, 1/10))

ggsave("./Analysis/Plots/acc_all_plots.pdf",
       acc_all_plots,
       width = 12,
       height = 6.67)

# RT and accuracy in single panel
# Save legend as grob
cond_legend <- gtable_filter(ggplot_gtable(ggplot_build(acc_nogo_plot + 
                                                          theme(legend.title = element_text(size = rel(1.1)),
                                                                legend.text = element_text(size = rel(0.9))))), "guide-box")

# Visualize
grid.arrange(RT_by_condition + labs(tag = "A", x = "") +
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(1.1)),
                     axis.title = element_text(size = rel(1.3)),
                     plot.title = element_text(size = rel(1.4)),
                     plot.tag = element_text(size = rel(1.5))),
             acc_go_plot + labs(tag = "B") +
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(1.1)),
                     axis.title = element_text(size = rel(1.3)),
                     plot.title = element_text(size = rel(1.4)),
                     plot.tag = element_text(size = rel(1.5))),
             acc_nogo_plot + labs(tag = "C", x = "") + 
               theme(legend.position = "none",
                     axis.text = element_text(size = rel(1.1)),
                     axis.title = element_text(size = rel(1.3)),
                     plot.title = element_text(size = rel(1.4)),
                     plot.tag = element_text(size = rel(1.5))) +
               labs(y = ""),
             cond_legend,
             nrow = 1,
             widths = c(5/17, 5/17, 5/17, 2/17))


rt_acc_plots <- arrangeGrob(RT_by_condition + labs(tag = "A", x = "") +
                              theme(legend.position = "none",
                                    axis.text = element_text(size = rel(1.1)),
                                    axis.title = element_text(size = rel(1.3)),
                                    plot.title = element_text(size = rel(1.4)),
                                    plot.tag = element_text(size = rel(1.5))),
                            acc_go_plot + labs(tag = "B") +
                              theme(legend.position = "none",
                                    axis.text = element_text(size = rel(1.1)),
                                    axis.title = element_text(size = rel(1.3)),
                                    plot.title = element_text(size = rel(1.4)),
                                    plot.tag = element_text(size = rel(1.5))),
                            acc_nogo_plot + labs(tag = "C", x = "") + 
                              theme(legend.position = "none",
                                    axis.text = element_text(size = rel(1.1)),
                                    axis.title = element_text(size = rel(1.3)),
                                    plot.title = element_text(size = rel(1.4)),
                                    plot.tag = element_text(size = rel(1.5))) +
                              labs(y = ""),
                            cond_legend,
                            nrow = 1,
                            widths = c(5/17, 5/17, 5/17, 2/17))

ggsave("./Analysis/Plots/rt_acc_plots.pdf",
       rt_acc_plots,
       width = 20,#16,
       height = 8.33,#11.94,
       units = "cm")

# For poster
ggplot2::ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/RT_by_condition.tiff",
                RT_by_condition + labs(title = "Experiment 1") +
                  theme(plot.title = element_text(size = rel(3.0)),
                        axis.title = element_text(size = rel(2.6)),
                        axis.text = element_text(size = rel(2.3)),
                        legend.title = element_text(size = rel(2.4)),
                        legend.text = element_text(size = rel(2.2)),
                        legend.key = element_blank(),
                        legend.box.spacing = unit(0, "pt"),
                        legend.margin = margin(5.5, 5.5, 5.5, 1)),
                width = 35,
                height = 23.33,
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


check_sphericity(fpAnova)

is_norm <- check_normality(fpAnova)

plot(is_norm)

plot(is_norm, type = "qq")

nice(fpAnova)
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




logFPAnova <- aov_ez(id = "ID",
                     dv = "meanLogRT",
                     data = summaryData2,
                     within = c("foreperiod", "condition"))

logFPAnova

#======================= 2.2. Sequential effects ============================
#========== 2.2.1. RT ============

# 2.2.1. Anova for FP n-1
seqEffAnova <- aov_ez(id = "ID",
                  dv = "meanRT",
                  data = summaryData2,
                  within = c("foreperiod", "condition", "oneBackFP"),
                  anova_table = list(es = "pes"))

# Without no-go trials
seqEffAnova <- aov_ez(id = "ID",
                      dv = "meanRT",
                      data = summaryData3,
                      within = c("foreperiod", "condition", "oneBackFP"),
                      anova_table = list(es = "pes"))

check_sphericity(seqEffAnova)

is_norm <- check_normality(seqEffAnova)

nice(seqEffAnova)

# Pairwise comparisons for conditions wihtin levels of foreperiod
conditionEmmeans <- emmeans(seqEffAnova,
                     pairwise ~ condition|foreperiod,
                     adjust = "holm")

# fpEmmeansContrasts <- contrast(fpEmmeans[[1]],
#                                interaction=c('poly'),
#                                adjust='bonferroni')


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


#================= 2.2.3. N-1 trial type, condition and sequential effects =================

# 2.2.3.1. Anova with foreperiod, condition and n-1 trial type
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

# 2.2.3.2. Anova with foreperiod, n-1 trial type and sequential effects
seqEfftrialTypeAnova <- aov_ez(id = "ID",
                         dv = "meanRT",
                         data = summaryData,
                         within = c("foreperiod", "oneBackFP", "oneBacktrialType"))

nice(seqEfftrialTypeAnova,
     correction='none')

# 2.2.3.3. Anova with foreperiod, condition, n-1 trial type and sequential effects
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


#==================== 2.4. Foreperiod, condition and block ======================
blocklm <- lm(meanRT ~ foreperiod * counterbalance * block,
              data = summaryData)

anova(blocklm)
Anova(blocklm)


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
