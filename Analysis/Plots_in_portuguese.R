# Load necessary packages

# Read and process data
library(tidyverse)
library(broom)
library(magrittr)
library(data.table)

# Plotting
library(lattice)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(ggdist)
library(ggpubr)
library(viridis)

# Linear models
library(car)
library(codingMatrices)
library(modelr)
library(afex)
library(emmeans)
library(rtdists)

# Mixed modeling
library(lme4)
library(performance)

# Bayesian analysis
library(BayesFactor)
library(bayestestR)



# Save defaults
graphical_defaults <- par()
options_defaults <- options() 


# Prepare data 
source('./Analysis/Prepare_data_6.R')

# RT by FP and condition em portugues
RT_by_condition_port <- ggplot(data = summaryData2 %>% 
                                 group_by(ID, foreperiod, condition) %>% 
                                 summarise(meanRT = mean(meanRT)),
                               aes(x = foreperiod,
                                   y = meanRT,
                                   color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.1, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 3.9, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 3.7, width = 0.1, geom = "errorbar") + 
  labs(title = "Experimento 2",
       x = "FP (s)",
       y = "TR médio (s)",
       color = "Condição") +
  theme(plot.title = element_text(size = rel(2.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.6)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.position = "none") +
  scale_color_manual(values = c("orange", "blue"), labels = c("Externa", "Ação"))

ggplot2::ggsave("./Analysis/Plots/RT_by_condition_port.tiff",
                RT_by_condition_port,
                width = 25,
                height = 16.66,
                units = "cm")


# Sequential effects separated by FP n-1 em portugues
seqEff_by_oneback_port <- ggplot(data = summaryData2 %>%
                                   group_by(ID, foreperiod, condition, oneBackFP) %>%
                                   summarise(meanRT = mean(meanRT)),
                                 aes(x = foreperiod,
                                     y = meanRT,
                                     color=condition)) +
  geom_jitter(height = 0, width = 0.30, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 3.9, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 3.7, width = 0.1, geom = "errorbar") + 
  labs(title = "Experimento 2: efeitos sequenciais",
       x = expression("FP"[n]*" (s)"),
       y = "TR médio (s)",
       color = "Condição") +
  theme(plot.title = element_text(size = rel(2.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.0)),
        axis.title = element_text(size = rel(2.0)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(1.7)),
        legend.title = element_text(size = rel(1.8)),
        legend.text = element_text(size = rel(1.6)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1)) +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`1` = "FP[n-1] == 1.0",
                                      `1.6` = "FP[n-1] == 1.6",
                                      `2.2` = "FP[n-1] == 2.2",
                                      `2.8` = "FP[n-1] == 2.8"),
                                    default = label_parsed)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("Externa", "Ação"))

# Save
ggplot2::ggsave("./Analysis/Plots/seqEffs_exp2_port.tiff",
                seqEff_by_oneback_port,
                width = 35,
                height = 23.32,
                units = "cm")
