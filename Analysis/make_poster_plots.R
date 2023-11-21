

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

# RT by condition and fp
RT_by_condition <- ggplot(data = summaryData2 %>% 
                            group_by(ID, foreperiod, condition) %>% 
                            summarise(meanRT = mean(meanRT)),
                          aes(x = foreperiod,
                              y = meanRT,
                              color = condition)) +
  geom_jitter(height = 0, width = 0.15, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point") +
  stat_summary(fun = "mean", geom = "line", linewidth = 4.3, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 4.1, width = 0.1, geom = "errorbar") + 
  labs(title = "Experiment 1 (n = 35)",
       x = "",
       y = "Mean RT (s)",
       color = "Condition") +
  theme(plot.title = element_text(size = rel(2.8), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.position = "none") +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

# Save
ggplot2::ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/RT_by_condition_exp1.tiff",
                RT_by_condition,
                width = 25,
                height = 16.66,
                units = "cm")


# Sequential effects separated by FP n-1
seqEff_by_oneback <- ggplot(data = summaryData2 %>%
                              group_by(ID, foreperiod, condition, oneBackFP) %>%
                              summarise(meanRT = mean(meanRT)),
                            aes(x = foreperiod,
                                y = meanRT,
                                color=condition)) +
  geom_jitter(height = 0, width = 0.30, size = 3.5, alpha = 0.5) +
  stat_summary(fun = "mean", geom = "point", size = 1.5) +
  stat_summary(fun = "mean", geom = "line", linewidth = 4.3, aes(group=condition)) +
  stat_summary(fun.data = "mean_cl_boot", linewidth = 4.1, width = 0.1, geom = "errorbar") + 
  labs(title = "Experiment 1",
       x = expression("FP"[n]*" (s)"),
       y = "Mean RT (s)",
       color = "Condition") +
  theme(plot.title = element_text(size = rel(3.0), hjust = 0.5, margin = margin(t = 0, r = 0, b = 25, l = 0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = rel(2.8)),
        axis.title = element_text(size = rel(2.8)),
        axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 2.75, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        axis.text = element_text(size = rel(2.5)),
        legend.title = element_text(size = rel(2.6)),
        legend.text = element_text(size = rel(2.4)),
        legend.key = element_blank(),
        legend.box.spacing = unit(0, "pt"),
        legend.margin = margin(5.5, 5.5, 5.5, 1),
        legend.position = "none") +
  facet_wrap(~oneBackFP,
             labeller = as_labeller(c(`0.6` = "FP[n-1] == 0.6",
                                      `1.2` = "FP[n-1] == 1.2",
                                      `1.8` = "FP[n-1] == 1.8"),
                                    default = label_parsed)) +
  scale_color_manual(values = c("orange", "blue"), labels = c("External", "Action"))

# Save
ggplot2::ggsave("G:/My Drive/Post-doc/Eventos/TRF-3/Poster/seqEffs_exp1.tiff",
                seqEff_by_oneback,
                width = 35,
                height = 23.32,
                units = "cm")
