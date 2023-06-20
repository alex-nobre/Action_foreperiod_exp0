

#===================================================================#
# Rolling regression analysis
# Runs rolling regression separately for each counterbalancing order
# Because participants went through distinct experimental conditions
# in distinct orders, their learning trajectories are not comparable
# Thus, here we only extract coefficients for the first half of the 
# experiment, in which they have not been exposed to different conditions
# before

# Changes:
# Uses RT instead of inverse RT
#===================================================================#

# Load packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(zoo)
library(reshape2)
library(forcats)
library(lemon)
library(tidyr)

#source('./Analysis/Prepare_data_4.R')
source('./Analysis/plot_theme.R')
source("./Analysis/helper_functions.R")

#======================== 1.  Prepare data ================================
# Create dataset including trials with wrong responses,
# outliers and no-go trials. These will later be dealt with within the loop)

# Read data
data <- read_csv("./Analysis/dataActionFPAll.csv")

# Remove unnecessary columns
data <- data %>%
  dplyr::select(ID, Acc, condition, block, trialType,
                foreperiod, RT, counterbalance, 
                extFixationDuration, action_trigger.rt,
                oneBackFP, twoBackFP, oneBacktrialType, twoBacktrialType, 
                oneBackEffect) %>%
  mutate(foreperiod = foreperiod * 1000,
         RT = RT *1000,
         extFixationDuration = extFixationDuration * 1000,
         action_trigger.rt = action_trigger.rt * 1000,
         oneBackFP = oneBackFP * 1000,
         twoBackFP = twoBackFP * 1000,
         oneBackEffect = oneBackEffect * 1000)

# Coerce to factors
data$ID <- as.factor(data$ID)
data$condition <- data$condition %>%
  as.factor() %>%
  forcats::fct_relevel(c("external", "action"))
data$block <- as.factor(data$block)
data$trialType <- as.factor(data$trialType)
data$foreperiod <- as.factor(data$foreperiod)
data$counterbalance <- as.factor(data$counterbalance)
data$oneBackFP <- as.factor(data$oneBackFP)
data$twoBackFP <- as.factor(data$twoBackFP)
data$oneBacktrialType <- as.factor(data$oneBacktrialType)
data$twoBacktrialType <- as.factor(data$twoBacktrialType)


# Remove practice trials
data <- data %>%
  filter(condition != 'practice')


# Create column for difference between current and previous FP:
# If FPn-1 is shorter than the current FP, value is 1
# If FPn-1 is equal to the current FP, value is 1
# If FPn-1 is longer than the current FP, value is -1
# If the previous trial was a no-go trial, use this as separate category
oneBackFPDiff <- diff(as.numeric(as.character(data$foreperiod)))
oneBackFPDiff <- c(NA, oneBackFPDiff)
data$numOneBackFPDiff <- oneBackFPDiff


data <- data %>%
  mutate(oneBackFPDiff = case_when(oneBacktrialType=="no-go" ~ "no-go",
                                   numOneBackFPDiff>0 ~ "1",
                                   numOneBackFPDiff<0 ~ "-1",
                                   numOneBackFPDiff==0 ~ "0")
         ) %>%
  mutate(oneBackFPDiff = as.factor(oneBackFPDiff)) %>%
  mutate(oneBackFPDiff=forcats::fct_relevel(oneBackFPDiff, c("no-go", "-1", "0", "1")))

# Column for FP n-1 including no-go
data <- data %>%
  mutate(oneBackFPGo = factor(case_when(oneBacktrialType == "no-go" ~ "no-go",
                                 oneBackFP == "600" ~ "600",
                                 oneBackFP == "1200" ~ "1200",
                                 oneBackFP == "1800" ~ "1800")))

data$oneBackFPGo <- forcats::fct_relevel(data$oneBackFPGo, c("600", "1200", "1800", "no-go"))
  
# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter(!is.na(oneBackFP))

goData3 <- data
  
# Coerce foreperiod and FP n-1 back to numeric
goData3$numForeperiod <- as.numeric(as.character(goData3$foreperiod))
goData3$numOneBackFP <- as.numeric(as.character(goData3$oneBackFP))

# Transform RT to reduce skew
# goData3$logRT <- ifelse(!is.na(goData3$RT), log10(goData3$RT), NA) # log-transform
# goData3$invRT <- ifelse(!is.na(goData3$RT), 1/goData3$RT, NA)

goData3 <- goData3 %>%
  mutate(logRT=ifelse(!is.na(RT), log10(RT), NA),
         invRT=ifelse(!is.na(RT), 1/RT, NA))

# Flag outliers
goData3 <- goData3 %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  mutate(outlier = ifelse(abs(logRTzscore) > 3,TRUE,FALSE)) %>%
  ungroup()

#===================== 2. Functions and parameters =======================
regr_fun <- function(z) {
  # some conversion statements (rollapply does strange things to variable columns..)
  dd <- z  %>% as.data.frame
  # center predictors and make sure FP and FP n-1 are numeric (to get effect as a single coefficient)
  dd$numForeperiod %<>% as.character  %>% as.numeric %>% subtract(1200)
  dd$RT %<>% as.character  %>% as.numeric
  dd$numOneBackFP %<>% as.character  %>% as.numeric %>% subtract(1200)
  # Fit linear model to RTs in window as a function of FP and FP n-1
  mod <- lm(RT ~ numForeperiod*numOneBackFP, data=dd, na.action = 'na.omit')
  return(coef(mod))
}


# use this to use different winsizes
# win = 11
win = 39
#win = 59

# Get number of trials by subject to define low and high limits
ntrials <- goData3 %>%
  group_by(ID) %>%
  summarise(n()) %>%
  ungroup() %>%
  rename(ntrials = `n()`) %>%
  select(ntrials)

# Check that all subs have the same number of trials
length(unique(ntrials$ntrials))==1

ntrials <- unique(ntrials$ntrials)

# Limits of rolling window for rollapply and for indices
lowlim <- win %/% 2 + 1# 1 + win/2
highlim <- ntrials - (win %/% 2) # (300-win/2)+1
  
# Function to get condition values within individual participant's data
# Works only for even values of rows
getCond <- function(dataset) {
  nrowsdataset <- nrow(dataset)
  startlim <- lowlim
  endlim <- highlim
  cond <- dataset %>%
    slice(startlim:endlim) %>%
    select(condition)
  return(cond)
}

#===================== 3. Get action coefs =======================
# Get coefficients for the first half of the experimento for participants
# in the "action-external" counterbalancing order

# Filter participants by counterbalancing order
goData4 <- goData3 %>%
  filter(counterbalance=='action-external')

## Apply lm with rolling window get beta-timecourses
goData4  %>% mutate(RT = ifelse(Acc==1,RT,NA), 
                    RT = ifelse(outlier==FALSE, RT, NA),
                    RT = ifelse(!is.na(numOneBackFP), RT, NA))  %>% 
  filter(condition=='action') %>%
  select(-condition) %>%
  group_by(ID) %>%
  mutate(trial_idx=seq(n()))  %>% # Create column for trial index
  # Fit regression function across rolling windows
  do( rollapply(., width = win,
                FUN = regr_fun,
                by.column = FALSE, align = "center",partial=FALSE)  %>% 
        data.frame  %>% 
        tibble::rownames_to_column(var='idx') ) -> bigcoeff

# Rename columns
names(bigcoeff) <- c('ID','trial','intercept','fpMain','oneBackFPMain','interaction')
bigcoeff$trial  %<>% as.numeric  %>% add(win %/% 2) # add win/2 size to correct trial index

# Reshape data to "long" format by ID and trial index and subtract first value from each coefficient
bigcoeffAction <- bigcoeff %>%
  melt(id.vars=c('ID','trial')) %>% 
  group_by(ID, variable) %>%
  mutate(value = value - value[1]) %>%
  ungroup()

colnames(bigcoeffAction) <- c('ID',c('trial','variable'),'dv')

# Compute means and confidence intervals for each time point
bigcoeffAction <- bigcoeffAction %>%
  group_by(trial, variable) %>%
  summarise(meanDV = mean(dv, na.rm = TRUE),
            ciUpper = mean(dv, na.rm=TRUE) + 1.96 * (sd(dv, na.rm = TRUE)/sqrt(n())),
            ciLower = mean(dv, na.rm=TRUE) - 1.96 * (sd(dv, na.rm = TRUE)/sqrt(n()))) %>% 
  data.frame

# Plot coefficients by predictor
cfPlotAction <- bigcoeffAction %>%
  ggplot(aes(x=trial,
             y=meanDV)) +
  geom_hline(yintercept=0.0, linewidth=lsz*.25, color='black') + 
  geom_ribbon(aes(ymin = ciLower, ymax = ciUpper, fill = variable), alpha = 0.2, size = 0) +
  geom_line(aes(color=variable), linewidth=lsz) + 
  mytheme +
  theme(legend.position='none', legend.justification = c(1,0)) + 
  colbetas + fillbetas + 
  facet_rep_grid(variable~., scales = "free") + 
  xlim(-1,300) +
  labs(x = "Trial",
       y = "beta coefficient",
       title = paste('Action beta coefs/window size =', win + 1))


# Lines to separate blocks
for (i in seq(4)) {
  idx = (i-1)* 75
  cfPlotAction <- cfPlotAction + geom_vline(xintercept = idx, color='black', linetype='dotted', linewidth=0.15*lsz)
}
cfPlotAction

# Save the plot
ggsave(paste("./Plots/Rolling_regression/Split_dataset/action_beta_split_rt_", win+1, ".tiff", sep = "")
       ,width = 7.5, height=12, units = 'cm')

#===================== 4. Get external coefs =======================
# Now we get coefficients for the second half of the experiment for participants
# in the "external-action" counterbalancing order

# Filter participants by counterbalancing order
goData5 <- goData3 %>%
  filter(counterbalance=='external-action')

## Apply lm with rolling window get beta-timecourses
goData5  %>% mutate(RT = ifelse(Acc==1,RT,NA), 
                    RT = ifelse(outlier==FALSE, RT, NA),
                    RT = ifelse(!is.na(numOneBackFP), RT, NA))  %>% 
  filter(condition=='external') %>%
  select(-condition) %>%
  group_by(ID) %>%
  mutate(trial_idx=seq(n()))  %>% # Create column for trial index
  # Fit regression function across rolling windows
  do( rollapply(., width = win,
                FUN = regr_fun,
                by.column = FALSE, align = "center",partial=FALSE)  %>% 
        data.frame  %>% 
        tibble::rownames_to_column(var='idx') ) -> bigcoeff

# Rename columns
names(bigcoeff) <- c('ID','trial','intercept','fpMain','oneBackFPMain','interaction')
bigcoeff$trial  %<>% as.numeric  %>% add(win %/% 2) # add win/2 size to correct trial index


# Reshape data to "long" format by ID and trial index and subtract first value from each coefficient
bigcoeffExternal <- bigcoeff %>%
  melt(id.vars=c('ID','trial')) %>% 
  group_by(ID, variable) %>%
  mutate(value = value - value[1]) %>%
  ungroup()

colnames(bigcoeffExternal) <- c('ID',c('trial','variable'),'dv')

# Compute means and confidence intervals for each time point
bigcoeffExternal <- bigcoeffExternal %>%
  group_by(trial, variable) %>%
  summarise(meanDV = mean(dv, na.rm = TRUE),
            ciUpper = mean(dv, na.rm=TRUE) + 1.96 * (sd(dv, na.rm = TRUE)/sqrt(n())),
            ciLower = mean(dv, na.rm=TRUE) - 1.96 * (sd(dv, na.rm = TRUE)/sqrt(n()))) %>% 
  data.frame


# Plot coefficients by predictor
cfPlotExternal <- bigcoeffExternal %>%
  ggplot(aes(x=trial,
             y=meanDV)) +
  geom_hline(yintercept=0.0, linewidth=lsz*.25, color='black') + 
  geom_ribbon(aes(ymin = ciLower, ymax = ciUpper, fill = variable), alpha = 0.2, size = 0) +
  geom_line(aes(color=variable), linewidth=lsz) +
  mytheme +
  theme(legend.position='none', legend.justification = c(1,0)) + 
  colbetas + fillbetas + 
  facet_rep_grid(variable~., scales = "free") + 
  xlim(-1,300) +
  labs(x = "Trial", 
       y = "beta coefficient", 
       title = paste('External beta coefs/window size =', win + 1))

# Lines to separate blocks
for (i in seq(4)) {
  idx = (i-1)* 75
  cfPlotExternal <- cfPlotExternal + geom_vline(xintercept = idx, color='black',linetype='dotted', linewidth=0.15*lsz)
}
cfPlotExternal

## Save the plot
ggsave(paste("./Plots/Rolling_regression/Split_dataset/external_beta_split_rt_", win+1, ".tiff", sep = "")
       ,width = 7.5, height=12, units = 'cm')


####### For the permutation test, the python packages mne, pandas and numpy are used:
bigcoeffAction  %>% write.csv('./Analysis/lm_coeff_action_split_rt.csv')
bigcoeffExternal  %>% write.csv('./Analysis/lm_coeff_external_split_rt.csv')

############################################################################################

### add line segments for significance (Computed using python-mne, and subsequently stored in lib/ folder!)
sign_data = read.csv('lib/permtest_sign.csv')  %>% 
  filter(winwidth==win)  %>% 
  mutate( istart = tstart+win/2, 
          iend   = tend+win/2)
cfplt <- cfplt + geom_segment(
  aes(x=istart, xend=iend, color=variable), y= -0.25, yend=-0.25, size=0.5*lsz, data=sign_data)
cfplt


