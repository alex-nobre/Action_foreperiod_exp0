

#=================================================================#
# Rolling regression analysis
#=================================================================#

source('./Analysis/Prepare_data_4.R')
source('./Analysis/plot_theme.R')

# Load packages
library(zoo)
library(reshape2)
library(forcats)
library(lemon)


#======================== 1.  Prepare data ================================
# Import data in right format (i.e., including trials with wrong responses,
# outliers and no-go trials. These will later be dealt with within the loop)
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

#===================== 2. Rolling regression =======================
# 2.1. Functions and parameters
# Functions (adapted from Kruijne et al., 2021) and varying window
# sizes to run rolling regression

# Function to run regression model at each time point
# regr_fun <- function(z) {
#   # some conversion statements (rollapply does strange things to variable columns..)
#   dd <- z  %>% as.data.frame
#   # center
#   dd$numForeperiod %<>% as.character  %>% as.numeric %>% subtract(1.2)
#   dd$invRT %<>% as.character  %>% as.numeric
#   dd$condition %<>% as.factor
#   dd$condition  = C(dd$condition, contr.sum)
#   mod <- lm(invRT ~ numForeperiod*condition, data=dd, na.action = 'na.omit')
#   return(coef(mod))
# }

regr_fun <- function(z) {
  # some conversion statements (rollapply does strange things to variable columns..)
  dd <- z  %>% as.data.frame
  # center
  dd$numForeperiod %<>% as.character  %>% as.numeric %>% subtract(1.2)
  dd$invRT %<>% as.character  %>% as.numeric
  dd$numOneBackFP %<>% as.character  %>% as.numeric %>% subtract(1.2)
  mod <- lm(invRT ~ numForeperiod*numOneBackFP, data=dd, na.action = 'na.omit')
  return(coef(mod))
}

# use this to use different winsizes
# win = 11
win = 39
# win = 59

# Limits of rolling window for rollapply and for indices
lowlim <- win %/% 2 + 1# 1 + win/2
highlim <- 300 - (win %/% 2) # (300-win/2)+1
  
# Function to get condition values within individual participant's data
# Works only for even values of rows
getCondIdx <- function(dd) {
  ddnrows <- nrow(dd)
  startlim <- lowlim
  endlim <- highlim
  condIdx <- dd %>%
    #mutate(trial_idx=seq(n())) %>%
    slice(startlim:endlim) %>%
    #slice((win/2):(ddnrows-(win/2))) %>%
    select(condition)
  return(condIdx)
}


## Apply lm with rolling window to get beta-timecourses
goData3  %>% mutate(invRT = ifelse(Acc==1,invRT,NA), 
                    invRT = ifelse(outlier==FALSE, invRT, NA),
                    invRT = ifelse(!is.na(numOneBackFP), invRT, NA))  %>% 
  group_by(ID) %>%
  mutate(trial_idx=seq(n()))  %>% 
  do( rollapply(., width = win,
                FUN = regr_fun,
                by.column = FALSE, align = "center",partial=FALSE)  %>% 
        data.frame  %>% 
        tibble::rownames_to_column(var='idx') ) -> bigcoeff

# Get condition values for each trial and bind to data frame with beta coefficients
cond_idx <- goData3 %>%
  group_by(ID) %>%
  getCondIdx() %>%
  ungroup()

bigcoeff$condition <- cond_idx$condition

names(bigcoeff) <- c('ID','idx','intercept','fpMain','oneBackFPMain','interaction','condition')
bigcoeff$idx  %<>% as.numeric  %>% add(win %/% 2) # add win/2 size to correct trial index

# 2.2. Plot coefficients
# Here, we split the beta-coefficient dataset by condition. Although it is a variable of
# interest for the comparisons, we cannot include it in the model because it is fixed 
# across trials within a block, so the model cannot be fitted.

# 2.2.1. Action trials
# bigcoeffAction <- bigcoeff %>%
#   filter(condition=='action') %>%
#   select(-condition)

bigcoeffAction <- bigcoeff %>%
  group_by(ID) %>%
  mutate(intercept = intercept - intercept[1],
         fpMain = fpMain - fpMain[1],
         oneBackFPMain = oneBackFPMain - oneBackFPMain[1],
         interaction = interaction - interaction[1]) %>%
  ungroup() %>%
  filter(condition=='action') %>%
  select(-condition) %>%
  melt(id.vars=c('ID','idx')) %>%
  group_by(ID,variable) %>% 
  select_at(c('ID',c('idx','variable'),'value'))

#cfMeansAction <- bigcoeffAction  %>%
# bigcoeffAction <- bigcoeff %>%
#   melt(id.vars=c('ID','idx')) %>%
#   group_by(ID,variable) %>% 
#   mutate(value=value-value[1]) %>%
#   filter(condition=='action') %>%
#   select(-condition) %>%
#   select_at(c('ID',c('idx','variable'),'value'))

colnames(bigcoeffAction) <- c('group',c('idx','variable'),'dv')

bigcoeffAction <- bigcoeffAction %>%
  group_by_at(c('idx','variable')) %>%
  summarize(dv = mean(dv, na.rm=TRUE))  %>% data.frame

rownames(bigcoeffAction) <- NULL
colnames(bigcoeffAction) <- c(c('idx','variable'),'value')

head(bigcoeffAction)

cfPlotAction <- bigcoeffAction %>%
  ggplot(aes(x=idx,
             y=value)) +
  geom_hline(yintercept=0.0, size=lsz*.25, color='black') + 
  geom_line(aes(color=variable), size=lsz) + mytheme +
  theme(legend.position='none', legend.justification = c(1,0)) + 
  colbetas + fillbetas + facet_rep_grid(variable~.) + ylim(-.40,.40) + xlim(-1,300) +
  xlab("Trial") + ylab("beta coefficient") +
  labs(title='Action beta coefs - windows size = 40')

for (i in seq(4)) {
  idx = (i-1)* 75
  cfPlotAction <- cfPlotAction + geom_vline(xintercept = idx, color='black',linetype='dotted', size=0.15*lsz)
}
cfPlotAction

## Save the plot with one of the following:
ggsave('./Plots/action_beta_traces_40.svg',width = 7.5, height=12, units = 'cm')
ggsave('./Plots/action_beta_traces_12.svg',width = 7.5, height=12, units = 'cm')
ggsave('./Plots/action_beta_traces_60.svg',width = 7.5, height=12, units = 'cm')

# 2.2.2. External trials
bigcoeffExternal <- bigcoeff %>%
  group_by(ID) %>%
  mutate(intercept = intercept - intercept[1],
         fpMain = fpMain - fpMain[1],
         oneBackFPMain = oneBackFPMain - oneBackFPMain[1],
         interaction = interaction - interaction[1]) %>%
  ungroup() %>%
  filter(condition=='external') %>%
  select(-condition) %>%
  melt(id.vars=c('ID','idx')) %>%
  group_by(ID,variable) %>% 
  select_at(c('ID',c('idx','variable'),'value'))

# bigcoeffExternal <- bigcoeff %>%
#   filter(condition=='external') %>%
#   select(-condition)
# 
# cfMeansExternal <- bigcoeffExternal  %>%
#   melt(id.vars=c('ID','idx')) %>%
#   group_by(ID,variable) %>% 
#   mutate(value=value-value[1]) %>%
#   select_at(c('ID',c('idx','variable'),'value'))

colnames(bigcoeffExternal) <- c('group',c('idx','variable'),'dv')

bigcoeffExternal <- bigcoeffExternal %>%
  group_by_at(c('idx','variable')) %>%
  summarize(dv = mean(dv, na.rm=TRUE))  %>% data.frame

rownames(bigcoeffExternal) <- NULL
colnames(bigcoeffExternal) <- c(c('idx','variable'),'value')

head(bigcoeffExternal)

cfPlotExternal <- bigcoeffExternal %>%
  ggplot(aes(x=idx,
             y=value)) +
  geom_hline(yintercept=0.0, size=lsz*.25, color='black') + 
  geom_line(aes(color=variable), size=lsz) + mytheme +
  theme(legend.position='none', legend.justification = c(1,0)) + 
  colbetas + fillbetas + facet_rep_grid(variable~.) + ylim(-.40,.40) + xlim(-1,300) +
  xlab("Trial") + ylab("beta coefficient") +
  labs(title='External beta coefs - windows size = 40')

for (i in seq(4)) {
  idx = (i-1)* 75
  cfPlotExternal <- cfPlotExternal + geom_vline(xintercept = idx, color='black',linetype='dotted', size=0.15*lsz)
}
cfPlotExternal

## Save the plot with one of the following:
ggsave('./Plots/external_beta_traces_40.svg',width = 7.5, height=12, units = 'cm')
ggsave('./Plots/external_beta_traces_12.svg',width = 7.5, height=12, units = 'cm')
ggsave('./Plots/external_beta_traces_60.svg',width = 7.5, height=12, units = 'cm')

#================== Export coefficients for permutation tests in python =====================
bigcoeff  %>% 
  filter(condition=='action') %>%
  select(-condition) %>%
  write.csv('./Analysis/lm_coeff_action.csv')

bigcoeff %>%
  filter(condition=='external') %>%
  select(-condition) %>%
  write.csv('./Analysis/lm_coeff_external.csv')

############################################################################################

### add line segments for significance (Computed using python-mne, and subsequently stored in lib/ folder!)
sign_data = read.csv('lib/permtest_sign.csv')  %>% 
  filter(winwidth==win)  %>% 
  mutate( istart = tstart+win/2, 
          iend   = tend+win/2)
cfplt <- cfplt + geom_segment(
  aes(x=istart, xend=iend, color=variable), y= -0.25, yend=-0.25, size=0.5*lsz, data=sign_data)
cfplt

