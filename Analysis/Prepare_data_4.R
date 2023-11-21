
#====================================================================================#
# Create dataset with numerical variables in seconds instead of milliseconds
#====================================================================================#

# Load necessary packages
library(readr)
library(ggplot2)
library(magrittr)
library(dplyr)
library(forcats)

source("./Analysis/helper_functions.R")

# Read data
data <- read_csv("./Analysis/dataActionFPAll.csv")

# Remove unnecessary columns
data <- data %>%
  dplyr::select(ID, Acc, condition, block, trialType,
                foreperiod, RT, counterbalance, 
                extFixationDuration, action_trigger.rt,
                oneBackFP, twoBackFP, oneBacktrialType, twoBacktrialType, 
                oneBackEffect)

# Coerce to factors
data <- data %>%
  mutate(across(c(ID, condition, block, trialType, foreperiod, counterbalance, oneBackFP, twoBackFP, oneBacktrialType, twoBacktrialType), as_factor))

data$condition <- data$condition %>%
  fct_relevel(c("external", "action"))

# Column for testing position based on counterbalancing
data <- data %>%
  mutate(testPos = case_when((condition == 'action' & counterbalance == 'action-external') ~ '1',
                             (condition == 'external' & counterbalance == 'action-external') ~ '2',
                             (condition == 'action' & counterbalance == 'external-action') ~ '2',
                             (condition == 'external' & counterbalance == 'external-action') ~ '1')) %>%
  mutate(testPos = as.factor(testPos))

# Remove practice trials
data <- data %>%
  filter(condition != 'practice')

# Create column for trial number
data <- data %>%
  group_by(ID) %>%
  mutate(trial = seq(1,n())) %>%
  ungroup() %>%
  group_by(ID, block) %>%
  mutate(trial_bl = seq(1, n())) %>%
  ungroup()

# Create numeric versions of foreperiod and FP n-1 (for computation of other variables)
data$numForeperiod <- as.numeric(as.character(data$foreperiod))
data$numOneBackFP <- as.numeric(as.character(data$oneBackFP))

# Create columns for numerical and categorical difference between current and previous FP, and for 
# whether the previous FP is longer than the current:
# If FPn-1 is shorter than the current FP, value is 0
# If FPn-1 is equal to the current FP, value is 0
# If FPn-1 is longer than the current FP, value is 1
data <- data %>%
  group_by(ID, block) %>%
  mutate(numOneBackFPDiff = c(NA, diff(numForeperiod))) %>%
  mutate(oneBackFPDiff = as.factor(numOneBackFPDiff), # to summarise RTs according to value of FP n-1
         squaredNumOneBackFPDiff = numOneBackFPDiff^2) %>% # quadratic term for difference between FP and FP n-1
  mutate(prevFPLonger = case_when(oneBacktrialType=="no-go" ~ "no-go",
                                  numOneBackFPDiff>0 ~ "0",
                                  numOneBackFPDiff<0 ~ "1",
                                  numOneBackFPDiff==0 ~ "0")
  ) %>%
  mutate(prevFPLonger = as.factor(prevFPLonger)) %>%
  mutate(prevFPLonger = fct_relevel(prevFPLonger, c("0", "1"))) %>%
  ungroup()

# Column for FP n-1 including no-go
data <- data %>%
  mutate(oneBackFPGo = factor(case_when(oneBacktrialType == "no-go" ~ "no-go",
                                        oneBackFP == "0.6" ~ "0.6",
                                        oneBackFP == "1.2" ~ "1.2",
                                        oneBackFP == "1.8" ~ "1.8")))

data$oneBackFPGo <- forcats::fct_relevel(data$oneBackFPGo, c("0.6", "1.2", "1.8", "no-go"))

# Remove trials without n-1 FP values (i.e., first of each block)
data <- data %>%
  filter(!is.na(oneBackFP))

# Save data with error trials to assess accuracy
dataAcc <- data

# Split between go and no-go trials for acc
dataAccGo <- dataAcc %>%
  filter(trialType == "go")

dataAccNoGo <- dataAcc %>%
  filter(trialType == "no-go")

# Keep only go trials with correct responses to analyze RT
goData <- data %>%
  filter(trialType == 'go', !is.na(RT), Acc == 1)

# Coerce foreperiod and FP n-1 back to numeric
goData$numForeperiod <- as.numeric(as.character(goData$foreperiod))
goData$numOneBackFP <- as.numeric(as.character(goData$oneBackFP))

dataAcc$numForeperiod <- as.numeric(as.character(dataAcc$foreperiod))
dataAcc$numOneBackFP <- as.numeric(as.character(dataAcc$oneBackFP))

dataAccGo$numForeperiod <- as.numeric(as.character(dataAccGo$foreperiod))
dataAccGo$numOneBackFP <- as.numeric(as.character(dataAccGo$oneBackFP))

dataAccNoGo$numForeperiod <- as.numeric(as.character(dataAccNoGo$foreperiod))
dataAccNoGo$numOneBackFP <- as.numeric(as.character(dataAccNoGo$oneBackFP))

# Create variable for error rate
dataAcc$Error <- abs(dataAcc$Acc - 1)
dataAccGo$Error <- abs(dataAccGo$Acc - 1)
dataAccNoGo$Error <- abs(dataAccNoGo$Acc - 1)

# Create log10 of continuous independent variables
goData$numLogFP <- log10(goData$numForeperiod)
goData$logFP <- log10(goData$numForeperiod)
goData$logOneBackFP <- log10(goData$numOneBackFP)

dataAcc$numLogFP <- log10(dataAcc$numForeperiod)
dataAcc$logFP <- log10(dataAcc$numForeperiod)
dataAcc$logOneBackFP <- log10(dataAcc$numOneBackFP)

dataAccGo$numLogFP <- log10(dataAccGo$numForeperiod)
dataAccGo$logFP <- log10(dataAccGo$numForeperiod)
dataAccGo$logOneBackFP <- log10(dataAccGo$numOneBackFP)

dataAccNoGo$numLogFP <- log10(dataAccNoGo$numForeperiod)
dataAccNoGo$logFP <- log10(dataAccNoGo$numForeperiod)
dataAccNoGo$logOneBackFP <- log10(dataAccNoGo$numOneBackFP)

dataAcc$squaredNumForeperiod <-  (dataAcc$numForeperiod)^2
dataAcc$squaredNumOneBackFP <- (dataAcc$numOneBackFP)^2
dataAcc$scaledNumForeperiod <-  scale(dataAcc$numForeperiod, scale = FALSE)[,1]
dataAcc$squaredScaledNumForeperiod <- (dataAcc$scaledNumForeperiod)^2
dataAcc$scaledNumOneBackFP <- scale(dataAcc$numOneBackFP, scale = FALSE)[,1]

dataAccGo$squaredNumForeperiod <-  (dataAccGo$numForeperiod)^2
dataAccGo$squaredNumOneBackFP <- (dataAccGo$numOneBackFP)^2
dataAccGo$scaledNumForeperiod <-  scale(dataAccGo$numForeperiod, scale = FALSE)[,1]
dataAccGo$squaredScaledNumForeperiod <- (dataAccGo$scaledNumForeperiod)^2
dataAccGo$scaledNumOneBackFP <- scale(dataAccGo$numOneBackFP, scale = FALSE)[,1]

dataAccNoGo$squaredNumForeperiod <-  (dataAccNoGo$numForeperiod)^2
dataAccNoGo$squaredNumOneBackFP <- (dataAccNoGo$numOneBackFP)^2
dataAccNoGo$scaledNumForeperiod <-  scale(dataAccNoGo$numForeperiod, scale = FALSE)[,1]
dataAccNoGo$squaredScaledNumForeperiod <- (dataAccNoGo$scaledNumForeperiod)^2
dataAccNoGo$scaledNumOneBackFP <- scale(dataAccNoGo$numOneBackFP, scale = FALSE)[,1]

# Factor version of accuracy/error rate
dataAcc$acc_result <- as.factor(dataAcc$Acc)
dataAccGo$acc_result <- as.factor(dataAccGo$Acc)
dataAccNoGo$acc_result <- as.factor(dataAccNoGo$Acc)

dataAcc$error_result <- as.factor(dataAcc$Error)
dataAccGo$error_result <- as.factor(dataAccGo$Error)
dataAccNoGo$error_result <- as.factor(dataAccNoGo$Error)

# Remove extreme values
ntrials_before_extrem <- nrow(goData)
goData <- goData %>%
  filter(RT < 1.0) %>%
  filter(RT > 0.15)
ntrials_after_extrem <- nrow(goData)

# Create quadratic values of continuous predictors
goData$squaredNumForeperiod <-  (goData$numForeperiod)^2
goData$squaredNumOneBackFP <- (goData$numOneBackFP)^2
goData$scaledNumForeperiod <-  scale(goData$numForeperiod, scale = FALSE)[,1]
goData$squaredScaledNumForeperiod <- (goData$scaledNumForeperiod)^2
goData$scaledNumOneBackFP <- scale(goData$numOneBackFP, scale = FALSE)[,1]

# Transform RT to reduce skew
goData$logRT <- ifelse(!is.na(goData$RT), log10(goData$RT), NA) # log-transform
goData$invRT <- ifelse(!is.na(goData$RT), 1/goData$RT, NA)

# RT Trimming
goData2 <- goData %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  filter(abs(logRTzscore) < 3) %>%
  ungroup()

# Recompute scaled variables (which were de-centered by trimming)
goData2$squaredNumForeperiod <-  (goData2$numForeperiod)^2
goData2$squaredNumOneBackFP <- (goData2$numOneBackFP)^2
goData2$scaledNumForeperiod <-  scale(goData2$numForeperiod, scale = FALSE)[,1]
goData2$squaredScaledNumForeperiod <- (goData2$scaledNumForeperiod)^2
goData2$scaledNumOneBackFP <- scale(goData2$numOneBackFP, scale = FALSE)[,1]

# Remove trials where trial type n-1 = no-go + RT trimming
goData3 <- goData %>%
  filter(oneBacktrialType == "go") %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  filter(abs(logRTzscore) < 3) %>%
  ungroup()

# Relevel oneBackFPDiff and prevFPLonger without no-go trials
goData3$oneBackFPDiff <- fct_drop(goData3$oneBackFPDiff)
goData3$prevFPLonger <- fct_drop(goData3$prevFPLonger)

# Compute scaled variables without no-go trials
goData3$squaredNumForeperiod <-  (goData3$numForeperiod)^2
goData3$squaredNumOneBackFP <- (goData3$numOneBackFP)^2
goData3$scaledNumForeperiod <-  scale(goData3$numForeperiod, scale = FALSE)[,1]
goData3$squaredScaledNumForeperiod <- (goData3$scaledNumForeperiod)^2
goData3$scaledNumOneBackFP <- scale(goData3$numOneBackFP, scale = FALSE)[,1]

# No RT trimming
goData <- goData %>%
  group_by(ID) %>%
  mutate(RTzscore=ifelse(!is.na(RT), compute_zscore(RT), NA),
         logRTzscore=ifelse(!is.na(RT), compute_zscore(logRT), NA)) %>%
  #filter(abs(logRTzscore) < 3) %>%
  ungroup()

###############################################################
# Add delay data
###############################################################
delayData <- read_csv("./Analysis/delayDataAll.csv") %>%
  mutate(across(c(ID, condition), as_factor)) %>%
  select(-condition)

goData <- inner_join(goData, delayData, by = c("trial", "ID"))
goData2 <- inner_join(goData2, delayData, by = c("trial", "ID"))
goData3 <- inner_join(goData3, delayData, by = c("trial", "ID"))
dataAcc <- inner_join(dataAcc, delayData, by = c("trial", "ID"))
dataAccGo <- inner_join(dataAccGo, delayData, by = c("trial", "ID"))
dataAccNoGo <- inner_join(dataAccNoGo, delayData, by = c("trial", "ID"))

goData <- goData %>%
  mutate(corRT = RT - delay)
goData2 <- goData2 %>%
  mutate(corRT = RT - delay)
goData3 <- goData3 %>%
  mutate(corRT = RT - delay)
dataAcc <- dataAcc %>%
  mutate(corRT = RT - delay)
dataAccGo <- dataAccGo %>%
  mutate(corRT = RT - delay)
dataAccNoGo <- dataAccNoGo %>%
  mutate(corRT = RT - delay)
################################################################

# Average data
summaryData <- goData %>%
  group_by(ID,foreperiod,logFP,condition,
           oneBackFP,twoBackFP,oneBackFPDiff,
           oneBacktrialType, oneBackFPGo,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect),
            meanCorRT = mean(corRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numOneBackFPDiff = as.numeric(as.character(oneBackFPDiff)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         squaredNumOneBackFPDiff = numOneBackFPDiff^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1],
         scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1],
         squaredScaledNumOneBackFPDiff = scaledNumOneBackFPDiff^2)

summaryData2 <- goData2 %>%
  group_by(ID,foreperiod,logFP,condition,
           oneBackFP,twoBackFP,oneBackFPDiff,
           oneBacktrialType, oneBackFPGo,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect),
            meanCorRT = mean(corRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numOneBackFPDiff = as.numeric(as.character(oneBackFPDiff)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         squaredNumOneBackFPDiff = numOneBackFPDiff^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1],
         scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1],
         squaredScaledNumOneBackFPDiff = scaledNumOneBackFPDiff^2)

summaryData3 <- goData3 %>%
  group_by(ID,foreperiod,logFP,condition,
           oneBackFP,twoBackFP,oneBackFPDiff,
           block,counterbalance) %>%
  summarise(meanRT = mean(RT),
            meanLogRT = mean(logRT),
            meanRTzscore = mean(RTzscore),
            meanInvRT = mean(invRT),
            meanSeqEff = mean(oneBackEffect),
            meanCorRT = mean(corRT)) %>%
  ungroup() %>%
  mutate(numForeperiod=as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP)),
         numOneBackFPDiff = as.numeric(as.character(oneBackFPDiff)),
         numLogFP = as.numeric(as.character(logFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         squaredNumOneBackFPDiff = numOneBackFPDiff^2,
         squaredNumLogFP = numLogFP^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1],
         scaledNumOneBackFPDiff = scale(numOneBackFPDiff)[,1],
         squaredScaledNumOneBackFPDiff = scaledNumOneBackFPDiff^2)

# Including all trials for accuracy analysis
summaryDataAcc <- dataAcc %>%
  group_by(ID,foreperiod,condition,
           oneBackFP, oneBacktrialType, oneBackFPGo) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc),
            errorRate = mean(Error),
            meanSeqEff = mean(oneBackEffect),
            meanCorRT = mean(corRT)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

summaryDataAccGo <- dataAccGo %>%
  group_by(ID,foreperiod,condition,
           oneBackFP, oneBacktrialType, oneBackFPGo) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc),
            errorRate = mean(Error),
            meanSeqEff = mean(oneBackEffect),
            meanCorRT = mean(corRT)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

summaryDataAccNoGo <- dataAccNoGo %>%
  group_by(ID,foreperiod,condition,
           oneBackFP, oneBacktrialType, oneBackFPGo) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc),
            errorRate = mean(Error),
            meanSeqEff = mean(oneBackEffect),
            meanCorRT = mean(corRT)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

#write_csv(goData, "./Analysis/goData.csv")
#write_csv(goData2, "./Analysis/goData2.csv"
#write_csv(goData2, "./Analysis/goData3.csv")
#write_csv(goData2, "./Analysis/dataAcc.csv")
#write_csv(goData2, "./Analysis/dataAccGo.csv")
#write_csv(goData2, "./Analysis/dataAccNoGo.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")
#write_csv(summaryData2, "./Analysis/summaryData3.csv")
#write_csv(summaryData2, "./Analysis/summaryDataAcc.csv")
#write_csv(summaryData2, "./Analysis/summaryDataAccGo.csv")
#write_csv(summaryData2, "./Analysis/summaryDataAccNoGo.csv")
