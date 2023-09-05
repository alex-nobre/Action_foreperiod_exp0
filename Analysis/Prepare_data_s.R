
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
  filter(!is.na(oneBackFP), !is.na(twoBackFP))

# Save data with error trials to assess accuracy
dataAll <- data

# Keep only go trials with correct responses to analyze RT
goData <- data %>%
  filter(trialType == 'go', !is.na(RT), Acc == 1)

# Coerce foreperiod and FP n-1 back to numeric
goData$numForeperiod <- as.numeric(as.character(goData$foreperiod))
goData$numOneBackFP <- as.numeric(as.character(goData$oneBackFP))

dataAll$numForeperiod <- as.numeric(as.character(dataAll$foreperiod))
dataAll$numOneBackFP <- as.numeric(as.character(dataAll$oneBackFP))

# Create log10 of continuous independent variables
goData$numLogFP <- log10(goData$numForeperiod)
goData$logFP <- log10(goData$numForeperiod)
goData$logOneBackFP <- log10(goData$numOneBackFP)

dataAll$numLogFP <- log10(dataAll$numForeperiod)
dataAll$logFP <- log10(dataAll$numForeperiod)
dataAll$logOneBackFP <- log10(dataAll$numOneBackFP)

dataAll$squaredNumForeperiod <-  (dataAll$numForeperiod)^2
dataAll$squaredNumOneBackFP <- (dataAll$numOneBackFP)^2
dataAll$scaledNumForeperiod <-  scale(dataAll$numForeperiod, scale = FALSE)[,1]
dataAll$squaredScaledNumForeperiod <- (dataAll$scaledNumForeperiod)^2
dataAll$scaledNumOneBackFP <- scale(dataAll$numOneBackFP, scale = FALSE)[,1]

# Factor version of accuracy
goData$acc_result <- as.factor(goData$Acc)

dataAll$acc_result <- as.factor(dataAll$Acc)

# Remove extreme values
goData <- goData %>%
  filter(RT < 1.0) %>%
  filter(RT > 0.15)

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
            meanSeqEff = mean(oneBackEffect)) %>%
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
            meanSeqEff = mean(oneBackEffect)) %>%
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
            meanSeqEff = mean(oneBackEffect)) %>%
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
summaryDataAll <- dataAll %>%
  group_by(ID,foreperiod,condition,
           oneBackFP, oneBacktrialType, oneBackFPGo) %>%
  summarise(meanRT = mean(RT),
            meanAcc = mean(Acc),
            meanSeqEff = mean(oneBackEffect)) %>%
  ungroup() %>%
  mutate(numForeperiod = as.numeric(as.character(foreperiod)),
         numOneBackFP = as.numeric(as.character(oneBackFP))) %>%
  mutate(squaredNumForeperiod = numForeperiod^2,
         scaledNumForeperiod = scale(numForeperiod)[,1],
         squaredScaledNumForeperiod = scaledNumForeperiod^2,
         scaledNumOneBackFP = scale(numOneBackFP)[,1])

#write_csv(goData, "./Analysis/goData.csv")
#write_csv(goData2, "./Analysis/goData2.csv")
#write_csv(summaryData, "./Analysis/summaryData.csv")
#write_csv(summaryData2, "./Analysis/summaryData2.csv")
