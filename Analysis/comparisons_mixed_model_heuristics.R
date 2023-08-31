

# Load necessary packages
library(tidyverse)
library(magrittr)
library(lattice)
library(afex)
library(emmeans)
library(lme4)
library(car)
library(data.table)
library(buildmer)
library(tidyr)
library(janitor)
library(broom.mixed)


# Save defaults
graphical_defaults <- par()
options_defaults <- options() 

# Read data
source("./Analysis/Prepare_data_4.R")

source("./Analysis/mixed_model_fitting_code.R")

#================================= 1. Buildmer x complete model ===================================

buildmerMod <- mixed(formula = RT ~  1 + numForeperiod + oneBackFP + numForeperiod:oneBackFP + 
                         condition + numForeperiod:condition + oneBackFP:condition + 
                         numForeperiod:oneBackFP:condition + (1 + condition | ID),
                       data=goData2,
                       control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                       progress = TRUE,
                       expand_re = TRUE,
                       method =  'KR',
                       REML=TRUE,
                       return = "merMod")

summary(buildmerMod)


fullMod <- mixed(formula = RT ~  numForeperiod * oneBackFP * condition +
                   (1 + numForeperiod * oneBackFP * condition | ID),
                 data=goData2,
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod")

anova(buildmerMod)
anova(fullMod)


#================================= 2. Buildmer x keep it maximal ===================================

# Keep it maximal
kimMod <- mixed(formula = RT ~  numForeperiod * oneBackFP * condition +
                   (1 + numForeperiod * oneBackFP * condition | ID),
                 data=goData2,
                 control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                 progress = TRUE,
                 expand_re = TRUE,
                 method =  'KR',
                 REML=TRUE,
                 return = "merMod")

# Remove correlations
kimMod <- mixed(formula = RT ~  numForeperiod * oneBackFP * condition +
                  (1 + numForeperiod * oneBackFP * condition || ID),
                data=goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'KR',
                REML=TRUE,
                return = "merMod")

# Remove interactions
kimMod <- mixed(formula = RT ~  numForeperiod * oneBackFP * condition +
                  (1 + numForeperiod + oneBackFP + condition || ID),
                data=goData2,
                control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
                progress = TRUE,
                expand_re = TRUE,
                method =  'KR',
                REML=TRUE,
                return = "merMod")

summary(kimMod)
anova(kimMod)

#================================= 3. Buildmer x alternative heuristic ===================================

ahMod <- mixed(formula = RT ~  numForeperiod * oneBackFP * condition +
                 (1 + numForeperiod * oneBackFP * condition || ID),
               data=goData2,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'KR',
               REML=TRUE,
               return = "merMod")

rmvdSlopes <- pare.random.slopes(ahMod)

ahMod <- mixed(formula = RT ~  numForeperiod * oneBackFP * condition +
                 (1 + condition || ID),
               data=goData2,
               control = lmerControl(optimizer = c("bobyqa"),optCtrl=list(maxfun=2e5),calc.derivs = FALSE),
               progress = TRUE,
               expand_re = TRUE,
               method =  'KR',
               REML=TRUE,
               return = "merMod")

anova(ahMod)
