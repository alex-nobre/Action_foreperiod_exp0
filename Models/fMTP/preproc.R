library(tools)
library(readr)
library(MASS)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(lubridate)
library(magrittr)
library(lme4)
library(BayesFactor)
library(knitr)
library(MuMIn)
library(lmerTest)
library(gridExtra)
library(emmeans)
library(lemon)
library(zoo)
library(broom)
library(tidyr)
library(xtable)

setwd("G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/data")

subfiles <- list.files()

lapply(subfiles, read_csv)
for (i in 1:length(subfiles)) assign(subfiles[i], read_csv(subfiles[i]))

preproc <- function(dat) {

	dat  %<>% select(logfile, subject_nr, Block_nr, blocktype, 
		ITI, image=im_basename,  cue_type, foreperiod, 
		correct, response, correct_response, RT=response_time, 
		longcue, shortcue, response_notice_anything, response_n_FPs, response_kb_mcpairing)
	dat %<>% mutate( sub_id = logfile  %>% file_path_sans_ext %>% basename ) 
	dat  %<>% mutate( exp_half = as.integer(Block_nr > 5) )

	dat$cue_type  %<>% as.character
	dat$longcue   %<>% as.character
	dat  %<>% mutate( 
		distribution_assoc=ifelse(cue_type==longcue,'anti-exp','exponential') )

	# dat$exp_half = as.integer(dat$)
	# find the M and SD for outlier removal -- using only correct, not_practice trials
    dat$outlier = FALSE
    dat$logRT <- log(dat$RT)
    M = dat  %>% filter(correct==1)  %>% .$logRT  %>% mean
    S = dat  %>% filter(correct==1)  %>% .$logRT  %>% sd
    dat$zlogRT <- (dat$logRT - M)/S
    dat$outlier <- abs(dat$zlogRT) > 3.0 
	return( dat )
}



# ProjectTemplate loads all raw datafiles as variables in the workspace; 
# mget can get variables by name
# so: fetch datasets by variable name, preprocess them, and merge into an all_dat table
patt <- '[1-9]+'
sub.dsets <- ls(pattern=patt)
sub.dsets <- mget(sub.dsets)
all_dat <- sub.dsets %>% ldply(preproc)#, progress='text')

summaryData <- all_data %>%
  

cache('all_dat')

#### some grand-average statistics:
gave_stats = all_dat  %>% filter(Block_nr > 0)  %>% group_by(sub_id)  %>% 
	summarize(acc = mean(correct), RT = mean(RT[correct==TRUE]) )  %>% 
	mutate(zRT = scale(RT)  %>%  as.numeric)

print("participants with low accuracies:")
acctab <- gave_stats  %>% filter(acc < .95)  %T>%  print   

print("participants with slow RT-responses:")
slowpptab <- gave_stats  %>% filter(abs(zRT)  > 2.5)  %T>%  print 

#### mark outliers:
all_dat   %<>% group_by(sub_id)  %>% mutate(
		bad_sub = (sub_id %in% acctab$sub_id) | (sub_id %in% slowpptab$sub_id)     )

# proportion trials discarded with a cutoff of 800ms:
d<- all_dat  %>% filter(Block_nr > 0, correct==1)  %>%  
	group_by(sub_id)  %>%
	summarize(P_SlowCorrect = mean(RT > 800) )
# d$sub_id = factor(d$sub_id, 
# 						levels = sort(d$P_SlowCorrect, decreasing=TRUE, index.return=TRUE)$ix )
d  %>% 
	ggplot(aes(sub_id, P_SlowCorrect)) + geom_bar(stat='identity') + theme_bw() + ggtitle("fraction trials discarded")


