

library(tidyr)
actionFPData <- summaryData %>%
  select(-c(twoBackFP,oneBackFPDiff,oneBacktrialType,block,counterbalance,meanSeqEff,logMeanRT,meanRTzscore)) %>%
  group_by(ID,foreperiod,condition,oneBackFP) %>%
  summarise(meanRT=mean(meanRT))%>%
  ungroup()%>%
  pivot_wider(names_from=c(condition,foreperiod,oneBackFP),values_from=meanRT)

write_csv(actionFPData,file='G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento 1/Analysis/ActionFPData.csv')
