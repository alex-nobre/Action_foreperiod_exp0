

# Calculate sd point-by-point in (time)series

data <- read.csv('./Analysis/lm_coeff_external_split.csv')

str(data)
data$ID = as.factor(data$ID)
data$idx = as.factor(data$idx)

library(ggplot2)

ggplot(data, aes(x=idx, y= intercept, group = ID, color = ID))+
    geom_line()
#z[, c("a", "b")]

data_wide <- data %>%
  pivot_wider(names_from = variable,
              values_from = value)

datatemp = data[,c('ID','intercept', 'idx')]


#Mean
data_intercept <- data.frame('intercept' =  colMeans(wide_intercept[,2:ncol(wide_intercept)]))

#CI
wide_intercept <- reshape(datatemp, v.names = 'intercept', idvar = 'ID', 
                          timevar = 'idx', direction = 'wide')
data_intercept$CI =  1.96*(apply(wide_intercept[,2:ncol(wide_intercept)], 
                                 2, sd, na.rm = TRUE)/sqrt(NROW(wide_intercept)))
data_intercept$idx = unique(datatemp$idx)


ggplot(data_intercept, aes(x=idx, y = intercept))+
    geom_line(color = 'red')+
    geom_ribbon(aes(ymin = intercept-CI, ymax=intercept+CI), alpha = 0.5)



