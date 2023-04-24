
#=====================================================================================================================#
# Functions sourced by other scripts:
# 1. 'not in' operator
# 2. Compute z-score
#=====================================================================================================================#

# 1. "not in" operator
'%!in%' <- function(x,y)!('%in%'(x,y))


# 2. Compute z-score
compute_zscore <- function(x) {
  z_score <- (x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE)
  return(z_score)
}

### CM-SE:
# 3. baguley (2012)'s function to compute 95% conf intervals
# (from Kruijne et al., 2021)
cm.ci <- function(data.frame, conf.level = 0.95, difference = TRUE) {
  #cousineau-morey within-subject CIs
  k = ncol(data.frame)
  if (difference == TRUE) 
    diff.factor = 2^0.5/2
  else diff.factor = 1
  n <- nrow(data.frame)
  df.stack <- stack(data.frame)
  index <- rep(1:n, k)
  p.means <- tapply(df.stack$values, index, mean)
  norm.df <- data.frame - p.means + (sum(data.frame)/(n * k))
  t.mat <- matrix(, k, 1)
  mean.mat <- matrix(, k, 1)
  for (i in 1:k) t.mat[i, ] <- t.test(norm.df[i])$statistic[1]
  for (i in 1:k) mean.mat[i, ] <- colMeans(norm.df[i])
  c.factor <- (k/(k - 1))^0.5
  moe.mat <- mean.mat/t.mat * qt(1 - (1 - conf.level)/2, n - 1) * c.factor * 
    diff.factor
  ci.mat <- matrix(, k, 2)
  dimnames(ci.mat) <- list(names(data.frame), c("lower", "upper"))
  for (i in 1:k) {
    ci.mat[i, 1] <- mean.mat[i] - moe.mat[i]
    ci.mat[i, 2] <- mean.mat[i] + moe.mat[i]
  }
  ci.mat
}

# Same as above,
cm.ci.lt <- function(df, dv, group='subject_nr'){ # does cm.ci for long-table format
  nm <-  names(df)
  cvs <- nm[-which(nm==group | nm == dv)]
  md <- melt(df,measure.var=dv) # df w/o group variable
  
  wt<-dcast(md, 
            formula=as.formula(paste(group,paste(c(cvs,'variable'), collapse='+'),sep='~')),
            fun.aggregate=mean)
  ci.mat <- cm.ci(wt[,-1])
  ret.df <- data.frame(ci.mat)
  cvdf<-ldply(strsplit(rownames(ci.mat),'_'))
  ret.df <- cbind(cvdf[,-ncol(cvdf)], ret.df)
  colnames(ret.df)[1:length(cvs)] <- cvs
  return(ret.df)
}

# 4. Plotting function from Kruijne et al. (2021)
plottab <- function(dat, gv, dv='RT', group='subject_nr'){
  d  <-  dat  %>% select_at(c(group, gv, dv))
  colnames(d) <- c('group',gv, 'dv')
  # compute grand average:
  mtab <- d  %>% group_by_at(gv)  %>% 
    summarize(dv = mean(dv, na.rm=TRUE))  %>% data.frame
  #citab <- d  %>% cm.ci.lt(.,dv='dv', group='group')
  
  #aggrtab <- cbind(mtab, citab[,c(ncol(citab)-1, ncol(citab) )])
  aggrtab <- mtab
  rownames(aggrtab) <- NULL
  #colnames(aggrtab) <- c(gv, dv, 'lower','upper')
  colnames(aggrtab) <- c(gv, dv)
  return(aggrtab)
}
