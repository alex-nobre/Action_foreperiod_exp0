## load relevant packages and identify directory where sample file is stored
library(lmerTest)

#file_dir <- "~/Downloads/"   # the directory where the sample file (gloves_Exp1a_usable_data.csv) is stored


## define functions for subsequent use

compute.lmerTest.dfs                      <- function(fitted.model, formula, data = NULL, REML = TRUE, control = lmerControl(), 
                                                      start = NULL, verbose = 0L, subset, weights, na.action, 
                                                      offset, contrasts = NULL, devFunOnly = FALSE) {
  # this function is equivalent to lmerTest::lmer(), but with one key difference:
  # it does not actually fit a model. it takes a fitted lme4::lmer() model
  # (fitted.model) plus all of the necessary arguments for the model call, then
  # runs the df-computing/significance testing functions on the fitted model.
  
  #' @rawNamespace
  #' if(getRversion() >= "3.3.0") {
  #'   importFrom("stats", sigma)
  #' } else {
  #'   export(sigma)
  #' }
  #'
  if(getRversion() < "3.3") {
    sigma <- function(object, ...) UseMethod("sigma")
    
    sigma.merMod <- function (object, ...)
    {
      dc <- object@devcomp
      dd <- dc$dims
      if (dd[["useSc"]])
        dc$cmp[[if (dd[["REML"]])
          "sigmaREML"
          else "sigmaML"]]
      else 1
    }
  }
  
  ##############################################
  ######## as_lmerModLT()
  ##############################################
  as_lmerModLT <- function(model, devfun, tol=1e-8) {
    is_reml <- getME(model, "is_REML")
    # Coerce 'lme4-model' to 'lmerModLmerTest':
    res <- as(model, "lmerModLmerTest")
    # Set relevant slots of the new model object:
    res@sigma <- sigma(model)
    res@vcov_beta <- as.matrix(vcov(model))
    varpar_opt <- unname(c(res@theta, res@sigma))
    # Compute Hessian:
    h <- numDeriv::hessian(func=devfun_vp, x=varpar_opt, devfun=devfun,
                           reml=is_reml)
    # Eigen decompose the Hessian:
    eig_h <- eigen(h, symmetric=TRUE)
    evals <- eig_h$values
    neg <- evals < -tol
    pos <- evals > tol
    zero <- evals > -tol & evals < tol
    if(sum(neg) > 0) { # negative eigenvalues
      eval_chr <- if(sum(neg) > 1) "eigenvalues" else "eigenvalue"
      evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
      warning(sprintf("Model failed to converge with %d negative %s: %s",
                      sum(neg), eval_chr, evals_num), call.=FALSE)
    }
    # Note: we warn about negative AND zero eigenvalues:
    if(sum(zero) > 0) { # some eigenvalues are zero
      eval_chr <- if(sum(zero) > 1) "eigenvalues" else "eigenvalue"
      evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
      warning(sprintf("Model may not have converged with %d %s close to zero: %s",
                      sum(zero), eval_chr, evals_num))
    }
    # Compute vcov(varpar):
    pos <- eig_h$values > tol
    q <- sum(pos)
    # Using the Moore-Penrose generalized inverse for h:
    h_inv <- with(eig_h, {
      vectors[, pos, drop=FALSE] %*% diag(1/values[pos], nrow=q) %*%
        t(vectors[, pos, drop=FALSE]) })
    res@vcov_varpar <- 2 * h_inv # vcov(varpar)
    # Compute Jacobian of cov(beta) for each varpar and save in list:
    Jac <- numDeriv::jacobian(func=get_covbeta, x=varpar_opt, devfun=devfun)
    res@Jac_list <- lapply(1:ncol(Jac), function(i)
      array(Jac[, i], dim=rep(length(res@beta), 2))) # k-list of jacobian matrices
    res
  }
  
  ##############################################
  ######## as_lmerModLmerTest()
  ##############################################
  #' Coerce lmerMod Objects to lmerModLmerTest
  #'
  #' Coercing an lme4::lmer model-object (of class 'lmerMod') to a model-object
  #' of class 'lmerModLmerTest' involves computing the covariance
  #' matrix of the variance parameters and the gradient (Jacobian) of cov(beta)
  #' with respect to the variance parameters.
  #'
  #' @param model and lmer model-object (of class 'lmerMod') -- the result of a
  #' call to \code{lme4::lmer()}
  #' @param tol tolerance for determining of eigenvalues are negative, zero or
  #' positive
  #'
  #' @return an object of class \code{'lmerModLmerTest'} which sets the following
  #' slots:
  #' \item{vcov_varpar}{the asymptotic covariance matrix of the variance parameters
  #' (theta, sigma).}
  #' \item{Jac_list}{list of Jacobian matrices; gradients of vcov(beta) with
  #' respect to the variance parameters.}
  #' \item{vcov_beta}{the asymptotic covariance matrix of the fixed-effect
  #' regression parameters (beta; vcov(beta)).}
  #' \item{sigma}{the residual standard deviation.}
  #'
  #' @seealso the class definition in \code{\link{lmerModLmerTest}}) and
  #' \code{\link{lmer}}
  #'
  #' @importFrom numDeriv hessian jacobian
  #' @importFrom stats vcov update
  #' @importFrom lme4 getME
  #'
  #' @author Rune Haubo B. Christensen
  #' @export
  #'
  #' @examples
  #' m <- lme4::lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
  #' bm <- as_lmerModLmerTest(m)
  #' slotNames(bm)
  #'
  as_lmerModLmerTest <- function(model, tol=1e-8) {
    if(!inherits(model, "lmerMod"))
      stop("model not of class 'lmerMod': cannot coerce to class 'lmerModLmerTest")
    # Get devfun:
    # 'Tricks' to ensure that we get the data to construct devfun even when
    # lmerTest is not attached or called inside a function:
    mc <- getCall(model)
    args <- c(as.list(mc), devFunOnly=TRUE)
    # if 'control' is not set we suppress potential message about rank deficient X
    # when evaulating devfun:
    if(!"control" %in% names(as.list(mc)))
      args$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
    Call <- as.call(c(list(quote(lme4::lmer)), args[-1]))
    ff <- environment(formula(model))
    pf <- parent.frame()  ## save parent frame in case we need it
    sf <- sys.frames()[[1]]
    ff2 <- environment(model)
    devfun <- tryCatch(eval(Call, envir=pf),
                       error=function(e) {
                         tryCatch(eval(Call, envir=ff),
                                  error=function(e) {
                                    tryCatch(eval(Call, envir=ff2),
                                             error=function(e) {
                                               tryCatch(eval(Call, envir=sf),
                                                        error=function(e) {
                                                          "error" })})})})
    if((is.character(devfun) && devfun == "error") ||
       !is.function(devfun) || names(formals(devfun)[1]) != "theta")
      stop("Unable to extract deviance function from model fit")
    as_lmerModLT(model, devfun, tol=tol)
  }
  
  
  ##############################################
  ######## devfun_vp()
  ##############################################
  #' Compute Deviance of an LMM as a Function of Variance Parameters
  #'
  #' This function is used for extracting the asymptotic variance-covariance matrix
  #'   of the variance parameters.
  #'
  #' @param varpar variance parameters; \code{varpar = c(theta, sigma)}.
  #' @param devfun deviance function as a function of theta only.
  #' @param reml if \code{TRUE} the REML deviance is computed;
  #'   if \code{FALSE}, the ML deviance is computed.
  #'
  #' @return the REML or ML deviance.
  #' @author Rune Haubo B. Christensen
  #' @keywords internal
  devfun_vp <- function(varpar, devfun, reml) {
    nvarpar <- length(varpar)
    sigma2 <- varpar[nvarpar]^2
    theta <- varpar[-nvarpar]
    df_envir <- environment(devfun)
    devfun(theta) # Evaluate deviance function at varpar
    n <- nrow(df_envir$pp$V)
    # Compute deviance for ML:
    dev <- df_envir$pp$ldL2() + (df_envir$resp$wrss() + df_envir$pp$sqrL(1))/sigma2 +
      n * log(2 * pi * sigma2)
    if(!reml) return(dev)
    # Adjust if REML is used:
    RX <- df_envir$pp$RX() # X'V^{-1}X ~ crossprod(RX^{-1}) = cov(beta)^{-1} / sigma^2
    dev + 2*c(determinant(RX)$modulus) - ncol(RX) * log(2 * pi * sigma2)
  }
  
  ##############################################
  ######## get_covbeta()
  ##############################################
  #' Compute cov(beta) as a Function of varpar of an LMM
  #'
  #' At the optimum cov(beta) is available as vcov(lmer-model). This function
  #' computes cov(beta) at non (RE)ML estimates of \code{varpar}.
  #'
  #' @inheritParams devfun_vp
  #'
  #' @return cov(beta) at supplied varpar values.
  #' @author Rune Haubo B. Christensen
  #' @keywords internal
  get_covbeta <- function(varpar, devfun) {
    nvarpar <- length(varpar)
    sigma <- varpar[nvarpar] # residual std.dev.
    theta <- varpar[-nvarpar] # ranef var-par
    devfun(theta) # evaluate REML or ML deviance 'criterion'
    df_envir <- environment(devfun) # extract model environment
    sigma^2 * tcrossprod(df_envir$pp$RXi()) # vcov(beta)
  }
  
  ##############################################
  ######## update.lmerModLmerTest()
  ##############################################
  ## We need our own update method for lmerModLmerTest objects because relying on
  ## lme4::update.merMod will sometimes return an object of class "lmerMod"
  ## instead of "lmerModLmerTest". This for instance happened if formula was a
  ## character vector, e.g.:
  ##   form <- "Informed.liking ~ Product+Information+
  ##   (1|Consumer) + (1|Product:Consumer) + (1|Information:Consumer)"
  ##   m <- lmer(form, data=ham)
  ##   class(m)                        # "lmerModLmerTest"
  ##   class(update(m, ~.- Product))   # "lmerMod"
  ## in versions < 3.0-1.9002.
  ##
  #' @importFrom stats getCall update.formula
  #' @export
  #' @keywords internal
  update.lmerModLmerTest <- function(object, formula., ..., evaluate = TRUE) {
    if(is.null(call <- getCall(object)))
      stop("object should contain a 'call' component")
    extras <- match.call(expand.dots = FALSE)$...
    if(!missing(formula.))
      call$formula <- update.formula(formula(object), formula.)
    if(length(extras) > 0) {
      existing <- !is.na(match(names(extras), names(call)))
      for(a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if(any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }
    if(evaluate) {
      ff <- environment(formula(object))
      pf <- parent.frame()
      sf <- sys.frames()[[1]]
      res <- tryCatch(eval(call, envir = ff), error = function(e) {
        tryCatch(eval(call, envir = sf), error = function(e) {
          eval(call, pf)
        })
      })
      # 'res' may be "lmerMod" instead of "lmerModLmerTest" in which case we
      # coerce to "lmerModLmerTest":
      if(inherits(res, "lmerMod") && !inherits(res, "lmerModLmerTest"))
        as_lmerModLmerTest(res) else res
    } else call
  }
  
  
  
  
  mc <- match.call()
  mc[["fitted.model"]] <- NULL   # remove added arg
  orig_call <- mc
  mc[[1L]] <- quote(lme4::lmer)
  #model <- eval.parent(mc)   # don't actually fit the model!
  model <-as(fitted.model, "lmerModLmerTest")   # convert fitted.model
  if (devFunOnly) 
    return(model)
  args <- as.list(mc)
  args$devFunOnly <- TRUE
  if (!"control" %in% names(as.list(mc))) 
    args$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
  Call <- as.call(c(list(quote(lme4::lmer)), args[-1]))
  devfun <- eval.parent(Call)
  res <- as_lmerModLT(model, devfun)
  res@call <- orig_call
  return(res)
}

pare.random.slopes                        <- function(fitted.lmer, prop.var.threshold=.01, consider.dropping.intercepts=FALSE) {
  # this function takes a fitted LMER (fitted.lmer) and, for each random factor, identifies the
  # number of random intercept/slopes n that account for less variance than a specified
  # threshold (prop.var.threshold; .01 by default) using the rePCA() function.
  # then, it identifies the n random intercept/slopes that allegedly account for the lowest
  # variance. it then creates and returns a list containing two objects:
  #    $keep is a string containing all supported random slopes, which can be dropped into
  #       a model's formula in place of all random effects. if consider.dropping.intercepts=TRUE,
  #       intercepts are eligible to be dropped during this process; otherwise, they will always be retained.
  #    $remove is a string containing the random effects that were dropped (if any).
  # 
  # this function will never remove a random factor entirely as it only determines
  # whether or not to drop an intercept/slope based on its *relative* share of
  # its random factor's variance. thus, it is not guaranteed to return a formula
  # that will converge.
  # 
  # NOTE: **This function does not work as intended if the model includes random correlations!!**
  #       It has also not been tested for optimizers other than 'bobyqa'.
  
  # check to make sure there are no random correlations
  variance.groups <- as.data.frame(summary(fitted.lmer)[["varcor"]])[["grp"]]
  
  rand.corrs.absent = length(variance.groups) == length(unique(variance.groups))
  
  if( rand.corrs.absent ) {
    
    intercept_string       <- "(Intercept)"
    
    # is.na() selects vars from the table as opposed to correlations
    fitted.lmer.varCor     <- subset(as.data.frame(summary(fitted.lmer)[["varcor"]]), is.na(var2))
    
    rePCA_summary          <- summary(rePCA(fitted.lmer))
    factorNames            <- names(rePCA_summary)
    
    allRandomTermsToKeep   <- ""
    allRandomTermsToRemove <- ""
    
    for( nextFactorName in factorNames ) {
      factor.num <- which(factorNames == nextFactorName)
      
      if( length(factor.num) == 1 ) {
        nextPropVars <- rePCA_summary[[nextFactorName]][["importance"]][2,]
        
        if(!all(is.na(nextPropVars))) {
          numSlopesToRemove <- sum(nextPropVars < prop.var.threshold)
        } else {   # triggers if there is a lone random intercept
          numSlopesToRemove <- 0
        }
        
        rowNumsOfFactor <- which(grepl(paste0(nextFactorName, "\\.*"), fitted.lmer.varCor[["grp"]]))
        factorVarCors <- fitted.lmer.varCor[rowNumsOfFactor,]
        nameSlopesAll <- factorVarCors[,"var1"]
        
        if( numSlopesToRemove > 0 ) {
          rowNumsOfSlopesToRemove <- order(factorVarCors[,"vcov"])[seq(numSlopesToRemove)]
          nameSlopesToRemove <- factorVarCors[rowNumsOfSlopesToRemove,"var1"]
          
          # if required to retain random intercept, drop that (if present) from list of slopes to remove
          if(!consider.dropping.intercepts & (intercept_string %in% nameSlopesToRemove)) {
            nameSlopesToRemove <- nameSlopesToRemove[nameSlopesToRemove != intercept_string]
            numSlopesToRemove <- numSlopesToRemove - 1
          }
          
          nameSlopesToKeep <- nameSlopesAll[!nameSlopesAll %in% nameSlopesToRemove]
        } else {
          nameSlopesToKeep <- nameSlopesAll
          nameSlopesToRemove <- NULL
        }
        
        # do this before adding 0 if needed
        corrTerm <- ifelse(length(nameSlopesToKeep) > 1, '||', '|')   # use single bar if only one random term left
        
        if( intercept_string %in% nameSlopesToKeep ) {
          nameSlopesToKeep[nameSlopesToKeep == intercept_string] <- "1"   # replace "(Intercept)" if present
        } else {
          nameSlopesToKeep <- c("0", nameSlopesToKeep)   # exclude intercept otherwise
        }
        
        randomTermToKeep <- paste0('(', paste(nameSlopesToKeep, collapse=' + '), ' ', corrTerm, ' ', nextFactorName, ')')
        
        if( nchar(allRandomTermsToKeep) == 0 ) {
          allRandomTermsToKeep <- randomTermToKeep
        } else {
          allRandomTermsToKeep <- paste0(allRandomTermsToKeep, " + ", randomTermToKeep)
        }
        
        
        if( numSlopesToRemove > 0 ) {
          if( intercept_string %in% nameSlopesToRemove ) {
            nameSlopesToRemove[nameSlopesToRemove == intercept_string] <- "1"
          }
          
          randomTermToRemove <- paste0('(', paste(nameSlopesToRemove, collapse=' + '), ' || ', nextFactorName, ')')
          
          if( nchar(allRandomTermsToRemove) == 0 ) {
            allRandomTermsToRemove <- randomTermToRemove
          } else {
            allRandomTermsToRemove <- paste0(allRandomTermsToRemove, " + ", randomTermToRemove)
          }
        }
      } else {
        stop("All random factors must be uniquely named!")
      }
    }
    
    # return
    return(list(keep=allRandomTermsToKeep, remove=allRandomTermsToRemove))
    
  } else {
    
    # identify correlated random slopes
    rand.comps.with.corr.slopes <- names(summary(as.factor(variance.groups))[summary(as.factor(variance.groups)) != 1])
    corr.slopes <- subset(as.data.frame(summary(fitted.lmer)[["varcor"]])[variance.groups %in% rand.comps.with.corr.slopes,], is.na(var2))
    rand.comps.with.corr.slopes <- unique(corr.slopes[["grp"]])   # same as above, but now in the same order as the function call
    
    # construct output to identify random slopes to user (and align spaces for easy reading)
    corr.slopes.output <- vector(mode="character", length=length(rand.comps.with.corr.slopes))
    max.uniqueComp.name.length <- max(sapply(rand.comps.with.corr.slopes, nchar))
    for( nextUniqueCompIndex in 1:length(rand.comps.with.corr.slopes) ) {
      nextUniqueComp <- rand.comps.with.corr.slopes[nextUniqueCompIndex]
      extra.spaces <- paste(rep(" ",max.uniqueComp.name.length - nchar(nextUniqueComp)), collapse="")
      corr.slopes.output[nextUniqueCompIndex] <- paste("   ", nextUniqueComp, ":   ", extra.spaces, paste(subset(corr.slopes, grp==nextUniqueComp)[["var1"]], collapse=", "), sep="")
    }
    
    output.text <- paste(
      "",
      "It appears this (g)lmer object has correlations between random slopes.",
      "The method used by this function to identify the slopes accounting for",
      "the least variance does not work properly when those correlations are",
      "present. (Note that this includes correlations between different levels",
      "of a single variable when it is specified as a nominal variable in the",
      "random effects; such variables must be explicitly recoded as numeric",
      "variables.)",
      "",
      "The correlated random slopes identified are shown below:",
      paste(corr.slopes.output, collapse="\n"),
      sep="\n")
    
    stop(output.text)
  }
}

numericize.factors                        <- function(model.components, component.to.numericize, num.col.suffix=".num") {
  # This function takes a list of model components and an associated data frame, then does the following:
  #    (1) All nominal predictors in the random effects terms (or whichever term is specified by component.to.numericize)
  #        will be converted to variables explicitly coded as numeric. (This is formally equivalent to leaving them as
  #        factors, but it makes it easier to remove random effects in pieces to facilitate model selection; furthermore,
  #        the || notation to remove random correlations does not work properly unless all random effects are coded this way.)
  #    (2) Columns representing the numeric predictors are added to the data frame, which can then be used to fit models with
  #        the given (fixed and) random effects.
  
  # model.components must be a list object that contains the following:
  #    $data:     the df containing data to be analyzed
  #    $formula:  must contain the following components as sub-lists:
  #       $dv:    the dependent variable
  #       $fixef: the fixed effects
  #       $ranef: must contain 1+ factors as sub-lists; each one must consist of the (intercept &) slopes to vary by that random factor
  #
  # component.to.numericize ("fixef" or "ranef") tells the model which component to convert to numeric
  #
  # num.col.suffix (by default, '.num') is a string that will be appended to the
  # end of the name of every factor column that is numericized using this function
  
  # kludge, but it works!
  model.components[["formula"]][["fixef"]] <- list(fixef=model.components[["formula"]][["fixef"]])
  
  variable.names.in.effects <- NULL
  
  # clumsily extract variable names from effects components
  for(next.effect.component in names(model.components[["formula"]][[component.to.numericize]])) {
    # list of characters that can be part of variable names
    variable.name.chars <- c(LETTERS, letters, 0:9, "_", ".")
    # split all characters, identify non-variable name characters, and use them to identify character groupings
    next.effect.component.split <- strsplit(model.components[["formula"]][[component.to.numericize]][[next.effect.component]], "")[[1]]
    variable.char.idx <- which(next.effect.component.split %in% variable.name.chars)
    variable.char.groupings <- cumsum(c(1, abs(variable.char.idx[-length(variable.char.idx)] - variable.char.idx[-1]) > 1))
    intercept.implicitly.included <- 1   # assume yes until otherwise specified
    # for each grouping, reconstruct the term and add it to the list if (1) it's not a number, and (2) it's unique (hasn't already been added)
    for(next.variable.char.grouping in unique(variable.char.groupings)) {
      next.variable.name <- paste(next.effect.component.split[variable.char.idx[variable.char.groupings==next.variable.char.grouping]], collapse="")
      # don't include 0 or 1 (absence/presence of intercept)
      # (suppress expected warnings due to NA coercion of factor names)
      if(suppressWarnings(is.na(as.numeric(next.variable.name)))) {
        if(!next.variable.name %in% variable.names.in.effects) {
          # confirm that the variable name is a column in the dataset
          if(next.variable.name %in% names(model.components[["data"]])) {
            variable.names.in.effects <- c(variable.names.in.effects, next.variable.name)
          } else {
            stop(paste0("Variable name '", next.variable.name, "' detected in effects component '", next.effect.component, "', but no matching column name found in data frame."))
          }
        }
      } else {   # it's a number, so (regardless of whether it's a 0 or 1) the intercept is not implicitly included
        intercept.implicitly.included <- 0
      }
    }
    
    if(intercept.implicitly.included) {   # if intercept is implicitly included, make it explicit
      # add "1 +" to the front (I think this should work,)
      model.components[["formula"]][[component.to.numericize]][[next.effect.component]] <- paste0("1 + ", model.components[["formula"]][[component.to.numericize]][[next.effect.component]])
    }
  }
  
  # sort variable names in order of length to make next step slightly more efficient
  if(length(variable.names.in.effects) > 0) {   # don't do anything else if the only random terms specified are numeric (thanks to Idan Blank for catching this)
    variable.name.lengths <- unlist(lapply(variable.names.in.effects, nchar))
    variable.names.in.effects <- variable.names.in.effects[order(variable.name.lengths)]
    # check to make sure that no variable name is contained entirely within any other variable name (which would make substitution, as accomplished here, impossible)
    if(length(variable.names.in.effects)>1) {
      for(next.variable.name.idx in 1:(length(variable.names.in.effects)-1)) {
        match.idx <- grepl(variable.names.in.effects[next.variable.name.idx], variable.names.in.effects[(next.variable.name.idx+1):length(variable.names.in.effects)])
        if(any(match.idx)) {
          stop(paste0("Variable name '", variable.names.in.effects[next.variable.name.idx], "' is contained within variable names '", paste(variable.names.in.effects[which(match.idx)], collapse="', '"), "', making substitution difficult."))
        }
      }
    }
    
    # initialize data frame of variable substitutions
    # (this object is not actually returned, but could be packaged with output for increased transparency)
    variable.substitutions <- as.data.frame(cbind(variable_name_old=variable.names.in.effects, variable_name_new=variable.names.in.effects))
    
    for(next.variable.name in variable.names.in.effects) {
      # check to make sure it's a factor; if not, don't change name
      if(is.factor(model.components[["data"]][[next.variable.name]])) {
        
        # make sure no columns already exist in df that match the naming pattern
        existing.matching.cols <- grepl(paste0("^", next.variable.name, num.col.suffix, "[0-9]+$"), names(model.components[["data"]]))
        if(any(existing.matching.cols)) {
          stop(paste0("data frame already contains columns with names matching pattern: ", paste(names(model.components[["data"]])[existing.matching.cols], collapse=", ")))
        }
        
        # create numeric columns
        next.variable.as.numeric <- model.matrix(~ model.components[["data"]][[next.variable.name]])
        # remove intercept column
        next.variable.as.numeric <- as.data.frame(next.variable.as.numeric[,2:ncol(next.variable.as.numeric)])
        # rename columns: if only 1 column, append ".num" (or num.col.suffix); if 2+ columns, append column #s as well (.num1, .num2, etc.)
        if(ncol(next.variable.as.numeric)>1) {
          names(next.variable.as.numeric) <- paste0(next.variable.name, num.col.suffix, 1:ncol(next.variable.as.numeric))
        } else {
          names(next.variable.as.numeric) <- paste0(next.variable.name, num.col.suffix                                  )
        }
        # attach columns to df
        model.components[["data"]] <- cbind(model.components[["data"]], next.variable.as.numeric)
        # update substitution list
        formula.substitution <- names(next.variable.as.numeric)
        if(length(formula.substitution) > 1) {
          formula.substitution <- paste0("(", paste(formula.substitution, collapse=" + "), ")")
        }
        variable.substitutions[(variable.substitutions[["variable_name_old"]]==next.variable.name),"variable_name_new"] <- formula.substitution
        
        # substitute numeric variables for factor names in effects
        for(next.effect.component in names(model.components[["formula"]][[component.to.numericize]])) {
          model.components[["formula"]][[component.to.numericize]][[next.effect.component]] <- gsub(next.variable.name, formula.substitution, model.components[["formula"]][[component.to.numericize]][[next.effect.component]])
        }
      }
    }
  }
  
  model.components[["formula"]][["fixef"]] <- model.components[["formula"]][["fixef"]][["fixef"]]
  
  return(model.components)
}

fit.lmer                                  <- function(model.components, return.all.models=FALSE, model.optimizer='bobyqa', verbosity=0, prop.var.threshold=.01, consider.dropping.intercepts=FALSE, num.col.suffix=".num", estimate.denom.df=TRUE) {
  # This function takes a list of model components and a data frame, then adopts the following approach to fit a converging model:
  #    (1) First, the model is fit with a maximal random effects structure (all random intercepts and slopes, and correlations between random effects).
  #    (2) If that does not converge, correlations between random effects are removed.
  #    (3) If that does not converge, all random effects that account for < X% (by default, X=1%) of the variance of their respective random factors are simultaneously removed.
  #        This step is repeated until the model converges (or until it would not change the random effects structure further, at which point the function gives up).
  
  # model.components: A list object that contains the following:
  #    $data:     the df containing data to be analyzed
  #    $formula:  must contain the following components as sub-lists:
  #       $dv:    the dependent variable
  #       $fixef: the fixed effects
  #       $ranef: must contain 1+ factors as sub-lists; each one must consist of the (intercept &) slopes to vary by that random factor
  #
  # return.all.models: If FALSE, this function will return a (g)lmerMod/lmerModLmerTest object. If TRUE, it will return a list object that contains the following:
  #    $components: Same as model.components above, but with several changes:
  #       (1) All nominal predictors in each index in $components$formula$ranef (the random effects) will be replaced by numeric predictors
  #       (2) Columns representing the numeric predictors will be added to $components$data (the data frame)
  #    $models: A list object that contains one sub-list for each model fit when this function was called: $full, $noRanCorrs, $reduced.1, $reduced.2, ... $reduced.n.
  #             (Note that $models$noRanCorrs may not exist if $models$full did not contain any random correlations.)
  #             Models are added to the list in the order in which they were fit; thus, if the returned object is called 'model.lmers',
  #             model.lmers$models[[length(model.lmers$models)]] will always return information about the last model fit (which most likely converged).
  #             Each of these objects will have three sub-objects:
  #       $formula:     The formula of the model that was fit
  #       $model:       The (g)lmerMod object corresponding to the model that was fit
  #       $convergence: A convergence code indicating whether the model converged (1) or not (0)
  #
  # model.optimizer: Set to 'bobyqa' by default.
  #
  # verbosity: Set to '0' (no output) by default. Does not affect model fitting.
  #
  # prop.var.threshold: The minimum proportion of variance that a random intercept/slope must account for (relative to its random factor)
  #                     to be retained during model selection; .01 (1%) by default. See pare.random.slopes() for more information.
  #
  # consider.dropping.intercepts: Whether or not random intercepts are eligible to be dropped during model selection.
  #                               See pare.random.slopes() for more information.
  #
  # num.col.suffix: A string that will be appended to the end of the name of every factor column that is numericized; '.num' by default.
  #
  # estimate.denom.df: Whether or not to use lmerTest to estimate denominator df of a converging model when appropriate. TRUE by default.
  
  
  # replace nominal variables in random effects with numeric equivalents
  model.components <- numericize.factors(model.components, 'ranef', num.col.suffix)
  
  # initialize all.models with model components & data (with numeric columns added)
  all.models <- list(components=model.components)
  
  # construct model formula from components
  formula.full   <- paste0(model.components[["formula"]][["dv"]], " ~ ", model.components[["formula"]][["fixef"]])
  for(next.ranef.component in names(model.components[["formula"]][["ranef"]])) {
    formula.full <- paste0(formula.full, " + (", model.components[["formula"]][["ranef"]][[next.ranef.component]], " | ", next.ranef.component, ")")
  }
  
  # determine whether to run lmer() or glmer()
  dv.is_binary <- all(model.components[["data"]][[model.components[["formula"]][["dv"]]]] %in% c(0,1))
  
  # run full model within try/catch statement in case lme4 throws an error before attempting to fit the model
  attempt.to.fit.model <- function(model.formula.as.text) {
    tryCatch({
      # use lme4 (g)lmer functions (denominator df are only computed on converging model to save time)
      if(dv.is_binary) {
        next.lmer <- lme4::glmer(formula(model.formula.as.text), data=model.components[["data"]], family=binomial, control=glmerControl(optimizer=model.optimizer), verbose=verbosity)
      } else {
        next.lmer <-  lme4::lmer(formula(model.formula.as.text), data=model.components[["data"]],                  control= lmerControl(optimizer=model.optimizer), verbose=verbosity)
      }
      return(next.lmer)
    }, error=function(cond) {   # may throw an error (e.g., if the number of random effects > the number of observations)
      return(NULL)
    })
  }
  
  lmer.full <- attempt.to.fit.model(formula.full)
  
  # if model fitted, check convergence; no model fitted is equivalent to non-convergence (and is likely due to overparamerization)
  if( typeof(lmer.full)=="S4" ) {
    most.recent.lmer <- lmer.full
    lmer.convergence <- length(most.recent.lmer@optinfo[["conv"]][["lme4"]])==0
  } else {
    lmer.convergence <- 0
  }
  
  all.models[["models"]][["full"]] <- list(formula=formula.full, model=lmer.full, convergence=lmer.convergence)
  
  if( !lmer.convergence ) {
    # if not, remove random correlations:
    
    model.components.noRanCorrs <- model.components
    formula.noRanCorrs <- paste0(model.components[["formula"]][["dv"]], " ~ ", model.components[["formula"]][["fixef"]])
    
    # replace "|" with "||" for all random factors with 2+ terms
    for(next.ranef.component in names(model.components.noRanCorrs[["formula"]][["ranef"]])) {
      # list of characters that can be part of variable names
      variable.name.chars <- c(LETTERS, letters, 0:9, "_", ".")
      # split all characters, identify non-variable name characters, and use them to identify character groupings
      next.ranef.component.split <- strsplit(model.components.noRanCorrs[["formula"]][["ranef"]][[next.ranef.component]], "")[[1]]
      variable.char.idx <- which(next.ranef.component.split %in% variable.name.chars)
      variable.char.groupings <- cumsum(c(1, abs(variable.char.idx[-length(variable.char.idx)] - variable.char.idx[-1]) > 1))
      next.ranef.component.formula <- paste0("(", model.components[["formula"]][["ranef"]][[next.ranef.component]], " | ", next.ranef.component, ")")
      if(length(unique(variable.char.groupings)) > 1) {   # more than one random term
        # NOTE: this will fail to identify an implicitly specified intercept, but numericize.ranef.factors should have made that explicit
        next.ranef.component.formula <- gsub("\\|", "\\||", next.ranef.component.formula)
      }
      formula.noRanCorrs <- paste0(formula.noRanCorrs, " + ", next.ranef.component.formula)
    }
    
    # run lmer unless formula hasn't changed (e.g., because model only has random intercepts -- although in that case, this function will not help...)
    if(formula.full!=formula.noRanCorrs) {
      lmer.noRanCorrs <- attempt.to.fit.model(formula.noRanCorrs)
    } else {
      lmer.noRanCorrs <- lmer.full
    }
    
    if( typeof(lmer.noRanCorrs)=="S4" ) {
      most.recent.lmer <- lmer.noRanCorrs
      lmer.convergence <- length(most.recent.lmer@optinfo[["conv"]][["lme4"]])==0
      all.models[["models"]][["noRanCorrs"]] <- list(formula=formula.noRanCorrs, model=lmer.noRanCorrs, convergence=lmer.convergence)
    } else {
      stop(cat("lme4 was unable to attempt to fit a model -- an error was thrown prior to fitting.",
                  "\nThis is probably because the model is mis-specified and/or because there are more",
                  "\nrandom effects terms than there are observations (check the console window for any errors)."))
    }
  }
  
  if( !lmer.convergence ) {
    num.reduced.models.run <- 1
    
    loop.aborted <- 0
    
    while(lmer.convergence==0 & loop.aborted==0) {   # remove recursively (though it is rare for >1 iteration to be needed)
      # if not, remove all random slopes that account for < 1% of the variance of their respective random factors
      
      random.slopes.to.keep.and.remove <- pare.random.slopes(most.recent.lmer, prop.var.threshold, consider.dropping.intercepts=consider.dropping.intercepts)
      
      if(nchar(random.slopes.to.keep.and.remove[["remove"]]) > 0) {   # only proceed if 1+ random slopes are being removed
        # remove those random slopes and re-fit the model
        formula.noRanCorrs.reduced <- paste0(model.components[["formula"]][["dv"]], " ~ ", model.components[["formula"]][["fixef"]], " + ", random.slopes.to.keep.and.remove[["keep"]])
        
        lmer.noRanCorrs.reduced <- attempt.to.fit.model(formula.noRanCorrs.reduced)
        most.recent.lmer <- lmer.noRanCorrs.reduced
        lmer.convergence <- length(most.recent.lmer@optinfo[["conv"]][["lme4"]])==0
        
        all.models[["models"]][[paste0("reduced.", num.reduced.models.run)]] <- list(formula=formula.noRanCorrs.reduced, model=lmer.noRanCorrs.reduced, convergence=lmer.convergence)
        
        num.reduced.models.run <- num.reduced.models.run + 1
      } else {   # no change in model formula between successive iterations; abort
        loop.aborted <- 1
      }
    }
    
    if(loop.aborted) {
      rePCA.results <- summary(rePCA(most.recent.lmer))
      rand.factors  <- names(rePCA.results)
      
      rePCA.output <- ""
      
      for(next.rand.factor.idx in 1:length(rand.factors)) {
        if(next.rand.factor.idx>1) {
          next.rePCA.output <- "\n   "
        } else {
          next.rePCA.output <- "   "
        }
        
        next.rePCA.output <- paste0(next.rePCA.output, rand.factors[next.rand.factor.idx], ": ")
        
        # conversion to data frame prevents an error when there is only one col
        next.rePCA.results <- as.data.frame(rePCA.results[[rand.factors[next.rand.factor.idx]]][["importance"]][c("Standard deviation","Proportion of Variance"),])
        for(next.rePCA.factor.idx in 1:ncol(next.rePCA.results)) {
          if(next.rePCA.factor.idx==1) {
            next.rePCA.factor.output <- ""
          } else {
            next.rePCA.factor.output <- ", "
          }
          next.rePCA.factor.output <- paste0(next.rePCA.factor.output, signif(next.rePCA.results["Standard deviation",next.rePCA.factor.idx],4), " (", round(100*next.rePCA.results["Proportion of Variance",next.rePCA.factor.idx],1), "%)")
          next.rePCA.output <- paste0(next.rePCA.output, next.rePCA.factor.output)
        }
        
        rePCA.output <- paste0(rePCA.output, next.rePCA.output)
      }
      
      intercepts.output.text <- ifelse(consider.dropping.intercepts, "", "random intercepts and/or ")
      
      cat(paste0("Model did not converge despite recursive removal of random effects components accounting for < ", 100*prop.var.threshold, "% of variance of their respective random factors.",
                 "\nConsider increasing that threshold and/or removing ", intercepts.output.text, "random factors altogether.",
                 "\nThe most recent pared-down effects structure that did not converge:",
                 "\n   ", random.slopes.to.keep.and.remove[["keep"]],
                 "\nsummary(rePCA(that.lmer)) is shown here. For each random factor, the SD accounted for by each intercept/slope (ranked) is shown, and the corresponding % of variance is shown in parentheses:",
                 "\n", rePCA.output
      ))
      # return list of models no matter what
      return(all.models)
    }
  }
  
  # compute denominator df using lmerTest functions if (1) specified via arguments, and (2) model has a continuous DV
  if(estimate.denom.df & !dv.is_binary) {
    all.models[["models"]][[length(all.models[["models"]])]][["model"]] <- compute.lmerTest.dfs(all.models[["models"]][[length(all.models[["models"]])]][["model"]], formula(all.models[["models"]][[length(all.models[["models"]])]][["model"]]), data=model.components[["data"]], control=lmerControl(optimizer=model.optimizer), verbose=verbosity)
  }
  
  if(return.all.models) {
    return(all.models)   # return list of all model objects, convergence statuses, etc.
  } else {
    return(all.models[["models"]][[length(all.models[["models"]])]][["model"]])   # return converging model only
  }
}


# ## prepare data and fit model ##
# 
# # import sample data
# sample.data <- read.delim(paste0(file_dir, "gloves_Exp1a_usable_data.csv"), sep=",")   # change to directory where data file was downloaded
# 
# # things that must be done manually before fitting the model:
# #    centering numeric variables
# #    ordering factor levels and setting contrast weightings for each level
# #    determining the correct initial model specification
# sample.data[["item_type"]]                    <- factor(sample.data[["item_type"]], c("control", "covid"))
# contrasts(sample.data[["item_type"]])         <- -contr.sum(2)/2
# sample.data[["interference_type"]]            <- factor(sample.data[["interference_type"]], c("noise", "cough"))
# contrasts(sample.data[["interference_type"]]) <- -contr.sum(2)/2
# 
# # create object, model.components, that stores:
# #    $data:     the df containing data to be analyzed
# #    $formula:  must contain the following components as sub-lists:
# #       $dv:    the dependent variable
# #       $fixef: the fixed effects
# #       $ranef: must contain 1+ factors as sub-lists; each one must consist of the (intercept &) slopes to vary by the random factor corresponding to its name
# model.components                            <- list()
# model.components$data                       <- sample.data
# model.components$formula$dv                 <- "match_target_related"
# model.components$formula$fixef              <- "item_type * interference_type"   # intercepts can be specified either implicitly or explicitly ("1 + ")
# # (better/more robust to use [[]] notation (vs. $ notation) to store info about
# # random factors, which are column names that can contain spaces)
# model.components$formula$ranef[["subj_id"]] <- "item_type * interference_type"   # make sure to change names of random factors when using your own data
# model.components$formula$ranef[["item_id"]] <- "item_type * interference_type"   # there should be one line for each random factor
# 
# # fit converging model using default parameters except that output should contain information about all models fit (not only the final one)
# model.lmers <- fit.lmer(model.components, return.all.models=TRUE)   # will likely take a few minutes
# 
# 
# ## inspect output ##
# 
# # random effects terms used in full-model formula after columns were numericized and any implicit intercepts were made explicit
# model.lmers$components$formula$ranef
# 
# # vector of names, each corresponding to a model that was fit
# names(model.lmers$models)
# 
# # information about first (full) model fit
# # can inspect other models by replacing "full" with "noRanCorrs", "reduced.1", etc.
# model.lmers$models$full$formula       # conversely, model.lmers$models[[1]]$formula, because it was the first model fit
# model.lmers$models$full$model         # full model; however, if NULL, that means the model was never fit because lme4 threw an error (due to, e.g., having more random effects than observations)
# model.lmers$models$full$convergence   # convergence status (TRUE=converged, FALSE=non-converged)
# 
# # if a model converged, it's always the last model fit
# converging.model.info <- model.lmers$models[[length(model.lmers$models)]]
# converging.model.info$convergence   # first, confirm convergence! if FALSE, there should have been text printed to the console window explaining where the process got tripped up
# converging.model.info$formula       # formula of converging model
# converging.model.info$model         # converging model (with denominator df estimated, if DV is continuous and estimate.denom.df==TRUE)
# 
# # one-line method to extract converging model (equivalent to running fit.lmer() with return.all.models=FALSE)
# converging.model <- model.lmers$models[[length(model.lmers$models)]]$model
