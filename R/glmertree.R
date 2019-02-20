utils::globalVariables(c(".tree", ".ranef", ".weights", ".cluster"))

lmertree <- function(formula, data, weights = NULL, cluster = NULL,
                     ranefstart = NULL, offset = NULL, joint = TRUE,
                     abstol = 0.001, maxit = 100L, dfsplit = TRUE, 
                     verbose = FALSE, plot = FALSE, 
                     lmer.control = lmerControl(), ...)
{
  ## remember call
  cl <- match.call()
  
  ## check if data is complet, print warning if not: 
  if (nrow(data) != sum(stats::complete.cases(data))) {
    warning("data contains missing values, note that listwise deletion will be employed.", immediate. = TRUE) 
  }
  
  ## process offset:
  q_offset <- substitute(offset)
  if (!is.null(q_offset)) {
    offset <- eval(q_offset, data)
  }
  
  ## process cluster:
  q_cluster <- substitute(cluster)
  if (!is.null(q_cluster)) {
    data$.cluster <- eval(q_cluster, data)
    if (length(eval(q_cluster, data)) != nrow(data)) {
      warning("Variable lengths differ for 'cluster' and 'data'.", immediate. = TRUE)
    }
  } 
  if (!is.null(q_cluster) && 
      !inherits(data$.cluster, c("numeric", "character", "factor", "integer"))) {
    warning("Argument 'cluster' should specify an object of class numeric, factor or character, or NULL.", immediate. = TRUE)
  } 
  
  # process weights:
  q_weights <- substitute(weights)
  if (!is.null(q_weights)) {
    data$.weights <- eval(q_weights, data)
    weights <- eval(q_weights, data)
    ## Note: Assigning weights and offset as variables outside data prevents lmer from yielding an error
  } else {
    data$.weights <- rep(1L, times = nrow(data))
  }
  
  
  ## formula processing (full, tree, random)
  ff <- Formula::as.Formula(formula)
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 3L))
  if(length(attr(ff, "rhs")[[2L]]) == 1L) {
    rf <- (. ~ (1 | id))[[3L]]
    rf[[2L]][[3L]] <- attr(ff, "rhs")[[2L]]
    attr(ff, "rhs")[[2L]] <- rf
  }
  if(joint) {
    rf <- formula(ff, lhs = 1L, rhs = 1L)
    rf <- update(rf, . ~ .tree / .)
    rf <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
                  lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
  } else {
    rf <- formula(ff, lhs = 1L, rhs = 2L)
  }
  
  ## initialization
  iteration <- 0L
  data$.ranef <- if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  } else if (ranefstart && length(ranefstart) == 1L) {
    ## generate ranefstart from lme null model: 
    predict(lmer(formula(ff, lhs = 1L, rhs = 2L),
                 data = data, weights = .weights, 
                 offset = offset, control = lmer.control),
            newdata = data)
  } else {
    ranefstart  
  }
  continue <- TRUE
  oldloglik <- -Inf
  
  ## iterate between lmer and lmtree estimation
  while (continue) {
    
    ## if offset was specified, add it to .ranef:
    if (!is.null(offset)) {
      data$.ranef <- data$.ranef + offset
    }
    
    iteration <- iteration + 1L
    
    ## fit tree
    if (is.null(q_cluster)) {
      tree <- lmtree(tf, data = data, offset = .ranef, weights = .weights, 
                     dfsplit = dfsplit, ...)
    } else {
      tree <- lmtree(tf, data = data, offset = .ranef, weights = .weights, 
                     cluster = .cluster, dfsplit = dfsplit, ...)
    }
    
    if (plot) plot(tree)
    
    data$.tree <- if (joint) {
      factor(predict(tree, newdata = data, type = "node"))
    } else {
      predict(tree, newdata = data, type = "response") 
      ## note that these predictions already include offset
    }
    
    ## fit lmer
    if (joint) {
      ## estimate full lmer model but force all coefficients from the
      ## .tree (and the overall intercept) to zero for the prediction
      if (length(tree) == 1L) {
        ## If tree of depth 1 was grown, (g)lmer model should not include interactions:
        rf.alt <- formula(ff, lhs = 1L, rhs = 1L)
        rf.alt <- formula(Formula::as.Formula(rf.alt, formula(ff, lhs = 0L, rhs = 2L)),
                          lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
        if (is.null(offset)) {
          lme <- lmer(rf.alt, data = data, weights = .weights, 
                      control = lmer.control)          
        } else {
          lme <- lmer(rf.alt, data = data, weights = .weights, 
                      offset = offset, control = lmer.control)
        }
      } else { # a tree was grown
        if (is.null(offset)) {
          lme <- lmer(rf, data = data, weights = .weights,
                      control = lmer.control)
        } else {
          lme <- lmer(rf, data = data, weights = .weights, 
                      offset = offset, control = lmer.control)          
        }
      }
      b <- structure(lme@beta, .Names = names(fixef(lme)))
      b[substr(names(b), 1L, 5L) %in% c("(Inte", ".tree")] <- 0
      data$.ranef <- suppressWarnings(suppressMessages(predict(lme, newdata = data, newparams = list(beta = b))))
    } else {
      ## estimate only a partial lmer model using the .tree fitted
      ## values as an offset
      lme <- lmer(rf, data = data, offset = .tree, weights = .weights, control = lmer.control)
      data$.ranef <- predict(lme, newdata = data)
      ## note that because newdata is specified, predict.merMod will not include offset in predictions
    }
    
    ## iteration information
    newloglik <- logLik(lme)    
    continue <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if (verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    lmer = lme,
    ranef = ranef(lme), 
    varcorr = VarCorr(lme),
    variance = attr(VarCorr(lme),"sc")^2, 
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    ranefstart = ranefstart, 
    abstol = abstol,
    mob.control = list(...),
    lmer.control = lmer.control,
    joint = joint
  )
  class(result) <- "lmertree"
  return(result)
}


glmertree <- function(formula, data, family = "binomial", weights = NULL,
                      cluster = NULL, ranefstart = NULL, offset = NULL,
                      joint = TRUE, abstol = 0.001, maxit = 100L,  
                      dfsplit = TRUE, verbose = FALSE, plot = FALSE, 
                      glmer.control = glmerControl(), ...)
{
  ## remember call
  cl <- match.call()
  
  ## check if data is complet, print warning if not: 
  if (nrow(data) != sum(stats::complete.cases(data))) {
    warning("data contains missing values, note that listwise deletion will be employed.", immediate. = TRUE) 
  }
  
  ## process offset:
  q_offset <- substitute(offset)
  if (!is.null(q_offset)) {
    offset <- eval(q_offset, data)
  }
  
  ## process cluster:
  q_cluster <- substitute(cluster)
  if (!is.null(q_cluster)) {
    data$.cluster <- eval(q_cluster, data)
    if (length(eval(q_cluster, data)) != nrow(data)) {
      warning("Variable lengths differ for 'cluster' and 'data'.", immediate. = TRUE)
    }
  } 
  if (!is.null(q_cluster) && 
      !inherits(data$.cluster, c("numeric", "character", "factor", "integer"))) {
    warning("Argument 'cluster' should specify an object of class numeric, factor or character, or NULL.", immediate. = TRUE)
  } 
  
  # process weights:
  q_weights <- substitute(weights)
  if (!is.null(q_weights)) {
    data$.weights <- eval(q_weights, data)
    weights <- eval(q_weights, data)
    ## Note: Assigning weights and offset as variables outside data prevents lmer from yielding an error.
  } else {
    data$.weights <- rep(1L, times = nrow(data))
  }
  
  ## formula processing (full, tree, random)
  ff <- Formula::as.Formula(formula)
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 3L))
  if (length(attr(ff, "rhs")[[2L]]) == 1L) {
    rf <- (. ~ (1 | id))[[3L]]
    rf[[2L]][[3L]] <- attr(ff, "rhs")[[2L]]
    attr(ff, "rhs")[[2L]] <- rf
  }
  if (joint) {
    rf <- formula(ff, lhs = 1L, rhs = 1L)
    rf <- update(rf, . ~ .tree / .)
    rf <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
                  lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
  } else {
    rf <- formula(ff, lhs = 1L, rhs = 2L)
  }
  
  ## initialization
  iteration <- 0L
  data$.ranef <- if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  } else if (ranefstart && length(ranefstart) == 1) {
    ## generate ranefstart from lme null model: 
    predict(glmer(formula(ff, lhs = 1L, rhs = 2L),
                  data = data, weights = .weights,
                  offset = offset, family = family, control = glmer.control),
            newdata = data, type = "link")
  } else {
    ranefstart  
  }
  continue <- TRUE
  oldloglik <- -Inf
  
  ## iterate between glmer and glmtree estimation
  while (continue) {
    
    iteration <- iteration + 1L
    
    ## if offset was specified, add it to .ranef:
    if (!is.null(offset)) {
      data$.ranef <- data$.ranef + offset
    }
    
    ## fit tree
    if (is.null(q_cluster)) {
      tree <- glmtree(tf, data = data, family = family, offset = .ranef, 
                      weights = .weights, dfsplit = dfsplit, ...)
    } else {
      tree <- glmtree(tf, data = data, family = family, offset = .ranef, 
                      weights = .weights, cluster = .cluster, 
                      dfsplit = dfsplit, ...)
    }
    
    if (plot) plot(tree)
    
    data$.tree <- if (joint) {
      factor(predict(tree, newdata = data, type = "node"))
    } else {
      predict(tree, newdata = data, type = "link")
      ## note that these predictions already include the offset
    }
    
    ## fit glmer
    if (joint) {
      ## estimate full glmer model but force all coefficients from the
      ## .tree (and the overall intercept) to zero for the prediction
      if (length(tree) == 1L) {
        ## If tree of depth 1 was grown, (g)lmer model should not include interactions:
        rf.alt <- formula(ff, lhs = 1L, rhs = 1L)
        rf.alt <- formula(Formula::as.Formula(rf.alt, formula(ff, lhs = 0L, rhs = 2L)),
                          lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
        if (is.null(offset)) {
          glme <- glmer(rf.alt, data = data, family = family, 
                        weights = .weights, control = glmer.control)
        } else {
          glme <- glmer(rf.alt, data = data, family = family, 
                        weights = .weights, offset = offset, 
                        control = glmer.control)          
        }
      } else {
        if (is.null(offset)) {
          glme <- glmer(rf, data = data, family = family, weights = .weights,
                        control = glmer.control)
        } else {
          glme <- glmer(rf, data = data, family = family, weights = .weights,
                        offset = offset, control = glmer.control)
        }
      }       
      b <- structure(glme@beta, .Names = names(fixef(glme)))
      b[substr(names(b), 1L, 5L) %in% c("(Inte", ".tree")] <- 0
      data$.ranef <- suppressWarnings(suppressMessages(
        predict(glme, newdata = data, type = "link", newparams = list(beta = b))))
    } else {
      glme <- glmer(rf, data = data, family = family, offset = .tree, 
                    weights = .weights, control = glmer.control)
      data$.ranef <- predict(glme, newdata = data, type = "link")
      ## note that because newdata is specified, predict.merMod will not include offset in predictions
    }
    ## iteration information
    newloglik <- logLik(glme)    
    continue <- (newloglik - oldloglik > abstol) & (iteration < maxit) 
    oldloglik <- newloglik
    if (verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    glmer = glme,
    ranef = ranef(glme), 
    varcorr = VarCorr(glme),
    variance = attr(VarCorr(glme),"sc")^2, 
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    ranefstart = ranefstart, 
    abstol = abstol,
    mob.control = list(...),
    glmer.control = glmer.control,
    joint = joint
  )
  class(result) <- "glmertree"
  return(result)
}

coef.lmertree <- coef.glmertree <- function(object, ...) {
  coef(object$tree, ...)
}

VarCorr.glmertree <- function(object, ...) {
  VarCorr(object$glmer)
}

VarCorr.lmertree <- function(object, ...) {
  VarCorr(object$lmer)
}

plot.lmertree <- plot.glmertree <- function(x, which = "all", ask = TRUE, ...) {    
  if (which != "ranef") {
    plot(x$tree, ...)
  }
  if (which != "tree") {
    if (which == "all" && ask == TRUE) {
      orig_devAsk <- devAskNewPage()
      devAskNewPage(ask = TRUE)
    }
    if (requireNamespace("lattice")) {
      print(lattice::dotplot(ranef(x$lmer, condVar = TRUE), main = TRUE))
    }
    if (which == "all" && ask == TRUE) {grDevices::devAskNewPage(ask = orig_devAsk)}
  }
}

plot.glmertree <- function(x, plotranef = FALSE, which = "all", ask = TRUE, ...) {
  if (which != "ranef") {
    plot(x$tree, ...)
  }
  if (which != "tree") {
    if (which == "all" && ask == TRUE) {
      orig_devAsk <- devAskNewPage()
      devAskNewPage(ask = TRUE)
    }
    if (requireNamespace("lattice")) {
      print(lattice::dotplot(ranef(x$glmer, condVar = TRUE), main = TRUE))
    }
    if (which == "all" && ask == TRUE) {grDevices::devAskNewPage(ask = orig_devAsk)}
  }
}

residuals.lmertree <- resid.lmertree <- function(object, type = NULL, scaled = FALSE, ...) {    
  if(object$joint) {
    if(is.null(type)) {
      resids <- residuals(object$lmer, scaled = scaled)
    } else {
      resids <- residuals(object$lmer, type = type, scaled = scaled)      
    }
  } else {
    resids <- object$data[, all.vars(object$formula[[2]])] - predict(
      object, newdata = NULL)
    if(scaled) {
      resids <- scale(resids, center = FALSE, scale = TRUE)
    }
  }
  return(resids)
}

residuals.glmertree <- resid.glmertree <- function(object, type = NULL, scaled = FALSE, ...) {    
  if(object$joint) {
    if(is.null(type)) {
      resids <- residuals(object$glmer, scaled = scaled)
    } else {
      resids <- residuals(object$glmer, type = type, scaled = scaled)      
    }
  } else {
    stop("To obtain residuals please fit the model with the default glmertree(..., joint = TRUE)")
  } 
  return(resids)
}

ranef.lmertree <- ranef.glmertree <- function(object, ...) {
  object$ranef
}

logLik.lmertree <- logLik.glmertree <- function(object, dfsplit = NULL, ...)
{
  if(is.null(dfsplit)) dfsplit <- object$dfsplit
  dfsplit <- as.integer(dfsplit) * (length(object$tree) - length(nodeids(object$tree, terminal = TRUE)))
  structure(object$loglik, df = object$df + dfsplit, nobs = object$nobs, class = "logLik")
}

model.frame.lmertree <- model.frame.glmertree <- function(formula, ...) {
  mf <- model.frame(formula$tree, ...)
  mf[["(offset)"]] <- mf[["(weights)"]] <- NULL
  dc <- attr(attr(mf, "terms"), "dataClasses")
  dc <- dc[!(names(dc) %in% c("(offset)", "(weights)"))]
  attr(attr(mf, "terms"), "dataClasses") <- dc
  return(mf)
}

terms.lmertree <- terms.glmertree <- function(x, ...) {
  terms(x$tree, ...)
}

as.party.lmertree <- as.party.glmertree <- function(obj, ...) {
  obj$tree
}

print.lmertree <- function(x, title = "Linear mixed model tree", ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  if(x$joint & length(fixef(x$lmer)[-grep(".tree", names(fixef(x$lmer)))]) > 1L) {
    cat("\nLinear fixed effects (from lmer model):\n")
    print(fixef(x$lmer)[-c(1L, grep(".tree", names(fixef(x$lmer))))])
  }
  if(!x$joint & length(fixef(x$lmer)) > 1L) {
    cat("\nLinear fixed effects (from lmer model):\n")
    print(fixef(x$lmer)[-1L])  
  }  
  invisible(x)
}

print.glmertree <- function(x, title = "Generalized linear mixed model tree", ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  if(x$joint & length(fixef(x$glmer)[-grep(".tree", names(fixef(x$glmer)))]) > 1L) {
    cat("\nLinear fixed effects (from glmer model):\n")
    print(fixef(x$glmer)[-c(1L, grep(".tree", names(fixef(x$glmer))))])
  }
  if(!x$joint & length(fixef(x$glmer)) > 1L) {
    cat("\nLinear fixed effects (from glmer model):\n")
    print(fixef(x$glmer)[-1L])  
  }  
  invisible(x)
}

predict.lmertree <- function(object, newdata = NULL, type = "response", 
                             re.form = NULL, ...) { 
  if(is.null(newdata)) {
    newdata <- object$data
  }
  if(type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    if(object$joint) {
      newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
      newdata$.tree <- factor(newdata$.tree)
      levels(newdata$.tree) <- levels(object$data$.tree)
      predict(object$lmer, newdata = newdata, type = type, re.form = re.form, ...)
    } else {
      newdata$.ranef <- predict(object$lmer, newdata = newdata, re.form = re.form, ...)
      predict(object$tree, newdata = newdata, type = type)
    }
  }
}

predict.glmertree <- function(object, newdata = NULL, type = "response", 
                              re.form = NULL, ...) { 
  if(is.null(newdata)) {
    newdata <- object$data
  }
  if(type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    if(object$joint) {
      newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
      newdata$.tree <- factor(newdata$.tree)
      levels(newdata$.tree) <- levels(object$data$.tree)
      predict(object$glmer, newdata = newdata, type = type, re.form = re.form, 
              ...)
    } else {
      newdata$.ranef <- predict(object$glmer, newdata = newdata, type = "link", 
                                re.form = re.form, ...)
      predict(object$tree, newdata = newdata, type = type)
    }
  }
}