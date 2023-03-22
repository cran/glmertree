utils::globalVariables(c(".tree", ".ranef", ".weights", ".cluster"))

betamertree <- function(formula, data, family = NULL, weights = NULL,
                        cluster = NULL, ranefstart = NULL, offset = NULL,
                        REML = TRUE, joint = TRUE, abstol = 0.001, maxit = 100L,  
                        dfsplit = TRUE, verbose = FALSE, plot = FALSE, 
                        glmmTMB.control = glmmTMB::glmmTMBControl(), ...)
{
  ## remember call
  cl <- match.call()
  
  ## TODO: check if data is complete
  
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
  ff <- Formula::as.Formula(formula) ## full formula
  tf <- formula(ff, lhs = 1L, rhs = c(1L, 3L)) ## tree formula
  tf <- Formula::as.Formula(tf)
  pf <- formula(tf, lhs = 0L, rhs = 2L) ## partition formula
  tf <- formula(tf, lhs = 1L, rhs = 1L)
  
  if (length(attr(ff, "rhs")[[2L]]) == 1L) {
    rf <- (. ~ (1 | id))[[3L]]
    rf[[2L]][[3L]] <- attr(ff, "rhs")[[2L]]
    attr(ff, "rhs")[[2L]] <- rf
  }
  #if (joint) {
    rf <- formula(ff, lhs = 1L, rhs = 1L)
    rf <- update(rf, . ~ .tree / .)
    rf <- formula(Formula::as.Formula(rf, formula(ff, lhs = 0L, rhs = 2L)),
                  lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
  #} else {
  #  rf <- formula(ff, lhs = 1L, rhs = 2L)
  #}
  
  ## TODO: process data
  # if (data_has_missings) {
  #   old_N <- nrow(data)
  #   data <- data[complete.cases(data[ , all_vars]), ]
  #   warning(paste0("New sample size is N = ", nrow(data), " (old sample size was N = ", old_N, ")."))
  # }
  
  ## initialization
  iteration <- 0L
  data$.ranef <- #if (is.null(ranefstart)) {
    rep(0, times = dim(data)[1L])
  #} else if (ranefstart && length(ranefstart) == 1) {
    ### generate ranefstart from lme null model: 
    #predict(glmer(formula(ff, lhs = 1L, rhs = 2L),
    #              data = data, weights = .weights, nAGQ = nAGQ,
    #              offset = offset, family = family, control = glmer.control),
    #        newdata = data, type = "link")
  #} else {
  #  ranefstart  
  #}
  continue <- TRUE
  oldloglik <- c(-Inf, - Inf) # last element is the oldest
  
  ## iterate between glmer and glmtree estimation
  while (continue) {
    
    iteration <- iteration + 1L
    
    ## if offset was specified, add it to .ranef:
    if (!is.null(cl$offset)) {
      data$.ranef <- data$.ranef + offset
    }
    
    ## fit tree
    #if (is.null(q_cluster)) {
      tree <- betareg::betatree(tf, partition = pf, data = data,  
                       offset = .ranef, link = "logit", link.phi = "log",
                       weights = .weights, dfsplit = dfsplit, ...)
    #} else {
    #  tree <- glmtree(tf, data = data, family = family, offset = .ranef, 
    #                  weights = .weights, cluster = .cluster, 
    #                  dfsplit = dfsplit, ...)
    #}
    
    if (plot) plot(tree, main = paste("Iteration", iteration))
    
    data$.tree <- #if (joint) {
      factor(predict(tree, newdata = data, type = "node"))
    #} else {
    #  predict(tree, newdata = data, type = "link")
    #  ## note that these predictions already include the offset
    #}
    
    ## fit glmmTMB
    #if (joint) {
      ## estimate full glmmTMB model but force all coefficients from the
      ## .tree (and the overall intercept) to zero for the prediction
      if (length(tree) == 1L) {
        ## If tree of depth 1 was grown, (g)lmer model should not include interactions:
        rf.alt <- formula(ff, lhs = 1L, rhs = 1L)
        rf.alt <- formula(Formula::as.Formula(rf.alt, formula(ff, lhs = 0L, rhs = 2L)),
                          lhs = 1L, rhs = c(1L, 2L), collapse = TRUE)
        if (is.null(offset)) {
          glmmTMB <- glmmTMB::glmmTMB(rf.alt, data = data, family = glmmTMB::beta_family, 
                             weights = .weights, control = glmmTMB.control,
                             REML = REML)
        } else {
          glmmTMB <- glmmTMB::glmmTMB(rf.alt, data = data, family = glmmTMB::beta_family, 
                             weights = .weights, offset = offset, 
                             control = glmmTMB.control, REML = REML)          
        }
      } else {
        if (is.null(offset)) {
          glmmTMB <- glmmTMB::glmmTMB(rf, data = data, family = glmmTMB::beta_family, 
                             weights = .weights, control = glmmTMB.control,
                             REML = REML)
        } else {
          glmmTMB <- glmmTMB::glmmTMB(rf, data = data, family = glmmTMB::beta_family, 
                             weights = .weights, offset = offset, 
                             control = glmmTMB.control, REML = REML)
        }
      }
    
    
    ## generate random-effects predictions
    data$.ranef <- predict(glmmTMB, re.form = NULL, type = "link") - ## mixed-eff predictions
        predict(glmmTMB, re.form = NA, type = "link") ## minus fixed-eff predictions
      
      
    #} else {
    #  glme <- glmer(rf, data = data, family = family, offset = .tree, 
    #                nAGQ = nAGQ, weights = .weights, control = glmer.control)
    #  data$.ranef <- predict(glme, newdata = data, type = "link")
    #  ## note that because newdata is specified, predict.merMod will not include offset in predictions
    #}
  
    ## iteration information
    newloglik <- logLik(glmmTMB)    
    continue <- (abs(newloglik - oldloglik[1]) > abstol) & (iteration < maxit) 
    if (continue & (abs(newloglik - oldloglik[2]) < abstol)) {
      if (newloglik > oldloglik[1]) continue <- FALSE
    }
    oldloglik[2] <- oldloglik[1]
    oldloglik[1] <- newloglik
    if (verbose) print(newloglik)
  }
  
  ## collect results
  result <- list(
    formula = formula,
    call = cl,
    tree = tree,
    glmmTMB = glmmTMB,
    ranef = ranef(glmmTMB), 
    varcorr = VarCorr(glmmTMB),
    variance = attr(VarCorr(glmmTMB),"sc")^2, 
    data = data,
    nobs = nrow(data),
    loglik = as.numeric(newloglik),
    df = attr(newloglik, "df"),
    dfsplit = dfsplit,
    iterations = iteration, 
    maxit = maxit,
    #ranefstart = ranefstart, 
    abstol = abstol,
    mob.control = list(...),
    glmmTMB.control = glmmTMB.control,
    joint = joint
  )
  class(result) <- "betamertree"
  return(result)
}



fixef.betamertree <- coef.betamertree <- 
  function(object, which = "tree", drop = FALSE, ...) {
  if (which == "tree") { 
    coefs <- coef(object$tree, drop = FALSE)
    coefs <- coefs[ , -which(colnames(coefs) == "(phi)_(Intercept)"), drop = FALSE]
    if (object$joint) { ## overwrite tree coefs with those (g)lmer:
      glmmTMB_fixef <- fixef(object[["glmmTMB"]])$cond
      if (nrow(coefs) > 1L) {
        ## Add the intercept to intercepts of all other nodes:
        glmmTMB_fixef[paste0(".tree", rownames(coefs)[-1L])] <- 
          glmmTMB_fixef[paste0(".tree", rownames(coefs)[-1L])] + glmmTMB_fixef[1L]
        ## Change name of intercept to first terminal node:
        names(glmmTMB_fixef)[1L] <- paste0(".tree", rownames(coefs)[1L])
        local_coefs <- glmmTMB_fixef[grepl(".tree", names(glmmTMB_fixef))]
        for (i in rownames(coefs)) {
          coef_names <- paste0(".tree", i)
          if (ncol(coefs) > 1L) {
            coef_names <- c(coef_names, paste0(coef_names, ":", 
                                               colnames(coefs)[-1]))
          }
          coefs[i, ] <- local_coefs[coef_names]
        }
      } else {
        coefs[1L, ] <- glmmTMB_fixef[colnames(coefs)]
      }
    }
  } else if (which == "global") {
    coefs <- fixef(object[["glmmTMB"]])
    coefs <- coefs[-which(names(coefs) == "(Intercept)")]
    coefs <- coefs[!grepl(".tree", names(coefs))]
  }
  if (drop) return(drop(coefs)) else return(coefs)
}



VarCorr.betamertree <- function(x, ...) {
  VarCorr(x[["glmmTMB"]], ...)
}



plot.betamertree <- function(x, which = "tree", type = "extended", 
  tp_args = list(), drop_terminal = TRUE, terminal_panel = NULL, ...) {
  
  ## Keep it simple, just plot the fitted betatree
  if (type == "extended") {
    plot(x$tree, type = type)
  } else if (type == "simple") {
    ## overwrite tree coefs with lmer coefs if joint estimation was used
    if (x$joint) {
      node_ids <- unique(x$tree$fitted[["(fitted)"]])
      gt_node <- as.list(x$tree$node)
      coefs <- coef(x, which = "tree")
      for (i in node_ids) {
        gt_node[[i]]$info$coefficients <- coefs[as.character(i), ]
      }
      x$tree$node <- as.partynode(gt_node)
    }
    FUN <- simple_terminal_func
    plot(x$tree, drop_terminal = drop_terminal,
         terminal_panel = node_terminal_glmertree, 
         tp_args = list(FUN = FUN, align = "right"), ...)
  } 
}




residuals.betamertree <- function(object, ...) {    
  if (object$joint) {
    resids <- residuals(object$glmmTMB, ...)
  } else {
    stop("To obtain residuals please fit the model with the default betamertree(..., joint = TRUE)")
  } 
  return(resids)
}

ranef.betamertree <- function(object, ...) {
  ranef(object$glmmTMB, ...)
}

logLik.betamertree <- function(object, dfsplit = NULL, ...) {
  if (is.null(dfsplit)) dfsplit <- object$dfsplit
  dfsplit <- as.integer(dfsplit) * 
    (length(object$tree) - length(nodeids(object$tree, terminal = TRUE)))
  structure(object$loglik, df = object$df + dfsplit, nobs = object$nobs, 
            class = "logLik")
}


model.frame.betamertree <- function(formula, ...) {
  mf <- model.frame(formula$tree, ...)
  mf[["(offset)"]] <- mf[["(weights)"]] <- NULL
  dc <- attr(attr(mf, "terms"), "dataClasses")
  dc <- dc[!(names(dc) %in% c("(offset)", "(weights)"))]
  attr(attr(mf, "terms"), "dataClasses") <- dc
  return(mf)
}


terms.betamertree <- function(x, ...) {
  terms(x$tree, ...)
}


as.party.betamertree <- function(obj, ...) {
  obj$tree
}


print.betamertree <- function(x, title = "Beta mixed-effects regression tree", 
                            ...) {
  print(x$tree, title = title, ...)
  cat("\nRandom effects:\n")
  print(x$ranef)
  cat("\nFixed effects (from glmmTMB model):\n")
  print(fixef(x$glmmTMB)[-c(1L, grep(".tree", names(fixef(x$glmmTMB))))])
  invisible(x)
}



predict.betamertree <- function(object, newdata = NULL, type = "response", 
                                re.form = NULL, ...) { 
  if (is.null(newdata)) {
    newdata <- object$data
  }
  if (type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    #if (object$joint) {
      newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
      newdata$.tree <- factor(newdata$.tree,
                              levels = levels(object$data$.tree))
      newdata$.weights <- 1L
      predict(object$glmmTMB, newdata = newdata, type = type, re.form = re.form, 
              ...)
    #} else {
    #  newdata$.ranef <- predict(object$glmer, newdata = newdata, type = "link", 
    #                            re.form = re.form, ...)
    #  predict(object$tree, newdata = newdata, type = type)
    #}
  }
}
