utils::globalVariables(c(".tree", ".ranef", ".weights", ".cluster"))

lmertree <- function(formula, data, weights = NULL, cluster = NULL,
                     ranefstart = NULL, offset = NULL, joint = TRUE,
                     abstol = 0.001, maxit = 100L, dfsplit = TRUE, 
                     verbose = FALSE, plot = FALSE, REML = TRUE,
                     lmer.control = lmerControl(), ...)
{
  ## remember call
  cl <- match.call()
  
  ## check if data is complet, print warning if not: 
  if (nrow(data) != sum(stats::complete.cases(data))) {
    warning("'data' contains missing values, note that listwise deletion will be employed.", immediate. = TRUE) 
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
  } else if (ranefstart && length(ranefstart) == 1L) {
    ## generate ranefstart from lme null model: 
    predict(lmer(formula(ff, lhs = 1L, rhs = 2L),
                 data = data, weights = .weights, REML = REML,
                 offset = offset, control = lmer.control),
            newdata = data)
  } else {
    ranefstart  
  }
  continue <- TRUE
  oldloglik <- c(-Inf, - Inf) # last element is the oldest
  
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
    
    if (plot) plot(tree, main = paste("Iteration", iteration))
    
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
                      REML = REML, control = lmer.control)          
        } else {
          lme <- lmer(rf.alt, data = data, weights = .weights, 
                      REML = REML, offset = offset, control = lmer.control)
        }
      } else { # a tree was grown
        if (is.null(offset)) {
          lme <- lmer(rf, data = data, weights = .weights,
                      REML = REML, control = lmer.control)
        } else {
          lme <- lmer(rf, data = data, weights = .weights, 
                      REML = REML, offset = offset, control = lmer.control)          
        }
      }
      b <- structure(lme@beta, .Names = names(fixef(lme)))
      b[substr(names(b), 1L, 5L) %in% c("(Inte", ".tree")] <- 0
      data$.ranef <- suppressWarnings(suppressMessages(
        predict(lme, newdata = data, newparams = list(beta = b))))
    } else {
      ## estimate only a partial lmer model using the .tree fitted
      ## values as an offset
      lme <- lmer(rf, data = data, offset = .tree, weights = .weights, 
                  REML = REML, control = lmer.control)
      data$.ranef <- predict(lme, newdata = data)
      ## note that because newdata is specified, predict.merMod will 
      ## not include offset in predictions
    }
    
    ## iteration information
    newloglik <- logLik(lme)    
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
                      nAGQ = 1L, glmer.control = glmerControl(), ...)
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
                  data = data, weights = .weights, nAGQ = nAGQ,
                  offset = offset, family = family, control = glmer.control),
            newdata = data, type = "link")
  } else {
    ranefstart  
  }
  continue <- TRUE
  oldloglik <- c(-Inf, - Inf) # last element is the oldest
  
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
    
    if (plot) plot(tree, main = paste("Iteration", iteration))
    
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
          glme <- glmer(rf.alt, data = data, family = family, nAGQ = nAGQ,
                        weights = .weights, control = glmer.control)
        } else {
          glme <- glmer(rf.alt, data = data, family = family, nAGQ = nAGQ,
                        weights = .weights, offset = offset, 
                        control = glmer.control)          
        }
      } else {
        if (is.null(offset)) {
          glme <- glmer(rf, data = data, family = family, weights = .weights,
                        nAGQ = nAGQ, control = glmer.control)
        } else {
          glme <- glmer(rf, data = data, family = family, weights = .weights,
                        nAGQ = nAGQ, offset = offset, control = glmer.control)
        }
      }       
      b <- structure(glme@beta, .Names = names(fixef(glme)))
      b[substr(names(b), 1L, 5L) %in% c("(Inte", ".tree")] <- 0
      data$.ranef <- suppressWarnings(suppressMessages(
        predict(glme, newdata = data, type = "link", newparams = list(beta = b))))
    } else {
      glme <- glmer(rf, data = data, family = family, offset = .tree, 
                    nAGQ = nAGQ, weights = .weights, control = glmer.control)
      data$.ranef <- predict(glme, newdata = data, type = "link")
      ## note that because newdata is specified, predict.merMod will not include offset in predictions
    }
    ## iteration information
    newloglik <- logLik(glme)    
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



fixef.lmertree <- coef.lmertree <- fixef.glmertree <- 
  coef.glmertree <- function(object, which = "tree", 
                             drop = FALSE, ...) {
  merMod_type <- ifelse(inherits(object, "lmertree"), "lmer", "glmer")
  if (which == "tree") { 
    coefs <- coef(object$tree, drop = FALSE)
    if (object$joint) { ## overwrite tree coefs with those (g)lmer:
      lmer_fixef <- fixef(object[[merMod_type]])
      if (nrow(coefs) > 1L) {
        ## Add the intercept to intercepts of all other nodes:
        lmer_fixef[paste0(".tree", rownames(coefs)[-1L])] <- 
          lmer_fixef[paste0(".tree", rownames(coefs)[-1L])] + lmer_fixef[1L]
        ## Change name of intercept to first terminal node:
        names(lmer_fixef)[1L] <- paste0(".tree", rownames(coefs)[1L])
        local_coefs <- lmer_fixef[grepl(".tree", names(lmer_fixef))]
        for (i in rownames(coefs)) {
          coef_names <- paste0(".tree", i)
          if (ncol(coefs) > 1L) {
            coef_names <- c(coef_names, paste0(coef_names, ":", 
                                               colnames(coefs)[-1]))
          }
          coefs[i, ] <- local_coefs[coef_names]
        }
      } else {
        coefs[1, ] <- lmer_fixef[colnames(coefs)]
      }
    }
  } else if (which == "global") {
    coefs <- fixef(object[[merMod_type]])
    coefs <- coefs[-which(names(coefs) == "(Intercept)")]
    coefs <- coefs[!grepl(".tree", names(coefs))]
  }
  if (drop) return(drop(coefs)) else return(coefs)
}



VarCorr.lmertree <- VarCorr.glmertree <- function(object, ...) {
  merMod_type <- ifelse(inherits(object, "lmertree"), "lmer", "glmer")
  VarCorr(object[[merMod_type]], ...)
}


# plot.lmertree <- plot.glmertree <- function(x, which = "all", ask = TRUE, 
#                                             type = "extended", ...) {    
# 
#   merMod_type <- ifelse(inherits(x, "lmertree"), "lmer", "glmer")  
#   if (which != "ranef") {
#     if (type == "extended") {
#       plot(x$tree, ...)
#     } else if (type == "simple") {
#       plot(x$tree, terminal_panel = node_terminal, tp_args = list(
#         FUN = simple_terminal_func, align = "right"))
#     }
#   }
#   if (which != "tree") {
#     if (which == "all" && ask == TRUE) {
#       orig_devAsk <- devAskNewPage()
#       devAskNewPage(ask = TRUE)
#     }
#     if (requireNamespace("lattice")) {
#       print(lattice::dotplot(ranef(x[[merMod_type]], condVar = TRUE), 
#                              main = TRUE))
#     }
#     if (which == "all" && ask == TRUE) {
#       grDevices::devAskNewPage(ask = orig_devAsk)
#     }
#   }
# }


## Helper function for extracting standard errors from merMod objects:
get_merMod_SEs <- function(object, which = "tree", global_intercept = TRUE) {
  ##
  ## object: should be of class "(g)lmertree"
  ## 
  ## which: should be "tree" (then a matrix in the shape of coef(object) will
  ##    be returned), or "global" (then a vector with SEs will be returned). 
  ##
  ## global_intercept: specifies whether(default) treatment contrast
  ##    coding was used in estimating the fixed effects in the (g)lmer 
  ##    model. If TRUE, standard errors for the intercepts of the 
  ##    non-first terminal nodes will be corrected for these intercepts 
  ##    being a sum of two coefficients.
  ##
  merMod_type <- ifelse(inherits(object, "lmertree"), "lmer", "glmer")
  ## Prepare matrix for returning SEs:
  SEs <- coef(object, which = which)
  SEs[is.numeric(SEs)] <- NA
  ## Get covariance matrix of parameter estimates:
  vc <- as.matrix(lme4::vcov.merMod(object[[merMod_type]]))
  colnames(vc)[1] <- rownames(vc)[1] <- paste0(".tree", rownames(SEs)[1])
  lmer_SEs <- sqrt(diag(vc))
  if (which == "tree") {
    lmer_SEs <- lmer_SEs[grepl(".tree", names(lmer_SEs))]
    for (i in rownames(SEs)) {
      coef_names <- paste0(".tree", i)
      if (ncol(SEs) > 1L) {
        coef_names <- c(coef_names, 
                        paste0(coef_names, ":", colnames(SEs)[-1]))
      }
      SEs[i, ] <- lmer_SEs[coef_names]
      if (global_intercept) {
        ## SE of the sum of two coefficients is given by 
        ##    sqrt( se(coef_1)^2 + se(coef_2)^2  + 2*cov(coef_1, coef_2) )
        if (i != rownames(SEs)[1L]) {
          SEs[i , "(Intercept)"] <- sqrt(
            SEs[1L, "(Intercept)"]^2 + SEs[i, "(Intercept)"]^2 + 
              2*vc[1L, rownames(vc) == paste0(".tree", i)])
        }
      }
    }
  } else SEs <- lmer_SEs[!grepl(".tree", names(lmer_SEs))]
  return(SEs)
}


plot.glmertree <- plot.lmertree <- function(x, which = "all", # c("tree", "ranef", "tree.coef", "global.coef", "all") 
                                              ask = TRUE,
                                              type = "extended", 
                                              observed = TRUE, 
                                              fitted = "combined", 
                                              tp_args = list(), 
                                              drop_terminal = TRUE, 
                                              terminal_panel = NULL, ...) {
  
  merMod_type <- ifelse(inherits(x, "lmertree"), "lmer", "glmer")
  ff <- Formula::as.Formula(x$formula)
  if (which %in% c("tree", "all", "tree.coef")) {
    if (attr(ff, "rhs")[[1]] == 1L && which != "tree.coef") {
      plot(x$tree, type = type)
    } else if (type == "simple") {
      FUN <- simple_terminal_func
      plot(x$tree, drop_terminal = drop_terminal,
           terminal_panel = node_terminal_glmertree, 
           tp_args = list(FUN = FUN, align = "right"))
    } else if (type == "extended") {
      mf <- model.frame(x)
      if (length(attr(ff, "rhs")[[2L]]) > 1L) {
        global_lm_vars <- all.vars(nobars(x$formula[[3]][[2]][[3]]))
        global_lm_vars <- global_lm_vars[!global_lm_vars %in% 
                                           names(ranef(x[[merMod_type]]))]
        if (length(global_lm_vars) > 0L) {
          warning("Global fixed effects were specified and/or estimated, but will not be plotted.")
        }
      } else {
        global_lm_vars <- character()
      }
      local_lm_vars <- all.vars(x$formula[[3]][[2]][[2]])
      lm_vars <- unique(c(global_lm_vars, local_lm_vars))
      if (which == "tree.coef") {
        coefs <- as.data.frame(fixef(x, which = "tree", drop = FALSE))
        long_coefs <- data.frame(utils::stack(coefs), node = rep(
          paste("node", rownames(coefs)), times = ncol(coefs)),
          stringsAsFactors = FALSE)
        if (x$joint) {
          local_SEs <- get_merMod_SEs(x, which = "tree")
          long_local_SEs <- utils::stack(as.data.frame(local_SEs))
          long_coefs$upper <- long_coefs$values + 1.96*long_local_SEs$values
          long_coefs$lower <- long_coefs$values - 1.96*long_local_SEs$values
        } else {
          long_coefs$upper <- long_coefs$lower <- long_coefs$values
        }
        panel <- function(x, y, lx, ux, subscripts, pch = 16) {
          lattice::panel.abline(h = unique(y), col = "grey")
          lattice::panel.arrows(lx[subscripts], y, ux[subscripts], y, 
                                col = 'black', length = 0, unit = "native", 
                                angle = 90, code = 3)
          lattice::panel.xyplot(x, y, pch = pch)
        }
        prepanel <- function(x, lx, ux, subscripts, ...) {
          list(xlim = range(x[subscripts], ux[subscripts], lx[subscripts], 
                            finite = TRUE))
        }
        print(lattice::dotplot(node ~ values | ind, data = long_coefs,
                      lx = long_coefs$lower, ux = long_coefs$upper,
                      prepanel = prepanel, panel = panel, as.table = TRUE,
                      scales = list(x = list(relation = "free")),
                      xlab = "Estimated coefficients", 
                      main = "Fixed effects from tree"))
      } else {
        if (length(local_lm_vars > 0L)) {
          if (is.null(tp_args$which)) {
            vars_to_plot <- 1L:length(local_lm_vars)
          } else {
            vars_to_plot <- tp_args$which
          }
          if (is.numeric(vars_to_plot)) vars_to_plot <- local_lm_vars[vars_to_plot]
          node_ids <- x$tree$fitted[["(fitted)"]]
          re.form <- if (fitted == "marginal") NA else NULL
          if (fitted == "marginal") {
            fitted_values <- matrix(nrow = length(node_ids), 
                                    ncol = length(vars_to_plot),
                                    dimnames = list(rownames(mf), vars_to_plot))
            for (varname in vars_to_plot) {
              newdata <- x$data
              remaining_lm_vars <- lm_vars[lm_vars != varname]
              if (length(lm_vars > 1L)) {
                for (i in remaining_lm_vars) {
                  if (class(x$data[, i]) == "numeric") {
                    ## set all values to mean value:
                    newdata[, i] <- mean(x$data[, i])
                  } else if (class(x$data[, i]) == "factor") {
                    ## set all values to most common class:
                    ux <- unique(x$data[, i])
                    newdata[, i] <- ux[which.max(tabulate(match(x$data[, i], ux)))]
                  }
                }
              }
              fitted_values[, varname] <- predict(x, newdata = newdata, 
                                                  re.form = re.form)
            }
          } else {
            fitted_values <- predict(x, newdata = x$data, re.form = re.form)
          }
          lt_node <- as.list(x$tree$node)
          for (i in unique(node_ids)) {
            if (fitted == "marginal") {
              lt_node[[i]]$info$fitted <- fitted_values[node_ids == i, ]  
            } else {
              lt_node[[i]]$info$fitted <- fitted_values[node_ids == i]          
            }
          }
          x$tree$node <- as.partynode(lt_node)
        }
        tp_args <- append(tp_args, values = list(
          ranef = ifelse(fitted == "marginal", "constant", "varying"), 
          fixef = ifelse(fitted == "marginal", "constant", "varying"), 
          fitmean = ifelse(fitted == "none", FALSE, TRUE),
          observed = observed))
        if (is.null(terminal_panel)) terminal_panel <- node_glmertree
        plot(x$tree, terminal_panel = terminal_panel, 
             drop_terminal = drop_terminal, tp_args = tp_args, ...)
      }
    }
  }
  if (which == "global.coef") {
    coefs <- fixef(x, which = "global")
    if (x$joint) {
      local_SEs <- get_merMod_SEs(x, which = "global")
      upper <- coefs + 1.96*local_SEs
      lower <- coefs - 1.96*local_SEs
    } else {
      upper <- lower <- coefs
    }
    panel <- function(x, y, lx, ux, subscripts, pch = 16) {
      lattice::panel.abline(h = unique(y), col = "grey")
      lattice::panel.arrows(lx[subscripts], y, ux[subscripts], y, 
                            col = 'black', length = 0, unit = "native", 
                            angle = 90, code = 3)
      lattice::panel.xyplot(x, y, pch = pch)
    }
    prepanel <- function(x, lx, ux, subscripts, ...) {
      list(xlim = range(x[subscripts], ux[subscripts], lx[subscripts], 
                        finite = TRUE))
    }
    print(lattice::dotplot(coefs, lx = lower, ux = upper,
                           prepanel = prepanel, panel = panel, as.table = TRUE,
                           scales = list(x = list(relation = "free")),
                           xlab = "Estimated coefficients", 
                           main = "Global fixed effects"))
  }
  if (which %in% c("all", "ranef")) {
    if (which == "all" && ask == TRUE) {
      orig_devAsk <- devAskNewPage()
      devAskNewPage(ask = TRUE)
    }
    if (requireNamespace("lattice")) {
      if (inherits(x, "lmertree")) {
        ## TODO: allow additional arguments to be passed to lattice dotplot?
        print(lattice::dotplot(ranef(x$lmer, condVar = TRUE), 
                               main = TRUE))
      } else {
        ## TODO: allow additional arguments to be passed to lattice dotplot?
        print(lattice::dotplot(ranef(x$glmer, condVar = TRUE), 
                               main = TRUE))
      }
    }
    if (which == "all" && ask == TRUE) {
      grDevices::devAskNewPage(ask = orig_devAsk)
    }
  }
}




residuals.lmertree <- resid.lmertree <- function(object, type = NULL, 
                                                 scaled = FALSE, ...) {    
  if (object$joint) {
    if (is.null(type)) {
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

residuals.glmertree <- resid.glmertree <- function(object, type = NULL, 
                                                   scaled = FALSE, ...) {    
  if (object$joint) {
    if (is.null(type)) {
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

logLik.lmertree <- logLik.glmertree <- function(object, dfsplit = NULL, ...) {
  if (is.null(dfsplit)) dfsplit <- object$dfsplit
  dfsplit <- as.integer(dfsplit) * 
    (length(object$tree) - length(nodeids(object$tree, terminal = TRUE)))
  structure(object$loglik, df = object$df + dfsplit, nobs = object$nobs, 
            class = "logLik")
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


print.glmertree <- function(x, title = "Generalized linear mixed model tree", 
                            ...) {
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
  if (is.null(newdata)) {
    newdata <- object$data
  }
  if (type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    if (object$joint) {
      newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
      newdata$.tree <- factor(newdata$.tree, 
                              levels = levels(object$data$.tree))
      predict(object$lmer, newdata = newdata, type = type, re.form = re.form, ...)
    } else {
      newdata$.ranef <- predict(object$lmer, newdata = newdata, re.form = re.form, ...)
      predict(object$tree, newdata = newdata, type = type)
    }
  }
}


predict.glmertree <- function(object, newdata = NULL, type = "response", 
                              re.form = NULL, ...) { 
  if (is.null(newdata)) {
    newdata <- object$data
  }
  if (type == "node") {
    predict(object$tree, newdata = newdata, type = "node")
  } else {
    if (object$joint) {
      newdata$.tree <- predict(object$tree, newdata = newdata, type = "node")
      newdata$.tree <- factor(newdata$.tree,
                              levels = levels(object$data$.tree))
      predict(object$glmer, newdata = newdata, type = type, re.form = re.form, 
              ...)
    } else {
      newdata$.ranef <- predict(object$glmer, newdata = newdata, type = "link", 
                                re.form = re.form, ...)
      predict(object$tree, newdata = newdata, type = type)
    }
  }
}
