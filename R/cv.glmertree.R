## CV function for GLMM trees
cv.lmertree <- cv.glmertree <- function(tree, newdata, reference = NULL, 
                                        omit.intercept = FALSE, ...) {
  
  merMod_type <- ifelse(inherits(tree, "lmertree"), "lmer", "glmer")
  
  ## Check if suitable value for reference has been specified, if so relevel
  if (!is.null(reference)) {
    if (!reference %in% c(0, as.numeric(levels(tree$data$.tree)))) {
      stop(paste0("Specified reference level should be a value in ", 
                   c(0, as.numeric(levels(tree$data$.tree))))) 
    } else {
      tree$data$.tree <- relevel(tree$data$.tree, ref = as.character(reference))
    }
  }
  
  ## print warning if no joint estimation was used
  if (!tree$joint) warning("Original lmertree was not fit using joint estimation.")
  
  ## Extract formula
  ff <- if (merMod_type == "lmer") {
    attr(tree$lmer, "call")$formula } else { attr(tree$glmer, "call")$formula 
  }
  
  ## Get node indicators for newdata
  newdata$.tree <- factor(predict(tree, newdata = newdata, type = "node"),
                          levels = levels(tree$data$.tree))
  
  ## Replace data
  omit_col_ids <- which(names(tree$tree$data) %in% c("(weights)", "(offset)", "(cluster)"))
  tree$tree$data <- newdata[names(tree$tree$data)[(1:ncol(tree$tree$data))[-omit_col_ids]]]
  
  tree$data <- newdata
  
  ## Replace (g)lmer
  if (omit.intercept) ff <- update(ff, . ~ . + 0)
  
  if (merMod_type == "lmer") {
    lmer.control <- if (!is.null(tree$call$lmer.control)) {
      lmerControl(tree$call$lmer.control)
    } else {
      lmerControl()
    }
    REML <- ifelse(!is.null(tree$call$REML), tree$call$REML, TRUE)
  } else {
    glmer.control <- if (!is.null(tree$call$glmer.control)) {
      glmerControl(tree$call$glmer.control) 
    } else { 
      glmerControl()
    }
    nAGQ <- ifelse(!is.null(tree$call$nAGQ), tree$call$nAGQ, 1L)
    family <- ifelse(!is.null(tree$call$family), tree$call$family, binomial)
  }
  
  if (merMod_type == "lmer") {
    tree$lmer <- lmer(ff, data = newdata, REML = REML, control = lmer.control)
  } else {
    tree$glmer <- glmer(ff, data = newdata, nAGQ = nAGQ, family = family,
                        control = glmer.control)
  }
  
  ## Replace stats in tree
  
  ## Extract node info
  tree_node <- as.list(tree$tree$node)
  nobs <- table(newdata$.tree)  
  coefs <- coef(tree, which = "tree")
  
  ## Replace node info
  for (i in length(tree_node):1L) {
    
    tree_node[[i]]$info$objfun <- NULL
    tree_node[[i]]$info$object <- NULL
    
    if (!is.null(tree_node[[i]]$kids)) { ## this is a splitting node
      
      tree_node[[i]]$info$test <- NULL
      tree_node[[i]]$info$p.value <- NULL
      tree_node[[i]]$info$coefficients <- NULL
      tree_node[[i]]$info$nobs <- sum(nobs[as.character(tree_node[[i]]$kids)]) 
      
    } else { ## this is a terminal node
      
      ## Replace
      tree_node[[i]]$info$nobs <- nobs[[as.character(i)]]
      
      if (is.null(reference) && (!omit.intercept)) {
        tree_node[[i]]$info$coefficients <- coefs[as.character(i), ]
      } else {
        ## TODO; This is slightly more involved and should be changed in (g)lmertree coef function
      }
      
    }
  }
  tree$tree$node <- as.partynode(tree_node)
  
  ## Assign class to result and return
  class(tree) <- if (inherits(tree, "lmertree")) {
    c("cv.lmertree", class(tree))
  } else {
    c("cv.glmertree", class(tree))
  }
  return(tree)
  
}



summary.cv.lmertree <- summary.cv.glmertree <- function(object, ...) {
  merMod_type <- ifelse(inherits(object, "lmertree"), "lmer", "glmer")
  summary(object[[merMod_type]], ...)
}
