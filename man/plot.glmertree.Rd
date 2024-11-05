\name{plot.glmertree}
\alias{plot.lmertree}
\alias{plot.glmertree}

\title{Plotting (Generalized) Linear Mixed Model Trees}

\description{
  \code{plot} method for \code{(g)lmertree} objects.
}

\usage{
\method{plot}{lmertree}(x, which = "all", nodesize_level = 1L, 
    cluster = NULL, ask = TRUE, type = "extended", 
    observed = TRUE, fitted = "combined", tp_args = list(), 
    drop_terminal = TRUE, terminal_panel = NULL, dotplot_args = list(), \dots)
\method{plot}{glmertree}(x, which = "all", nodesize_level = 1L, 
    cluster = NULL, ask = TRUE, type = "extended", 
    observed = TRUE, fitted = "combined", tp_args = list(), 
    drop_terminal = TRUE, terminal_panel = NULL, dotplot_args = list(), \dots)
}

\arguments{
  \item{x}{an object of class \code{lmertree} or \code{glmertree}.}
  \item{which}{character; \code{"all"} (default), \code{"tree"}, 
      \code{"random"}, \code{"tree.coef"} of \code{"growth"}. Specifies whether the tree, random          
      effects, or both should be plotted. \code{"growth"} should only be used
      in longitudinal models, it yields a tree with growth curves for each of the
      subjects in the tree nodes, instead of individual datapoints, and a 
      thick red curve for the estimated node-specific fixed effectsm representing the
      average trajectory within the terminal node.
      \code{"tree.coef"}           
      yields caterpillar plots of the estimated fixed-effects coefficients
      in every terminal node of the tree, but omits the tree structure 
      (see Details).} 
  \item{nodesize_level}{numeric. At which grouping level should sample size 
      printed above each terminal node be computed? Defaults to 1, which
      is the lowest level of observation. If a value of 2 is specified,
      sample size at the cluster level will be printed above each terminal node.
      This only works if \code{x} (the (g)lmertree) was fitted using 
      the \code{cluster} argument. Alternatively, a character vector of length
      one can be supplied, which gives the name of the grouping indicator
      in the \code{data.frame} used to fit \code{x}.}
  \item{cluster}{vector of cluster ids. Only used if \code{which = "growth"}.
      Need not be specified if clustered covariances were used for partitioning 
      (i.e., argument \code{cluster} was specified). If \code{cluster} was not
      specified in the call to functions \code{(g)lmertree}, this argument 
      should be specified for the plotting function to identify which
      individual observations belong to the same subject.}
  \item{ask}{logical. Should user be asked for input, before a new figure
      is drawn?}
  \item{type}{character; \code{"extended"} (default) or \code{"simple"}.
      \code{"extended"} yields a plotted tree with observed data and/or
      fitted means plotted in the terminal nodes; \code{"simple"} yields a 
      plotted tree with the value of fixed and/or random effects coefficients 
      reported in the terminal nodes.}
  \item{observed}{logical. Should observed datapoints be plotted in the tree?
      Defaults to \code{TRUE}, \code{FALSE} is only supported for objects of 
      class \code{lmertree}, not of class\code{glmertree}.}
  \item{fitted}{character. \code{"combined"} (default), \code{"marginal"} 
      or \code{"none"}. Specifies whether and how fitted values should be 
      computed and visualized. Only used when predictor variables for the 
      node-specific (G)LMs were specified. If \code{"combined"}, fitted 
      values will computed, using observed values of the remaining 
      (random and fixed-effects) predictor variables, which can yield
      very wiggly curves. If \code{"marginal"}, fitted values will be calculated, 
      fixing all remaining predictor variables (with random and/or fixed effects) 
      at the observed sample mean (or majority class).} 
  \item{tp_args}{list of arguments to be passed to panel generating function
      \code{node_glmertree}. See arguments \code{node_bivplot} in 
      \code{\link[partykit]{panelfunctions}}.}
  \item{drop_terminal}{logical. Should all terminal nodes be plotted at the 
      bottom of the plot?}
  \item{terminal_panel}{an optional panel generating function to be passed to 
      \code{plot.party()}, but will most likely be ignored. For passing 
      arguments to the panel generating functions, use argument
      \code{tp_args}. For using a custom panel generating function,
      see Details.}
  \item{dotplot_args}{Optional list of additional arguments to be passed to 
  \code{\link[lattice]{dotplot}}. Only relevant when random- or fixed-effects plots
    are requested by specifying \code{which} as \code{"ranef"}, \code{"all"}, or 
    \code{"ranef"}.} 
  \item{\dots}{Additional arguments to be passed to \code{plot.party()}.
      See \code{\link[partykit]{party-plot}} documentation for details.}
}

  
\details{
  If the node-specific model of the \code{(g)lmertree} object specified by 
  argument \code{x} is an intercept-only model, observed data distributions
  will be plotted in the terminal nodes of the tree (using 
  \code{\link[partykit]{node_barplot}} (for categorical responses) or 
  \code{\link[partykit]{node_boxplot}} (for numerical responses). Otherwise, 
  fitted values will be plotted, in addition to observed datapoints, using a
  function taking similar arguments as \code{\link[partykit]{node_bivplot}}.
  
  Exceptions:
  
  If \code{fitted = "marginal"}, fitted values will be plotted by assuming the
  mean (continuous predictors) or mode (categorical predictors) for all 
  predictor variables, except the variable on the x-axis of the current plot.
  
  If \code{which = "growth"}, individual growth curves will be plotted as thin grey
  lines in the terminal nodes, while the node-specific fixed effect will be plotted
  on top of that as a thicker red curve. 
  
  If \code{which = "tree.coef"}), caterpillar plot(s) are created for the local 
  (node-specific) fixed effects. These depict the estimated fixed-effects 
  coefficients with 95\% confidence intervals, but note that these CIs do not 
  account for the searching of the tree structure and are therefore likely 
  too narrow. There is currently no way to adjust CIs for searching of the tree 
  structure, but the CIs can be useful to obtain an indication of the variability 
  of the coefficient estimates, not for statistical significance testing. 
  
  If \code{which = "ranef"} or \code{"all"}, caterpillar plot(s) for the random 
  effect(s) created, depicting the predicted random effects with 95\% confidence 
  intervals. See also \code{\link[lme4]{ranef}} for more info. Note that the CIs
  do not account for the searching of the tree structure and may be too narrow.
  
  If users want to specify custom panel generating functions, it might be best to
  not use the plotting method for (g)lmertrees. Instead, extract the (g)lmtree from
  the fitted (g)lmertree object (which is a list containing, amongst others, an element 
  \code{$tree}). On this tree, most of the customization options from 
  \code{\link[partykit]{party-plot}} can then be applied.
  
  The code is still under development and might change in future versions.
}


\references{
Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2018).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
    with Generalized Linear Mixed-Effects Model Trees}.
  \emph{Behavior Research Methods}, \bold{50}(5), 2016--2034.
  \doi{10.3758/s13428-017-0971-x}
  
Fokkema M, Zeileis A (2024).
  \dQuote{Subgroup Detection in Linear Growth Curve Models with
    Generalized Linear Mixed Model (GLMM) Trees.}
  \emph{Behavior Research Methods}, \bold{56}(7), 6759--6780. 
  \doi{10.3758/s13428-024-02389-1}
}


\examples{
\donttest{
## load artificial example data
data("DepressionDemo", package = "glmertree")

## fit linear regression LMM tree for continuous outcome
lt <- lmertree(depression ~ treatment + age | cluster | anxiety + duration,
  data = DepressionDemo)
plot(lt)
plot(lt, type = "simple")
plot(lt, which = "tree", fitted = "combined")
plot(lt, which = "tree", fitted = "none")
plot(lt, which = "tree", observed = FALSE)
plot(lt, which = "tree.coef")
plot(lt, which = "ranef")

## fit logistic regression GLMM tree for binary outcome
gt <- glmertree(depression_bin ~ treatment + age | cluster | 
  anxiety + duration, data = DepressionDemo)
plot(gt)  
plot(gt, type = "simple")
plot(gt, which = "tree", fitted = "combined")
plot(gt, which = "tree", fitted = "none")
plot(gt, which = "tree.coef")
plot(gt, which = "ranef")
}
}

\seealso{\code{\link{lmertree}}, \code{\link{glmertree}}, 
    \code{\link[partykit]{party-plot}}.}

\keyword{hplot}
