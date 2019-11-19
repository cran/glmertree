\name{plot.glmertree}
\alias{plot.lmertree}
\alias{plot.glmertree}

\encoding{latin1}

\title{Plotting (Generalized) Linear Mixed Model Trees}

\description{
  \code{plot} method for \code{(g)lmertree} objects.
}

\usage{
\method{plot}{lmertree}(x, which = "all", ask = TRUE, type = "extended", 
    observed = TRUE, fitted = "combined", tp_args = list(), 
    drop_terminal = TRUE, terminal_panel = NULL, \dots)
\method{plot}{glmertree}(x, which = "all", ask = TRUE, type = "extended", 
    observed = TRUE, fitted = "combined", tp_args = list(), 
    drop_terminal = TRUE, terminal_panel = NULL, \dots)
}

\arguments{
  \item{x}{an object of class \code{lmertree} or \code{glmertree}.}
  \item{which}{character; \code{"all"} (default), \code{"tree"}, 
      \code{"random"} or \code{"tree.coef"}. Specifies whether, tree, random          
      effects, or both should be plotted. Alternatively, \code{"tree.coef"}           
      yields caterpillar plots of the estimated fixed-effects coefficients
      in every terminal node of the tree, but omits the tree structure 
      (see Details).}
  \item{ask}{logical. Should user be asked for input, before a new figure
      is drawn?}
  \item{type}{character; \code{"extended"} (default) or \code{"simple"}.
      \code{type = "extended"} yields a plotted tree with observed data and/or
      fitted means plotted in the terminal nodes; \code{"simple"} yields a 
      plottedtree with the value of fixed and/or random effects coefficients 
      reported in the terminal nodes.}
  \item{observed}{logical. Should observed datapoints be plotted in the tree?
      Defaults to \code{TRUE}, \code{FALSE} is only supported for objects of 
      class \code{lmertree}, not of class\code{glmertree}.}
  \item{fitted}{character. \code{"combined"} (default), \code{"marginal"} 
      or \code{"none"}. Specifies whether and how fitted values should be 
      computed and visualized. Only used when predictor variables for the 
      node-specific (G)LMs were specified. If \code{"combined"}, fitted 
      values will computed, based on the observed values of the remaining 
      (random and fixed-effects) predictor variables, and their estimated
      effects. If \code{"marginal"}, fitted values will be calculated, 
      keeping all remaining variables (with random and/or fixed effects) 
      fixed at their (population and sample) means (or majority class).} 
  \item{tp_args}{list of arguments to be passed to panel generating function
      \code{node_glmertree}. See arguments \code{node_bivplot} in 
      \code{\link[partykit]{panelfunctions}}.}
  \item{drop_terminal}{logical. Should all terminal nodes be plotted at the 
      bottom of the plot?}
  \item{terminal_panel}{an optional panel function to be passed to 
      \code{plot.party()}. See \code{\link[partykit]{party-plot}} 
      documentation for details.}
  \item{\dots}{Additional arguments to be passed to \code{plot.party()}.
      See \code{\link[partykit]{party-plot}} documentation for details.}
}

  
\details{
  The caterpillar plot(s) for the local (node-specific) fixed effects (created 
  when \code{which = "tree.coef"}) depict the estimated fixed-effects 
  coefficients with 95\% confidence intervals, but these CIs do not account for
  the searching of the tree structure and are therefore likely too narrow. 
  There is currently no way to adjust CIs for searching of the tree structure,
  but the CIs can be useful to obtain an indication of the variability 
  of the coefficient estimates, not for statistical significance testing. 
  
  The caterpillar plot(s) for the random effect (created if \code{which = "ranef"} 
  or \code{"all"}) depict the predicted random effects with 95\% confidence 
  intervals. See also \code{\link[lme4]{ranef}}.
  
  The code is still under development and might change in future versions.
}


\references{
Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2018).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
  with Generalized Linear Mixed-Effects Model Trees}. Behavior Research 
  Methods, 50(5), 2016-2034. \url{https://doi.org/10.3758/s13428-017-0971-x}
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