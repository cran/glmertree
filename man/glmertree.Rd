\name{glmertree}
\alias{lmertree}
\alias{glmertree}

\alias{print.lmertree}
\alias{plot.lmertree}
\alias{coef.lmertree}
\alias{ranef.lmertree}
\alias{logLik.lmertree}
\alias{model.frame.lmertree}
\alias{predict.lmertree}
\alias{terms.lmertree}
\alias{as.party.lmertree}
\alias{residuals.lmertree}
\alias{resid.lmertree}

\alias{print.glmertree}
\alias{plot.glmertree}
\alias{coef.glmertree}
\alias{ranef.glmertree}
\alias{logLik.glmertree}
\alias{model.frame.glmertree}
\alias{predict.glmertree}
\alias{terms.glmertree}
\alias{as.party.glmertree}
\alias{residuals.glmertree}
\alias{resid.glmertree}

\encoding{latin1}

\title{(Generalized) Linear Mixed Model Trees}

\description{
  Model-based recursive partitioning based on (generalized)
  linear mixed models.
}

\usage{
lmertree(formula, data, weights = NULL,
  ranefstart = NULL, abstol = 0.001, maxit = 100, 
  joint = TRUE, dfsplit = TRUE, verbose = FALSE, plot = FALSE,
  lmer.control = lmerControl(), \dots)

glmertree(formula, data, family = "binomial", weights = NULL,
  ranefstart = NULL, abstol = 0.001, maxit = 100, 
  joint = TRUE, dfsplit = TRUE, verbose = FALSE, plot = FALSE,
  glmer.control = glmerControl(), \dots)
}

\arguments{
  \item{formula}{formula specifying the response variable and a three-part
    right-hand-side describing the regressors, random effects, and
    partitioning variables, respectively. For details see below.}
  \item{data}{data.frame to be used for estimating the model tree.}
  \item{family}{family specification for \code{glmtree} and \code{glmer}.
    See \code{\link[stats]{glm}} documentation for families.}
  \item{weights}{numeric. An optional numeric vector of weights. (Note that
    this is passed with standard evaluation, i.e., it is not enough to pass
    the name of a column in \code{data}.)}
  \item{ranefstart}{numeric. A vector of length \code{nrow(data)}, to be used as
    an offset in estimation of the first tree. \code{NULL} by default,
    which results in a zero offset initialization.}
  \item{abstol}{numeric. The convergence criterion used for estimation of the model.
    When the difference in log-likelihoods of the random-effects model from two
    consecutive iterations is smaller than \code{abstol}, estimation of the
    model tree has converged.} 
  \item{maxit}{numeric. The maximum number of iterations to be performed in
    estimation of the model tree.}
  \item{joint}{logical. Should the fixed effects from the tree be (re-)estimated
    jointly along with the random effects?}
  \item{dfsplit}{logical or numeric. \code{as.integer(dfsplit)} is the
    degrees of freedom per selected split employed when extracting
    the log-likelihood.}
  \item{verbose}{Should the log-likelihood value of the estimated
    random-effects model be printed for every iteration of the estimation?} 
  \item{plot}{Should the tree be plotted at every iteration of the estimation?
    Note that selecting this option slows down execution of the function.}
  \item{lmer.control, glmer.control}{list. An optional list with control
    parameters to be passed to \code{lmer()} and \code{glmer()}, respectively.
    See \code{\link[lme4]{lmerControl}} and \code{\link[lme4]{glmerControl}} for details.} 
  \item{\dots}{Additional arguments to be passed to \code{lmtree()} or \code{glmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  
\details{
  (G)LMM trees learn a tree where each terminal node is associated with
  different regression coefficients while adjusting for global random effects
  (such as a random intercept). This allows for detection of subgroup-specific
  fixed effects, keeping the random effects constant throughout the tree.
  The estimation algorithm iterates between (1) estimation of the tree given
  an offset of random effects, and (2) estimation of a random-effects model
  given the tree structure. See Fokkema et al. (2015) for a detailed introduction.
  
  To specify all variables in the model a \code{formula} such as
  \code{y ~ x1 + x2 | random | z1 + z2 + z3} is used, where \code{y} is the
  response, \code{x1} and \code{x2} are the regressors in every node of the
  tree, \code{random} is the random effect, and \code{z1} to \code{z3} are
  the partitioning variables considered for growing the tree. If \code{random}
  is only a single variable such as \code{id} a random intercept with respect
  to \code{id} is used. Alternatively, it may be an explicit random-effects
  formula such as \code{(1 | id)} or a more complicated formula. (Note that
  in the latter case, the brackets are necessary to protect the pipes in the
  random effects formulation.)

  In the random-effects model from step (2), two strategies are available:
  Either the fitted values from the tree can be supplied as an offset
  (\code{joint = FALSE}) so that only the random effects are estimated.
  Or the fixed effects are (re-)estimated along with the random effects
  using a nesting factor with nodes from the tree (\code{joint = TRUE}).
  In the former case, the estimation of each random-effects model is typically
  faster but more iterations are required.

  The code is still under development and might change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{tree}{The final \code{lmtree}/\code{glmtree}.}
  \item{lmer}{The final \code{lmer} random-effects model.}
  \item{ranef}{The corresponding random effects of \code{lmer}.} 
  \item{varcorr}{The corresponding \code{VarCorr(lmer)}.}
  \item{variance}{The corresponding \code{attr(VarCorr(lmer), "sc")^2}.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.ranef} and \code{.tree}
    from the last iteration.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{iterations}{The number of iterations used to estimate the \code{lmertree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{ranefstart}{The random effects used as an offset, as specified with
    the \code{ranefstart} argument.}
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{randomformula}{The formula as specified with the \code{randomformula} argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{mob.control}{A list containing control parameters passed to
    \code{lmtree()}, as specified with \dots.}
  \item{lmer.control}{A list containing control parameters passed to
    \code{lmer()}, as specified in the \code{lmer.control} argument.}
  \item{joint}{Whether the fixed effects from the tree were (re-)estimated jointly along 
    with the random effects, specified with the \code{joint} argument.}  
}

\references{
Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2015).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
  with Generalized Linear Mixed-Effects Model Trees}.
  Working Paper 2015-10. Working Papers in Economics and Statistics,
  Research Platform Empirical and Experimental Economics, Universität Innsbruck.
  \url{http://EconPapers.RePEc.org/RePEc:inn:wpaper:2015-10}
}

\seealso{\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\examples{
## artificial example data
data("DepressionDemo", package = "glmertree")

## fit normal linear regression LMM tree for continuous outcome
lt <- lmertree(depression ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo)
print(lt)
plot(lt, which = "all") # default behavior, which may also be "tree" or "ranef" 
coef(lt)
ranef(lt)
predict(lt, type = "response") # default behavior, type may also be "node" 
residuals(lt)

\donttest{
## fit logistic regression GLMM tree for binary outcome
gt <- glmertree(depression_bin ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo)
print(gt)
plot(gt, which = "all") # default behavior, which may also be "tree" or "ranef" 
coef(gt)
ranef(gt)
predict(gt, type = "response") # default behavior, type may also be "node" or "link"
residuals(gt)
}
}

\keyword{tree}
