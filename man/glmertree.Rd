\name{glmertree}
\alias{lmertree}
\alias{glmertree}

\alias{print.lmertree}
\alias{ranef.lmertree}
\alias{logLik.lmertree}
\alias{model.frame.lmertree}
\alias{predict.lmertree}
\alias{terms.lmertree}
\alias{as.party.lmertree}
\alias{residuals.lmertree}
\alias{resid.lmertree}
\alias{VarCorr.lmertree}

\alias{print.glmertree}
\alias{ranef.glmertree}
\alias{logLik.glmertree}
\alias{model.frame.glmertree}
\alias{predict.glmertree}
\alias{terms.glmertree}
\alias{as.party.glmertree}
\alias{residuals.glmertree}
\alias{resid.glmertree}
\alias{VarCorr.glmertree}

\title{(Generalized) Linear Mixed Model Trees}

\description{
  Model-based recursive partitioning based on (generalized)
  linear mixed models.
}

\usage{
lmertree(formula, data, weights = NULL, cluster = NULL, 
  ranefstart = NULL, offset = NULL, joint = TRUE, 
  abstol = 0.001, maxit = 100, dfsplit = TRUE, verbose = FALSE, 
  plot = FALSE, REML = TRUE, lmer.control = lmerControl(), \dots)

glmertree(formula, data, family = "binomial", weights = NULL,
  cluster = NULL, ranefstart = NULL, offset = NULL, joint = TRUE,
  abstol = 0.001, maxit = 100, dfsplit = TRUE, verbose = FALSE, 
  plot = FALSE, nAGQ = 1L, glmer.control = glmerControl(), \dots)
}

\arguments{
  \item{formula}{formula specifying the response variable and a three-part
    right-hand-side describing the regressors, random effects, and
    partitioning variables, respectively. For details see below.}
  \item{data}{data.frame to be used for estimating the model tree.}
  \item{family}{family specification for \code{glmtree} and \code{glmer}.
    See \code{\link[stats]{glm}} documentation for families.}
  \item{weights}{numeric. An optional numeric vector of weights. Can be a 
    name of a column in data or a vector of length \code{nrow(data)}.}
  \item{cluster}{optional vector of cluster IDs to be employed for clustered 
    covariances in the parameter stability tests. Can be a name of a column 
    in \code{data} or a vector of length \code{nrow(data)}. If 
    \code{cluster = NULL} (the default), observation-level covariances 
    are employed in the parameter stability tests. If 
    partitioning variables are measured on the cluster level, this can
    be accounted for by specifying the name of the cluster argument here, 
    as a result cluster-level covariances will be employed in the parameter 
    stability tests.}
  \item{ranefstart}{\code{NULL} (the default), \code{TRUE}, or a numeric 
    vector of length \code{nrow(data)}. Specifies the offset to be used in 
    estimation of the first tree. \code{NULL} by default, yielding a zero offset 
    initialization. If \code{ranefstart = TRUE} is specified, the random effects 
    will be estimated first and the first tree will be grown using the 
    random-effects predictions as an offset.}
  \item{offset}{optional numeric vector to be included in the linear predictor 
    with a coeffcient of one. Note that \code{offset} can be a name of a column 
    in \code{data} or a a numeric vector of length \code{nrow(data)}.}
  \item{joint}{logical. Should the fixed effects from the tree be (re-)estimated
    jointly along with the random effects?}
  \item{abstol}{numeric. The convergence criterion used for estimation of the model.
    When the difference in log-likelihoods of the random-effects model from two
    consecutive iterations is smaller than \code{abstol}, estimation of the
    model tree has converged.} 
  \item{maxit}{numeric. The maximum number of iterations to be performed in
    estimation of the model tree.}
  \item{dfsplit}{logical or numeric. \code{as.integer(dfsplit)} is the
    degrees of freedom per selected split employed when extracting
    the log-likelihood.}
  \item{verbose}{Should the log-likelihood value of the estimated
    random-effects model be printed for every iteration of the estimation?} 
  \item{plot}{Should the tree be plotted at every iteration of the estimation?
    Note that selecting this option slows down execution of the function.}
  \item{REML}{logical scalar. Should the fixed-effects estimates be chosen to 
    optimize the REML criterion (as opposed to the log-likelihood)? Will be 
    passed to funtion \code{lmer()}. See \code{\link[lme4]{lmer}} for details.} 
  \item{nAGQ}{integer scalar. Specifies the number of points per axis for evaluating 
    the adaptive Gauss-Hermite approximation to the log-likelihood, to be passed
    to function \code{glmer()}. See \code{\link[lme4]{glmer}} for details.} 
  \item{lmer.control, glmer.control}{list. An optional list with control
    parameters to be passed to \code{lmer()} or \code{glmer()}, respectively.
    See \code{\link[lme4]{lmerControl}} for details.} 
  \item{\dots}{Additional arguments to be passed to \code{lmtree()} or \code{glmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  
\details{
  (G)LMM trees learn a tree where each terminal node is associated with
  different fixed-effects regression coefficients while adjusting for global 
  random effects (such as a random intercept). This allows for detection of 
  subgroups with different fixed-effects parameter estimates, keeping the 
  random effects constant throughout the tree (i.e., random effects are
  estimated globally). The estimation algorithm iterates between (1) estimation 
  of the tree given an offset of random effects, and (2) estimation of the 
  random effects given the tree structure. See Fokkema et al. (2018) for 
  a detailed introduction.
  
  To specify all variables in the model a \code{formula} such as
  \code{y ~ x1 + x2 | random | z1 + z2 + z3} is used, where \code{y} is the
  response, \code{x1} and \code{x2} are the regressors in every node of the
  tree, \code{random} is the random effects, and \code{z1} to \code{z3} are
  the partitioning variables considered for growing the tree. If \code{random}
  is only a single variable such as \code{id} a random intercept with respect
  to \code{id} is used. Alternatively, it may be an explicit random-effects
  formula such as \code{(1 | id)} or a more complicated formula such as
  \code{((1+time) | id)}. (Note that in the latter two formulas, the brackets 
  are necessary to protect the pipes in the random-effects formulation.)

  In the random-effects model from step (2), two strategies are available:
  Either the fitted values from the tree can be supplied as an offset
  (\code{joint = FALSE}) so that only the random effects are estimated.
  Or the fixed effects are (re-)estimated along with the random effects
  using a nesting factor with nodes from the tree (\code{joint = TRUE}).
  In the former case, the estimation of each random-effects model is typically
  faster, but more iterations are required.

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
Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2018).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
    with Generalized Linear Mixed-Effects Model Trees}.
  \emph{Behavior Research Methods}, \bold{50}(5), 2016--2034.
  \doi{10.3758/s13428-017-0971-x}
  
Fokkema M, Edbrooke-Childs J, Wolpert M (2021).
  \dQuote{Generalized Linear Mixed-Model (GLMM) Trees: A Flexible
    Decision-Tree Method for Multilevel and Longitudinal Data.} 
  \emph{Psychotherapy Research}, \bold{31}(3), 329--341.
  \doi{10.1080/10503307.2020.1785037}
  
Fokkema M, Zeileis A (2024).
  \dQuote{Subgroup Detection in Linear Growth Curve Models with
    Generalized Linear Mixed Model (GLMM) Trees.}
  \emph{Behavior Research Methods}, \bold{56}(7), 6759--6780. 
  \doi{10.3758/s13428-024-02389-1}
}

\seealso{\code{\link{plot.lmertree}}, \code{\link{plot.glmertree}},
  \code{\link{cv.lmertree}}, \code{\link{cv.glmertree}},
  \code{\link{GrowthCurveDemo}},
  \code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link[partykit]{lmtree}}, \code{\link[partykit]{glmtree}}}

\examples{
\donttest{
## artificial example data
data("DepressionDemo", package = "glmertree")

## fit normal linear regression LMM tree for continuous outcome
lt <- lmertree(depression ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo)
print(lt)
plot(lt, which = "all") # default behavior, may also be "tree" or "ranef" 
coef(lt)
ranef(lt)
predict(lt, type = "response") # default behavior, may also be "node"
predict(lt, re.form = NA) # excludes random effects, see ?lme4::predict.merMod
residuals(lt)
VarCorr(lt) # see lme4::VarCorr


## fit logistic regression GLMM tree for binary outcome
gt <- glmertree(depression_bin ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo)
print(gt)
plot(gt, which = "all") # default behavior, may also be "tree" or "ranef" 
coef(gt)
ranef(gt)
predict(gt, type = "response") # default behavior, may also be "node" or "link"
predict(gt, re.form = NA) # excludes random effects, see ?lme4::predict.merMod
residuals(gt)
VarCorr(gt) # see lme4::VarCorr

## Alternative specification for binomial family: no. of successes and failures
DepressionDemo$failures <- as.numeric(DepressionDemo$depression_bin) - 1
DepressionDemo$successes <- 1 - DepressionDemo$failures
gt <- glmertree(cbind(failures, successes) ~ treatment | cluster | age + anxiety + duration,
  data = DepressionDemo, ytype = "matrix") ## see also ?partykit::mob_control
}
}
\keyword{tree}
