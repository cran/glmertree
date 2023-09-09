\name{betamertree}
\alias{betamertree}
\alias{print.betamertree}
\alias{ranef.betamertree}
\alias{logLik.betamertree}
\alias{model.frame.betamertree}
\alias{predict.betamertree}
\alias{terms.betamertree}
\alias{as.party.betamertree}
\alias{residuals.betamertree}
\alias{resid.betamertree}
\alias{VarCorr.betamertree}

\encoding{latin1}

\title{Beta Mixed-Effects Regression Trees}

\description{
  Model-based recursive partitioning based on mixed-effects beta
  regression.
}

\usage{
betamertree(formula, data, family = NULL, weights = NULL, cluster = NULL, 
  ranefstart = NULL, offset = NULL, REML = TRUE, joint = TRUE, 
  abstol = 0.001, maxit = 100, dfsplit = TRUE, verbose = FALSE, 
  plot = FALSE, glmmTMB.control = glmmTMB::glmmTMBControl(), \dots)
}

\arguments{
  \item{formula}{formula specifying the response variable and a three-part
    right-hand-side describing the regressors, random effects, and
    partitioning variables, respectively. For details see below.}
  \item{data}{data.frame to be used for estimating the model tree.}
  \item{family}{currently not used. The default beta distribution parameterization 
    of package \code{betareg} is used, see also \code{?glmmTMB::beta_family.}}
  \item{weights}{numeric. An optional numeric vector of weights. Can be a 
    name of a column in data or a vector of length \code{nrow(data)}.}
  \item{cluster}{currently not used.}
  \item{ranefstart}{currently not used.}
  \item{offset}{optional numeric vector to be included in the linear predictor 
    with a coeffcient of one. Note that \code{offset} can be a name of a column 
    in \code{data} or a a numeric vector of length \code{nrow(data)}.}
  \item{joint}{currently not used. Fixed effects from the tree are
    always (re-)estimated jointly along with the random effects.}
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
    passed to funtion \code{glmmTMB()}. See \code{\link[glmmTMB]{glmmTMB}} for details.} 
  \item{glmmTMB.control}{list. An optional list with control
    parameters to be passed to \code{glmmTMB()}.
    See \code{\link[glmmTMB]{glmmTMBControl}} for details.} 
  \item{\dots}{Additional arguments to be passed to \code{lmtree()} or \code{glmtree()}.
    See \code{\link[partykit]{mob_control}} documentation for details.}
}
  
\details{
  Experimental function.
  
  Function betamertree aims to learn a tree where each terminal node is associated with
  different fixed-effects regression coefficients, while adjusting for global 
  random effects (such as a random intercept). It is a generalization of the ideas
  underlying function \code{glmertree}, to allow for detection of 
  subgroups with different fixed-effects parameter estimates, keeping the 
  random effects constant throughout the tree (i.e., random effects are
  estimated globally). The estimation algorithm iterates between (1) estimation 
  of the tree given an offset of random effects, and (2) estimation of the 
  random effects given the tree structure. See Fokkema et al. (2018) for 
  a detailed description.
  
  Where \code{glmertree} uses function \code{glmtree} from package partykit to 
  find the subgroups, and function \code{glmer} from package lme4 to estimate
  the mixed-effects model, \code{betamertree} uses function \code{betatree} from
  package betareg to find the subgroups, and function \code{glmmTMB} from 
  package package glmmTMB to estimate the mixed-effects model. 

  The code is experimental and will change in future versions.
}

\value{
The function returns a list with the following objects:
  \item{tree}{The final \code{betatree}.}
  \item{glmmTMB}{The final \code{glmmTMB} random-effects model.}
  \item{ranef}{The corresponding random effects of \code{glmmTMB}.} 
  \item{varcorr}{The corresponding \code{VarCorr(glmmTMB)}.}
  \item{variance}{The corresponding \code{attr(VarCorr(glmmTMB), "sc")^2}.}
  \item{data}{The dataset specified with the \code{data} argument
    including added auxiliary variables \code{.ranef} and \code{.tree}
    from the last iteration.}
  \item{loglik}{The log-likelihood value of the last iteration.}
  \item{iterations}{The number of iterations used to estimate the \code{betamertree}.} 
  \item{maxit}{The maximum number of iterations specified with the \code{maxit} argument.}
  \item{ranefstart}{The random effects used as an offset, as specified with
    the \code{ranefstart} argument.}
  \item{formula}{The formula as specified with the \code{formula} argument.}
  \item{randomformula}{The formula as specified with the \code{randomformula} argument.}
  \item{abstol}{The prespecified value for the change in log-likelihood to evaluate
    convergence, as specified with the \code{abstol} argument.} 
  \item{mob.control}{A list containing control parameters passed to
    \code{betatree()}, as specified with \dots.}
  \item{glmmTMB.control}{A list containing control parameters passed to
    \code{glmmTMB()}, as specified in the \code{control} argument of
    function \code{glmmTMB}.}
  \item{joint}{Whether the fixed effects from the tree were (re-)estimated jointly along 
    with the random effects, specified with the \code{joint} argument.}  
}

\references{
Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2018).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
  with Generalized Linear Mixed-Effects Model Trees}. Behavior Research 
  Methods, 50(5), 2016-2034. \doi{10.3758/s13428-017-0971-x}
}

\seealso{\code{\link[glmmTMB]{glmmTMB}}, 
\code{\link[betareg]{betatree}}}

\examples{
\donttest{
if (require("betareg") && require("glmmTMB")) {
## load example data
data("ReadingSkills", package = "betareg")
## add random noise (not associated with reading scores)
set.seed(1071)
ReadingSkills$x1 <- rnorm(nrow(ReadingSkills))
ReadingSkills$x2 <- runif(nrow(ReadingSkills))
ReadingSkills$x3 <- factor(rnorm(nrow(ReadingSkills)) > 0)
ReadingSkills$gr <- factor(rep(letters[1:5], length.out = nrow(ReadingSkills))) 

## Fit beta mixed-effects regression tree 
betamer_form <- accuracy ~ iq | gr | dyslexia + x1 + x2 + x3
bmertree <- betamertree(betamer_form, data = ReadingSkills, minsize = 10)
VarCorr(bmertree)
fixef(bmertree)
coef(bmertree)
plot(bmertree)
plot(bmertree, type = "simple")
predict(bmertree, newdata = ReadingSkills[1:5,])
predict(bmertree) ## see ?predict.glmmmTMB for other arguments that can be passed
residuals(bmertree) ## see ?residuals.glmmmTMB for other arguments that can be passed
}
}
}
\keyword{tree}