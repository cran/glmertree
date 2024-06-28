\name{GrowthCurveDemo}
\alias{GrowthCurveDemo}

\title{Artificial dataset for partitioning of linear growth curve models}

\description{
Artificial dataset to illustrate fitting of LMM trees with growth curve models  
in the terminal nodes. 
}

\usage{
data("GrowthCurveDemo")
}

\format{
A data frame containing 1250 repeated observations on 250 persons. x1 - x8 are
time-invariant partitioning variables. Thus, they are measurements on the person
(i.e., cluster) level, not on the individual observation level.
\describe{
  \item{person}{numeric. Indicator linking repeated measurements to persons.}
  \item{time}{factor. Indicator for timepoint.}
  \item{y}{numeric. Response variable.}
  \item{x1}{numeric. Potential partitioning variable.}
  \item{x2}{numeric. Potential partitioning variable.}
  \item{x3}{numeric. Potential partitioning variable.}
  \item{x4}{numeric. Potential partitioning variable.}
  \item{x5}{numeric. Potential partitioning variable.}
  \item{x6}{numeric. Potential partitioning variable.}
  \item{x7}{numeric. Potential partitioning variable.}
  \item{x8}{numeric. Potential partitioning variable.}
}}

\references{
Fokkema M & Zeileis A (2024). Subgroup detection in linear growth curve models with 
  generalized linear mixed model (GLMM) trees. Behavior Research Methods. 
  \doi{10.3758/s13428-024-02389-1}
  }

\details{Data were generated so that \code{x1}, \code{x2} and \code{x3} are
true partitioning variables, \code{x4} through \code{x8} are noise
variables. The (potential) partitioning variables are time invariant.
Time-varying covariates can also be included in the model. For partitioning
growth curves these should probably not be potential partitioning variables,
as this could result in observations from the same person ending up in 
different terminal nodes. Thus, time-varying covariates are probably
best included as predictors in the node-specific regression model. E.g.: 
\code{y ~ time + timevarying_cov | person | x1 + x2 + x3 + x4}.}

\seealso{\code{\link{lmertree}}, \code{\link{glmertree}}}

\examples{
\donttest{
data("GrowthCurveDemo", package = "glmertree")
head(GrowthCurveDemo)

## Fit LMM tree with a random intercept w.r.t. person:
form <- y ~ time | person | x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
lt.default <- lmertree(form, data = GrowthCurveDemo)
plot(lt.default, which = "tree") ## yields too large tree
VarCorr(lt.default)

## Account for measurement level of the partitioning variables:
lt.cluster <- lmertree(form, cluster = person, data = GrowthCurveDemo)
plot(lt.cluster, which = "tree") ## yields correct tree
plot(lt.cluster, which = "growth") ## plot individual growth curves not datapoints 
coef(lt.cluster) ## node-specific fixed effects
VarCorr(lt.cluster) ## with smaller trees random effects explain more variance

## Fit LMM tree with random intercept and random slope of time w.r.t. person:
form.s <- y ~ time | (1 + time | person) | x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
lt.s.cluster <- lmertree(form.s, cluster = person, data = GrowthCurveDemo)
plot(lt.s.cluster, which = "tree") ## same tree as before
coef(lt.cluster) 
VarCorr(lt.s.cluster)
}
}

\keyword{datasets}