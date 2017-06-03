\name{DepressionDemo}
\alias{DepressionDemo}

\title{Artificial depression treatment dataset}

\description{
Simulated dataset of a randomized clinical trial (N = 150) to illustrate fitting 
of (G)LMM trees.
}

\usage{
data("DepressionDemo")
}

\format{
A data frame containing 150 observations on 6 variables:
\describe{
  \item{depression}{numeric. Continuous treatment outcome variable 
    (range: 3-16, M = 9.12, SD = 2.66).}
  \item{treatment}{factor. Binary treatment variable.}
  \item{cluster}{factor. Indicator for cluster with 10 levels.}
  \item{age}{numeric. Continuous partitioning variable (range: 18-69, 
    M = 45, SD = 9.56).}
  \item{anxiety}{numeric. Continuous partitioning variable (range: 3-18, 
    M = 10.26, SD = 3.05).}
  \item{duration}{numeric. Continuous partitioning variable (range: 1-17, 
    M = 6.97, SD = 2.90).}
  \item{depression_bin}{factor. Binarized treatment outcome variable 
    (0 = recovered, 1 = not recovered.}
}}


\details{
The data were generated such that the duration and anxiety covariates 
characterized three subgroups with differences in treatment effects. The 
\code{cluster} variable was used to introduce a random intercept that should 
be accounted for. The treatment outcome is an index of depressive symptomatology. 
}

\seealso{\code{\link{lmertree}}}

\examples{
\donttest{
data("DepressionDemo", package = "glmertree")
summary(DepressionDemo)
lt <- lmertree(depression ~ treatment | cluster | anxiety + duration + age, 
        data = DepressionDemo)
plot(lt)
gt <- glmertree(depression_bin ~ treatment | cluster | anxiety + duration + age, 
        data = DepressionDemo)
plot(gt)
}
}

\keyword{datasets}