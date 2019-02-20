\name{MHserviceDemo}
\alias{MHserviceDemo}

\title{Artificial mental-health service outcomes dataset}

\description{
Artificial dataset of treatment outcomes (N = 3739) of 13 mental-health
services to illustrate fitting of (G)LMM trees with constant fits in 
terminal nodes. 
}

\usage{
data("MHserviceDemo")
}

\format{
A data frame containing 3739 observations on 8 variables:
\describe{
  \item{age}{numeric. Variable representing age in years (range: 4.8 - 
    23.6, M = 11.46).}
  \item{impact}{numeric. Continuous variable representing severity of 
    and impairment due to mental-health problems at baseline. Higher 
    values indicate higher severity and impairment.}
  \item{gender}{factor. Indicator for gender.}
  \item{emotional}{factor. Indicator for presence of emotional disorder 
    at baseline.}
  \item{autism}{factor. Indicator for presence of autistic disorder at 
    baseline.}
  \item{conduct}{factor. Indicator for mental-health service provider.}
  \item{cluster_id}{factor. Binarized treatment outcome variable 
    (0 = recovered, 1 = not recovered.}
  \item{outcome}{numeric. Variable representing treatment outcome as
    measured by a total mental-health difficulties score assessed about
    6 months after baseline, corrected for the baseline assessment. Higher
    values indicate poorer outcome.}
}}



\details{
Dataset was modelled after Edbrooke-Childs et al. (2017), who analyzed a 
sample of $N = 3,739$ young people who received treatment at one of 13 
mental-health service providers in the UK. Note that the data were 
artificially generated and do not reflect actual empirical findings.
}

\seealso{\code{\link{lmertree}}, \code{\link{glmertree}}}

\examples{
\donttest{
data("MHserviceDemo", package = "glmertree")
summary(MHserviceDemo)
lt <- lmertree(outcome ~ 1 | cluster_id | age + gender + emotional + 
               autism + impact + conduct, data = MHserviceDemo)
plot(lt)

gt <- glmertree(factor(outcome > 0) ~ 1 | cluster_id | age + gender + 
                emotional + autism + impact + conduct, 
                data = MHserviceDemo, family = "binomial")
plot(gt)
}
}

\keyword{datasets}