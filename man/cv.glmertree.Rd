\name{cv.glmertree}
\alias{cv.glmertree}
\alias{cv.lmertree}

\encoding{latin1}

\title{Cross Validation of (Generalized) Linear Mixed Model Trees}

\description{
  Performs cross-validation of a model-based recursive partition based on (generalized)
  linear mixed models. Using the tree or subgroup structure estimated from a training dataset,
  the full mixed-effects model parameters are re-estimated using a new set of test observations,
  providing valid computation of standard errors and valid inference. The approach is inspired
  by Athey & Imbens (2016), and "enables the construction of valid confidence intervals [...] whereby 
  one sample is used to construct the partition and another to estimate [...] effects for each 
  subpopulation."
}

\usage{
cv.lmertree(tree, newdata, reference = NULL, omit.intercept = FALSE, \dots)

cv.glmertree(tree, newdata, reference = NULL, omit.intercept = FALSE, \dots)
}

\arguments{
  \item{tree}{An object of class \code{lmertree} or \code{glmertree} that was fitted on 
  a set of training data.}
  \item{newdata}{A \code{data.frame} containing a new set of observations on the same
  variables that were used to fit \code{tree}.}
  \item{reference}{Numeric or character scalar, indicating the number of the terminal node of 
which the intercept should be taken as a reference for intercepts in all other nodes. 
If \code{NULL}, the default of taking the first terminal node's intercept as the reference category
will be used. If the interest is in testing significance of differences between the different
nodes intercepts, this can be overruled by specifying the number of the terminal node
that should be used as the reference category.}
  \item{omit.intercept}{Logical scalar, indicating whether the intercept should be omitted from the model.
The default (\code{FALSE}) includes the intercept of the first terminal node as the intercept and allows
for significance testing of the differences between the first and the other terminal node's intercepts.
Specifying \code{TRUE} will test the value of each terminal node's intercept against zero.}
  \item{\dots}{Not currently used.}
}

\details{The approach is inspired
  by Athey & Imbens (2016), and "enables the construction of valid confidence intervals [...] whereby 
  one sample is used to construct the partition and another to estimate [...] effects for each 
  subpopulation."}

\value{An object of with classes \code{lmertree} and \code{cv.lmertree}, or \code{glmertree} and \code{cv.glmertree}. It is the original (g)lmertree specified by argument \code{tree}, but the parametric model model estimated based on the data specified by argument \code{newdata}. The default S3 methods for classes \code{lmertree} and \code{glmertree} can be used to inspect the results: \code{plot}, \code{predict}, \code{coef}, \code{fixef}, \code{ranef} and \code{VarCorr}. In addition, there is a dedicated \code{summary} method for classes \code{cv.lmertree} and \code{cv.glmertree}, which prints valid parameter estimates and standard errors, resulting from \code{\link[lme4]{summary.merMod}}. For objects of clas \code{cv.lmertree}, hypothesis tests (i.e., p-values) can be obtained by loading package \code{lmerTest} PRIOR to loading package(s) \code{glmertree} (and \code{lme4}), see examples.
}

\references{

Athey S & Imbens G (2016). Recursive partitioning for heterogeneous causal effects. 
  \dQuote{Proceedings of the National Academy of Sciences, 113(27), 7353-7360.}
  \doi{10.1073/pnas.1510489113}

Fokkema M, Smits N, Zeileis A, Hothorn T, Kelderman H (2018).
  \dQuote{Detecting Treatment-Subgroup Interactions in Clustered Data
  with Generalized Linear Mixed-Effects Model Trees}. Behavior Research 
  Methods, 50(5), 2016-2034. \doi{10.3758/s13428-017-0971-x}
  
Fokkema M, Edbrooke-Childs J & Wolpert M (2021). \dQuote{Generalized linear mixed-model 
  (GLMM) trees: A flexible decision-tree method for multilevel and longitudinal data.} 
  Psychotherapy Research, 31(3), 329-341. \doi{10.1080/10503307.2020.1785037}
  
Fokkema M & Zeileis A (2024). Subgroup detection in linear growth curve models with 
  generalized linear mixed model (GLMM) trees. Behavior Research Methods. 
  \doi{10.3758/s13428-024-02389-1}
}

\examples{

\donttest{
require("lmerTest") ## load BEFORE lme4 and glmertree to obtain hypothesis tests / p-values

## Create artificial training and test datasets
set.seed(42)
train <- sample(1:nrow(DepressionDemo), size = 200, replace = TRUE)
test <- sample(1:nrow(DepressionDemo), size = 200, replace = TRUE)

## Fit tree on training data
tree1 <- lmertree(depression ~ treatment | cluster | age + anxiety + duration,
                 data = DepressionDemo[train, ])
                 
## Obtain honest estimates of parameters and standard errors using test data
tree2 <- cv.lmertree(tree1, newdata = DepressionDemo[test, ])
tree3 <- cv.lmertree(tree1, newdata = DepressionDemo[test, ], 
                     reference = 7, omit.intercept = TRUE)

summary(tree2)
summary(tree3)

coef(tree1)
coef(tree2)
coef(tree3)

plot(tree1, which = "tree")
plot(tree2, which = "tree")
plot(tree3, which = "tree")

predict(tree1, newdata = DepressionDemo[1:5, ])
predict(tree2, newdata = DepressionDemo[1:5, ])
}}


\keyword{tree}

\seealso{\code{\link[lme4]{lmer}}, \code{\link[lme4]{glmer}},
  \code{\link{lmertree}}, \code{\link{glmertree}}, 
  \code{\link[lme4]{summary.merMod}}}