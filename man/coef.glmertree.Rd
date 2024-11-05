\name{coef.glmertree}
\alias{coef.lmertree}
\alias{coef.glmertree}
\alias{fixef.lmertree}
\alias{fixef.glmertree}

\title{Obtaining Fixed-Effects Coefficient Estimates of (Generalized) Linear Mixed Model Trees}

\description{
  \code{coef} and \code{fixef} methods for \code{(g)lmertree} objects.
}

\usage{
\method{coef}{lmertree}(object, which = "tree", drop = FALSE, \dots)
\method{fixef}{lmertree}(object, which = "tree", drop = FALSE, \dots)
\method{coef}{glmertree}(object, which = "tree", drop = FALSE, \dots)
\method{fixef}{glmertree}(object, which = "tree", drop = FALSE, \dots)
}

\arguments{
  \item{object}{an object of class \code{lmertree} or \code{glmertree}.}
  \item{which}{character; \code{"tree"} (default) or \code{"global"}. Specifies 
    whether local (tree) or global fixed-effects estimates should be returned.}
  \item{drop}{logical. Only used when \code{which = "tree"}; delete the 
    dimensions of the resulting array if it has only one level?}
  \item{\dots}{Additional arguments, curretnly not used.} 
}

  
\details{
  The code is still under development and might change in future versions.
}

\value{
  If \code{type = "local"}, returns a matrix of estimated local fixed-effects
  coefficients, with a row for every terminal node and a column for every 
  fixed effect. If \code{type = "global"}, returns a numeric vector of 
  estimated global fixed-effects coefficients.
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

## fit LMM tree with local fixed effects only
lt <- lmertree(depression ~ treatment + age | cluster | anxiety + duration,
  data = DepressionDemo)
coef(lt)

## fit LMM tree including both local and global fixed effect
lt <- lmertree(depression ~ treatment | (age + (1|cluster)) | anxiety + duration,
  data = DepressionDemo)
coef(lt, which = "tree") # default behaviour
coef(lt, which = "global")


## fit GLMM tree with local fixed effects only
gt <- glmertree(depression_bin ~ treatment | cluster | 
  age + anxiety + duration, data = DepressionDemo)
coef(gt)

## fit GLMM tree including both local and global fixed effect
gt <- glmertree(depression_bin ~ treatment | (age + (1|cluster)) | 
  anxiety + duration, data = DepressionDemo)
coef(gt, which = "tree") # default behaviour
coef(gt, which = "global")
}
}

\seealso{\code{\link{lmertree}}, \code{\link{glmertree}}, 
    \code{\link[partykit]{party-plot}}.}
