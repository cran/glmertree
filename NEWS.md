# glmertree 0.2-6

* Paper on linear growth curve modeling with `lmertree()` published:
  Fokkema M, Zeileis A (2024). "Subgroup Detection in Linear Growth
  Curve Models with Generalized Linear Mixed Model (GLMM) Trees."
  _Behavior Research Methods_, **56**(7), 6759-6780.
  [doi:10.3758/s13428-024-02389-1](https://doi.org/10.3758/s13428-024-02389-1)

* Reduce precision in tests in order to avoid NOTEs in CRAN checks.


# glmertree 0.2-5

* Added functions `cv.lmertree()` and `cv.glmertree()`, implementing honest estimation of
  model coefficients per Athey & Imbens (2016)

* Fixed bug in `coef()` method `betamertree()`.


# glmertree 0.2-4

* Additional arguments to `lattice::dotplot` can now be passed to 
  `plot.lmertree()`/`plot.glmertree()`.
  
* Dedicated plotting function for partitioning growth curves has
  been added, see `which = "growth"` in `?plot.lmertree` and
  `?plot.glmertree`. 


# glmertree 0.2-3

* Minor bug fixes.

* Experimental function `betamertree()` added, for partitioning 
  mixed-effects beta regression trees.


# glmertree 0.2-1

* Bug fix for missing data: Listwise deletion over all predictor,
  response, and splitting variables is now applied. Otherwise, 
  global and local parts of (g)lmertree model could be estimated 
  based on a different sample size or an error could occur.

* Fixed bug: If dot `.` was used in `formula` for specifying partitioning 
  variables, tree fitting would pick up tree structure from previous
  iteration. Fixed now.

* Fixed bug: When `type = "simple"` in call to function `plot.lmertree()`
  and `plot.glmertree()` is specified, further arguments are now passed
  correctly. 


# glmertree 0.2-0

* Bug fix in `predict()` method: Erroneous allocation of new observations
  to tree nodes could occur, which has been fixed now. 

* New plotting method for (g)lmertrees. Arguments `which`, `type`, and `fitted`
  now support wider range of plots. E.g., `fitted` argument provides 
  different ways to compute fitted values in terminal panel plots. 
  E.g., `type` argument now supports the plotting of fixed-effects 
  coefficients with standard-error bars.

* Two artificial example datasets added to illustate fitting of constant 
  fits (`MHserviceDemo`) and growth-curve models (`GrowthCurveDemo`) in
  terminal nodes. 

* New sections added to vignette, illustrating how mixed-effects 
  regression trees with constant fits and growth curve models in terminal
  nodes can be fitted.


# glmertree 0.1-2

* Functions `glmertree()` and `lmertree()` now take offset arguments.

* Argument `ranefstart` of functions `glmertree()` and `lmertree()` can now be 
  set to `TRUE`. As a result, the random effects will be estimated
  before fitting the tree in the first iteration. This may yield better 
  results when random effects are expected to be substantial.

* Functions `glmertree()` and `lmertree()` now take `cluster` argument, an optional 
  vector with cluster IDs to be employed for clustered covariances in the parameter 
  stability tests.

* Bugs fixed: Argument `dfsplit` is now passed correctly to tree fitting functions,
  additional arguments are now passed correctly to `lmer()` and `glmer()`.

* Arguments `lmer.control` in `lmertree()` and `glmer.control` in `glmertree()`
  are now actually passed to `lmer()` and `glmer()` internally.


# glmertree 0.1-1

* First CRAN release of the `glmertree` package for fitting generalized linear
  mixed-effects model trees in R. For an introduction to the underlying methods see
 
  Fokkema, Smits, Zeileis, Hothorn, Kelderman (2015). Detecting
  Treatment-Subgroup Interactions in Clustered Data with Generalized
  Linear Mixed-Effects Model Trees. Working Paper 2015-10. Working Papers
  in Economics and Statistics, Research Platform Empirical and Experimental
  Economics, Universitaet Innsbruck.
  URL <https://EconPapers.RePEc.org/RePEc:inn:wpaper:2015-10>

  The package is under development on R-Forge at
  <https://R-Forge.R-project.org/projects/partykit/>
