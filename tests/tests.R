library("glmertree")

## Tests are rounded to three decimals to avoid warnings with noLD tests 
## on CRAN since new lme4 version of Feb. 6 2019

## load artificial example data
data("DepressionDemo", package = "glmertree")

## fit and test lmertree:
lt <- lmertree(depression ~ treatment | cluster | age + anxiety + duration,
               data = DepressionDemo)
round(coef(lt), digits = 3L)
round(ranef(lt)$cluster, digits = 3L)
round(predict(lt, type = "response"), digits = 3L)
round(predict(lt, re.form = NA), digits = 3L)
round(residuals(lt), digits = 3L)
round(attr(VarCorr(lt)$cluster, "stddev"), digits = 3L)


## fit and test glmertree:
gt <- glmertree(depression_bin ~ treatment | cluster | age + anxiety + duration,
                data = DepressionDemo)
round(coef(gt), digits = 3L)
round(ranef(gt)$cluster, digits = 3L)
round(predict(gt, type = "response"), digits = 3L)
round(predict(gt, re.form = NA), digits = 3L)
round(residuals(gt), digits = 3L)
round(attr(VarCorr(gt)$cluster, "stddev"), digits = 3L)