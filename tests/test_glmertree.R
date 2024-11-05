library("glmertree")
options(width = 70, prompt = "R> ", continue = "+  ", digits = 4)

data("DepressionDemo", package = "glmertree")
summary(DepressionDemo)

lmm_tree <- lmertree(depression ~ treatment | cluster | age + duration + anxiety,
  data = DepressionDemo)
plot(lmm_tree, which = "tree")
plot(lmm_tree, which = "ranef")
coef(lmm_tree)
ranef(lmm_tree)$cluster
predict(lmm_tree, newdata = DepressionDemo[1:7, ])
predict(lmm_tree, newdata = DepressionDemo[1:7, -3], re.form = NA)
residuals(lmm_tree)[1:10]
predict(lmm_tree)[1:10]

glmm_tree <- glmertree(depression_bin ~ treatment | cluster | age + duration + anxiety,
  data = DepressionDemo, family = binomial)
plot(glmm_tree, which = "tree")
plot(glmm_tree, which = "ranef")
coef(glmm_tree)
ranef(glmm_tree)$cluster
predict(glmm_tree, newdata = DepressionDemo[1:7, ])
predict(glmm_tree, newdata = DepressionDemo[1:7, -3], re.form = NA)
residuals(glmm_tree)[1:10]
predict(glmm_tree)[1:10]
