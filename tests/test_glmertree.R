library("glmertree")
options(width = 70, prompt = "R> ", continue = "+  ")
data("DepressionDemo", package = "glmertree")
summary(DepressionDemo)
lmm_tree <- lmertree(depression ~ treatment | cluster | 
                       age + duration + anxiety, data = DepressionDemo)
plot(lmm_tree, which = "tree")
plot(lmm_tree, which = "ranef")
round(coef(lmm_tree), digits = 6)
round(ranef(lmm_tree)$cluster, digits = 6)
formatC(predict(lmm_tree, newdata = DepressionDemo[1:7, ]), format = "f", 
        digits = 7)
formatC(predict(lmm_tree, newdata = DepressionDemo[1:7, -3], re.form = NA),
        format = "f", digits = 7)
formatC(residuals(lmm_tree)[1:10], format = "f", digits = 7)
formatC(predict(lmm_tree)[1:10], format = "f", digits = 7)

glmm_tree <- glmertree(depression_bin ~ treatment | 
                         cluster | age + duration + anxiety, data = DepressionDemo, 
                       family = binomial)
plot(glmm_tree, which = "tree")
plot(glmm_tree, which = "ranef")
round(coef(glmm_tree), digits = 6)
round(ranef(glmm_tree)$cluster, digits = 6)
formatC(predict(glmm_tree, newdata = DepressionDemo[1:7, ]), format = "f", 
        digits = 7)
formatC(predict(glmm_tree, newdata = DepressionDemo[1:7, -3], re.form = NA), 
        format = "f", digits = 7)
formatC(residuals(glmm_tree)[1:10], format = "f", digits = 7)
formatC(predict(glmm_tree)[1:10], format = "f", digits = 7)