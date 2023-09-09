### R code from vignette source 'glmertree.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("glmertree")
options(width = 70, prompt = "R> ", continue = "+  ")


###################################################
### code chunk number 2: glmertree.Rnw:61-62 (eval = FALSE)
###################################################
## install.packages("glmertree")


###################################################
### code chunk number 3: glmertree.Rnw:67-68 (eval = FALSE)
###################################################
## install.packages("glmertree", repos = "http://R-Forge.R-project.org")


###################################################
### code chunk number 4: glmertree.Rnw:73-74 (eval = FALSE)
###################################################
## library("glmertree")


###################################################
### code chunk number 5: glmertree.Rnw:92-94
###################################################
data("MHserviceDemo", package = "glmertree")
summary(MHserviceDemo)


###################################################
### code chunk number 6: glmertree.Rnw:99-101
###################################################
MH_tree <- lmertree(outcome ~ 1 | cluster_id | age + gender + emotional + 
                      autism + impact + conduct, data = MHserviceDemo)


###################################################
### code chunk number 7: glmertree.Rnw:108-111 (eval = FALSE)
###################################################
## MH_tree2 <- lmertree(outcome ~ 1 | age + (1 | cluster_id) | gender + 
##                        emotional + autism + impact + conduct, 
##                      data = MHserviceDemo)


###################################################
### code chunk number 8: glmertree.Rnw:118-122 (eval = FALSE)
###################################################
## MHserviceDemo$outcome_bin <- factor(MHserviceDemo$outcome > 0)
## MH_gtree <- glmertree(outcome_bin ~ 1 | cluster_id | age + gender + 
##                         emotional + autism + impact + conduct, 
##                       data = MHserviceDemo, family = "binomial")


###################################################
### code chunk number 9: glmertree.Rnw:127-128 (eval = FALSE)
###################################################
## plot(MH_tree)


###################################################
### code chunk number 10: glmertree.Rnw:135-136
###################################################
plot(MH_tree, which = "tree")


###################################################
### code chunk number 11: glmertree.Rnw:153-154
###################################################
plot(MH_tree, which = "ranef")


###################################################
### code chunk number 12: glmertree.Rnw:162-167 (eval = FALSE)
###################################################
## print(MH_tree)
## coef(MH_tree)
## fixef(MH_tree)
## ranef(MH_tree)
## VarCorr(MH_tree)


###################################################
### code chunk number 13: glmertree.Rnw:172-173
###################################################
predict(MH_tree, newdata = MHserviceDemo[1:10,])


###################################################
### code chunk number 14: glmertree.Rnw:178-179
###################################################
predict(MH_tree, newdata = MHserviceDemo[1:10, -7], re.form = NA)


###################################################
### code chunk number 15: glmertree.Rnw:187-191 (eval = FALSE)
###################################################
## resids <- residuals(MH_tree)
## preds <- predict(MH_tree)
## plot(MHserviceDemo$cluster_id, resids)
## scatter.smooth(preds, resids)


###################################################
### code chunk number 16: glmertree.Rnw:197-200
###################################################
resids <- residuals(MH_tree)
preds <- predict(MH_tree)
plot(MHserviceDemo$cluster_id, resids, xlab = "Cluster", ylab = "Residuals")


###################################################
### code chunk number 17: glmertree.Rnw:202-203
###################################################
scatter.smooth(preds, resids, xlab = "Predicted values", ylab = "Residuals")


###################################################
### code chunk number 18: glmertree.Rnw:217-219
###################################################
data("DepressionDemo", package = "glmertree")
summary(DepressionDemo)


###################################################
### code chunk number 19: glmertree.Rnw:226-228
###################################################
lmmt <- lmertree(depression ~ treatment | cluster | age + 
                   duration + anxiety, data = DepressionDemo)


###################################################
### code chunk number 20: glmertree.Rnw:235-237 (eval = FALSE)
###################################################
## depression ~ treatment | (age + (1 + age | cluster)) | age + 
##   duration + anxiety


###################################################
### code chunk number 21: glmertree.Rnw:244-245
###################################################
depression ~ treatment | (1|center/area) | age + duration + anxiety


###################################################
### code chunk number 22: glmertree.Rnw:250-251 (eval = FALSE)
###################################################
## plot(lmmt)


###################################################
### code chunk number 23: glmertree.Rnw:258-259
###################################################
plot(lmmt, which = "tree")


###################################################
### code chunk number 24: glmertree.Rnw:272-273
###################################################
plot(lmmt, which = "ranef")


###################################################
### code chunk number 25: glmertree.Rnw:282-283 (eval = FALSE)
###################################################
## plot(lmmt, which = "tree.coef")


###################################################
### code chunk number 26: glmertree.Rnw:289-290
###################################################
plot(lmmt, which = "tree.coef")


###################################################
### code chunk number 27: glmertree.Rnw:299-304 (eval = FALSE)
###################################################
## print(lmmt)
## coef(lmmt)
## fixef(lmmt)
## ranef(lmmt)
## VarCorr(lmmt)


###################################################
### code chunk number 28: glmertree.Rnw:309-310
###################################################
predict(lmmt, newdata = DepressionDemo[1:7,])


###################################################
### code chunk number 29: glmertree.Rnw:315-316
###################################################
predict(lmmt, newdata = DepressionDemo[1:7, -3], re.form = NA)


###################################################
### code chunk number 30: glmertree.Rnw:324-328 (eval = FALSE)
###################################################
## resids <- residuals(lmmt)
## preds <- predict(lmmt)
## plot(factor(DepressionDemo$cluster), resids)
## scatter.smooth(preds, resids)


###################################################
### code chunk number 31: glmertree.Rnw:334-337
###################################################
resids <- residuals(lmmt)
preds <- predict(lmmt)
plot(factor(DepressionDemo$cluster), resids, xlab = "Cluster", ylab = "Residuals")


###################################################
### code chunk number 32: glmertree.Rnw:339-340
###################################################
scatter.smooth(preds, resids, xlab = "Predicted values", ylab = "Residuals")


###################################################
### code chunk number 33: glmertree.Rnw:347-349
###################################################
fligner.test(resids ~ DepressionDemo$cluster)
bartlett.test(resids ~ DepressionDemo$cluster)


###################################################
### code chunk number 34: glmertree.Rnw:363-366
###################################################
data("GrowthCurveDemo", package = "glmertree")
dim(GrowthCurveDemo)
names(GrowthCurveDemo)


###################################################
### code chunk number 35: glmertree.Rnw:375-377
###################################################
gc_tree <- lmertree(y ~ time | person | x1 + x2 + x3 + x4 + x5 + x6 + 
                      x7 + x8, cluster = person, data = GrowthCurveDemo)


###################################################
### code chunk number 36: glmertree.Rnw:382-383
###################################################
width(gc_tree$tree)


###################################################
### code chunk number 37: glmertree.Rnw:388-391
###################################################
gc_obs_tree <- lmertree(y ~ time | person | x1 + x2 + x3 + x4 + x5 + 
                          x6 + x7 + x8, data = GrowthCurveDemo)
width(gc_obs_tree$tree)


###################################################
### code chunk number 38: glmertree.Rnw:396-397 (eval = FALSE)
###################################################
## plot(gc_tree, which = "tree", fitted = "marginal")


###################################################
### code chunk number 39: glmertree.Rnw:404-405 (eval = FALSE)
###################################################
## plot(gc_tree, which = "growth")


###################################################
### code chunk number 40: glmertree.Rnw:411-412
###################################################
plot(gc_tree, which = "growth")


###################################################
### code chunk number 41: glmertree.Rnw:421-422 (eval = FALSE)
###################################################
## plot(gc_tree, type = "simple", which = "tree")


###################################################
### code chunk number 42: glmertree.Rnw:428-429
###################################################
plot(gc_tree, type = "simple", which = "tree")


###################################################
### code chunk number 43: glmertree.Rnw:436-437 (eval = FALSE)
###################################################
## plot(gc_tree, which = "tree.coef")


###################################################
### code chunk number 44: glmertree.Rnw:443-444
###################################################
plot(gc_tree, which = "tree.coef")


###################################################
### code chunk number 45: glmertree.Rnw:453-455
###################################################
varcor <- VarCorr(gc_tree)
varcor


###################################################
### code chunk number 46: glmertree.Rnw:460-464
###################################################
res_var <- attr(varcor, "sc")^2
int_var <- as.numeric(varcor$person)
ICC <- int_var / (res_var + int_var)
ICC


###################################################
### code chunk number 47: glmertree.Rnw:474-477
###################################################
form_s <- formula(paste0("y ~ time | (1 + time | person) | ", 
                         paste0("x", 1:8, collapse = " + ")))
form_s


###################################################
### code chunk number 48: glmertree.Rnw:482-483
###################################################
gc_tree_s <- lmertree(form_s, cluster = person, data = GrowthCurveDemo)


###################################################
### code chunk number 49: glmertree.Rnw:488-489
###################################################
VarCorr(gc_tree_s)


###################################################
### code chunk number 50: appendix1
###################################################
set.seed(123)
treatment <- rbinom(n = 150, size = 1, prob = .5)
duration <- round(rnorm(150, mean = 7, sd = 3))
anxiety <- round(rnorm(150, mean = 10, sd = 3))
age <- round(rnorm(150, mean = 45, sd = 10))
error <- rnorm(150, 0, 2)


###################################################
### code chunk number 51: appendix2
###################################################
cluster <- error + rnorm(150, 0, 6)
rand_int <- sort(rep(rnorm(10, 0, 1), each = 15))
rand_int[order(cluster)] <- rand_int 
error <- error - rand_int
cluster[order(cluster)] <- rep(1:10, each = 15)


###################################################
### code chunk number 52: appendix3
###################################################
node3t1 <- ifelse(duration <= 8 & anxiety <= 10 & treatment == 0, -2, 0)
node3t2 <- ifelse(duration <= 8 & anxiety <= 10 & treatment == 1, 2, 0)
node5t1 <- ifelse(duration > 8 & treatment == 0, 2.5, 0)
node5t2 <- ifelse(duration > 8 & treatment == 1, -2.5, 0)


###################################################
### code chunk number 53: appendix4
###################################################
depression <- round(9 + node3t1 + node3t2 + node5t1 + node5t2 +
  .4 * treatment + error + rand_int)
depression_bin <- factor(as.numeric(depression > 9))


###################################################
### code chunk number 54: appendix5
###################################################
treatment <- factor(treatment, labels = c("Treatment 1", "Treatment 2"))
DepressionDemo <- data.frame(depression, treatment, cluster,
  age, anxiety, duration, depression_bin)


