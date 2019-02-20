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
### code chunk number 5: glmertree.Rnw:84-86
###################################################
data("MHserviceDemo", package = "glmertree")
summary(MHserviceDemo)


###################################################
### code chunk number 6: glmertree.Rnw:93-95
###################################################
lmmt <- lmertree(outcome ~ 1 | cluster_id | age + gender + emotional + 
                 autism + impact + conduct, data = MHserviceDemo)


###################################################
### code chunk number 7: glmertree.Rnw:102-105 (eval = FALSE)
###################################################
## glmmt <- glmertree(factor(outcome > 0) ~ 1 | cluster_id | age + gender + 
##                   emotional + autism + impact + conduct, 
##                 data = MHserviceDemo, family = "binomial")


###################################################
### code chunk number 8: glmertree.Rnw:110-111 (eval = FALSE)
###################################################
## plot(lmmt)


###################################################
### code chunk number 9: glmertree.Rnw:117-118
###################################################
plot(lmmt, which = "tree")


###################################################
### code chunk number 10: glmertree.Rnw:126-127
###################################################
plot(lmmt, which = "ranef")


###################################################
### code chunk number 11: glmertree.Rnw:142-145 (eval = FALSE)
###################################################
## print(lmmt)
## coef(lmmt)
## ranef(lmmt)


###################################################
### code chunk number 12: glmertree.Rnw:150-151
###################################################
predict(lmmt, newdata = MHserviceDemo[1:10,])


###################################################
### code chunk number 13: glmertree.Rnw:156-157
###################################################
predict(lmmt, newdata = MHserviceDemo[1:10, -7], re.form = NA)


###################################################
### code chunk number 14: glmertree.Rnw:165-169 (eval = FALSE)
###################################################
## resids <- residuals(lmmt)
## preds <- predict(lmmt)
## plot(MHserviceDemo$cluster_id, resids)
## scatter.smooth(preds, resids)


###################################################
### code chunk number 15: glmertree.Rnw:175-178
###################################################
resids <- residuals(lmmt)
preds <- predict(lmmt)
plot(MHserviceDemo$cluster_id, resids, xlab = "Cluster", ylab = "Residuals")


###################################################
### code chunk number 16: glmertree.Rnw:180-181
###################################################
scatter.smooth(preds, resids, xlab = "Predicted values", ylab = "Residuals")


###################################################
### code chunk number 17: glmertree.Rnw:194-196
###################################################
data("DepressionDemo", package = "glmertree")
summary(DepressionDemo)


###################################################
### code chunk number 18: glmertree.Rnw:203-205
###################################################
lmm_tree <- lmertree(depression ~ treatment | cluster |
  age + duration + anxiety, data = DepressionDemo)


###################################################
### code chunk number 19: glmertree.Rnw:212-214
###################################################
glmm_tree <- glmertree(depression_bin ~ treatment | cluster |
  age + duration + anxiety, data = DepressionDemo, family = binomial)


###################################################
### code chunk number 20: glmertree.Rnw:219-220 (eval = FALSE)
###################################################
## plot(lmm_tree)


###################################################
### code chunk number 21: glmertree.Rnw:226-227
###################################################
plot(lmm_tree, which = "tree")


###################################################
### code chunk number 22: glmertree.Rnw:235-236
###################################################
plot(lmm_tree, which = "ranef")


###################################################
### code chunk number 23: glmertree.Rnw:251-254
###################################################
print(lmm_tree)
coef(lmm_tree)
ranef(lmm_tree)


###################################################
### code chunk number 24: glmertree.Rnw:259-260
###################################################
predict(lmm_tree, newdata = DepressionDemo[1:7,])


###################################################
### code chunk number 25: glmertree.Rnw:265-266
###################################################
predict(lmm_tree, newdata = DepressionDemo[1:7, -3], re.form = NA)


###################################################
### code chunk number 26: glmertree.Rnw:274-278 (eval = FALSE)
###################################################
## resids <- residuals(lmm_tree)
## preds <- predict(lmm_tree)
## plot(factor(DepressionDemo$cluster), resids)
## scatter.smooth(preds, resids)


###################################################
### code chunk number 27: glmertree.Rnw:284-287
###################################################
resids <- residuals(lmm_tree)
preds <- predict(lmm_tree)
plot(factor(DepressionDemo$cluster), resids, xlab = "Cluster", ylab = "Residuals")


###################################################
### code chunk number 28: glmertree.Rnw:289-290
###################################################
scatter.smooth(preds, resids, xlab = "Predicted values", ylab = "Residuals")


###################################################
### code chunk number 29: appendix1
###################################################
set.seed(123)
treatment <- rbinom(n = 150, size = 1, prob = .5)
duration <- round(rnorm(150, mean = 7, sd = 3))
anxiety <- round(rnorm(150, mean = 10, sd = 3))
age <- round(rnorm(150, mean = 45, sd = 10))
error <- rnorm(150, 0, 2)


###################################################
### code chunk number 30: appendix2
###################################################
cluster <- error + rnorm(150, 0, 6)
rand_int <- sort(rep(rnorm(10, 0, 1), each = 15))
rand_int[order(cluster)] <- rand_int 
error <- error - rand_int
cluster[order(cluster)] <- rep(1:10, each = 15)


###################################################
### code chunk number 31: appendix3
###################################################
node3t1 <- ifelse(duration <= 8 & anxiety <= 10 & treatment == 0, -2, 0)
node3t2 <- ifelse(duration <= 8 & anxiety <= 10 & treatment == 1, 2, 0)
node5t1 <- ifelse(duration > 8 & treatment == 0, 2.5, 0)
node5t2 <- ifelse(duration > 8 & treatment == 1, -2.5, 0)


###################################################
### code chunk number 32: appendix4
###################################################
depression <- round(9 + node3t1 + node3t2 + node5t1 + node5t2 +
  .4 * treatment + error + rand_int)
depression_bin <- factor(as.numeric(depression > 9))


###################################################
### code chunk number 33: appendix5
###################################################
treatment <- factor(treatment, labels = c("Treatment 1", "Treatment 2"))
DepressionDemo <- data.frame(depression, treatment, cluster,
  age, anxiety, duration, depression_bin)


