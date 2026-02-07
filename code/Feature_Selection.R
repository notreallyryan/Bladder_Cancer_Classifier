library(dplyr)
library(ggplot2)
library(glmnet)
library(survival)

#loading training data ONLY.
train <- readRDS("data/train.rds")

#some parameter values
VAR_CUTOFF <- 0

#FEATURE SELECTION
#Let's first get rid of low variance genes.
vars <- apply(train$exprs, 2, var)
train$exprs <- train$exprs[, vars > VAR_CUTOFF, drop = FALSE]

#how many genes are important to survival? using LASSO regression to reduce 
#non important genes to zero.
y <- Surv(train$RFS_time, train$Recurrence)

#caret was being annoying so i just swapped to base glmnet
set.seed(67)
cvfit <- cv.glmnet(
  train$exprs, y,
  family = "cox",
  alpha = 1,
  nfolds = 10
)

#Use the best case found by cv, we can fit it to a LASSO model to find best genes
#related to survival.
fit <- glmnet(train$exprs, y, family = "cox", alpha = 1, lambda = cvfit$lambda.min)
coef_mat <- coef(fit)
survival_genes <- which(coef_mat != 0)
final_genes_list <- rownames(coef_mat)[survival_genes]

#saving the top genes to be used in the model separately
write.csv(final_genes_list,
          file = "data/final_gene_list.csv")

#clear environment
rm(list=ls())
