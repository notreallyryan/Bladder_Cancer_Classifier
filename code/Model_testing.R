library(dplyr)
library(survival)
library(caret)
library(ggplot2)
library(survminer)
library(pROC)

#loading in data
train_data <- readRDS("data/train.rds")
test_data <- readRDS("data/test.rds")
M_data <- readRDS("data/Microarray.rds")
final_gene_list <- read.csv("data/final_gene_list.csv")$x

#keeping only the genes in final_gene_list
train_data$exprs <- train_data$exprs[,colnames(train_data$exprs) %in% final_gene_list]
test_data$exprs <- test_data$exprs[,colnames(test_data$exprs) %in% final_gene_list]
M_data$exprs <- M_data$exprs[,colnames(M_data$exprs) %in% final_gene_list]

#To ensure model learns shape instead of values, we Z score the expression values.
#This way model is applicable to the differently scaled microarray data. 
train_mean <- apply(train_data$exprs, 2, mean)
train_sd <- apply(train_data$exprs, 2, sd)

train_data$exprs <- scale(train_data$exprs)
test_data$exprs <-scale(test_data$exprs, train_mean, train_sd)
M_data$exprs <- scale(M_data$exprs)

#making final X and Y dataframes
to_drop <- c("Progression", "PFS_time.", "Recurrence", "RFS_time", "FUtime_days.")

X <- train_data %>%
  select(-all_of(c("UROMOL.ID", to_drop)))
X <- cbind(X, X$exprs)
X <- subset(X, select = -c(exprs))

test_X <- test_data %>%
  select(-all_of(c("UROMOL.ID", to_drop)))
test_X <- cbind(test_X, test_X$exprs)
test_X <- subset(test_X, select = -c(exprs))

M_X <- M_data %>%
  select(-all_of(c("knowles_ID", to_drop)))
M_X <- cbind(M_X, M_X$exprs)
M_X <- subset(M_X, select = -c(exprs))

#Using a Support Vector Machine approach with polynomial kernel. 
set.seed(67) #DO NOT CHANGE THIS.
control <- trainControl(method = "cv",
                        number = 10,
                        classProbs = TRUE)

svmPoly <- train(label ~ .,
                 data = X,
                 method = "svmPoly",
                 trControl = control)

#adding prediction labels back to main dataframe.
train_data$preds <- predict(svmPoly, X)
test_data$preds <- predict(svmPoly, test_X)
M_data$preds <- predict(svmPoly, M_X)

#Getting cm objects
train_cm <- confusionMatrix(train_data$preds, train_data$label)
test_cm <- confusionMatrix(test_data$preds, test_data$label)
M_cm <- confusionMatrix(M_data$preds, M_data$label)

#saving confusion matrix
write.table(train_cm$table, "results/train_cm.csv")
write.table(test_cm$table, "results/test_cm.csv")
write.table(M_cm$table, "results/Microarray_cm.csv")

#Accuracy, Precision, Recall, and F1 can be obtained from the confusion Matrix objects.
Accuracy <- c(train_cm$overall[['Accuracy']] , test_cm$overall[['Accuracy']], M_cm$overall[['Accuracy']])
Precision <- c(train_cm$byClass[['Precision']] , test_cm$byClass[['Precision']], M_cm$byClass[['Precision']])
Recall <- c(train_cm$byClass[['Recall']] , test_cm$byClass[['Recall']], M_cm$byClass[['Recall']])
F1 <- c(train_cm$byClass[['F1']] , test_cm$byClass[['F1']], M_cm$byClass[['F1']])

#calculating AUC values
train_AUC <- auc(train_data$label, predict(svmPoly, X, type = "prob")[,1])
test_AUC <- auc(test_data$label, predict(svmPoly, test_X, type = "prob")[,1])
M_AUC <- auc(M_data$label, predict(svmPoly, M_X, type = "prob")[,1])
AUC <- c(train_AUC, test_AUC, M_AUC)

#making dataframe for saving
score_df <- data.frame(Accuracy, Precision, Recall, F1, AUC)
rownames(score_df) <- c("Train RNA-seq", "Test RNA-seq", "Microarray")
write.csv(score_df,
          file = "results/scores.csv")

#Plotting Survival curves by predicted class
#saving as pdf because ggsave will not work with ggsurvplot for some reason.
train_fit <- surv_fit(Surv(RFS_time, Recurrence) ~ preds, data = train_data)
train_plot <- ggsurvplot(train_fit,
           title = "Training Data Survival Plot",
           xlab = "Time in Months")
pdf("results/train_survplot.pdf")
print(train_plot, newpage = FALSE)
dev.off()

test_fit <- surv_fit(Surv(RFS_time, Recurrence) ~ preds, data = test_data)
test_plot <- ggsurvplot(test_fit, 
           title = "Testing Data Survival Plot",
           xlab = "Time in Months")
pdf("results/test_survplot.pdf")
print(test_plot, newpage = FALSE)
dev.off()

M_fit <- surv_fit(Surv(RFS_time, Recurrence) ~ preds, data = M_data)
M_plot <- ggsurvplot(M_fit, 
           title = "Microarray Data Survival Plot",
           xlab = "Time in Months")
pdf("results/M_survplot.pdf")
print(M_plot, newpage = FALSE)
dev.off()

#clear environment
rm(list=ls())