library(dplyr)

RNAseq <- readRDS("data/UROMOL_TaLG.teachingcohort.rds")
Microarray <- readRDS("data/knowles_matched_TaLG_final.rds")

# DATA CLEANING
# - If Recurrence is NA, delete row. 
# - If Recurrence does occur, but RFS is NA, assume median of class
# - If Recurrence does not occur but RFS is NA, set RFS to 0

# As these values are being predicted, they won't be input into the model so 
# there is no reason to save them.

#train data cleaning
RNAseq <- RNAseq[(!is.na(RNAseq$Recurrence)),]
RNAseq$RFS_time[is.na(RNAseq$RFS_time) & RNAseq$Recurrence == 1] = 
  median(RNAseq$RFS_time, na.rm = TRUE)
RNAseq$RFS_time[is.na(RNAseq$RFS_time) & RNAseq$Recurrence == 0] = 0

#test data cleaning
Microarray <- Microarray[(!is.na(Microarray$Recurrence)),]
Microarray$RFS_time[is.na(Microarray$RFS_time) & Microarray$Recurrence == 1] = 
  median(Microarray$RFS_time, na.rm = TRUE)
Microarray$RFS_time[is.na(Microarray$RFS_time) & Microarray$Recurrence == 0] = 0

# ASSIGNING LABELS
# Based on the data exploration we will use the following grouping:
# HIGH RISK: Recurrence within 21 days
# LOW RISK: no Recurrence or Recurrence beyond 21 days.

RNA_labels <- with(RNAseq, 
                     ifelse((Recurrence & RFS_time <= 36), 
                            "HIGH", "LOW"))

M_labels <- with(Microarray,
                    ifelse(Recurrence & RFS_time <= 36, 
                           "HIGH", "LOW"))

RNAseq$label <- as.factor(RNA_labels)
Microarray$label <- as.factor(M_labels)

#KEEPING ONLY THE SHARED GENES
#no point in making the model if it is highly focused on features only found in
#the training dataset.
common_genes <- intersect(colnames(RNAseq$exprs), colnames(Microarray$exprs))
RNAseq$exprs <- RNAseq$exprs[,colnames(RNAseq$exprs) %in% common_genes]
Microarray$exprs <- Microarray$exprs[,colnames(Microarray$exprs) %in% common_genes]

#A lot of demographic/clinical variables are NA - we need to fill these
#assume "null" case scenario
RNAseq$Smoking[is.na(RNAseq$Smoking)] <- "Never"
RNAseq$Tumor.size[is.na(RNAseq$Tumor.size)] <- "< 3 cm"
RNAseq$Age[is.na(RNAseq$Age)] <- median(RNAseq$Age, na.rm = TRUE)

#notably the microarray dataset is missing smoking, size, and incident information
#in each case, assume the best scenario - never smoker, <3cm, and incident tumor.
Microarray$Smoking <- "Never"
Microarray$Tumor.size <- "< 3 cm"
Microarray$Incident.tumor <- "Yes"

#remove features that will not be used.(Stage, Grade, EAU, UROMOL classes)
RNAseq <- subset(RNAseq, select = -c(Tumor.stage, Tumor.grade, EAU.risk, UROMOL2021.classification))
Microarray <- subset(Microarray, select = -c(Tumor.stage, Tumor.grade, UROMOL2021.classification))

#convert all character columns to factors. This will help down the line.
RNAseq <- RNAseq %>%
  mutate(across(where(is.character), as.factor))
Microarray <- Microarray %>%
  mutate(across(where(is.character), as.factor))

#convert BCG to factor, and age to a numeric
RNAseq$BCG <- as.factor(RNAseq$BCG)
Microarray$BCG <- as.factor(Microarray$BCG)

RNAseq$Age <- as.numeric(RNAseq$Age)
Microarray$Age <- as.numeric(Microarray$Age)

#now we need to randomly split the training data into training and test
#This is especially important, our testing set is purely microarray data. 
#we also want to see how well the model does on unseen RNA-seq data!
set.seed(67)
split <- sample(c(TRUE, FALSE), nrow(RNAseq), replace=TRUE, prob=c(0.7,0.3))
train <- RNAseq[split, ]
test <- RNAseq[!split, ]

#SAVE DATASETS
saveRDS(train, "data/train.rds")
saveRDS(test, "data/test.rds")
saveRDS(Microarray, "data/Microarray.rds")

#clear environment
rm(list=ls())