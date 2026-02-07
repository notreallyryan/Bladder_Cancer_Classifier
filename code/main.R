#literally just runs the entire project

#clean data
source("code/Data_Processing.R")

#Filter for most important genomic features
source("code/Feature_Selection.R")

#Train Model and make figures
source("code/Model_testing.R")