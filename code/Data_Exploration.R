library(ggplot2)
library(dplyr)

training_data <-readRDS("data/UROMOL_TaLG.teachingcohort.rds")

#Plotting PFS with EUA groupings
ggplot(data = training_data, aes(x=UROMOL2021.classification, y=PFS_time.)) + geom_boxplot()

#Plotting RFS with EUA groupings
ggplot(data = training_data, aes(x=UROMOL2021.classification, y=RFS_time)) + geom_boxplot()

#Plotting RFS with EUA groupings
ggplot(data = training_data, aes(x=EAU.risk, y=FUtime_days.)) + geom_boxplot()

#Well that wasn't very helpful.

#how many NA values are there in the Recurrence section?
length(which(is.na(training_data$Recurrence)))
length(which(is.na(training_data$RFS_time)))

#What about progression?
length(which(is.na(training_data$Progression)))
length(which(is.na(training_data$PFS_time.)))

#how many people have progression?
length(which(training_data$Progression == 1))

#filling in missing progression values with zero
training_data$Progression[is.na(training_data$Progression)] = 0

#removing those that do not have recurrence data:
training_data <- training_data[(!is.na(training_data$Recurrence)),]

#plotting distribution of RFS times
RFS_data <- training_data[training_data$Recurrence == 1, ]
ggplot(data=RFS_data, aes(x=RFS_time)) + geom_histogram()
#a large amount of patients experience recurrence within 21 days.

#making labels:
#rules:
# progression at any time or recurrence within 21 days = positive case
# no recurrence or recurrence occurs outside of 21 days = negative case

training_data$group <- with(training_data, ifelse(Progression | (Recurrence & RFS_time <= 42), "HIGH RISK", "LOW RISK"))
training_data <- training_data %>% select(group, everything())

#plotting number of LOW and HIGH risk cases
ggplot(data=training_data, aes(x=group)) + geom_bar()


training_data %>% count(group)
