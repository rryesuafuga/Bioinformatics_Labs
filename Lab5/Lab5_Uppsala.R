
##-----------------------------Lab5------------------------------

# 1
#install.packages("devtools")
#library(devtools)
#install_github("mategarb/R.ROSETTA")
library(R.ROSETTA)

data("autcon")
View(autcon)


# 2
## 35 features
# number of objects in each class
summary(autcon$decision)


# 3
autconDefault = rosetta(autcon)
table_rule <- autconDefault$main

autconDefault$quality

## 3a 10 cross validations
## 3b Johnson
## 3c Equal Frequency 3 bins
## 3d  0.821818 mean ...., meadian

## 3e 
table_rule[1:3,]
nrow(table_rule) ##number of rules

nrow(table_rule[which(table_rule$PVAL<= 0.05),])


# 4
# saveLineByLine(table_rule, "table_rule.txt", 
#               discrete=FALSE, filterByPval=FALSE, pval=0.05)


# 5
## big red line strong connection
## big bands
## many connections
## check table of rules
## rules CUT_COND (decision)
## significant ---> accuracy, support, number of rules