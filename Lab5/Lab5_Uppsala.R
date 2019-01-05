
##-----------------------------Lab5------------------------------

# 1
#install.packages("devtools")
#library(devtools)
#install_github("mategarb/R.ROSETTA")
library(R.ROSETTA)

data("autcon")
View(autcon)


# 2
ncol(autcon)-1
## 35 features

# number of objects in each class
summary(autcon$decision)


# 3
autconDefault = rosetta(autcon)
table_rule <- autconDefault$main

autconDefault$quality

## 3a 10 cross validations (from Rosetta documentation)
## According to Wikipedia, cross validation is any of the many model validation
## techniques for assessing how the results of a statistical analysis will
## generalize to an independent data set.

## 3b Johnson (from Rosetta documentation)
## A reduct is a minimal subset of attributes such that it preserves
## indiscernibility between objects

## 3c Equal Frequency 3 bins (from Rosetta documentation)
## Discritization is when you classify values obtained from measuring a 
## phenomenon into distinct states. For example, recording temperature values 
## ranging from 34.5 to 37 degrees centigrade as "Fever", "Normal" and
## "Hyperthermia.

## 3d  0.821818 mean, 0.8 median

## 3e 
table_rule[1:3,]
head(table_rule, 10)
nrow(table_rule) ##number of rules

nrow(table_rule[which(table_rule$PVAL<= 0.05),])
head(table_rule[which(table_rule$PVAL<= 0.05),], 3)
y <- table_rule[which(table_rule$PVAL<= 0.05),]
y <- y[order(y$PVAL),]
head(y, 3)

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