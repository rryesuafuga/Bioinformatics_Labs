#install.packages("devtools")
#library(devtools)
#install_github("mategarb/R.ROSETTA")
#library(R.ROSETTA)
View(autcon)
data("autcon")
decision <- table(autcon$decision)
autconDefault = rosetta(autcon)
autconDefault$main
table(autconDefault$main$DECISION)
autconDefault$quality
#CV = 10
#reducer = "Johnson
#discreteMethod = "EqualFrequency"
#discreteParam = 3
#Mean accurancy = 0.821818
autconDefault$main[1:3,]
nrow(autconDefault$main[which(autconDefault$main$PVAL <0.05),])
saveLineByLine(autconDefault$main, "rules.txt")
save.image(file = "Rosetta_Lab5.RData")
load("Rosetta_Lab5.RData")


