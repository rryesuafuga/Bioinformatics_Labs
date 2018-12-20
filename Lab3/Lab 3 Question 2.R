

# ------------------ Question 2 Lab 3 ----------------------


# Question 2

## 2.1
library(ade4)
data(carni70)
View(carni70)
str(carni70)
names(carni70)

carni70_1 <- carni70[[1]]
carni70_2 <- carni70[[2]]

library(tidyverse)

hist(carni70_2$size, xlab="size trait")
hist(carni70_2$range, xlab="range trait")
boxplot(carni70_2$size, main="Size")
boxplot(carni70_2$range, main="Range")
ggplot()+aes(x=names(carni70_2)[1], y=carni70_2$size)+
  geom_boxplot()+
  geom_boxplot(aes(x=names(carni70_2[2]), y=carni70_2$range))


### generally, they have small sizes
### Ursus_arctos has the biggest size


rownames(carni70_2[which.max(carni70_2[,"size"]), ])
rownames(carni70_2[which.min(carni70_2[,"size"]), ])
rownames(carni70_2[which.max(carni70_2[,"range"]), ])
rownames(carni70_2[which.min(carni70_2[,"range"]), ])



library(ape)
tree_phylo <- ape::read.tree(text=carni70_1)
plot(tree_phylo)


## 2.2

### 2.2.1 Both traits evolve as independent Brownian motions
library(mvMORPH)
mvMORPH::mvBM(tree_phylo, carni70$tab, model="BM1", 
             param = list(constraint = "diagonal"))


### 2.2.2 The traits evolve as a correlated Brownian motion
library(ouch)
tree_ouch <- ouch::ape2ouch(tree_phylo, branch.lengths = tree_phylo$edge.length)
library(mvSLOUCH)
bm_cor <- mvSLOUCH::BrownianMotionModel(tree_ouch, data = as.matrix(carni70$tab))


### 2.2.3 independent Ornstein Uhlenbeck processes
mvMORPH::mvOU(tree_phylo, data=carni70$tab$size, model = c("OU1"), 
              diagnostic = FALSE, echo = TRUE)
mvMORPH::mvOU(tree_phylo, data=carni70$tab$range, model = c("OU1"), 
              diagnostic = FALSE, echo = TRUE)

# mvslouchModel(tree_ouch, 
#               data = as.matrix(cbind(carni70$tab$size,carni70$tab$range)))



### 2.2.4 traits evolve as a bivariate Ornstein-Uhlenbeck process
mvMORPH::mvOU(tree_phylo, data=carni70$tab, model = c("OU1"), 
               diagnostic = TRUE, echo = TRUE)


### 2.2.5 size evolves as a Brownian motion and range as an Ornstein Uhlenbeck
# mvslouchModel(tree_ouch, kY=1, 
#               data = as.matrix(cbind(carni70$tab$range,carni70$tab$size)))



model1=c("Both traits evolve as independent Brownian motions", 1186.475, 1186.771)

model2=c("Both traits as a correlated Brownian motion", 1323.473, 1323.921)

model3=c("both traits as independent Ornstein Uhlenbeck processes", 1158.5068, 1159.2341)

model4=c("both traits as bivariate Ornstein-Uhlenbeck processes", 
         1161.213, 1162.312) 

model5=c("size as Brownian motion and range as Ornstein Uhlenbeck", 
         1295.603, 1296.235)

comparison_matrix <- rbind(model1, model2, model3, model4, model5)
comparison_matrix <- as.data.frame(comparison_matrix)
comparison_matrix[,2] <- as.numeric(as.character(comparison_matrix[,2]))
comparison_matrix[,3] <- as.numeric(as.character(comparison_matrix[,3]))
colnames(comparison_matrix) <- c("model", "AIC", "AIC.c")
rownames(comparison_matrix) <- NULL

print(comparison_matrix[,1:2])
