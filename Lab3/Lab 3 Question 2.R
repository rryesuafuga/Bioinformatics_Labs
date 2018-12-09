

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


# is.binary.tree(tree_phylo)
# tree2 <- multi2di(tree_phylo)
# size2 <- pic(carni70_2$size, tree2)
# range2 <- pic(carni70_2$range, tree2)
# z <- lm(size2 ~ range2 - 1)    # the "-1" forces line through origin
# summary(z)
# correlationz <- sqrt(summary(z)$r.squared)
# ## no correlation
# bm.tree = fastBM(tree_phylo, a=0, sig2=1.0, internal = TRUE)
# phenogram(tree_phylo, bm.tree, spread.labels = TRUE)
# vcv.phylo(tree_phylo, cor=TRUE)

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



### 2.2.4 traits evolve as a bivariate Ornstein-Uhlenbeck process
mvMORPH::mvOU(tree_phylo, data=carni70$tab, model = c("OU1"), 
               diagnostic = TRUE, echo = TRUE)


### 2.2.5
mvslouchModel(tree_ouch, kY=1, 
              data = as.matrix(cbind(carni70$tab$range,carni70$tab$size)))
