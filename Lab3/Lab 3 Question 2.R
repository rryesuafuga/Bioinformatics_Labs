

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

ggplot(data = carni70_2)+
  layer(mapping = aes(x=size, y=range), geom = "point", stat = "identity",
        position = "identity")

ggplot(data = carni70_2)+
  layer(mapping = aes(x=rownames(carni70_2), y=size), geom = "bar", 
        stat = "identity", position = "identity")+coord_flip()

ggplot(data = carni70_2)+
  layer(mapping = aes(x=rownames(carni70_2), y=range), geom = "bar", 
        stat = "identity", position = "identity")

### generally, they have small sizes
### Ursus_arctos has the biggest size
### package meant for multivariate data analysis.

library(ape)
tree <- read.tree(text=carni70_1)
plot(tree)



## 2.2.1


is.binary.tree(tree)
tree2 <- multi2di(tree)
size2 <- pic(carni70_2$size, tree2)
range2 <- pic(carni70_2$range, tree2)
z <- lm(size2 ~ range2 - 1)    # the "-1" forces line through origin
summary(z)
correlationz <- sqrt(summary(z)$r.squared)
## no correlation
bm.tree = fastBM(tree, a=0, sig2=1.0, internal = TRUE)
phenogram(tree, bm.tree, spread.labels = TRUE)
vcv.phylo(tree, cor=TRUE)


### 2.2.1 Both traits evolve as independent Brownian motions
library(mvMORPH)
# Fitting the models
bm_i <- mvBM(tree, carni70$tab, model="BM1", param = list(constraint = "diagonal"))


### 2.2.2 The traits evolve as a correlated Brownian motion
library(ouch)
tree_ouch <- ouch::ape2ouch(tree, branch.lengths = tree$edge.length)
library(mvSLOUCH)
bm_cor <- mvSLOUCH::BrownianMotionModel(tree_ouch, data = as.matrix(carni70$tab))


### 2.2.3 independent Ornstein Uhlenbeck processes
mvOU(tree, data=carni70$tab$size, model = c("OU1"), diagnostic = FALSE, echo = TRUE)
mvOU(tree, data=carni70$tab$range, model = c("OU1"), diagnostic = FALSE, echo = TRUE)



### 2.2.4 traits evolve as a bivariate Ornstein{Uhlenbeck process
mvOU(tree, data=carni70$tab, model = c("OU1"), diagnostic = TRUE, echo = TRUE)


### 2.2.5
