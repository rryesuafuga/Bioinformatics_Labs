

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



## 2.2

### 2.2a Both traits evolve as independent Brownian motions
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

library(mvMORPH)
# Simulate a trait evolving by brownian motion on the tree
#trait<-rTraitCont(tree)
# Fitting the models
#sizebm <- mvBM(tree, carni70_2$size, model="BM1", method="pic")
#rangebm <- mvBM(tree, carni70_2$range, model="BM1", method="pic")


library(ouch)
# treeo <- convert(tree)
s <- ape2ouch(tree, scale = TRUE, branch.lengths = tree$edge.length)
bm.size <- ouch::brown(carni70_1$size, s)


### 2.2b The traits evolve as a correlated Brownian motion



### 2.2c independent Ornstein Uhlenbeck processes





