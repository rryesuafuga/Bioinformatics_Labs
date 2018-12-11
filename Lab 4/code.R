#For install the requirement package for this lab
install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")
BiocManager::install("simpleaffy", version = "3.8")
BiocManager::install("affyPLM", version = "3.8")
BiocManager::install("hgu133plus2.db", version = "3.8")
install.packages("RColorBrewer")


#Start running the code
##################

library(GEOquery)
#Pay attention at the working directory before running this line
#download data
x = getGEOSuppFiles("GSE20986")
#show list of downloaded file
x
#unzip to "data" folder
untar("GSE20986/GSE20986_RAW.tar", exdir = "data")
#cels is the list of gz file
cels = list.files("data/", pattern = "[gz]")
#unzip file
sapply(paste("data", cels, sep = "/"), gunzip)
phenodata = matrix(rep(list.files("data"), 2), ncol =2)
class(phenodata)

phenodata <- as.data.frame(phenodata)
colnames(phenodata) <- c("Name", "FileName")
phenodata$Targets <- c("iris", 
                       "retina", 
                       "retina", 
                       "iris", 
                       "retina", 
                       "iris", 
                       "choroid", 
                       "choroid", 
                       "choroid", 
                       "huvec", 
                       "huvec", 
                       "huvec")
#Write the list of downloaded content to a file
write.table(phenodata, "data/phenodata.txt", quote = F, sep = "\t", row.names = F)

library(simpleaffy)
#Using read.affy function to read.. 
celfiles <- read.affy(covdesc = "phenodata.txt", path = "data")
boxplot(celfiles)

library(RColorBrewer)
cols = brewer.pal(8, "Set1")
eset <- exprs(celfiles)
samples <- celfiles$Targets
colnames(eset)

colnames(eset) <- samples
boxplot(celfiles, col = cols, las = 2) #las=2 make the axis labels horizontal
distance <- dist(t(eset), method = "maximum")
clusters <- hclust(distance)
plot(clusters)


require(simpleaffy)
require(affyPLM)
celfiles.gcrma = gcrma(celfiles)

par(mfrow=c(1,2))
boxplot(celfiles.gcrma, col = cols, las = 2, main = "Post-Normalization") 
boxplot(celfiles, col = cols, las = 2, main = "Pre-Normalization")

dev.off()

distance <- dist(t(exprs(celfiles.gcrma)), method = "maximum")
clusters <- hclust(distance)
plot(clusters)

library(limma)
phenodata

samples <- as.factor(samples)
design <- model.matrix(~0+samples)
colnames(design)

colnames(design) <- c("choroid", "huvec", "iris", "retina")
design

contrast.matrix = makeContrasts(
  huvec_choroid = huvec - choroid, 
  huvec_retina = huvec - retina, 
  huvec_iris = huvec - iris, 
  levels = design)

fit = lmFit(celfiles.gcrma, design)
huvec_fit <- contrasts.fit(fit, contrast.matrix)
huvec_ebay <- eBayes(huvec_fit)

library(hgu133plus2.db)
library(annotate)

probenames.list <- rownames(topTable(huvec_ebay, number = 100000))
getsymbols <- getSYMBOL(probenames.list, "hgu133plus2")
results <- topTable(huvec_ebay, number = 100000, coef = "huvec_choroid")
results <- cbind(results, getsymbols)
summary(results)

results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
table(results$threshold)

library(ggplot2)
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
  geom_point() + 
  scale_color_manual(values = c("black", "red", "green"), 
                     labels = c("Not Significant", "Upregulated", "Downregulated"), 
                     name = "Key/Legend")

volcano + 
  geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5), aes(x = logFC, y = -1*log10(adj.P.Val), colour = threshold, label = getsymbols)  )


###
#Question 2
iris <- eset[,1]
retina <- eset[,2]
choroid <- eset[,7]
huvec <- eset[,10]

plot(x=huvec ,y=iris,xlab="huvec",ylab="iris", main="Scatterplot of raw data") 
plot(x=huvec ,y=retina,xlab="huvec",ylab="retina", main="Scatterplot of raw data") 
plot(x=huvec ,y=choroid,xlab="huvec",ylab="choroid", main="Scatterplot of raw data") 

plot(x=log(huvec) ,y=log(iris),xlab="huvec",ylab="iris", main="Scatterplot of raw data") 
plot(x=log(huvec) ,y=log(retina),xlab="huvec",ylab="retina", main="Scatterplot of raw data") 
plot(x=log(huvec) ,y=log(choroid),xlab="huvec",ylab="choroid", main="Scatterplot of raw data") 
##log
#huvec_iris <- pairwise.comparison(celfiles,"Targets",c("huvec","iris"))
#huvec_retina <-pairwise.comparison(celfiles,"Targets",c("huvec","retina"))
#huvec_choroid <-pairwise.comparison(celfiles,"Targets",c("huvec","choroid"))
#plot(huvec_iris)
#plot(huvec_retina)
#plot(huvec_choroid)

## MA plot
iris_huvec <- eset[,c(1,10)]
retina_huvec <- eset[,c(2,10)]
chronoid_huvec <- eset[,c(7,10)]
library(affy)
#inspired by wiki

ma.plot( rowMeans(log2(iris_huvec)), log2(iris_huvec[, 1])-log2(iris_huvec[, 2]), cex=1 )
ma.plot( rowMeans(log2(retina_huvec)), log2(retina_huvec[, 1])-log2(retina_huvec[, 2]), cex=1 )
ma.plot( rowMeans(log2(chronoid_huvec)), log2(chronoid_huvec[, 1])-log2(chronoid_huvec[, 2]), cex=1 )


###################################### normalized #######################################

##?? how to normalize??
eset_norm <- normalize.quantiles(eset)

iris <- eset_norm[,1]
retina <- eset_norm[,2]
choroid <- eset_norm[,7]
huvec <- eset_norm[,10]

plot(x=huvec ,y=iris,xlab="huvec",ylab="iris", main="Scatterplot of raw data") 
plot(x=huvec ,y=retina,xlab="huvec",ylab="retina", main="Scatterplot of raw data") 
plot(x=huvec ,y=choroid,xlab="huvec",ylab="choroid", main="Scatterplot of raw data") 


plot(x=log(huvec) ,y=log(iris),xlab="huvec",ylab="iris", main="Scatterplot of raw data") 
plot(x=log(huvec) ,y=log(retina),xlab="huvec",ylab="retina", main="Scatterplot of raw data") 
plot(x=log(huvec) ,y=log(choroid),xlab="huvec",ylab="choroid", main="Scatterplot of raw data") 

#huvec_iris <- pairwise.comparison(celfiles.gcrma,"Targets",c("huvec","iris"))
#huvec_retina <-pairwise.comparison(celfiles.gcrma,"Targets",c("huvec","retina"))
#huvec_choroid <-pairwise.comparison(celfiles.gcrma,"Targets",c("huvec","choroid"))
#plot(huvec_iris)
#plot(huvec_retina)
#plot(huvec_choroid)




library(preprocessCore)

#do a quantile normalization
norm_iris_huvec <- normalize.quantiles(iris_huvec)
norm_retina_huvec <- normalize.quantiles(retina_huvec)
norm_chronoid_huvec <- normalize.quantiles(chronoid_huvec)



##normalized
ma.plot( rowMeans(log2(norm_iris_huvec)), log2(norm_iris_huvec[, 1])-log2(norm_iris_huvec[, 2]), cex=1 )
ma.plot( rowMeans(log2(norm_retina_huvec)), log2(norm_retina_huvec[, 1])-log2(norm_retina_huvec[, 2]), cex=1 ) 
ma.plot( rowMeans(log2(norm_chronoid_huvec)), log2(norm_chronoid_huvec[, 1])-log2(norm_chronoid_huvec[, 2]), cex=1 ) 


###############################

results <- topTable(huvec_ebay, number = 100000, coef = "huvec_retina")
results <- cbind(results, getsymbols)
summary(results)

results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
table(results$threshold)

library(ggplot2)
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
  geom_point() + 
  scale_color_manual(values = c("black", "red", "green"), 
                     labels = c("Not Significant", "Upregulated", "Downregulated"), 
                     name = "Key/Legend")

volcano + 
  geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5), aes(x = logFC, y = -1*log10(adj.P.Val), colour = threshold, label = getsymbols)  )


results <- topTable(huvec_ebay, number = 100000, coef = "huvec_iris")
results <- cbind(results, getsymbols)
summary(results)

results$threshold <- "1"
a <- subset(results, adj.P.Val < 0.05 & logFC > 5)
results[rownames(a), "threshold"] <- "2"
b <- subset(results, adj.P.Val < 0.05 & logFC < -5)
results[rownames(b), "threshold"] <- "3"
table(results$threshold)

library(ggplot2)
volcano <- ggplot(data = results, 
                  aes(x = logFC, y = -1*log10(adj.P.Val), 
                      colour = threshold, 
                      label = getsymbols))

volcano <- volcano + 
  geom_point() + 
  scale_color_manual(values = c("black", "red", "green"), 
                     labels = c("Not Significant", "Upregulated", "Downregulated"), 
                     name = "Key/Legend")

volcano + 
  geom_text(data = subset(results, logFC > 5 & -1*log10(adj.P.Val) > 5), aes(x = logFC, y = -1*log10(adj.P.Val), colour = threshold, label = getsymbols)  )
