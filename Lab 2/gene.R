library(ape)
library(seqinr)
library(phangorn)
library(markovchain)
library(msa)

#0
library(ape)
## Gene bank accession numbers taken from http://www.jcsantosresearch.org/Class_2014_Spring_Comparative/pdf/week_2/Jan_13_15_2015_GenBank_part_2.pdf
lizards_accession_numbers <- c("JF806202", "HM161150", "FJ356743", "JF806205", 
                               "JQ073190", "GU457971", "FJ356741", "JF806207",
                               "JF806210", "AY662592", "AY662591", "FJ356748",       
                               "JN112660", "AY662594", "JN112661", "HQ876437", 
                               "HQ876434", "AY662590", "FJ356740", "JF806214", 
                               "JQ073188", "FJ356749", "JQ073189", "JF806216", 
                               "AY662598", "JN112653", "JF806204", "FJ356747", 
                               "FJ356744", "HQ876440", "JN112651", "JF806215",
                               "JF806209") 
lizards_sequences<-ape::read.GenBank(lizards_accession_numbers)
print(lizards_sequences)
ape::write.dna(lizards_sequences, file ="lizard_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = "", colw = 10)

#1.1
clean <- function(template_gene){
  nucleotide <- c("a", "c", "g", "t")
  for (i in 1:length(template_gene)) {
    #Remove the " " that created when reading a file
    template_gene[[i]] <- template_gene[[i]][template_gene[[i]]!= " "]
    #Remove the character that not nucleotide (eg: name of species...)
    template_gene[[i]] <- template_gene[[i]][match(template_gene[[i]], nucleotide)]
  }
  return(template_gene)
}

# Read and clean 
lizards_sequences <-  read.fasta("lizard_seqs.fasta")
lizards_sequences <- clean(lizards_sequences)

#Simulate AI gen function
simulate_gene <- function(template_gene)
{
  ai_gene <- list()
  gene_num <- length(template_gene)
  nucleotide <- c("a", "t", "g", "c")
  
  #Scan all gene and get some information
  for (i in 1:gene_num) {
    
    template_sequence <- template_gene[[i]]
    #get leng and base compotision of gene
    this_leng <- length(template_sequence)
    this_compotision = seqinr::count(template_sequence,1)/this_leng
    
    #generate a new sequence base on sample function
    this_sequence <- sample(nucleotide, size=this_leng ,prob = this_compotision, replace = TRUE)
    #print(this_sequence)
    
    #add to list
    ai_gene[i] <- list(this_sequence)
  }
  
  #write to a file
  ape::write.dna(ai_gene, file ="AI_gene.fasta", format = "fasta", colsep =" ")
  return("Created an AI gene and saved in file: AI_gene.fasta")
}


simulate_gene(lizards_sequences)
ai_gene_1.1 <- read.fasta("AI_gene.fasta")
ai_gene_1.1 <- clean(ai_gene_1.1)


#1.2
tree <- rtree(length(lizards_sequences))
plot(tree)

#the rates matrix
rates <- matrix(0, ncol = 4, nrow = 4)
rownames(rates) <- c("a", "c", "g", "t")
colnames(rates) <- c("a", "c", "g", "t")

#fill value to rates matrix
for (i in 1:4) {
  rate = runif(3, 0.22, 0.28)
  for (j in 1:4) {
    if (j==4) 
    {
      rates[i,j]= 1- sum(rate)
    }else rates[i,j] = rate[j]
  }
}

#print the rates
rates
#create the ai_gene 2
ai_gene_1.2 <- phangorn::simSeq(tree, l = 1000, Q=rates , type = "DNA")

#rename 
for (i in 1:length(ai_gene_1.2)){
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == 1] = "a"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == 2] = "c"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == 3] = "g"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == 4] = "t"
}

ape::write.dna(ai_gene_1.2, file ="AI_gene2.fasta", format = "fasta", colsep ="")

#2.1
ai_gene_1.2 <- read.fasta("AI_gene2.fasta")

for (i in 1:length(lizards_sequences)) {
  cat(paste("For the sequences number: ", i , "\n"))
  print("The composition of lizards dataset:")
  print(round(count(lizards_sequences[[i]],2)/length(lizards_sequences[[i]]), 4))
  
  print("The composition of AI_gene1 dataset:")
  print(round(count(ai_gene_1.1[[i]],2)/length(ai_gene_1.1[[i]]), 4))
  
  print("The composition of Ai_gene2 dataset:")
  print(round(count(ai_gene_1.2[[i]],2)/length(ai_gene_1.2[[i]]), 4))
  
  cat("\n")
}


#2.2
library(markovchain)
markovchainFit(lizards_sequences)
markovchainFit(ai_gene_1.1)
markovchainFit(ai_gene_1.2)

#The fitHigherOrder function work only with a list, not list of list
#So, my idea is take some random sequences, see the order and then make the conclusion 
fitHigherOrder(lizards_sequences[[1]])


#2.3
library(msa)
#vignette("msa")
#example("msa")
#Remember that the fasta file in here should not have any space (colsep ="" when write)
mySequenceFile <- system.file("examples", "lizard_seqs.fasta", package="msa")
mySequences <- readAAStringSet(mySequenceFile)
mySequences


real_alignseq<- msaConvert(real_align, type="seqinr::alignment") 
dist_real <- as.matrix(dist.alignment(real_alignseq, "identity"))

plot_ly(x=colnames(dist_real), y=rownames(dist_real), 
        z=dist_real, type="heatmap", colors =colorRamp(c("yellow", "red")))%>% 
  layout(title = "Heatmap of the Lizard sequences")

#Do the same thing with 2 other sequences. 
#You can see that the value in heatmap of AI sequences is quite low
# It means that the Ai gene has no connections (because it randomly created) while the original gene is highly connected


#3.1
treeUPGMA <- upgma(dist_real)
treeNJ <- upgma(dist_real)
layout(matrix(c(1,2), 2, 1), height=c(1,2))
par(mar = c(0,0,2,0)+ 0.1)
plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")

#haven't done the last question
