library(ape)
library(seqinr)

#1.1
lizards_sequences <-  read.fasta("lizard_seqs.fasta")

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
    this_sequence <- sample(nucleotide, prob = this_compotision, replace = TRUE)
    print(this_sequence)
    
    #add to list
    ai_gene[i] <- list(this_sequence)
  }
  
  #write to a file
  ape::write.dna(ai_gene, file ="AI_gene.fasta", format = "fasta", colsep =" ")
  return("Created an AI gene and saved in file: AI_gene.fasta")
}


simulate_gene(lizards_sequences)
ai_gene_1.1 <- read.fasta("AI_gene.fasta")

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
ai_gene_1.2 <- simSeq(tree, l = 1000, Q=rates , type = "DNA")

#rename 
for (i in 1:length(ai_gene_1.2)){
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "1"] = "a"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "2"] = "c"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "3"] = "g"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "4"] = "t"
}

ape::write.dna(ai_gene_1.2, file ="AI_gene2.fasta", format = "fasta", colsep =" ")

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


