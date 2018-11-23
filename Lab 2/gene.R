library(ape)
library(seqinr)

#1.1
lizards_sequences <-  read.fasta("lizard_seqs.fasta")

simulate_gene <- function(template_gene){
  ai_gene <- list()
  gene_num <- length(template_gene)
  nucleotide <- c("a", "t", "g", "c")
  
  #Scan all gene and get some information
  for (i in 1:gene_num) {
    template_sequence <- template_gene[[i]]
    #get leng and base compotision of gene
    this_leng <- length(template_sequence)
    this_compotision = seqinr::count(template_sequence,1)/this_leng
    this_sequence <- c()
    
    # create a number of each nucleotide based on the template_sequence 
    for (j in 1:4) {
      letter <- sample(nucleotide[j], size = this_compotision[[j]]*this_leng, replace = TRUE)
      this_sequence <- c(this_sequence, letter)
    }
    #random order
    rand <- sample(length(this_sequence))
    this_sequence <- this_sequence[rand]
    
    #add to list
    ai_gene[i] <- list(this_sequence)
  }
  
  #write to a file
  ape::write.dna(ai_gene, file ="AI_gene.fasta", format = "fasta", colsep =" ")
  return("Created an AI gene and saved in file: AI_gene.fasta")
}

simuate_gene(lizards_sequences)

#1.2
tree <- rtree(length(lizards_sequences))
plot(tree)


rates <- matrix(0, ncol = 4, nrow = 4)
rownames(rates) <- c("a", "c", "g", "t")
colnames(rates) <- c("a", "c", "g", "t")

for (i in 1:4) {
  rate = runif(3, 0.2, 0.3)
  for (j in 1:4) {
    if (j==4) 
      {
      rates[i,j]= 1- sum(rate)
    }else rates[i,j] = rate[j]
  }
}

rates

ai_gene_1.2 = list()
for (i in 1:length(lizards_sequences)) {
  this_leng <- length(lizards_sequences[[i]])
  ai_sequence2 <- simSeq(tree, l = this_leng, Q=rates , type = "DNA")
  ai_gene_1.2[i] <- ai_sequence2
}

for (i in 1:length(ai_gene_1.2)){
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "1"] = "a"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "2"] = "c"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "3"] = "g"
  ai_gene_1.2[[i]][ai_gene_1.2[[i]] == "4"] = "t"
}

ape::write.dna(ai_gene_1.2, file ="AI_gene2.fasta", format = "fasta", colsep =" ")



# 2.2
