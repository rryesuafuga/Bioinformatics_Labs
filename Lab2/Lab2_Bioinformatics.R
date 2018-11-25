
library(ape)
library(seqinr)


lizard_seq_seqinr <- read.fasta(file = "lizard_seqs.fasta", seqtype = "DNA",
                                as.string = TRUE, forceDNAtolower = FALSE)

# Assignment 1

## Question 1.1
artificialseq <- list()
set.seed(12345)
for(i in 1:length(lizard_seq_seqinr)){
  
  unsampled_vector <- unlist(strsplit(lizard_seq_seqinr[[i]], " ", fixed = TRUE))
  artificialseq[[i]] <- sample(unsampled_vector, replace = TRUE,
                               prob = NULL)
  artificialseq[[i]] <- paste(artificialseq[[i]], sep=" ", collapse = " ")
  names(artificialseq)[i] <- paste0("lizard",i)
  
}

# need to change writing method
ape::write.dna(artificialseq, file ="artificial_lizard_seqs.fasta", format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)



## Question 1.2
second_artificialseq <- list()
for(i in 1:33){
  
  second_artificialseq[[i]] <- sample(lizards_sequences[[i]], 
                               size=length(lizards_sequences[[i]]))
  
}

tree <- rtree(n = 33)
plot(tree)
