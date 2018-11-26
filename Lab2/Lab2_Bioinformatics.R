
library(ape)
library(seqinr)
library(phangorn)
library(markovchain)


lizard_seq_seqinr <- read.fasta(file = "lizard_seqs.fasta", seqtype = "DNA",
                                as.string = TRUE, forceDNAtolower = FALSE)

# Assignment 1

## Question 1.1

### Martin's method
lizards_sequences <- read.fasta("lizard_seqs.fasta")
lizard_artificial <- lizards_sequences
nucleotide <- c("a", "c", "g", "t")
for(i in 1:33){
  length_s <- length(lizards_sequences[[i]])
  composition = seqinr::count(lizards_sequences[[i]],1)/length_s
  lizard_artificial[[i]] <- sample(nucleotide, size = length_s, 
                                   prob=composition, replace=TRUE)
  names(lizard_artificial)[i] <- paste0("lizard",i)
  
}
print(lizard_artificial)

ape::write.dna(lizard_artificial, file ="lizard_artificial.fasta", 
               format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)



### My method
artificialseq <- list()
set.seed(12345)
for(i in 1:length(lizard_seq_seqinr)){
  
  unsampled_vector <- unlist(strsplit(lizard_seq_seqinr[[i]], " ", fixed = TRUE))
  artificialseq[[i]] <- sample(unsampled_vector, replace = TRUE,
                               prob = NULL)
  artificialseq[[i]] <- paste(artificialseq[[i]], sep=" ", collapse = " ")
  names(artificialseq)[i] <- paste0("lizard",i)
  
}


write.fasta(sequences=artificialseq,
            names=names(artificialseq),
            file.out="artificial_lizard_seqs.fasta")


### Comment on base composition



## Question 1.2
lizard_artificial_2 <- lizards_sequences
nucleotide <- c("a", "c", "g", "t")
for(i in 1:33){
  length_s <- length(lizards_sequences[[i]])
  composition = seqinr::count(lizards_sequences[[i]],1)/length_s
  lizard_artificial_2[[i]] <- sample(nucleotide, size = length_s, 
                                   prob=composition, replace=TRUE)
  names(lizard_artificial_2)[i] <- paste0("lizard",i)

}
print(lizard_artificial_2)

ape::write.dna(lizard_artificial_2, file ="lizard_artificial_2.fasta", 
               format = "fasta", append =FALSE, nbcol = 6, colsep = " ", colw = 10)



tree <- rtree(n = 33)
plot(tree)
l2 <- length(unlist(lizards_sequences))
bc2 <- seqinr::count(unlist(lizards_sequences),1)/length(unlist(lizards_sequences))
Q2 <- matrix(c(0.9,0.05,0.03,0.02,0.05,0.87,0.07,0.01,
               0.03,0.07,0.89,0.01,0.02,0.01,0.01,0.96), 
             ncol = 4, nrow = 4)
from_tree <- simSeq(tree, l = l2, type="DNA", bf=as.vector(bc2), Q=Q2)

write.fasta(sequences=from_tree,
            names=names(lizard_artificial_2),
            file.out="from_tree.fasta")

### will come back to this later




# Assignment 2

## 2.1

ape::base.freq(as.DNAbin(lizard_artificial), freq = FALSE, all = FALSE)
GC.content(as.DNAbin(lizard_artificial))

ape::base.freq(as.DNAbin(lizard_artificial_2), freq = FALSE, all = FALSE)
GC.content(as.DNAbin(lizard_artificial_2))

### How do we know whether to choose Homosapien or anything else?
### The results are not similar
### There are many stop codons marked by stars
### Yes. There are many stop codons in the true sequence


## 2.2

lizards_sequences_markov <- 
  markovchainListFit(data=lizards_sequences, name = "lizards_sequences_markov")

artificial_markov <- 
  markovchainListFit(data=lizard_artificial, name = "lizard_artificial")

artificial_markov_2 <- 
  markovchainListFit(data=lizard_artificial_2, name = "lizard_artificial_2")
### I get a very large list



