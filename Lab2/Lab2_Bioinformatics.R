
library(ape)
library(seqinr)
library(phangorn)
library(markovchain)
library(msa)



# Assignment 1

## Question 1.1


clean <- function(template_gene){
  for (i in 1:length(template_gene)) {
    #Remove the " " that created when reading a file
    template_gene[[i]] <- template_gene[[i]][template_gene[[i]]!= " "]
  }
  return(template_gene)
}

### Martin's method
lizards_sequences <- seqinr::read.fasta("lizard_seqs.fasta")
lizards_sequences <- clean(lizards_sequences)
lizard_artificial <- lizards_sequences
nucleotide <- c("a", "c", "g", "t")
for(i in 1:33){
  length_s <- length(lizards_sequences[[i]])
  composition = seqinr::count(lizards_sequences[[i]],1)/length_s
  lizard_artificial[[i]] <- sample(nucleotide, size = length_s, 
                                   prob=composition, replace=TRUE)
  names(lizard_artificial)[i] <- paste0("lizard",i)
  
}
# print(lizard_artificial)

ape::write.dna(lizard_artificial, file ="lizard_artificial.fasta", 
        format = "fasta", append =FALSE, nbcol = 6, colsep = "", colw = 10)

ape::base.freq(as.DNAbin(lizard_artificial), freq = FALSE, all = FALSE)


### My method
# lizard_seq_seqinr <- read.fasta(file = "lizard_seqs.fasta", seqtype = "DNA",
#                                 as.string = TRUE, forceDNAtolower = FALSE)
# 
# artificialseq <- list()
# set.seed(12345)
# for(i in 1:length(lizard_seq_seqinr)){
#   
#   unsampled_vector <- unlist(strsplit(lizard_seq_seqinr[[i]], " ", fixed = TRUE))
#   artificialseq[[i]] <- sample(unsampled_vector, replace = TRUE,
#                                prob = NULL)
#   artificialseq[[i]] <- paste(artificialseq[[i]], sep=" ", collapse = " ")
#   names(artificialseq)[i] <- paste0("lizard",i)
#   
# }
# 
# write.fasta(sequences=artificialseq,
#             names=names(artificialseq),
#             file.out="artificial_lizard_seqs.fasta")




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
# print(lizard_artificial_2)

ape::write.dna(lizard_artificial_2, file ="lizard_artificial_2.fasta", 
               format = "fasta", append =FALSE, nbcol = 6, colsep = "", colw = 10)

ape::base.freq(as.DNAbin(lizard_artificial_2), freq = FALSE, all = FALSE)



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
            file.out="from_tree.fasta", colsep ="")

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

markovchainFit(data=lizards_sequences)
markovchainFit(data=lizard_artificial)
markovchainFit(data=unlist(lizard_artificial_2))



## 2.3 

real_align <- msaClustalW("lizard_seqs.fasta",type="dna") 
real_alignseq<- msaConvert(real_align, type="seqinr::alignment") 
dist_real <- as.matrix(dist.alignment(real_alignseq, "identity"))
heatmap(dist_real)

artificial_align <- msaClustalW("lizard_artificial.fasta",type="dna") 
artificial_alignseq<- msaConvert(real_align, type="seqinr::alignment") 
dist_artificial <- as.matrix(dist.alignment(real_alignseq, "identity"))
heatmap(dist_artificial)



# Assignment 3

## 3.1

# dna <- as.DNAbin(real_alignseq)
# D <- dist.dna(dna, model="TN93")
# F = function(x) njs(x)
# D[is.nan(D)] <- 10
# tree1<-F(D)
# v1 <- boot.phylo(phy =  tree1, x = dna,
#                  FUN = F, quiet = T,
#                  1000)
# plot(tree1)
# v1
# dna <- as.DNAbin(ai1.1_alignseq)
# D <- dist.dna(dna, model="TN93")
# D[is.nan(D)] <- 10
# tree2<-F(D)
# v2 <- boot.phylo(phy =  tree2, x = dna,
#                  FUN = F, quiet = T,
#                  1000)
# plot(tree2)
# v2
# dna <- as.DNAbin(ai1.2_alignseq)
# D <- dist.dna(dna, model="TN93")
# D[is.nan(D)] <- 10
# tree3<-F(D)
# v3 <- boot.phylo(phy =  tree3, x = dna,
#                  FUN = F, quiet = T,
#                  1000)
# plot(tree3)
# v3
# 
