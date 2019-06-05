library(bio3d)
alignment <- read.fasta("alignment_trim.fst")
str(alignment)
align_seq <- seqidentity(alignment, normalize=TRUE, similarity=FALSE, ncore=1, nseg.scale=1)

# Plot a heat map with clustering dendogram
heatmap(align_seq,symm = TRUE,margins = c(12,12))


## BLAST 
con_align <- consensus(alignment)
blast <- blast.pdb(con_align$seq)

##annotation
protein1 <- pdb.annotate("1YVL")
protein2 <- pdb.annotate("1BF5")
protein3 <- pdb.annotate("3CWG")
