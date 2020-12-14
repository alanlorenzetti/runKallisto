# alorenzetti 20200519
# this script is going to take supplied annotation files
# and extract sequences

# args
args = commandArgs(trailingOnly = T)

# temp paths
# outputdir = "~/gdrive/runKallisto"
# genefile = "~/gdrive/dlsm/de_analysis/misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3"
# isfile = "~/gdrive/dlsm/de_analysis/misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3"
# genomefile = "~/gdrive/dlsm/misc/Hsalinarum.fa"
# bside = "y"

# def paths
outputdir = args[1]
genefile = args[2]
isfile = args[3]
genomefile = args[4]
bside = args[5]

# removing / from the end of outputdir
outputdir = sub("/$","",outputdir)

# fasta output
fastaoutput = paste0(outputdir, "/seqs.fa")

# there will be a pseudoantisense file if required
# this set is the reverse complement of each gene
if(bside == "y"){bsideoutput = paste0(outputdir, "/seqs_pseudoAS.fa")}

# packages
library(pacman)

packs = c("tidyverse", "BSgenome",
          "rtracklayer", "Biostrings")

p_load(char=packs)

# loading annotations and parsing
is = rtracklayer::import(isfile, format = "gff")
is$locus_tag = is$Name

gene = rtracklayer::import(genefile, format = "gff")
gene = subset(gene, type == "gene")

# replacing transposase gene loci by IS annotation
#gene = gene[!gene %over% is]

# merging gene and is annotations
#annot = c(gene, is)

# not using is features this time
annot = gene

# loading genome
geno = readDNAStringSet(genomefile)
names(geno) = names(geno) %>% sub(" .*$","",.)

# getting seqs and parsing appropriate names
annotSeqs = getSeq(geno, annot)
names(annotSeqs) = paste0(annot$locus_tag, "|",
                          seqnames(annot),":",
                          start(annot),"-",
                          end(annot),":",
                          strand(annot))

# in case pseudoantisenses are required
if(bside == "y"){
  annotPsiAS = invertStrand(annot)
  annotPsiASseqs = getSeq(geno, annotPsiAS)
  names(annotPsiASseqs) = paste0(annotPsiAS$locus_tag, "_pseudoAS", "|",
                                 seqnames(annotPsiAS),":",
                                 start(annotPsiAS),"-",
                                 end(annotPsiAS),":",
                                 strand(annotPsiAS))
  
  # saving pseudoantisenses to file
  writeXStringSet(annotPsiASseqs, bsideoutput, format = "fasta")
}

# saving seqs to file
writeXStringSet(annotSeqs, fastaoutput, format = "fasta")
