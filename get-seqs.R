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

# def paths
outputdir = args[1]
genefile = args[2]
isfile = args[3]
genomefile = args[4]

# removing / from the end of outputdir
outputdir = sub("/$","",outputdir)

# fasta output
fastaoutput = paste0(outputdir, "/seqs.fa")

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

# saving seqs to file
writeXStringSet(annotSeqs, fastaoutput, format = "fasta")