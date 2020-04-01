# alorenzetti 20191206

# loading libraries
library(pacman)
packs = c("tidyverse")

p_load(char=packs)

# parsing args
args = commandArgs(trailingOnly = T)

rawdir = args[1]
#rawdir = "raw"

outputdir = args[2]
#outputdir = "output"

# getting fastq prefixes
prefixes = list.files(path=rawdir) %>% 
  sub("_R[12].fastq.gz","",.) %>%
  unique()

# parsing counts to single table
tableest = NULL
tabletpm = NULL
flag = "init"

for(i in prefixes){
  if(flag == "init"){
    df = read_delim(paste0(outputdir, "/", i, "/", "abundance.tsv"), delim="\t")
    estCounts = df %>% select(-tpm)
    tpm = df %>% select(-est_counts)
    
    names(estCounts)[4] = i
    names(tpm)[4] = i
    
    tableest = bind_cols(tableest, estCounts)
    tabletpm = bind_cols(tabletpm, tpm)
    
    flag = "end"
  } else {
    df = read_delim(paste0(outputdir, "/", i, "/", "abundance.tsv"), delim="\t")
    
    estCounts = df %>% select(est_counts)
    tpm = df %>% select(tpm)
    
    names(estCounts) = i
    names(tpm) = i
    
    tableest = bind_cols(tableest, estCounts)
    tabletpm = bind_cols(tabletpm, tpm)
  }
}

# saving counts
write_delim(tableest, paste0(outputdir, "/", "tableEstCounts.tsv"), delim="\t")
write_delim(tabletpm, paste0(outputdir, "/", "tableTpm.tsv"), delim="\t")