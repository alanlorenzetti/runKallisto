# alorenzetti April 2020

# this script will take the kmer test outputs and plot

# loading libs
library(tidyverse)
theme_set(theme_bw())

# reading results file for riboseq exp
cols=c("kmerSize","libName","isTrimmed","totalReads","psAlnReads")
ribo = read_delim("kmerResultsRiboSeq.txt", delim = " ", col_names = F)
colnames(ribo) = cols

# creating percentage attribute in relation to nonTrimmed total
nonTrimmed = ribo %>% filter(isTrimmed == "nonTrimmed")
trimmed = ribo %>% filter(isTrimmed == "trimmed")

nonTrimmed$pct = nonTrimmed$psAlnReads/nonTrimmed$totalReads
trimmed$pct = trimmed$psAlnReads/nonTrimmed$totalReads

# binding tibbles
ribo = bind_rows(nonTrimmed, trimmed)

# plotting
ggplot(ribo, aes(x=as.character(kmerSize), y=pct)) +
  geom_boxplot() +
  facet_grid(. ~ isTrimmed)

# doing it for total rna libraries

# reading results file for riboseq exp
cols=c("kmerSize","libName","isTrimmed","totalReads","psAlnReads")
ribo = read_delim("kmerResultsTotalRNA.txt", delim = " ", col_names = F)
colnames(ribo) = cols

# creating percentage attribute in relation to nonTrimmed total
nonTrimmed = ribo %>% filter(isTrimmed == "nonTrimmed")
trimmed = ribo %>% filter(isTrimmed == "trimmed")

nonTrimmed$pct = nonTrimmed$psAlnReads/nonTrimmed$totalReads
trimmed$pct = trimmed$psAlnReads/nonTrimmed$totalReads

# binding tibbles
ribo = bind_rows(nonTrimmed, trimmed)

# plotting
ggplot(ribo, aes(x=as.character(kmerSize), y=pct)) +
  geom_boxplot() +
  facet_grid(. ~ isTrimmed)
