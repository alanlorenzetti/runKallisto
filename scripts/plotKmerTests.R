# alorenzetti 202004

# this script will take the kmer test outputs and plot distributions

# Loading libs ################################################ 
# pacman is a nice package manager; make it easier to load and install packages
if(!require("pacman")){install.packages("pacman"); library("pacman")}

# required packs
packs = "tidyverse"

# loading packs
p_load(char = packs)

# setting theme for ggplot2
theme_set(theme_bw())

# Processing starts here ################################################ 
# reading results file
cols=c("kmerSize","libName","isTrimmed","totalReads","psAlnReads")
ribo = read_delim("results/kmerResults.txt", delim = " ", col_names = F)
colnames(ribo) = cols

# creating percentage attribute in relation to nonTrimmed total
nonTrimmed = ribo %>% filter(isTrimmed == "nonTrimmed")
trimmed = ribo %>% filter(isTrimmed == "trimmed")

nonTrimmed$pct = nonTrimmed$psAlnReads/nonTrimmed$totalReads
trimmed$pct = trimmed$psAlnReads/nonTrimmed$totalReads

# binding tibbles
ribo = bind_rows(nonTrimmed, trimmed)

# plotting
p = ggplot(ribo, aes(x=as.character(kmerSize), y=pct)) +
    geom_boxplot() +
    facet_grid(. ~ isTrimmed)

# saving
ggsave(plot = plot,
       filename = "results/kmerResults.png",
       width = 5,
       height = 5,
       dpi = 300,
       units = "in")