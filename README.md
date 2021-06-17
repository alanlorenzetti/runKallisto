# Transcriptome quantification using kallisto

This set of scripts was conceived to pseudoalign Illumina RNA-Seq reads to a prokaryotic non-redundant reference transcriptome. Briefly, `run-kallisto.sh` will perform:

* **Generation of a non-redundant reference transcriptome**: It will take a reference genome (fasta) and its correspondent annotation (gff3), extract the sequences, and them perform redundancy removal using cd-hit; or it could accept a transcriptome fasta file as well.
* **Read trimming and mapping**: Reads will be trimmed using Trimmomatic and then mapped to the previously generated non-redundant reference transcriptome using Kallisto.
* **Output parsing**: Individual Kallisto count tables (one per library) will be parsed and made available in two bulky files: i. containing TPMs count for each transcript for each library. ii. containing raw estimated counts for each transcript for each sample.

A few other scripts are provided to perform additional tests, but up to the moment they have not been documented and should still be considered for development purposes only.
 
