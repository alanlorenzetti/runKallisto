#!/bin/bash

for k in 11 13 15 17 19 21 23 25 27 ; do bash run-kallisto.sh -t 12 -o outputTrimmed$k -g ../misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3 -i ../misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 -G ../misc/Hsalinarum.fa -a ../misc/adap.fa -s y -p n -k $k -T y -B y > run-kallistoTrimmed$k.log 2> run-kallistoTrimmed$k.err; done
for k in 11 13 15 17 19 21 23 25 27 ; do bash run-kallisto.sh -t 12 -o output$k -g ../misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3 -i ../misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 -G ../misc/Hsalinarum.fa -a ../misc/adap.fa -s y -p n -k $k -T n -B y > run-kallisto$k.log 2> run-kallisto$k.err; done
