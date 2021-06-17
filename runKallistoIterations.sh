#!/bin/bash

for k in 11 13 15 17 19 21 23 25 27 ; do bash run-kallisto.sh -t 12 -o outputTrimmed$k -g ../misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3 -i ../misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 -G ~/dlsm/misc/Hsalinarum_nrtx.fa -x n -a ../misc/adap.fa -s n -p y -k $k -T y -B n > run-kallistoTrimmed$k.log 2> run-kallistoTrimmed$k.err; done
for k in 11 13 15 17 19 21 23 25 27 ; do bash run-kallisto.sh -t 12 -o output$k -g ../misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3 -i ../misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 -G ~/dlsm/misc/Hsalinarum_nrtx.fa -x n -a ../misc/adap.fa -s n -p y -k $k -T n -B n > run-kallisto$k.log 2> run-kallisto$k.err; done
