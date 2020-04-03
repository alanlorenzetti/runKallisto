#!/bin/bash

kmers="15 17 19 21 23 25 27"
prefixes=`ls raw/*.fastq.gz | sed 's/^raw\///;s/_R[12].*$//' | sort | uniq`

for k in $kmers ; do
  for p in $prefixes ; do
    line=`grep '\[quant\] processed' output${k}/${p}/kallistoQuant.log`
    total=`echo $line | perl -pe 's/^.*processed (.*) reads,.*$/\1/;s/,//g'`
    psaln=`echo $line | perl -pe 's/^.*reads, (.*) reads pseudoaligned.*$/\1/;s/,//g'`
    pct=`echo $psaln/$total | bc`

    lineTrimmed=`grep '\[quant\] processed' outputTrimmed$k/${p}/kallistoQuant.log`
    totalTrimmed=`echo $line | perl -pe 's/^.*processed (.*) reads,.*$/\1/;s/,//g'`
    psalnTrimmed=`echo $line | perl -pe 's/^.*reads, (.*) reads pseudoaligned.*$/\1/;s/,//g'`
    pctTrimmed=`echo $psalnTrimmed/$totalTrimmed | bc`

    echo $total $psaln $pct $totalTrimmed $psalnTrimmed $pctTrimmed
  done
done
