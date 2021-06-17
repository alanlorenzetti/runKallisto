#!/bin/bash

kmers="11 13 15 17 19 21 23 25 27"
prefixes=`ls raw/*.fastq.gz | sed 's/^raw\///;s/_R[12].*$//' | sort | uniq`

# creating directory to store aligned reads
# per kmer size tested
rm -r results

if [[ ! -d results ]] ; then mkdir results ; fi

touch results/kmerResults.txt

# generating aligned reads file
for k in $kmers ; do
  for p in $prefixes ; do
    line=`grep '\[quant\] processed' output${k}/${p}/kallistoQuant.log`
    total=`echo $line | perl -pe 's/^.*processed (.*) reads,.*$/\1/;s/,//g'`
    psaln=`echo $line | perl -pe 's/^.*reads, (.*) reads pseudoaligned.*$/\1/;s/,//g'`
    pct=`echo $psaln/$total | bc`

    lineTrimmed=`grep '\[quant\] processed' outputTrimmed$k/${p}/kallistoQuant.log`
    totalTrimmed=`echo $lineTrimmed | perl -pe 's/^.*processed (.*) reads,.*$/\1/;s/,//g'`
    psalnTrimmed=`echo $lineTrimmed | perl -pe 's/^.*reads, (.*) reads pseudoaligned.*$/\1/;s/,//g'`
    pctTrimmed=`echo $psalnTrimmed/$totalTrimmed | bc`

    echo "$k $p nonTrimmed $total $psaln" >> results/kmerResults.txt
    echo "$k $p trimmed $totalTrimmed $psalnTrimmed" >> results/kmerResults.txt

  done
done
