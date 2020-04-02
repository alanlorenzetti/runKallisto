#!/bin/bash

# alorenzetti 201912
# this script is intended as an extension for frtc pipeline
# after the frtc run, you can use it to quantify transcripts
# using the pseudoalignment approach

version=0.2
lastupdate=20191208

####################################
# PRINT FUNCTIONS
####################################


helpCond () {
  echo -e "Usage:\nbash run-kallisto.sh -t <numberOfThreads> -o <outputDir> -g <GFF3annotationFile> -i <ISannotationFile> -G <GenomeFile> -s <strandedExp> -p <PEorSE> -f <readSize>\n"

  echo -e "Example:\nbash run-kallisto.sh -t 20 -o output -g ~/dlsm/de_analysis/misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3 -i ~/dlsm/de_analysis/misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 -G ~/dlsm/misc/Hsalinarum.fa -s y -p n -f 75"
}

helpFull () {
  echo -e "full help details"
}

####################################
# ARGUMENT VARIABLES
####################################
# t = number of threads
# o = output dir
# g = gene annotation
# i = is annotation
# G = genome file
# s = invert strand [yn]
# p = paired-end [yn]
# f = fragment-size (read size) [int]

# argument parser
while getopts 't:o:g:G:i:s:p:f:h' OPTION ; do
  case $OPTION in
    t)
      threads=$OPTARG
      ;;
    o)
      outputdir=$OPTARG
      ;;
    g)
      geneannot=$OPTARG
      ;;
    i)
      isannot=$OPTARG
      ;;
    G)
      genome=$OPTARG
      ;;
    s)
      invertstrand=$OPTARG
      ;;
    p)
      pairedend=$OPTARG
      ;;
    f)
      readSize=$OPTARG
      ;;
    h)
      helpCond >&2
      exit 0
      ;;
    ?)
      echo "usage is blablalba" >&2
      exit 1
      ;;
  esac
done

if [[ "$inverstrand" == "y" ]] ; then
  strandflag="--rf-stranded" ; else
  strandflag="--fr-stranded"
fi

# getting readSize and computing ~ten percent sd
readSd=`echo "$readSize*0.1" | bc`
readSd=`printf %0.f $readSd`

if [[ "$pairedend" == "y" ]] ; then
  pairedflag="" ; else
  pairedflag="--single -l $readSize -s $readSd"
fi

# number of threads to run the applications
if [ "$threads" -gt 99 ] ; then echo >&2 "Threads argument value can not be greater than 99" ; exit 1 ; fi

####################################
# HARD CODED VARIABLES
####################################
rawdir=raw
prefixes=`ls $rawdir/*.fastq.gz | sed 's/^raw\///;s/_R[12].*$//' | sort | uniq`

####################################
# PROGRAM STAMP
####################################
echo "Kallisto for frtc
Version: $version
Last update: $lastupdate"

####################################
# CALL AND DATE
####################################
dateAndTime=`date`
echo "$dateAndTime"
echo "Call: $0 $@"

####################################
# CHECKING DEPENDENCIES
####################################
echo "Checking dependencies..."

# checking files
if [ -z "$prefixes" ] ; then echo >&2 "FASTQ files not found. Aborting" ; exit 1 ; fi

# checking programs
for i in kallisto cd-hit; do
        command -v $i > /dev/null >&1 || { echo >&2 "$i is not installed. Aborting" ; exit 1; }
done

R --slave -e 'if(!require("pacman", quietly=T)){quit(save="no", status=1)}else{quit(save="no", status=0)}'
if [ $? == 1 ] ; then echo >&2 "pacman package is not installed. Aborting" ; exit 1 ; fi

# checking scripts
if [ ! -e get-seqs.R ]  ; then echo >&2 "Missing get-seqs.R script. Aborting" ; exit 1 ; fi

echo "Done!"

####################################
# PROCESSING STARTS HERE
####################################

# creating outputdir if it doesnt exist
if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; else rm -r $outputdir ; mkdir $outputdir ; fi

# creating fasta file to input in cd-hit
R -q -f get-seqs.R --args $outputdir $geneannot $isannot $genome > $outputdir/get-seqs.log 2>&1

echo "Fasta file generation done!"

# running cd-hit using the fasta file created above
cd-hit -i $outputdir/seqs.fa \
       -o $outputdir/cdhit-output.fa \
       -c 0.95 \
       -T $threads \
       -M 10000 \
       -sc 1 > $outputdir/cdhit-output.log 2>&1

echo "CD-HIT done!"

# creating kallisto index
kallisto index -i $outputdir/kallistoidx $outputdir/cdhit-output.fa > $outputdir/kallisto-index.log 2>&1

# creating count tables using kallisto
for i in $prefixes ; do

  echo "Processing $i..."

  if [[ "$pairedend" == "y" ]]; then
    R2=$rawdir/${i}_R2.fastq.gz; else
    R2=""
  fi

  kallisto quant -i $outputdir/kallistoidx \
                 -o $outputdir/$i \
                 --plaintext \
                 $strandflag \
                 -t $threads \
                 $pairedflag \
                 $rawdir/${i}_R1.fastq.gz $R2 > $outputdir/$i/kallistoQuant.log 2>&1
  echo "$i done!"
done

# parsing count tables
R -q -f parse-counts.R --args $rawdir $outputdir > $outputdir/parse-counts.log 2>&1
