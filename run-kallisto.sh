#!/bin/bash

# alorenzetti
# this script is intended to run pseudoalignment
# to quantify transcripts using kallisto

version=0.3.1
lastupdate=20201130

####################################
# PRINT FUNCTIONS
####################################

helpCond () {
  echo -e "Usage:\nbash run-kallisto.sh -t <numberOfThreads> -o <outputDir> -g <GFF3annotationFile> -i <ISannotationFile> -G <GenomeFile> -a <adapterFile> -s <invertStrand> -p <PEorSE> -k <kmersize> -T <trimfiles> -B <bside>\n"

  echo -e "Example:\nbash run-kallisto.sh -t 20 -o output -g ~/dlsm/de_analysis/misc/Hsalinarum-gene-annotation-pfeiffer2019.gff3 -i ~/dlsm/de_analysis/misc/Hsalinarum-IS-annotation-intact-pfeiffer2019.gff3 -G ~/dlsm/misc/Hsalinarum.fa -a ~/dlsm/misc/adap.fa -s y -p n -k 21 -T n -B n"
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
# a = adapter file
# s = invert strand [yn]
# p = paired-end [yn]
# k = kmer size for kallisto index (only odd numbers)
# T = trim files before pseudoalignment
# B = quantify opposite strand too (simplified antisense quantification) [yn]
# h = print help

# argument parser
while getopts 't:o:g:i:G:a:s:p:k:T:B:h' OPTION ; do
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
    a)
      adapter=$OPTARG
      ;;
    s)
      invertstrand=$OPTARG
      ;;
    p)
      pairedend=$OPTARG
      ;;
    k)
      kmer=$OPTARG
      ;;
    T)
      trimming=$OPTARG
      ;;
    B)
      bside=$OPTARG
      ;;
    h)
      helpCond >&2
      exit 0
      ;;
    ?)
      echo "Please, run 'bash run-kallisto.sh -h' for help." >&2
      exit 1
      ;;
  esac
done

# number of threads to run the applications
if [ "$threads" -gt 99 ] ; then echo >&2 "Threads argument value can not be greater than 99" ; exit 1 ; fi

####################################
# HARD CODED VARIABLES
####################################
rawdir=raw
trimmeddir=trimmed
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
for i in kallisto cd-hit perl; do
        command -v $i > /dev/null >&1 || { echo >&2 "$i is not installed. Aborting" ; exit 1; }
done

if [ ! -e /opt/Trimmomatic-0.39/trimmomatic-0.39.jar ] ; then echo >&2 "trimmomatic is not installed. Aborting" ; exit 1; fi

R --slave -e 'if(!require("pacman", quietly=T)){quit(save="no", status=1)}else{quit(save="no", status=0)}'
if [ $? == 1 ] ; then echo >&2 "pacman package is not installed. Aborting" ; exit 1 ; fi

# checking scripts
if [ ! -e get-seqs.R ]  ; then echo >&2 "Missing get-seqs.R script. Aborting" ; exit 1 ; fi

echo "Done!"

####################################
# PROCESSING STARTS HERE
####################################

# creating outputdir
if [[ ! -d $outputdir ]] ; then mkdir $outputdir ; else rm -r $outputdir ; mkdir $outputdir ; fi

# creating fasta file to input in cd-hit
R -q -f get-seqs.R --args $outputdir $geneannot $isannot $genome $bside > $outputdir/get-seqs.log 2>&1

echo "Fasta file generation done!"

# running cd-hit using the fasta file created above
cd-hit -i $outputdir/seqs.fa \
       -o $outputdir/cdhit-output.fa \
       -c 0.95 \
       -T $threads \
       -M 10000 \
       -sc 1 > $outputdir/cdhit-output.log 2>&1

# in case one allowed pseudoAntisense quantification
# we should include them in the ref. transcriptome file
# fasta linearization taken from https://www.biostars.org/p/9262/
if [[ "$bside" == "y" ]] ; then
  grep ">" $outputdir/cdhit-output.fa | perl -pe 's/^>(.*)\|NC.*$/\1/' > $outputdir/nrSeqs.fa.tmp

  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $outputdir/seqs_pseudoAS.fa |\
  perl -pe 's/^>(.*)_pseudoAS.*$/\1/' > $outputdir/seqs_pseudoAS.fa.tmp

  grep -x -A 1 -f $outputdir/nrSeqs.fa.tmp $outputdir/seqs_pseudoAS.fa.tmp |\
  perl -pe 's/^(VNG.*)$/>\1_pseudoAS/' |\
  sed '/^--/d;/^$/d' > $outputdir/seqs_pseudoAS_NR.fa

  awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' $outputdir/cdhit-output.fa |\
  sed '/^$/d' > $outputdir/seqsSense.fa

  cat $outputdir/seqsSense.fa $outputdir/seqs_pseudoAS_NR.fa > $outputdir/cdhit-output_plus_pseudoAS_NR.fa

  rm $outputdir/nrSeqs.fa.tmp $outputdir/seqs_pseudoAS.fa.tmp
fi

echo "CD-HIT done!"

# trimming files
if [[ "$trimming" == "y" ]]; then
  if [ ! -d $trimmeddir ] ; then
    mkdir $trimmeddir

    if [[ "$pairedend" == "y" ]] ; then
      for prefix in $prefixes ; do
        echo "Trimming $prefix"
        R1=$rawdir/$prefix"_R1.fastq.gz"
        R2=$rawdir/$prefix"_R2.fastq.gz"
        outpairedR1=$trimmeddir/$prefix"-paired_R1.fastq.gz"
        outpairedR2=$trimmeddir/$prefix"-paired_R2.fastq.gz"
        outunpairedR1=$trimmeddir/$prefix"-unpaired_R1.fastq.gz"
        outunpairedR2=$trimmeddir/$prefix"-unpaired_R2.fastq.gz"
        logfile=$trimmeddir/$prefix".log"

        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
        -threads $threads \
        $R1 $R2 \
        $outpairedR1 $outunpairedR1 \
        $outpairedR2 $outunpairedR2 \
        ILLUMINACLIP:$miscdir/adap.fa:1:30:10 \
        SLIDINGWINDOW:4:30 \
        MINLEN:16 > $logfile 2>&1
      done
    else
      for prefix in $prefixes ; do
        echo "Trimming $prefix"
        R1=$rawdir/$prefix"_R1.fastq.gz"
        outunpairedR1=$trimmeddir/$prefix"-unpaired_R1.fastq.gz"
        logfile=$trimmeddir/$prefix".log"

        java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar SE \
        -threads $threads \
        $R1 \
        $outunpairedR1 \
        ILLUMINACLIP:${adapter}:1:30:10 \
        SLIDINGWINDOW:4:30 \
        MINLEN:16 > $logfile 2>&1
      done
    fi
  fi
fi

# creating kallisto index
if [[ "$bside" == "y" ]] ; then
  kallisto index -k $kmer -i $outputdir/kallistoidx $outputdir/cdhit-output_plus_pseudoAS_NR.fa > $outputdir/kallisto-index.log 2>&1
else
  kallisto index -k $kmer -i $outputdir/kallistoidx $outputdir/cdhit-output.fa > $outputdir/kallisto-index.log 2>&1
fi

# creating count tables using kallisto
for i in $prefixes ; do
  echo "Processing $i..."

  # creating outputdir
  if [[ ! -d $outputdir/$i ]] ; then mkdir $outputdir/$i ; else rm -r $outputdir/$i ; mkdir $outputdir/$i ; fi

  if [[ "$trimming" == "y" ]]; then
    if [[ "$pairedend" == "y" ]] ; then
        inputfastq=$trimmeddir/${i}"-paired_R1.fastq.gz"
    else
        inputfastq=$trimmeddir/${i}"-unpaired_R1.fastq.gz"
    fi
  else
    inputfastq=$rawdir/${i}"_R1.fastq.gz"
  fi

  # computing readSize and standardDev
  # awk code copied from https://www.biostars.org/p/243552/
  # credits for pierreLindenbaum and guillermo.luque.ds
  meansd=`zcat $inputfastq | awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{m=t/n;printf("%0.f-%0.f\n",m,sqrt(sq/n-m*m));}'`
  mean=${meansd/-*/}
  sd=${meansd/*-/}
  if [[ $sd -eq 0 ]]; then
    sd=1
  fi

  # preparing option according to
  # strandness
  if [[ "$invertstrand" == "y" ]] ; then
    strandflag="--rf-stranded"
  else
    strandflag="--fr-stranded"
  fi

  if [[ "$pairedend" == "y" ]]; then
    R2=${inputfastq/_R1.fastq.gz/_R2.fastq.gz}
    pairedflag=""
  else
    R2=""
    pairedflag="--single -l $mean -s $sd"
  fi

  echo "Read mean size: $mean; SD: $sd"

  # running kallisto quant
  kallisto quant -i $outputdir/kallistoidx \
                 -o $outputdir/$i \
                 --plaintext \
                 $strandflag \
                 -t $threads \
                 $pairedflag \
                 $inputfastq $R2 > $outputdir/$i/kallistoQuant.log 2>&1
  echo "$i done!"
done

# parsing count tables
R -q -f parse-counts.R --args $rawdir $outputdir > $outputdir/parse-counts.log 2>&1
