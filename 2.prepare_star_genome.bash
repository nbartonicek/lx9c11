#!/bin/bash

#load all modules

module load gi/samtools/1.0
module load nenbar/star/2.4.0d
#number of cores
ncore=8

#project directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"

#genome directory
genome="slfn"
genomesDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/$genome/"
#annotationDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/sequin"
starLine="STAR --runMode genomeGenerate --genomeSAindexNbases 7 --sjdbOverhang 49 --genomeDir $genomesDir  --genomeFastaFiles $genomesDir/$genome.fa --runThreadN $ncore"
#qsub -N star_hg19 -b y -cwd -j y -R y -pe smp $ncore -V $starLine 

qsub -N mm10 -b y -cwd -j y -R y -pe smp $ncore -V $starLine
