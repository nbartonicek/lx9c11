#!/bin/bash


#module load gi/bwa/0.7.9a
#cd $GENOMES/bwa/mm10_BAC
#bwa index -p mm10_BAC -a bwtsw mm10_BAC.fa 

module load gi/gcc/4.8.2
module load gi/samtools/1.0
module load gi/novosort/precompiled/1.03.08

#directory hierarchy
#raw files directory
numcores=30
tag="-P TumourProgression"
#tag=""

homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/claudia"
resultsDir="$projectDir/project_results/"

genomeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/bwa/mm10/mm10.fa"
#extension of the files to be used

#scripts directory
scriptsPath="/share/ClusterShare/biodata/contrib/nenbar/projects/claudia/scripts/integration"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

i=0
        
inPath="$projectDir/project_results/bwa_mm10"
inPath="/share/ScratchGeneral/nenbar/projects/claudia/project_results/bwa_mm10"
outPath="$projectDir/project_results/svaba_mm10vsWT"
        #log and command files for bsub
logPath="logs"
        #make the directory structure   
mkdir -p $outPath
mkdir -p $logPath

wtFile="/share/ClusterShare/biodata/contrib/nenbar/projects/claudia/project_results/bwa_mm10/H35HGCCX2_8_190819_FD09251145_Mus-musculus__R_190816_CECKIN_DNA_M001.sorted.bam"


files=`ls $inPath/*.bam`
for file in ${files[@]};do
	uniqueID=`basename $file | sed s/.bam//`
	mkdir -p $outPath/$uniqueID
		svaba_line="/share/ScratchGeneral/nenbar/local/lib/svaba/bin/svaba run -t $file -n $wtFile -p $numcores -L 6 -I -a somatic_run -G $genomeFile"
		#svaba_line="/share/ScratchGeneral/nenbar/local/lib/svaba/bin/svaba run -t $file -p $numcores -L 6 -I -a germline_run -G $genomeFile -D /share/ClusterShare/biodata/contrib/nenbar/projects/claudia/annotation/vcf/mgp.v5.indels.pass.chr.sort.vcf"
	qsub -N $uniqueID -b y -wd $outPath/$uniqueID -j y -R y -pe smp $numcores $tag -q short.q -V $svaba_line
done 

