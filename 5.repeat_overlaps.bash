#!/bin/bash

module load briglo/R/3.4.2 
#module load gi/bcftools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/claudia"
resultsDir="$projectDir/project_results/"

#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/claudia/scripts/repeats"
logDir=$projectDir"/scripts/repeats/logs"
mkdir -p $logDir

projectnames=( "lx9" "NOD" "B6_CVB4" "CVB4_D4" )
projectnames=( "B6_CVB4" )
for projectname in ${projectnames[@]};do
	inPath="/share/ScratchGeneral/nenbar/projects/claudia/project_results/$projectname.starUnique"
	files=`ls $inPath/**/*.sorted.bam.bam | grep FD02704009`
	for file in ${files[@]}; do
	        sampleName=`basename $file | sed s/.genomeAligned.sorted.bam.bam//`
	        echo $sampleName
	        qsub -q long.q -b y -j y -N repeats.$sampleName -wd $logDir -pe smp 12 -V "R --vanilla <$scriptsPath/1a.repeat_overlap_parallel_lx9.R --projectname=$projectname --UniqueID=$sampleName"
	       
	done;
done;
