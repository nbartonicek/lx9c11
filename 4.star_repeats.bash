#!/bin/bash

numcores=10
module load nenbar/star/2.4.0d
module load gi/samtools/1.0
#directory hierarchy
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
#temphomedir="/share/ClusterScratch/nenbar"
projectDir="$homedir/projects/claudia"
#tempprojectDir="$temphomedir/capseq"

resultsDir="$projectDir/project_results/"
genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10"
#extension of the files to be used
inExt="gz"

#scripts directory
scriptsPath="$homedir/projects/claudia/scripts/QC"
logDir=$projectDir"/scripts/QC/logs"

#name to append to projectname and create a folder
inType="trimgalore"

projectnames=( "erv" )

for projectname in ${projectnames[@]}; do

    #out directory
    inPath="$projectDir/project_results/$projectname.trimgalore"

	outPath="$resultsDir/$projectname.star_unique/"
    #log and command files for bsub
    logPath="logs"
    commandPath="commands"
    #make the directory structure   
    mkdir -p $outPath
    mkdir -p $logPath
    mkdir -p $commandPath

    rm -f $commandFile

    subs=0

    #get the name of the script for the logs
    scriptName=`basename $0`
    i=0
    echo $inPath
	
	############## fetch file names ##############

	i=0   
	files=( $(ls $inPath/*.fq.gz) )
	for file in ${files[@]};do
	        echo The file used is: $file
	        filesTotal[i]=$file;
	        let i++;
	done;


	j=0
	echo ${#filesTotal[@]}
	while [ $j -lt ${#filesTotal[@]} ]; do

		inFile1=${filesTotal[$j]}
		inFile2=${filesTotal[$(($j+1))]}
        uniqueID=`basename $inFile1 | sed s/_R1.fq.gz//`
        name=$uniqueID
        outDir=$outPath/$uniqueID/
		mkdir -p $outDir
	        
		echo $name
		#echo $command_line


		starJobName="star."$name
		samSortJobName="samSort"$name
		bamJobName="bam."$name
		sortJobName="sort."$name
		
		indexJobName="index."$name
		indexStatsJobName="indexstats."$name
		outSam=$outDir"Aligned.out.sam"
		outSortedSam=$outDir"Aligned.sorted.sam"
		outBam=$outDir"$name.bam"
		outSortedBam=$outDir"$name.sorted.bam"

       star_line="STAR --runMode alignReads \
        --readFilesCommand zcat \
        --genomeDir $genomeDir \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD\
        --outFilterMultimapNmax 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1500000 \
        --alignMatesGapMax 1500000 \
        --alignSJoverhangMin 6 \
        --alignSJDBoverhangMin 1 \
        --readFilesIn $inFile1 $inFile2 \
        --outFileNamePrefix $outDir \
        --runThreadN $numcores \
        --outFilterMatchNmin 101"
		bam_line="samtools view -m 16G -@ $numcores -h -S $outSam -b -o $outBam"
		samtools_line="samtools sort -m 16G -@ $numcores $outBam $outDir$name.sorted"

        echo $star_line
        #echo "samtools index $outSortedBam"
        #echo "samtools idxstats $outSortedBam | head -n1 | cut -f3 >$outDir$name.tmp;"
		qsub -N $starJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $star_line 
		qsub -N $bamJobName -hold_jid $starJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $bam_line
		qsub -N $sortJobName -hold_jid $bamJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $samtools_line 
		qsub -N $indexJobName -hold_jid $sortJobName -b y -wd $logDir -j y -R y -pe smp 1 -V "samtools index $outSortedBam;"
		#qsub -N $indexStatsJobName -hold_jid $indexJobName -b y -wd $logDir -j y -R y -pe smp 1 -V "samtools idxstats $outSortedBam | head -n1 | cut -f3 >$outDir$name.tmp;"

		j=$(($j+2))


	done;
done; 

