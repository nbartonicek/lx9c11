#!/bin/bash


#module load gi/bwa/0.7.9a
#cd $GENOMES/bwa/mm10_BAC
#bwa index -p mm10_BAC -a bwtsw mm10_BAC.fa 

module load gi/bwa/0.7.8
module load gi/gcc/4.8.2
module load gi/samtools/1.0
module load gi/novosort/precompiled/1.03.08

#directory hierarchy
#raw files directory
numcores=36
tag="-P TumourProgression"
#tag=""

homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/claudia"
resultsDir="$projectDir/project_results/"

genomeFile="/share/ClusterShare/biodata/contrib/nenbar/genomes/bwa/mm10/mm10.fa"
#extension of the files to be used
inExt="fastq.gz"


#scripts directory
scriptsPath="/share/ClusterShare/biodata/contrib/nenbar/projects/claudia/scripts/integration"
logDir=$scriptsPath"/logs"
mkdir -p $logDir

i=0


        
        inPath=$projectDir"/raw_files/DNA/"
	inPath="/share/ScratchGeneral/nenbar/projects/claudia/raw_files/DNA/WT"
	outPath="$projectDir/project_results/bwa_mm10"
        #log and command files for bsub
        logPath="logs"
        commandPath="commands"
        #make the directory structure   
        mkdir -p $outPath
        mkdir -p $logPath

        subs=0

        #get the name of the script for the logs
        scriptName=`basename $0`
        
        echo $inPath
	files=`ls $inPath/*.fastq.gz`
	for file in ${files[@]};do
                        echo The file used is: $file
                        filesTotal[i]=$file;
                        let i++;
        done 

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

        #dir=`echo ${filesTotal[$j]}`
        #files=`ls $inPath/$dir/*.$inExt`
	
        inFile1=${filesTotal[$j]}
        inFile2=${filesTotal[$(($j+1))]}
        uniqueID=`basename $inFile1 | sed s/_R1.fastq.gz//`
        name=$uniqueID
        outDir=$outPath/$name/
	mkdir -p $outDir
        echo $name
	#echo $command_line


	bwaJobName="bwa."$name
	samSortJobName="samSort"$name
	bamJobName="bam."$name
	sortJobName="sort."$name
	
	indexJobName="index."$name
	indexStatsJobName="indexstats."$name
	outSam=$outDir"Aligned.out.sam"
	outSortedSam=$outDir"Aligned.sorted.sam"
	outBam=$outDir"$name.bam"
	outSortedBam=$outDir"Aligned.sortedByCoord.out.bam"

	#bwa_line="bwa mem $genomeFile -t $numcores $inFile1 $inFile2 >$outDir$name.sam"
        bwa_line="bwa mem $genomeFile -t $numcores $inFile1 $inFile2 | samtools view -b -u -S - |novosort -m 7G -c $numcores -i -o $outDir$name.sorted.bam -"
        #| novosort -c $numcores -i -o $outDir$name.sorted.bam"
        #samtools_line="~/local/lib/samtools/bin/samtools view -m 16G -@ $numcores -b -u -S - >$outDir$name.bam"
        #$outDir$name.sam >$outDir$name.bam"
        #"| novosort -c $numcores -i -o "$outDir$name.sam" - "
        

#samtools view -m 8G -@ 1 -u -S - |


        
        #echo $bwa_line
        qsub -N $bwaJobName -hold_jid bwa_mm10 -b y -wd $logDir -j y -R y -pe smp $numcores $tag -q short.q -V $bwa_line 
	#qsub -N $bamJobName -hold_jid $bwaJobName -b y -wd $logDir -j y -R y -pe smp $numcores  -P DSGClinicalGenomics -q short.q -V $bwa_line
	#qsub -N $sortJobName -hold_jid $bamJobName -b y -wd $logDir -j y -R y -pe smp $numcores -V $samtools_line 

        j=$(($j+2))


done;
#        bwa_line="/home/nenbar/local/lib/bwa/bin/Linux_x86_64/bwa --genomeDir $genomeDir \
#                --readFilesCommand zcat \
#                --runMode alignReads \
#                --readFilesIn $inFile1 $inFile2 \
#                --outFilterMultimapNmax 1 \
#                --outFileNamePrefix $outDir \
#                --runThreadN $numcores \
#                --sjdbOverhang 149 \
#                --outFilterType BySJout \
#                --outSAMattributes NH HI AS NM MD\
#                --outFilterMismatchNmax 999 \
#                --outFilterMismatchNoverReadLmax 0.04 \
#                --alignIntronMin 20 \
#                --alignIntronMax 1 \
#                --alignMatesGapMax 1500000 \
#                --alignSJoverhangMin -1 \
#                --alignSJDBoverhangMin 5 \
#                --outSAMtype BAM Unsorted \
#                --chimSegmentMin 20 \
#                --chimJunctionOverhangMin 20 \
#                --chimScoreMin 1 \
#                --chimScoreDropMax 30 \
#
#                --chimScoreSeparation 10" 
#
