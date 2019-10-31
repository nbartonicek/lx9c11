#!/bin/bash

#module load gi/star/2.3.0e
module load gi/samtools/1.0
module load gi/novosort/precompiled/1.03.08
module load nenbar/star/2.4.0d
module load borgue/rsem/1.2.26
numcores=15
tag="-P DSGClinicalGenomics"


#directory hierarchy
#raw files directory
homedir="/share/ClusterShare/biodata/contrib/nenbar"
projectDir="$homedir/projects/claudia"
resultsDir="$projectDir/project_results/"

genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/slfn"
#genomeDir="/share/ClusterShare/biodata/contrib/nenbar/genomes/star/mm10_sequin"
#extension of the files to be used
inExt="fq.gz"

#scripts directory
scriptsPath="/share/ClusterShare/biodata/contrib/nenbar/projects/claudia/scripts/QC"
logDir=$scriptsPath"/logs"

#name to append to projectname and create a folder
inType="trimgalore"

projectnames=( "lx9" "NOD" "B6_CVB4" "CVB4_D4" "B6_CVB4" )
#projectnames=( "B6_CVB4" )

for projectname in ${projectnames[@]}; do


        #out directory
        
        inPath="$homedir/projects/claudia/project_results/$projectname.$inType/"

	outPath="/share/ScratchGeneral/nenbar/projects/claudia/project_results/$projectname.star.slfn"
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
	files=`ls $inPath | grep CBDF7ANXX_1_170727_FD02700503_Mouse_GTGAAA_R_170726_CLALOE_RNA_M001`
        for file in ${files[@]};do
                        echo The file used is: $file
                        filesTotal[i]=$file;
                        let i++;
        done 
done;

j=0
echo ${#filesTotal[@]}
while [ $j -lt ${#filesTotal[@]} ]; do

        dir=`echo ${filesTotal[$j]}`
        files=`ls $inPath/$dir/*.$inExt`
	
	inFile1=${files[0]}
	inFile2=${files[1]}
        uniqueID=`basename $dir`
        name=$uniqueID
        outDir=$outPath/$uniqueID/
	mkdir -p $outDir
        rsemDir="$homedir/projects/claudia/project_results/$projectname.rsem/"
        rsem_index="/share/ClusterShare/biodata/contrib/nenbar/genomes/rsem/mm10-gencode.vM9/genome"
        mkdir -p $rsemDir

        echo $name
	#echo $command_line


	starJobName="star."$name
	samSortJobName="samSort"$name
	bamJobName="bam."$name
	sortJobName="sort."$name
	filterJobName="filter."$name
	indexJobName="index."$name
	indexStatsJobName="indexstats."$name
        rsemJobName="rsem."$name
	outSam=$outDir"Aligned.out.sam"
	outSortedSam=$outDir"Aligned.sorted.sam"
	outBam=$outDir"$name.bam"
	outSortedBam=$outDir"$name.sorted.bam"
        outTranscriptomeBam=$outDir"Aligned.toTranscriptome.out.bam"
        outFilteredBam=$outDir"Aligned.toTranscriptome.filtered.out.bam"
	#star_line="/home/nenbar/local/lib/STAR-STAR_2.4.0i/source/STAR --genomeDir $genomeDir --runMode alignReads --readFilesIn $inFile1 $inFile2 --outFileNamePrefix $outDir --runThreadN 4 --outSAMattributes Standard --outSAMstrandField intronMotif --sjdbOverhang 99" 
	
	star_line="STAR \
                --runMode alignReads \
		--genomeDir $genomeDir \
		--readFilesIn $inFile1 $inFile2 \
                --outFileNamePrefix $outDir \
		--runThreadN $numcores \
		--sjdbOverhang 100 \
		--readFilesCommand zcat \
        --outFilterType BySJout \
        --outSAMattributes NH HI AS NM MD\
        --outFilterMultimapNmax 20 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1500000 \
        --alignMatesGapMax 1500000 \
        --alignSJoverhangMin 6 \
        --alignSJDBoverhangMin 1 \
        --quantMode TranscriptomeSAM \
        --outFilterMatchNmin 101"


        #rsem_line="/share/ClusterShare/biodata/contrib/nenbar/projects/MRTA/dnanexus/applets/quantification/rsem-quantify-bam-0.5.v1.2.18/resources/usr/bin/rsem-calculate-expression \
        rsem_line="rsem-calculate-expression \
                --paired-end \
                --bam \
                --no-bam-output \
                --seed 12345 -p 8 \
                --forward-prob 0 \
                $outSortedBam.bam \
                $rsem_index  \
                $rsemDir$name"
        echo $index_line

        filter_line="samtools view -m 16G -@ $numcores $outTranscriptomeBam -f 3 -b > $outFilteredBam"
        sortname_line="novosort -n -m 16G -c $numcores $outFilteredBam >$outSortedBam.bam"
        index_line="samtools index $outSortedBam.bam"
        
        #--outStd BAM_Quant | samtools view -f 3 -u - | novosort -c $numcores -n -m 16G -o $outSortedBam -"
	#bam_line="samtools view -m 16G -h -S $outSam -b -o $outBam"
	#samtools_line="samtools sort -m 16G $outBam $outDir$name.sorted"

        #echo $star_line
        qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores $tag -V"
        #echo "samtools index $outSortedBam"
        #echo "samtools idxstats $outSortedBam | head -n1 | cut -f3 >$outDir$name.tmp;"
	
        $qsubLine -N $starJobName -hold_jid trimgalore  $star_line 
        qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcores -l h_vmem=17G,mem_requested=16G $tag -V"

        $qsubLine -N $filterJobName -hold_jid $starJobName $filter_line
	$qsubLine -N $sortJobName -hold_jid $filterJobName $sortname_line 
        qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 1 $tag -V"
	$qsubLine -N $indexJobName -hold_jid $sortJobName $index_line
	#qsubLine="qsub -q short.q -b y -wd $logDir -j y -R y -pe smp 8 $tag -V"
        #$qsubLine -N $rsemJobName -hold_jid $indexJobName $rsem_line
        
        j=$(($j+1))


done;
#        --outWigType wiggle \
# \
#        --outStd BAM_Quant | samtools view -f 3 -u - | novosort -n -m 16G -o output.bam -
#qsub -q short.q -b y -wd $logDir -j y -R y -pe smp $numcoresSmall $tag -V
#        --outSAMtype BAM Unsorted \
#-outStd BAM_Quant | samtools view -f 3 -u - | novosort -c $numcores -n -m 16G -o $outSortedBam -"

