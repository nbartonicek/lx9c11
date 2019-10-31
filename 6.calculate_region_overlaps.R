
library(GenomicRanges)
library(ShortRead)
library(R.utils)
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)

timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84

UniqueID="CBDF7ANXX_170221_FD02704009_Mouse_ACAGTG_R_170206_CLALOE_RNA_M001"
projectname="B6_CVB4"
args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["UniqueID"]])){UniqueID = args$UniqueID} 
if (!is.null(args[["projectname"]])){projectname = args$projectname} 



homedir="/share/ScratchGeneral/nenbar"
inPath=paste0(homedir,"/projects/claudia/project_results/",projectname,".starUnique/",UniqueID,"/")


######## directory structure #######
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0("/share/ClusterShare/biodata/contrib/nenbar/projects/claudia/annotation/API/")
robjectsDir = paste(resultsDir,"/all.Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
sizeDir=paste(resultsDir,"/libSize/",sep="")

outPath=paste0(resultsDir,"/repeatOverlap/")
system(paste("mkdir",outPath))
system(paste("mkdir",sizeDir))

chrs=seqlengths(Mmusculus)[!grepl("_",names(seqlengths(Mmusculus)))]
#chrs=seqlengths(Mmusculus)

#names(chrs)=gsub("chr","",names(chrs))
#names(chrs)[25]="MT"
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

#Load both type and class lists
load(paste0(annotationDir,"class.Rdata"))
grClass<-class
load(paste0(annotationDir,"type.Rdata"))
grType<-type

results<-list()

inFile=list.files(inPath,pattern=paste0(UniqueID,".genomeAligned.sorted.bam.bam$"), full.names=T)
if(!length(inFile)==0){
    #Load the bam file and find overlaps with repeats

    sampleName=basename(inFile)
    sampleName=gsub(".sorted.bam.bam","",sampleName)
    outFileClass=paste0(outPath,sampleName,".Rdata")
    cat(".")
    ##############   Read Bam files     #####################
    #chromosome lengths

    #Instead of loading the whole file, will load just the specific bit
    #chrRanges<-IRanges(start=start(gr),width=as.integer(width(gr)),names=as.character(seqnames(gr)))
    #chrRangesL=split(chrRanges,1:length(chrRanges))
    #chrRangeList <- do.call(RangesList, chrRangesL)

    what <- c("qname","rname","strand","pos","qwidth")
    flag <- scanBamFlag(isUnmappedQuery=FALSE)
    #repeatGR<-unlist(grClass)
    #param <- ScanBamParam(what=what, which=repeatGR,flag=flag)

    #load in the BAM file, the count for each of the reads is in the name of the read, after "c"
    
    param <- ScanBamParam(what=what, which=gr,flag=flag)
    system.time(reads<-readGAlignmentPairs(inFile,param=param))
    cat("loaded the reads\n")
    libSize=length(reads)
    save(libSize,file=paste0(sizeDir,sampleName,".Rdata"))
    
     
    countsType=countOverlaps(grType,reads)
    countsClass=countOverlaps(grClass,reads)
    cat("counted the overlaps\n")

    outFileClass=paste0(outPath,UniqueID,"_class.Rdata")
    outFileType=paste0(outPath,UniqueID,"_type.Rdata")
    cat(outFileClass)
    save(countsClass,file=outFileClass)
    save(countsType,file=outFileType)

}    








