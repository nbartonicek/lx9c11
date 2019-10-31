library(GenomicRanges)
library(annotatr)
library(rtracklayer)
library(ShortRead)
library("BSgenome.Mmusculus.UCSC.mm10")

projectDir="/share/ClusterShare/biodata/contrib/nenbar/projects/claudia/"
inPath=paste0(projectDir,"project_results/svaba_mm10vsWT")

samples<-list.files(inPath)

liftoverMm9<-function(gr1){
  gr.mm9=gr1
  ch = import.chain(paste0(projectDir,"/annotation/UCSC/mm9ToMm10.over.chain"))
  seqlevelsStyle(gr.mm9) = "UCSC"  # necessary
  gr.mm10 = liftOver(gr.mm9, ch)
  gr.mm10  = unlist(gr.mm10 )
  genome(gr.mm10) = "mm10"
  return(gr.mm10)
}

results<-GRangesList()
resultsVCF<-list()
for(sampleName in samples){
	inDir<-paste0(inPath,"/",sampleName)
	inFile<-paste0(inDir,"/somatic_run.svaba.somatic.indel.vcf")
	data=read.table(inFile)
	gr<-GRanges(seqnames=data$V1,IRanges(start=data$V2,width=nchar(as.character(data$V4))))
	results[[sampleName]]=gr
	resultsVCF[[sampleName]]=data
}

#do the filtering
#how many total: 1297+1372=2669, but only 1739 unique in total
grReduce<-reduce(unlist(results))

#how many overlapping between two samples: 904
gr1<-results[[1]]
gr2<-results[[2]]

mat<-findOverlaps(gr1,gr2)
gr<-gr1[unique(queryHits(mat))]

#how many in each length category
dataShort<-resultsVCF[[1]][unique(queryHits(mat)),]

#insertions vs deletions
dataShort$l1<-nchar(as.character(dataShort$V4))
dataShort$l2<-nchar(as.character(dataShort$V5))

#length of 
dataShort$indel<-apply(dataShort[,c("l1","l2")],1,function(x){ifelse(x[1]>x[2],"deletion","insertion")})
dataShort$width<-apply(dataShort[,c("l1","l2")],1,function(x){abs(x[1]-x[2])})

#how many of those overlap each of the functional categories
#for inserts vs deletions
#for length: 1-5, 5-10, 10-50,50-100,>100
dataShort$widthCategory<-cut(dataShort$width,c(0,5,10,50,100,200))

#for functional category
if(!file.exists("annotations.Rdata")){
	annots=builtin_annotations()
	annots = annots[grepl("mm10",annots)]
	annots = annots[!grepl("cpg",annots)]

# Build the annotations (a single GRanges object)
	annotations = build_annotations(genome = 'mm10', annotations = annots)
	save(annotations,file="annotations.Rdata")
} else {load("annotations.Rdata")}

# Intersect the regions we read in with the annotations
gr_annotated = annotate_regions(
    regions = gr,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
# A GRanges object is returned
df_gr_annotated = data.frame(gr_annotated)
gr_annotated$annotation=df_gr_annotated$annot.type

mat<-findOverlaps(gr,gr_annotated)
grAnnotatedL<-split(gr_annotated[subjectHits(mat)],queryHits(mat))
annotation<-sapply(grAnnotatedL,function(x){paste0(unique(x$annotation),collapse=",")})
dataShort$annotation<-annotation
temp<-dataShort

dataShort$annotation[grepl("mm10_genes_intergenic",dataShort$annotation)]<-"intergenic"
dataShort$annotation[grepl("mm10_lncrna_gencode",dataShort$annotation)&grepl("mm10_genes_exons",dataShort$annotation)]<-"lncRNA_exons"
dataShort$annotation[grepl("mm10_lncrna_gencode",dataShort$annotation)&grepl("mm10_genes_introns",dataShort$annotation)]<-"lncRNA_introns"
dataShort$annotation[grepl("mm10_genes_exons",dataShort$annotation)&grepl("mm10_genes_3UTRs",dataShort$annotation)]<-"UTR_exons"
dataShort$annotation[grepl("mm10_lncrna_gencode",dataShort$annotation)]<-"lncRNA_introns"
dataShort$annotation[grepl("mm10_genes_exons",dataShort$annotation)]<-"protein_coding_exons"
dataShort$annotation[grepl("mm10_enhancers_fantom",dataShort$annotation)]<-"enhancers"
dataShort$annotation[grepl("mm10_genes_intronexonboundaries",dataShort$annotation)]<-"intronexonboundaries"
dataShort$annotation[grepl("mm10_genes_exonintronboundaries",dataShort$annotation)]<-"intronexonboundaries"
dataShort$annotation[grepl("mm10_genes_promoters",dataShort$annotation)]<-"promoters"
dataShort$annotation[grepl("mm10_genes_introns",dataShort$annotation)]<-"protein_coding_introns"


############## for overlap with Lx9
lx9df=read.table("../../annotation/API/lx9.bed")
lx9<-GRanges(seqnames=lx9df$V1,IRanges(lx9df$V2,lx9df$V3))

mat<-findOverlaps(gr,lx9)
dataShort$lx9Overlap<-0
dataShort$lx9Overlap[unique(queryHits(mat))]<-1


############## for overlap with offtarget site
data=read.table("lx9_offtargets.txt",header=F)
grOfftarget<-GRanges(seqnames=data$V2,IRanges(start=data$V4,width=1))
grOfftarget<-resize(grOfftarget,50,fix="center")
grMm10<-liftoverMm9(grOfftarget)

mat<-findOverlaps(gr,grMm10)
dataShort$offtarget=0
dataShort$offtarget[unique(queryHits(mat))]<-1

############## for overlap with low complexity 
a=DNAStringSet(dataShort$V4)
dataShort$lDust<-dustyScore(DNAStringSet(dataShort$V4))
dataShort$rDust<-dustyScore(DNAStringSet(dataShort$V5))

dataShort$lowComplexity<-dataShort$lDust|dataShort$rDust>20


table(dataShort[,c("annotation","lowComplexity")])
table(dataShort[,c("widthCategory","lowComplexity")])

dataShortGR<-GRanges(seqnames=dataShort$V1,IRanges(start=dataShort$V2,width=1))
dataShortGR<-resize(dataShortGR,20,fix="center")
seqs<-getSeq(Mmusculus,dataShortGR)
dataShort$lowComplexityGenome<-dustyScore(seqs)>20
dataShort$LC<-dataShort$lowComplexityGenome|dataShort$lowComplexity

dataShortNonLC<-dataShort[dataShort$LC==FALSE,]

table(dataShortNonLC[,c("annotation","lowComplexity")])
table(dataShortNonLC[,c("widthCategory","lowComplexity")])











