
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(pheatmap)
library(GenomicRanges)
library(DESeq)
library(ChIPpeakAnno)
library(EnsDb.Mmusculus.v79)
library(edgeR)
library(org.Mm.eg.db)

library("BSgenome.Mmusculus.UCSC.mm10")

homedir="/share/ClusterShare/biodata/contrib/nenbar"

#enable for Rstudio
#homedir="../../../../"
timeStamp <- format(Sys.time(), "%Y_%m_%d")


######## directory structure #######
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures/")
annotationDir=paste0(projectDir,"/annotation/API/")

load("../../annotation/API/repeats.Rdata")

#Cell type
repeatCategories<-unique(repeatsGR$type[grepl("LTR",repeatsGR$class)])

#1. Fetch the differentially expressed LTRs



for(repeats in "type"){
	#Load the data
	inDir=paste0(projectDir,"/project_results/erv.repeatOverlap/")
	inFiles=list.files(inDir,full.names=T)
	inFiles=inFiles[grepl(repeats,inFiles)]
	results<-list()
	if(repeats=="class"){
		for (file in inFiles){
			sampleName=basename(file)
			sampleName=gsub(".Rdata","",sampleName)
			load(file)
			results[[sampleName]]<-countsClass
		}
	} else {
		for (file in inFiles){
			sampleName=basename(file)
			sampleName=gsub(".Rdata","",sampleName)
			load(file)
			results[[sampleName]]<-countsType
		}
	}
	#load the library sizes
	sampleSizeDir<-paste0(projectDir,"/project_results/erv.libSize/")
	sizeFiles<-list.files(sampleSizeDir,full.names=T)
	sizeList<-list()
	for(sizeFile in sizeFiles){
		sampleName=basename(sizeFile)
		sampleName=gsub(".Rdata","",sampleName)
		load(sizeFile)
		sizeList[[sampleName]]<-libSize
	}

	data=do.call("cbind",results)
	data=data[,c(2,1,3,4,6,5)]
	temp=data
	#data=data[row.names(data) %in% repeatCategories,]
	colnames(data)=gsub(paste0("_",repeats),"",colnames(data))
	total=unlist(sizeList)
	total=total[colnames(data)]
	data=rbind(data,total=total)

	#Filter out small categories.
	include<-apply(data,1,function(x){sum(x>=50)>=2})
	data<-data[include,]
	data[data==0]=1
	dataNorm<-apply(data,2,function(x){x*1000000/x["total"]})
	dataNorm<-dataNorm[-dim(dataNorm)[1],]
	dataNorm<-log10(dataNorm)

	######### perform differential expression 
	data=temp
	grouping <- colnames(data)[1:6]
	grouping <- factor(grouping)

	countTable=data[,1:6]
	condition=grouping
	cds = newCountDataSet( countTable, condition )
	cds = estimateSizeFactors( cds )
	cds = estimateDispersions( cds, method='blind',sharingMode="fit-only" )

	results<-list()
	for(i in c(2,4,6)){
	    res = nbinomTest( cds, colnames(data)[i-1],colnames(data)[i])
	    sampleName=gsub("_.*","",colnames(data)[i])
	    results[[sampleName]]<-res
	}


	diffExprs<-do.call(cbind,results)
	colnames(diffExprs)<-paste0(rep(names(results),each=8),"_",colnames(res))
	colnames(diffExprs)[1]<-"id"
	diffExprs<-diffExprs[,-grep("_id",colnames(diffExprs))]
	diffExprsLTRs=diffExprs[diffExprs$id %in% repeatCategories,]
	row.names(diffExprsLTRs)=diffExprsLTRs$id

	write.table(diffExprsLTRs,"diffExprsLTRs.txt",row.names=F,quote=F,sep="\t")

	diffExprsLTRsB6=diffExprsLTRs[diffExprsLTRs$B6_padj<0.1&!is.na(diffExprsLTRs$B6_padj),1]
	diffExprsLTRsNod=diffExprsLTRs[diffExprsLTRs$Nod_padj<0.1&!is.na(diffExprsLTRs$Nod_padj),1]
	diffExprsLTRsSJL=diffExprsLTRs[diffExprsLTRs$SJL_padj<0.1&!is.na(diffExprsLTRs$SJL_padj),1]

	ids<-unique(c(diffExprsLTRsB6,diffExprsLTRsNod,diffExprsLTRsSJL))
}
#2. Find what is their location
diffExprsLTRs<-repeatsGR[repeatsGR$type %in% ids]

#3. Find their distance in relationship to the non-diff exprs LTRs to genes/diff.exprs.genes
annoData <- genes(EnsDb.Mmusculus.v79)
seqlevelsStyle(diffExprsLTRs) <- seqlevelsStyle(annoData)
anno <- annotatePeakInBatch(diffExprsLTRs, AnnotationData=annoData)

nonDiffIds<-repeatCategories[!repeatCategories %in% ids]
nodiffExprsLTRs<-repeatsGR[repeatsGR$type %in% nonDiffIds]
seqlevelsStyle(nodiffExprsLTRs) <- seqlevelsStyle(annoData)
annoND <- annotatePeakInBatch(nodiffExprsLTRs, AnnotationData=annoData)

#all genes


pdf("diffExprs_LTRs_distance_types_all_genes.pdf",width=12,height=8)
par(mfrow=c(1,2))
pie1(table(anno$insideFeature),main="Differentially Expressed LTRs")
pie1(table(annoND$insideFeature),main="All LTRs")
dev.off()


#only differentially expressed genes
diffGenes<-read.table("../RNAseq/all_inf.txt",header=T,sep="\t")
annoDiff<-annoData[annoData$gene_name %in% row.names(diffGenes)]
anno <- annotatePeakInBatch(diffExprsLTRs, AnnotationData=annoDiff)
annoND <- annotatePeakInBatch(nodiffExprsLTRs, AnnotationData=annoDiff)



pdf("diffExprs_LTRs_distance_types_diff_genes.pdf",width=12,height=8)
par(mfrow=c(1,2))
pie1(table(anno$insideFeature),main="Differentially Expressed LTRs")
pie1(table(annoND$insideFeature),main="All LTRs")
dev.off()

############### load the protein coding data

homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../../../"
inPath=paste0(homedir,"/projects/claudia/project_results/erv.rsem/")


######## directory structure #######
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
projectname="erv"
robjectsDir = paste(resultsDir,"/erv.Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
sizeDir=paste(resultsDir,"/erv.libSize/",sep="")

outPath=paste0(resultsDir,"/erv.repeatOverlap/")
system(paste("mkdir",outPath))
system(paste("mkdir",sizeDir))

chrs=seqlengths(Mmusculus)[!grepl("_",names(seqlengths(Mmusculus)))]
#names(chrs)=gsub("chr","",names(chrs))
#names(chrs)[25]="MT"
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

#Load both type and class lists
#load(paste0(annotationDir,"class.Rdata"))
#grClass<-class
#load(paste0(annotationDir,"type.Rdata"))
#grType<-type



#########################################
########## 1. Load in the files
#########################################

geneFiles<-list.files(inPath,pattern="genes",full.names=T)

geneResults<-list()

for(file in geneFiles){
  sampleName<-basename(file)
  sampleName<-gsub(".transcriptome.*","",sampleName)
  cat(sampleName)
  cat("\n")
  data<-read.table(file,header=T)  
  geneResults[[sampleName]]<-as.integer(data$expected_count)
}
temp<-geneResults

df<-as.data.frame(geneResults)
row.names(df)<-data$gene_id



#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]


treatment <- rep(c("inflammed", "control"), 3)
type <- c(rep(c("B6", "Nod", "SJL"), each=2))

design <- model.matrix(~type + treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=df[1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt <- glmLRT(fit)




temp=as.data.frame(topTags(lrt,n=40000))


#4. Look at the distribution of distance and expression in vs expression of top over/underexpressed genes
#Take annoData for low and high, genes. And then find those that are in proximity of both categories. 
#Then find their expression and 
diffGenesUp=diffGenes[diffGenes$logFC>0,]
diffGenesUp=diffGenesUp[1:500,]
diffGenesUpIDs<-row.names(diffGenesUp)
diffGenesDown=diffGenes[diffGenes$logFC<0,]
diffGenesDown=diffGenesDown[1:65,]
diffGenesDownIDs<-row.names(diffGenesDown)

diffGenesNonchanging=temp[temp$FDR>0.95,]
diffGenesNonchangingIDs=row.names(diffGenesNonchanging)
#take random 500 genes
diffGenesNonchangingIDs=sample(diffGenesNonchangingIDs,1000)
#Change from ensembl to symbol
egENSEMBL <- toTable(org.Mm.egENSEMBL)
m <- match(diffGenesNonchangingIDs, egENSEMBL$ensembl_id)
entrezGene<-egENSEMBL$gene_id[m]
entrezGene=entrezGene[!is.na(entrezGene)]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(entrezGene, egSYMBOL$gene_id)
symbol<-egSYMBOL$symbol[m]
diffGenesNonchangingIDsSymbol<-symbol[!is.na(symbol)]




annoDiffUp<-annoData[annoData$gene_name %in% row.names(diffGenesUp)]
annoDiffDown<-annoData[annoData$gene_name %in% row.names(diffGenesDown)]
annoDiffNC<-annoData[annoData$gene_name %in% diffGenesNonchangingIDsSymbol]

annoUp <- annotatePeakInBatch(diffExprsLTRs, AnnotationData=annoDiffUp)
annoNDUp <- annotatePeakInBatch(nodiffExprsLTRs, AnnotationData=annoDiffUp)
annoDown <- annotatePeakInBatch(diffExprsLTRs, AnnotationData=annoDiffDown)
annoNDDown <- annotatePeakInBatch(nodiffExprsLTRs, AnnotationData=annoDiffDown)
annoNC <- annotatePeakInBatch(diffExprsLTRs, AnnotationData=annoDiffNC)
annoNDNC <- annotatePeakInBatch(nodiffExprsLTRs, AnnotationData=annoDiffNC)

results<-list()
results[["diffLTR_upGenes"]]<-annoUp
results[["diffLTR_downGenes"]]<-annoDown
results[["nondiffLTR_upGenes"]]<-annoNDUp
results[["nondiffLTR_downGenes"]]<-annoNDDown
results[["diffLTR_control_Genes"]]<-annoNC
results[["nondiffLTRcontrol_Genes"]]<-annoNDNC

distances<-list()
expressions<-list()
for(sampleName in names(results)){
	distances[[sampleName]]<-results[[sampleName]]$shortestDistance
#	expressions[[sampleName]]<-results[[sampleName]]$shortestDistance
}

dataDistances<-unlist(distances)
dataLength<-sapply(distances,length)
dataM<-data.frame(distances=dataDistances,type=rep(names(distances),dataLength))
dataM<-dataM[!is.na(dataM$distances),]
dataM$distances<-log10(dataM$distances)
dataM$type<-factor(dataM$type,levels=c(names(results)))
pdf("distances_of_diffExprsLtrs_to_100diffExprs_genes.pdf",width=12,height=8)
p<-ggplot(dataM,aes(type,distances))
p<-p+geom_boxplot()
p
dev.off()
pdf("distances_of_diffExprsLtrs_to_100diffExprs_genes_violin.pdf",width=12,height=8)
p<-ggplot(dataM,aes(type,distances))
p<-p+geom_violin()
p
dev.off()

dataM<-data.frame(distances=dataDistances,type=rep(names(distances),dataLength))
dataM<-dataM[!is.na(dataM$distances),]
dataM$type<-factor(dataM$type,levels=c(names(results)))
pdf("distances_of_diffExprsLtrs_to_100diffExprs_genes_violin_nonlog.pdf",width=12,height=8)
p<-ggplot(dataM,aes(type,distances))
p<-p+geom_violin()
p
dev.off()


#pvalues
dataM<-data.frame(distances=dataDistances,type=rep(names(distances),dataLength))
dataM<-dataM[!is.na(dataM$distances),]
forKS<-split(dataM$distances,dataM$type)
ks.test(forKS[["diffLTR_upGenes"]],forKS[["diffLTR_control_Genes"]])$p.value
ks.test(forKS[["diffLTR_downGenes"]],forKS[["diffLTR_control_Genes"]])$p.value
ks.test(forKS[["nondiffLTR_upGenes"]],forKS[["nondiffLTRcontrol_Genes"]])$p.value
ks.test(forKS[["nondiffLTR_downGenes"]],forKS[["nondiffLTRcontrol_Genes"]])$p.value
dataM$distances<-log10(dataM$distances)
dataM<-dataM[!is.na(dataM$distances),]

forKS<-split(dataM$distances,dataM$type)
ks.test(forKS[["diffLTR_upGenes"]],forKS[["diffLTR_control_Genes"]])$p.value
ks.test(forKS[["diffLTR_downGenes"]],forKS[["diffLTR_control_Genes"]])$p.value
ks.test(forKS[["nondiffLTR_upGenes"]],forKS[["nondiffLTRcontrol_Genes"]])$p.value
ks.test(forKS[["nondiffLTR_downGenes"]],forKS[["nondiffLTRcontrol_Genes"]])$p.value
ks.test(forKS[["diffLTR_control_Genes"]],forKS[["nondiffLTRcontrol_Genes"]])$p.value

forKS<-split(dataM$distances,dataM$type)
wilcox.test(forKS[["diffLTR_upGenes"]],forKS[["diffLTR_control_Genes"]],alternative = "greater")$p.value
wilcox.test(forKS[["diffLTR_downGenes"]],forKS[["diffLTR_control_Genes"]],alternative = "greater")$p.value
wilcox.test(forKS[["nondiffLTR_upGenes"]],forKS[["nondiffLTRcontrol_Genes"]],alternative = "greater")$p.value
wilcox.test(forKS[["nondiffLTR_downGenes"]],forKS[["nondiffLTRcontrol_Genes"]],alternative = "greater")$p.value
wilcox.test(forKS[["diffLTR_control_Genes"]],forKS[["nondiffLTRcontrol_Genes"]],alternative = "greater")$p.value
wilcox.test(forKS[["diffLTR_upGenes"]],forKS[["diffLTR_downGenes"]],alternative = "less")$p.value
wilcox.test(forKS[["nondiffLTR_upGenes"]],forKS[["nondiffLTR_downGenes"]],alternative = "less")$p.value



genePromoters<-promoters(annoDiffUp)
diffOverlap<-sum(countOverlaps(diffExprsLTRs,genePromoters)>0)/length(diffExprsLTRs)
nodiffOverlap<-sum(countOverlaps(nodiffExprsLTRs,genePromoters)>0)/length(nodiffExprsLTRs)


genePromoters<-promoters(annoDiffDown)
diffOverlap<-sum(countOverlaps(diffExprsLTRs,genePromoters)>0)/length(diffExprsLTRs)
nodiffOverlap<-sum(countOverlaps(nodiffExprsLTRs,genePromoters)>0)/length(nodiffExprsLTRs)


repeatNames<-c("RLTR45-int","MER21C","ERVB2_1-I_MM-int")
diffGenesUp=diffGenes[diffGenes$logFC>0,]
diffGenesUp=diffGenesUp[1:500,]
diffGenesUpIDs<-row.names(diffGenesUp)
diffGenesDown=diffGenes[diffGenes$logFC<0,]
diffGenesDown=diffGenesDown[1:65,]
diffGenesDownIDs<-row.names(diffGenesDown)
annoNodif<-annoData[annoData$gene_id %in% diffGenesNonchangingIDs]

for(repeatName in repeatNames){

	#check what percent overlaps up promoters, down promoters, randoms
	#get locations
	cat("Repeat type: ")
	cat(repeatName)
	cat("\n")

	genePromotersUp<-promoters(annoDiffUp)
	genePromotersDown<-promoters(annoDiffDown)
	genePromotersDown<-promoters(annoNodif)
	
	diffOverlapUp<-sum(countOverlaps(genePromotersUp,diffExprsLTRs[diffExprsLTRs$type %in% repeatName])>0)/length(genePromotersUp)
	diffOverlapDown<-sum(countOverlaps(annoDiffDown,diffExprsLTRs[diffExprsLTRs$type %in% repeatName])>0)/length(annoDiffDown)
	diffOverlapNoDif<-sum(countOverlaps(genePromotersUp,diffExprsLTRs[diffExprsLTRs$type %in% repeatName])>0)/length(genePromotersUp)
	cat("Percent of UP genes with the repeat in promoter\n")
	cat(diffOverlapUp)
	cat("\n")
	cat("Percent of DOWN genes with the repeat in promoter\n")
	cat(diffOverlapDown)	
	cat("\n")
	cat("Control (random 500 nodiff genes)\n")	
	cat(diffOverlapNoDif)
	cat("\n")
}









