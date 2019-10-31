
library(GenomicRanges)
#library(ShortRead)
#library(R.utils)
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
#library(RUVSeq)
library(org.Mm.eg.db)
library(DESeq)
library(ChIPseeker)
library(ChIPpeakAnno)
library(ggrepel)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)
library(plyr)
library(ggrastr)
library(pheatmap)
library(JASPAR2018)
library(TFBSTools)
library(fgsea)
library(ggrepel)
library(clusterProfiler)

timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84


#setwd("~/contrib/nenbar/projects/claudia/scripts/RNAseq")
#homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../../../"



######## directory structure #######
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
robjectsDir = paste(resultsDir,"/Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
sizeDir=paste(resultsDir,"/libSize/",sep="")

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
sampleNames=c( "lx9_wt", "lx9_CVB4", "B6_wt", "B6_CVB4" )
geneResults<-list()

inPath=paste0(homedir,"projects/claudia/project_results/rsem/")
countFile<-paste0(homedir,"projects/claudia/project_results/Robjects/geneResults.Rdata")
#if(!file.exists(countFile)){
  for(file in list.files(inPath,full.names=T)){
  
    sampleName<-basename(file)
    sampleName<-gsub("_Mouse.*","",sampleName)
    sampleName<-gsub(".*_","",sampleName)
  
    fileSequin<-gsub("rsem","sequin",file)
    cat(sampleName)
    cat("\n")
    data<-read.table(file,header=T)  
    dataS<-read.table(fileSequin,header=T)  
  
    geneResults[[sampleName]]<-c(as.integer(data$expected_count),as.integer(dataS$expected_count))
    temp<-geneResults
  }
  df<-as.data.frame(geneResults)
  row.names(df)<-c(as.character(data$gene_id),as.character(dataS$gene_id))
  save(df,file=countFile)
#} else {load(countFile)}

#df<-as.data.frame(geneResults)
#row.names(df)<-c(as.character(data$gene_id),as.character(dataS$gene_id))

#load the repetitive elements
repeatResults<-list()

inPathRepeats=paste0(homedir,"/projects/claudia/project_results/repeatOverlap/")
repeatFiles<-list.files(inPathRepeats,pattern="type",full.names=T)
for(file in repeatFiles){
  sampleName<-basename(file)
  #sampleName<-gsub("_type.Rdata","",sampleName)
  sampleName<-gsub("_Mouse.*","",sampleName)
  sampleName<-gsub(".*_","",sampleName)
  cat(sampleName)
  cat("\n")
  load(file)
  repeatResults[[sampleName]]<-countsType
}
dfRepeat<-do.call("cbind",repeatResults)
dfRepeat<-dfRepeat[,colnames(dfRepeat) %in% colnames(df)]
df<-rbind(df,dfRepeat)

temp=df

ids<-read.table("ids.txt")
newNames=ids[match(names(geneResults),ids$V1),"V2"]
#newNames=newNames[!is.na(newNames)]

dfS=df[,7:21]
newNames=newNames[7:21]
colnames(dfS)=newNames
dfS=dfS[,c(1:6,10:15)]
colnames(dfS)[10:12]=c("LX9_CVB4_1","LX9_CVB4_2","LX9_CVB4_3")
df=dfS

#-dim(dfRepeat)[1]-1):dim(df)[1]

#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]
df$EntrezGene[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]<-row.names(df)[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]
df$symbol[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]<-row.names(df)[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]

#eliminate duplicated symbols
o <- order(rowSums(df[,1:12]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$symbol)
df<-df[!d,]

#eliminate lowly expressed
#include<-apply(df[,1:12],1,function(x){sum(x>=10)>=2})
#df<-df[include,]
#include<-apply(df[,1:12],1,function(x){sum(x==0)>=6})
#df<-df[!include,]
#
include<-apply(df[,1:12],1,function(x){sum(x>=5)>=2})
df<-df[include,]
#include<-apply(df[,1:12],1,function(x){sum(x==0)>=3})
#df<-df[!include,]



df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

temp=df

#effect of infection on B6
dfLINE<-df[,c(1:12)]
type <- c(rep(c("B6","LX9"), each=6))
treatment <- rep(rep(c("CVB4","WT"), each=3),2)
design <- model.matrix(~type+type:treatment)
expr <- DGEList(counts=dfLINE[1:12])
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_B6_CVB4 <- glmLRT(fit)
logFCs_B6_CVB4<-as.data.frame(topTags(lrt_B6_CVB4,n=20000))
logFCs_B6_CVB4=logFCs_B6_CVB4[!grepl("^R1_",row.names(logFCs_B6_CVB4)),]
logFCs_B6_CVB4=logFCs_B6_CVB4[!grepl("^R2_",row.names(logFCs_B6_CVB4)),]

write.csv(logFCs_B6_CVB4,"B6_CVB4_pvalue_repeats.csv")

plus_B6_CVB4<-logFCs_B6_CVB4[logFCs_B6_CVB4$logFC<=-1,]
minus_B6_CVB4<-logFCs_B6_CVB4[logFCs_B6_CVB4$logFC>=1,]
significantIDs_B6_CVB4<-c(row.names(plus_B6_CVB4[1:25,]),row.names(minus_B6_CVB4[1:25,]))
fittedValues_B6_CVB4<-fitted.values(lrt_B6_CVB4)
fittedValues100_B6_CVB4<-fittedValues_B6_CVB4[row.names(fittedValues_B6_CVB4) %in% significantIDs_B6_CVB4,]

fittedValues50_B6_CVB4=fittedValues100_B6_CVB4
fittedValues50_B6_CVB4<-fittedValues50_B6_CVB4[,c(6,5,4,3,2,1)]
colsums<-apply(fittedValues50_B6_CVB4[,4:6],1,median)
fittedValues50_B6_CVB4<-fittedValues50_B6_CVB4[rev(order(colsums)),]
pheatmap(log1p(fittedValues50_B6_CVB4),cex=0.7,cluster_rows = F,cluster_cols=F)

write.csv(fittedValues_B6_CVB4,"B6_CVB4_all_repeats.csv")
write.csv(fittedValues100_B6_CVB4,"B6_CVB4_topUpDown200_repeats.csv")

repeats_B6_CVB4<-logFCs_B6_CVB4[row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat),]
write.csv(repeats_B6_CVB4,"logFC_B6_CVB4_repeats_only.csv")
break()

repeats_B6_CVB4[row.names(repeats_B6_CVB4) %in% c("Lx9","L1_Mus3"),]
fittedValues_B6_CVB4[row.names(fittedValues_B6_CVB4) %in% c("Lx9","L1_Mus3"),]
df[row.names(df) %in% c("Lx9","L1_Mus3"),]

repeats_Lx9CVB4_vs_CVB4[row.names(repeats_Lx9CVB4_vs_CVB4) %in% c("Lx9","L1_Mus3"),]
fittedValues_LINE[row.names(fittedValues_LINE) %in% c("Lx9","L1_Mus3"),]

#kegg analysis

#take a table without repeats
df<-dfS

#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]
df$EntrezGene[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]<-row.names(df)[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]
df$symbol[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]<-row.names(df)[c(dim(df)[1]-76-dim(dfRepeat)[1]):dim(df)[1]]

#eliminate duplicated symbols
o <- order(rowSums(df[,1:12]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$symbol)
df<-df[!d,]

#change the id to entrez
o <- order(rowSums(df[,1:12]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$EntrezGene)
df<-df[!d,]
df<-df[!is.na(df$EntrezGene),]
row.names(df)<-df$EntrezGene
df<-df[!grepl("^R",row.names(df)),]


dfLINE<-df[,c(1:12)]
type <- c(rep(c("B6","LX9"), each=6))
treatment <- rep(rep(c("CVB4","WT"), each=3),2)
design <- model.matrix(~type+type:treatment)
expr <- DGEList(counts=dfLINE[1:12])
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_KEGG <- glmLRT(fit)




kegLx9wt <- kegga(lrt_KEGG, species="Mm")
topKEGG(kegLx9wt)
write.table(topKEGG(kegLx9wt,number=Inf),"KEGG_categories_lx9wt.txt",row.names=T,quote=F,sep="\t")




#volcano plot
genes_B6_CVB4<-logFCs_B6_CVB4[!(row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat)),]
thresholdDE<-genes_B6_CVB4$FDR<=0.001
genes_B6_CVB4$threshold <- thresholdDE 
genes_B6_CVB4_ordered <- genes_B6_CVB4[order(genes_B6_CVB4$FDR), ] 
genes_B6_CVB4_ordered$genelabels <- ""
genes_B6_CVB4_ordered$genelabels[1:50] <- rownames(genes_B6_CVB4_ordered)[1:50]
#repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[repeats_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf(paste0(imageDir,"/volcano_logFCs_Lx9CVB4_vs_CVB4_rasterized.pdf"),width=3,height=3, useDingbats=FALSE)
ggplot(genes_B6_CVB4_ordered) +
  geom_point_rast(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("Lx9+CVB4 vs Lx9 Genes") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels),size=1.5,segment.size = 0.05) +
  xlim(-15,15) +
  scale_color_manual(values=c("#4363AE","#A51F22")) +
  theme_bw(base_size=3)  

dev.off()



#make the lines under the volcano

logFCs_B6_CVB4<-as.data.frame(topTags(lrt_B6_CVB4,n=40000))
genes_B6_CVB4<-logFCs_B6_CVB4[!(row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat)),]

rank<-1:dim(genes_B6_CVB4)[1]
names(rank)<-row.names(genes_B6_CVB4[order(genes_B6_CVB4$logFC),])

thresholdDE<-genes_B6_CVB4$FDR<=0.001
genes_B6_CVB4$threshold <- thresholdDE 

deGenes<-row.names(genes_B6_CVB4)[genes_B6_CVB4$threshold]
pathways<-list()
pathways[["slfnGenes"]]<-deGenes[grepl("Slfn",deGenes)]
pathways[["ifnGenes"]]<-deGenes[grepl("Irf",deGenes)]

#cxclGenes"
cxclGenes<-bitr_kegg("mmu04060", "Path", "ncbi-geneid", "mmu")
m <- match(cxclGenes[,2], egSYMBOL$gene_id)
cxclGenes$symbol<-egSYMBOL$symbol[m]
pathways[["cxclGenes"]]<-unique(cxclGenes$symbol)

#RIG-I-LIKE receptor signaling pathway"
rigGenes<-bitr_kegg("mmu04622", "Path", "ncbi-geneid", "mmu")
m <- match(rigGenes[,2], egSYMBOL$gene_id)
rigGenes$symbol<-egSYMBOL$symbol[m]
pathways[["RIGGenes"]]<-unique(rigGenes$symbol)


fgseaRes <- fgsea(pathways = pathways, 
                  stats = rank,
                  minSize=5, maxSize=500,
                  nperm=100000)

plotGseaTable(pathways[1:4], rank, fgseaRes, gseaParam = 0.5)


sampleNames<-c("slfnGenes","ifnGenes","cxclGenes","RIGGenes")
for(sampleName in sampleNames){
  pathway=pathways[[sampleName]]
  cat(sampleName)


  pathwayGenes=pathway
  stats=rank
  gseaParam=1
  rnk <- rank(stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  
  leRanks<-data.frame(rank=rnk[names(rnk) %in% pathwayGenes],genes=names(rnk)[names(rnk) %in% pathwayGenes])
  
  temp=toPlot
  merged=merge(toPlot,leRanks,by.x="x",by.y="rank",all.x=T)
  #merged$genes[is.na(merged$genes)]=""
  toPlot=merged
  dup=duplicated(toPlot$genes)
  toPlot$genes[dup]=NA
  pdf(paste0("~/Desktop/",sampleName,"_GSEA2.pdf"),width=18,height=4)
  p<-plotEnrichment(pathways[[sampleName]],stats=statsAdj) + labs(title=sampleName)
  print(p)
  dev.off()
  pdf(paste0("~/Desktop/",sampleName,"_GSEA.pdf"),width=9,height=2)
  g <- ggplot(toPlot, aes(x = x, y = y)) + 
    geom_point(color = "green", size = 0.1) + 
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
    geom_hline(yintercept = 0, colour = "black") + 
    geom_line(color = "green") + 
    theme_bw() + 
    ylim(c(-1.5,1.5)) +
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x, 
                                                               y = -diff*4, xend = x, yend = diff*4), size = 0.4) + 
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score") + 
    geom_text_repel(data=subset(toPlot,!is.na(genes)),
                    angle        = 0,
                    segment.size  = 0.8,
                    segment.colour = "gray80",
                    segment.alpha = 0.5,
                    direction     = "y",
                    ylim=c(0.5,1.5),size=4,
                    #      xlim=(c(0,n+1)),
                    aes(x = x, y = diff/2, 
                        label = genes))  
  print(g)
  dev.off()
  
}











data=read.table("KEGG_categories_lx9wt.txt",header=T,sep="\t",stringsAsFactors=F)
data$Pathway<-gsub(" pathway","",data$Pathway)
dfU<-data[order(data$P.Up),][1:30,]
dfD<-data[order(data$P.Down),][1:12,]
dfU<-dfU[!(dfU$Pathway %in% c("Osteoclast differentiation","Pathways in cancer")),][1:12,]

df<-dfU
df$ratio<-df$Up/(df$Up+df$Down)
df<-df[,c(1,2,5,7)]
df$P.Up<--log10(df$P.Up)
df$Pathway<-factor(df$Pathway,levels=rev(df$Pathway))
df1<-df
colnames(df1)[3]<-"P"
#pdf(paste0(imageDir,"/lx9_pathway_final_up.pdf"),width=4.8,height=2.8, useDingbats=FALSE)
pdf(paste0(imageDir,"/lx9_pathway_final_up.pdf"),width=6,height=6, useDingbats=FALSE)
p<-ggplot(df1,aes(P,Pathway))
p<-p+geom_point(aes(size=N,colour=ratio))
p<-p+xlim(c(0,40))
p<-p+scale_size(limits=c(0,1500),range=c(1,10),breaks=c(-100,0,25,50,100,200,500,1000,1500))
p<-p+scale_colour_gradient2(limits=c(0.5,1),midpoint = 0.5,low="#4363AE",high="#A51F22")
p<-p+theme(panel.background = element_rect(fill = "white"),panel.grid.major = element_line(colour = "gray90",size=0.2))
p
dev.off()

df<-dfD
df$ratio<-df$Up/(df$Up+df$Down)
df<-df[,c(1,2,6,7)]
df$P.Down<--log10(df$P.Down)
df$Pathway<-factor(df$Pathway,levels=rev(df$Pathway))
df2=df
colnames(df2)[3]<-"P"
pdf(paste0(imageDir,"/lx9_pathway_final_down.pdf"),width=5.5,height=6, useDingbats=FALSE)
p<-ggplot(df2,aes(P,Pathway))
p<-p+geom_point(aes(size=N,colour=ratio))
p<-p+xlim(c(0,75))
p<-p+scale_size(limits=c(-100,1500),range=c(1,10),breaks=c(-100,0,25,50,100,500,1000,1500))
p<-p+scale_colour_gradient2(limits=c(0,1),midpoint = 0.5,low="#4363AE",high="#A51F22")
p<-p+theme(panel.background = element_rect(fill = "white"),panel.grid.major = element_line(colour = "gray90",size=0.2))
p
dev.off()
#same for GO categories
data=read.table("GO_categories_B6CVB4_vs_CVB4_clean.txt",header=T,sep="\t",stringsAsFactors=F)
data=data[data$Ont=="BP",]
dfU<-data[order(data$P.Up),][1:30,]
dfD<-data[order(data$P.Down),][1:10,]
dfU<-dfU[1:10,]

df<-dfU
df$ratio<-df$Up/df$N
df<-df[,c(1,3,6,8)]
df$P.Up<--log10(df$P.Up)
df$Term<-factor(df$Term,levels=rev(df$Term))
df1<-df
colnames(df1)[3]<-"P"
p<-ggplot(df1,aes(P,Term))
p<-p+geom_point(aes(size=N,colour=ratio))
p<-p+xlim(c(0,100))
p<-p+scale_size(limits=c(-100,1000),range=c(1,10),breaks=c(-100,0,25,50,100,200,500,1000))
p<-p+scale_colour_gradient2(limits=c(0,1),midpoint = 0.5,low="darkblue",high="darkred")
p<-p+theme(panel.background = element_rect(fill = "white"),panel.grid.major = element_line(colour = "gray90",size=0.2))
p


df<-dfD
df$ratio<-df$Up/df$N
df<-df[,c(1,3,7,8)]
df$P.Down<--log10(df$P.Down)
df$Term<-factor(df$Term,levels=rev(df$Term))
df2=df
colnames(df2)[3]<-"P"
p<-ggplot(df2,aes(P,Term))
p<-p+geom_point(aes(size=N,colour=ratio))
p<-p+xlim(c(0,25))
p<-p+scale_size(limits=c(-100,1000),range=c(1,10),breaks=c(-100,0,25,50,100,200,500,1000))
p<-p+scale_colour_gradient2(limits=c(0,0.75),midpoint = 0.5,low="darkblue",high="darkred")
p<-p+theme(panel.background = element_rect(fill = "white"),panel.grid.major = element_line(colour = "gray90",size=0.2))
p




#volcano plot
repeatAnnotation<-read.table("../../annotation/API/unique_categories.txt",header=F,sep="\t")
colnames(repeatAnnotation)<-c("id","class")
repeats_B6_CVB4<-logFCs_B6_CVB4[row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat),]

thresholdDE<-repeats_B6_CVB4$FDR<=0.20
repeats_B6_CVB4$threshold <- thresholdDE 
repeats_B6_CVB4_ordered <- repeats_B6_CVB4[order(repeats_B6_CVB4$FDR), ] 
#repeats_B6_CVB4_ordered$genelabels <- ""
#repeats_B6_CVB4_ordered$genelabels[1:25] <- rownames(repeats_B6_CVB4_ordered)[1:25]
repeats_B6_CVB4_ordered$genelabels <- rownames(repeats_B6_CVB4_ordered)

merged<-merge(repeats_B6_CVB4_ordered,repeatAnnotation,by.x="genelabels",by.y="id",all.x=T)
repeats_B6_CVB4_ordered<-merged
repeats_B6_CVB4_ordered<-repeats_B6_CVB4_ordered[order(repeats_B6_CVB4_ordered$FDR),]
repeats_B6_CVB4_orderedShort<-repeats_B6_CVB4_ordered[repeats_B6_CVB4_ordered$class %in% c("LTRs","Type I Transposons/LINE"),]
repeats_B6_CVB4_orderedShort$repeatlabels <- ""

repeats_B6_CVB4_orderedShort$repeatlabels[repeats_B6_CVB4_orderedShort$class=="LTRs"][1:10]<-repeats_B6_CVB4_orderedShort$genelabels[repeats_B6_CVB4_orderedShort$class=="LTRs"][1:10]
repeats_B6_CVB4_orderedShort$repeatlabels[repeats_B6_CVB4_orderedShort$class=="Type I Transposons/LINE"][1:10]<-repeats_B6_CVB4_orderedShort$genelabels[repeats_B6_CVB4_orderedShort$class=="Type I Transposons/LINE"][1:10]
repeats_B6_CVB4_orderedShort$repeatlabels[repeats_B6_CVB4_orderedShort$genelabels=="Lx9"]<-"Lx9"
#repeats_B6_CVB4_orderedShort$FC<-ifelse(repeats_B6_CVB4_orderedShort$logFC>0,exp(repeats_B6_CVB4_orderedShort$logFC,-exp(repeats_B6_CVB4_orderedShort)))
#repeats_B6_CVB4_orderedShort$repeatlabels[1:25] <- repeats_B6_CVB4_orderedShort$genelabels[1:25]

#repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[repeats_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf(paste0(imageDir,"/volcano_repeats_logFCs_B6_CVB4_vs_CVB4_rasterized_byType.pdf"),width=12,height=8)
ggplot(repeats_B6_CVB4_orderedShort) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("B6+CVB4 vs B6 Repeats") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  facet_wrap(~class,ncol=2)+
  scale_x_continuous(breaks = round(seq(-7, 7, by = 1),1)) +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = repeatlabels),min.segment.length=0.1) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()


repeats_B6_CVB4<-logFCs_B6_CVB4[row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat),]

#plus_B6_CVB4<-repeats_B6_CVB4[repeats_B6_CVB4$logFC<=-1,]
#minus_B6_CVB4<-repeats_B6_CVB4[repeats_B6_CVB4$logFC>=1,]
#significantIDs_B6_CVB4<-c(row.names(plus_B6_CVB4[1:10,]),row.names(minus_B6_CVB4[1:10,]))
#fittedValues_B6_CVB4<-fitted.values(lrt_B6_CVB4)
#fittedValues100_B6_CVB4<-fittedValues_B6_CVB4[row.names(fittedValues_B6_CVB4) %in% significantIDs_B6_CVB4,]


fittedValues50_B6_CVB4<-fittedValues100_B6_CVB4[,c(6,5,4,3,2,1)]
#colsums<-apply(fittedValues50_B6_CVB4[,4:6],1,median)
#fittedValues50_B6_CVB4<-fittedValues50_B6_CVB4[rev(order(colsums)),]
pheatmap(log1p(fittedValues50_B6_CVB4),cex=0.7,cluster_rows = T,cluster_cols=F)






# Figure 1b. distance

# 1. get promoters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=10000, downstream=0)
#promoterTest <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
geneid <- mapIds(org.Mm.eg.db, names(promoter), "SYMBOL","ENTREZID")
promoter=promoter[!duplicated(geneid)]
names(promoter)<-geneid[!duplicated(geneid)]
promoterDE<-promoter[names(promoter) %in% row.names(minus_B6_CVB4[1:51,])]
strand(promoterDE)="*"
promoterDE_plus<-promoterDE[strand(promoterDE)=="+"]
promoterDE_minus<-promoterDE[strand(promoterDE)=="-"]

# 2. Get repetitive elements distances
#load all LTRs, split into DE vs non DE
#load(paste0(annotationDir,"/class.Rdata"))
LTRs<-class[grep("LTR",names(class))]
LINEs<-class[grep("LINE",names(class))]
LTRs<-unlist(LTRs)
LINEs<-unlist(LINEs)

mat<-findOverlaps(promoterDE,LTRs)
LTRs<-LTRs[unique(subjectHits(mat))]
names(LTRs)<-1:length(LTRs)
mat<-findOverlaps(promoterDE,LINEs)
LINEs<-LINEs[unique(subjectHits(mat))]
names(LINEs)<-1:length(LINEs)

#options:
#1. find all that oeverlap within 10kb
#2. for each find distance (Chippeakanno)
#3. plot the distance while retaining info

LTRs$ids<-paste0("LTR_",1:length(LTRs))
LINEs$ids<-paste0("LINE_",1:length(LINEs))

seqsLTRs<-as.character(getSeq(Mmusculus, LTRs))
seqsLINEs<-as.character(getSeq(Mmusculus, LINEs))

mappLTRS<-list()
for(i in 1:length(seqsLTRs)){
  cat(".")
  load(paste0(robjectsDir,"LTR_",i,".Rdata"))
  mappLTRS[[i]]<-mappSum
}
LTRs$mapp<-unlist(mappLTRS)

mappLINEs<-list()
for(i in 1:length(seqsLINEs)){
  cat(".")
  load(paste0(robjectsDir,"LINE_",i,".Rdata"))
  mappLINEs[[i]]<-mappSum
}
LINEs$mapp<-unlist(mappLINEs)

peakAnnoLTRs <- annotatePeak(LTRs, tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnnoLINEs <- annotatePeak(LINEs, tssRegion=c(-10000, 10000), TxDb=txdb, annoDb="org.Mm.eg.db")

#make data frame with element, distance, gene
dfLTR<-as.data.frame(values(peakAnnoLTRs@anno)[c("type","distanceToTSS","SYMBOL","ids","mapp")])
dfLINE<-as.data.frame(values(peakAnnoLINEs@anno)[c("type","distanceToTSS","SYMBOL","ids","mapp")])

df<-rbind(dfLTR,dfLINE)
df<-df[df$distanceToTSS<0,]
df<-df[df$SYMBOL %in% row.names(minus_B6_CVB4[1:51,]),]
df<-df[order(df$distanceToTSS,decreasing=T),]

dup<-duplicated(df$SYMBOL)
symbols<-df$SYMBOL[!dup]
symbols<-symbols[!is.na(symbols)]

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "IRF1"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
irf1=PFMatrixList[["MA0050.1"]]

opts <- list()
opts[["species"]] <- 9606
opts[["name"]] <- "STAT3"
#opts[["type"]] <- "SELEX"
opts[["all_versions"]] <- TRUE
PFMatrixList <- getMatrixSet(JASPAR2018, opts)
stat1=PFMatrixList[["	MA0144.2"]]

#now for each element find if there is a match in it
LINEsShort<-LINEs[LINEs$ids %in% df$ids]
LTRsShort<-LTRs[LTRs$ids %in% df$ids]
gr<-c(LINEsShort,LTRsShort)

seqs<-getSeq(Mmusculus, gr)
gr$irf<-sapply(seqs,function(x){length(matchPWM(as.matrix(irf1),x))>0})
gr$stat<-sapply(seqs,function(x){length(matchPWM(as.matrix(stat1),x))>0})
gr$score<-gr$irf+gr$stat

df$TFscore<-gr$score[match(df$ids,gr$ids)]

df$SYMBOL<-factor(df$SYMBOL, levels=rev(symbols))
rm(mapp)
df$mapp<-log10(df$mapp)
df$TFscore <- as.factor(df$TFscore)
pdf(paste0(imageDir,"/LTR_LINE_distance2genes.pdf"),width=8,height=8)
p<-ggplot(df,aes(SYMBOL,distanceToTSS))
p<-p+geom_point(aes(colour=TFscore))
p<-p+coord_flip()
p<-p+scale_color_manual(values=c("gray30", "orange", "red"))
p
dev.off()



LINEsLx9=LINEs[LINEs$type=="Lx9"]
mat<-findOverlaps(promoterDE,LINEsLx9)
LINEsLx9<-LINEsLx9[unique(subjectHits(mat))]
LINEsLx9$type
promoterDE[unique(queryHits(mat))]


# 4. get their locations
mat<-findOverlaps(LTRs,promoters_DE)
LTRsOverlap<-LTRs[unique(queryHits(mat))]







deLTRs<-row.names(repeats_Lx9CVB4_vs_CVB4)[repeats_Lx9CVB4_vs_CVB4$FDR<0.01]

LTRsDE=LTRsOverlap[LTRsOverlap$type %in% deLTRs]
LTRsNonDE=LTRsOverlap[!(LTRsOverlap$type %in% deLTRs)]

allLTRs<-c(LTRsDE,LTRsNonDE)
values(allLTRs)$type=rep(c("DE","nonDE"),times=c(length(LTRsDE),length(LTRsNonDE)))

tssDE=allLTRs

sampleName="UP"
gr=resize(promoters_UP,1,fix="center")
gr$score=1
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score",
                         extend = 2500, mean_mode = "w0", w = 10)

sampleName="DOWN"
gr=resize(promoters_DOWN,1,fix="center")
gr$score=1
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score",
                         extend = 2500, mean_mode = "w0", w = 10)


pdf(paste0(imageDir,"/DE_genes_DE_LTRs_rasterized.pdf"),width=8,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="UP",column_title = "UP")+
  EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="DOWN",column_title = "DOWN")
dev.off()











# 1. get DE
n=500
significantIDs_LINE<-c(row.names(plus_Lx9CVB4_vs_CVB4[1:n,]),row.names(minus_Lx9CVB4_vs_CVB4[1:n,]))


# 2. get promoters
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")), upstream=2500, downstream=2500)

ids_DE<-df$EntrezGene[df$symbol %in% significantIDs_LINE]
ids_Pos<-df$EntrezGene[df$symbol %in% row.names(plus_Lx9CVB4_vs_CVB4[1:n,])]
ids_Neg<-df$EntrezGene[df$symbol %in% row.names(minus_Lx9CVB4_vs_CVB4[1:n,])]


promoters_DE<-promoter[promoter$gene_id %in% ids_DE]
promoters_UP<-promoter[promoter$gene_id %in% ids_Pos]
promoters_DOWN<-promoter[promoter$gene_id %in% ids_Neg]

# 3. get LTRs that are DE / nonDE
repeats_Lx9CVB4_vs_CVB4<-logFCs_Lx9CVB4_vs_CVB4[row.names(logFCs_Lx9CVB4_vs_CVB4) %in% row.names(dfRepeat),]

#load all LTRs, split into DE vs non DE
load(paste0(annotationDir,"/class.Rdata"))
LTRs<-class[grep("LTR",names(class))]
LTRs<-unlist(LTRs)

# 4. get their locations
mat<-findOverlaps(LTRs,promoters_DE)
LTRsOverlap<-LTRs[unique(queryHits(mat))]

# 5. plot
deLTRs<-row.names(repeats_Lx9CVB4_vs_CVB4)[repeats_Lx9CVB4_vs_CVB4$FDR<0.01]

LTRsDE=LTRsOverlap[LTRsOverlap$type %in% deLTRs]
LTRsNonDE=LTRsOverlap[!(LTRsOverlap$type %in% deLTRs)]

allLTRs<-c(LTRsDE,LTRsNonDE)
values(allLTRs)$type=rep(c("DE","nonDE"),times=c(length(LTRsDE),length(LTRsNonDE)))

tssDE=allLTRs

sampleName="UP"
gr=resize(promoters_UP,1,fix="center")
gr$score=1
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score",
    extend = 2500, mean_mode = "w0", w = 10)

sampleName="DOWN"
gr=resize(promoters_DOWN,1,fix="center")
gr$score=1
mat1 = normalizeToMatrix(gr, tssDE, value_column = "score",
    extend = 2500, mean_mode = "w0", w = 10)


pdf(paste0(imageDir,"/DE_genes_DE_LTRs_rasterized.pdf"),width=8,height=8)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="UP",column_title = "UP")+
EnrichedHeatmap(mat1,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, col = c("white", "darkblue"),name ="DOWN",column_title = "DOWN")
dev.off()



sampleName="All DE promoters"
gr=resize(promoters_DE,1,fix="center")
gr$score=1
mat0 = normalizeToMatrix(gr, tssDE, value_column = "score",
    extend = 5000, mean_mode = "w0", w = 30)

pdf(paste0(imageDir,"/DE_all_genes_DE_LTRs_rasterized.pdf"),width=4,height=6)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="DE Gene Promoters",column_title = "DE Gene Promoters")

dev.off()



########### all genes

sampleName="All  promoters"
gr=resize(promoter,1,fix="center")
gr$score=1

mat<-findOverlaps(LTRs,promoter)
LTRsOverlap<-LTRs[unique(queryHits(mat))]

# 5. plot
deLTRs<-row.names(repeats_Lx9CVB4_vs_CVB4)[repeats_Lx9CVB4_vs_CVB4$FDR<0.01]

LTRsDE=LTRsOverlap[LTRsOverlap$type %in% deLTRs]
LTRsNonDE=LTRsOverlap[!(LTRsOverlap$type %in% deLTRs)]

allLTRs<-c(LTRsDE,LTRsNonDE)
values(allLTRs)$type=rep(c("DE","nonDE"),times=c(length(LTRsDE),length(LTRsNonDE)))

tssDE=allLTRs




mat0 = normalizeToMatrix(gr, tssDE, value_column = "score",
    extend = 2500, mean_mode = "w0", w = 30)

pdf(paste0(imageDir,"/all_genes_DE_LTRs_rasterized.pdf"),width=4,height=6)
#lgd = Legend(at = c(unique(tssDE$type)), title = "Tamoxifen sensitivity", 
#    type = "lines", legend_gp = gpar(col = c("darkblue","darkred")))
EnrichedHeatmap(mat0,   top_annotation = HeatmapAnnotation(line = anno_enriched(gp = gpar(col = c("darkblue","darkred")))), use_raster = TRUE, split = tssDE$type, col = c("white", "darkblue"),name ="All Gene Promoters",column_title = "All Gene Promoters")

dev.off()




#volcano plot
thresholdDE<-repeats_Lx9CVB4_vs_CVB4$FDR>=0
repeats_Lx9CVB4_vs_CVB4$threshold <- thresholdDE 
repeats_Lx9CVB4_vs_CVB4_ordered <- repeats_Lx9CVB4_vs_CVB4[order(repeats_Lx9CVB4_vs_CVB4$FDR), ] 
repeats_Lx9CVB4_vs_CVB4_ordered$genelabels <- ""
repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[1:50] <- rownames(repeats_Lx9CVB4_vs_CVB4_ordered)[1:50]
#repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[repeats_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf(paste0(imageDir,"/volcano_repeats_logFCs_Lx9CVB4_vs_CVB4.pdf"),width=12,height=8)
ggplot(repeats_Lx9CVB4_vs_CVB4_ordered) +
  geom_point_rast(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("Influence of lx9 knockout on CVB4 response Repeats") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  xlim(-12,12) +
  ylim(0,30) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()






