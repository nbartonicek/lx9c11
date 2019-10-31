
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
library(clusterProfiler)
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
timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84


#setwd("~/contrib1/nenbar/projects/claudia/scripts/RNAseq")
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
if(!file.exists(countFile)){
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
} else {load(countFile)}

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
df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol
#eliminate lowly expressed
include<-apply(df[,1:12],1,function(x){sum(x>=10)>=2})
df<-df[include,]
include<-apply(df[,1:12],1,function(x){sum(x==0)>=6})
df<-df[!include,]

df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

temp=df

#effect of infection on B6
dfL<-df[,c(6:1)]
treatment <- rep(c("B6_CVB4","B6"), each=3)
design <- model.matrix(~treatment)
expr <- DGEList(counts=dfL[1:6])
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
significantIDs_B6_CVB4<-c(row.names(plus_B6_CVB4[1:50,]),row.names(minus_B6_CVB4[1:50,]))
fittedValues_B6_CVB4<-fitted.values(lrt_B6_CVB4)
fittedValues100_B6_CVB4<-fittedValues_B6_CVB4[row.names(fittedValues_B6_CVB4) %in% significantIDs_B6_CVB4,]

fittedValues50_B6_CVB4=fittedValues100_B6_CVB4
fittedValues50_B6_CVB4<-fittedValues50_B6_CVB4[,c(6,5,4,3,2,1)]
colsums<-apply(fittedValues50_B6_CVB4[,4:6],1,median)
fittedValues50_B6_CVB4<-fittedValues50_B6_CVB4[rev(order(colsums)),]
pdf("supFig1.pdf",width=8,height=12)
pheatmap(log1p(fittedValues50_B6_CVB4),cex=0.9,cluster_rows = F,cluster_cols=F)
dev.off()
write.csv(fittedValues_B6_CVB4,"B6_CVB4_all_repeats.csv")
write.csv(fittedValues100_B6_CVB4,"B6_CVB4_topUpDown200_repeats.csv")

repeats_B6_CVB4<-logFCs_B6_CVB4[row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat),]
write.csv(repeats_B6_CVB4,"logFC_B6_CVB4_repeats_only.csv")

repeats_B6_CVB4[row.names(repeats_B6_CVB4) %in% c("Lx9","L1_Mus3"),]
fittedValues_B6_CVB4[row.names(fittedValues_B6_CVB4) %in% c("Lx9","L1_Mus3"),]
df[row.names(df) %in% c("Lx9","L1_Mus3"),]

repeats_Lx9CVB4_vs_CVB4[row.names(repeats_Lx9CVB4_vs_CVB4) %in% c("Lx9","L1_Mus3"),]
fittedValues_LINE[row.names(fittedValues_LINE) %in% c("Lx9","L1_Mus3"),]


#volcano plot
genes_B6_CVB4<-logFCs_B6_CVB4[!(row.names(logFCs_B6_CVB4) %in% row.names(dfRepeat)),]
thresholdDE<-genes_B6_CVB4$FDR<=0.001
genes_B6_CVB4$threshold <- thresholdDE 
genes_B6_CVB4_ordered <- genes_B6_CVB4[order(genes_B6_CVB4$FDR), ] 
genes_B6_CVB4_ordered$genelabels <- ""
genes_B6_CVB4_ordered$genelabels[1:100] <- rownames(genes_B6_CVB4_ordered)[1:100]
#repeats_Lx9CVB4_vs_CVB4_ordered$genelabels[repeats_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf(paste0(imageDir,"/volcano_genes_logFCs_B6_CVB4_vs_CVB4_rasterized.pdf"),width=12,height=8)
ggplot(genes_B6_CVB4_ordered) +
  geom_point_rast(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("B6+CVB4 vs CVB4 Genes") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  xlim(-12,12) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()

data=read.table("KEGG_categories_B6CVB4_vs_B6.txt",header=T,sep="\t",stringsAsFactors=F)
dfU<-data[order(data$P.Up),][1:30,]
dfD<-data[order(data$P.Down),][1:10,]
dfU<-dfU[!(dfU$Pathway %in% c("Osteoclast differentiation","Pathways in cancer")),][1:10,]

df<-dfU
df$ratio<-df$Up/df$N
df<-df[,c(1,2,5,7)]
df$P.Up<--log10(df$P.Up)
df$Pathway<-factor(df$Pathway,levels=rev(df$Pathway))
df1<-df
colnames(df1)[3]<-"P"
pdf(paste0(imageDir,"/BL6_CVB4_pathway_final_up.pdf"),width=6,height=6, useDingbats=FALSE)
p<-ggplot(df1,aes(P,Pathway))
p<-p+geom_point(aes(size=N,colour=ratio))
p<-p+xlim(c(0,25))
p<-p+scale_size(limits=c(-100,1000),range=c(1,10),breaks=c(-100,0,25,50,100,200,500,1000))
p<-p+scale_colour_gradient2(limits=c(0.25,0.75),midpoint = 0.5,low="#4363AE",high="#A51F22")
p<-p+theme(panel.background = element_rect(fill = "white"),panel.grid.major = element_line(colour = "gray90",size=0.2))
p
dev.off()

df<-dfD
df$ratio<-df$Up/df$N
df<-df[,c(1,2,6,7)]
df$P.Down<--log10(df$P.Down)
df$Pathway<-factor(df$Pathway,levels=rev(df$Pathway))
df2=df
colnames(df2)[3]<-"P"
pdf(paste0(imageDir,"/BL6_CVB4_pathway_final_down.pdf"),width=5.5,height=6, useDingbats=FALSE)
p<-ggplot(df2,aes(P,Pathway))
p<-p+geom_point(aes(size=N,colour=ratio))
p<-p+xlim(c(0,25))
p<-p+scale_size(limits=c(-100,1000),range=c(1,10),breaks=c(-100,0,25,50,100,200,500,1000))
p<-p+scale_colour_gradient2(limits=c(0,0.75),midpoint = 0.5,low="darkblue",high="darkred")
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


#get sequences of promoters for overexpressed genes
logFCs_B6_CVB4<-read.table("B6_CVB4_pvalue_repeats.csv",sep=",",header=T,row.names=1)
minus_B6_CVB4<-logFCs_B6_CVB4[logFCs_B6_CVB4$logFC>=1,]

genes<-row.names(minus_B6_CVB4)[1:500]
#TSS distribution

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
promoter <- promoters(genes(txdb,columns=c("gene_id")))
gns<-bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
promoterSubset=promoter[promoter$gene_id %in% gns$ENTREZID]

seq = BSgenome::getSeq(BSgenome.Mmusculus.UCSC.mm10, promoterSubset)
names(seq) = paste0("SEQUENCEALL_", seq_along(seq))
Biostrings::writeXStringSet(seq, paste0(projectDir,"/project_results/fasta/B6_CVB4_top500_upregulated_promoters.fa"))



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
pdf(paste0(imageDir,"/volcano_repeats_logFCs_B6_CVB4_vs_CVB4_rasterized_byType.pdf"),width=5,height=2, useDingbats=FALSE)
ggplot(repeats_B6_CVB4_orderedShort) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold),size=0.6) +
  ggtitle("B6+CVB4 vs B6 Repeats") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  facet_wrap(~class,ncol=2)+
  scale_color_manual(values=c("#4363AE","#A51F22"))+
  scale_x_continuous(breaks = round(seq(-7, 7, by = 1),1)) +
  theme_bw(base_size=6) +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = repeatlabels),segment.size = 0.5,min.segment.length=0.1,size=2) 

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
promoterDE<-promoter[names(promoter) %in% row.names(minus_B6_CVB4[1:500,])]
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
#this was changed from 40
df<-df[df$SYMBOL %in% row.names(minus_B6_CVB4[1:500,]),]
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
stat1=PFMatrixList[["MA0144.2"]]

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

dfL<-split(df,df$SYMBOL)
dfMin<-lapply(dfL,function(x){min(x$distanceToTSS)})
names(dfMin)<-names(dfL)
dfMin<-sort(unlist(dfMin))
min25<-names(dfMin)[1:25]
df<-df[df$SYMBOL %in% min25,]
pdf(paste0(imageDir,"/LTR_LINE_distance2genes.pdf"),width=4,height=4, useDingbats=FALSE)
p<-ggplot(df,aes(SYMBOL,distanceToTSS))
p<-p+geom_point(aes(colour=TFscore))
p<-p+coord_flip()
p<-p+scale_color_manual(values=c("gray30", "orange", "red"))
p<-p+theme_bw()
p
dev.off()


#do the enrichment of TFs
promoterDE<-promoter[names(promoter) %in% row.names(minus_B6_CVB4[1:500,])]

#total universe (previously the genome)
chrGR<-promoterDE
strand(chrGR)<-"*"
result<-list()
oddsRatio<-function(region,regionDatabase){

        #for a given repeat in region, what are the odds it has a TF
        totalNuc<-sum(width(region))
        reducedTFDatabase=reduce(regionDatabase)

        #mat<-findOverlaps(region,reducedsnpDatabase)
        #total hits is the length of intersect
        strand(region)="*"
        strand(reducedTFDatabase)="*"

        totalHits<-sum(width(intersect(region,reducedTFDatabase)))

        #totalHits<-length(unique(subjectHits(mat)))
        firstRatio<-totalHits/(totalNuc-totalHits)

        nonRegion<-sum(as.numeric(width(chrGR)))-totalNuc
        genomeHits<-sum(as.numeric(width(reducedTFDatabase)))-totalHits
        secondRatio<-genomeHits/(nonRegion-genomeHits)

        odds<-firstRatio/secondRatio
        se<-sqrt(1/totalNuc+1/(totalNuc-totalHits)+1/(nonRegion-genomeHits)+1/genomeHits)
        #cat(odds)
        #cat("\n")
        #cat(se)
        #cat("\n")
        result<-list()
        result[[1]]=odds
        result[[2]]=se


        #p-value
        SNPenrichment <-
     matrix(c(as.numeric(totalHits), as.numeric(totalNuc), as.numeric(genomeHits), as.numeric(nonRegion)),
            nrow = 2,
            dimnames = list(Regions = c("WithSNP", "Total"),
                            Genome = c("WithSNP", "Total")))
        result[[3]]=chisq.test(SNPenrichment)$p.value

        return(result)

}
repeats=gr
chrGRPlus<-chrGR
strand(chrGRPlus)<-"+"
seqs<-getSeq(Mmusculus, chrGRPlus)

irfGR<-lapply(seqs,function(x){matchPWM(as.matrix(irf1),x)})
statGR<-lapply(seqs,function(x){matchPWM(as.matrix(stat1),x)})

irfGRL<-GRangesList()
statGRL<-GRangesList()

for(gene in names(irfGR)){
  cat(gene)
  startsI<-start(irfGR[[gene]])
  endsI<-end(irfGR[[gene]])
  startsS<-start(statGR[[gene]])
  endsS<-end(statGR[[gene]]) 
  grNewI<-GRanges(seqnames=seqnames(chrGRPlus[gene]),IRanges(start=start(chrGRPlus[gene])+startsI,end=start(chrGRPlus[gene])+endsI))
  irfGRL[[gene]]<-grNewI
  grNewS<-GRanges(seqnames=seqnames(chrGRPlus[gene]),IRanges(start=start(chrGRPlus[gene])+startsS,end=start(chrGRPlus[gene])+endsS))
  statGRL[[gene]]<-grNewS
}

TFs<-reduce(c(unlist(irfGRL),unlist(statGRL)))
repeats<-reduce(gr)
unlist(oddsRatio(repeats,TFs))






























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






