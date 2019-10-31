
library(GenomicRanges)
library(ShortRead)
#library(R.utils)
#biocLite("BSgenome.Mmusculus.UCSC.mm10")
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
library(RUVSeq)
library(org.Mm.eg.db)
library(DESeq)
#library(Go.db)
library(ggrepel)
library(ggrastr)
timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84



#homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../../.."

inPath=paste0(homedir,"/projects/claudia/project_results/rsem/")


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

ids<-read.table("ids.txt")
newNames=ids[match(names(geneResults),ids$V1),"V2"]
#newNames=newNames[!is.na(newNames)]

geneResults=geneResults[7:21]
newNames=newNames[7:21]
names(geneResults)=newNames
geneResults=geneResults[c(1:6,10:15)]
df<-as.data.frame(geneResults)
row.names(df)<-c(as.character(data$gene_id),as.character(dataS$gene_id))
colnames(df)[10:12]=c("LX9_CVB4_1","LX9_CVB4_2","LX9_CVB4_3")

#########################################
########## 1. Load in the files
#########################################



cols<-rep(brewer.pal(6,"Set1"),each=3)
pdf("../project_results/figures/samples_pca.pdf",width=12,height=8)
pca<-princomp(df)
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, colnames(df),pos = 1)
dev.off()

pdf("../project_results/figures/samples_pca_components.pdf",width=12,height=8)
plot(pca)
dev.off()


#trying out mds
group=gl(4, 3, labels = c("B6", "B6_CVB4", "LX9", "LX9_CVB4"))
y <- DGEList(df,group=group)
y <- calcNormFactors(y)

# without labels, indexes of samples are plotted.
col <- as.numeric(group)
pdf("../project_results/figures/samples_MDS.pdf",width=12,height=8)
mds <- plotMDS(y, top=200, col=col,labels=group)
dev.off()
#Annotate with symbols, aggregate

#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]
df$EntrezGene[c(dim(df)[1]-76):dim(df)[1]]<-as.character(dataS$gene_id)

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]
df$symbol[c(dim(df)[1]-76):dim(df)[1]]<-as.character(dataS$gene_id)

#eliminate duplicated symbols
o <- order(rowSums(df[,1:12]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$symbol)
df<-df[!d,]

#add the gene length
#data$gene_id=gsub("\\..*","",data$gene_id)
#df_merged<-merge(df,data[,c(1,4)],by.x=0,by.y="gene_id")
#df<-df_merged
#row.names(df)<-df[,1]
#df<-df[,-1]
#eliminate lowly expressed
include<-apply(df[,1:12],1,function(x){sum(x>=5)>=2})
df<-df[include,]
include<-apply(df[,1:12],1,function(x){sum(x==0)>=3})
df<-df[!include,]

df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

#RUVseq
#x <- as.factor(rep(c("B6", "B6_CVB4", "LX9", "LX9_CVB4"), each=3))
#set <- newSeqExpressionSet(as.matrix(df[,1:12]),
#                           phenoData = data.frame(x, row.names=colnames(df[,1:12])))
#set
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#EDASeq::plotPCA(set, col=colors[x], cex=1.2)
#
#set <- betweenLaneNormalization(set, which="full")
#plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#EDASeq::plotPCA(set, col=colors[x], cex=1.2)
#
##use sequins to normalize
#spikes<-as.character(row.names(df)[grepl("R1|R2_",row.names(df))])
#set1 <- RUVg(set, spikes, k=1)
#pData(set1)
#plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
#EDASeq::plotPCA(set1, col=colors[x], cex=1.2)
#doesn't work as well


#change the id to entrez
#o <- order(rowSums(df[,1:12]), decreasing=TRUE)
#df=df[o,]
#d<-duplicated(df$EntrezGene)
#df<-df[!d,]
#row.names(df)<-df$EntrezGene


#provide the EC numbers
#egENZYME <- toTable(org.Mm.egENZYME)
#m <- match(row.names(df), egENZYME$gene_id)
#
#df[,"EC number"]<-egENZYME$ec_number[m]

cols<-rep(brewer.pal(6,"Paired"),each=3)
pdf("../project_results/figures/pca_clean_all.pdf",width=12,height=8)
pca<-princomp(df[1:12])
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, colnames(df)[1:12],pos = 1)
dev.off()

pdf("../project_results/figures/pca_components_clean_all.pdf",width=12,height=8)
plot(pca)
dev.off()

temp=df
#Differential expression of LX9 vs wt
dfL<-df[,c(1:3,7:9)]
treatment <- rep(c("LX9","WT"), each=3)
design <- model.matrix(~treatment)
expr <- DGEList(counts=dfL[1:6])
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_lx9wt <- glmLRT(fit)
logFCs_lx9wt<-as.data.frame(topTags(lrt_lx9wt,n=20000))
logFCs_lx9wt=logFCs_lx9wt[!grepl("^R1_",row.names(logFCs_lx9wt)),]
logFCs_lx9wt=logFCs_lx9wt[!grepl("^R2_",row.names(logFCs_lx9wt)),]

write.csv(logFCs_lx9wt,"logFCs_FC_pvalue.csv")

plus<-logFCs_lx9wt[logFCs_lx9wt$logFC<=-1,]
minus<-logFCs_lx9wt[logFCs_lx9wt$logFC>=1,]
significantIDs_lx9wt<-c(row.names(plus[1:100,]),row.names(minus[1:100,]))
fittedValues<-fitted.values(lrt_lx9wt)
fittedValues100_lx9wt<-fittedValues[row.names(fittedValues) %in% significantIDs_lx9wt,]

write.csv(fittedValues,"lx9wt_all.csv")
write.csv(fittedValues100_lx9,"lx9wt_topUpDown200.csv")

fittedValues100_lx9wt=fittedValues100_lx9wt[c(1:50,151:200),]
pdf("heatmap_lx9wt.pdf",width=8,height=8)
pheatmap(log10(fittedValues100_lx9wt),cex=0.7)
dev.off()

fittedValues100_lx9wt<-fittedValues[row.names(fittedValues) %in% immune_genes,]

pdf("heatmap_lx9wt_immune_genes.pdf",width=8,height=8)
pheatmap(log10(fittedValues100_lx9wt),cex=0.7)
dev.off()



#First do the differential expression of LX9

dfL<-df[,7:12]
treatment <- rep(c("CVB4","WT"), each=3)
#type <- c(rep(c("B6", "LX9"), each=6))
#treatment <- rep(rep(c("WT", "CVB4"), each=3),2)
#type <- c(rep(c("B6", "LX9"), each=6))

#design <- model.matrix(~0.+treatment)
design <- model.matrix(~treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=dfL[1:6])
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_lx9 <- glmLRT(fit)



logFCs<-as.data.frame(topTags(lrt_lx9,n=20000))
write.csv(logFCs,"lx9_all_FC_pvalue.csv")

#results<-list()
#for(i in 1:10){logFCS_short=logFCs[logFCs$FDR<10^-i,];results[[i]]=dim(logFCS_short[logFCS_short$logFC>=1,])[1]}
#plot(1:10,unlist(results))
#volcano plot
thresholdDE<-logFCs$FDR<1e-40
logFCs$threshold <- thresholdDE 

pdf("lx9_vs_CVB4.pdf",width=12,height=8)
ggplot(logFCs) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("Influence of CVB4 on lx9 knockouts") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  
dev.off()



logFCs<-logFCs[logFCs$FDR<0.00001,]
plus<-logFCs[logFCs$logFC<=-1,]
minus<-logFCs[logFCs$logFC>=1,]
significantIDs_lx9<-c(row.names(plus[1:100,]),row.names(minus[1:100,]))
fittedValues<-fitted.values(lrt_lx9)
fittedValues100_lx9<-fittedValues[row.names(fittedValues) %in% significantIDs_lx9,]

write.csv(fittedValues,"lx9_all.csv")
write.csv(fittedValues100_lx9,"lx9_topUpDown200.csv")


######### same for B6

dfB<-df[,1:6]
treatment <- rep(c("CVB4","WT"), each=3)
#type <- c(rep(c("B6", "LX9"), each=6))
#treatment <- rep(rep(c("WT", "CVB4"), each=3),2)
#type <- c(rep(c("B6", "LX9"), each=6))

#design <- model.matrix(~0.+treatment)
design <- model.matrix(~treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=dfB[1:6])
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_B6 <- glmLRT(fit)


logFCs_B6<-as.data.frame(topTags(lrt_B6,n=20000))
write.csv(logFCs_B6,"B6_all_FC_pvalue.csv")

thresholdDE<-logFCs_B6$FDR<1e-40
logFCs_B6$threshold <- thresholdDE 
ggplot(logFCs_B6) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("Influence of CVB4 on lx9 knockouts") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  #scale_y_continuous(limits = c(0,50)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  




logFCs_B6<-logFCs_B6[logFCs_B6$FDR<0.00001,]
plus<-logFCs_B6[logFCs_B6$logFC<=-1,]
minus<-logFCs_B6[logFCs_B6$logFC>=1,]
significantIDs_B6<-c(row.names(plus[1:100,]),row.names(minus[1:100,]))
fittedValues_B6<-fitted.values(lrt_B6)
fittedValues100_B6<-fittedValues[row.names(fittedValues) %in% significantIDs_B6,]
significantIDs_B6[significantIDs_B6 %in% significantIDs_B6]


write.csv(fittedValues_B6,"B6_all.csv")
write.csv(fittedValues100_B6,"B6_topUpDown200.csv")

#effect of deletion of LINE on CVB4 infection

dfLINE<-df[,c(1:12)]
#treatment <- rep(c("CVB4", "lx9_CVB4"), each=3)
type <- c(rep(c("B6","LX9"), each=6))
treatment <- rep(rep(c("CVB4","WT"), each=3),2)
#type <- c(rep(c("B6", "LX9"), each=6))

#design <- model.matrix(~0.+treatment)
design <- model.matrix(~type+type:treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=dfLINE)
expr <- calcNormFactors(expr,method="TMM")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_LINE <- glmLRT(fit)


logFCs_Lx9CVB4_vs_CVB4<-as.data.frame(topTags(lrt_LINE,n=40000))
write.csv(logFCs_Lx9CVB4_vs_CVB4,"Lx9CVB4_vs_CVB4_pvalue.csv")

logFCs_Lx9CVB4_vs_CVB4<-logFCs_Lx9CVB4_vs_CVB4[logFCs_Lx9CVB4_vs_CVB4$FDR<0.00001,]
plus<-logFCs_Lx9CVB4_vs_CVB4[logFCs_Lx9CVB4_vs_CVB4$logFC<=-1,]
minus<-logFCs_Lx9CVB4_vs_CVB4[logFCs_Lx9CVB4_vs_CVB4$logFC>=1,]
significantIDs_LINE<-c(row.names(plus[1:100,]),row.names(minus[1:100,]))
fittedValues_LINE<-fitted.values(lrt_LINE)
fittedValues100_LINE<-fittedValues_LINE[row.names(fittedValues_LINE) %in% significantIDs_LINE,]
#significantIDs_LINE[significantIDs_LINE %in% significantIDs_LINE]

fittedValues100_LINE=fittedValues100_LINE[c(1:50,151:200),]

pheatmap(log10(fittedValues100_LINE),cex=0.7)


write.csv(fittedValues_LINE,"Lx9CVB4_vs_CVB4_all.csv")
write.csv(fittedValues100_LINE,"Lx9CVB4_vs_CVB4_topUpDown200.csv")








##change the id to entrez
o <- order(rowSums(df[,1:12]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$EntrezGene)
df<-df[!d,]
row.names(df)<-df$EntrezGene

dfL<-df[,c(1:3,7:9)]
treatment <- rep(c("LX9","WT"), each=3)
design <- model.matrix(~treatment)
expr <- DGEList(counts=dfL[1:6])
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_lx9wt <- glmLRT(fit)



#treatment <- rep(c("WT", "CVB4"), each=3)
#design <- model.matrix(~treatment)
#
#dfL<-df[,7:12]
#expr <- DGEList(counts=dfL[1:6])
#expr <- calcNormFactors(expr,method="RLE")
#expr <- estimateDisp(expr,design)
#fit <- glmFit(expr,design)
#lrt_lx9 <- glmLRT(fit)
#
#
#dfB<-df[,1:6]
#expr <- DGEList(counts=dfB[1:6])
#expr <- calcNormFactors(expr,method="RLE")
#expr <- estimateDisp(expr,design)
#fit <- glmFit(expr,design)
#lrt_B6 <- glmLRT(fit)
#
#
dfLINE<-df[,c(1:12)]
type <- c(rep(c("B6","LX9"), each=6))
treatment <- rep(rep(c("CVB4","WT"), each=3),2)
design <- model.matrix(~type+type:treatment)

expr <- DGEList(counts=dfLINE)
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_LINE <- glmLRT(fit)
#
##lx9wt
goLx9wt <- goana(lrt_lx9wt, species="Mm")
topGO(goLx9wt)
write.table(topGO(goLx9wt,number=Inf),"GO_categories_lx9wt.txt",row.names=T,quote=F,sep="\t")

kegLx9wt <- kegga(lrt_lx9wt, species="Mm")
topKEGG(kegLx9wt)
write.table(topKEGG(kegLx9wt,number=Inf),"KEGG_categories_lx9wt.txt",row.names=T,quote=F,sep="\t")

###lx9
#goLx9 <- goana(lrt_lx9, species="Mm")
#topGO(goLx9)
#write.table(topGO(goLx9,number=Inf),"GO_categories_lx9.txt",row.names=T,quote=F,sep="\t")
#
#kegLx9 <- kegga(lrt_lx9, species="Mm")
#topKEGG(kegLx9)
#write.table(topKEGG(kegLx9,number=Inf),"KEGG_categories_lx9.txt",row.names=T,quote=F,sep="\t")
#
##B6
#goB6 <- goana(lrt_B6, species="Mm")
#topGO(goB6)
#write.table(topGO(goB6,number=Inf),"GO_categories_B6.txt",row.names=T,quote=F,sep="\t")
#
#kegB6 <- kegga(lrt_B6, species="Mm")
#topKEGG(kegB6)
#write.table(topKEGG(kegB6,number=Inf),"KEGG_categories_B6.txt",row.names=T,quote=F,sep="\t")
#
#
##LINE
goLINE <- goana(lrt_LINE, species="Mm")
topGO(goLINE)
write.table(topGO(goLINE,number=Inf),"GO_categories_Lx9CVB4_vs_CVB4.txt",row.names=T,quote=F,sep="\t")

kegLINE <- kegga(lrt_LINE, species="Mm")
topKEGG(kegLINE)
write.table(topKEGG(kegLINE,number=Inf),"KEGG_categories_Lx9CVB4_vs_CVB4.txt",row.names=T,quote=F,sep="\t")
#
#
#
#
#
#volcano plot
thresholdDE<-logFCs_Lx9CVB4_vs_CVB4$FDR>=0
logFCs_Lx9CVB4_vs_CVB4$threshold <- thresholdDE 
logFCs_Lx9CVB4_vs_CVB4_ordered <- logFCs_Lx9CVB4_vs_CVB4[order(logFCs_Lx9CVB4_vs_CVB4$FDR), ] 
logFCs_Lx9CVB4_vs_CVB4_ordered$genelabels <- ""
logFCs_Lx9CVB4_vs_CVB4_ordered$genelabels[1:50] <- rownames(logFCs_Lx9CVB4_vs_CVB4_ordered)[1:50]
logFCs_Lx9CVB4_vs_CVB4_ordered$genelabels[logFCs_Lx9CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf("volcano_logFCs_Lx9CVB4_vs_CVB4.pdf",width=12,height=8)
ggplot(logFCs_Lx9CVB4_vs_CVB4_ordered) +
  geom_point_rast(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("Influence of lx9 knockout on CVB4 response") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  xlim(-12,12) +
  ylim(0,70) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()

immune_genes<-logFCs_Lx9CVB4_vs_CVB4_ordered$genelabels
thresholdDE<-logFCs_lx9wt$FDR>=0
logFCs_lx9wt$threshold <- thresholdDE 
logFCs_lx9wt_ordered <- logFCs_lx9wt[order(logFCs_lx9wt$FDR), ] 
logFCs_lx9wt_ordered$genelabels <- ""

logFCs_lx9wt_ordered$genelabels[row.names(logFCs_lx9wt_ordered) %in% immune_genes] = row.names(logFCs_lx9wt_ordered)[row.names(logFCs_lx9wt_ordered) %in% immune_genes]

#logFCs_lx9wt_ordered$genelabels[1:100] <- rownames(logFCs_lx9wt_ordered)[1:100]
#logFCs_lx9wt_ordered$genelabels[logFCs_lx9wt_ordered$logFC<0]<-""
pdf("volcano_logFCs_Lx9_vs_wt_selectedGenes.pdf",width=12,height=8)
ggplot(logFCs_lx9wt_ordered) +
  geom_point(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("Influence of lx9 knockout on gene expression") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  xlim(-12,12) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25))) 
dev.off()



########## expression of BL6 inf vs uninf ############


dfLINE<-df[,c(1:6)]
treatment <- rep(c("CVB4","WT"), each=3)
design <- model.matrix(~treatment)

expr <- DGEList(counts=dfLINE)
expr <- calcNormFactors(expr,method="RLE")
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt_LINE <- glmLRT(fit)
#

logFCs_B6CVB4_vs_CVB4<-as.data.frame(topTags(lrt_LINE,n=40000))
write.csv(logFCs_B6CVB4_vs_CVB4,"B6CVB4_vs_CVB4_pvalue.csv")

logFCs_B6CVB4_vs_CVB4<-logFCs_B6CVB4_vs_CVB4[logFCs_B6CVB4_vs_CVB4$FDR<0.0001,]
plus<-logFCs_B6CVB4_vs_CVB4[logFCs_B6CVB4_vs_CVB4$logFC<=-1,]
minus<-logFCs_B6CVB4_vs_CVB4[logFCs_B6CVB4_vs_CVB4$logFC>=1,]
significantIDs_LINE<-c(row.names(plus[1:100,]),row.names(minus[1:100,]))
fittedValues_LINE<-fitted.values(lrt_LINE)
fittedValues100_LINE<-fittedValues_LINE[row.names(fittedValues_LINE) %in% significantIDs_LINE,]
#significantIDs_LINE[significantIDs_LINE %in% significantIDs_LINE]

#volcano plot
thresholdDE<-logFCs_B6CVB4_vs_CVB4$FDR>=0
logFCs_B6CVB4_vs_CVB4$threshold <- thresholdDE 
logFCs_B6CVB4_vs_CVB4_ordered <- logFCs_B6CVB4_vs_CVB4[order(logFCs_B6CVB4_vs_CVB4$FDR), ] 
logFCs_B6CVB4_vs_CVB4_ordered$genelabels <- ""
logFCs_B6CVB4_vs_CVB4_ordered$genelabels[1:50] <- rownames(logFCs_B6CVB4_vs_CVB4_ordered)[1:50]
#logFCs_B6CVB4_vs_CVB4_ordered$genelabels[logFCs_B6CVB4_vs_CVB4_ordered$logFC<0]<-""
pdf("volcano_logFCs_B6CVB4_vs_CVB4.pdf",width=12,height=8)
ggplot(logFCs_B6CVB4_vs_CVB4_ordered) +
  geom_point_rast(aes(x=logFC, y=-log10(FDR), colour=threshold)) +
  ggtitle("B6 response to CVB4") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  geom_text_repel(aes(x = logFC, y = -log10(FDR), label = genelabels)) +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))  

dev.off()

goLINE <- goana(lrt_LINE, species="Mm")
topGO(goLINE)
write.table(topGO(goLINE,number=Inf),"GO_categories_B6CVB4_vs_CVB4.txt",row.names=T,quote=F,sep="\t")

kegLINE <- kegga(lrt_LINE, species="Mm")
topKEGG(kegLINE)
write.table(topKEGG(kegLINE,number=Inf),"KEGG_categories_B6CVB4_vs_CVB4.txt",row.names=T,quote=F,sep="\t")

