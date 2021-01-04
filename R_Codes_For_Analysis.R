#=======================================================================================
#Batch Normalization
#=======================================================================================
source("https://bioconductor.org/biocLite.R")
biocLite("sva")
biocLite("limma")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("sva", version = "3.8")


library(sva)
library(limma)
setwd("")  #set work dictionary

rt=read.table("merge2.txt",sep="\t",header=T,check.names=F) 
#  "merge2.txt" represents the expression matrix merged from TPM matrix (obtained from shell script using fastq data from SRP118628/GSE104140)
# and expression matrix from GSE28829. The "VLOOKUP" function in excel was used to annoate probe matrix of GSE28829 and merge expression matrices
# of these 2 datasets to get "merge2.txt" ("merge()" function in R can also be used)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)  
batchType=c(rep(1,29),rep(2,32))  
modType=c(rep("Advanced",16),rep("Eealy",13),rep("Advanced",19),rep("Early",13)) 
#batchType anD modType represents batch information.
mod = model.matrix(~as.factor(modType))
outTab=ComBat(data, batchType, mod, par.prior=TRUE)
outTab=rbind(geneNames=colnames(outTab),outTab)
write.table(outTab,file="normalize.txt",sep="\t",quote=F,col.names=F)





#=======================================================================================
#DEG analysis
#=======================================================================================
source("http://bioconductor.org/biocLite.R")
biocLite("limma")

logFoldChange=1
adjustP=0.05

library(limma)
setwd("")
rt=read.table("Final_Gene_Matrix.txt",sep="\t",header=T,check.names=F,row.names=1)
#After the "normalize.txt" was obtained from batch normalization, we changed its file name into "Final_Gene_Matrix.txt"
Targets=read.table("Final_Clinical_Info.txt",sep="\t",header=T,row.names=1,check.names=F)
#"Final_Clinical_Info.txt" include "Sample_ID" (Begin with GSM or SRR) and "Disease_Group" (Grouping information) obtained from GEO (For GSE28829) and SRA database (For SRP118628).
# The grouping information labeled each sample with advanced plaque(1) or early plaque(0) for further analysis.


#rt=as.matrix(rt)
#rownames(rt)=rt[,1]
#exp=rt[,2:ncol(rt)]
#dimnames=list(rownames(exp),colnames(exp))
#rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames) 

#differential
targets <-Targets
design <- model.matrix(~Disease_Group,data=targets)
colnames(design)

fit <- lmFit(rt,design)
fit <- eBayes(fit)


allDiff=topTable(fit,coef="Disease_Group",adjust='fdr',number=200000)
write.table(allDiff,file="limmaTab.xls",sep="\t",quote=F)

#write table
diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffSig,file="diff.xls",sep="\t",quote=F)
diffUp <- allDiff[with(allDiff, (logFC>logFoldChange & adj.P.Val < adjustP )), ]
write.table(diffUp,file="up.xls",sep="\t",quote=F)
diffDown <- allDiff[with(allDiff, (logFC<(-logFoldChange) & adj.P.Val < adjustP )), ]
write.table(diffDown,file="down.xls",sep="\t",quote=F)

#write expression level of diff gene
hmExp=rt[as.vector(diffSig[,1]),]
diffExp=rbind(id=colnames(hmExp),hmExp)
write.table(diffExp,file="diffExp.txt",sep="\t",quote=F,col.names=F)

#volcano
pdf(file="vol.pdf")
xMax=max(-log10(allDiff$adj.P.Val))
yMax=max(abs(allDiff$logFC))
plot(-log10(allDiff$adj.P.Val), allDiff$logFC, xlab="-log10(adj.P.Val)",ylab="logFC",
     main="Volcano", xlim=c(0,xMax),ylim=c(-yMax,yMax),yaxs="i",pch=20, cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC>logFoldChange)
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="red",cex=0.8)
diffSub=subset(allDiff, adj.P.Val<adjustP & logFC<(-logFoldChange))
points(-log10(diffSub$adj.P.Val), diffSub$logFC, pch=20, col="green",cex=0.8)
abline(h=0,lty=2,lwd=3)
dev.off()





#=======================================================================================
#Heatmap
#=======================================================================================
install.packages("pheatmap")
install.packages("Rcpp")
setwd("")     

rt=read.table("New_Data_Matrix.txt",sep="\t",header=T,check.names=F)
#We extracted the expression data of DEGs using "VLOOKUP" function in excel the matrix with DEG expression data were named "New_Data_Matrix.txt")
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
#rt=log2(rt+1)
#rt[rt>15]=15

library(pheatmap)
Type=c(rep("Advanced Plaque",35),rep("Early Plaque",26)) 
names(Type)=colnames(rt)
ann=cbind(Type)
ann=as.data.frame(ann)

tiff(file="heatmap.tiff",
     width = 20,            
     height =25,            
     units ="cm",
     compression="lzw",
     bg="white",
     res=600)
pheatmap(rt, annotation=ann, 
         color = colorRampPalette(c("green", "black", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         scale="row",
         show_rownames = F,
         show_colnames = F,
         fontsize_row=0.000001,
         fontsize_col=0.000001)
dev.off()





#=======================================================================================
#Functional enrichment analysis
#=======================================================================================
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

#Preprocessing
setwd("")
Gene_For_Annotation <- read.table("diff_list.txt",sep="\t",header=T,check.names=F)
#"diff_list.txt" had only one column "Gene", which contained the gene symbols for DEGs. In this analysis, the functional
# enrichment analysis for genes of WGCNA gene module was conducted. The gene list was replaced by gene symbols of each module.
Symbol <- Gene_For_Annotation$Gene
EntrezID <- bitr(Symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
EntrezID <- EntrezID$ENTREZID
head(EntrezID)

#Gene Ontology

#GO analysis
GO <- enrichGO(gene = EntrezID,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable = T)
write.table(GO,file="GO.txt",sep="\t",quote=F,row.names = F)

#Barplot
tiff(file="GO_barplot.tiff",width = 26,height = 20,units ="cm",compression="lzw",bg="white",res=600)
barplot(GO, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

#Dotplot
tiff(file="GO_dotplot.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
dotplot(GO, showCategory = 20)
dev.off()

#KEGG pathway

#KEGG pathway analysis
KEGG_Pathway <- enrichKEGG(gene = EntrezID, 
                           organism = "hsa", 
                           pvalueCutoff =0.05, 
                           qvalueCutoff =0.05)
write.table(KEGG_Pathway,file="KEGG_Pathway.txt",sep="\t",quote=F,row.names = F)

#Barplot
tiff(file="KEGG_barplot.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
barplot(KEGG_Pathway, drop = TRUE, showCategory = 20)
dev.off()

#Dotplot
tiff(file="KEGG_dotplot.tiff",width = 20,height = 12,units ="cm",compression="lzw",bg="white",res=600)
dotplot(KEGG_Pathway, showCategory = 20)
dev.off()





#=======================================================================================
#WGCNA
#=======================================================================================
source("http://bioconductor.org/biocLite.R")
biocLite("WGCNA")
biocLite("stringr")
install.packages("WGCNA")
install.packages("stringr")
install.packages(c("AnnotationDbi","impute","GO.db","preprocessCore"))
install.packages("GO.db")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("preprocessCore")


#=========================================================
setwd('')
# The following setting is important, do not omit. 
options(stringsAsFactors = FALSE); 
expro=read.csv('Advanced_Plaque.txt',sep = '\t',row.names = 1) 
#"Advanced_Plaque.txt" was obtained by extracting expression data of advanced plaque group from "Final_Gene_Matrix.txt"
dim(expro)
#m.vars=apply(expro,1,var) 
#expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4]),]
#dim(expro.upper) 
#write.table(expro.upper,file="geneInput_variancetop0.25.txt",sep='\t ',quote=F,row.names=T)
datExpr0=as.data.frame(t(expro)); 
library(WGCNA)
gsg = goodSamplesGenes(datExpr0, verbose = 3); 
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed: 
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[! gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[! gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data: datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


#=========================================================
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9); 
par(cex = 0.45);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 60, col = "red"); 
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 57, minSize = 10) 
table(clust)

# clust 1 contains the samples we want to keep. 
keepSamples = (clust==1)  
datExpr = datExpr0[keepSamples, ] 
dim(datExpr0)
dim(datExpr)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


#=========================================================
# Re-cluster samples
dim(datExpr)
sampleTree2 = hclust(dist(datExpr), method = "average")
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9); 
par(cex = 0.45);
par(mar = c(0,4,2,0))
plot(sampleTree2, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)


#=========================================================
# Load the WGCNA package 
library(WGCNA)
# The following setting is important, do not omit. 
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above. 
enableWGCNAThreads()
# Load the data saved in the first part  
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables.
lnames 


#=========================================================
# Choose a set of soft-thresholding powers  
powers = c(c(1:10), seq(from = 12, to=20, by=2)) # Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5) 
par(mfrow = c(1,2)); 
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft- thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[, 2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[, 2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h 
abline(h=0.85,col="red") 
Spower=16
# Mean connectivity as a function of the soft-thresholding power 
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# here we define the adjacency matrix using soft thresholding with beta=6
ADJ1=abs(cor(datExpr,use="p"))^Spower
# When you have relatively few genes (<5000) use the following code 
k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code 
k=softConnectivity(datE=datExpr,power=Spower)
# Plot a histogram of k and a scale free topology plot 
sizeGrWindow(10,5)
par(mfrow=c(1,2)) 
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")


#=========================================================
softPower = Spower;
adjacency = adjacency(datExpr, power = softPower)


#=========================================================
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


#=========================================================
geneTree = hclust(as.dist(dissTOM), method = "average"); 
# Plot the resulting clustering tree (dendrogram) sizeGrWindow(12,12)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);


#=========================================================
minModuleSize = 100;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro =FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods)
dynamicColors = labels2colors(dynamicMods) 
table(dynamicColors)
sizeGrWindow(8,12)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#=========================================================
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes # Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average"); 
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes", 
     xlab = "", sub = "")


#=========================================================
MEDissThres = 0.25 
# abline=0.25 
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors 
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


#=========================================================
sizeGrWindow(12, 9)
#pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6) 
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
#dev.off()


#=========================================================
# Rename to moduleColors 
moduleColors = mergedColors
colorOrder = c("grey", standardColors(50)); 
moduleLabels = match(moduleColors, colorOrder)-1; 
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts 
save(MEs, moduleLabels, moduleColors, geneTree, file = "G-02-networkConstruction-StepByStep.RData")


#=========================================================
# Define numbers of genes and samples 
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
# names (colors) of the modules 
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p")); 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep=""); 
names(MMPvalue) = paste("p.MM", modNames, sep="");


#=========================================================
# Create the starting data frame
geneInfo0 = data.frame(geneSymbol=rownames(geneModuleMembership),                 
                       moduleColor = moduleColors, 
                       geneModuleMembership, 
                       MMPvalue)
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor);
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")


#=========================================================
# Select module 
table(moduleColors)
module = "yellow";
# Select module probes 
probes = names(datExpr)
inModule = (moduleColors==module); 
modProbes = probes[inModule];
IMConn = softConnectivity(datExpr[, modProbes],power=9); 
dat1=datExpr[inModule]
datExp_IMConn <-data.frame(IMConn,t(dat1)) 
datExp_IMConn=data.frame(datExp_IMConn)
write.table(datExp_IMConn,
            file =paste("Intramodule_connectivity-",module," .txt"),sep='\t')
nTop = 30;
top = (rank(-IMConn) <= nTop)
dat2=t(datExp_IMConn) 
dat2<-data.frame(dat2) 
dat3<-dat2[top]
dat3<-t(dat3)
dat3<-data.frame(dat3) 
write.table(dat3,
            file = paste("Intramodule_connectivity-",module,"- top30.txt"),sep='\t')


#=========================================================
setwd('');  # Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit. 
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part 
lnames = load(file = "GBM-01-dataInput.RData");
lnames = load(file = "G-01-dataInput.RData");
#The variable lnames contains the names of loaded variables. 
lnames
# Load network data saved in the second part.
lnames = load(file = "GBM-02-networkConstruction-StepByStep.RData"); 
lnames = load(file = "G-02-networkConstruction-StepByStep.RData"); 
lnames


#=========================================================
# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power = Spower); 
table(moduleColors)
# Select mbodules
modules = c("yellow"); 
# Select module probes probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules)); 
modProbes = probes[inModule];
modGenes = modProbes;
# Select the corresponding Topological Overlap 
modTOM = TOM[inModule, inModule]; 
dimnames(modTOM) = list(modProbes, modProbes)


#=========================================================
  nTop = 10000;
IMConn = softConnectivity(datExpr[, modProbes]); 
top = (rank(-IMConn) <= nTop)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM[top, top],
                               edgeFile = paste("CytoscapeInput- edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput- nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE, threshold = 0.02, nodeNames = modProbes, altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])





#=======================================================================================
#Lasso regression and LDA analysis
#=======================================================================================
library(car) #package to calculate Variance Inflation Factor
library(corrplot) #correlation plots
library(leaps) #best subsets regression
library(glmnet) #allows ridge regression, LASSO and elastic net
library(caret) #this will help identify the appropriate parameters
carotictr <- read.table("Gene100.csv",header=TRUE,row.names="Tag",sep=",") 
## "gene100.csv" is the table which contains the expression of most related 100 genes in MCODE and can be extracted from Final_Gene_Matrix.txt.
carotrain <- data.frame(carotictr)
x <- as.matrix(carotrain[, 2:100])
y <- carotrain[, 1]

###lasso

lasso <- glmnet(x, y, family = "binomial", alpha = 1)
print(lasso)

plot(lasso, xvar = "lambda", label = TRUE)
lasso.coef <- predict(lasso, s = 0.0041, type = "coefficients")
lasso.coef


#### lasso.cv
set.seed(317)lasso.cv = cv.glmnet(x, y, nfolds = 10)
plot(lasso.cv)
lasso.cv$lambda.min #minimum
lasso.cv$lambda.1se #one standard error away
coef(lasso.cv, s = "lambda.min")
##the minimum of the lambda was selected
#####
traintrian <- read.table("gene15.csv",header=TRUE,row.names="Tag",sep=",") 
## "gene15.csv" containts the expression of 15 genes that were selected by LvASSO and can also be extracted from Final_Gene_Matrix.txt
trian2<-data.frame(traintrian) 
train=sample(1:nrow(trian2)*3/4)
train_set = carotrain[train,]
test_set = carotrain[-(train),]
##3/4 samples was included in the train group and 1/4 samples were in test group
###
library(MASS)
lda.fit<-lda(Disease_Group~. ,data=train_set)
plot(lda.fit)
lda.pred<-predict(lda.fit,test_set)
confusionMatrix<-table(lda.pred$class,test_set$Disease_Group)
confusionMatrix
mean(test_set$Disease_Group!=lda.pred$class)
###LDA was used to verify the classification

install.packages(nnet)
require(nnet)
multi.fitting<-multinom(Disease_Group~.,data=train_set)
logistics.pred = predict(multi.fitting,newdata=test_set,"class")
##
confusionMatrix<-table(logistics.pred,test_set$Disease_Group)
confusionMatrix
mean(test_set$Disease_Group!=logistics.pred)
## the model was verified in the test group