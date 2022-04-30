library(pcadapt)
library(ggplot2)
library(qqman)
library(plyr)
library(dplyr)
library(RColorBrewer)
library(stringr)
library(gridExtra)

path_to_file <- "C:\\Users\\eviec\\OneDrive\\Desktop\\Plink\\GSE76645_INDONESIA_SUMBA_TIMOR_mind05_noTIMOR.bed" 
data <- read.pcadapt(path_to_file, type = "bed") 

x <- pcadapt(input = data, K = 20)
plot(x, option="screeplot")

setwd("C:/Users/eviec/OneDrive/Documents/Archaeology/RGenetics")
villages=read.csv("villages.csv", header=FALSE) #Read in list of ID's and villages


#First PCA
poplist.names=villages$V4 
plot(x, option = "scores", pop = poplist.names)
#ggsave(file="PCA1.tiff")
dev.off()


#Read in new data after removing related individuals - KING trim
path_to_file <- "C:\\Users\\eviec\\OneDrive\\Documents\\Archaeology\\RGenetics\\GSE76645_INDONESIA_SUMBA_TIMOR_mind05_noTIMOR_norelatives.bed" 
data <- read.pcadapt(path_to_file, type = "bed") #new data: no Timor individuals and no relatives


#Change to K=6 after initial score plot
x2 <- pcadapt(input = data, K = 6) 
plot(x2, option="screeplot")
#ggsave(file="Scree_1_a.tiff", width=5, height=5, dpi=300)
dev.off()

#Changing 'villages' list to match trimmed data
#Took the output list of individuals to keep from king, put it in .csv to read in
villages_unrelated_king=read.csv("villages_after_king_1.csv", header=FALSE)
summary(villages_unrelated_king) #210 (correct number after removing 7 Timor individuals)
villages_after_king=subset(villages,villages$V2%in%villages_unrelated_king$V1)
summary(villages_after_king) #210
poplist.names=villages_after_king$V4 #for PCA


#Second PCA - after removing relatives
#singular.values is #"the vector containing the K ordered squared root of the proportion of variance explained by each PC"
variance.explained=(x2$singular.values)^2 
variance.explained #PC1:  0.009640499, PC2: 0.009036479
sum(variance.explained)#0.05093707 approx. 5% altogether
plot(x2, option = "scores", pop = poplist.names)+
  labs(title="Sumba Principle Component Analysis", x="PC1 (0.96%)", y="PC2 (0.9%)", color="Villages")+
  scale_color_brewer(palette="Dark2")
#ggsave(file="PCA_1_c.tiff", width=7, height=5, dpi=300)
dev.off()

#PC's 3 and 4
plot(x2, option = "scores", i=3, j=4, pop = poplist.names)+
  labs(title="Sumba Principle Component Analysis", x="PC3 (0.85%)", y="PC4 (0.82%)", color="Villages")+
  scale_color_brewer(palette="Dark2")
#ggsave(file="PCA_1_d.tiff", width=7, height=5, dpi=300)
dev.off()


#Benjamini-Hochberg procedure - outliers before linkage disequilibrium clumping 
padj <- p.adjust(x2$pvalues,method="BH") 
alpha <- 0.1
outliers <- which(padj < alpha)
length(outliers) #500 - qvalue method had same number of outliers


#Read in SNP data
SNPs=read.table("GSE76645_INDONESIA_SUMBA_TIMOR_mind05_noTIMOR_norelatives.bim", header=FALSE, sep="", quote= " ")
SNPs$index=seq.int(nrow(SNPs))
summary(SNPs) #index to 1:713246


#SNPs before LD
SNPs.b4.LD=data.frame(SNPs)
tracemem(SNPs.b4.LD)==tracemem(SNPs)
SNPs.b4.LD$pval=x2$pvalues
summary(SNPs.b4.LD$pval) 
SNPs.b4.LD$logP=(-log10(SNPs.b4.LD$pval))
SNPs.b4.LD.no.NA=na.exclude(SNPs.b4.LD)
nrow(SNPs.b4.LD.no.NA)#545644


#LD clumping step 1: GGplot of loadings before LD
SNPs.b4.LD$loadings.PCA1=x2$loadings[,1]
#SNPs.b4.LD$V1=as.factor(SNPs.b4.LD$V1) - for chrom colour
p1=ggplot(SNPs.b4.LD, aes(index,loadings.PCA1))+
  geom_point(size=0.8)+ #aes(color=V1)
  scale_y_continuous(limits=c(-0.01,0.01))+
  labs(x="Index", y="Loadings PC1") #color="Chromosome"
#ggsave(filename="Loadings_before_LD.tiff", width=7, height=5, units="in",dpi=300)
dev.off()


#LD clumping step 2: Remove SNP's with r2 higher than 0.1 in blocks of 200 
res <- pcadapt(data, K = 6, LD.clumping = list(size = 200, thr = 0.1)) 
(res$singular.values)^2
plot(res, option = "screeplot")
#ggsave(file="Scree_2.tiff", width=5, height=5)
dev.off()


#GGplot of loadings after clumping
SNPs$loadings.PCA1=res$loadings[,1]
p2=ggplot(SNPs, aes(index,loadings.PCA1))+
  geom_point(size=0.8)
scale_y_continuous(limits=c(-0.025,0.025))+
  labs(x="Index", y="Loadings PC1")
#ggsave(filename="Loadings_after_LD.tiff", width=7, height=5, units="in",dpi=300)
dev.off()

#GGplot of loadings after clumping in colour
SNPs$loadings.PCA1=res$loadings[,1]
SNPs$V1=as.factor(SNPs$V1) # convert to factor for colour
p3=ggplot(SNPs, aes(index,loadings.PCA1))+
  geom_point(size=0.8)+ aes(color=V1)+
  scale_y_continuous(limits=c(-0.025,0.025))+
  labs(x="Index", y="Loadings PC1", color="Chromosome")+
  guides(colour = guide_legend(override.aes = list(size=4)))+
  scale_color_manual(values=c("#ec3db6","#54e06e","#1e46c2", "#a8c100", "#9c77ff","#019e31", "#e80750",  "#007629",   "#ff9cee", "#a08b00", "#0055b0","#fb5c33",  "#01a89c", "#902477", "#a9d380", "#92303d", "#98baff","#9a5500","#7d4b71", "#e5c261","#c27c91", "#638151","#ff9167"), labels=c(1:22,"X"))
#ggsave(filename="Loadings_after_LD_colour_1.tiff", width=7, height=5, units="in",dpi=300)
dev.off()
p3

grid.arrange(p1, p2, nrow=2)
g=arrangeGrob(p1,p2,ncol=2)
#ggsave(filename="Both_loadings_col.tiff", g, width=14, height=6, units="in",dpi=300)
dev.off()

grid.arrange(p1, p3, nrow=2)
g2=arrangeGrob(p1,p3,nrow=2)
#ggsave(filename="Both_loadings_row_colour.tiff", g2, width=14, height=15, units="in",dpi=300)
dev.off()


#Benjamini-Hochberg - how many outliers AFTER LD clumping
padj <- p.adjust(res$pvalues,method="BH") #length(padj) 713246
alpha <- 0.1
outliers.after.LD <- which(padj < alpha) #returns position of SNPs in the bed/bim files
length(outliers.after.LD) #Benjamini used as more conservative, 401 outliers when qq gave 500


#Manhattan plot: highlighing outliers after LD
SNPs$pval=res$pvalues
summary(SNPs$pval)#167602 NA's: P-values with maf<0.05%
SNPs$Outliers=(SNPs$index%in%outliers.after.LD) #add TRUE/FALSE outlier column to SNP table
SNPs$logP=(-log10(SNPs$pval))#column added to SNPs NOT SNP.outliers
SNPs.no.NA=na.exclude(SNPs)#exclude NA's 
nrow(SNPs.no.NA)#69252
summary(SNPs)
summary(SNPs.no.NA)
SNP.outliers=subset(SNPs, SNPs$index%in%outliers.after.LD) #subset of outlier SNPs
highlight.outliers=SNP.outliers$V2
summary(SNP.outliers)


#Manhattan without outliers highlighted
ggplot(SNPs.no.NA, aes(V4, logP))+
  geom_point(size=0.8, shape=20)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x="Chromosome", y="-log10 P-value", title="Manhattan Plot")


#Manhattan with outliers highlighted
cols1=c("TRUE"="#E7298A", "FALSE"="#1B9E77")
ggplot(SNPs.no.NA, aes(V4, logP))+
  geom_point(size=0.8, shape=20, aes(colour=outliers))+
  facet_grid(.~V1, scales="free_x", switch="x")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(x="Chromosome", y="-log10 P-value", title="Manhattan Plot", colour="SNPs")+
  scale_color_manual(labels=c("Outlier", "Inlier"), values=cols1)
#ggsave("Manhattan_after_LD_1_c.tiff", width=7, height=5, units="in",dpi=300)
dev.off()



#Read in gene list
genelist=read.table("BiomartProteinCodingGeneList.bed", header=FALSE, sep="", quote=" ")
genelist=as.data.frame(genelist) #20327
summary(genelist)


#Function: Find SNPs within genes, and the nearest upstream/downstream gene for each SNP
in_genes_x=(function(x){
  x[12]=NA
  x[13]=NA
  x[14]=NA
  x[15]=NA
  x[16]=NA
  SNP_pos=as.numeric(x[4][[1]]) #as.numeric to make sure its treating cols numerically
  genelist_chr=subset(genelist, str_trim(x[1][[1]])==(str_remove(genelist[,1], "chr")))
  test=subset(genelist_chr, SNP_pos>=as.numeric(genelist_chr[,2]) & SNP_pos<=as.numeric(genelist_chr[,3]))
  if(nrow(test)>0){
    
    x[12]=as.matrix(test)[1,4] #SNPs in genes 
  }
  test=subset(genelist_chr, SNP_pos>as.numeric(genelist_chr[,3])) #Nearest upstream gene
  if(nrow(test)>0){
    x[13]=as.matrix(test)[nrow(test),4]
    x[14]=SNP_pos-as.numeric(as.matrix(test)[nrow(test),3])
  }
  test=subset(genelist_chr, SNP_pos<as.numeric(genelist_chr[,2])) #Nearest downstream gene
  if(nrow(test)>0){
    x[15]=as.matrix(test)[1,4]
    x[16]=as.numeric(as.matrix(test)[1,2])-SNP_pos
  }
  return(x)
})


#Output: 
out=apply(X=SNP.outliers, FUN=in_genes_x, MARGIN=1)#apply with Margin=1 to go row-wise
out1=matrix(out, nrow=ncol(out), ncol=nrow(out), byrow=TRUE)
out1
colnames(out1)=c("Chromosome", "rsID", "V3", "Position", "V5", "V6", "Index", "PC1 loadings","Pval", "Outlier","-log10 Pval" ,"Gene", "GeneLH", "LH Gap", "GeneRH", "RH Gap")
head(out1)
#write.csv(out1, file="SNPs and genes.csv")


which(out1[,12]=="AK2") #will return the index number of SNP if it is an outlier
which(out1[,12]=="SLC4A1")#None in SAO gene
which(out1[,12]=="SCN5A")#None in Brugada gene
unique(out1[,12])#116 genes with outlier SNPs inside them
as.data.frame(out1)


out1[order(as.numeric(out1[,9])),] 

#All SNPs
#write.csv(out1[order(as.numeric(out1[,9])),], file="SNPs_and_genes.csv")

#Top 10 smallest p-value SNPs overall 
smallest.p.SNPs=out1[order(as.numeric(out1[,9])),]
smallest.p.SNPs[1:10,]
SNP.outliers[order(SNP.outliers[,9]),] #check
tab.2=smallest.p.SNPs[1:10,c(1,2,4,9,12,13,14,15,16)]
tab.2
#write.csv(tab.2, file="smallest.p.SNPs.csv")


#Top 10 smallest P-value SNPs WITHIN genes
SNPs.in.genes=subset(out1, !is.na(out1[,12]))
summary(SNPs.in.genes)
nrow(SNPs.in.genes) #170 - previously 109
SNPs.in.genes=SNPs.in.genes[order(as.numeric(SNPs.in.genes[,9])),]
SNPs.in.genes[1:10,]
tab.1=SNPs.in.genes[1:10,c(1,2,4,9,12)]
tab.1
#write.csv(tab.1, file="SNPs.in.genes.csv")
#write.csv(SNPs.in.genes, file="All.SNPs.in.genes.csv")


#Genes with more than 1 SNP
dups=subset(data.frame(table(SNPs.in.genes[,12])), Freq>1)
nrow(dups) #31 Genes with more than 1 SNP
table(dups[,2]==2)#15
table(dups[,2]==3)#12
table(dups[,2]==4)#2
table(dups[,2]==5)#1
table(dups[,2]==6)#1
sum(dups[,2])#85 SNPs total
colnames(dups)=c("Gene", "Frequency")
dups=dups[order(dups[,2]),]
#write.csv(dups, file="duplicated.genes.csv")
which(SNPs.in.genes[1:10,12]%in%dups[,1])#1 5 6 7 8
which(smallest.p.SNPs[1:10,12]%in%dups[,1])#2,10 -- AcSS3 (4!) and FAM167A 

