# Manisha Munasinghe - Last Updated: 9/17/20
# Script used for:
#   - Differential Expression Analysis
# For More Details See:
#   - deseq2_personal_func.R
#   - 


.libPaths('/workdir/mam737/Rlibs')
#Set this to wherever you store your Rlibs#

required_packages <- c("geneplotter","ggplot2","plyr","dplyr","DESeq2","gplots",
	"RColorBrewer","stringr","topGO","genefilter","biomaRt","EDASeq",
	"fdrtool","limma","reshape2","pheatmap","TxDb.Dmelanogaster.UCSC.dm6.ensGene",
	"gridExtra","goseq","sva","treemap","RcisTarget")
lapply(required_packages,library,character.only=TRUE)

#sessionInfo() 
# Output Is Listed in R.sessionInfo.txt


#Seed set for reproducibility
set.seed(1994)

## Relevant Directories ##
#setwd('/workdir/mam737/mitoY_rnaseq_final/')
setwd('/fs/cbsuclarkfs1/storage/mam737/mitoY_rnaseq_final')
source('./deseq2_personal_func.R')
fastqDir <- file.path('./raw_files')
htseqDir <- file.path('./htseq_counts')
## Relevant Directories ##

############################### CREATE METADATA ###############################
fastq <- list.files(fastqDir, pattern = '*.fastq')
sampleNO <- str_sub(fastq, 39, 43)
mito <- c(rep("B38",12),rep("B39",12),rep("I02",12),rep("N01",12),rep("N02",12), rep("N23",12))
Ycsome <- c(rep(c("B04","B04","B11","B11","N03","N03","N07","N07","ZH23","ZH23","ZW139","ZW139"),6))
mito_Y <- paste(mito,"-",Ycsome,sep='')
rep <- c(rep(c("A","B"),36))
batch <- c(rep('B1',36),rep(c("B1","B2"),18))
libraryName = paste(mito,"-",Ycsome,"-",sampleNO,sep='')
countFile = list.files(htseqDir, pattern = '*.counts')

metadata <- data.frame(sampleNO=sampleNO, batch=batch,mito=mito, Ycsome=Ycsome, mito_Y=mito_Y, rep=rep,libraryName=libraryName, fastq = fastq, countFile=countFile)
metadata <- lapply(metadata, as.character)
############################### CREATE METADATA ###############################

############################### CREATE DESEQ2 OBJ #############################
## Create DESeq2 Count Table ##
sampleTable <- data.frame(sampleName = metadata$libraryName, fileName = metadata$countFile, batch=metadata$batch,mito= metadata$mito, Ycsome = metadata$Ycsome, mito_Y = metadata$mito_Y,rep=metadata$rep ,sampleNO=metadata$sampleNO,fastq = metadata$fastq)
DESeq2Table <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = htseqDir, design = ~ batch + mito + Ycsome)
print("Number of Genes from HTSEQ Count")
rowData(DESeq2Table)
############################### CREATE DESEQ2 OBJ #############################

############################### QC & FILTERING ################################
raw_counts <- counts(DESeq2Table)##Raw Counts Matrix
#Remove rows where all values are zero
idx.nz <- apply(raw_counts[,-1],1,function(x) !all(x==0))
print("Number of genes with non-zero counts across samples")
sum(idx.nz)
DESeq2Table<- DESeq2Table[idx.nz,]

### Normalization Method
DESeq2Table <- estimateSizeFactors(DESeq2Table)
#sizeFactors(DESeq2Table)
norm_counts <- counts(DESeq2Table, normalized = TRUE)#<- normalized data
norm_counts2 <- norm_counts[rowSums(norm_counts > 20) >=36,]
DESeq2Table  <- DESeq2Table[rowSums(norm_counts > 20) >=36,]
dim(norm_counts2)

#If you want to visualize the impact of
# normalizing run the following chunk of code below

#norm_plots <- generate_normalization_plots(raw_counts,norm_counts,norm_counts2)
#pdf("./Rplots/final/normalization_plots.pdf")
#for (i in seq(1, length(norm_plots), 1)) {
#	grid.arrange(grobs=norm_plots[i],ncol=1,nrow=1)
#}
#dev.off()

############################### QC & FILTERING ################################

############################### Incorporate SVA ###############################

## SVA - For Y & Mito Tests
mod1 <- model.matrix(~batch + mito + Ycsome, colData(DESeq2Table))
mod0 <- model.matrix(~1, colData(DESeq2Table))

svseq <- svaseq(norm_counts2,mod1,mod0)
sva_DESeq2Table <- DESeq2Table
sva_DESeq2Table$SV1 <- svseq$sv[,1]
sva_DESeq2Table$SV2 <- svseq$sv[,2]
sva_DESeq2Table$SV3 <- svseq$sv[,3]
sva_DESeq2Table$SV4 <- svseq$sv[,4]
sva_DESeq2Table$SV5 <- svseq$sv[,5]
sva_DESeq2Table$SV6 <- svseq$sv[,6]
sva_DESeq2Table$SV7 <- svseq$sv[,7]
sva_DESeq2Table$SV8 <- svseq$sv[,8]
design(sva_DESeq2Table) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + batch+ mito + Ycsome

## SVA - For Mito:Y Test

mY_DESeq2Table <- DESeq2Table
design(mY_DESeq2Table) <- ~ batch + mito + Ycsome + mito:Ycsome
mod1 <- model.matrix(~batch + mito + Ycsome + mito:Ycsome, colData(mY_DESeq2Table))
mod0 <- model.matrix(~1, colData(mY_DESeq2Table))
mY_svseq <- svaseq(norm_counts2,mod1,mod0)
mY_sva_DESeq2Table <- mY_DESeq2Table
mY_sva_DESeq2Table$SV1 <- mY_svseq$sv[,1]
mY_sva_DESeq2Table$SV2 <- mY_svseq$sv[,2]
mY_sva_DESeq2Table$SV3 <- mY_svseq$sv[,3]
mY_sva_DESeq2Table$SV4 <- mY_svseq$sv[,4]
mY_sva_DESeq2Table$SV5 <- mY_svseq$sv[,5]
mY_sva_DESeq2Table$SV6 <- mY_svseq$sv[,6]
design(mY_sva_DESeq2Table) <- ~ SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + batch + mito + Ycsome + mito:Ycsome


#If you want to visualize results of PCA analysis
# run the following

rld_blind_SVA <- rlog(sva_DESeq2Table, blind=TRUE)

pca_plots <- generate_pca_plots(rld_blind_SVA)

# Regress out surrogate variables from data to see impact on PCA

sv_mat <- cbind(svseq$sv[,1],svseq$sv[,2],svseq$sv[,3],svseq$sv[,4],svseq$sv[,5],svseq$sv[,6],svseq$sv[,7],svseq$sv[,8])
regress_out_svs_mat <- removeBatchEffect(assay(rld_blind_SVA),covariates=sv_mat)
rv <- rowVars(regress_out_svs_mat)
select <- order(rv, decreasing=TRUE)[seq_len(min(500,length(rv)))]
pca <-prcomp(t(regress_out_svs_mat[select,]))
percentVar <- round(100*(pca$sdev^2/sum(pca$sdev^2)))
d <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2],PC3=pca$x[,3],PC4=pca$x[,4],PC5=pca$x[,5],PC6=pca$x[,6],mito=mito,Ycsome=Ycsome,percentVar=percentVar)
pca_plots[['SVA_adj']] <- ggplot(d, aes(x=PC1,y=PC2, color=mito,shape=Ycsome))+geom_point(size=3)+xlab(paste0("PC1: ",percentVar[1],"% Variance"))+ylab(paste0("PC2: ", percentVar[2],"% Variance"))+ggtitle('PCA: Mito-Y Surrogate Variables Regressed Out')+coord_fixed()+scale_color_manual(values=c("red4","red2","darkorange3",'darkorange','darkgoldenrod1','gold1'))

#pdf('./Rplots/final/pca_plots.pdf')
#for (i in seq(1, length(pca_plots), 1)) {
#	grid.arrange(grobs=pca_plots[i],ncol=1,nrow=1)
#}
#dev.off()

############################### Incorporate SVA ###############################

############################### RUN DESeq2 DE Analysis ########################

# Test for differences between the Y haplotypes
# design = ~ SV1+...+SV8 + batch + mito + Ycsome vs ~ SV1+...+SV8 + batch + mito

# RUN DESeq to test for DE
sva_dds_Y <- DESeq(sva_DESeq2Table,test='LRT',reduced=~SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + batch + mito)
# Extract the results w/ Multiple Test Correction FDR = 0.05
sva_filtered_res_Y <- results(sva_dds_Y,alpha=0.05)
# Order the results by Adjusted p-value
sva_ordered_filtered_res_Y <- sva_filtered_res_Y[order(sva_filtered_res_Y$padj),]
# Prints the Number of Significant Hits
sum(sva_ordered_filtered_res_Y$padj < 0.05,  na.rm=TRUE)
# Subset the data to only include significant hits
sva_Y_sig_hits_padj <- sva_ordered_filtered_res_Y[which(sva_ordered_filtered_res_Y$padj < 0.05),]

#Output results of Y DE analysis rank ordered by p-value:
#write.csv(data.frame(ID=rownames(sva_ordered_filtered_res_Y),raw_pvalue = sva_ordered_filtered_res_Y$pvalue,adj_pvalue = sva_ordered_filtered_res_Y$padj),'./Rplots/final/ordered_Y_results.csv')

# Test for differences between the mito haplogroups
# design = ~ SV1+...+SV8 + batch + mito + Ycsome vs ~ SV1+...+SV8 + batch + Ycsome

# RUN DESeq to test for DE
sva_dds_mito <- DESeq(sva_DESeq2Table,test='LRT',reduced=~SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + SV7 + SV8 + batch + Ycsome)
# Extract the results w/ Multiple Test Correction FDR = 0.05
sva_filtered_res_mito <- results(sva_dds_mito,alpha=0.05)
# Order the results by Adjusted p-value
sva_ordered_filtered_res_mito <- sva_filtered_res_mito[order(sva_filtered_res_mito$padj),]
# Prints the Number of Significant Hits
sum(sva_filtered_res_mito$padj < 0.05,  na.rm=TRUE)
# Subset the data to only include significant hits
sva_mito_sig_hits_padj <- sva_filtered_res_mito[which(sva_filtered_res_mito$padj < 0.05),]

#Output results of mito DE analysis rank ordered by p-value:
#write.csv(data.frame(ID=rownames(sva_ordered_filtered_res_mito),raw_pvalue = sva_ordered_filtered_res_mito$pvalue,adj_pvalue = sva_ordered_filtered_res_mito$padj),'./Rplots/final/ordered_mito_results.csv')


# Test for differences between the mito:Y haplotypes
# design = ~ SV1+...+SV6 + batch + mito + Ycsome vs ~ SV1+...+SV6 + batch + mito + Ycsome

# RUN DESeq to test for DE
sva_dds_mitoY <- DESeq(mY_sva_DESeq2Table,test='LRT',reduced=~SV1 + SV2 + SV3 + SV4 + SV5 + SV6 + batch + mito + Ycsome)
# Extract the results w/ Multiple Test Correction FDR = 0.05
sva_filtered_res_mitoY <- results(sva_dds_mitoY,alpha=0.05)
# Order the results by Adjusted p-value
sva_ordered_filtered_res_mitoY <- sva_filtered_res_mitoY[order(sva_filtered_res_mitoY$padj),]
# Prints the Number of Significant Hits
sum(sva_filtered_res_mitoY$padj < 0.05,  na.rm=TRUE)
# Subset the data to only include significant hits
sva_mitoY_sig_hits_padj <- sva_filtered_res_mitoY[which(sva_filtered_res_mitoY$padj < 0.05),]

#Output results of mito:Y DE analysis rank ordered by p-value:
#write.csv(data.frame(ID=rownames(sva_ordered_filtered_res_mitoY),raw_pvalue = sva_ordered_filtered_res_mitoY$pvalue,adj_pvalue = sva_ordered_filtered_res_mitoY$padj),'./Rplots/final/ordered_mitoY_results.csv')

############################### RUN DESeq2 DE Analysis ########################

############################### RUN GO Analysis ###############################

#Create the GO Database
txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
txsByGene <- transcriptsBy(txdb,'gene')
lengthData <- median(width(txsByGene))
dmel_ensembl_vBDGP6 = useDataset("dmelanogaster_gene_ensembl",mart=useMart ("ensembl"))

# Test for enriched GO terms among Y sensitive hits
# Extract all over- and under-represented GO Cats
sva_Y_go <- go_analysis(sva_dds_Y,dmel_ensembl_vBDGP6)
# Create data.frame of just enriched GO Cats pvals
enriched_Y_go.df <- data.frame(category=sva_Y_go$category,over_represented_pvalue=sva_Y_go$over_represented_pvalue,numDEInCat=sva_Y_go$numDEInCat,numInCat=sva_Y_go$numInCat,term=sva_Y_go$term)
enriched_Y_go.df$over_represented_padj <- p.adjust(enriched_Y_go.df$over_represented_pvalue,method='BH')
# Extract significant GO cats
eY_go_plot.df <- data.frame(category=enriched_Y_go.df[enriched_Y_go.df$over_represented_padj < 0.05,]$category,term=enriched_Y_go.df[enriched_Y_go.df$over_represented_padj < 0.05,]$term,numDEInCat=enriched_Y_go.df[enriched_Y_go.df$over_represented_padj < 0.05,]$numDEInCat,numInCat=enriched_Y_go.df[enriched_Y_go.df$over_represented_padj < 0.05,]$numInCat,padj=enriched_Y_go.df[enriched_Y_go.df$over_represented_padj < 0.05,]$over_represented_padj)

# Create data.frame of just underrep Go cats
underrep_Y_go.df <- data.frame(category=sva_Y_go$category,under_represented_pvalue=sva_Y_go$under_represented_pvalue,numDEInCat=sva_Y_go$numDEInCat,numInCat=sva_Y_go$numInCat,term=sva_Y_go$term)
underrep_Y_go.df$under_represented_padj <- p.adjust(underrep_Y_go.df$under_represented_pvalue,method='BH')
# Extract significant GO cats
uY_go_plot.df <- data.frame(category=underrep_Y_go.df[underrep_Y_go.df$under_represented_padj < 0.05,]$category,term=underrep_Y_go.df[underrep_Y_go.df$under_represented_padj < 0.05,]$term,numDEInCat=(-1*(underrep_Y_go.df[underrep_Y_go.df$under_represented_padj < 0.05,]$numDEInCat)),numInCat=underrep_Y_go.df[underrep_Y_go.df$under_represented_padj < 0.05,]$numInCat,padj=underrep_Y_go.df[underrep_Y_go.df$under_represented_padj < 0.05,]$under_represented_padj)

# Merge all significant GO cats into one data.frame
Y_go_plot.df <- rbind(eY_go_plot.df,uY_go_plot.df)
Y_go_plot.df$term <- factor(Y_go_plot.df$term,levels=Y_go_plot.df$term[order(Y_go_plot.df$numDEInCat)])

#write.csv(Y_go_plot.df,"./Rplots/final/go_Y_results.csv")

# Visualize Results in a Lollipop Plot
#pdf('./Rplots/final/Ysva_golol.pdf')
#ggplot(Y_go_plot.df,aes(x=reorder(term,numDEInCat),y=numDEInCat,label=numDEInCat))+geom_segment(aes(x=term,xend=term,y=0,yend=numDEInCat),color='grey')+geom_point(color='darkorange3',size=5.5)+theme_light() +geom_text(colour='white',size=2)+ coord_flip() + labs(title='Significant GO Categories',x='Number of DE Genes in Category',y='GO term')
#dev.off()


# Test for enriched GO terms among mito sensitive hits
# Extract all over- and under-represented GO Cats
sva_mito_go <- go_analysis(sva_dds_mito,dmel_ensembl_vBDGP6)
# Create data.frame of just enriched GO Cats pvals

enriched_mito_go.df <- data.frame(category=sva_mito_go$category,over_represented_pvalue=sva_mito_go$over_represented_pvalue,numDEInCat=sva_mito_go$numDEInCat,numInCat=sva_mito_go$numInCat,term=sva_mito_go$term)
enriched_mito_go.df$over_represented_padj <- p.adjust(enriched_mito_go.df$over_represented_pvalue,method='BH')
# Extract significant GO cats
emito_go_plot.df <- data.frame(category=enriched_mito_go.df[enriched_mito_go.df$over_represented_padj < 0.05,]$category,term=enriched_mito_go.df[enriched_mito_go.df$over_represented_padj < 0.05,]$term,numDEInCat=enriched_mito_go.df[enriched_mito_go.df$over_represented_padj < 0.05,]$numDEInCat,numInCat=enriched_mito_go.df[enriched_mito_go.df$over_represented_padj < 0.05,]$numInCat,padj=enriched_mito_go.df[enriched_mito_go.df$over_represented_padj < 0.05,]$over_represented_padj)


# Create data.frame of just underrep Go cats
underrep_mito_go.df <- data.frame(category=sva_mito_go$category,under_represented_pvalue=sva_mito_go$under_represented_pvalue,numDEInCat=sva_mito_go$numDEInCat,numInCat=sva_mito_go$numInCat,term=sva_mito_go$term)
underrep_mito_go.df$under_represented_padj <- p.adjust(underrep_mito_go.df$under_represented_pvalue,method='BH')
# Extract significant GO cats
umito_go_plot.df <- data.frame(category=underrep_mito_go.df[underrep_mito_go.df$under_represented_padj < 0.05,]$category,term=underrep_mito_go.df[underrep_mito_go.df$under_represented_padj < 0.05,]$term,numDEInCat=(-1*(underrep_mito_go.df[underrep_mito_go.df$under_represented_padj < 0.05,]$numDEInCat)),numInCat=underrep_mito_go.df[underrep_mito_go.df$under_represented_padj < 0.05,]$numInCat,padj=underrep_mito_go.df[underrep_mito_go.df$under_represented_padj < 0.05,]$under_represented_padj)

# Merge all significant GO cats into one data.frame
mito_go_plot.df <- rbind(emito_go_plot.df,umito_go_plot.df)
mito_go_plot.df$term <- factor(mito_go_plot.df$term,levels=mito_go_plot.df$term[order(mito_go_plot.df$numDEInCat)])
# Rename one cat for plotting purposes
levels(mito_go_plot.df$term)[levels(mito_go_plot.df$term)=='oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen'] <- 'oxidoreductase activity, acting on ...'


#write.csv(mito_go_plot.df,"./Rplots/final/go_mito_results.csv")

# Visualize Results in a Lollipop Plot
#pdf('./Rplots/final/mito1_sva_golol.pdf')
#ggplot(mito_go_plot.df,aes(x=reorder(term,numDEInCat),y=numDEInCat,label=numDEInCat))+geom_segment(aes(x=term,xend=term,y=0,yend=numDEInCat),color='grey')+geom_point(color='dodgerblue4',size=5.5)+theme_light() +geom_text(colour='white',size=2)+ coord_flip() + labs(title='Significant GO Categories',x='Number of DE Genes in Category',y='GO term')
#dev.off()



# Test for enriched GO terms among mito:Y sensitive hits
# Extract all over- and under-represented GO Cats
sva_mitoY_go <- go_analysis(sva_dds_mitoY,dmel_ensembl_vBDGP6)
# Create data.frame of just enriched GO Cats pvals
enriched_mitoY_go.df <- data.frame(category=sva_mitoY_go$category,over_represented_pvalue=sva_mitoY_go$over_represented_pvalue,numDEInCat=sva_mitoY_go$numDEInCat,term=sva_mitoY_go$term)
enriched_mitoY_go.df$over_represented_padj <- p.adjust(enriched_mitoY_go.df$over_represented_pvalue,method='BH')
# Extract significant GO cats
emitoY_go_plot.df <- data.frame(category=enriched_mitoY_go.df[enriched_mitoY_go.df$over_represented_padj < 0.05,]$category,term=enriched_mitoY_go.df[enriched_mitoY_go.df$over_represented_padj < 0.05,]$term,numDEInCat=enriched_mitoY_go.df[enriched_mitoY_go.df$over_represented_padj < 0.05,]$numDEInCat)

# Create data.frame of just underrep Go cats
underrep_mitoY_go.df <- data.frame(category=sva_mitoY_go$category,under_represented_pvalue=sva_mitoY_go$under_represented_pvalue,numDEInCat=sva_mitoY_go$numDEInCat,term=sva_mitoY_go$term)
underrep_mitoY_go.df$under_represented_padj <- p.adjust(underrep_mitoY_go.df$under_represented_pvalue,method='BH')
# Extract significant GO cats
umitoY_go_plot.df <- data.frame(category=underrep_mitoY_go.df[underrep_mitoY_go.df$under_represented_padj < 0.05,]$category,term=underrep_mitoY_go.df[underrep_mitoY_go.df$under_represented_padj < 0.05,]$term,numDEInCat=(-1*(underrep_mitoY_go.df[underrep_mitoY_go.df$under_represented_padj < 0.05,]$numDEInCat)))

# If you want to extract a list of genes within a DE category
#   you can run the following; just adjust appropriately
#	it will output a vector of all genes in that category

#genes_in_dego(dmel_ensembl_vBDGP6,sva_dds_Y,'GO:0007601')

############################### RUN GO Analysis ###############################

############################### RUN Tissue Bias Analysis ######################

# Create the FlyAtlas Tissue Expression Database
FlyAtlas_FPKM <- read.delim("/workdir/mam737/ref_seqs/FlyAtlas/FlyAtlas_FPKM.txt", stringsAsFactors=FALSE)
# Subset as we are only interested in male tissue
male_FlyAtlas_FPKM <- FlyAtlas_FPKM[,grepl("Male",names(FlyAtlas_FPKM))]
male_FlyAtlas_FPKM <- cbind(ID=FlyAtlas_FPKM$ID,male_FlyAtlas_FPKM)
male_FlyAtlas_FPKM <- cbind(male_FlyAtlas_FPKM,sum_fpkm=rowSums(male_FlyAtlas_FPKM[,-1]))

# Test for biased expression among Y sensitive hits
# Extract any DE genes that are biased for expression in the testis
Y_testis_enrich <- tissue_enrich(male_FlyAtlas_FPKM,'Testis',sva_Y_sig_hits_padj)
te_gene_names <- c("CG30487","CG13494","CG1324","Mst98Ca","CG10407","CG11106","CG34241","CG12860","Tengl1","CG11379","Adgf-A2","CG31473","CR43264","CG13843","CG15025","gdl","RpS5b","CG44044")
Y_testis_enrich <- cbind(name=te_gene_names,Y_testis_enrich)

# Visualize Results in Vertical Bar Plot
#pdf("./Rplots/final/Y_testis_bp.pdf")
#ggplot(Y_testis_enrich,aes(x=reorder(name,enrich_val),y=enrich_val,fill=-log10(padj)))+geom_bar(position='dodge',stat='identity')+coord_flip()+scale_fill_gradient(low='darkorange3',high='darkorange')+ylab("Testis Enrichment Value")+xlab("Gene Name") + labs(color="-log10(P-Adj)")
#dev.off()

# Extract any DE genes that are biased for expression in the accessory glands
Y_acp_enrich <- tissue_enrich(male_FlyAtlas_FPKM,'Accessory',sva_Y_sig_hits_padj)
acp_gene_names <- c("CR43978","Dup99B")
Y_acp_enrich <- cbind(name=acp_gene_names,Y_acp_enrich)

# Visualize Results in Vertical Bar Plot
#pdf("./Rplots/final/Y_acp_bp.pdf")
#ggplot(Y_acp_enrich,aes(x=reorder(name,enrich_val),y=enrich_val,fill=-log10(padj)))+geom_bar(position='dodge',stat='identity')+coord_flip()+scale_fill_gradient(low='darkorange3',high='darkorange')+ylab("Accessory Gland Enrichment Value")+xlab("Gene Name") + labs(color="-log10(P-Adj)")
#dev.off()


# Test for biased expression among mito sensitive hits
# Extract any DE genes that are biased for expression in the testis
mito_testis_enrich <- tissue_enrich(male_FlyAtlas_FPKM,'Testis',sva_mito_sig_hits_padj)
te_gene_names <- c("ProtB","CG42841","CG3124","CG43288","CG12853","ssp5","CG12699","CG13326","CR43810","CG3330","Atg8b","CG13442","CG10750","CG15483","CG31294","CG31802","Tim17a1","CG9101","CG4714","CG31029","CG31870","CG7202","CG16782","CG32643","CR44370","CR44308","CG31407","CG17377","CG43059","CG17717","CG4095","CG31538","CG13884","fzr2","CG42288","CG34300","CG34167","CG43659","CG11362","CG1722","CG8043","CG3408","CG16825","Ote","CG30334")
mito_testis_enrich <- cbind(name=te_gene_names,mito_testis_enrich)

# Visualize Results in Vertical Bar Plot
#pdf("./Rplots/final/mito_testis_bp.pdf")
#ggplot(mito_testis_enrich,aes(x=reorder(name,enrich_val),y=enrich_val,fill=-log10(padj)))+geom_bar(position='dodge',stat='identity')+coord_flip()+scale_fill_gradient(low='dodgerblue4',high='deepskyblue2')+ylab("Testis Enrichment Value")+xlab("Gene Name") + labs(color="-log10(P-Adj)")
#dev.off()

# Extract any DE genes that are biased for expression in the accessory glands
mito_acp_enrich <- tissue_enrich(male_FlyAtlas_FPKM,'Accessory',sva_mito_sig_hits_padj)
acp_gene_names <- c("lectin-46Ca","Sems","CR43633","CG14034","Acp53C14a","CG9029","CR43715","Acp33A","CG44477","CG43407","CR43994","CR44805","CG43829","CR44744","CG42869","CG43788","CG44139","CR45800","CaBP1","sel","Gp93","CG2918")
mito_acp_enrich <- cbind(name=acp_gene_names,mito_acp_enrich)

# Visualize Results in Vertical Bar Plot
#pdf("./Rplots/final/mito_acp_bp.pdf")
#ggplot(mito_acp_enrich,aes(x=reorder(name,enrich_val),y=enrich_val,fill=-log10(padj)))+geom_bar(position='dodge',stat='identity')+coord_flip()+scale_fill_gradient(low='dodgerblue4',high='deepskyblue2')+ylab("Accessory Gland Enrichment Value")+xlab("Gene Name") + labs(color="-log10(P-Adj)")
#dev.off()


# Test for biased expression among mito:Y sensitive hits
# Extract any DE genes that are biased for expression in the testis
mitoY_testis_enrich <- tissue_enrich(male_FlyAtlas_FPKM,'Testis',sva_mitoY_sig_hits_padj)

# Extract any DE genes that are biased for expression in the accessory glands
mitoY_acp_enrich <- tissue_enrich(male_FlyAtlas_FPKM,'Accessory',sva_mitoY_sig_hits_padj)
acp_gene_names <- c("CR43404")
mitoY_acp_enrich <- cbind(name=acp_gene_names,mitoY_acp_enrich)

############################### RUN Tissue Bias Analysis ######################

############################### RUN TF Binding Motif Analysis ######################

### Need to turn rownames of significant hits into gene symbol
### geneSet <- c("")

motifRankings <-importRankings("./dm6-5kb-upstream-full-tx-11species.mc8nr.feather")
data(motifAnnotations_dmel)

Y_sighits_names <- rownames(sva_Y_sig_hits_padj)
mito_sighits_names <- rownames(sva_mito_sig_hits_padj)
mitoY_sighits_names <- rownames(sva_mitoY_sig_hits_padj)

Y_geneSet <-
mito_GeneSet <-
mitoY_GeneSet <-

Y_geneLists <- list(geneSetname=Y_geneSet)
Y_motifs_AUC <- calcAUC(Y_geneLists,motifRankings)
Y_motifEnrichmentTable <- addMotifAnnotation(Y_motifs_AUC,motifAnnot=motifAnnotations_dmel)
T_TFbm_results <- motifEnrichmentTable[,-"TF_lowConf",with=FALSE]


mito_geneLists <- list(geneSetname=mito_geneSet)
mito_motifs_AUC <- calcAUC(mito_geneLists,motifRankings)
mito_motifEnrichmentTable <- addMotifAnnotation(mito_motifs_AUC,motifAnnot=motifAnnotations_dmel)
mito_TFbm_results <- motifEnrichmentTable[,-"TF_lowConf",with=FALSE]

mitoY_geneLists <- list(geneSetname=mitoY_geneSet)
mitoY_motifs_AUC <- calcAUC(mitoY_geneLists,motifRankings)
mitoY_motifEnrichmentTable <- addMotifAnnotation(mitoY_motifs_AUC,motifAnnot=motifAnnotations_dmel)
mitoT_TFbm_results <- motifEnrichmentTable[,-"TF_lowConf",with=FALSE]



############################### Result Visualizations #########################

### Count Plots ###
# These are the 3 plots in Figure 3 of the paper
#  it can be adjusted to any genes of your interest
#  just pass the FB gene ID, FB gene name, and whether you want it oriented
#  such that Y haplotype is on the x (Ycsome) or the mito haplogroup (mito)

plotlist <- list()

plotlist[[1]] <- sva_reg_joint_count_plots(regress_out_svs_mat,'FBgn0250832','Dup99B','Ycsome')
plotlist[[2]] <- sva_reg_joint_count_plots(regress_out_svs_mat,'FBgn0013675','mt:COII','mito')
plotlist[[3]] <- sva_reg_joint_count_plots(regress_out_svs_mat,'FBgn0004623','GB76C','Ycsome')

#pdf('./Rplots/final/fig3_counts.pdf',width=50,height=25)
#do.call("grid.arrange",c(plotlist,ncol=3))
#dev.off()

### Count Plots ###

### Heatmaps ###
# This will generate a heatmap for a collection of genes
# You must pass it a 2 vectors, one of the genes of interests FBgn IDS and one of the gene names
# NOTE: The heatmap simply takes the rld transformed counts, not the counts with the SVs regressed out

# The following commands are for all figures in the paper


mitoY_genes <- c("FBgn0000053","FBgn0028396","FBgn0029507","FBgn0030310","FBgn0031176","FBgn0032382","FBgn0032587","FBgn0033301","FBgn0033782","FBgn0033926","FBgn0033981","FBgn0035343","FBgn0035355","FBgn0035715","FBgn0035770","FBgn0036024","FBgn0036619","FBgn0036857","FBgn0038858","FBgn0040001","FBgn0041094","FBgn0041184","FBgn0044812","FBgn0051075","FBgn0051974","FBgn0085195","FBgn0086254","FBgn0259951","FBgn0263323")
mitoY_gene_names <- c("ade3","TotA","Tsp42Ed","PGRP-SA","CG1678","Mal-B2","CG5953","CG12780","sug","Arc1","Cyp6a21","CG16762","CG16985","CG10103","pst","CG18180","Cpr72Ec","CG9629","CG5793","FASN3","scyl","Socs36E","TotC","CG31075","CG31974","CG34166","CG6084","Sfp24Ba","CR43404")


pdf('./Rplots/final/fig4_mitoY_heatmap.pdf',width=48, height=32)
generate_heatmap(mitoY_genes,mitoY_gene_names,rld_blind_SVA,"mito")
dev.off()

metabolic_genes <- c("FBgn0002719","FBgn0001248","FBgn0004057","FBgn0002570","FBgn0002569","FBgn0002571","FBgn0033294","FBgn0033296","FBgn0033297","FBgn0032382")
metabolic_gene_names <- c("Men","Idh","Zw","Mal-A1","Mal-A2","Mal-A3","Mal-A4","Mal-A7","Mal-A8","Mal-B2")


pdf('./Rplots/final/metabolic_heatmap.pdf',width=48, height=32)
generate_heatmap(metabolic_genes,metabolic_gene_names,rld_blind_SVA,"mito")
dev.off()

vis_genes <- c("FBgn0000120","FBgn0000121","FBgn0002936","FBgn0003250","FBgn0004623","FBgn0267252")
vis_gene_names <- c("Arr1","Arr2","ninaA","Rh4","Gbeta76C","Ggamma30A")


pdf('./Rplots/final/vision_heatmap.pdf',width=48, height=32)
generate_heatmap(vis_genes,vis_gene_names,rld_blind_SVA,"Ycsome")
dev.off()