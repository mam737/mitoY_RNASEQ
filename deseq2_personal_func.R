#The following function generates 
#  a series of plots visualizing effects of
#  normalization and filtering thresholds
generate_normalization_plots <- function(raw_counts_matrix,norm_counts_matrix,norm_counts2_matrix) {
	plotList <- list()

	pseudoRawCounts <- log2(raw_counts_matrix +1 )
	melt_raw_counts <- melt(raw_counts_matrix)
	melt_pseudoRawCounts <- melt(pseudoRawCounts)

	pseudoNormCounts <- log2(norm_counts_matrix +1 )
	melt_norm_counts <- melt(norm_counts_matrix)
	melt_pseudoNormCounts <- melt(pseudoNormCounts)

	pseudoNormCounts2 <- log2(norm_counts2_matrix + 1)
	melt_norm_counts2 <- melt(norm_counts2_matrix)
	melt_pseudoNormCounts2 <- melt(pseudoNormCounts2)	

	plotList[['a']] <- ggplot(melt_raw_counts, aes(x=value)) + geom_histogram(breaks=seq(0,400000,by=5000)) + ggtitle("Histogram of Raw HTSEQ Counts") + labs(x='HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['b']] <- ggplot(melt_pseudoRawCounts, aes(x=value)) + geom_histogram(breaks=seq(0,20,by=1)) + ggtitle("Histogram of -log2(Raw HTSEQ Counts)") + labs(x='HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['c']] <- ggplot(melt_norm_counts, aes(x=value)) + geom_histogram(breaks=seq(0,400000,by=10000)) + ggtitle("Histogram of Normalized HTSEQ Counts") + labs(x='Normalized HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['d']] <- ggplot(melt_pseudoNormCounts, aes(x=value)) + geom_histogram(breaks=seq(0,20,by=1)) + ggtitle("Histogram of -log2(Normalized HTSEQ Counts)") + labs(x='Normalzied HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['e']] <- ggplot(melt_pseudoRawCounts, aes(x=Var2,y=value,fill=Var2))+geom_boxplot()+ggtitle('Boxplot Distribution of Normalized Counts in Each Sample')+labs(x='Sample',y='log2(Normalized Count+1)')+scale_color_brewer(palette='Spectral')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='none')
	plotList[['f']] <- ggplot(melt_norm_counts2, aes(x=value)) + geom_histogram(breaks=seq(0,400000,by=10000)) + ggtitle("Histogram of Normalized HTSEQ Counts: At least 36 samples have norm counts > 20") + labs(x='Normalized HTSEQ Counts ',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['g']] <- ggplot(melt_pseudoNormCounts2, aes(x=value)) + geom_histogram(breaks=seq(0,20,by=1)) + ggtitle("Histogram of -log2(Normalized HTSEQ Counts): At least 36 samples have norm counts >20") + labs(x='Normalzied HTSEQ Counts',y='Frequency') + theme(plot.title=element_text(hjust=0.5))
	plotList[['h']] <- ggplot(melt_pseudoNormCounts2, aes(x=Var2,y=value,fill=Var2))+geom_boxplot()+ggtitle('Boxplot Distribution of Normalized Counts in Each Sample: At least 36 samples have norm counts >20 ')+labs(x='Sample',y='log2(Normalized Count+1)')+scale_color_brewer(palette='Spectral')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),legend.position='none')

	return(plotList)
}

generate_pca_plots <- function(rld_object) {
	plotList <- list()

	pcaData <- plotPCA(rld_object, intgroup=c('mito','Ycsome','rep','batch'), returnData=TRUE)
	percentVar <- round(100*attr(pcaData,"percentVar"))

	plotList[['batch']] <- ggplot(pcaData, aes(x=PC1,y=PC2, color=batch))+geom_point(size=3)+xlab(paste0("PC1: ",percentVar[1],"% Variance"))+ylab(paste0("PC2: ", percentVar[2],"% Variance"))+ggtitle('PCA: Batch')+coord_fixed()
	plotList[['mito']] <- ggplot(pcaData, aes(x=PC1,y=PC2, color=mito))+geom_point(size=3)+xlab(paste0("PC1: ",percentVar[1],"% Variance"))+ylab(paste0("PC2: ", percentVar[2],"% Variance"))+ggtitle('PCA: Mito')+coord_fixed()+scale_color_manual(values=c("red4","red2","darkorange3",'darkorange','darkgoldenrod1','gold1'))
	plotList[['Ycsome']] <- ggplot(pcaData, aes(x=PC1,y=PC2, color=Ycsome))+geom_point(size=3)+xlab(paste0("PC1: ",percentVar[1],"% Variance"))+ylab(paste0("PC2: ", percentVar[2],"% Variance"))+ggtitle('PCA: Y')+coord_fixed()+scale_color_manual(values=c("darkgreen","green3","dodgerblue4",'deepskyblue2','purple4','mediumpurple2'))
	plotList[['mitoY']] <- ggplot(pcaData, aes(x=PC1,y=PC2, color=mito,shape=Ycsome))+geom_point(size=3)+xlab(paste0("PC1: ",percentVar[1],"% Variance"))+ylab(paste0("PC2: ", percentVar[2],"% Variance"))+ggtitle('PCA: Mito-Y')+coord_fixed()+scale_color_manual(values=c("red4","red2","darkorange3",'darkorange','darkgoldenrod1','gold1'))

	return(plotList)
}


#The following function runs 
#  a GO analysis
#  to identify over- and under-represented GO terms in the DE results
go_analysis <- function(dds_obj,data_set) {
	fdr.threshold <- 0.05
	go_res <- results(dds_obj,independentFiltering=FALSE)
	assayed.genes <- rownames(go_res)
	de.genes <- rownames(go_res)[which(go_res$padj < fdr.threshold)]
	gene.vectors <- as.integer(assayed.genes%in%de.genes)
	names(gene.vectors) <- assayed.genes

	dmel_ensembl_vBDGP6 = data_set

	present_genes <- subset(lengthData,names(lengthData)%in%assayed.genes)
	missing_genes <- subset(assayed.genes, !(assayed.genes%in%names(lengthData)))
	test <- rep(NA,length(missing_genes))
	names(test) <- missing_genes
	missing_genes <- test

	gene.vectors <- gene.vectors[order(names(gene.vectors))]
	final_length_set <- lengthData[order(names(gene.vectors))]
	GOmap = getBM (filters = "ensembl_gene_id", attributes =c("ensembl_gene_id", "go_id"), values = names(gene.vectors), mart = dmel_ensembl_vBDGP6)

	pwf <- nullp(gene.vectors,bias.data=final_length_set)
	go <- goseq(pwf,gene2cat=GOmap)

	return(go)
}

genes_in_dego <- function(data_set,dds_obj,go_id){
	fdr.threshold <- 0.05
	go_res <- results(dds_obj,independentFiltering=FALSE)
	assayed.genes <- rownames(go_res)
	de.genes <- rownames(go_res)[which(go_res$padj < fdr.threshold)]
	gene.vectors <- as.integer(assayed.genes%in%de.genes)
	names(gene.vectors) <- assayed.genes
	dmel_ensembl_vBDGP6 = data_set

	present_genes <- subset(lengthData,names(lengthData)%in%assayed.genes)
	missing_genes <- subset(assayed.genes, !(assayed.genes%in%names(lengthData)))
	test <- rep(NA,length(missing_genes))
	names(test) <- missing_genes
	missing_genes <- test

	gene.vectors <- gene.vectors[order(names(gene.vectors))]
	final_length_set <- lengthData[order(names(gene.vectors))]
	GOmap = getBM (filters = "ensembl_gene_id", attributes =c("ensembl_gene_id", "go_id"), values = names(gene.vectors), mart = dmel_ensembl_vBDGP6)


	genes_in_go <- GOmap[which(GOmap$go_id==go_id),]
	de_genes_in_go.df <- genes_in_go[which(genes_in_go$ensembl_gene_id %in% de.genes),]
	de_genes.vec <- de_genes_in_go.df$ensembl_gene_id
	return(de_genes.vec)

}

tissue_enrich <- function(FlyAtlas_enrich.df,tissue,dds_res) {
	list_of_sig_genes <- rownames(dds_res)
	padj_info.df <- data.frame(ID=rownames(dds_res),padj=dds_res$padj)

	FlyAtlas_enrich.df <- FlyAtlas_enrich.df[!(FlyAtlas_enrich.df$sum_fpkm==0),]

	tissue_fpkm <- FlyAtlas_enrich.df[,grepl(tissue,names(FlyAtlas_enrich.df))]
	enrich_val <- tissue_fpkm/FlyAtlas_enrich.df$sum_fpkm

	tissue_enrich.df <- data.frame(ID=FlyAtlas_enrich.df$ID,enrich_val=enrich_val)

	sig_genes_enrich <- tissue_enrich.df[which(tissue_enrich.df$ID %in% list_of_sig_genes),]

	

	#Return Top Enriched Genes
	top_sig_genes_enrich <- sig_genes_enrich[sig_genes_enrich$enrich_val > 0.5,]
	top_sig_genes_enrich <- merge(top_sig_genes_enrich,padj_info.df,by='ID')
	top_sig_genes_enrich <- top_sig_genes_enrich[order(-top_sig_genes_enrich$enrich_val),]
	return(top_sig_genes_enrich)
}

# The following function generates
#   a heatmap for a given set of genes 
generate_heatmap <- function(vec_of_genes, vec_of_gene_names,rld_obj,orientation) {
	specified_rows <- match(vec_of_genes, row.names(rld_obj))
	specified_phenotype_mat <- assay(rld_obj)[specified_rows,]
	specified_phenotype_mat <- specified_phenotype_mat - rowMeans(specified_phenotype_mat)

	breaksList = seq(-1.0,1.0,length.out=100)
	
	if (orientation=='mito') {
		pheatmap(specified_phenotype_mat,cluster_cols=FALSE,labels_row=vec_of_gene_names,breaks=breaksList,fontsize_row=20,fontsize_column=15)
	}

	if (orientation =='Ycsome') {
		Yor_specified_phenotype_mat <- specified_phenotype_mat[,y_ord]
		pheatmap(Yor_specified_phenotype_mat,cluster_cols=FALSE,labels_row=vis_gene_names,breaks=breaksList,fontsize_row=20,fontsize_column=15)
	}

}

sva_reg_joint_count_plots <- function(regress_mat,gene_id,gene_name,orientation) {

	gene_counts <- subset(regress_mat,rownames(regress_mat)==gene_id)
	d <- cbind(t(gene_counts),mito,Ycsome,rep)
	colnames(d)[1] <- 'count'
	d <- as.data.frame(d)
	d$count <- as.numeric(as.character(d$count))
	y_title <- paste("Read Counts")

	if (orientation=='Ycsome') {
		temp <- ggplot(d,aes(x=Ycsome,y=count))+
			geom_boxplot(width=0.4)+
			scale_y_log10(breaks=seq(0,20,by=0.5))+
			theme(legend.position="none",axis.title=element_text(size=60),axis.text.x=element_text(size=55,angle=90),axis.text.y=element_text(size=45),panel.background=element_rect(fill='white',colour='grey50'),panel.grid.minor = element_line())+
			geom_point(aes(color=mito,shape=rep),position=position_jitter(w=0.05,h=0),size=8,alpha=0.75)+
			labs(x='Y Haplotype',y=y_title)+
			scale_color_manual(name='Mitochondrial Haplogroup',values=c("red4","red2","darkorange3",'darkorange','darkgoldenrod1','gold1'))+
			scale_shape(name='Replicate')
	}

	if (orientation=='mito') {
		temp <- ggplot(d,aes(x=mito,y=count))+
			geom_boxplot(width=0.4)+
			scale_y_log10(breaks=seq(0,20,by=0.5))+
			theme(legend.position="none",axis.title=element_text(size=60),axis.text.x=element_text(size=55,angle=90),axis.text.y=element_text(size=45),panel.background=element_rect(fill='white',colour='grey50'),panel.grid.minor = element_line())+
			geom_point(aes(color=Ycsome,shape=rep),position=position_jitter(w=0.05,h=0),size=8,alpha=0.75)+
			labs(x='Mitochondrial Haplotype',y=y_title)+
			scale_color_manual(name='Y Haplotype',values=c("darkgreen","green3","dodgerblue4",'deepskyblue2','purple4','mediumpurple2'))+
			scale_shape(name='Replicate') 
	}
	return(temp)

}
