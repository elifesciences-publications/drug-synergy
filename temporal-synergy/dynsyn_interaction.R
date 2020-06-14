#
#!classifying treatment dependent temporal behaviour of expression
#B. Losic (2018)
#! 
#
#
library(limma)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(Mfuzz)
library(clusterSim)
library(DESeq2)
library(plyr)
library(dplyr)
#
#
setwd('~/analysis/projects/dynsyn/dynsyn')
#
#read in count data
gcounts <- read.csv('counts_bojan.csv')
rownames(gcounts) <- gcounts$X
gcounts$X <- NULL
#
phenotype <- read.csv('pheno_bojan.csv')
rownames(phenotype) <- phenotype$X
phenotype$X <- NULL
#make factor level time too
phenotype$Time_unordered <- factor(phenotype$Time)
#
identical(rownames(phenotype), colnames(gcounts))
#gene annotations
ann_g <- read.table( 'Homo_sapiens.GRCh37.75.gene.gtf', header=F)
sum(rownames(gcounts) %in% ann_g$V13)
#
#
#specify a model with time and treatment and their interaction (Treatment + Time + Treatment:Time). Each coefficeint of the regression term
#is basically from a LRT of the full interaction model with the bare model without the interaction term: Treatment + Time. 
#note that this specifies a model for increasing logFC with increasing Time, which is reasonable to first approx
#Note here that interaction simply means that we are testing difference of SLOPES between treaments across time. Thus there 
#are 6 (0,3,6,9,12,24) - 1 = 5 degrees of freedom in the design matrix. The reference is taken to be DMSO by default factor ordering.
xmm_full <- model.matrix(~0 + Treatment*Time, data = phenotype)
colnames(xmm_full) <- make.names(colnames(xmm_full))
#!edgeR operations
dge <- DGEList( counts = gcounts)
dge <- dge[rowSums(cpm(dge$counts)> .5) > 5,,keep.lib.sizes = F]
dge <- calcNormFactors(dge)
#voom transformations with just treatment distinctions
v <- voom(dge, model.matrix(~0 + Treatment, data = phenotype), plot = T)
##fit
fit_full <- eBayes(lmFit(v, xmm_full, robust = T))
summary(decideTests(fit_full))
#huge interaction effects, but most are small. Note again that interactions are automatically measured with respect to DMSO. 
#make contrasts
con <- makeContrasts(seeTM = -TreatmentT - TreatmentM + TreatmentDMSO + TreatmentTM, #so-called SYN conditions
	                 seeTW = -TreatmentT - TreatmentW + TreatmentDMSO + TreatmentTW,
	                 seeMW = -TreatmentM - TreatmentW + TreatmentDMSO + TreatmentMW,
				     transient_genesTM = TreatmentTM.Time - (TreatmentT.Time + TreatmentM.Time)/2,  #TM specific temporal evolution
					 transient_genesTW = TreatmentTW.Time - (TreatmentT.Time + TreatmentW.Time)/2,  #TW specific temporal evolution
					 transient_genesMW = TreatmentMW.Time - (TreatmentM.Time + TreatmentW.Time)/2,  #MW specific temporal evolution
					 time = Time, levels =  xmm_full)
#compute the effect sizes and FDR of these specific contrasts
fit_con <- eBayes(contrasts.fit(fit_full, con))
summary(decideTests(fit_con))
###############################
#compute empitifical Bayes moderated t-test p-values relative to a minimum effect size threshold. By definition there will be many "slopes" different
#but we want the strongest set of differences. 
fit_con_treat <- treat( fit_con, lfc = .05)
summary(decideTests(fit_con_treat))
#significantly reduced effect burden, but still appreciable, especially with TW, which has no synergistic burden. 
#TM interaction effect
topTreat(fit_con_treat, coef = 4, num = 10, p.value = .05 , sort.by = 'logFC')
#TW interaction effect
topTreat(fit_con_treat, coef = 5, num = 10, p.value = .05 , sort.by = 'AveExpr')
#TM SEE effect
topTreat(fit_con_treat, coef = 1, num = Inf, p.value = .05 )
###############################
#!
#let us collect the full results into a data frame
conmat <- colnames(con)
contrast_list <- NULL
for (i in 1:length(conmat)){
#
				#contrast_list[[i]] <- topTable( fit_con, coef = conmat[i], p.value = 1, num = Inf)
				contrast_list[[i]] <- topTreat( fit_con_treat, coef = conmat[i], p.value = 1, num = Inf)
				contrast_list[[i]]$gene <- factor(rownames(contrast_list[[i]]))
contrast_list[[i]]$sig <- factor(ifelse( -log10(contrast_list[[i]]$adj.P.Val) > 
					as.numeric(quantile(-log10(contrast_list[[i]]$adj.P.Val), probs = .995)), 'sig', 'no') )
contrast_list[[i]]$FDR_sig <- factor(ifelse( -log10(contrast_list[[i]]$adj.P.Val) < 1.3 , 'no', 'sig') )

#
							}
names(contrast_list) <- conmat
#
contrast_df <- do.call(rbind, contrast_list)
contrast_df$contrast <- factor(gsub('\\..*', '', rownames(contrast_df)))
contrast_df <- contrast_df[order(contrast_df$adj.P.Val),]
#
#unique to interaction:
#decideTests(fit_con)[, c(1,4)][decideTests(fit_con)[, c(1,4)][,1] == 0 & 
#							decideTests(fit_con)[, c(1,4)][,2] != 0,]

#volcano plot
global_volcano <- ggplot(contrast_df[contrast_df$contrast != 'time',], 
			aes(logFC, -log10(adj.P.Val))) + geom_point() +
			facet_wrap(~contrast, scales = 'free') +
			geom_text_repel(data = contrast_df[contrast_df$contrast != 'time' & 
							contrast_df$sig == 'sig',],aes(label = gene), size =2) + geom_hline(yintercept = 1.3, color = 'red')
global_volcano
#MA plot
global_MA <- ggplot(contrast_df[contrast_df$contrast != 'time',], 
			aes(AveExpr, logFC, alpha = -log10(adj.P.Val))) + geom_point() +
			facet_wrap(~contrast, scales = 'free') +
			geom_text_repel(data = contrast_df[contrast_df$contrast != 'time' & 
							contrast_df$sig == 'sig',],aes(label = gene), size =2, color = 'red', alpha = 1) 
global_MA	

pdf('global_volcano_MA_synergistic_interaction.pdf', height = 10, width = 12)
global_volcano
global_MA
dev.off()

#
#
#now let us make a fixed effect using DESeq2 with an explicit LRT test between the interacting and bare model
#specify full model. This will include ALL slope differences and not specify any ordering in time. This is then the FULL 
#subset of interactioin effects. The voom/limma model above is a SUBSET of these. However we wish to cluster all possible 
#patterns first and see in which clusters the monotonic model results fit in. 
cycl_expr_dds <- DESeqDataSetFromMatrix(countData = gcounts, colData = phenotype, 
					design = ~ Treatment + Time_unordered + Treatment:Time_unordered)
#take LRT with bare model 
cycl_expr_deseq <- DESeq(cycl_expr_dds, test = "LRT", reduced = ~ Treatment + Time_unordered)

########################################################################################################################################
########################################################################################################################################
#clumsily extract list of results
cycl_expr_table <- as.data.frame(results(cycl_expr_deseq,name = "TreatmentTM.Time_unordered24"))
cycl_expr_table$log2FoldChange <- NULL
#cycl_expr_table$lfcSE <- NULL
#cycl_expr_table$baseMean <- NULL
#take the full point wise estimates with respect to DMSO; pick only TM for now 
a <- as.data.frame(results(cycl_expr_deseq,name = "TreatmentTM.Time_unordered3")[,2:3])
colnames(a) <- c('log2FoldChange_TreatmentTM.Time_unordered3','lfcSE_TreatmentTM.Time_unordered3') 
b <- as.data.frame(results(cycl_expr_deseq,name = "TreatmentTM.Time_unordered6")[,2:3])
colnames(b) <- c('log2FoldChange_TreatmentTM.Time_unordered6','lfcSE_TreatmentTM.Time_unordered6') 
c <- as.data.frame(results(cycl_expr_deseq,name = "TreatmentTM.Time_unordered9")[,2:3])
colnames(c) <- c('log2FoldChange_TreatmentTM.Time_unordered9','lfcSE_TreatmentTM.Time_unordered9') 
d <- as.data.frame(results(cycl_expr_deseq,name = "TreatmentTM.Time_unordered12")[,2:3])
colnames(d) <- c('log2FoldChange_TreatmentTM.Time_unordered12','lfcSE_TreatmentTM.Time_unordered12')
e <- as.data.frame(results(cycl_expr_deseq,name = "TreatmentTM.Time_unordered24")[,2:3])
colnames(e) <- c('log2FoldChange_TreatmentTM.Time_unordered24','lfcSE_TreatmentTM.Time_unordered24')
#data merger
cycl_expr_table <- merge( cycl_expr_table, a, by = 0)
rownames(cycl_expr_table) <- cycl_expr_table$Row.names
cycl_expr_table$Row.names <- NULL
cycl_expr_table <- merge( cycl_expr_table, b, by = 0)
rownames(cycl_expr_table) <- cycl_expr_table$Row.names
cycl_expr_table$Row.names <- NULL
cycl_expr_table <- merge( cycl_expr_table, c, by = 0)
rownames(cycl_expr_table) <- cycl_expr_table$Row.names
cycl_expr_table <- merge( cycl_expr_table, d, by = 0)
rownames(cycl_expr_table) <- cycl_expr_table$Row.names
cycl_expr_table$Row.names <- NULL
cycl_expr_table <- merge( cycl_expr_table, e, by = 0)
rownames(cycl_expr_table) <- cycl_expr_table$Row.names
#
cycl_expr_table$Row.names <- NULL
cycl_expr_table$GeneSymbol <- factor(rownames(cycl_expr_table))
#cycl_expr_table <- merge( cycl_expr_table, gannots, by = 'GeneSymbol')
#cycl_expr_table <- cycl_expr_table[!is.na(cycl_expr_table$padj),]
cycl_expr_table$test_type <- 'treatment_time_interaction'
cycl_expr_table$significance <- factor(ifelse( cycl_expr_table$padj < .05, 'FDR < 5%', 'FDR >= 5%' ))
cycl_expr_table <- cycl_expr_table[order(cycl_expr_table$padj),]
########################################################################################################################################
########################################################################################################################################
####visualize some genes
#extract specific gene normalized counts
geneid <- 'ABCA9'
data <- plotCounts(cycl_expr_deseq, geneid, 
                   intgroup=c("Time_unordered","Treatment"), returnData=TRUE)
#reset to numeric for smoothing
data$Time_numeric <- as.numeric(as.character(data$Time_unordered))
#plot cycle
ggplot(data, aes(x = Time_numeric, y = count, color = Treatment, fill = Treatment)) + geom_point() + 
							stat_smooth(se = F, method = 'loess') + scale_y_log10() + ggtitle(geneid)


#! cluster trajectories
#data melting, merging
norm_counts <- melt(log2(counts(cycl_expr_deseq, normalized = T)+1e-100))
colnames(norm_counts) <- c("GeneSymbol", "Sample", "log2_norm_counts")
#extract log2cpm-s of genes that show significant evidence of different temporal trajectories between the disease and control groups at FDR < 5% as given by the LRT test
norm_counts_cycling <- norm_counts[which(!is.na(match(norm_counts$GeneSymbol, cycl_expr_table[cycl_expr_table$padj < .05 & 
												!is.na(cycl_expr_table$padj),]$GeneSymbol))),]
norm_counts_cycling <- merge(norm_counts_cycling, phenotype,  all.x = T, by.x = 'Sample', by.y= 0)
#generate trajectories across timepoints for each gene by calculating the median of the log2-normalized count values for each timepoint and gene in each group (Tremant vs DMSO)
median_trajectories <- dcast(norm_counts_cycling[,c(2,3,4,11)], GeneSymbol + Treatment ~ Time_unordered, median, 
																	value.var = "log2_norm_counts")
#remove zero variance trajectories and add gene annotations
library(matrixStats)
median_trajectories <- median_trajectories[rowVars(as.matrix(median_trajectories[,-c(1,2)])) != 0,]
#annot_median_trajectories <- merge(ann_g, median_trajectories,  all.y = T, by.y = 'GeneSymbol', by.x = 'V13')
#cluster disease trajectories using Mfuzz
median_TM_trajectories <- median_trajectories[median_trajectories$Treatment == "TM",]
#
clust_trajectories <- median_TM_trajectories[,-c(2)]
rownames(clust_trajectories) <- clust_trajectories$GeneSymbol
clust_trajectories$GeneSymbol <- NULL
#clust_trajectories <- clust_trajectories[,-1]

#first perform hclust to know how many clusters mfuzz needs to identify
hier_clust <- hclust(as.dist(1-cor(t(clust_trajectories), method = "spearman")), method = "complete")
pdf("disease_trajectories_hclust_dendrogam.pdf")
plot(hier_clust)
abline(h = 0.5, col = "red")
dev.off()
clusters <- cutree(hier_clust, h = 0.5)
cluster_count <- max(clusters)
#5 main clusters noted
#standardize trajectories by calculating z-scores
scaled_clust_trajectories <- data.Normalization(clust_trajectories, type = "n1", normalization = "row")
#mfuzz : cluster the trajectories
expr_set_trajectories <- ExpressionSet(assayData = as.matrix(scaled_clust_trajectories))
m1 <- mestimate(expr_set_trajectories)
#mfuzz_clust <- mfuzz(expr_set_trajectories, centers = cluster_count, m = m1)
mfuzz_clust <- mfuzz(expr_set_trajectories, centers = 5, m = m1)
#extra convenience conversions
#extract cluster memberships, assign gene names and turn cluster numbers into letter
mfuzz_clusters <- as.data.frame(mfuzz_clust$cluster)
mfuzz_clusters$Gene <- rownames(scaled_clust_trajectories)
colnames(mfuzz_clusters) <- c("Cluster", "GeneSymbol")
mfuzz_clusters$Cluster <- factor(chartr('12345', 'ABCDE', mfuzz_clusters$Cluster))
mfuzz_clusters$Cluster <- factor(chartr( 'ABCDE', 'DECBA',mfuzz_clusters$Cluster))
#merge scaled median dss trajectories, cluster memberships and gene annotations for visualizations
TM_trajectories_clusters <- merge(mfuzz_clusters, scaled_clust_trajectories, by = 0)
TM_trajectories_clusters <- TM_trajectories_clusters[, c(3,2, 4:9)]
colnames(TM_trajectories_clusters)[1] <- 'GeneSymbol'
#
TM_scaled_trajectories_long <- melt(TM_trajectories_clusters, id.vars  = c("GeneSymbol", "Cluster"), 
	 														measure.vars = c("0", "3", "6", "9", "12", "24"))
colnames(TM_scaled_trajectories_long)[c(3,4)] <- c("Time_unordered", "scaled_median_expr") 
#ok now finally plot the full expression trajectories
trajectoryplot <- ggplot(TM_scaled_trajectories_long, aes( Time_unordered, scaled_median_expr, group = GeneSymbol)) + 
							geom_line(aes(color = Cluster)) + 
							stat_smooth(aes(group = Cluster), method = "loess", se = F, color = "black") + facet_wrap(~Cluster, scales = "free") +
									ggtitle( 'Treatment specific expression trajectories') + ylab( 'scaled median expression')+
									theme(text = element_text(size = 11, face = "bold")) + xlab('Time (hours)')
trajectoryplot 
#merge to create some master tables
clusters_annots_cycling <- merge(cycl_expr_table[cycl_expr_table$padj < .05 & !is.na(cycl_expr_table$padj),], 
											 mfuzz_clusters, by = 'GeneSymbol')
#!
pdf('Fig_S3a.pdf', height = 10, width = 12)
trajectoryplot
dev.off()


#plot by gene bio-type -- mostly protein coding, not meaningful in this case because cell-type is fixed
# clusters_annots_cycling <- merge(cycl_expr_table[cycl_expr_table$padj < .05 & !is.na(cycl_expr_table$padj),], 
# 											merge( mfuzz_clusters, ann_g, by.x = 'GeneSymbol', by.y = 'V13'), by = 'GeneSymbol')
# rownames(clusters_annots_cycling) <- clusters_annots_cycling$GeneSymbol
# #				
# trajectory_cluster_plot<-	ggplot(clusters_annots_cycling, aes( V19, Cluster, color = log(baseMean), size = -log10(padj))) + geom_jitter() + 
# 			scale_colour_gradient(limits = c(min(log(clusters_annots_cycling$baseMean)), 
# 														max(log(clusters_annots_cycling$baseMean))), low = "lightgreen", high = "red") + 
# 			geom_text_repel(data = clusters_annots_cycling[clusters_annots_cycling$GeneBiotype != "protein_coding",], 
# 			aes(GeneBiotype, Cluster, label = GeneSymbol), color = "black", size = 3) + 
# 			theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), text = element_text(size = 15, face = "bold"))+
# 			ggtitle( 'Gene type by trajectory cluster ')
# trajectory_cluster_plot


#too large, but reasonable plot with some pruning: not now. 

#create one volcanoplot for all slopes of interaction given same pdj spectrum of LRT
# ints <- as.matrix( c('TreatmentTM.Time_unordered3', 'TreatmentTM.Time_unordered6', 'TreatmentTM.Time_unordered9',
# 							'TreatmentTM.Time_unordered12','TreatmentTM.Time_unordered24') )
# intcont <- NULL
# for (i in 1:nrow(ints)){
# 	 intcont[[i]] <- as.data.frame(results(cycl_expr_deseq, name = ints[i]))
# 	 intcont[[i]]$interaction <- ints[i]
# 	 					}
# intcont <- do.call(rbind, intcont)
# intcont$interaction <- factor( intcont$interaction)
# intcont$GeneSymbol <- factor(rownames(intcont))
# dat <- intcont[!is.na(intcont$padj),]
# #
# volcano_interaction_plot <- ggplot(dat, aes( log2FoldChange, -log10(padj), color = interaction) ) + 
# 											geom_point(aes(alpha = -log10(padj + 1e-10), shape = interaction)) + 
# 											geom_hline(yintercept = -log10(0.05), color = 'red' ) + 
# 			#geom_text_repel( data = dat[-log10(dat$padj) > 4 | ( abs(dat$log2FoldChange) > 8  & -log10(dat$padj) > -log10(.05)),], 
# 			#														aes(label = GeneSymbol), size = 3)+
# 			 ggtitle( 'TM specific time dependent genes') +
# 			 theme(text = element_text(size = 15, face = "bold"))
# #volcano_interaction_plot
# pdf('trajectory_volcanoplots.pdf', height = 12, width = 15)
# #multiplot(volcano_interaction_plot,trajectoryplot,trajectory_cluster_plot, cols = 3 )
# multiplot(trajectoryplot,trajectory_cluster_plot, cols = 2 )
# trajectoryplot
# trajectory_cluster_plot
# dev.off()


###
####visualize some genes
#what are the most top transient genes for each called cluster? 
ddply(clusters_annots_cycling, .(Cluster), function(x) head(x[order(x$padj[x$padj != 0]),], 3) )
#misc
head(contrast_df)
head(contrast_df[contrast_df$contrast != 'time' & (contrast_df$logFC) > 1.0 & contrast_df$adj.P.Val < .05,],20)
(contrast_df[contrast_df$contrast != 'time' & abs(contrast_df$logFC) > .2 & contrast_df$adj.P.Val < .05 & contrast_df$contrast != 'seeTM',])
#misc


#!#!
#filter the total interaction set with that of the monotonic thresholded interaction results from voom/limma model: 
summary(clusters_annots_cycling[rownames(topTreat(fit_con_treat, coef = 4, num = Inf, p.value = .05 , sort.by = 'logFC')),])
#total_filtered_interaction_TM <- clusters_annots_cycling[rownames(topTreat(fit_con_treat, coef = 4, num = Inf, p.value = .05 , sort.by = 'logFC')),]
total_filtered_interaction_TM <- merge( clusters_annots_cycling, topTreat(fit_con_treat, coef = 4, num = Inf, p.value = .05 , sort.by = 'logFC'), by = 0)
total_filtered_interaction_TM <- total_filtered_interaction_TM[order(total_filtered_interaction_TM$adj.P.Val),]
#total_filtered_interaction_TW <- clusters_annots_cycling[rownames(topTreat(fit_con_treat, coef = 5, num = Inf, p.value = .05 , sort.by = 'logFC')),]
total_filtered_interaction_TW <- merge( clusters_annots_cycling, topTreat(fit_con_treat, coef = 5, num = Inf, p.value = .05 , sort.by = 'logFC'), by = 0)

#!
TMdf <- total_filtered_interaction_TM[, c(2,21:26)]
rownames(TMdf) <- TMdf$GeneSymbol
TMdf$GeneSymbol <- NULL
write.csv( TMdf, 'Supplementary_table_2.csv', row.names = T, col.names = T, quote = F)

#!what are the best expressed TM specific genes 
tophits <- ddply(total_filtered_interaction_TM, .(Cluster), function(x) x[which.min(x$adj.P.Val),] )
tophits <- ddply(total_filtered_interaction_TM, .(Cluster), function(x) x[which.max(x$baseMean),] )
tophits <- ddply(total_filtered_interaction_TW, .(Cluster), function(x) x[which.min(x$adj.P.Val),] )
#straight topsig rank
tophits <- head(total_filtered_interaction_TM)

#tophits <- ddply(total_filtered_interaction_TW, .(Cluster), function(x) x[which.min(x$adj.P.Val),] )

#!#!
#!#!#!#plot some genes
rownames(clusters_annots_cycling) <- clusters_annots_cycling$GeneSymbol
#extract specific gene normalized counts: note subtlety: for loops have no separate variable scope so we need to use local to wrap the for block
p <- list()
for (j in 1:nrow(tophits))local({
	j <- j
geneid <- as.character(tophits$GeneSymbol[j])
clusters_annots_cycling[geneid,]
clusterid <- as.character(clusters_annots_cycling[geneid,]$Cluster)
data <- plotCounts(cycl_expr_deseq, geneid, 
                   intgroup=c("Time_unordered","Treatment"), returnData=TRUE)
#reset to numeric for smoothing
data$Time <- as.numeric(as.character(data$Time))
#plot cycle
p1 <- ggplot(data, aes(x = Time, y = count, color = Treatment, fill = Treatment)) + geom_point() + 
							stat_smooth(se = F, method = 'loess') + scale_y_log10() + ggtitle( paste(geneid, ': cluster ', clusterid, sep = '') ) + 
							scale_y_continuous(breaks = round(seq(min(data$count), max(data$count), by = floor(max(data$count)/5)),1)) #+ 
							#scale_x_discrete( labels = c("5", "12", "17", "36" ) )
print(p1)
p[[j]] <<- p1
###
								})

pdf('transientTW_genes_expression_rank2.pdf', height = 8, width = 12)
multiplot( p[[1]],p[[4]], p[[2]],p[[5]],p[[3]], cols = 3)
dev.off()

#!normalized to DMSO

p <- list()
for (j in 1:nrow(tophits))local({
	j <- j
geneid <- as.character(tophits$GeneSymbol[j])
clusters_annots_cycling[geneid,]
clusterid <- as.character(clusters_annots_cycling[geneid,]$Cluster)
data <- plotCounts(cycl_expr_deseq, geneid, 
                   intgroup=c("Time_unordered","Treatment"), returnData=TRUE)
data$replicate <- factor(gsub('.*_', '', rownames(data)))
data2 <- data[data$Treatment == 'DMSO',]
colnames(data2)[1] <- 'DMSO_count'
data <- merge( data, data2[, c('replicate', 'DMSO_count')], by = 'replicate')
data$DMSO_relative_count <- log10(data$count/data$DMSO_count)
#reset to numeric for smoothing
data$Time <- as.numeric(as.character(data$Time))
#plot cycle
p1 <- ggplot(data[data$Treatment %in% c('T', 'M', 'TM'),], aes(x = Time, y = DMSO_relative_count, color = Treatment, fill = Treatment)) + geom_point() + 
							stat_smooth(se = F, method = 'loess') + ggtitle( paste(geneid, ': cluster ', clusterid, sep = '') ) + 
							scale_y_continuous(breaks = round(seq(min(data$DMSO_relative_count), max(data$DMSO_relative_count), 
																								by = (max(data$DMSO_relative_count)/5)),1)) + 
							xlab( 'Time (hours)') + ylab( 'log10(cpm) - log10(cpm_DMSO)' )																	 #+ 
							#scale_x_discrete( labels = c("5", "12", "17", "36" ) )
print(p1)
p[[j]] <<- p1
###
								})

pdf('Fig_S3b.pdf', height = 8, width = 12)
multiplot( p[[1]],p[[4]], p[[2]],p[[5]],p[[3]], cols = 3)
dev.off()


##################################################################################################################
##################################################################################################################
#!clustering via mfuzz for just 925 TM specific genes
# tmig <- topTreat(fit_con_treat, coef = 4, num = Inf, p.value = .05 , sort.by = 'logFC') 
# scaled_clust_trajectories_tmig <- scaled_clust_trajectories[rownames(scaled_clust_trajectories) %in% rownames(tmig),]
# ##
# expr_set_trajectories_tmig <- ExpressionSet(assayData = as.matrix(scaled_clust_trajectories_tmig))
# m1_tmig <- mestimate(expr_set_trajectories_tmig)
# #mfuzz_clust <- mfuzz(expr_set_trajectories, centers = cluster_count, m = m1)
# mfuzz_clust_tmig <- mfuzz(expr_set_trajectories_tmig, centers = 4, m = m1_tmig)
# #extra convenience conversions
# #extract cluster memberships, assign gene names and turn cluster numbers into letter
# mfuzz_clusters_tmig <- as.data.frame(mfuzz_clust_tmig$cluster)
# mfuzz_clusters_tmig$Gene <- rownames(scaled_clust_trajectories_tmig)
# colnames(mfuzz_clusters_tmig) <- c("Cluster", "GeneSymbol")
# mfuzz_clusters_tmig$Cluster <- factor(chartr('12345', 'ABCDE', mfuzz_clusters_tmig$Cluster))
# mfuzz_clusters_tmig$Cluster <- factor(chartr( 'ABCDE', 'DECBA',mfuzz_clusters_tmig$Cluster))
# #merge scaled median dss trajectories, cluster memberships and gene annotations for visualizations
# TM_trajectories_clusters_tmig <- merge(mfuzz_clusters_tmig, scaled_clust_trajectories_tmig, by = 0)
# TM_trajectories_clusters_tmig <- TM_trajectories_clusters_tmig[, c(3,2, 4:9)]
# colnames(TM_trajectories_clusters_tmig)[1] <- 'GeneSymbol'
# #
# TM_scaled_trajectories_tmig_long <- melt(TM_trajectories_clusters_tmig, id.vars  = c("GeneSymbol", "Cluster"), 
# 	 														measure.vars = c("0", "3", "6", "9", "12", "24"))
# colnames(TM_scaled_trajectories_tmig_long)[c(3,4)] <- c("Time_unordered", "scaled_median_expr") 
# #ok now finally plot the full expression trajectories
# trajectoryplot_tmig <- ggplot(TM_scaled_trajectories_tmig_long, aes( Time_unordered, scaled_median_expr, group = GeneSymbol)) + 
# 							geom_line(aes(color = Cluster)) + 
# 							stat_smooth(aes(group = Cluster), method = "loess", se = F, color = "black") + facet_wrap(~Cluster, scales = "free") +
# 									ggtitle( 'TM specific expression trajectories') + ylab( 'scaled median expression')+
# 									theme(text = element_text(size = 15, face = "bold"))
# trajectoryplot_tmig
##################################################################################################################
##################################################################################################################
#!clustering via mfuzz for just 925 TM specific genes


# geneid <- 'FASN'
# clusters_annots_cycling[geneid,]
# clusterid <- clusters_annots_cycling[geneid,]$Cluster
# data <- plotCounts(cycl_expr_deseq, geneid, 
#                    intgroup=c("Time_unordered","Treatment"), returnData=TRUE)
# #reset to numeric for smoothing
# data$Time <- as.numeric(as.character(data$Time))
# #plot cycle
# clustgene <- ggplot(data, aes(x = Time, y = count, color = Treatment, fill = Treatment)) + geom_point() + 
# 							stat_smooth(se = F, method = 'loess') + scale_y_log10() + ggtitle( paste(geneid, ': cluster ', clusterid, sep = '') ) + 
# 							scale_y_continuous(breaks = round(seq(min(data$count), max(data$count), by = floor(max(data$count)/5)),1)) #+ 
# clustgene							#scale_x_
# summary(clusters_annots_cycling[rownames(topTreat(fit_con_treat, coef = 4, num = Inf, p.value = .05 , sort.by = 'logFC')),])


# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
