#B. Losic 
#Given exon level counts (with unique rownames), return differentially expressed exons along with estimates of genes that are differentially spliced as a whole
#adjust model xmm, and count thresholds prior to running
#USAGE
#module load R/3.2.0
#Rscript DEX_DSG_dynsyn.R exon_raw_counts phenotype exon_hg19_annotation


args <- commandArgs(TRUE)
library(limma)
library(edgeR)
ecounts <- read.table( args[1], row.names = 1 , header=T)
phenotype <- read.table( args[2], row.names =1 , header=T)
annotation <- read.table( args[3], header=T)

#21 nondiagonal combinations amongst 6 members

paircombos <- as.matrix(t(combn( unique(phenotype$combo), 2)))

for (i in 1:dim(paircombos)[1]){

#pairs = c("TMX_MFL", "MFL_singlet")
pairs <- as.character(paircombos[i,])
subset_idx <- rownames((phenotype[ phenotype$combo == pairs[1] | phenotype$combo == pairs[2], ]))
nphenotype <-  phenotype[subset_idx, ]
nphenotype$combo <- factor( nphenotype$combo, pairs)

#Construct DGEList to conveniently encode annotation and library normalization
dge_exon <- DGEList( counts = ecounts[, rownames(nphenotype) ], genes = annotation )
#mild filter
dge_exon <- dge_exon[ rowSums(cpm(dge_exon) > 1) > 9,,keep.lib.sizes=FALSE]
dge_exon <- calcNormFactors(dge_exon)

#simple testing the combo contrasts 

xmm <- model.matrix( ~as.factor(combo), data = nphenotype)

#voom
v_exon <- voom( dge_exon, xmm )

##Differentially expressed exons

#Threshold lfc 0

fx <- lmFit(v_exon, xmm)
fit_DE_exons <- eBayes( fx, robust = T)

#Threshold lfc 0.5

fit_DE_exons_treat <- treat( fx, lfc = 0.5)

#Threshold lfc 1.0

fit_DE_exons_treat2 <- treat( fx, lfc = 1.0)

#Write DE exon results relative to coefficient of interest

write.table( topTable( fit_DE_exons, coef = 2, num = Inf, p.value = 0.05), paste( pairs[1], pairs[2], "DE_lmfit", args[1], args[3], sep = "_"), row.names = T, col.names = T, sep = "\t", quote = F)
write.table( topTreat( fit_DE_exons_treat, coef = 2, num = Inf, p.value = 0.05), paste( pairs[1], pairs[2], "DE_treat_lfc_0.5", args[1], args[3],  sep = "_"), row.names = T, col.names = T, sep = "\t", quote = F)
write.table( topTreat( fit_DE_exons_treat2, coef = 2, num = Inf, p.value = 0.05), paste( pairs[1], pairs[2], "DE_treat_lfc_1.0", args[1], args[3],  sep = "_"), row.names = T, col.names = T, sep = "\t", quote = F)

##############Whole-gene statistics for differential exon usage 

#geneids <- gsub( ":.*", "", rownames(dge_exon$counts))
#exonids <- rownames(dge_exon$counts)
#ex <- diffSplice( fx, geneid = geneids, exonid = exonids)

ex <- diffSplice( fx, geneid="Geneid", exonid = "Start" )

#Write splicing results

#t-test for each exon's differences between all other exons of a given gene
write.table( topSplice(ex, coef=2, test="t", number=Inf, FDR=.1), paste( pairs[1], pairs[2], args[1], "exons_intragene_t_test_FDR_10", sep = "_"), row.names = T, col.names = T, sep = "\t", quote = F)
#F-test at gene level for DS given several exons having signal
write.table( topSplice(ex, coef=2, test="F", number=Inf, FDR=.1), paste( pairs[1], pairs[2], args[1], "genes_F_test_splicing_FDR_10", sep = "_"), row.names = T, col.names = T, sep = "\t", quote = F)
#simes-test at gene level for DS given rare exons having signal
write.table( topSplice(ex, coef=2, test="simes", number=Inf, FDR=.1), paste( pairs[1], pairs[2], args[1], "genes_simes_test_rare_splicing_FDR_10", sep = "_"), row.names = T, col.names = T, sep = "\t", quote = F)


##############################################################

volcanoplot2 <- function (fit, coef = 1, highlight = 0, names = fit$genes$ID,
    xlab = "Log Fold Change", ylab = "Log Odds", pch = 16, cex = 0.35,
    ...)
{
    if (!is(fit, "MArrayLM"))
        stop("fit must be an MArrayLM")
    if (is.null(fit$lods))  
        stop("No B-statistics found, perhaps eBayes() not yet run")
    x <- as.matrix(fit$coef)[, coef]
    y <- as.matrix(fit$lods)[, coef]
    plot(x, y, xlab = xlab, ylab = ylab, pch = pch, cex = cex,
        ...)
    if (highlight > 0) {
        if (is.null(names))
            names <- 1:length(x)
        names <- as.character(names)
        o <- order(x, decreasing = TRUE)
        on <- order(-x, decreasing=TRUE)
        i <- o[1:highlight]
        text(x[i], y[i], labels = names[i],
            cex = 0.5, col = "blue", pos=3)
        i2 <- on[1:round(highlight)]
        text(x[i2], y[i2], labels = names[i2],
            cex = 0.5, col = "red", pos = 3)
    }
    invisible()
}

######



#plots

pdf(paste(pairs[1], pairs[2], args[1], "exon_summary.pdf", sep = "_"))
#pdf('exon_summary.pdf')
v_exon <- voom( dge_exon, xmm, plot = T)
plotMDS(v_exon, col = as.numeric(nphenotype$combo), labels = colnames( v_exon), cex = 0.6, main = paste( pairs[1], pairs[2], "v$E", sep = "_"),cex.main=0.7,cex.lab=0.8 )
volcanoplot2(fit_DE_exons, coef = 2, highlight = 50, names = rownames(fit_DE_exons) )
top_DS_genes <- as.matrix(topSplice( ex, coef = 2, num = 50, FDR = .05, test = "F")$Geneid)
for (i in 1:dim(top_DS_genes)[1]){ 
plotSplice( ex, coef = 2, geneid = top_DS_genes[i], genecol = "Geneid", FDR = 0.05) 
plotExons( fit_DE_exons, genecolname = "Geneid", exoncolname = "Start", coef = 2, geneid = top_DS_genes[i], FDR = 0.05) 
				 }
#top_DE_genes <- as.matrix( topTable( fit_DE_exons, coef = 2, num = 20, p.value = .05)$Geneid)
#for (i in 1:dim(top_DE_genes)[1]){ plotExons( fit_DE_exons, genecolname = "Geneid", exoncolname = "Start", coef = 2, geneid = top_DE_genes[i], FDR = 0.05) } 				
dev.off()

#plotExons(fit_DE_exons, coef = 2, genecolname = "Geneid", exoncolname="Start", geneid = "Igfbp5")
#plotSplice( ex, coef = 2, geneid = "Igfbp5", genecol = "Geneid")


save.image(paste( pairs[1], pairs[2], args[1], 'splice_voom.RData', sep = "_"))

}
