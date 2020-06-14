library(edgeR)
library(limma)




lncap_counts=read.table("data/lncap_ucsc_raw_counts",header = TRUE)
lncap_annotation=read.csv("data/lncap_anno.csv",row.names = 1)

lncap_count=lncap_counts[,7:150]

rownames(lncap_count)=as.character(lncap_counts[,1])

y_lncap=DGEList(lncap_count)
ids_lncap=substr(colnames(lncap_counts)[7:150],59,63)
colnames(y_lncap)=ids_lncap

is_expr=rowSums(cpm(y_lncap) > 0.5) >=3


y_lncap <- y_lncap[is_expr,,keep.lib.sizes=FALSE]

y_lncap <- calcNormFactors(y_lncap)


lncap_annotation=lncap_annotation[colnames(y_lncap),]

lncap_design=model.matrix( ~0 + Contrast,data = lncap_annotation)
colnames(lncap_design)=substr(colnames(lncap_design),9,nchar(colnames(lncap_design)))
rownames(lncap_design)=colnames(y_lncap)

lncap_voom <- voom(y_lncap,lncap_design, plot = TRUE)
lncap_fit=lmFit(lncap_voom,lncap_design)


conds=c("TM","T","M","W","TW","MW","TMW")
time=c(0,3,6,9,12,24)
time=as.character(time)

for( i in conds)
{
  for(j in time)
  {
    cond=paste0(i,"_",j)
    control=paste0("DMSO_",j)
    contrast.matrix1 <- makeContrasts(contrasts=paste0(cond,"-",control),levels=lncap_design)
    fit2 <- contrasts.fit(lncap_fit, contrast.matrix1)
    fit2 <-eBayes(fit2)
    dummy=topTable(fit2,  adjust="BH",p.value=1,n=Inf)
    dummyname=paste0("lncap_limma_eb_",cond)
    assign(dummyname,dummy) 
    rm(cond)
  }
  rm(control)
}


rm(time,conds)