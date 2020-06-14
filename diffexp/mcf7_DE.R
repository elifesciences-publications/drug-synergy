load("data/counts_all.Rdata")
pheno_mcf7=readRDS("data/pheno_mcf7.RDS")

counts_all_n48=counts_all[,which(pheno_mcf7[,2]!="48")]
pheno_mcf7_n48=pheno_mcf7[which(pheno_mcf7[,2]!="48"),]
pheno_mcf7_n48=data.frame(pheno_mcf7_n48) 
design_n48=model.matrix(~ 0 + Contrast, data=pheno_mcf7_n48)
colnames(design_n48)=substr(colnames(design_n48),9,nchar(colnames(design_n48)))
y_MCF7<-DGEList(counts=counts_all_n48)
isexpr <- rowSums(counts_all_n48[,which(pheno_mcf7_n48[,1]=="DMSO")] > 1) >=3

y_MCF7<-y_MCF7[isexpr,,keep.lib.sizes=FALSE]

y_MCF7 <- calcNormFactors(y_MCF7)
voom_data <- voom(y_MCF7, design_n48, plot = TRUE)

fit=lmFit(voom_data,design_n48)






############ limma
samples_for_diff=setdiff(colnames(design_n48),
                         paste0("DMSO_",c(0,3,6,9,12,24)))

for( i in samples_for_diff)
{  
  control=paste0("DMSO_",substr(i,gregexpr("_",i)[[1]][1]+1,nchar(i)))
  contrast.matrix1 <- makeContrasts(contrasts=paste0(i,"-",control),levels=design_n48)
  fit2 <- contrasts.fit(fit, contrast.matrix1)
  fit2 <-eBayes(fit2)
  dummy=topTable(fit2,  adjust="BH",p.value=1,n=Inf)
  dummyname=paste0("limma_eb_",i)
  assign(dummyname,dummy)
  l=l+1
  rm(fit2,control,contrast.matrix1,dummyname,dummy)
}

