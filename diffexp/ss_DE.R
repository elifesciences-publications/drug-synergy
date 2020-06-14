ss_count <- read.delim("C:/Users/mehmeteren/OneDrive/Drug_synergy_R_studio/self_synergy/self_synergy_fcounts_ucsc.mtrx", comment.char="#", stringsAsFactors=FALSE)
rownames(ss_count)=ss_count[,1]
ndim=dim(ss_count)[2]
ss_count=ss_count[,7:ndim]
rm(ndim)
y_ss=DGEList(ss_count)
#is_expr=rowSums(cpm(y_ss) > 0.5) >=3


y_ss <- y_ss[rownames(y_MCF7),,keep.lib.sizes=FALSE]

y_ss <- calcNormFactors(y_ss)
#rm(is_expr)


ids_ss=substr(colnames(ss_count),61,65)
colnames(y_ss)=ids_ss
colnames(ss_count)=colnames(y_ss)

y_ss=y_ss[,rownames(ss_phenotype)]
ss_count=ss_count[,rownames(ss_phenotype)]

ss_phenotype <- read.delim("C:/Users/mehmeteren/OneDrive/Drug_synergy_R_studio/self_synergy/pheno_ss.txt", stringsAsFactors=FALSE)
ss_viability<-read.delim("C:/Users/mehmeteren/OneDrive/Drug_synergy_R_studio/self_synergy/dose_ss.txt", stringsAsFactors=FALSE)

rownames(ss_phenotype)=ss_phenotype[,1]

nphen=length(colnames(ss_phenotype))
ss_phenotype=cbind(ss_phenotype,rep(c("A","B"),9))
colnames(ss_phenotype)[nphen+1]="Replicate"


ss_phenotype=cbind(ss_phenotype,paste0(ss_phenotype$Drug,"_",
                                       ss_phenotype$Concentration..uM.,"_",ss_phenotype$Replicate))
colnames(ss_phenotype)[nphen+2]="Sample"
ss_phenotype=cbind(ss_phenotype,paste0(ss_phenotype$Drug,"_",
                                       ss_phenotype$Concentration..uM.))
colnames(ss_phenotype)[nphen+3]="Contrast"

colnames(y_ss)=ss_phenotype[colnames(y_ss),"Sample"]




ss_design=model.matrix( ~0 + Contrast,data = ss_phenotype)
rownames(ss_design)=colnames(y_ss)
colnames(ss_design)=substr(colnames(ss_design),9,nchar(colnames(ss_design)))


ss_voom <- voom(y_ss,ss_design, plot = TRUE)
ss_fit=lmFit(ss_voom,ss_design)
## Combine and Plot
##
###

conds1=c(paste0("DMSO_24",c("A","B","C")),
         paste0("T_24",c("A","B","C")),
         paste0("M_24",c("A","B","C")))

conds=colnames(counts_all_n48)
y_dummy=cbind(counts_all_n48[,conds],ss_count)
y_dummy=DGEList(y_dummy)
is_expr=rowSums(cpm(y_dummy) > 0.5) >=3


y_dummy <- y_dummy[is_expr,,keep.lib.sizes=FALSE]

y_dummy <- calcNormFactors(y_dummy)
#colnames(y_dummy)=c(colnames(y_MCF7),colnames(y_ss))




samples=setdiff(ss_phenotype$Contrast,"DMSO")

for( i in samples)
{
  
  cond=i
  control="DMSO"
  contrast.matrix1 <- makeContrasts(contrasts=paste0(cond,"-",control),levels=ss_design)
  fit2 <- contrasts.fit(ss_fit, contrast.matrix1)
  fit2 <-eBayes(fit2)
  dummy=topTable(fit2,  adjust="BH",p.value=1,n=Inf)
  dummyname=paste0("ss_limma_eb_",cond)
  assign(dummyname,dummy) 
  rm(cond,control,dummy,dummyname,contrast.matrix1,fit2)
}

