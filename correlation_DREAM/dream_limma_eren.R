##4/3/19 limma on dream

#limma with eren's expression table
setwd('~/Google_Drive/ibm_synergy/dream/limma_new')
#comp = read.csv('dream_treatments_24hrs_IC20_DMSO.csv',header=TRUE)

dream_list = list('Methotrexate', 'Blebbistatin', 'Camptothecin',
               'Doxorubicin hydrochloride', 'Cycloheximide',
              'Aclacinomycin A', 'Etoposide', 'Geldanamycin', 'Mitomycin C',
              'Rapamycin', 'Trichostatin A', 'Vincristine', 'Monastrol',
              'H_7, Dihydrochloride')
log2data = read.csv('../eren_correlation/log2data.csv',row.names='X')

#stuff needed for each comparison
dmso = as.character(c(list('DM_24_'),paste0('DM_24_.',1:7)))
design = model.matrix(~-1+factor(c(1,1,1,2,2,2,2,2,2,2,2)))
colnames(design) = c('t','c')

#   ## Run limma
for(i in seq_along(dream_list)) {
  drug = as.character(paste0(substr(dream_list[[i]],1,2),'_24_I20',c('','.1','.2')))
  dfsub <- log2data[, c(drug,dmso)]
  fit <- limma::lmFit(dfsub, design) # Fit a linear model for each gene based on the given series of arrays
  contrast.matrix <- limma::makeContrasts(contrasts="t-c", levels=design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2) # Computes moderated t-statistics and log-odds of differential expression by empirical Bayes shrinkage of the standard errors towards a common value.
  limmaDF <- limma::topTable(fit2, coef=1, adjust.method='BH', sort.by="none", number=Inf)

  #write.csv(limmaDF,file=paste('eren/limma_',dream_list[[i]],'.csv',sep=''))	
	
}
 

