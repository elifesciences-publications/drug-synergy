#check kmeans cluster lists against all gene sets of selected libraries
#nov 4 2016

import pandas as pd
from itertools import product
from scipy import stats
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mult




def fishers(row):
    de = sets[row.dxset]
    s = set(cand.ix[row.candset].dropna()).intersection(allRNAseq)
    notde = allRNAseq - de
    nots = allRNAseq - s
    ##contingency table: 
            ## [[diffex&inset,diffex&notinset,
            ## [notdiffex&inset,notdiffex&notinset]]
    if len(de)>=minlist and len(s)>=minlist and len(de.intersection(s))>=minov:
        return stats.fisher_exact([[len(de.intersection(s)),len(de.intersection(nots))], \
                [len(notde.intersection(s)),len(notde.intersection(nots))]])[1]
    else:
        return None


##set and dysregulated genes list minimum length
minlist = 3
##set minimum overlap 
minov = 2



finalsets = [
            #cell cycle
            'mitotic cell cycle (GO:0000278)', \
            
            #ER
            'BHAT_ESR1_TARGETS_NOT_VIA_AKT1_UP', 'BHAT_ESR1_TARGETS_NOT_VIA_AKT1_DN', \

            #growth
             'negative regulation of growth (GO:0045926)',

            #nuclease
            'positive regulation of nuclease activity (GO:0032075)', \
            'regulation of nuclease activity (GO:0032069)', \

            #unfolded protein
            'cellular response to external stimulus (GO:0071496)', \
            
             'activation of signaling protein activity involved in unfolded protein response (GO:0006987)', \
            'response to unfolded protein (GO:0006986)', \
            'response to topologically incorrect protein (GO:0035966)', \
            'cellular response to topologically incorrect protein (GO:0035967)', \
             'cellular response to unfolded protein (GO:0034620)', \
            'endoplasmic reticulum unfolded protein response (GO:0030968)', \
            'ER-nucleus signaling pathway (GO:0006984)', \
            'response to endoplasmic reticulum stress (GO:0034976)', \

            #cell death
            'intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress (GO:0070059)', \
            'intrinsic apoptotic signaling pathway (GO:0097193)', \
            'apoptotic signaling pathway (GO:0097190)',  \
            
             #immune
            'regulation of cytokine production (GO:0001817)', \
            'toll-like receptor signaling pathway (GO:0002224)', \
            'MyD88-independent toll-like receptor signaling pathway (GO:0002756)', \
             'TRIF-dependent toll-like receptor signaling pathway (GO:0035666)', \
            'toll-like receptor 3 signaling pathway (GO:0034138)', \
            'innate immune response-activating signal transduction (GO:0002758)', \
            'activation of innate immune response (GO:0002218)', \

            #signaling
            'pattern recognition receptor signaling pathway (GO:0002221)', \
            'transforming growth factor beta receptor signaling pathway (GO:0007179)', \

            #phosphorylation
             'negative regulation of phosphorylation (GO:0042326)', \
            'positive regulation of kinase activity (GO:0033674)', \

            #transcription
            'regulation of sequence-specific DNA binding transcription factor activity (GO:0051090)', \
            'transcription from RNA polymerase II promoter (GO:0006366)', \
            
            #autophagy
            'autophagy (GO:0006914)', \
            'positive regulation of autophagy (GO:0010508)', \
            
            #metabolism
            'cofactor metabolic process (GO:0051186)', \
            'generation of precursor metabolites and energy (GO:0006091)', \
             'carbohydrate derivative biosynthetic process (GO:1901137)', \
            'organophosphate biosynthetic process (GO:0090407)', \
            'sulfur compound metabolic process (GO:0006790)', \
            'cellular aldehyde metabolic process (GO:0006081)', \
            'cellular carbohydrate metabolic process (GO:0044262)', \
            'cellular modified amino acid metabolic process (GO:0006575)', \
            
            ]




#collect clusters
sets = {}
clusters = [1,2,3,4]
df = pd.read_csv('T-M-TM_k_4_bestKc.tsv',header=0,sep='\t')
for cluster in clusters:
    sets[cluster] = set(df[df.cluster ==cluster].gene)
    

#background set
allRNAseq = set(pd.read_csv('../gene_expression/diff_limma_december/limma_t0_M_0.txt',header=0,
sep='\t',index_col='1').index)



#collect gene set libraries of interest
go = pd.read_csv('../libraries_from_enrichr/GO_Biological_Process.txt',sep='\t',
                 names = list(range(10000)))

##kegg = pd.read_csv('../libraries_from_enrichr/KEGG_2015.txt',sep='\t',
##                   names = list(range(300)))

er = pd.read_csv('genesets.gmt',sep='\t',names = list(range(300)))

cand = pd.concat([go,er]).set_index(0)
cand = cand[cand.index.isin(finalsets)]



rslt = pd.DataFrame(list(product(sets.keys(),cand.index)), columns=['dxset', 'candset'])
rslt['pval'] = rslt.apply(fishers,axis=1)
nan = rslt[np.isnan(rslt.pval)]
nan['fdr'] = np.nan
notnan = rslt[~np.isnan(rslt.pval)]
notnan['fdr'] = mult(notnan.pval,method='fdr_bh')[1]
rslt = pd.concat([nan,notnan])
rslt.reset_index(drop=True,inplace=True)

fdrs = rslt.pivot('candset','dxset','fdr')
print(fdrs.shape)
fdrs.dropna(how='all',inplace=True)
print(fdrs.shape)
fdrs.fillna(1.,inplace=True)
fdrs = fdrs[fdrs[fdrs<0.05].any(axis=1)]
print(fdrs.shape)
fdrs.to_csv('new_candidate_geneset_4clusters_fdr.tsv',sep='\t')
logged = -np.log10(fdrs)
#thomas sent me a renumbering of the clusters
logged = logged.T[finalsets].T[[1,4,2,3]]
logged.to_csv('new_candidate_geneset_4clusters_fdr_log_.tsv',sep='\t')

