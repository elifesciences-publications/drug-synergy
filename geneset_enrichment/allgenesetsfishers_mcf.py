#check diff exp lists against all gene sets of selected libraries
#oct 15 2016

import pandas as pd
from itertools import product
from scipy import stats
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mult

##finallogged = pd.read_csv('hit_geneset_enrichment_fdr_log.tsv',sep='\t',index_col = 'candset')#
##print('read final file from previous run')

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


##set p value cutoff for diff exp
pval = 1E-18
    
txs = ['TM','M','T','MW','TW','W'] #
times = ['3','6','9','12','24']
dirs = ['up','down']

denames = []
for tx in txs:
    for time in times:
        denames.append('../gene_expression/diff_limma_december/limma_t0_'+tx+'_'+time+'.txt')

dfs = {}
sets = {}
for name in denames:
    df = pd.read_csv(name,header=0,sep='\t',index_col='1')
    dfs['up_'+name.split('t0')[1].split('_')[1]+'_'+
        name.split('t0')[1].split('_')[2].strip('.txt')] = df[(df['adj.P.Val'] <pval)&(df['logFC']>0)].index
    dfs['down_'+name.split('t0')[1].split('_')[1]+'_'+
        name.split('t0')[1].split('_')[2].strip('.txt')] = df[(df['adj.P.Val'] <pval)&(df['logFC']<0)].index

    sets['up_'+name.split('t0')[1].split('_')[1]+'_'+
        name.split('t0')[1].split('_')[2].strip('.txt')] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']>0)].index)
    sets['down_'+name.split('t0')[1].split('_')[1]+'_'+
        name.split('t0')[1].split('_')[2].strip('.txt')] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']<0)].index)


for time in times:
    for di in dirs:
        sets[di+'_TUM_'+time] = sets[di+'_T_'+time].union(sets[di+'_M_'+time])
        sets[di+'_TUW_'+time] = sets[di+'_T_'+time].union(sets[di+'_W_'+time])
        sets[di+'_MUW_'+time] = sets[di+'_M_'+time].union(sets[di+'_W_'+time])


allRNAseq = set(pd.read_csv('../gene_expression/diff_limma_december/limma_t0_M_0.txt',header=0,
sep='\t',index_col='1').index)

##set and dysregulated genes list minimum length
minlist = 3
##set minimum overlap 
minov = 2


#collect gene set libraries of interest
go = pd.read_csv('../libraries_from_enrichr/GO_Biological_Process.txt',sep='\t',
                 names = list(range(10000)))

##kegg = pd.read_csv('../libraries_from_enrichr/KEGG_2015.txt',sep='\t',
##                   names = list(range(300)))

er = pd.read_csv('genesets.gmt',sep='\t',names = list(range(300)))

cand = pd.concat([go,er]).set_index(0)



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
fdrs.to_csv('all_geneset_enrichment_fdr.tsv',sep='\t')
logged = -np.log10(fdrs)
logged.to_csv('all_geneset_enrichment_fdr_log_.tsv',sep='\t')


upcols= []
downcols = []
uupcols = []
udowncols = []
tmcols = []
for col in fdrs.columns:
    if any(st in col for st in ['_TM_','_T_','_M_','_TUM_']):
        tmcols.append(col)
    if '_TM_' in col:
        if 'up_' in col:
            upcols.append(col)
        else:
            downcols.append(col)
    elif '_TUM_'in col:
        if 'up_' in col:
            uupcols.append(col)
        else:
            udowncols.append(col)
tmup = fdrs[upcols]
tmdown = fdrs[downcols]
uup = fdrs[uupcols]
udown = fdrs[udowncols]
tm = fdrs[tmcols]

#tmup[tmup[tmup<0.00001].any(axis=1)]
#tmup[(tmup[tmup<0.00001].any(axis=1))&((tmup<0.01).all(axis=1))] --> 23 processes
#tmdown[(tmdown[tmdown<0.00001].any(axis=1))&((tmdown<0.01).all(axis=1))]

##unexpected gene sets
sig = fdrs[(((tmup[tmup<0.00001].any(axis=1))&((tmup<0.01).all(axis=1)))&
     ~((uup[uup<0.00001].any(axis=1))&((uup<0.01).all(axis=1))))| #synergistic
           (((tmdown[tmdown<0.00001].any(axis=1))&((tmdown<0.01).all(axis=1)))&
     ~((udown[udown<0.00001].any(axis=1))&((udown<0.01).all(axis=1))))| #synergistic
     ((uup[uup<0.00001].any(axis=1))&((uup<0.01).all(axis=1)))|#additive
           ((udown[udown<0.00001].any(axis=1))&((udown<0.01).all(axis=1)))] #additive

##drug targets
targets = 'autophagy|BHAT_ESR1_TARGETS|estrogen'
tar = fdrs[fdrs.index.str.contains(targets)&((tmup[tmup<0.01].any(axis=1))|(tmdown[tmdown<0.01].any(axis=1)))]

##hallmarks of cancer
#selected from each high level GO term associated with a hallmark of cancer (weinberg update), and 1-2 levels of its children connected by "is a" 
candidatesets = [
'GO:0001525', #angiogenesis
'GO:0001837', #emt
'GO:0006955', #immune response
'GO:0045087|GO:0006959|GO:0002250|GO:0002418', #innate humoral adaptive | response to tumor
'GO:0008219|GO:0012501|GO:0070265|GO:0036473', #cell death | programmed cell death| necrotic cell death | cell death in response to oxidative stress  
'GO:0097190|GO:0097191|GO:0097193|GO:0070059|GO:0006915|GO:0043276', #apoptotic signaling pathway | extrinsic | intrinsic | response to ER stress | apoptotic process | anoikis
'GO:0001896|GO:0070269|GO:0097468|GO:0097707|GO:0048102|GO:0097300', #autolysis|pyroptosis|response to ROS|ferroptosis|cornification|autophagy|programmed necrosis
'GO:0007569|GO:0001302|GO:0090399|GO:0090398|GO:0001300', #cell aging | replicative cell aging | replicative senscence | cellular senescence | chronological cell aging
'GO:0008283|GO:0050673', #cell proliferation | epithelial cell proliferation
'GO:0007049|GO:0000278', #cell cycle | mitotic cell cycle
'GO:0040007|GO:0045927|GO:0045926', #growth positive negative
'GO:0006954|GO:0002544|GO:0002437|GO:0090594|GO:0002526', ##immune reponse | chronic | antigen| wounding | acute
'GO:0008152|GO:0044237',  #metabolism | cellular | children
'GO:0034641|GO:0006575|GO:0009850|GO:0072592|GO:0017144|GO:0006081|GO:0090345|GO:0006793|GO:0044107|GO:0044260|\
GO:0097354|GO:0043446|GO:0044248|GO:0035383|GO:0044249|GO:0046950|GO:0009404|GO:0046483|GO:0006805|GO:0018942|\
GO:0043449|GO:0006091|GO:0010191|GO:0072593|GO:0051186|GO:0031323|GO:0008218|GO:0042133|GO:0044255|GO:0006790|\
GO:0006730|GO:0006725|GO:0006413|GO:0006082|GO:0042810|GO:0006127|GO:0042180|GO:2001290|GO:0001887|GO:0034754|\
GO:0044262|GO:0043452']

#GO genesets
sets = ''
for i in candidatesets:
    sets +=i
    sets +='|'
sets = sets.strip('|')

###kegg genesets
##k = 'angiogenesis|EMT|epithelial to mesenchymal transition|immune|inflammation|cell death|apoptosis|autophagy|proliferation|senescence|growth|metabolism|invasion|metastasis'

c = fdrs[(fdrs.index.str.contains(sets))&((tmup[tmup<0.01].any(axis=1))|(tmdown[tmdown<0.01].any(axis=1)))]

final = pd.concat([sig,tar,c]).drop_duplicates()

final.to_csv('hit_geneset_enrichment_fdr.tsv',sep='\t')
finallogged = -np.log10(final)
finallogged.to_csv('hit_geneset_enrichment_fdr_log.tsv',sep='\t')


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


finallogged = finallogged.T[finalsets].T

tmcols = ['down_M_3', 'down_M_6', 'down_M_9', 'down_M_12', 'down_M_24',
        'down_T_3', 'down_T_6', 'down_T_9', 'down_T_12', 'down_T_24',
        'down_TUM_3', 'down_TUM_6','down_TUM_9', 'down_TUM_12', 'down_TUM_24',
        'down_TM_3', 'down_TM_6', 'down_TM_9','down_TM_12', 'down_TM_24',
        
        'up_M_3', 'up_M_6', 'up_M_9', 'up_M_12', 'up_M_24',
        'up_T_3', 'up_T_6', 'up_T_9', 'up_T_12', 'up_T_24',
        'up_TUM_3', 'up_TUM_6', 'up_TUM_9', 'up_TUM_12', 'up_TUM_24', 
        'up_TM_3', 'up_TM_6', 'up_TM_9','up_TM_12', 'up_TM_24']

finallogged[tmcols].to_csv('TM.tsv',sep='\t')

mwcols = ['down_M_3', 'down_M_6', 'down_M_9', 'down_M_12', 'down_M_24',
          'down_W_3', 'down_W_6', 'down_W_9', 'down_W_12', 'down_W_24',
        'down_MUW_3', 'down_MUW_6','down_MUW_9', 'down_MUW_12', 'down_MUW_24',
        'down_MW_3', 'down_MW_6', 'down_MW_9','down_MW_12', 'down_MW_24',
        
        'up_M_3', 'up_M_6', 'up_M_9', 'up_M_12', 'up_M_24',
          'up_W_3', 'up_W_6', 'up_W_9', 'up_W_12', 'up_W_24',
        'up_MUW_3', 'up_MUW_6', 'up_MUW_9', 'up_MUW_12', 'up_MUW_24', 
        'up_MW_3', 'up_MW_6', 'up_MW_9','up_MW_12', 'up_MW_24']

finallogged[mwcols].to_csv('MW.tsv',sep='\t')

twcols = [
        'down_T_3', 'down_T_6', 'down_T_9', 'down_T_12', 'down_T_24',
        'down_W_3', 'down_W_6', 'down_W_9', 'down_W_12', 'down_W_24',
        'down_TUW_3', 'down_TUW_6','down_TUW_9', 'down_TUW_12', 'down_TUW_24',
        'down_TW_3', 'down_TW_6', 'down_TW_9','down_TW_12', 'down_TW_24',
        
        'up_T_3', 'up_T_6', 'up_T_9', 'up_T_12', 'up_T_24',
        'up_W_3', 'up_W_6', 'up_W_9', 'up_W_12', 'up_W_24',
        'up_TUW_3', 'up_TUW_6', 'up_TUW_9', 'up_TUW_12', 'up_TUW_24', 
        'up_TW_3', 'up_TW_6', 'up_TW_9','up_TW_12', 'up_TW_24']

finallogged[twcols].to_csv('TW.tsv',sep='\t')

#excluded: things not relevant to breast cancer cells i.e. photosynthesis, neuron death
#'positive regulation of neuron apoptotic process (GO:0043525)', 'positive regulation of neuron death (GO:1901216)', \

