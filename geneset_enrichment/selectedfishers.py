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
pval = 1E-5

##p val cutoff for combos
oldpval = 1E-18
    
sstxs = ['M_2.5','M_5','M_10','M_15','T_5','T_10','T_20','T_25']
dirs = ['up','down']

dfs = {}
sets = {}
for tx in sstxs:
    df = pd.read_csv('../gene_expression/self_synergy_data/ss_limma_eb_'+tx+'.csv',
                     header=0,index_col='Unnamed: 0')
    dfs['up_'+tx] = df[(df['adj.P.Val'] <pval)&(df['logFC']>0)]
    dfs['down_'+tx] = df[(df['adj.P.Val'] <pval)&(df['logFC']<0)]

    sets['up_'+tx] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']>0)].index)
    sets['down_'+tx] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']<0)].index)

for di in dirs:
    sets[di+'_T5uT20'] = sets[di+'_T_5'].union(sets[di+'_T_20'])
    sets[di+'_M5uM10'] = sets[di+'_M_5'].union(sets[di+'_M_10'])

ctxs = ['M','T','TM']
hours = ['0','3','6','9','12','24']
for hour in hours:
    for tx in ctxs:
        df = pd.read_csv('../gene_expression/diff_limma_december/limma_t0_'+tx+'_'+hour+'.txt',
                         header=0,sep='\t',index_col='1')
        dfs['up_'+tx+'_'+hour] = df[(df['adj.P.Val'] <oldpval)&(df['logFC']>0)]
        dfs['down_'+tx+'_'+hour] = df[(df['adj.P.Val'] <oldpval)&(df['logFC']<0)]
        sets['up_'+tx+'_'+hour] = set(df[(df['adj.P.Val'] <oldpval)&(df['logFC']>0)].index)
        sets['down_'+tx+'_'+hour] = set(df[(df['adj.P.Val'] <oldpval)&(df['logFC']<0)].index)

for hour in hours:
    for di in dirs:
        sets[di+'_TUM_'+hour] = sets[di+'_T_'+hour]|sets[di+'_M_'+hour]
        
txs = sstxs[:3]+['M5uM10']+sstxs[3:7]+['T5uT20']+sstxs[7:]+\
[tx+'_'+hour for tx in ['M','T'] for hour in hours]+['TUM_'+hour for hour in hours]+\
['TM_'+hour for hour in hours]

allRNAseq = set(pd.read_csv('../gene_expression/self_synergy_data/ss_limma_eb_M_2.5.csv',header=0,index_col='Unnamed: 0').index)

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

cand = cand[cand.index.isin(finalsets)]

#calculate fishers
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
#fdrs = fdrs[fdrs[fdrs<0.05].any(axis=1)]
print(fdrs.shape)
#fdrs.to_csv('all_geneset_enrichment_fdr.tsv',sep='\t')
logged = -np.log10(fdrs)
#logged.to_csv('all_geneset_enrichment_fdr_log_.tsv',sep='\t')


cols=[t + '_' + h for t in ['T','M','TUM','TM'] for h in [
        '3','6','9','12','24']] + \
        ['T_5','T_10','T_20','T5uT20','T_25','M_10','M5uM10','M_15']
        
y = [x/255 for x in [249, 255, 60]]
b = [x/255 for x in [86, 112, 255]]
r = [x/255 for x in [255, 9, 60]]
g = [x/255 for x in [0,255,79]]
o = [x/255 for x in [255,170,0]]
m = [x/255 for x in [255,0,255]]

def make_rgb_transparent(color, alpha):
    rgb = colors.colorConverter.to_rgb(color)
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, (1,1,1))]
    
clist = [make_rgb_transparent(b,0.8)]*5+[make_rgb_transparent(r,2/3)]*5+\
    [make_rgb_transparent(b,0.8)]*5+[m]*5 +\
    [make_rgb_transparent(b,a) for a in [x/25 for x in [5,10,20,5,25]]] +\
    [make_rgb_transparent(r,a) for a in [x/15 for x in [10,5,15]]]


cols = ['down_'+tx for tx in cols]+['up_'+tx for tx in cols]

lut = dict(zip(cols,clist*2))

mat = logged.reindex(index=finalsets[::-1])[cols]

hm = sns.clustermap(mat,
    cmap='hot_r',col_cluster=False,
                    row_cluster=False,yticklabels=True,xticklabels=True,
                    cbar_kws={'label':'-log10(fdr)'},
                    vmax=6,vmin=0,figsize=(11,6),col_colors=mat.columns.map(lut))
hm.ax_heatmap.set_xticklabels([col.replace('_',' ') for col in mat.columns],
                               fontsize=6)
fig = hm.fig
#fig.set_xticklabels(ontsize=5)
#fig
#.xticks(fontsize=2)
fig.subplots_adjust(left=0.01,bottom=0.15,right=0.5,top=0.95)
fig.savefig('combo_self_synergy_fishers.pdf')


####all this code is about finding synergistic processes.
###For now just check the processes identified in teh TM combo previously.
##upcols= []
##downcols = []
##uupcols = []
##udowncols = []
##tmcols = []
##for col in fdrs.columns:
##    if any(st in col for st in ['_TM_','_T_','_M_','_TUM_']):
##        tmcols.append(col)
##    if '_TM_' in col:
##        if 'up_' in col:
##            upcols.append(col)
##        else:
##            downcols.append(col)
##    elif '_TUM_'in col:
##        if 'up_' in col:
##            uupcols.append(col)
##        else:
##            udowncols.append(col)
##tmup = fdrs[upcols]
##tmdown = fdrs[downcols]
##uup = fdrs[uupcols]
##udown = fdrs[udowncols]
##tm = fdrs[tmcols]
##
###tmup[tmup[tmup<0.00001].any(axis=1)]
###tmup[(tmup[tmup<0.00001].any(axis=1))&((tmup<0.01).all(axis=1))] --> 23 processes
###tmdown[(tmdown[tmdown<0.00001].any(axis=1))&((tmdown<0.01).all(axis=1))]
##
####unexpected gene sets
##sig = fdrs[(((tmup[tmup<0.00001].any(axis=1))&((tmup<0.01).all(axis=1)))&
##     ~((uup[uup<0.00001].any(axis=1))&((uup<0.01).all(axis=1))))| #synergistic
##           (((tmdown[tmdown<0.00001].any(axis=1))&((tmdown<0.01).all(axis=1)))&
##     ~((udown[udown<0.00001].any(axis=1))&((udown<0.01).all(axis=1))))| #synergistic
##     ((uup[uup<0.00001].any(axis=1))&((uup<0.01).all(axis=1)))|#additive
##           ((udown[udown<0.00001].any(axis=1))&((udown<0.01).all(axis=1)))] #additive
##
####drug targets
##targets = 'autophagy|BHAT_ESR1_TARGETS|estrogen'
##tar = fdrs[fdrs.index.str.contains(targets)&((tmup[tmup<0.01].any(axis=1))|(tmdown[tmdown<0.01].any(axis=1)))]
##
####hallmarks of cancer
###selected from each high level GO term associated with a hallmark of cancer (weinberg update), and 1-2 levels of its children connected by "is a" 
##candidatesets = [
##'GO:0001525', #angiogenesis
##'GO:0001837', #emt
##'GO:0006955', #immune response
##'GO:0045087|GO:0006959|GO:0002250|GO:0002418', #innate humoral adaptive | response to tumor
##'GO:0008219|GO:0012501|GO:0070265|GO:0036473', #cell death | programmed cell death| necrotic cell death | cell death in response to oxidative stress  
##'GO:0097190|GO:0097191|GO:0097193|GO:0070059|GO:0006915|GO:0043276', #apoptotic signaling pathway | extrinsic | intrinsic | response to ER stress | apoptotic process | anoikis
##'GO:0001896|GO:0070269|GO:0097468|GO:0097707|GO:0048102|GO:0097300', #autolysis|pyroptosis|response to ROS|ferroptosis|cornification|autophagy|programmed necrosis
##'GO:0007569|GO:0001302|GO:0090399|GO:0090398|GO:0001300', #cell aging | replicative cell aging | replicative senscence | cellular senescence | chronological cell aging
##'GO:0008283|GO:0050673', #cell proliferation | epithelial cell proliferation
##'GO:0007049|GO:0000278', #cell cycle | mitotic cell cycle
##'GO:0040007|GO:0045927|GO:0045926', #growth positive negative
##'GO:0006954|GO:0002544|GO:0002437|GO:0090594|GO:0002526', ##immune reponse | chronic | antigen| wounding | acute
##'GO:0008152|GO:0044237',  #metabolism | cellular | children
##'GO:0034641|GO:0006575|GO:0009850|GO:0072592|GO:0017144|GO:0006081|GO:0090345|GO:0006793|GO:0044107|GO:0044260|\
##GO:0097354|GO:0043446|GO:0044248|GO:0035383|GO:0044249|GO:0046950|GO:0009404|GO:0046483|GO:0006805|GO:0018942|\
##GO:0043449|GO:0006091|GO:0010191|GO:0072593|GO:0051186|GO:0031323|GO:0008218|GO:0042133|GO:0044255|GO:0006790|\
##GO:0006730|GO:0006725|GO:0006413|GO:0006082|GO:0042810|GO:0006127|GO:0042180|GO:2001290|GO:0001887|GO:0034754|\
##GO:0044262|GO:0043452']
##
###GO genesets
##sets = ''
##for i in candidatesets:
##    sets +=i
##    sets +='|'
##sets = sets.strip('|')
##
#####kegg genesets
####k = 'angiogenesis|EMT|epithelial to mesenchymal transition|immune|inflammation|cell death|apoptosis|autophagy|proliferation|senescence|growth|metabolism|invasion|metastasis'
##
##c = fdrs[(fdrs.index.str.contains(sets))&((tmup[tmup<0.01].any(axis=1))|(tmdown[tmdown<0.01].any(axis=1)))]
##
##final = pd.concat([sig,tar,c]).drop_duplicates()
##
##final.to_csv('hit_geneset_enrichment_fdr.tsv',sep='\t')
##finallogged = -np.log10(final)
##finallogged.to_csv('hit_geneset_enrichment_fdr_log.tsv',sep='\t')

#
#
#
#TM = pd.read_csv('../pathways_unbiased_newest/TM.tsv',sep='\t',
#                 index_col='candset',header=0,usecols=['candset',
#                                                       'down_M_24','down_T_24',
#                                                       'down_TUM_24','down_TM_24',
#                                                       'up_M_24','up_T_24',
#                                                       'up_TUM_24','up_TM_24'])
#
#finallogged = logged.T[finalsets].T
#
#finallogged= finallogged.merge(TM,how='outer',left_index=True,right_index=True)
#
#mcols = ['down_M_2.5', 'down_M_5', 'down_M_10', 'down_M5uM10', 'down_M_15', 
#         'down_M_24','down_TUM_24', 'down_TM_24',
#
#        'up_M_2.5', 'up_M_5', 'up_M_10', 'up_M5uM10', 'up_M_15', 
#        'up_M_24','up_TUM_24', 'up_TM_24']
#
#finallogged[mcols].to_csv('M_allcandidates.tsv',sep='\t')
#
#tcols = ['down_T_5', 'down_T_10', 'down_T_20', 'down_T5uT20', 'down_T_25', 
#         'down_T_24','down_TUM_24', 'down_TM_24',
#        
#        'up_T_5', 'up_T_10', 'up_T_20', 'up_T5uT20', 'up_T_25', 
#        'up_T_24','up_TUM_24', 'up_TM_24']
#
#finallogged[tcols].to_csv('T_allcandidates.tsv',sep='\t')
#
##excluded: things not relevant to breast cancer cells i.e. photosynthesis, neuron death
##'positive regulation of neuron apoptotic process (GO:0043525)', 'positive regulation of neuron death (GO:1901216)', \
#
