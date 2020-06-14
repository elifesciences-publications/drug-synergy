#check diff exp lists against all gene sets of selected libraries
#oct 15 2016

import pandas as pd
from itertools import product
from scipy import stats
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mult
from matplotlib import colors
import seaborn as sns
import matplotlib.pyplot as plt

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


#candidate gene sets for molecular function and cellular component
#from lysosomotropism paper
candsets = ['pre-autophagosomal structure','cytoplasmic vesicle',
            'mitochondrial matrix','MCM complex','mitochondrial nucleoid',
            'smooth endoplasmic reticulum','endosome membrane',
            'autophagic vacuole membrane','COPII vesicle coat',
            'Fanconi anaemia nuclear complex','ER to Golgi transport vesicle',
            'mitochondrial intermembrane space','pernuclear region of cytoplasm',
            'nuclear speck']
#the gene term 'cytoplasmic vesicle' is not in this library. probably Avi's lab
#removed it because it's too big -- looks like it has 6433 human genes
#instead the ones I can include are cytoplasmic vesicle membrane and part

setstr = ''
for i in candsets:
    setstr += '|'
    setstr += i
    
setstr = setstr[1:]

candsets = ['lysosome (GO:0005764)', 'nucleolus (GO:0005730)',
            'Golgi membrane (GO:0000139)','nuclear pore (GO:0005643)',
            'spliceosomal complex (GO:0005681)','endosome (GO:0005768)']

#lysosomotropism transcriptional targets
reg = pd.read_csv('../network_new/MCF7_Woo_et_al-Li_et_al_genelevel_noregulatorsastargets_20160126_JELD.txt',
                  sep='\t')
#of interest are 'TFEB' (too few targets in our network) and 'TEF3'
ls = []

for txtype in ['Transcriptional Activator','Transcriptional Repressor']:
    tmp = ['TFEB_TFE3_'+txtype+'s','']
    for regulator in ['TFEB','TFE3']:
        tmp+=reg[(reg.Regulator==regulator)&(
                        reg.Type==txtype)].Target.tolist()
    ls.append(tmp)
tfs = pd.DataFrame(ls).set_index(0).drop(1,axis=1)

candsets += tfs.index.tolist()

#PLD signature
pld = pd.read_excel('../lysosomotopism/PLD_signature.xlsx',sheet_name='Sheet1')
ls = []
for col in pld.columns:
    ls.append( ['PLD_'+col.split(' ')[1].split('R')[0],'']+pld[col].str.split(' - ').str[0].tolist())
pld = pd.DataFrame(ls).set_index(0).drop(1,axis=1)

candsets+=pld.index.tolist()

name = 'GO_Cellular_Component_2015'
key = 'cc2015'


#collect gene set libraries of interest
go = pd.read_csv('../libraries_from_enrichr/'+name+'.txt',sep='\t',
                 names = list(range(10000)),index_col=0)
go.drop(1,inplace=True,axis=1)

go = pd.concat([go,tfs,pld])


finalsets = [#
        'ER to Golgi transport vesicle (GO:0030134)',
        'Golgi membrane (GO:0000139)',
        
        
        #ER
        'endoplasmic reticulum membrane (GO:0005789)',
        'mitochondrial matrix (GO:0005759)',
        'extracellular vesicular exosome (GO:0070062)',
    
        
        'lysosome (GO:0005764)',
        'cytoplasmic vesicle membrane (GO:0030659)',
        'cytoplasmic vesicle part (GO:0044433)', 
        
        'pre-autophagosomal structure (GO:0000407)',
        'pre-autophagosomal structure membrane (GO:0034045)',
        
        
        'COPII vesicle coat (GO:0030127)',
        'TFEB_TFE3_Transcriptional Activators',
        
        'nuclear pore (GO:0005643)',
        'spliceosomal complex (GO:0005681)',
        'nucleolus (GO:0005730)',
        'nucleoplasm (GO:0005654)',
        'MCM complex (GO:0042555)',
        
        'PLD_Down', 'PLD_Up', 
       ]
           
 

cand = go[go.index.isin(finalsets)]

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
fdrs.to_csv('all_geneset_enrichment_fdr_'+key+'.tsv',sep='\t')
logged = -np.log10(fdrs)
logged.to_csv('all_geneset_enrichment_fdr_log_'+key+'.tsv',sep='\t')


cols=[t + '_' + h for t in ['T','M','TUM','TM'] for h in [
        '3','6','9','12','24']] + \
        ['T_5','T_10','T_20','T5uT20','T_25','M_10','M5uM10','M_15']
        
y = [x/255 for x in [249, 255, 60]]
b = 'b'#[x/255 for x in [86, 112, 255]]
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

mat = logged.reindex(index=finalsets)[cols] #[::-1]

plt.rcParams['font.family'] = "Arial" #Narrow
plt.rcParams['font.size'] = 6

hm = sns.clustermap(mat,
    cmap='hot_r',col_cluster=False,
                    row_cluster=False,yticklabels=True,xticklabels=True,
                    cbar_kws={'label':'-log10(fdr)'},
                    vmax=50,#6
                    vmin=0,figsize=(11,6),col_colors=mat.columns.map(lut))
hm.ax_heatmap.set_xticklabels([col.replace('_',' ') for col in mat.columns],
                               fontsize=6)
fig = hm.fig
#fig.set_xticklabels(ontsize=5)
#fig
#.xticks(fontsize=2)
fig.subplots_adjust(left=0.01,bottom=0.5,right=0.7,top=0.95)
fig.savefig('combo_self_synergy_fishers_'+key+'_highvmax.pdf') #


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
