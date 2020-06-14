#check diff exp lists against all gene sets of selected libraries
#oct 15 2016

import pandas as pd
from itertools import product
from scipy import stats
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mult
import seaborn as sns

##finallogged = pd.read_csv('hit_geneset_enrichment_fdr_log.tsv',sep='\t',index_col = 'candset')#
##print('read final file from previous run')

def fishers(row):
    de = sets[row.dxset]
    s = set(go.ix[row.candset].dropna()).intersection(allRNAseq)
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


sns.reset_defaults()
sns.set(font='Arial',font_scale=0.5,color_codes=False)

##set p value cutoff for diff exp
pval = 1E-15
    
txs = ['T','M','TM','W','TW','MW'] #
times = ['3','6','9','12','24']
dirs = ['up','down']

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

dfs = {}
sets = {}
for tx in txs:
    for time in times:
        df = pd.read_csv('../gene_expression/LNCaP/LIMMA_DATA/DE_limma_'+tx+'_'+time+'.csv',
                         header=0,index_col='Unnamed: 0')
        dfs[tx+'_'+time] = df
        sets[tx+'_'+time] = set(df[df['adj.P.Val']<pval].index)
        sets['up_'+tx+'_'+time] = set(df[(df['adj.P.Val']<pval)&(df.logFC>0)].index)
        sets['down_'+tx+'_'+time] = set(df[(df['adj.P.Val']<pval)&(df.logFC<0)].index)

for time in times:
    for di in dirs:
        sets[di+'_TUM_'+time] = sets[di+'_T_'+time].union(sets[di+'_M_'+time])
        sets[di+'_TUW_'+time] = sets[di+'_T_'+time].union(sets[di+'_W_'+time])
        sets[di+'_MUW_'+time] = sets[di+'_M_'+time].union(sets[di+'_W_'+time])


allRNAseq = set(pd.read_csv('../gene_expression/LNCaP/LIMMA_DATA/DE_limma_M_3.csv',
                            header=0,index_col='Unnamed: 0').index)

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



rslt = pd.DataFrame(list(product(sets.keys(),go.index)), columns=['dxset', 'candset'])
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
fdrs.to_csv('all_geneset_enrichment_fdr_cc2015.tsv',sep='\t')
logged = -np.log10(fdrs)
logged.to_csv('all_geneset_enrichment_fdr_log_cc2015.tsv',sep='\t')


###read stuff back in after the first run
#fdrs = pd.read_csv('all_geneset_enrichment_fdr.tsv',sep='\t',
#                     index_col='candset')
#logged = pd.read_csv('all_geneset_enrichment_fdr_log_.tsv',sep='\t',
#                     index_col='candset')


final = pd.DataFrame()
for c in ['TM','TW','MW']:
    #set up combo
    upcols= ['up_'+c+'_'+time for time in times]
    downcols = ['down_'+c+'_'+time for time in times]
    uupcols = ['up_'+c[0]+'U'+c[1]+'_'+time for time in times]
    udowncols = ['down_'+c[0]+'U'+c[1]+'_'+time for time in times]
    ccols = ['up_'+c[0]+'_'+time for time in times]+\
        ['up_'+c[1]+'_'+time for time in times]+\
        ['down_'+c[0]+'_'+time for time in times]+\
        ['down_'+c[1]+'_'+time for time in times]+\
        upcols+downcols+uupcols+udowncols
    
    cup = fdrs[upcols]
    cdown = fdrs[downcols]
    uup = fdrs[uupcols]
    udown = fdrs[udowncols]
    com = fdrs[ccols]
    
    #I'm going to do something slightly different here than i did for biological process:
    #I'm going to include additivity in the candidate analysis
    #anything significnat in t,m,tum, or tm
    c = fdrs[(fdrs.index.str.contains(setstr)|fdrs.index.isin(candsets))&
             (fdrs[ccols][fdrs[ccols]<0.01]).any(axis=1)]
    
    ##unexpected gene sets
    sig = fdrs[(((cup[cup<0.00001].any(axis=1))&((cup<0.01).all(axis=1)))&
         ~((uup[uup<0.00001].any(axis=1))&((uup<0.01).all(axis=1))))| #synergistic
               (((cdown[cdown<0.00001].any(axis=1))&((cdown<0.01).all(axis=1)))&
         ~((udown[udown<0.00001].any(axis=1))&((udown<0.01).all(axis=1))))| #synergistic
         ((uup[uup<0.00001].any(axis=1))&((uup<0.01).all(axis=1)))|#additive
               ((udown[udown<0.00001].any(axis=1))&((udown<0.01).all(axis=1)))] #additive
    
    final = pd.concat([final,sig,c]).drop_duplicates()

final.to_csv('hit_geneset_enrichment_fdr_3combos_cc2015.tsv',sep='\t')
finallogged = -np.log10(final)
finallogged.to_csv('hit_geneset_enrichment_fdr_log_3combos_cc2015.tsv',sep='\t')

#finallogged = pd.read_csv('hit_geneset_enrichment_fdr_log_3combos.tsv',sep='\t',
#                          index_col='candset')
newtxs = ['T', 'M', 'TUM','TM', 'W', 'TUW','TW', 'MUW','MW']
finallogged = finallogged[['down_'+tx+'_'+time for tx in newtxs for time in times] +\
    ['up_'+tx+'_'+time for tx in newtxs for time in times]]
#seaborn heatmap
c = [el for ls in [[x]*5 for x in ['b','r','r','magenta', #'b','magenta'
                   'gold','gold','g','gold','orange']] for el in ls]*2 #'b','g','r','orange'


lut = dict(zip(finallogged.columns.tolist(),c))

hm=sns.clustermap(finallogged,col_cluster=False,cmap='hot_r',
                  vmax=50,#6
               yticklabels=True,cbar_kws={'label':'-log10(fdr)'},figsize=(12,6),
               col_colors=finallogged.columns.map(lut)) #row_cluster=False,

fig = hm.fig
fig.subplots_adjust(left=0.01,bottom=0.1,right=0.718,top=0.99)
fig.savefig('alltxs_cc2015_highvmax.pdf')


##use this line to look at the sets that showed up in MCF7
#finallogged = finallogged[finallogged.index.isin(finalsets)]

