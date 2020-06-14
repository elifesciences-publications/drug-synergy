#check diff exp lists against all gene sets of selected libraries
#oct 15 2016

import pandas as pd
from itertools import product
from scipy import stats
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests as mult
from matplotlib import colors
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
    

#collect diff exp sets
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


#background
allRNAseq = set(pd.read_csv('../gene_expression/diff_limma_december/limma_t0_M_0.txt',header=0,
sep='\t',index_col='1').index)

##set and dysregulated genes list minimum length
minlist = 3
##set minimum overlap 
minov = 2


d = {'cc2015':'GO_Cellular_Component_2015'#, 
     #'mf2015':'GO_Molecular_Function_2015'
     }


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

for key in d.keys():
    
    name = d[key]
    
    #collect gene set libraries of interest
    go = pd.read_csv('../libraries_from_enrichr/'+name+'.txt',sep='\t',
                     names = list(range(10000)),index_col=0)
    go.drop(1,inplace=True,axis=1)

    go = pd.concat([go,tfs,pld])
#
#    
    rslt = pd.DataFrame(list(product(sets.keys(),go.index)), columns=['dxset', 'candset'])
    rslt['pval'] = rslt.apply(fishers,axis=1)
    nan = rslt[np.isnan(rslt.pval)]
    nan['fdr'] = np.nan
    notnan = rslt[~np.isnan(rslt.pval)]
    notnan['fdr'] = mult(notnan.pval,method='fdr_bh')[1]
    rslt = pd.concat([nan,notnan])
    rslt.reset_index(drop=True,inplace=True)
#    
    fdrs = rslt.pivot('candset','dxset','fdr')
    print(fdrs.shape)
    fdrs.dropna(how='all',inplace=True)
    print(fdrs.shape)
    fdrs.fillna(1.,inplace=True)
    fdrs = fdrs[fdrs[fdrs<0.05].any(axis=1)]
    print(fdrs.shape)
    fdrs.to_csv('all_geneset_enrichment_fdr_'+key+'.tsv',sep='\t')
    logged = -np.log10(fdrs)
    logged.to_csv('all_geneset_enrichment_fdr_log_'+key+'.tsv',sep='\t')
    
    
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
    
    
    #these are the candidate gene sets for biological processes
    ###drug targets
    #targets = 'autophagy|BHAT_ESR1_TARGETS|estrogen'
    #tar = fdrs[fdrs.index.str.contains(targets)&((tmup[tmup<0.01].any(axis=1))|(tmdown[tmdown<0.01].any(axis=1)))]
    #
    ###hallmarks of cancer
    ##selected from each high level GO term associated with a hallmark of cancer (weinberg update), and 1-2 levels of its children connected by "is a" 
    #candidatesets = [
    #'GO:0001525', #angiogenesis
    #'GO:0001837', #emt
    #'GO:0006955', #immune response
    #'GO:0045087|GO:0006959|GO:0002250|GO:0002418', #innate humoral adaptive | response to tumor
    #'GO:0008219|GO:0012501|GO:0070265|GO:0036473', #cell death | programmed cell death| necrotic cell death | cell death in response to oxidative stress  
    #'GO:0097190|GO:0097191|GO:0097193|GO:0070059|GO:0006915|GO:0043276', #apoptotic signaling pathway | extrinsic | intrinsic | response to ER stress | apoptotic process | anoikis
    #'GO:0001896|GO:0070269|GO:0097468|GO:0097707|GO:0048102|GO:0097300', #autolysis|pyroptosis|response to ROS|ferroptosis|cornification|autophagy|programmed necrosis
    #'GO:0007569|GO:0001302|GO:0090399|GO:0090398|GO:0001300', #cell aging | replicative cell aging | replicative senscence | cellular senescence | chronological cell aging
    #'GO:0008283|GO:0050673', #cell proliferation | epithelial cell proliferation
    #'GO:0007049|GO:0000278', #cell cycle | mitotic cell cycle
    #'GO:0040007|GO:0045927|GO:0045926', #growth positive negative
    #'GO:0006954|GO:0002544|GO:0002437|GO:0090594|GO:0002526', ##immune reponse | chronic | antigen| wounding | acute
    #'GO:0008152|GO:0044237',  #metabolism | cellular | children
    #'GO:0034641|GO:0006575|GO:0009850|GO:0072592|GO:0017144|GO:0006081|GO:0090345|GO:0006793|GO:0044107|GO:0044260|\
    #GO:0097354|GO:0043446|GO:0044248|GO:0035383|GO:0044249|GO:0046950|GO:0009404|GO:0046483|GO:0006805|GO:0018942|\
    #GO:0043449|GO:0006091|GO:0010191|GO:0072593|GO:0051186|GO:0031323|GO:0008218|GO:0042133|GO:0044255|GO:0006790|\
    #GO:0006730|GO:0006725|GO:0006413|GO:0006082|GO:0042810|GO:0006127|GO:0042180|GO:2001290|GO:0001887|GO:0034754|\
    #GO:0044262|GO:0043452']
    #
    ##GO genesets
    #sets = ''
    #for i in candidatesets:
    #    sets +=i
    #    sets +='|'
    #sets = sets.strip('|')
    #
    ####kegg genesets
    ###k = 'angiogenesis|EMT|epithelial to mesenchymal transition|immune|inflammation|cell death|apoptosis|autophagy|proliferation|senescence|growth|metabolism|invasion|metastasis'
    
    
    #I'm going to do something slightly different here than i did for biological process:
    #I'm going to include additivity in the candidate analysis
    #anything significnat in t,m,tum, or tm
    c = fdrs[(fdrs.index.str.contains(setstr)|fdrs.index.isin(candsets))&
             (fdrs[tmcols][fdrs[tmcols]<0.01]).any(axis=1)]
    
    final = pd.concat([sig,c]).drop_duplicates()
    
    final.to_csv('hit_geneset_enrichment_fdr_'+key+'.tsv',sep='\t')
    finallogged = -np.log10(final)
    finallogged.to_csv('hit_geneset_enrichment_fdr_log_'+key+'.tsv',sep='\t')
    
#    finallogged = pd.read_csv('hit_geneset_enrichment_fdr_log_'+key+'.tsv',sep='\t')
    
    #order
    order = pd.DataFrame([finallogged[upcols].max(axis=1)-finallogged[uupcols].max(axis=1),
                          finallogged[downcols].max(axis=1)-finallogged[udowncols].max(
                                  axis=1)]).T.max(axis=1).sort_values()
    fo = finallogged+0.01
    order = pd.DataFrame([fo[upcols].max(axis=1)-fo[uupcols].max(
            axis=1)/fo[uupcols].max(axis=1),
                          fo[downcols].max(axis=1)-fo[udowncols].max(
                                  axis=1)/fo[udowncols].max(
                                  axis=1)]).T.max(axis=1).sort_values()

    cols=[t + '_' + h for t in ['T','M','TUM','TM','W','MUW','MW','TUW','TW'] for h in [
            '3','6','9','12','24']] 
            
    y = 'gold' #[x/255 for x in [249, 255, 60]]
    b = 'b' #[x/255 for x in [86, 112, 255]]
    r = 'r' #[x/255 for x in [255, 9, 60]]
    g = 'g' #[x/255 for x in [0,255,79]]
    o = 'orange' #[x/255 for x in [255,170,0]]
    m = [x/255 for x in [255,0,255]]
    
    def make_rgb_transparent(color, alpha):
        rgb = colors.colorConverter.to_rgb(color)
        return [alpha * c1 + (1 - alpha) * c2
                for (c1, c2) in zip(rgb, (1,1,1))]
        
    clist = [make_rgb_transparent(b,0.8)]*5+[make_rgb_transparent(r,0.8)]*5+\
        [make_rgb_transparent(b,0.8)]*5+[m]*5 + [make_rgb_transparent(y,0.8)]*5+\
        [make_rgb_transparent(r,0.8)]*5+[make_rgb_transparent(o,0.8)]*5+\
        [make_rgb_transparent(b,0.8)]*5+[make_rgb_transparent(g,0.8)]*5
    
    
    cols = ['down_'+tx for tx in cols]+['up_'+tx for tx in cols]
    
    lut = dict(zip(cols,clist*2))
    
    mat = finallogged.reindex(index=order.sort_values().index)[cols]
#    mat2 = finallogged[udowncols+downcols+uupcols+upcols]
#    
#    hm2 = sns.clustermap(mat2,
#        cmap='hot_r',col_cluster=False,
#                        row_cluster=True,yticklabels=True,xticklabels=True,
#                        cbar_kws={'label':'-log10(fdr)'},
#                        vmax=6,vmin=0,figsize=(11,6),col_colors=mat.columns.map(lut))
#    
    hm = sns.clustermap(mat,
        cmap='hot_r',col_cluster=False,
                        row_cluster=True,yticklabels=True,xticklabels=True,
                        cbar_kws={'label':'-log10(fdr)'},
                        vmax=6,vmin=0,figsize=(11,6),col_colors=mat.columns.map(lut))#vmax=6
    hm.ax_heatmap.set_xticklabels([col.replace('_',' ') for col in mat.columns],
                                   fontsize=6)
    fig = hm.fig
    #fig.set_xticklabels(ontsize=5)
    #fig
    #.xticks(fontsize=2)
    fig.subplots_adjust(left=0.01,bottom=0.15,right=0.62,top=0.95)
    fig.savefig(key+'_fishers.pdf')
    
    
    
    #tmcols = ['down_M_3', 'down_M_6', 'down_M_9', 'down_M_12', 'down_M_24',
    #        'down_T_3', 'down_T_6', 'down_T_9', 'down_T_12', 'down_T_24',
    #        'down_TUM_3', 'down_TUM_6','down_TUM_9', 'down_TUM_12', 'down_TUM_24',
    #        'down_TM_3', 'down_TM_6', 'down_TM_9','down_TM_12', 'down_TM_24',
    #        
    #        'up_M_3', 'up_M_6', 'up_M_9', 'up_M_12', 'up_M_24',
    #        'up_T_3', 'up_T_6', 'up_T_9', 'up_T_12', 'up_T_24',
    #        'up_TUM_3', 'up_TUM_6', 'up_TUM_9', 'up_TUM_12', 'up_TUM_24', 
    #        'up_TM_3', 'up_TM_6', 'up_TM_9','up_TM_12', 'up_TM_24']
    #
    #finallogged[tmcols].to_csv('TM.tsv',sep='\t')
    #
    #mwcols = ['down_M_3', 'down_M_6', 'down_M_9', 'down_M_12', 'down_M_24',
    #          'down_W_3', 'down_W_6', 'down_W_9', 'down_W_12', 'down_W_24',
    #        'down_MUW_3', 'down_MUW_6','down_MUW_9', 'down_MUW_12', 'down_MUW_24',
    #        'down_MW_3', 'down_MW_6', 'down_MW_9','down_MW_12', 'down_MW_24',
    #        
    #        'up_M_3', 'up_M_6', 'up_M_9', 'up_M_12', 'up_M_24',
    #          'up_W_3', 'up_W_6', 'up_W_9', 'up_W_12', 'up_W_24',
    #        'up_MUW_3', 'up_MUW_6', 'up_MUW_9', 'up_MUW_12', 'up_MUW_24', 
    #        'up_MW_3', 'up_MW_6', 'up_MW_9','up_MW_12', 'up_MW_24']
    #
    #finallogged[mwcols].to_csv('MW.tsv',sep='\t')
    #
    #twcols = [
    #        'down_T_3', 'down_T_6', 'down_T_9', 'down_T_12', 'down_T_24',
    #        'down_W_3', 'down_W_6', 'down_W_9', 'down_W_12', 'down_W_24',
    #        'down_TUW_3', 'down_TUW_6','down_TUW_9', 'down_TUW_12', 'down_TUW_24',
    #        'down_TW_3', 'down_TW_6', 'down_TW_9','down_TW_12', 'down_TW_24',
    #        
    #        'up_T_3', 'up_T_6', 'up_T_9', 'up_T_12', 'up_T_24',
    #        'up_W_3', 'up_W_6', 'up_W_9', 'up_W_12', 'up_W_24',
    #        'up_TUW_3', 'up_TUW_6', 'up_TUW_9', 'up_TUW_12', 'up_TUW_24', 
    #        'up_TW_3', 'up_TW_6', 'up_TW_9','up_TW_12', 'up_TW_24']
    #
    #finallogged[twcols].to_csv('TW.tsv',sep='\t')
    #
    ##excluded: things not relevant to breast cancer cells i.e. photosynthesis, neuron death
    ##'positive regulation of neuron apoptotic process (GO:0043525)', 'positive regulation of neuron death (GO:1901216)', \
    #

#def inter(row):
#    de = sets[row.dxset]
#    s = set(go.ix[row.candset].dropna()).intersection(allRNAseq)
#    notde = allRNAseq - de
#    nots = allRNAseq - s
#    ##contingency table: 
#            ## [[diffex&inset,diffex&notinset,
#            ## [notdiffex&inset,notdiffex&notinset]]
#    if len(de)>=minlist and len(s)>=minlist and len(de.intersection(s))>=minov:
#        return len(de.intersection(s))
#    else:
#        return None
#    
#interrslt = pd.DataFrame(list(product(sets.keys(),go.index)), columns=['dxset', 'candset'])
#interrslt['intersection'] = interrslt.apply(inter,axis=1)
#interrslt2 = interrslt.pivot('candset','dxset','intersection')