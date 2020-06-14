#synergistic/diff expressed gene sets, for each combo
#jul 8 2015

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn as venn
import csv
from matplotlib import colors
from matplotlib.colors import ColorConverter
import numpy as np
import matplotlib
from scipy import stats

plt.ioff()

#access color dict for color mixing
c = ColorConverter.colors

#using matplotlib-venn colormixing similar to venn.__venn3.compute_colors
def mix(col1,col2,col3=None):
    #input a color as string
    if col3==None:
        return venn._common.mix_colors(np.array(c[col1]),np.array(c[col2]))
    else:
        return venn._common.mix_colors(np.array(c[col1]),np.array(c[col2]),np.array(c[col3]))
    
    
def make_rgb_transparent(color, alpha):
    rgb = colors.colorConverter.to_rgb(color)
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, (1,1,1))]


def Fishers(set1,set2,minlist,minov):
    
    not1 = allRNAseq.difference(set1)
    not2 = allRNAseq.difference(set2)
    
    if len(set1)>=minlist and len(set2)>=minlist and len(set1.intersection(set2))>=minov:
        cont = [[len(set1.intersection(set2)),len(set1.intersection(not2))], 
                [len(not1.intersection(set2)),len(not1.intersection(not2))]]
        p = stats.fisher_exact(cont)[1]
    else:
        p = np.nan
    return p


##set fonts
font = {'fontname':'Arial','size':10}
#matplotlib.rc('font', **font)

##set old p value cutoff
pval_old = 1E-18

##set new p value cutoff
pval = 1E-5

##import data
txs = ['M_2.5','M_5','M_10','M_15','T_5','T_10','T_20','T_25']

allRNAseq_old = set(pd.read_csv('../diff_limma_december/limma_t0_M_0.txt',header=0,
                            sep='\t',index_col='1').index)
print('all genes old',len(allRNAseq_old))

allRNAseq = set(pd.read_csv('ss_limma_eb_M_2.5.csv',header=0,index_col='Unnamed: 0').index)
print('all genes ss',len(allRNAseq))

denames = []
for tx in txs:
    denames.append('ss_limma_eb_'+tx+'.csv')

dfs = {}

for name in denames:
    df = pd.read_csv(name,header=0,index_col='Unnamed: 0')
    dfs[name.strip('sslimmaeb_.csv')] = df

times = [0,3,6,9,12,24]
for drug in ['T','M','TM','W','TW','MW']:
    for time in times:
        df = pd.read_csv('../diff_limma_december/limma_t0_'+drug+'_'+str(time)+'.txt',header=0,
                                sep='\t',index_col='1')
        dfs[drug+'_'+str(time)+'_triplicate']=df
    
##testing p values
test_dict = {'M_10':[],'T_20':[]}
    
for p in [10**-i for i in range(3,15)]:
    for key in test_dict.keys():
        df = dfs[key]
        test_dict[key].append(df[df['adj.P.Val'] <p].shape[0])

fig, ax = plt.subplots()
ax.scatter(test_dict['M_10'], test_dict['T_20'])
for i in range(3,15):
    ax.annotate(str(i), (test_dict['M_10'][i-3], test_dict['T_20'][i-3]))
ax.plot(dfs['M_24_triplicate'][dfs['M_24_triplicate']['adj.P.Val']<pval_old].shape[0],
        dfs['T_24_triplicate'][dfs['T_24_triplicate']['adj.P.Val']<pval_old].shape[0],'ro')
ax.set_xlabel('DE genes in Mefloquine')
ax.set_ylabel('DE genes in Tamoxifen')
ax.text(x = 500,y=1000,s='adj.P.Val = 1E-n')
ax.text(x = 200,y=400,s='M_10,T_20 in triplicate',color='r')
#plt.show()

#apoptosis genes
apop = pd.read_csv('../../libraries_from_enrichr/GO_Biological_Process.txt',sep='\t',
                 names = list(range(400)),index_col=0,dtype=str)
apop = apop[apop.index.isin(['intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress (GO:0070059)', \
            'intrinsic apoptotic signaling pathway (GO:0097193)', \
            'apoptotic signaling pathway (GO:0097190)',\
            'generation of precursor metabolites and energy (GO:0006091)'])].dropna(axis=1,how='all')

#build sets
sets = {}
for key in dfs.keys():
    df = dfs[key]
    if 'triplicate' in key:
        sets[key] = set(df[df['adj.P.Val'] <pval_old].index)
    else:
        sets[key] = set(df[df['adj.P.Val'] <pval].index)
    sets[key+'_up'] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']>0)].index)
    sets[key+'_down'] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']<0)].index)

##bar graph of synergistic genes
plt.figure()
plt.bar(x=list(range(6)),height=[len(sets['M_10'].difference(sets['M_5'])),
         len((sets['M_15'].difference(sets['M_10'])).difference(sets['M_5'])),
         len(sets['T_10'].difference(sets['T_5'])),
         len(sets['T_20'].difference(sets['T_10'])),
         len((sets['T_25'].difference(sets['T_20'])).difference(sets['T_5'])),
         len((sets['TM_24_triplicate'].difference(sets['T_24_triplicate'])).difference(sets['M_24_triplicate']))],
        align='center',
        tick_label=['M_10','M15','T10','T20','T25','TM'])
plt.ylabel('Number of synergistic genes')
#plt.show()

#get differential expression venn diagrams

######all diff expressed genes
##plt.figure()
##v = venn.venn3([sets['M_2.5'],sets['M_5'],sets['M_10']],
##set_labels = ('Mefloquine_2.5','Mefloquine_5','Mefloquine_10'))
####for i in ['001','100','010','110','011','101','111']:
####      if v.get_label_by_id(i) != None:
####          v.get_label_by_id(i).set_text('')
####for text in v.set_labels:
####    text.set_fontsize(46)
####    text.set_family('Arial')
###plt.title('Differentially Expressed Genes at '+time+' Hours')
###plt.show()
##plt.savefig('venn/M_2.5_5_10.pdf')

plt.figure()
v = venn.venn2([sets['M_10'],sets['M_15']],
set_labels = ('Mefloquine_10','Mefloquine_15'))
##for i in ['001','100','010','110','011','101','111']:
##      if v.get_label_by_id(i) != None:
##          v.get_label_by_id(i).set_text('')
##for text in v.set_labels:
##    text.set_fontsize(46)
##    text.set_family('Arial')
#plt.title('Differentially Expressed Genes at '+time+' Hours')
#plt.show()
plt.savefig('venn/M_10_15.png')

plt.figure()
v = venn.venn3([sets['T_5'],sets['T_10'],sets['T_20']],
set_labels = ('Tamoxifen_5','Tamoxifen_10','Tamoxifen_20'))
##for i in ['001','100','010','110','011','101','111']:
##      if v.get_label_by_id(i) != None:
##          v.get_label_by_id(i).set_text('')
##for text in v.set_labels:
##    text.set_fontsize(46)
##    text.set_family('Arial')
#plt.title('Differentially Expressed Genes at '+time+' Hours')
#plt.show()
plt.savefig('venn/T_5_10_20.png')

plt.figure()
v = venn.venn3([sets['T_5'],sets['T_20'],sets['T_25']],
set_labels = ('Tamoxifen_5','Tamoxifen_20','Tamoxifen_25'))
##for i in ['001','100','010','110','011','101','111']:
##      if v.get_label_by_id(i) != None:
##          v.get_label_by_id(i).set_text('')
##for text in v.set_labels:
##    text.set_fontsize(46)
##    text.set_family('Arial')
#plt.title('Differentially Expressed Genes at '+time+' Hours')
#plt.show()
plt.savefig('venn/T_5_20_25.png')

###plot overlaps between DE genes
plt.figure()
v = venn.venn2([sets['M_10'],sets['M_24_triplicate']],
               set_labels=('Mefloquine_10_duplicate','Mefloquine_10_triplicate'))
plt.savefig('venn/M_duplicate_vs_triplicate.png')
plt.figure()
v = venn.venn2([sets['T_20'],sets['T_24_triplicate']],
               set_labels=('Tamoxifen_20_duplicate','Tamoxifen_20_triplicate'))
plt.savefig('venn/T_duplicate_vs_triplicate.png')


#plot overlaps between synergistic genes
plt.figure()
v = venn.venn3(
    [(sets['M_15'].difference(sets['M_10'])),
     (sets['T_25'].difference(sets['T_20'])).difference(
         sets['T_5']),(sets['TM_24_triplicate'].difference(sets['T_24_triplicate'])).difference(
        sets['M_24_triplicate']),],
    set_labels = ('T_25_synergistic','M_15_synergistic','TM_synergistic'))
##    pt = Fishers((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
##            sets['limma_t0_M_'+time]),
##                 (sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
##            sets['limma_t0_W_'+time]),
##                 3,2)
##    pm = Fishers((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
##            sets['limma_t0_M_'+time]),
##                 (sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
##             sets['limma_t0_W_'+time]),
##                 3,2)
##    pw = Fishers((sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
##            sets['limma_t0_W_'+time]),
##                 (sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
##             sets['limma_t0_W_'+time]),
##                 3,2)
##    if v.get_label_by_id('110') != None:
##        v.get_label_by_id('110').set_text(str(v.get_label_by_id('110')).split(',')[2].strip(")'")+'\n'+str(pt))
##        v.get_label_by_id('101').set_text(str(v.get_label_by_id('101')).split(',')[2].strip(")'")+'\n'+str(pm))
##        v.get_label_by_id('011').set_text(str(v.get_label_by_id('011')).split(',')[2].strip(")'")+'\n'+str(pw)) 
plt.title('Synergistic genes')
plt.savefig('venn/synergistic.png')
##

#get diff exp genes in dataframe for fischer's exact tests 
for tx in txs:
    df = pd.DataFrame(data=None,columns = ['gene','direction'])
    print(tx)
    for direc in ['up','down']:
        newdf = pd.DataFrame(data={'gene':list(sets[tx+'_'+direc]),'direction':direc})
        df = pd.concat([df,newdf])
    print(df.head())
    df.to_csv('diffexp_selfsynergy_'+tx+'_JELD.csv',index=False)
##


##filter genes that don't vary and calculate correlations
##always = sets['M_10'].copy()
##for key in sets.keys():
##    if any(i in key for i in ['_0_','2.5','M_5','M_10']):
##        always = always.intersection(sets[key])
##
##print('always',len(always)) #0
##
##never = allRNAseq.copy()
##for key in sets.keys():
##    never = never.difference(sets[key])
##
##print('never',len(never))#6237 genes that are never DE --these don't make a big difference so don't use this filter

def corr(dfs,pval=1.1):
    corr_dict = {}
    p_corr_dict = {}
    for time in times:
        for tx in ['TM','TW','MW']:
            s = set(dfs[tx[0]+'_'+str(time)+'_triplicate'][dfs[tx[0]+'_'+str(time)+'_triplicate']['adj.P.Val']<pval].index).union(
                set(dfs[tx[1]+'_'+str(time)+'_triplicate'][dfs[tx[1]+'_'+str(time)+'_triplicate']['adj.P.Val']<pval].index))
            corr_dict[tx+'_'+str(time)+'_triplicate'],p_corr_dict[tx+'_'+str(time)+'_triplicate'] = stats.pearsonr(
                dfs[tx[0]+'_'+str(time)+'_triplicate'][dfs[tx[0]+'_'+str(time)+'_triplicate'].index.isin(s)].sort_index().logFC,
                dfs[tx[1]+'_'+str(time)+'_triplicate'][dfs[tx[1]+'_'+str(time)+'_triplicate'].index.isin(s)].sort_index().logFC)
    combos = [['M_10','T_20','M10.T20']]
    combos += [[pair[0],pair[1],
                pair[0][:2]+str(int(float(pair[0].split('_')[1])+
                                    float(pair[1].split('_')[1])))] for pair \
               in [['M_2.5']*2,['M_5']*2,['M_5','M_10'],['T_5']*2,['T_10']*2,['T_5','T_20']]]
    for combo in combos:
        s = set(dfs[combo[0]][dfs[combo[0]]['adj.P.Val']<pval].index).union(
            set(dfs[combo[1]][dfs[combo[1]]['adj.P.Val']<pval].index))
        corr_dict[combo[2]],p_corr_dict[combo[2]] = stats.pearsonr(
                dfs[combo[0]][dfs[combo[0]].index.isin(s)].sort_index().logFC,
                dfs[combo[1]][dfs[combo[1]].index.isin(s)].sort_index().logFC)
    return (corr_dict,p_corr_dict)

corr_dict,p_corr_dict = corr(dfs)

##filtered = {}
##for key in dfs.keys():
##    filtered[key] = dfs[key][~dfs[key].index.isin(never)]
    
new_corr_dict,new_p_corr_dict = corr(dfs,0.1)

ind = np.arange(len(corr_dict.keys()))
width= 0.35
fig, ax = plt.subplots()
p1=ax.bar(ind,corr_dict.values(),width,color='r')
p2 = ax.bar(ind+width,new_corr_dict.values(),width,color='b')
ax.set_title('correlation between monotherapies')
ax.set_ylabel('spearman r of logFC')
ax.set_xticks(ind+width/2)
ax.set_xticklabels([i.strip('_triplicate') for i in new_corr_dict.keys()])
ax.legend((p1[0],p2[0]),('All genes','DE genes'))
ax.tick_params(axis='x', rotation=60,pad=0)
ax.margins(x=0)
plt.savefig('correlations.png')
#plt.show()

##viability data
ss_phen = pd.read_csv('ss_phenotype.csv') #not using this anymore
##this sheet is not ideal
ss = pd.read_excel('../../../drug synergy project/viability/Self Synergy/May 2016 MCF7 Viability.xlsx',
                   sheet_name='MCF7 DRC 24Hrs CTG',header=0)
ss.dropna(how='all',inplace=True)
ss = ss[ss.Plate.str.contains('MCF7 24Hrs InCell')]
ss.groupby(['Cmpd','Concentration'])[['Viability']].mean()
ss_viab = pd.concat([ss.groupby(['Cmpd','Concentration'])[['Viability']].mean()*100,
                     ss.groupby(['Cmpd','Concentration'])[['Viability']].std()/np.sqrt(1.5)*100],axis=1)
ss_viab.columns = ['Viability','Std']
ss_viab['index']=['DMSO','M_2.5','M_5','M_10','M_15','T_5','T_10','T_20','T_25','T_30']
ss_viab.set_index('index',inplace=True)
#ss_viab = ss_phen.groupby('Treatment_Time')[['InCell']].mean()#.rename(columns={'InCell':'InCell_mean'})
#ss_viab = ss_viab.merge(ss_phen.groupby('Treatment_Time')[['InCell']].std())
#ss_viab['Viability']=100*ss_viab.InCell/ss_viab.ix['DMSO','InCell']
ss_eob = pd.DataFrame([['M_5','M_10','M_15','T_10','T_20','T_25'],
			  ['M_2.5','M_5','M_5','T_5','T_10','T_5'],
			  ['M_2.5','M_5','M_10','T_5','T_10','T_20']],
			  index=['combo','Drug1','Drug2']).T
ss_eob['Drug 1 Viability'] = ss_eob.Drug1.apply(lambda x: ss_viab.ix[x,'Viability'])
ss_eob['Drug 2 Viability'] = ss_eob.Drug2.apply(lambda x: ss_viab.ix[x,'Viability'])
ss_eob['Observed Viability'] = ss_eob.combo.apply(lambda x: ss_viab.ix[x,'Viability'])

combo_eob = pd.read_excel(
    '../../../drug synergy project/viability/IBM MCF7 Combination Viability February 2015 v3.xlsx',
    sheet_name='Combinations')
combo_eob = combo_eob[(combo_eob.Assay=='InCell')&(combo_eob.Hours<48)]
combo_eob.Drug1 = combo_eob.Drug1.str[0]+'_'+combo_eob.Hours.astype(int).astype(str)+'_triplicate'
combo_eob.Drug2 = combo_eob.Drug2.str[0]+'_'+combo_eob.Hours.astype(int).astype(str)+'_triplicate'
combo_eob['combo']=combo_eob.Drug1.str[0]+combo_eob.Drug2.str[0]+'_'+combo_eob.Hours.astype(int).astype(str)+'_triplicate'
combo_viab = combo_eob[['combo','Observed Viability','Observed Inhibition StdDev']]
combo_viab.set_index('combo',inplace=True)
combo_viab.columns = ['Viability','Std']

combo_eob = combo_eob[['combo',
    'Hours','Drug1','Drug 1 Viability','Drug2','Drug 2 Viability','Observed Viability']]
eobs = pd.concat([ss_eob,combo_eob]).drop('Hours',axis=1)

#artificially add a data point where T at 24 hrs is replaced with T_20:
row = eobs.loc[38]
row['Drug1']=eobs.loc[5]['Drug2']
row['Drug 1 Viability']=eobs.loc[5]['Drug 2 Viability']
#row.combo = 'TM_24_test'
eobs.loc[65]=row
row['Drug2']= eobs.loc[2].Drug2
row['Drug 2 Viability']=eobs.loc[2]['Drug 2 Viability']
#row.combo = 'TM_24_test2'
eobs.loc[66]=row

eobs['EOB']= 100*(eobs['Drug 1 Viability']/100*eobs['Drug 2 Viability']/100 - eobs['Observed Viability']/100)
eobs['pearson_r']=eobs.combo.apply(lambda x: new_corr_dict[x]) #choose type of correlation here #new_corr_dict[x] <---this one is DEGs only
eobs['pearson_p']=eobs.combo.apply(lambda x: new_p_corr_dict[x])

def lensyn(row):
    if len(row.Drug1)<6:
        if len(sets[row.combo])>0:
            return pd.Series([
                    len((sets[row.combo].difference(sets[row.Drug1])).difference(sets[row.Drug2])),
                    len((sets[row.combo].difference(sets[row.Drug1])).difference(sets[row.Drug2]))/len(sets[row.combo])]+
    [len(((sets[row.combo].difference(sets[row.Drug1])).difference(sets[row.Drug2])).intersection(
            set(apop.loc[i].dropna().unique()))) for i in apop.index]+
    [-np.log(Fishers(sets[row.combo+'_down'] - sets[row.Drug1+'_down'] - sets[row.Drug2+'_down'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))]+
    [-np.log(Fishers(sets[row.combo+'_up'] - sets[row.Drug1+'_up'] - sets[row.Drug2+'_up'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))])
        else:
            return pd.Series([
                    len((sets[row.combo].difference(sets[row.Drug1])).difference(sets[row.Drug2])),0,0,
                    -np.log(Fishers(sets[row.combo+'_down'] - sets[row.Drug1+'_down'] - sets[row.Drug2+'_down'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2)),
                    -np.log(Fishers(sets[row.combo+'_up'] - sets[row.Drug1+'_up'] - sets[row.Drug2+'_up'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))])
    else:
        if len(sets[row.combo])>0:
            return pd.Series([
                    len((sets[row.combo].difference(sets[row.combo[0]+'_'+row.combo[3:]])).difference(
                sets[row.combo[1]+'_'+row.combo[3:]])),
                len((sets[row.combo].difference(sets[row.combo[0]+'_'+row.combo[3:]])).difference(
                sets[row.combo[1]+'_'+row.combo[3:]]))/len(sets[row.combo])]+
    [len(((sets[row.combo].difference(sets[row.combo[0]+'_'+row.combo[3:]])).difference(
                sets[row.combo[1]+'_'+row.combo[3:]])).intersection(
            set(apop.loc[i].dropna().unique()))) for i in apop.index]+
    [-np.log(Fishers(sets[row.combo+'_down'] - sets[row.combo[0]+'_'+row.combo[3:]+'_down'] - sets[row.combo[1]+'_'+row.combo[3:]+'_down'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))]+
    [-np.log(Fishers(sets[row.combo+'_up'] - sets[row.combo[0]+'_'+row.combo[3:]+'_up'] - sets[row.combo[1]+'_'+row.combo[3:]+'_up'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))])
        else:
            return pd.Series([
                    len((sets[row.combo].difference(sets[row.combo[0]+'_'+row.combo[3:]])).difference(
                sets[row.combo[1]+'_'+row.combo[3:]])),0,0,
                -np.log(Fishers(sets[row.combo+'_down'] - sets[row.combo[0]+'_'+row.combo[3:]+'_down'] - sets[row.combo[1]+'_'+row.combo[3:]+'_down'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2)),
                -np.log(Fishers(sets[row.combo+'_up'] - sets[row.combo[0]+'_'+row.combo[3:]+'_up'] - sets[row.combo[1]+'_'+row.combo[3:]+'_up'],
                   set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))])
            
eobs[['nsyn','psyn']+apop.index.tolist()+['energy-enrichment-down',
     'energy-enrichment-up']]=eobs.apply(lensyn,axis=1)

eobs.combo=eobs.combo.str.split('_triplicate').str[0]

eobs['nsyn_log'] = np.log(eobs.nsyn)
print(eobs[['combo','EOB','pearson_r','nsyn_log']])

plt.bar(x=range(eobs.shape[0]),height=eobs['Observed Viability'],tick_label=eobs.combo,
        color=['b','g','r','b','g','r']+(['b']*5+['r'])*3)
plt.xticks(rotation='vertical')
plt.ylabel('Viability')
plt.title('Viability of Combination Treatments')
plt.savefig('viabilities.png')

small = eobs[eobs.combo.isin(['M_15','T_25','TM_24'])]

fig,ax = plt.subplots()
ax.scatter(eobs.pearson_r,eobs.EOB)
ax.scatter(small.pearson_r,small.EOB,color='r')
for i, row in eobs.iterrows():
    ax.annotate(row.combo,(row.pearson_r,row.EOB))
plt.xlabel('Pearson r for DEGs')
plt.ylabel('Excess Over Bliss')
plt.title('EOB vs Correlation')
plt.savefig('eob_pearson_degs.png')

fig,ax = plt.subplots()
ax.scatter(eobs.nsyn,eobs.EOB)
ax.scatter(small.nsyn,small.EOB,color='r')
for i, row in eobs.iterrows():
    ax.annotate(row.combo,(row.nsyn,row.EOB))
plt.xlabel('Synergistic Genes')
plt.ylabel('Excess Over Bliss')
plt.title('EOB vs Synergistic Genes')
plt.savefig('syngenes_eob.png')

fig,ax = plt.subplots()
ax.scatter(eobs.pearson_r,eobs.nsyn)
ax.scatter(small.pearson_r,small.nsyn,color='r')
for i, row in eobs.iterrows():
    ax.annotate(row.combo,(row.pearson_r,row.nsyn))
plt.xlabel('Pearson r for DEGs')
plt.ylabel('Synergistic Genes')
plt.title('Synergistic Genes vs Correlation')
plt.savefig('syngenes_pearson_degs.png')

fig,ax = plt.subplots()
ax.scatter(eobs['Observed Viability'],eobs.EOB)
ax.scatter(small['Observed Viability'],small.EOB,color='r')
for i, row in eobs.iterrows():
    ax.annotate(row.combo,(row['Observed Viability'],row.EOB))
plt.xlabel('Observed Viability')
plt.ylabel('Excess Over Bliss')
plt.title('EOB vs Viability')
plt.savefig('eob_viability.png')

fig,ax = plt.subplots()
ax.scatter(eobs['Observed Viability'],eobs.nsyn)
ax.scatter(small['Observed Viability'],small.nsyn,color='r')
for i, row in eobs.iterrows():
    ax.annotate(row.combo,(row['Observed Viability'],row.nsyn))
plt.xlabel('Observed Viability')
plt.ylabel('Synergistic Genes')
plt.title('Synergistic Genes vs Viability')
plt.savefig('syngenes_viabilty.png')

eobs = eobs[eobs.Drug1.str.len()<6]
eobs['Hours']=24
eobs['cell_line']='MCF7'
eobs.set_index('combo',inplace=True)

eobs[['Hours','Drug1','Drug 1 Viability','Drug2','Drug 2 Viability',
      'Observed Viability','EOB','psyn','nsyn']+apop.index.tolist()+['pearson_r',
                                              'pearson_p','cell_line',
                                              'energy-enrichment-down',
                                              'energy-enrichment-up']].to_csv(
      'self_synergy_combination_degs0.1.csv')

##viability plot
viab = pd.concat([ss_viab,combo_viab])
viab.index=viab.index.str.strip('_triplicate')
b = viab.loc[['M_2.5', 'M_5', 'M_10', 'M_15', 'T_5', 'T_10', 'T_20', 'T_25','TM_24']]
alphas = [float(tx.split('_')[1])/25. for tx in b.index[:-1]]+[1.]
bcolors = ['r']*4+['b']*4+['magenta']

plt.close('all')
plt.figure(figsize=[8,5])
plt.bar(range(b.shape[0]),height=b.Viability,yerr=b.Std,
        color=[make_rgb_transparent(bcolors[i],alphas[i]) for i in range(len(bcolors))])
plt.xticks(range(b.shape[0]),[
        tx.replace('_',' ')+ r'$\mu$'+'M' for tx in b.index if not 'TM' in tx]+['TM'],**font)
plt.ylabel('Viability',**font)
plt.title('Viability of Treated MCF7 Cells at 24 hours',**font)
plt.show()
plt.savefig('self_synergy_viab_forpaper3.pdf')