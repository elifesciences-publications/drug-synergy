#synergistic/diff expressed gene sets, for each combo
#jul 8 2015

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn as venn
import csv
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
font = {'fontname':'Arial','size':19}
#matplotlib.rc('font', **font)

##set p value cutoff
pval = 1E-15

txs = ['M','T','W','TM','TW','MW','TMW']
times = ['0','3','6','9','12','24']

allRNAseq = set(pd.read_csv('LIMMA_DATA/DE_limma_M_0.csv',header=0,
                            index_col='Unnamed: 0').index)

denames = []
for tx in txs:
    for time in times:
        denames.append('LIMMA_DATA/DE_limma_'+tx+'_'+time+'.csv')

dfs = {}
sets = {}
for name in denames:
    df = pd.read_csv(name,header=0,index_col='Unnamed: 0')
    dfs[name.split('limma_')[1].strip('.csv')] = df
    p = []
    sets[name.split('limma_')[1].strip('.csv')] = set(df[df['adj.P.Val'] <pval].index)
    sets[name.split('limma_')[1].strip('.csv')+'_up'] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']>0)].index)
    sets[name.split('limma_')[1].strip('.csv')+'_down'] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']<0)].index)

tm = []
tw = []
mw = []

T1 = []
M1 = []
TM = []
TMnT = []
TMnM = []
TnM = []
TnMnTM = []

W2= []
T2 = []
TW = []
TWnW = []
TWnT = []
TnW = []
TWnTnW = []

W3= []
M3 = []
MW = []
MWnM = []
MWnW = []
MnW = []
MWnMnW = []

#viabilities
viab = pd.read_excel('../../../drug synergy project/viability/LNCap Combination Viability Results Nov 2015.xlsx')
viab = viab[(viab['Viability Assay']=='InCell')&(viab.Hours<48)][[
    'Hours','Drug1','Drug 1 Viability','Drug2','Drug 2 Viability','Observed Viability']]
viab['EOB']= 100*(viab['Drug 1 Viability']/100*viab['Drug 2 Viability']/100 - viab['Observed Viability']/100)
viab['combo']=viab.Drug1.str[0]+viab.Drug2.str[0]+'_'+viab.Hours.astype(int).astype(str)
viab.set_index('combo',inplace=True)
viab['psyn'] = np.nan
viab['nsyn'] = np.nan

#apoptosis genes
apop = pd.read_csv('../../libraries_from_enrichr/GO_Biological_Process.txt',sep='\t',
                 names = list(range(400)),index_col=0,dtype=str)
apop = apop[apop.index.isin(['intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress (GO:0070059)', \
            'intrinsic apoptotic signaling pathway (GO:0097193)', \
            'apoptotic signaling pathway (GO:0097190)',\
            'generation of precursor metabolites and energy (GO:0006091)'])].dropna(axis=1,how='all')

for i in apop.index:
    viab[i]=np.nan
viab['energy-enrichment'] = np.nan

#get differential expression venn diagrams and a line graph of percent synergistic
for time in times: 
    ######all diff expressed genes
    plt.figure()
    v = venn.venn3([sets['T_'+time],sets['M_'+time],sets['TM_'+time]],
    set_labels = ('Tamoxifen','Mefloquine','TM'))
    for i in ['001','100','010','110','011','101','111']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_text('')
    for text in v.set_labels:
        text.set_fontsize(46)
        text.set_family('Arial')
#    plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TM_'+time+'forpaper.pdf')
    viab.loc['TM_'+time,'energy-enrichment-down'] = -np.log(Fishers(sets['TM_'+time+'_down'] - sets['T_'+time+'_down'] - sets['M_'+time+'_down'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    viab.loc['TM_'+time,'energy-enrichment-up'] = -np.log(Fishers(sets['TM_'+time+'_up'] - sets['T_'+time+'_up'] - sets['M_'+time+'_up'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    if len(sets['TM_'+time]) >0:
        tm.append(len((sets['TM_'+time].difference(sets['T_'+time])).difference(
            sets['M_'+time]))/len(sets['TM_'+time]))
        viab.loc['TM_'+time,'psyn']=tm[-1]
        viab.loc['TM_'+time,'nsyn']=len((sets['TM_'+time].difference(sets['T_'+time])).difference(
            sets['M_'+time]))
        for i in apop.index:
            viab.loc['TM_'+time,i]=len(((sets['TM_'+time].difference(sets['T_'+time])).difference(
            sets['M_'+time])).intersection(set(apop.loc[i].dropna())))
    else:
        tm.append(0)
        viab.loc['TM_'+time,'psyn']=0
        viab.loc['TM_'+time,'nsyn']=0
        for i in apop.index:
            viab.loc['TM_'+time,i]=0
    T1.append(len(((sets['T_'+time]).difference(sets['M_'+time])).difference(
        sets['TM_'+time])))
    M1.append(len(((sets['M_'+time]).difference(sets['T_'+time])).difference(
        sets['TM_'+time])))
    TnM.append(len(((sets['T_'+time]).intersection(sets['M_'+time])).difference(
        sets['TM_'+time])))
    TnMnTM.append(len(((sets['T_'+time]).intersection(sets['M_'+time])).intersection(
        sets['TM_'+time])))
    TMnT.append(len(((sets['T_'+time]).intersection(sets['TM_'+time])).difference(
        sets['M_'+time])))
    TMnM.append(len(((sets['M_'+time]).intersection(sets['TM_'+time])).difference(
        sets['T_'+time])))
    TM.append(len((sets['TM_'+time].difference(sets['T_'+time])).difference(
        sets['M_'+time])))
    #
    plt.figure()
    v = venn.venn3([sets['T_'+time],sets['W_'+time],sets['TW_'+time]],
    set_labels = ('Tamoxifen','Withaferin','TW'))
    for i in ['001','100','010','110','011','101','111']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_text('')
    for text in v.set_labels:
        text.set_fontsize(46)
        text.set_family('Arial')
#    plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TW_'+time+'forpaper.pdf')
    viab.loc['TW_'+time,'energy-enrichment-down'] = -np.log(Fishers(sets['TW_'+time+'_down'] - sets['T_'+time+'_down'] - sets['W_'+time+'_down'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    viab.loc['TW_'+time,'energy-enrichment-up'] = -np.log(Fishers(sets['TW_'+time+'_up'] - sets['T_'+time+'_up'] - sets['W_'+time+'_up'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    if len(sets['TW_'+time]) >0:
        tw.append(len((sets['TW_'+time].difference(sets['T_'+time])).difference(
            sets['W_'+time]))/len(sets['TW_'+time]))
        viab.loc['TW_'+time,'psyn']=tw[-1]
        viab.loc['TW_'+time,'nsyn']=len((sets['TW_'+time].difference(sets['T_'+time])).difference(
            sets['W_'+time]))
        for i in apop.index:
            viab.loc['TW_'+time,i]=len(((sets['TW_'+time].difference(sets['T_'+time])).difference(
            sets['W_'+time])).intersection(set(apop.loc[i].dropna())))
    else:
        tw.append(0)
        viab.loc['TW_'+time,'psyn']=0
        viab.loc['TW_'+time,'nsyn']=0
        for i in apop.index:
            viab.loc['TW_'+time,i]=0
    W2.append(len(((sets['W_'+time]).difference(sets['T_'+time])).difference(
        sets['TW_'+time])))
    T2.append(len(((sets['T_'+time]).difference(sets['W_'+time])).difference(
        sets['TW_'+time])))
    TnW.append(len(((sets['T_'+time]).intersection(sets['W_'+time])).difference(
        sets['TW_'+time])))
    TWnTnW.append(len(((sets['T_'+time]).intersection(sets['W_'+time])).intersection(
        sets['TW_'+time])))
    TWnW.append(len(((sets['TW_'+time]).intersection(sets['W_'+time])).difference(
        sets['T_'+time])))
    TWnT.append(len(((sets['TW_'+time]).intersection(sets['T_'+time])).difference(
        sets['W_'+time])))
    TW.append(len((sets['TW_'+time].difference(sets['T_'+time])).difference(
        sets['W_'+time])))
    #
    plt.figure()
    v = venn.venn3([sets['M_'+time],sets['W_'+time],sets['MW_'+time]],
    set_labels = ('Mefloquine','Withaferin','MW'))
    for i in ['001','100','010','110','011','101','111']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_text('')
    for text in v.set_labels:
        text.set_fontsize(46)
        text.set_family('Arial')
#    plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/MW_'+time+'forpaper.pdf')
    viab.loc['MW_'+time,'energy-enrichment-down'] = -np.log(Fishers(sets['MW_'+time+'_down'] - sets['M_'+time+'_down'] - sets['W_'+time+'_down'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    viab.loc['MW_'+time,'energy-enrichment-up'] = -np.log(Fishers(sets['MW_'+time+'_up'] - sets['M_'+time+'_up'] - sets['W_'+time+'_up'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    if len(sets['MW_'+time]) >0:
        mw.append(len((sets['MW_'+time].difference(sets['M_'+time])).difference(
            sets['W_'+time]))/len(sets['MW_'+time]))
        viab.loc['MW_'+time,'psyn']=mw[-1]
        viab.loc['MW_'+time,'nsyn']=len((sets['MW_'+time].difference(sets['M_'+time])).difference(
            sets['W_'+time]))
        for i in apop.index:
            viab.loc['MW_'+time,i]=len(((sets['MW_'+time].difference(sets['M_'+time])).difference(
            sets['W_'+time])).intersection(set(apop.loc[i].dropna())))
    else:
        mw.append(0)
        viab.loc['MW_'+time,'psyn']=0
        viab.loc['MW_'+time,'nsyn']=0
        for i in apop.index:
            viab.loc['MW_'+time,i]=0
    W3.append(len(((sets['W_'+time]).difference(sets['M_'+time])).difference(
        sets['MW_'+time])))
    M3.append(len(((sets['M_'+time]).difference(sets['W_'+time])).difference(
        sets['MW_'+time])))
    MnW.append(len(((sets['W_'+time]).intersection(sets['M_'+time])).difference(
        sets['MW_'+time])))
    MWnMnW.append(len(((sets['W_'+time]).intersection(sets['M_'+time])).intersection(
        sets['MW_'+time])))
    MWnW.append(len(((sets['MW_'+time]).intersection(sets['W_'+time])).difference(
        sets['M_'+time])))
    MWnM.append(len(((sets['MW_'+time]).intersection(sets['M_'+time])).difference(
        sets['W_'+time])))
    MW.append(len((sets['MW_'+time].difference(sets['M_'+time])).difference(
        sets['W_'+time])))

    ##plot veen diagrams for triple comb
    plt.figure()
    v = venn.venn3([sets['TM_'+time],sets['MW_'+time],sets['TMW_'+time]],
    set_labels = ('TM','MW','TMW'))
##    for i in ['001','100','010','110','011','101','111']:
##          if v.get_label_by_id(i) != None:
##              v.get_label_by_id(i).set_text('')
##    for text in v.set_labels:
##        text.set_fontsize(46)
##        text.set_family('Arial')
    #plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TM_MW_TMW'+time+'.png')

    plt.figure()
    v = venn.venn3([sets['TM_'+time],sets['TW_'+time],sets['TMW_'+time]],
    set_labels = ('TM','TW','TMW'))
##    for i in ['001','100','010','110','011','101','111']:
##          if v.get_label_by_id(i) != None:
##              v.get_label_by_id(i).set_text('')
##    for text in v.set_labels:
##        text.set_fontsize(46)
##        text.set_family('Arial')
    #plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TM_TW_TMW'+time+'.png')

    plt.figure()
    v = venn.venn3([sets['TW_'+time],sets['MW_'+time],sets['TMW_'+time]],
    set_labels = ('TW','MW','TMW'))
##    for i in ['001','100','010','110','011','101','111']:
##          if v.get_label_by_id(i) != None:
##              v.get_label_by_id(i).set_text('')
##    for text in v.set_labels:
##        text.set_fontsize(46)
##        text.set_family('Arial')
    #plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TW_MW_TMW'+time+'.png')

#plot overlaps between synergistic genes
    plt.figure()
    v = venn.venn3(
        [(sets['TM_'+time].difference(sets['T_'+time])).difference(
            sets['M_'+time]),
        (sets['TW_'+time].difference(sets['T_'+time])).difference(
            sets['W_'+time]),
         (sets['MW_'+time].difference(sets['M_'+time])).difference(
             sets['W_'+time])],
        set_labels = ('TM','TW','MW'))
    pt = Fishers((sets['TM_'+time].difference(sets['T_'+time])).difference(
            sets['M_'+time]),
                 (sets['TW_'+time].difference(sets['T_'+time])).difference(
            sets['W_'+time]),
                 3,2)
    pm = Fishers((sets['TM_'+time].difference(sets['T_'+time])).difference(
            sets['M_'+time]),
                 (sets['MW_'+time].difference(sets['M_'+time])).difference(
             sets['W_'+time]),
                 3,2)
    pw = Fishers((sets['TW_'+time].difference(sets['T_'+time])).difference(
            sets['W_'+time]),
                 (sets['MW_'+time].difference(sets['M_'+time])).difference(
             sets['W_'+time]),
                 3,2)
    if v.get_label_by_id('110') != None:
        v.get_label_by_id('110').set_text(str(v.get_label_by_id('110')).split(',')[2].strip(")'")+'\n'+str(pt))
        v.get_label_by_id('101').set_text(str(v.get_label_by_id('101')).split(',')[2].strip(")'")+'\n'+str(pm))
        v.get_label_by_id('011').set_text(str(v.get_label_by_id('011')).split(',')[2].strip(")'")+'\n'+str(pw)) 
    plt.title('Synergistic genes in each combination at time '+time)
    plt.savefig('venn/synergistic_'+time+'.png')
    
    
    
#
#plt.figure()
#plt.plot(tm,'k',linewidth=4)
#plt.plot(tw, 'k--',linewidth=4)
#plt.plot(mw,'k:',linewidth=4)
#plt.xticks([0,1,2,3,4,5],times)
#plt.yticks([0,0.2,0.4,0.6,0.8,1.0],[0,20,40,60,80,100])
#plt.legend(['TM','TW','MW'],bbox_to_anchor=[1,0.6])
#plt.xlabel('Time')
#plt.ylabel('Percent of Differentially Expressed Genes')
#plt.title('Percent of DE Genes that are Synergistic')
#plt.savefig('percent_synergistic.png')
#

plt.figure()
plt.stackplot([0,3,6,9,12,24],T1,TnM,M1,TMnT,TnMnTM,TMnM,TM,linewidth=0,
              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
              #colors=('r',[.7,.35,.008],'g',[.7,.008,.7],[.008,.35,.7],[.4,.2,.4],'b'),alpha=0.4)

plt.xticks([0,3,6,9,12,24],times,**font)
plt.ylabel('Differentially Expressed Genes',**font)
plt.xlabel('Time (hours)',**font)
#plt.title('Differentially Expressed Genes in Combination TM',**font)
plt.xlim([0,24])
plt.ylim([0,6000])
plt.yticks(**font)
#plt.legend(['Tamoxifen','Tamoxifen-Mefloquine','Mefloquine','Tamoxifen-TM','Mefloquine-TM','Intersection','TM'])#,bbox_to_anchor=[0.5,1] 
plt.tight_layout()
plt.savefig('Genes_TM.pdf') #

plt.figure()
plt.stackplot([0,3,6,9,12,24],T2,TnW,W2,TWnT,TWnTnW,TWnW,TW,linewidth=0,
              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
plt.xticks([0,3,6,9,12,24],times,**font)
plt.ylabel('Differentially Expressed Genes',**font)
plt.xlabel('Time (hours)',**font)
#plt.title('Differentially Expressed Genes in Combination TW',**font)
plt.xlim([0,24])
plt.ylim([0,6000])
plt.yticks(**font)
#plt.legend(['Withaferin','Intersection','Withaferin-TW','TW'],bbox_to_anchor=[0.5,1])
plt.tight_layout()
plt.savefig('Genes_TW.pdf') #,dpi=300

plt.figure()
plt.stackplot([0,3,6,9,12,24],M3,MnW,W3,MWnM,MWnMnW,MWnW,MW,linewidth=0,
              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
plt.xticks([0,3,6,9,12,24],times,**font)
plt.ylabel('Differentially Expressed Genes',**font)
plt.xlabel('Time (hours)',**font)
#plt.title('Differentially Expressed Genes in Combination MW',**font)
plt.xlim([0,24])
plt.ylim([0,6000])
plt.yticks(**font)
#plt.legend(['Withaferin','Intersection','Withaferin-MW','MW'],bbox_to_anchor=[0.5,1])
plt.tight_layout()
plt.savefig('Genes_MW.pdf') #,dpi=300

#
###    
###    ######up genes
###    plt.figure()
###    v = venn.venn3([sets['diffexp_T_'+time+'_up'],sets['diffexp_M_'+time+'_up'],
###    sets['diffexp_TM_'+time+'_up']],set_labels = ('Tamoxifen','Mefloquine','TM'))
###    plt.title(time+' Upregulated')
###    plt.show()
###    plt.savefig('TM_'+time+'_up.png')
###    #
###    plt.figure()
###    v = venn.venn3([sets['diffexp_T_'+time+'_up'],sets['diffexp_W_'+time+'_up'],
###    sets['diffexp_TW_'+time+'_up']], set_labels = ('Tamoxifen','Witheraferin','TW'))
###    plt.title(time+' Upgregulated')
###    plt.show()
###    plt.savefig('TW_'+time+'_up.png')
###    #
###    plt.figure()
###    v = venn.venn3([sets['diffexp_M_'+time+'_up'],sets['diffexp_W_'+time+'_up'],
###    sets['diffexp_MW_'+time+'_up']],set_labels = ('Mefloquine','Witheraferin','MW'))
###    plt.title(time+' Upregulated')
###    plt.show()
###    plt.savefig('MW_'+time+'_up.png')
###    #
###    ######down genes
###    plt.figure()
###    v = venn.venn3([sets['diffexp_T_'+time+'_down'],sets['diffexp_M_'+time+'_down'],
###    sets['diffexp_TM_'+time+'_down']],set_labels = ('Tamoxifen','Mefloquine','TM'))
###    plt.title(time+' Downregulated')
###    plt.show()
###    plt.savefig('TM_'+time+'_down.png')
###    #
###    plt.figure()
###    v = venn.venn3([sets['diffexp_T_'+time+'_down'],sets['diffexp_W_'+time+'_down'],
###    sets['diffexp_TW_'+time+'_down']],set_labels = ('Tamoxifen','Witheraferin','TW'))
###    plt.title(time+' Downregulated')
###    plt.show()
###    plt.savefig('TW_'+time+'_down.png')
###    #
###    plt.figure()
###    v = venn.venn3([sets['diffexp_M_'+time+'_down'],sets['diffexp_W_'+time+'_down'],
###    sets['diffexp_MW_'+time+'_down']],set_labels = ('Mefloquine','Witheraferin','MW'))
###    plt.title(time+' Downregulated')
###    plt.show()
###    plt.savefig('MW_'+time+'_down.png')
#    
#
#get diff exp genes in dataframe for fischer's exact tests 
#print(len(allRNAseq))
diffls = []
for tx in txs:
    df = pd.DataFrame(data=None,columns = ['gene','time','direction'])
    #print(tx)
    for time in times:
        diffls+=list(sets[tx+'_'+time])
        #print(time)
        #print(len(sets[tx+'_'+time]-allRNAseq))
        for direc in ['up','down']:
            newdf = pd.DataFrame(data={'gene':list(sets[tx+'_'+time+'_'+direc]),
            'time':time,'direction':direc})
            df = pd.concat([df,newdf])
    #print(df.head())
    df.to_csv('diffexp_lncap_'+tx+'_alltimes_JELD.csv',index=False)
##
#
viab = viab[viab.Drug1!='Tamoxifen - Mefloquine']
plt.figure()
plt.scatter(viab.psyn,viab.EOB)
for i,row in viab.iterrows():
    plt.annotate(row.name,(row.psyn,row.EOB))
plt.xlabel('percent synergistic genes')
plt.ylabel('EOB')
plt.savefig('psyn vs eob lncap.png')

plt.figure()
plt.scatter(viab.nsyn,viab.EOB)
for i,row in viab.iterrows():
    plt.annotate(row.name,(row.nsyn,row.EOB))
plt.xlabel('number of synergistic genes')
plt.ylabel('EOB')
plt.savefig('nsyn vs eob lncap.png')

def corr(dfs,pval=1.1):
    corr_dict = {}
    p_corr_dict = {}
    for time in times:
        for tx in ['TM','TW','MW']:
            s = set(dfs[tx[0]+'_'+str(time)][dfs[tx[0]+'_'+str(time)]['adj.P.Val']<pval].index).union(
                set(dfs[tx[1]+'_'+str(time)][dfs[tx[1]+'_'+str(time)]['adj.P.Val']<pval].index))
            corr_dict[tx+'_'+str(time)],p_corr_dict[tx+'_'+str(time)] = stats.pearsonr(
                dfs[tx[0]+'_'+str(time)][dfs[tx[0]+'_'+str(time)].index.isin(s)].sort_index().logFC,
                dfs[tx[1]+'_'+str(time)][dfs[tx[1]+'_'+str(time)].index.isin(s)].sort_index().logFC)
    #combos = [['M_10','T_20','M10.T20']]
#    combos += [[pair[0],pair[1],
#                pair[0][:2]+str(int(float(pair[0].split('_')[1])+
#                                    float(pair[1].split('_')[1])))] for pair \
#               in [['M_2.5']*2,['M_5']*2,['M_5','M_10'],['T_5']*2,['T_10']*2,['T_5','T_20']]]
#    for combo in combos:
#        s = set(dfs[combo[0]][dfs[combo[0]]['adj.P.Val']<pval].index).union(
#            set(dfs[combo[1]][dfs[combo[1]]['adj.P.Val']<pval].index))
#        corr_dict[combo[2]],p_corr_dict[combo[2]] = stats.pearsonr(
#                dfs[combo[0]][dfs[combo[0]].index.isin(s)].sort_index().logFC,
#                dfs[combo[1]][dfs[combo[1]].index.isin(s)].sort_index().logFC)
    return (corr_dict,p_corr_dict)

corr_dict,p_corr_dict = corr(dfs,pval=0.1)
viab['combo'] = viab.index
viab['pearson_r'] = viab.combo.apply(lambda x: corr_dict[x])
viab['pearson_p'] = viab.combo.apply(lambda x: p_corr_dict[x])
viab.drop('combo',axis=1,inplace=True)

plt.figure()
plt.scatter(viab['Observed Viability'],viab.EOB)
for i, row in viab.iterrows():
    plt.annotate(row.name,(row['Observed Viability'],row.EOB))
plt.xlabel('Viability')
plt.ylabel('Excess Over Bliss')
plt.title('EOB vs Viability')
plt.savefig('eob_viability.png')

plt.figure()
plt.scatter(viab.pearson_r,viab.EOB)
for i, row in viab.iterrows():
    plt.annotate(row.name,(row.pearson_r,row.EOB))
plt.xlabel('Pearson r for DEGs')
plt.ylabel('Excess Over Bliss')
plt.title('EOB vs Correlation')
plt.savefig('eob_pearson_degs.png')

plt.figure()
plt.scatter(viab.pearson_r,viab.nsyn)
for i, row in viab.iterrows():
    plt.annotate(row.name,(row.pearson_r,row.nsyn))
plt.xlabel('Pearson r for DEGs')
plt.ylabel('Synergistic Genes')
plt.title('Synergistic Genes vs Correlation')
plt.savefig('nsyn_pearson_degs.png')
    
viab['cell_line'] = 'LNCaP'
viab.to_csv('LNCaP_combinations0.1.csv')

plt.close('all')