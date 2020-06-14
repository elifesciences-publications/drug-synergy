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
pval = 1E-18

txs = ['M','T','W','TM','TW','MW']
times = ['0','3','6','9','12','24']

allRNAseq = set(pd.read_csv('limma_t0_M_0.txt',header=0,
                            sep='\t',index_col='1').index)

denames = []
for tx in txs:
    for time in times:
        denames.append('limma_t0_'+tx+'_'+time+'.txt')

dfs = {}
sets = {}
for name in denames:
    df = pd.read_csv(name,header=0,sep='\t',index_col='1')
    dfs[name.strip('.txt')] = df
    p = []
    sets[name.strip('.txt')] = set(df[df['adj.P.Val'] <pval].index)
    sets[name.strip('.txt')+'_up'] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']>0)].index)
    sets[name.strip('.txt')+'_down'] = set(df[(df['adj.P.Val'] <pval)&(df['logFC']<0)].index)

for key in dfs.keys():
    print(key.split('t0_')[1])
    df = dfs[key]
    print(df.loc['ABCB1'][['logFC','t','adj.P.Val']])
    
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

#viability data
viab = pd.read_excel('../../../drug synergy project/viability/IBM MCF7 Combination Viability February 2015 v3.xlsx',
    sheet_name='Combinations')
viab = viab[(viab['Assay']=='InCell')&(viab.Hours<48)][[
    'Hours','Drug1','Drug 1 Viability','Drug2','Drug 2 Viability','Observed Viability']]
viab['EOB']= 100*(viab['Drug 1 Viability']/100*viab['Drug 2 Viability']/100 - viab['Observed Viability']/100)
viab['combo']=viab.Drug1.str[0]+viab.Drug2.str[0]+'_'+viab.Hours.astype(int).astype(str)
viab.set_index('combo',inplace=True)
viab['psyn'] = np.nan
viab['nsyn'] = np.nan

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
    v = venn.venn3([sets['limma_t0_T_'+time],sets['limma_t0_M_'+time],sets['limma_t0_TM_'+time]],
    set_labels = ('Tamoxifen','Mefloquine','TM'))
    for i in ['001','100','010','110','011','101','111']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_text('')
    for text in v.set_labels:
        text.set_fontsize(36)
        text.set_family('Arial')
    #plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TM_'+time+'forpaper.pdf')
    if len(sets['limma_t0_TM_'+time]) >0:
        tm.append(len((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_M_'+time]))/len(sets['limma_t0_TM_'+time]))
        viab.loc['TM_'+time,'psyn']=tm[-1]
        viab.loc['TM_'+time,'nsyn']=len((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_M_'+time]))
        for i in apop.index:
            viab.loc['TM_'+time,i]=len(((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_M_'+time])).intersection(set(apop.loc[i].dropna())))
        viab.loc['TM_'+time,'energy-enrichment-down'] = -np.log(Fishers(sets['limma_t0_TM_'+time+'_down'] - sets['limma_t0_T_'+time+'_down'] - sets['limma_t0_M_'+time+'_down'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
        viab.loc['TM_'+time,'energy-enrichment-up'] = -np.log(Fishers(sets['limma_t0_TM_'+time+'_up'] - sets['limma_t0_T_'+time+'_up'] - sets['limma_t0_M_'+time+'_up'],
            set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
            3,2))
    else:
        tm.append(0)
        viab.loc['TM_'+time,'psyn']=0
        viab.loc['TM_'+time,'nsyn']=0
        for i in apop.index:
            viab.loc['TM_'+time,i]=0
    T1.append(len(((sets['limma_t0_T_'+time]).difference(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_TM_'+time])))
    M1.append(len(((sets['limma_t0_M_'+time]).difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_TM_'+time])))
    TnM.append(len(((sets['limma_t0_T_'+time]).intersection(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_TM_'+time])))
    TnMnTM.append(len(((sets['limma_t0_T_'+time]).intersection(sets['limma_t0_M_'+time])).intersection(
        sets['limma_t0_TM_'+time])))
    TMnT.append(len(((sets['limma_t0_T_'+time]).intersection(sets['limma_t0_TM_'+time])).difference(
        sets['limma_t0_M_'+time])))
    TMnM.append(len(((sets['limma_t0_M_'+time]).intersection(sets['limma_t0_TM_'+time])).difference(
        sets['limma_t0_T_'+time])))
    TM.append(len((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_M_'+time])))
    #
    plt.figure()
    v = venn.venn3([sets['limma_t0_T_'+time],sets['limma_t0_W_'+time],sets['limma_t0_TW_'+time]],
    set_labels = ('Tamoxifen','Withaferin','TW'))
    for i in ['001','100','010','110','011','101','111']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_text('')
    for text in v.set_labels:
        text.set_fontsize(36)
        text.set_family('Arial')
    #plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/TW_'+time+'forpaper.pdf')
    tw.append(len((sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_W_'+time]))/len(sets['limma_t0_TW_'+time]))
    viab.loc['TW_'+time,'psyn']=tw[-1]
    viab.loc['TW_'+time,'nsyn']=len((sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_W_'+time]))
    for i in apop.index:
        viab.loc['TW_'+time,i]=len(((sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_W_'+time])).intersection(set(apop.loc[i].dropna())))
    viab.loc['TW_'+time,'energy-enrichment-down'] = -np.log(Fishers(sets['limma_t0_TW_'+time+'_down'] - sets['limma_t0_T_'+time+'_down'] - sets['limma_t0_W_'+time+'_down'],
        set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
        3,2))
    viab.loc['TW_'+time,'energy-enrichment-up'] = -np.log(Fishers(sets['limma_t0_TW_'+time+'_up'] - sets['limma_t0_T_'+time+'_up'] - sets['limma_t0_W_'+time+'_up'],
        set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
        3,2))
    W2.append(len(((sets['limma_t0_W_'+time]).difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_TW_'+time])))
    T2.append(len(((sets['limma_t0_T_'+time]).difference(sets['limma_t0_W_'+time])).difference(
        sets['limma_t0_TW_'+time])))
    TnW.append(len(((sets['limma_t0_T_'+time]).intersection(sets['limma_t0_W_'+time])).difference(
        sets['limma_t0_TW_'+time])))
    TWnTnW.append(len(((sets['limma_t0_T_'+time]).intersection(sets['limma_t0_W_'+time])).intersection(
        sets['limma_t0_TW_'+time])))
    TWnW.append(len(((sets['limma_t0_TW_'+time]).intersection(sets['limma_t0_W_'+time])).difference(
        sets['limma_t0_T_'+time])))
    TWnT.append(len(((sets['limma_t0_TW_'+time]).intersection(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_W_'+time])))
    TW.append(len((sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
        sets['limma_t0_W_'+time])))
    #
    plt.figure()
    v = venn.venn3([sets['limma_t0_M_'+time],sets['limma_t0_W_'+time],sets['limma_t0_MW_'+time]],
    set_labels = ('Mefloquine','Withaferin','MW'))
    for i in ['001','100','010','110','011','101','111']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_text('')
    for text in v.set_labels:
        text.set_fontsize(36)
        text.set_family('Arial')
    #plt.title('Differentially Expressed Genes at '+time+' Hours')
    #plt.show()
    plt.savefig('venn/MW_'+time+'forpaper.pdf')
    mw.append(len((sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_W_'+time]))/len(sets['limma_t0_MW_'+time]))
    viab.loc['MW_'+time,'psyn']=mw[-1]
    viab.loc['MW_'+time,'nsyn']=len((sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_W_'+time]))
    for i in apop.index:
        viab.loc['MW_'+time,i]=len(((sets['limma_t0_MW_'+time].difference(sets['limma_t0_W_'+time])).difference(
        sets['limma_t0_M_'+time])).intersection(set(apop.loc[i].dropna())))
    viab.loc['MW_'+time,'energy-enrichment-down'] = -np.log(Fishers(sets['limma_t0_MW_'+time+'_down'] - sets['limma_t0_W_'+time+'_down'] - sets['limma_t0_M_'+time+'_down'],
        set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
        3,2))
    viab.loc['MW_'+time,'energy-enrichment-up'] = -np.log(Fishers(sets['limma_t0_MW_'+time+'_up'] - sets['limma_t0_W_'+time+'_up'] - sets['limma_t0_M_'+time+'_up'],
        set(apop.loc['generation of precursor metabolites and energy (GO:0006091)'].dropna().unique()),
        3,2))
    W3.append(len(((sets['limma_t0_W_'+time]).difference(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_MW_'+time])))
    M3.append(len(((sets['limma_t0_M_'+time]).difference(sets['limma_t0_W_'+time])).difference(
        sets['limma_t0_MW_'+time])))
    MnW.append(len(((sets['limma_t0_W_'+time]).intersection(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_MW_'+time])))
    MWnMnW.append(len(((sets['limma_t0_W_'+time]).intersection(sets['limma_t0_M_'+time])).intersection(
        sets['limma_t0_MW_'+time])))
    MWnW.append(len(((sets['limma_t0_MW_'+time]).intersection(sets['limma_t0_W_'+time])).difference(
        sets['limma_t0_M_'+time])))
    MWnM.append(len(((sets['limma_t0_MW_'+time]).intersection(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_W_'+time])))
    MW.append(len((sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
        sets['limma_t0_W_'+time])))

#plot overlaps between synergistic genes
    plt.figure()
    v = venn.venn3(
        [(sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_M_'+time]),
        (sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_W_'+time]),
         (sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
             sets['limma_t0_W_'+time])],
        set_labels = ('TM','TW','MW'))
    pt = Fishers((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_M_'+time]),
                 (sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_W_'+time]),
                 3,2)
    pm = Fishers((sets['limma_t0_TM_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_M_'+time]),
                 (sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
             sets['limma_t0_W_'+time]),
                 3,2)
    pw = Fishers((sets['limma_t0_TW_'+time].difference(sets['limma_t0_T_'+time])).difference(
            sets['limma_t0_W_'+time]),
                 (sets['limma_t0_MW_'+time].difference(sets['limma_t0_M_'+time])).difference(
             sets['limma_t0_W_'+time]),
                 3,2)
    if v.get_label_by_id('110') != None:
        v.get_label_by_id('110').set_text(str(v.get_label_by_id('110')).split(',')[2].strip(")'")+'\n'+str(pt))
        v.get_label_by_id('101').set_text(str(v.get_label_by_id('101')).split(',')[2].strip(")'")+'\n'+str(pm))
        v.get_label_by_id('011').set_text(str(v.get_label_by_id('011')).split(',')[2].strip(")'")+'\n'+str(pw)) 
    plt.title('Synergistic genes in each combination at time '+time)
    plt.savefig('venn/synergistic_'+time+'.png')
    
    
#    
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
#
#plt.figure()
#plt.stackplot([0,3,6,9,12,24],T1,TnM,M1,TMnT,TnMnTM,TMnM,TM,linewidth=0,
#              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
#              #colors=('r',[.7,.35,.008],'g',[.7,.008,.7],[.008,.35,.7],[.4,.2,.4],'b'),alpha=0.4)
#
#plt.xticks([0,3,6,9,12,24],times,**font)
#plt.ylabel('Differentially Expressed Genes',**font)
#plt.xlabel('Time (hours)',**font)
##plt.title('Differentially Expressed Genes in Combination TM',**font)
#plt.xlim([0,24])
#plt.ylim([0,7000])
#plt.yticks(**font)
##plt.legend(['Tamoxifen','Tamoxifen-Mefloquine','Mefloquine','Tamoxifen-TM','Mefloquine-TM','Intersection','TM'])#,bbox_to_anchor=[0.5,1] 
#plt.savefig('Genes_TM.pdf') #,dpi=300
#
#plt.figure()
#plt.stackplot([0,3,6,9,12,24],T2,TnW,W2,TWnT,TWnTnW,TWnW,TW,linewidth=0,
#              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
#plt.xticks([0,3,6,9,12,24],times,**font)
#plt.ylabel('Differentially Expressed Genes',**font)
#plt.xlabel('Time (hours)',**font)
##plt.title('Differentially Expressed Genes in Combination TW',**font)
#plt.xlim([0,24])
#plt.ylim([0,7000])
#plt.yticks(**font)
##plt.legend(['Withaferin','Intersection','Withaferin-TW','TW'],bbox_to_anchor=[0.5,1])
#plt.savefig('Genes_TW.pdf') #,dpi=300
#
#plt.figure()
#plt.stackplot([0,3,6,9,12,24],M3,MnW,W3,MWnM,MWnMnW,MWnW,MW,linewidth=0,
#              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
#plt.xticks([0,3,6,9,12,24],times,**font)
#plt.ylabel('Differentially Expressed Genes',**font)
#plt.xlabel('Time (hours)',**font)
##plt.title('Differentially Expressed Genes in Combination MW',**font)
#plt.xlim([0,24])
#plt.ylim([0,7000])
#plt.yticks(**font)
##plt.legend(['Withaferin','Intersection','Withaferin-MW','MW'],bbox_to_anchor=[0.5,1])
#plt.savefig('Genes_MW.pdf') #,dpi=300
#

##    
##    ######up genes
##    plt.figure()
##    v = venn.venn3([sets['diffexp_T_'+time+'_up'],sets['diffexp_M_'+time+'_up'],
##    sets['diffexp_TM_'+time+'_up']],set_labels = ('Tamoxifen','Mefloquine','TM'))
##    plt.title(time+' Upregulated')
##    plt.show()
##    plt.savefig('TM_'+time+'_up.png')
##    #
##    plt.figure()
##    v = venn.venn3([sets['diffexp_T_'+time+'_up'],sets['diffexp_W_'+time+'_up'],
##    sets['diffexp_TW_'+time+'_up']], set_labels = ('Tamoxifen','Witheraferin','TW'))
##    plt.title(time+' Upgregulated')
##    plt.show()
##    plt.savefig('TW_'+time+'_up.png')
##    #
##    plt.figure()
##    v = venn.venn3([sets['diffexp_M_'+time+'_up'],sets['diffexp_W_'+time+'_up'],
##    sets['diffexp_MW_'+time+'_up']],set_labels = ('Mefloquine','Witheraferin','MW'))
##    plt.title(time+' Upregulated')
##    plt.show()
##    plt.savefig('MW_'+time+'_up.png')
##    #
##    ######down genes
##    plt.figure()
##    v = venn.venn3([sets['diffexp_T_'+time+'_down'],sets['diffexp_M_'+time+'_down'],
##    sets['diffexp_TM_'+time+'_down']],set_labels = ('Tamoxifen','Mefloquine','TM'))
##    plt.title(time+' Downregulated')
##    plt.show()
##    plt.savefig('TM_'+time+'_down.png')
##    #
##    plt.figure()
##    v = venn.venn3([sets['diffexp_T_'+time+'_down'],sets['diffexp_W_'+time+'_down'],
##    sets['diffexp_TW_'+time+'_down']],set_labels = ('Tamoxifen','Witheraferin','TW'))
##    plt.title(time+' Downregulated')
##    plt.show()
##    plt.savefig('TW_'+time+'_down.png')
##    #
##    plt.figure()
##    v = venn.venn3([sets['diffexp_M_'+time+'_down'],sets['diffexp_W_'+time+'_down'],
##    sets['diffexp_MW_'+time+'_down']],set_labels = ('Mefloquine','Witheraferin','MW'))
##    plt.title(time+' Downregulated')
##    plt.show()
##    plt.savefig('MW_'+time+'_down.png')
    

#turns out both of the withaferin combos are basically the same gene sets as withaferin alone
#mefloquine and tamoxifen gene sets alone are very small. 
#TM combo gene sets are quite large. 
#Want to avoid running fischer's exact test with small numbers and giving it the same
#weight as large gene sets if possible. 
#in addition, from Eren's data, there are few synergistic genes in either withaferin condition
#and a lot in TM. About half of the differentially expressed genes are synergistic (both up and down)
#and many more are NOT DIFFERENTIALLY EXPRESSED. Would like to know what biological effects 
#these are all having.

#gene sets of interest:
#with open('synergistic_gene_set_library_august2015.gmt','wb') as libfile:
#    libfile = csv.writer(libfile,delimiter='\t')
#    for time in times[1:]:
#        #In TM:
#        #synergistic & not differentially expressed
#        libfile.writerow(['syngenes_notdiffexp_TM_'+time,''] + \
#        list(sets['syngenes_TM_'+time] - sets['diffexp_TM_'+time]))
#        #synergistic & upregulated
#        libfile.writerow(['syngenes_up_TM_'+time,''] + \
#        list(sets['syngenes_TM_'+time] & sets['diffexp_TM_'+time+'_up']))
#        #synergistic & downregulated
#        libfile.writerow(['syngenes_down_TM_'+time,''] + \
#        list(sets['syngenes_TM_'+time] & sets['diffexp_TM_'+time+'_down']))
#        #upregulated but not synergistic
#        libfile.writerow(['upregulated_notsynerg_TM_'+time,''] + \
#        list(sets['diffexp_TM_'+time+'_up'] - sets['syngenes_TM_'+time]))
#        #downregulated but not synergistic
#        libfile.writerow(['downregulated_notsynerg_TM_'+time,''] + \
#        list(sets['diffexp_TM_'+time+'_down'] - sets['syngenes_TM_'+time]))
#        
#        #upregulated in witheraferin
#        libfile.writerow(['upregulated_W_'+time,''] + list(sets['diffexp_W_'+time+'_up']))
#        #downregulated in withaferin
#        libfile.writerow(['downregulated_W_'+time,''] + list(sets['diffexp_W_'+time+'_down']))
#        
#        #synergistic in TW (nearly all of these are not differentially expressed)
#        libfile.writerow(['syngenes_TW_'+time,'']+list(sets['syngenes_TW_'+time]))
##    
###get diff exp genes in dataframe for fischer's exact tests 
##for tx in txs:
##    df = pd.DataFrame(data=None,columns = ['gene','time','direction'])
##    print(tx)
##    for time in times:
##        for direc in ['up','down']:
##            newdf = pd.DataFrame(data={'gene':list(sets['limma_t0_'+ tx+'_'+time+'_'+direc]),
##            'time':time,'direction':direc})
##            df = pd.concat([df,newdf])
##    print(df.head())
##    df.to_csv('diffexp_'+tx+'_alltimes_JELD.csv',index=False)
####

plt.figure()
plt.scatter(viab.psyn,viab.EOB)
plt.xlabel('percent synergistic genes')
plt.ylabel('EOB')
plt.savefig('psyn vs eob mcf7.png')

def get_corr(row):
    pd.Series(stats.spearmanr(dfs['limma_t0_'+row.name.split('_')[0][0]+
                               '_'+row.name.split('_')[1]].logFC,
                     dfs['limma_t0_'+row.name.split('_')[0][1]+
                               '_'+row.name.split('_')[1]].logFC))
    row.name.split('_')[0][1]

def corr(dfs,pval=1.1):
    corr_dict = {}
    p_corr_dict = {}
    for time in times:
        for tx in ['TM','TW','MW']:
            s = set(dfs['limma_t0_'+tx[0]+'_'+str(time)][dfs['limma_t0_'+tx[0]+'_'+str(time)]['adj.P.Val']<pval].index).union(
                set(dfs['limma_t0_'+tx[1]+'_'+str(time)][dfs['limma_t0_'+tx[1]+'_'+str(time)]['adj.P.Val']<pval].index))
            corr_dict[tx+'_'+str(time)],p_corr_dict[tx+'_'+str(time)] = stats.pearsonr(
                dfs['limma_t0_'+tx[0]+'_'+str(time)][dfs['limma_t0_'+tx[0]+'_'+str(time)].index.isin(s)].sort_index().logFC,
                dfs['limma_t0_'+tx[1]+'_'+str(time)][dfs['limma_t0_'+tx[1]+'_'+str(time)].index.isin(s)].sort_index().logFC)
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

viab['cell_line'] = 'MCF7'
viab.to_csv('MCF7_combinations0.1.csv')