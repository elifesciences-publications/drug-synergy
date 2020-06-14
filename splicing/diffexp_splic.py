#diff expressed/spliced gene sets, for each combo
#may 25 2016

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn as venn
import csv
from matplotlib.colors import ColorConverter
import numpy as np


##set fonts
font = {'fontname':'Arial','size':19} #use size 14 for figure 4

#access color dict for color mixing
c = ColorConverter.colors

###colors
##y = [i/255. for i in [249, 255, 0]]#[249, 255, 60]
##b = [i/255. for i in [0, 0, 255]]#[86, 112, 255]

#using matplotlib-venn colormixing similar to venn.__venn3.compute_colors
def mix(col1,col2,col3=None):
    #input a color as string
    if col3==None:
        return venn._common.mix_colors(np.array(c[col1]),np.array(c[col2]))
    else:
        return venn._common.mix_colors(np.array(c[col1]),np.array(c[col2]),np.array(c[col3]))


##set p value cutoff
pval = 1E-18

txs = ['M','T','W','TM','TW','MW']
times = ['0','3','6','9','12','24']

sets = {}
denames = []
for tx in txs:
    for time in times:
        df = pd.read_csv('../gene_expression/diff_limma_december/limma_t0_'+tx+'_'+time+'.txt',header=0,sep='\t',index_col='1')
        sets[tx+'_'+time+'_diff'] = set(df[df['adj.P.Val'] <pval].index)

txd = {'M':'MFL_singlet','T':'TMX_singlet','W':'WFA_singlet','TM':'TMX_MFL','TW':'TMX_WFA','MW':'MFL_WFA'}

for tx in txs:
    for time in times:
        try: pd.read_csv('../Differential Splicing/'+
            txd[tx]+'_DMSO_singlet_'+time+'_dynsyn_exon_O_fcounts_O_flatten_exon_name.mtrx_genes_simes_test_rare_splicing_FDR_10',
            sep='\t')
        except: sets[tx+'_'+time+'_splic'] = set()
        else: 
            df = pd.read_csv('../Differential Splicing/'+
                txd[tx]+'_DMSO_singlet_'+time+'_dynsyn_exon_O_fcounts_O_flatten_exon_name.mtrx_genes_simes_test_rare_splicing_FDR_10',
                sep='\t')
            sets[tx+'_'+time+'_splic']  = set(df[df.FDR<0.05].V12)

for combo in ['MW','TW','TM']:           
    for time in times:
        sets[combo+'_'+time+'_exp_syn'] = (sets[combo+'_'+time+'_diff'].difference(sets[combo[0]+'_'+time+'_diff'])).difference(
            sets[combo[1]+'_'+time+'_diff'])
        sets[combo+'_'+time+'_splic_syn'] = (sets[combo+'_'+time+'_splic'].difference(sets[combo[0]+'_'+time+'_splic'])).difference(
            sets[combo[1]+'_'+time+'_splic'])

#diffexp vs splicing stackplots
for tx in txs:
    both= []
    splic = []
    diffexp = []
    
    for time in times: 
        diffexp.append(len((sets[tx+'_'+time+'_diff']).difference(sets[tx+'_'+time+'_splic'])))
        splic.append(len((sets[tx+'_'+time+'_splic']).difference(sets[tx+'_'+time+'_diff'])))
        both.append(len((sets[tx+'_'+time+'_splic']).intersection(sets[tx+'_'+time+'_diff'])))

    plt.figure()
    plt.stackplot([0,1,2,3,4,5],diffexp,both,splic,linewidth = 0.1,colors=('b','g','y'),alpha=0.4)
    plt.xticks([0,1,2,3,4,5],times,**font)
    plt.yticks(**font)
    plt.xlim([0,24])
    plt.ylabel('Genes',**font)
    plt.xlabel('Time (hours)',**font)
    plt.title('Differentially Expressed or Spliced Genes in '+tx,**font)
    plt.rc('font',family='Arial')
    #plt.legend(['Differentially Expressed','Both','Differentially Spliced'],bbox_to_anchor=[0.37,0.99],
               #fontsize='smaller')#
    #plt.savefig('Genes_'+tx+'.png',dpi=300)

    #venn diagram at time 3 - match figure 1
    plt.figure()
    v = venn.venn2([sets[tx+'_3_diff'],sets[tx+'_3_splic']],set_labels=('',''),
               set_colors = ('b','y'))
    v.get_label_by_id('10').set_text('Differentially\nExpressed')
    v.get_label_by_id('01').set_text('Differentially\nSpliced')
    if v.get_label_by_id('11') != None:
        v.get_label_by_id('11').set_text('')
    for i in ['01','10','11']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_fontsize(30) #30 for figure 4
              v.get_label_by_id(i).set_family('Arial')
    if v.get_patch_by_id('11') != None:
        v.get_patch_by_id('11').set_color('g')
        v.get_patch_by_id('11').set_alpha(0.3)
        v.get_patch_by_id('11').set_edgecolor('none')
##    for text in v.set_labels:
##        text.set_fontsize(56)
##        text.set_family('Arial')
    #plt.savefig('Venn_differential_'+tx+'.png')


for combo in ['MW','TW','TM']:
    tx = combo
    both= []
    splic = []
    diffexp = []
    
    for time in times: 
        diffexp.append(len((sets[tx+'_'+time+'_exp_syn']).difference(sets[tx+'_'+time+'_splic_syn'])))
        splic.append(len((sets[tx+'_'+time+'_splic_syn']).difference(sets[tx+'_'+time+'_exp_syn'])))
        both.append(len((sets[tx+'_'+time+'_splic_syn']).intersection(sets[tx+'_'+time+'_exp_syn'])))

    plt.figure()
    plt.stackplot([0,3,6,9,12,24],diffexp,both,splic,linewidth = 0.1,colors=('b','g','y'),alpha=0.4)
    plt.xticks([0,3,6,9,12,24],times,**font)
    plt.yticks(**font)
    plt.ylim([0,6000])
    plt.xlim([0,24])
    plt.ylabel('Genes',**font)
    plt.xlabel('Time (hours)',**font)
    plt.title('Synergistically Expressed or Spliced Genes in '+tx,**font)
    plt.rc('font',family='Arial')
    #plt.legend(['Synergistically Expressed','Both','Synergistically Spliced'],bbox_to_anchor=[0.38,0.99],
               #fontsize='smaller')#
    plt.savefig('Genes_synergistic_'+tx+'.pdf',dpi=300)
    
    #venn diagram at time 3 - match figure 1
    plt.figure()
    v = venn.venn2([sets[tx+'_3_exp_syn'],sets[tx+'_3_splic_syn']],set_labels=('',''),
                   set_colors = ('b','y'))
    v.get_label_by_id('10').set_text('Synergistically     \nExpressed')
    v.get_label_by_id('01').set_text('      Synergistically\nSpliced')
    if v.get_label_by_id('11') != None:
        v.get_label_by_id('11').set_text('')
    for i in ['01','10','11']:
          if v.get_label_by_id(i) != None:
              v.get_label_by_id(i).set_fontsize(38) #30 for figure 4
              v.get_label_by_id(i).set_family('Arial')
    if v.get_patch_by_id('11') != None:
        v.get_patch_by_id('11').set_color('g')
        v.get_patch_by_id('11').set_alpha(0.3)
        v.get_patch_by_id('11').set_edgecolor('none')
##    for text in v.set_labels:
##        text.set_fontsize(56)
##        text.set_family('Arial')
    plt.savefig('Venn_synergistic_'+tx+'.pdf')
