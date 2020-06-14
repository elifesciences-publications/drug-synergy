#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 15:23:19 2019

@author: jenniferlongdiaz
heatmaps for paper
"""

import pandas as pd
import seaborn as sns
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

sns.reset_defaults()

def make_rgb_transparent(color, alpha):
    rgb = colors.colorConverter.to_rgb(color)
    return [alpha * c1 + (1 - alpha) * c2
            for (c1, c2) in zip(rgb, (1,1,1))]
    
sns.set(font='Arial',color_codes=False)

scpm = pd.read_csv('../self_synergy_data/ss_log_cpm.csv',index_col=0,header=0)
scpm = 2**scpm
smeans = pd.DataFrame()
for tx in list(set([col[:-1] for col in scpm.columns])):
    smeans[tx] = scpm[[col for col in scpm.columns if tx in col]].mean(axis=1)
smeans.columns = [col.strip('_') for col in smeans.columns]

lcpm = pd.read_csv('../LNCaP/lncap_expression.csv',index_col=0,header=0)
lcpm = 2**lcpm
lmeans = pd.DataFrame()
for tx in list(set([col[:-1] for col in lcpm.columns])):
    lmeans[tx] = lcpm[[col for col in lcpm.columns if tx in col]].mean(axis=1)
lmeans = lmeans[[col for col in lmeans.columns if not 'TMW' in col]]

mcpm = pd.read_csv('../log_cpm_viab.csv',index_col=0,header=0,skiprows=[1])
mcpm = 2**mcpm
mmeans = pd.DataFrame()
for tx in list(set([col[:-1] for col in mcpm.columns])):
    mmeans[tx] = mcpm[[col for col in mcpm.columns if tx in col]].mean(axis=1)

#diffexp lists
dfs = dict()
ss_genes = []
ln_genes = []
mc_genes = []
for tx in smeans.columns:
    if 'DMSO' not in tx:
        df = pd.read_csv(
                '../self_synergy_data/diffexp_selfsynergy_'+tx+'_JELD.csv')
        dfs['ss_'+tx] = df
        ss_genes+=df.gene.tolist()
for tx in ['M','T','W','MW','TW','TM']:
    df = pd.read_csv('../LNCaP/diffexp_lncap_'+tx+'_alltimes_JELD.csv')
    dfs['ln_'+tx] = df
    ln_genes+=df.gene.tolist()
for tx in ['M','T','W','MW','TW','TM']:
    df = pd.read_csv('../diff_limma_december/diffexp_'+tx+'_alltimes_JELD.csv')
    dfs['mc_'+tx] = df
    mc_genes+=df.gene.tolist()

ss_genes = list(set(ss_genes))
ln_genes = list(set(ln_genes))
mc_genes = list(set(mc_genes))

##make some heatmaps
#ss
sstxs =['DMSO', 'M_2.5', 'M_5','M_10','M_15','T_5', 'T_10','T_20','T_25']
alphas = [1.]+[float(tx.split('_')[1])/25. for tx in sstxs[1:]]
bcolors = ['gray']+['r']*4+['b']*4
doses = [2.5,5,10,15,20,25]

lut = dict(zip(sstxs,[make_rgb_transparent(bcolors[i],alphas[i]) for i in range(len(bcolors))]))
#smeans = smeans.loc[ss_genes][sstxs]
smeans = smeans[sstxs]

cg = sns.clustermap(smeans,cmap='bwr',col_cluster=False,
               z_score=0,
               yticklabels=False,col_colors=smeans.columns.map(lut),
               cbar_kws={'label':'z-score(cpm)'},vmin=-2.5,vmax=2.5,
               figsize=(5,8))#)
cg.ax_row_dendrogram.set_visible(False)
# Use the dendrogram box to reposition the colour bar
dendro_box = cg.ax_row_dendrogram.get_position()
dendro_box.x0 = (dendro_box.x0 + 2*dendro_box.x1) / 3
dendro_box.x1 = (dendro_box.x0 + dendro_box.x1)/2

lbox = cg.ax_col_dendrogram.get_position()
lbox.x0 = dendro_box.x0*0.3
lbox.x1 = dendro_box.x1
lbox.y1 = dendro_box.y1*1.06
#dendro_box.y0 = dendro_box.y0*1.5
dendro_box.y1 = dendro_box.y1*0.6
lbox.y0 = dendro_box.y1
cg.cax.set_position(dendro_box)
# Move the ticks to the left (https://stackoverflow.com/a/36939552/1878788)
cg.cax.yaxis.set_ticks_position("left")
cg.cax.yaxis.set_label_position('left')
#color legend
cg.ax_col_dendrogram.set_position(lbox)
cg.ax_col_dendrogram.legend(
        [Line2D([0],[0],
                color='w',markerfacecolor=c,marker='s') for c in ['grey','r','b']]+\
        [Line2D([0],[0],
                color='w',markerfacecolor='k',marker='s',
                alpha = d/25.) for d in doses],
        ['DMSO','M','T']+[str(d)+ r'$\mu$'+'M' for d in doses],
        loc='upper left',frameon=False,title='Treatment')
#plt.savefig('ss_heatmap_all.pdf')
plt.savefig('ss_heatmap_all.png')
#plt.show()

##lncap
tcolors = {'DMSO':'gray','M':'r','T':'b','TM':'magenta',
          'W':'gold','TW':'green','MW':'orange'
          }
hours = ['0','3','6','9','12','24']
ltxs = [tx+'_'+hour for tx in tcolors.keys() for hour in hours]
lcolors = dict()
for tx in ltxs:
    lcolors[tx] = tcolors[tx.split('_')[0]]

#lmeans = lmeans.loc[ln_genes][ltxs]
lmeans = lmeans[ltxs]
cg = sns.clustermap(lmeans,cmap='bwr',col_cluster=False,
               z_score=0,
               yticklabels=False,col_colors=lmeans.columns.map(lcolors),
               cbar_kws={'label':'z-score(cpm)'},vmin=-2.5,vmax=2.5)#)
cg.ax_row_dendrogram.set_visible(False)
# Use the dendrogram box to reposition the colour bar
dendro_box = cg.ax_row_dendrogram.get_position()
dendro_box.x0 = (dendro_box.x0 + 2*dendro_box.x1) / 3
dendro_box.x1 = (dendro_box.x0 + dendro_box.x1)/2

lbox = cg.ax_col_dendrogram.get_position()
lbox.x0 = dendro_box.x0*0.75
lbox.x1 = dendro_box.x1
lbox.y1 = dendro_box.y1*1.05
#dendro_box.y0 = dendro_box.y0*1.5
dendro_box.y1 = dendro_box.y1*0.7
lbox.y0 = dendro_box.y1
cg.cax.set_position(dendro_box)
# Move the ticks to the left (https://stackoverflow.com/a/36939552/1878788)
cg.cax.yaxis.set_ticks_position("left")
cg.cax.yaxis.set_label_position('left')
#color legend
cg.ax_col_dendrogram.set_position(lbox)
cg.ax_col_dendrogram.legend(
        [Line2D([0],[0],
                color='w',markerfacecolor=value,marker='s') for value in tcolors.values()],
        list(tcolors.keys()),loc='upper left',frameon=False,title='Treatment')
#plt.savefig('ln_heatmap_all.pdf')
plt.savefig('ln_heatmap_all.png')
#plt.show()

##mcf7
#mmeans = mmeans.loc[mc_genes][ltxs]
mmeans = mmeans[ltxs]
cg = sns.clustermap(mmeans,cmap='bwr',col_cluster=False,
               z_score=0,
               yticklabels=False,col_colors=lmeans.columns.map(lcolors),
               cbar_kws={'label':'z-score(cpm)'},vmin=-2.5,vmax=2.5)#)
cg.ax_row_dendrogram.set_visible(False)
# Use the dendrogram box to reposition the colour bar
dendro_box = cg.ax_row_dendrogram.get_position()
dendro_box.x0 = (dendro_box.x0 + 2*dendro_box.x1) / 3
dendro_box.x1 = (dendro_box.x0 + dendro_box.x1)/2

lbox = cg.ax_col_dendrogram.get_position()
lbox.x0 = dendro_box.x0*0.75
lbox.x1 = dendro_box.x1
lbox.y1 = dendro_box.y1*1.05
#dendro_box.y0 = dendro_box.y0*1.5
dendro_box.y1 = dendro_box.y1*0.7
lbox.y0 = dendro_box.y1
cg.cax.set_position(dendro_box)
# Move the ticks to the left (https://stackoverflow.com/a/36939552/1878788)
cg.cax.yaxis.set_ticks_position("left")
cg.cax.yaxis.set_label_position('left')
#color legend
cg.ax_col_dendrogram.set_position(lbox)
cg.ax_col_dendrogram.legend(
        [Line2D([0],[0],
                color='w',markerfacecolor=value,marker='s') for value in tcolors.values()],
        list(tcolors.keys()),loc='upper left',frameon=False,title='Treatment')
#plt.savefig('mc_heatmap_all.pdf')
plt.savefig('mc_heatmap_all.png')
#plt.show()

