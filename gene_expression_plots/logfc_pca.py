#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 23:05:49 2019

@author: jenniferlongdiaz

pca of logfc
"""

import pandas as pd
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import os

##set fonts
font = {'fontname':'Arial','size':10}

dfs = dict()
ldfs = dict()

for file in [f for f in os.listdir('../self_synergy_data') if 'ss_limma_eb' in f]:
    dfs['self_'+file.split('limma_eb_')[1].split('.csv')[0]] = pd.read_csv(
            '../self_synergy_data/'+file,index_col=0)
for file in [f for f in os.listdir('../diff_limma_december/') if 'limma_t0' in f]:
    dfs['MCF7_'+file.split('limma_t0_')[1].split('.txt')[0]] = pd.read_csv(
            '../diff_limma_december/'+file,sep='\t',index_col='1')

for file in [f for f in os.listdir('../LNCaP/LIMMA_DATA/') if 'limma' in f]:
    if not 'TMW' in file:
        ldfs[file.split('DE_limma_')[1].split('.csv')[0]] = pd.read_csv(
                '../LNCaP/LIMMA_DATA/'+file,index_col=0)

mcf7 = pd.DataFrame()
lncap = pd.DataFrame()

for key in dfs.keys():
    mcf7[key] = dfs[key].logFC
for key in ldfs.keys():
    lncap[key] = ldfs[key].logFC
    
 
######MCF7#####
#pca
pca = PCA(2)
pca.fit(mcf7.T)
comp = pd.DataFrame(pca.transform(mcf7.T),
                    index=mcf7.columns)

#plotting 
comp['combo'] = comp.index.str.split('_').str[1]                  
colors = {'M':'r','T':'b','W':'gold',
          'TM':'magenta','MW':'orange','TW':'green',
          }
comp['color'] = comp.combo.map(colors) 

mcf = comp[comp.index.str.contains('MCF7')]
mcf['hours'] = mcf.index.str.split('_').str[2].astype(int)
mcf.sort_values(['combo','hours'],inplace=True)
#mcf['alpha'] = np.where(mcf.combo=='M',0.4,
#   np.where(mcf.combo=='T',0.8,
#            np.where(mcf.combo=='W',0.2,1)))
ss = comp[comp.index.str.contains('self')]
ss['dose'] = ss.index.str.replace('DMSO_','DMSO_25_').str.split('_').str[2].astype(float)
ss['alpha'] = ss.dose/25.
ss.sort_values(['dose','combo'],inplace=True)
comp = pd.concat([mcf,ss])
      
#make the scatterplot
fig,ax =plt.subplots(figsize=(6,7))
ax.scatter(mcf[0],mcf[1],color=mcf.color,s=mcf.hours+1)
for i,row in ss.iterrows():
    ax.scatter(row[0],row[1],color=row.color,alpha=row.alpha,marker='v')
ax.set_xlabel('First Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[0]*100,1))+'% of variance',**font)
ax.set_ylabel('Second Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[1]*100,1))+'% of variance',**font)
ax.set_title('Prinicpal Components Analysis of MCF7 RNA-Seq',**font)
# Shrink current axis by 10%
ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.8])
# Put a legend to the right of the current axis
ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
            [Line2D([0],[0],marker='o',color='w',
                   markerfacecolor='w',
                   markeredgecolor='k',
                   markersize=np.sqrt(sz)) for sz in (mcf.hours+1).unique()]+
    [Line2D([0],[0],marker=m,color='w',
                   markerfacecolor='k') for m in ['o','v']]+
    [Line2D([0],[0],alpha=d/25,marker='v',color='w',
                   markerfacecolor='k',
                   markeredgecolor='k') for d in ss.dose.unique()],
           list(colors.keys())+
           [str(int(x))+' hours' for x in mcf.hours.astype(str).unique()]+['Combos','Doses']+
           [x+' '+ r'$\mu$'+'M' for x in ss.dose.astype(str).unique()],
    loc='upper center',bbox_to_anchor=(0.5,-0.1),ncol=4,prop={'family':'Arial'})

plt.savefig('mcf7_logfc_pca.pdf')


####LNCaP####
#pca
pca = PCA(2)
pca.fit(lncap.T)
comp = pd.DataFrame(pca.transform(lncap.T),
                    index=lncap.columns)

#plotting 
comp['combo'] = comp.index.str.split('_').str[0] 
comp['color'] = comp.combo.map(colors) 
comp['hours'] = comp.index.str.split('_').str[1].astype(int)
comp.sort_values(['combo','hours'],inplace=True)

#scatterplot
fig,ax =plt.subplots()
ax.scatter(comp[0],comp[1],color=comp.color,s=comp.hours+1)
ax.set_xlabel('First Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[0]*100,1))+'% of variance',**font)
ax.set_ylabel('Second Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[1]*100,1))+'% of variance',**font)
ax.set_title('Prinicpal Components Analysis of LNCaP RNA-Seq',**font)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
            [Line2D([0],[0],marker='o',color='w',
                   markerfacecolor='w',
                   markeredgecolor='k',
                   markersize=np.sqrt(sz)) for sz in (comp.hours+1).unique()],
           list(colors.keys())+
           [str(int(x))+' hours' for x in comp.hours.astype(str).unique()],
    loc='lower left',bbox_to_anchor=(1,0.1),ncol=1,prop={'family':'Arial'})

plt.savefig('lncap_logfc_pca.pdf')
