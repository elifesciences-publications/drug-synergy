#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 12:57:35 2019

@author: jenniferlongdiaz
#combine MCF7 combos, self-synergy, and lncap combos before removing batch effects with limma
"""
import pandas as pd
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


##first combine the files
cpm = pd.read_csv('../log_cpm_viab.csv',index_col=0,header=0,skiprows=[1])
cpm.columns = ['MCF7_'+col for col in cpm.columns]
lcpm = pd.read_csv('../LNCaP/lncap_expression.csv',index_col=0,header=0)
#lcpm.columns = ['LNCaP_'+col for col in lcpm.columns]
scpm = pd.read_csv('../self_synergy_data/ss_log_cpm.csv',index_col=0,header=0)
scpm.columns = ['self_'+col for col in scpm.columns]

cpm = cpm.merge(scpm,how='inner',left_index=True,right_index=True)

cpm.to_csv('mcf7_log_cpm.csv')##use this file to remove batch effects in limma

rmeans = pd.DataFrame()
for tx in list(set([col[:-1] for col in cpm.columns])):
    rmeans[tx] = cpm[[col for col in cpm.columns if tx in col]].mean(axis=1)

rmeans.to_csv('mcf7_log_cpm_means.csv')
##import the cpm file with batch effects removed
cpm = pd.read_csv('mcf7_log_cpm_removedbatch.csv',index_col=0,header=0)
##import the means file with batch effects removed
means = pd.read_csv('mcf7_log_cpm_means_removedbatch.csv',index_col=0,header=0)

#means.columns = [col.split('_')[1]+'_'+col.split('_')[2] for col in means.columns]
#means.to_csv('mcf7_batchremoved_means.csv')

lmeans = pd.DataFrame()
for tx in list(set([col[:-1] for col in lcpm.columns])):
    lmeans[tx] = lcpm[[col for col in lcpm.columns if tx in col]].mean(axis=1)
lmeans = lmeans[[col for col in lmeans.columns if not 'TMW' in col]]
    
######MCF7#####
#pca
pca = PCA()
pca.fit(rmeans.T)
comp = pd.DataFrame(pca.transform(rmeans.T),
                    index=rmeans.columns)

#plotting 
comp['combo'] = comp.index.str.split('_').str[1]                  
    colors = {'DMSO':'k','M':'r','T':'b','W':'gold',
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
fig,ax =plt.subplots()
ax.scatter(mcf[0],mcf[1],color=mcf.color,s=mcf.hours+1)
for i,row in ss.iterrows():
    ax.scatter(row[0],row[1],color=row.color,alpha=row.alpha,marker='v')
ax.set_xlabel('First Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[0]*100,1))+'% of variance')
ax.set_ylabel('Second Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[1]*100,1))+'% of variance')
ax.set_title('Prinicpal Components Analysis of MCF7 RNA-Seq')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
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
           [x+' uM' for x in ss.dose.astype(str).unique()],
    loc='lower left',bbox_to_anchor=(1,-0.15),ncol=1)

plt.savefig('mcf7_pca.png')


####LNCaP####
#pca
pca = PCA(2)
pca.fit(lmeans.T)
comp = pd.DataFrame(pca.transform(lmeans.T),
                    index=lmeans.columns)

#plotting 
comp['combo'] = comp.index.str.split('_').str[0] 
comp['color'] = comp.combo.map(colors) 
comp['hours'] = comp.index.str.split('_').str[1].astype(int)
comp.sort_values(['combo','hours'],inplace=True)

#scatterplot
fig,ax =plt.subplots()
ax.scatter(comp[0],comp[1],color=comp.color,s=comp.hours+1)
ax.set_xlabel('First Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[0]*100,1))+'% of variance')
ax.set_ylabel('Second Prinicpal Component: '+ str(
        round(pca.explained_variance_ratio_[1]*100,1))+'% of variance')
ax.set_title('Prinicpal Components Analysis of LNCaP RNA-Seq')
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
    loc='lower left',bbox_to_anchor=(1,0.1),ncol=1)

plt.savefig('lncap_pca.pdf')

