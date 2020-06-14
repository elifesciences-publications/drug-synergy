#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 30 10:42:58 2018

@author: jenniferlongdiaz
"""

#eobs for mcf7 and lncap

import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy import stats

plt.ioff()

##set fonts
font = {'fontname':'Arial','size':10}

lncap = pd.read_csv('../LNCaP/LNCaP_combinations0.1.csv')
mcf7 = pd.read_csv('../diff_limma_december/MCF7_combinations0.1.csv')
ss = pd.read_csv('../self_synergy_data/self_synergy_combination_degs0.1.csv')
ss = ss[ss.combo!='TM_24']

both = pd.concat([lncap,mcf7,ss])
both = both[both.Hours>0]

both['lognsyn'] =np.log(both.nsyn)

for col in ['apoptotic signaling pathway (GO:0097190)',
            'intrinsic apoptotic signaling pathway (GO:0097193)',
            'intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress (GO:0070059)',
            'generation of precursor metabolites and energy (GO:0006091)']:
    both[col.split(' ')[-1]+'_synergistic'] = both[col]/both.nsyn

both['apoptosis-energy'] = both['apoptotic signaling pathway (GO:0097190)']+\
    both['generation of precursor metabolites and energy (GO:0006091)']

both['energy-enrichment']=both[['energy-enrichment-down','energy-enrichment-up']].max(axis=1).fillna(0)

colors = {'MW':'orange','TW':'green','TM':'magenta','M':'r','T':'b'}
#markers = {'LNCaP':'x','MCF7':'o'}

both['color']=both.apply(lambda row: colors[row.combo.split('_')[0]],axis=1)
#both['marker']=both.apply(lambda row: markers[row.cell_line],axis=1)

mcf = both[both.cell_line=='MCF7']
ss = mcf[mcf.Drug1.str.len()<6]
mcf = mcf[mcf.Drug1.str.len()>5]
ss['alpha'] = ss.combo.str.split('_').str[1].astype(int)/25
ln = both[both.cell_line=='LNCaP']

##corr vs nsyn
print('corr vs nsyn: ',stats.spearmanr(both.pearson_r,both.nsyn))

fig,ax = plt.subplots()
ax.scatter(mcf.pearson_r,mcf.nsyn,s=mcf.Hours.astype(int)+1,color=mcf.color,
            marker='o')
ax.scatter(ln.pearson_r,ln.nsyn,s=ln.Hours.astype(int)+1,color=ln.color,
            marker='x')
for i,row in ss.iterrows():
    ax.scatter(row.pearson_r,row.nsyn,s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
ax.set_xlabel('Correlation of Monotherapies')
ax.set_ylabel('Number of Synergistic Genes')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
           
            [Line2D([0],[0],marker='o',color='w',
                   markerfacecolor='w',
                   markeredgecolor='k',
                   markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
    [Line2D([0],[0],marker=m,color='w',
                   markerfacecolor='k') for m in ['o','x']]+
    [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                   markerfacecolor='k',
                   markeredgecolor='k') for dose in ss.combo.str.split('_').str[1].astype(int).unique()],
           list(colors.keys())+
           [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
           [x+' uM' for x in ss.combo.str.split('_').str[1].unique()],
    loc='lower left',bbox_to_anchor=(1,0),ncol=1) #markersize=15
plt.title('All In-House Data: Synergistic Genes vs Correlation')
plt.savefig('corr_vs_nsyn.png',dpi=300)

##eob plots
b = both[~both.combo.isin(ss.combo)]
print('corr vs eob spearman: ',stats.spearmanr(both.pearson_r,both.EOB))
print('corr vs eob combos only spearman: ', stats.spearmanr(b.pearson_r,b.EOB))
print('nsyn vs eob spearman: ',stats.spearmanr(both.nsyn,both.EOB))
print('lognsyn vs eob spearman: ',stats.spearmanr(both.lognsyn,both.EOB))

print('corr vs eob pearson: ',stats.pearsonr(both.pearson_r,both.EOB))
print('corr vs eob combos only spearman: ', stats.pearsonr(b.pearson_r,b.EOB))
print('nsyn vs eob pearson: ',stats.pearsonr(both.nsyn,both.EOB))
print('lognsyn vs eob pearson: ',stats.pearsonr(both.lognsyn,both.EOB))

tmp = ss[ss.combo!='M_5']
col = ['lognsyn','Number of Synergistic Genes (log)','Synergistic Genes(log)']
fig,ax = plt.subplots()
ax.scatter(mcf[col[0]],mcf.EOB,s=mcf.Hours.astype(int)+1,color=mcf.color,
            marker='o')
ax.scatter(ln[col[0]],ln.EOB,s=ln.Hours.astype(int)+1,color=ln.color,
            marker='x')
for i,row in tmp.iterrows():
    ax.scatter(row[col[0]],row.EOB,s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
#    if col in cols[-3:]:
#        ax.set_xlim(0,0.1)
ax.set_xlabel(col[1])
ax.set_ylabel('Excess Over Bliss')
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
           
            [Line2D([0],[0],marker='o',color='w',
                   markerfacecolor='w',
                   markeredgecolor='k',
                   markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
    [Line2D([0],[0],marker=m,color='w',
                   markerfacecolor='k',markeredgecolor='k') for m in ['o','x']]+
    [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                   markerfacecolor='k',
                   markeredgecolor='k') for dose in tmp.combo.str.split('_').str[1].astype(int).unique()],
           list(colors.keys())+
           [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
           [x+' uM' for x in tmp.combo.str.split('_').str[1].unique()],
    loc='lower left',bbox_to_anchor=(1,0),ncol=1) #markersize=15
plt.title('All In-House Data: EOB vs '+col[2])
plt.savefig(col[2].split(' ')[-1]+'_vs_EOB.png',dpi=300)

cols = [['pearson_r','Correlation of Monotherapies','Correlation'],
        ['nsyn','Number of Synergistic Genes','Synergistic Genes'],
       ['apoptotic signaling pathway (GO:0097190)','Number of Synergistic Apoptosis Genes','Synergistic Apoptosis'],
       ['intrinsic apoptotic signaling pathway (GO:0097193)',
        'Number of Synergistic Intrinsic Apoptosis Genes','Synergistic GO:0097193'],
       ['intrinsic apoptotic signaling pathway in response to endoplasmic reticulum stress (GO:0070059)',
       'Number of Synergistic Apoptosis-ER Stress Genes','Synergistic GO:0070059'],
       ['generation of precursor metabolites and energy (GO:0006091)',
       'Number of Synergistic Energy Genes','Synergistic GO: 0006091'],
        ['apoptosis-energy','apoptosis-energy','apoptosis-energy'],
        ['energy-enrichment','energy-enrichment','energy-enrichment']]

       
for col in cols:
    fig,ax = plt.subplots()
    ax.scatter(mcf[col[0]],mcf.EOB,s=mcf.Hours.astype(int)+1,color=mcf.color,
                marker='o')
    ax.scatter(ln[col[0]],ln.EOB,s=ln.Hours.astype(int)+1,color=ln.color,
                marker='x')
    for i,row in ss.iterrows():
        ax.scatter(row[col[0]],row.EOB,s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
#    if col in cols[-3:]:
#        ax.set_xlim(0,0.1)
    ax.set_xlabel(col[1])
    ax.set_ylabel('Excess Over Bliss')
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
               
                [Line2D([0],[0],marker='o',color='w',
                       markerfacecolor='w',
                       markeredgecolor='k',
                       markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
        [Line2D([0],[0],marker=m,color='w',
                       markerfacecolor='k',markeredgecolor='k') for m in ['o','x']]+
        [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                       markerfacecolor='k',
                       markeredgecolor='k') for dose in ss.combo.str.split('_').str[1].astype(int).unique()],
               list(colors.keys())+
               [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
               [x+' uM' for x in ss.combo.str.split('_').str[1].unique()],
        loc='lower left',bbox_to_anchor=(1,0),ncol=1) #markersize=15
    plt.title('All In-House Data: EOB vs '+col[2])
    plt.savefig(col[2].split(' ')[-1]+'_vs_EOB.png',dpi=300)

#viability plots
cols = [['pearson_r','Correlation of Monotherapies','Correlation'],
        ['nsyn','Number of Synergistic Genes','Synergistic Genes'],
        ['EOB','Excess Over Bliss','EOB']]
for col in cols:
    fig,ax = plt.subplots()
    ax.scatter(mcf['Observed Viability'],mcf[col[0]],s=mcf.Hours.astype(int)+1,color=mcf.color,
                marker='o')
    ax.scatter(ln['Observed Viability'],ln[col[0]],s=ln.Hours.astype(int)+1,color=ln.color,
                marker='x')
    for i,row in ss.iterrows():
        ax.scatter(row['Observed Viability'],row[col[0]],s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
    ax.set_xlabel('Viability')
    ax.set_ylabel(col[1])
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
               
                [Line2D([0],[0],marker='o',color='w',
                       markerfacecolor='w',
                       markeredgecolor='k',
                       markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
        [Line2D([0],[0],marker=m,color='w',
                       markerfacecolor='k',markeredgecolor='k') for m in ['o','x']]+
        [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                       markerfacecolor='k',
                       markeredgecolor='k') for dose in ss.combo.str.split('_').str[1].astype(int).unique()],
               list(colors.keys())+
               [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
               [x+' uM' for x in ss.combo.str.split('_').str[1].unique()],
        loc='lower left',bbox_to_anchor=(1,0),ncol=1) #markersize=15
    plt.title('All In-House Data: '+col[2]+' vs Viability')
    plt.savefig(col[0]+'_vs_viab.png',dpi=300)

    
cols = [['(GO:0097190)_synergistic','apoptosis/synergistic','(GO:0097190)_synergistic'],
       ['(GO:0097193)_synergistic','intrinsic-apoptosis/synergistic','(GO:0097193)_synergistic'],
       ['(GO:0070059)_synergistic','ER-apoptosis/synergistic','(GO:0070059)_synergistic'],
        ['(GO:0006091)_synergistic','energy/synergistic','(GO:0006091)_synergistic']]
    
for col in cols:
    fig,ax = plt.subplots()
    ax.scatter(mcf[col[0]],mcf.EOB,s=mcf.Hours.astype(int)+1,color=mcf.color,
                marker='o')
    ax.scatter(ln[col[0]],ln.EOB,s=ln.Hours.astype(int)+1,color=ln.color,
                marker='x')
    for i,row in ss.iterrows():
        ax.scatter(row[col[0]],row.EOB,s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
    ax.set_xlim(0,0.1)
    ax.set_xlabel(col[1])
    ax.set_ylabel('Excess Over Bliss')
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
               
                [Line2D([0],[0],marker='o',color='w',
                       markerfacecolor='w',
                       markeredgecolor='k',
                       markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
        [Line2D([0],[0],marker=m,color='w',
                       markerfacecolor='k',markeredgecolor='k') for m in ['o','x']]+
        [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                       markerfacecolor='k',
                       markeredgecolor='k') for dose in ss.combo.str.split('_').str[1].astype(int).unique()],
               list(colors.keys())+
               [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
               [x+' uM' for x in ss.combo.str.split('_').str[1].unique()],
        loc='lower left',bbox_to_anchor=(1,0),ncol=1) #markersize=15
    plt.title('All In-House Data: EOB vs '+col[2])
    plt.savefig(col[2].split(' ')[-1]+'_vs_EOB.png',dpi=300)
    
    


###PLOTS FOR PAPER
tmp = ss[ss.combo!='M_5']
col = ['lognsyn','Number of Synergistically Expressed Genes (log)','Synergistically Expressed Genes (log)']
fig,ax = plt.subplots()
ax.scatter(mcf[col[0]],mcf.EOB,s=mcf.Hours.astype(int)+1,color=mcf.color,
            marker='o')
ax.scatter(ln[col[0]],ln.EOB,s=ln.Hours.astype(int)+1,color=ln.color,
            marker='x')
for i,row in tmp.iterrows():
    ax.scatter(row[col[0]],row.EOB,s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
#    if col in cols[-3:]:
#        ax.set_xlim(0,0.1)
ax.set_xlabel(col[1],**font)
ax.set_ylabel('Excess Over Bliss',**font)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
           
            [Line2D([0],[0],marker='o',color='w',
                   markerfacecolor='w',
                   markeredgecolor='k',
                   markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
    [Line2D([0],[0],marker=m,color='w',
                   markerfacecolor='k',markeredgecolor='k') for m in ['o','x']]+
    [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                   markerfacecolor='k',
                   markeredgecolor='k') for dose in tmp.combo.str.split('_').str[1].astype(int).unique()],
           list(colors.keys())+
           [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
           [x+' uM' for x in tmp.combo.str.split('_').str[1].unique()],
    loc='lower left',bbox_to_anchor=(1,0),ncol=1,prop={'family':'Arial'}) #markersize=15
plt.title('Synergy vs '+col[2],**font)
plt.savefig(col[2].split(' ')[-1]+'_vs_EOB.pdf')



col = ['pearson_r','Pearson r of Differentially Expressed Genes','Correlation of Monotherapies']

fig,ax = plt.subplots()
ax.scatter(mcf[col[0]],mcf.EOB,s=mcf.Hours.astype(int)+1,color=mcf.color,
            marker='o')
ax.scatter(ln[col[0]],ln.EOB,s=ln.Hours.astype(int)+1,color=ln.color,
            marker='x')
for i,row in ss.iterrows():
    ax.scatter(row[col[0]],row.EOB,s=int(row.Hours)+1,color=row.color,alpha=row.alpha)
#    if col in cols[-3:]:
#        ax.set_xlim(0,0.1)
ax.set_xlabel(col[1],**font)
ax.set_ylabel('Excess Over Bliss',**font)
# Shrink current axis by 20%
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
# Put a legend to the right of the current axis
ax.legend([Line2D([0],[0],color=value) for value in colors.values()]+
           
            [Line2D([0],[0],marker='o',color='w',
                   markerfacecolor='w',
                   markeredgecolor='k',
                   markersize=np.sqrt(sz)) for sz in (both.Hours.astype(int)+1).unique()]+
    [Line2D([0],[0],marker=m,color='w',
                   markerfacecolor='k',markeredgecolor='k') for m in ['o','x']]+
    [Line2D([0],[0],alpha=dose/25,marker='o',color='w',
                   markerfacecolor='k',
                   markeredgecolor='k') for dose in ss.combo.str.split('_').str[1].astype(int).unique()],
           list(colors.keys())+
           [str(int(x))+' hours' for x in both.Hours.unique()]+['MCF7','LNCaP']+
           [x+' uM' for x in ss.combo.str.split('_').str[1].unique()],
    loc='lower left',bbox_to_anchor=(1,0),ncol=1,prop={'family':'Arial'}) #markersize=15
plt.title('Synergy vs '+col[2],**font)
plt.savefig(col[2].split(' ')[-1]+'_vs_EOB.pdf')
