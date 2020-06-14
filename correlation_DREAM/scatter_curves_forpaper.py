#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 01:22:26 2019

@author: jenniferlongdiaz

all the different learning methods on DREAM
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn import ensemble
from sklearn import metrics
from scipy import stats
import seaborn as sns
import os
from itertools import combinations

sns.reset_defaults()

#import EOBs
eob = pd.read_csv('../dream_raw/dream/NCI-DREAMSubchallenge2_GoldStandard.txt',sep='\t')
eob.index = eob['Cmpd A'].str.replace('-','_')+'-'+eob['Cmpd B'].str.replace('-','_')
eob.drop(['Cmpd A','Cmpd B'],axis=1,inplace=True)
eob.columns = ['EOB','SEM']
eob.EOB = -eob.EOB
eob['snr']=eob.EOB.abs()/eob.SEM
eob['synergistic'] = np.where((eob.EOB>0)&(eob.snr>2),1,0).astype(int)

#import digre
digre = pd.read_csv('../dream_raw/Team Predictions/Rank1.csv')
digre.index=digre['Drug Combination'].str.replace('-','_').str.replace(
        ' & ','-')
digre.drop('Compound pair with additive activity (IC36)',axis=0,inplace=True)
digre.Rank = digre.Rank.astype(int)
digre['Digre_inverted'] = 92-digre.Rank #the ranks are flipped for some reason

#set up dataframe for correlations
pairs = pd.DataFrame(index=eob.index)
pairs['drug1'] = pairs.index.str.split('-').str[0]
pairs['drug2'] = pairs.index.str.split('-').str[1]

#calculate correlations    
def get_corr_new(row,function,folder,param,p_mask):
    df1 = pd.read_csv('limma_new/'+folder+'/limma_'+row.drug1.replace(
            'H_7','H-7, Dihydrochloride').replace('Doxorubicin','Doxorubicin hydrochloride')+'.csv',
                      index_col='Unnamed: 0')
    df2 = pd.read_csv('limma_new/'+folder+'/limma_'+row.drug2.replace(
            'H_7','H-7, Dihydrochloride').replace('Doxorubicin','Doxorubicin hydrochloride')+'.csv',
                      index_col='Unnamed: 0').reindex(df1.index)
    s = set(df1[df1['adj.P.Val']<p_mask].index).union(
            set(df2[df2['adj.P.Val']<p_mask].index))
    return pd.Series(function(df1[df1.index.isin(s)][param],
                                  df2[df2.index.isin(s)][param]))

param='logFC'
pval = 0.1

pairs[['spearman_r_mask_new_probe','spearman_p_mask_new_probe']]=pairs.apply(get_corr_new,axis=1,
     function=stats.spearmanr,folder='probe',param=param,p_mask=pval) #'logFC'
pairs[['pearson_r_mask_new_probe','pearson_p_mask_new_probe']]=pairs.apply(get_corr_new,axis=1,
     function=stats.pearsonr,folder='probe',param=param,p_mask=pval)

pairs[['spearman_r_new_probe','spearman_p_new_probe']]=pairs.apply(get_corr_new,axis=1,
     function=stats.spearmanr,folder='probe',param=param,p_mask=1) #'logFC'
pairs[['pearson_r_new_probe','pearson_p_new_probe']]=pairs.apply(get_corr_new,axis=1,
     function=stats.pearsonr,folder='probe',param=param,p_mask=1)

pairs[['spearman_r_mask_new_gene','spearman_p_mask_new_gene']]=pairs.apply(get_corr_new,axis=1,
     function=stats.spearmanr,folder='gene',param=param,p_mask=pval) #'logFC'
pairs[['pearson_r_mask_new_gene','pearson_p_mask_new_gene']]=pairs.apply(get_corr_new,axis=1,
     function=stats.pearsonr,folder='gene',param=param,p_mask=pval)

pairs[['spearman_r_new_gene','spearman_p_new_gene']]=pairs.apply(get_corr_new,axis=1,
     function=stats.spearmanr,folder='gene',param=param,p_mask=1) #'logFC'
pairs[['pearson_r_new_gene','pearson_p_new_gene']]=pairs.apply(get_corr_new,axis=1,
     function=stats.pearsonr,folder='gene',param=param,p_mask=1)

#put everything together for plotting
pairs = pd.concat([pairs,digre[['Digre_inverted']],eob],axis=1)

corr = pairs.pearson_r_mask_new_gene

###FOR PAPER
pfont = {'fontname':'Arial','fontsize':10}
lfont = {'family':'Arial'} #,'fontsize':10
ifont = {'fontname':'Arial','fontsize':30}

i=0
for metric, label in [[corr,'Correlation'],[pairs.Digre_inverted,'DIGRE']]:#,'min_rf','max_rf']:
    print(metric)
    fpr,tpr,thresh, = metrics.roc_curve(pairs.synergistic,metric)
    print('AUROC',metrics.auc(fpr,tpr))
    plt.figure(12)
    p=plt.plot(fpr,tpr,label=label)
    plt.text(0.29,0.8-0.2*i,'AUC '+str(round(metrics.auc(fpr,tpr),2))[:4],
             color = p[-1].get_color(),**pfont)
    

    
    precision, recall, thresh = metrics.precision_recall_curve(pairs.synergistic,
                                                               metric)
    print('AUPR',metrics.auc(recall,precision))
    plt.figure(13)
    q = plt.plot(recall,precision,label=label)
    plt.text(0.5,0.45-0.2*i,'AUC '+str(round(metrics.auc(recall,precision),2))[:4],
             color = q[-1].get_color(),**pfont)
    
    i+=1

plt.figure(12)
plt.plot([0,1],[0,1],'k',label='Random')
plt.legend(prop=lfont)
plt.xticks(**pfont)
plt.yticks(**pfont)
plt.xlabel('False Positive Rate',**pfont)
plt.ylabel('True Positive Rate',**pfont)
plt.title('Synergistic Drug Pair Predictions: ROC',**pfont)
plt.savefig('ROCforpaper.pdf')

plt.figure(13)
plt.plot([0,1],[eob.synergistic.sum()/eob.shape[0]]*2,'k',label='Random')
plt.legend(prop=lfont)
plt.xticks(**pfont)
plt.yticks(**pfont)
plt.xlabel('Recall',**pfont)
plt.ylabel('Precision',**pfont)
plt.title('Synergistic Drug Pair Predictions: PR',**pfont)
plt.savefig('PRforpaper.pdf')

#remake the scatterplot
plt.figure()
plt.scatter(corr,pairs.EOB,c='k')
plt.xticks(**pfont)
plt.yticks(**pfont)
plt.xlabel('Gene Expression Correlation of Monotherapies',**pfont)
plt.ylabel('Excess over Bliss (EOB)',**pfont)
plt.title('Phenotypic Synergy vs Correlation of Monotherapies',**pfont)
plt.savefig('scatter_forpaper.pdf')

#remake the boxplot
plt.figure()
sns.boxplot(data=[corr[pairs.EOB<-2.5],corr[pairs.EOB>2.5]],color='gray')
sns.swarmplot(data = [corr[pairs.EOB<-2.5],corr[pairs.EOB>2.5]],color='k',size=10)
plt.xticks([0,1],['EOB<-2.5','EOB>2.5'],**ifont)
plt.ylim([-0.4,0.8])
plt.yticks(**ifont)
plt.ylabel('Correlation',**ifont)
plt.tight_layout()
plt.savefig('boxplot_inset.pdf')
