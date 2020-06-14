#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 16:48:07 2018

@author: jenniferlongdiaz
"""
##plot raw viability data

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import xlrd
import numpy as np

mpl.rcParams["errorbar.capsize"]  = 2

lncap = pd.read_excel('LnCap InCell DMSO Counts Viability_novonly.xlsx',
                                        sheet_name='DMSO')[['Hour',
                      'Measurement']].rename(columns={'Hour':'Hours'})
lncap['Treatment'] = 'DMSO'
l2 = pd.read_excel('LnCap InCell DMSO Counts Viability_novonly.xlsx',
                                        sheet_name='Drugs Alone')[['Hours',
                 'Compound','Measurement']].rename(columns={'Compound':'Treatment'})
l2 = l2[~l2.Treatment.str.contains(' - ')]
lncap = pd.concat([lncap,l2])
l2 = pd.read_excel('LnCap InCell DMSO Counts Viability_novonly.xlsx',
                                        sheet_name='Combinations')[['Hours',
                 'Drug1','Drug2','Measurement']]
l2 = l2[~l2.Drug1.str.contains(' - ')]
l2['Treatment'] = l2.Drug1.str[0]+l2.Drug2.str[0]
lncap = pd.concat([lncap,l2[['Hours','Treatment','Measurement']]])
lncap['Cell_line'] = 'LNCaP'


mcf = pd.read_excel('MCF7 InCell DMSO Counts Viability.xlsx',
                                        sheet_name='DMSO')[['Hour',
                      'Measurement']].rename(columns={'Hour':'Hours'})
mcf['Treatment'] = 'DMSO'
mcf = pd.concat([mcf,pd.read_excel('MCF7 InCell DMSO Counts Viability.xlsx',
                                        sheet_name='Drugs Alone')[['Hour',
                 'Compound','Measurement']].rename(columns={'Hour':'Hours',
                    'Compound':'Treatment'})])
m2 = pd.read_excel('MCF7 InCell DMSO Counts Viability.xlsx',
                                        sheet_name='Combination Wells')[['Hours',
                 'Drug1','Drug2','Measurement']]
m2['Treatment'] = m2.Drug1.str[0]+m2.Drug2.str[0]
mcf = pd.concat([mcf,m2[['Hours','Treatment','Measurement']]])
mcf['Cell_line'] = 'MCF7'

both = pd.concat([lncap,mcf])

#now normalize everything to time 0 
both['Fold_change']=both.apply(lambda row: row.Measurement/both[(both.Cell_line==row.Cell_line)&(
        both.Treatment==row.Treatment)&(both.Hours==0)].Measurement.mean(),axis=1)

both = both[both.Hours<48]
means = both.groupby(['Cell_line', 'Treatment','Hours'],as_index=False).mean()
stds = both.groupby(['Cell_line', 'Treatment','Hours']).std().reset_index()
stds.columns = [x+'_std' for x in stds.columns]
means = pd.concat([means,stds[['Measurement_std','Fold_change_std']]],axis=1)

colordict = {'DMSO':'k','Tamoxifen':'blue','Mefloquine':'r','TM':'magenta',
             'Withaferin A':'y','TW':'green','MW':'orange'}

l = means[means.Cell_line=='LNCaP']
m = means[means.Cell_line=='MCF7']

#plot the raw cell counts for everything
plt.figure()
ax = plt.subplot(111)
for treatment in colordict.keys():
    tm = m[m.Treatment==treatment] 
    tl = l[l.Treatment==treatment]
    ax.errorbar(tm.Hours,tm.Measurement,tm.Measurement_std,linestyle='-',
                c=colordict[treatment],label='MCF7 '+treatment)
    ax.errorbar(tl.Hours,tl.Measurement,tl.Measurement_std,linestyle='--',
                c=colordict[treatment],label='LNCaP '+treatment)
ax.set_xticks(means.Hours.unique())
ax.set_xticklabels(means.Hours.unique())
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Cell count')
box = ax.get_position()
ax.set_position([box.x0,box.y0,box.width*0.8,box.height])
#plt.legend(loc='best',bbox_to_anchor=(1,0.5))
ax.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.show()

#split lncap and mcf and plot fold change
plt.figure()
ax = plt.subplot(111)
for treatment in colordict.keys():
    tm = m[m.Treatment==treatment]
    ax.errorbar(tm.Hours,tm.Fold_change,tm.Fold_change_std,linestyle='-',
                c=colordict[treatment],label=treatment)
ax.set_xticks(means.Hours.unique())
ax.set_xticklabels(means.Hours.unique())
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Fold change over time 0')
ax.set_title('MCF7')
box = ax.get_position()
ax.set_position([box.x0,box.y0,box.width*0.8,box.height])
#plt.legend(loc='best',bbox_to_anchor=(1,0.5))
ax.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.show()

plt.figure()
ax = plt.subplot(111)
for treatment in colordict.keys():
    tl = l[l.Treatment==treatment]
    ax.errorbar(tl.Hours,tl.Fold_change,tl.Fold_change_std,linestyle='--',
                c=colordict[treatment],label=treatment)
ax.set_xticks(means.Hours.unique())  
ax.set_xticklabels(means.Hours.unique())
ax.set_xlabel('Time (hours)')
ax.set_ylabel('Fold change over time 0')
ax.set_title('LNCaP')
box = ax.get_position()
ax.set_position([box.x0,box.y0,box.width*0.8,box.height])
#plt.legend(loc='best',bbox_to_anchor=(1,0.5))
ax.legend(loc='center left',bbox_to_anchor=(1,0.5))
plt.show()


#remake the standard viability plots
dmso = lncap[lncap.Treatment=='DMSO'].groupby('Hours').mean()
lncaptx = lncap[lncap.Treatment!='DMSO']
lncaptx['Fold_Change_over_DMSO'] = lncaptx.apply(lambda row: row.Measurement/dmso.loc[row.Hours,
       'Measurement'],axis=1)*100

lnfc = lncaptx[['Treatment','Hours','Fold_Change_over_DMSO']].groupby(
        ['Treatment','Hours'],as_index=False).mean()
lnfcs = lncaptx[['Treatment','Hours','Fold_Change_over_DMSO']].groupby(
        ['Treatment','Hours']).std().reset_index()
lnfcs.rename(columns={'Fold_Change_over_DMSO':'Fold_Change_DMSO_STD'},inplace=True)
lnfcs['SEM'] = lnfcs.Fold_Change_DMSO_STD/np.sqrt(1.5) ## <--the files from Ron and Chuck have this factor in the STD
