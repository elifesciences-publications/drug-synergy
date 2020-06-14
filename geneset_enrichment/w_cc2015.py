#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:41:33 2019

@author: jenniferlongdiaz

make a new withaferin figure to reduce redundancy
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import colors

sns.reset_defaults()
sns.set(font='Arial',font_scale=0.5,color_codes=False)

w = pd.read_csv('all_geneset_enrichment_fdr_log_cc2015.tsv',sep='\t',index_col='candset')

w = w[['down_W_3', 'down_W_6', 'down_W_9', 'down_W_12', 'down_W_24',
               'down_T_3', 'down_T_6', 'down_T_9', 'down_T_12', 'down_T_24',
               'down_TUW_3', 'down_TUW_6', 'down_TUW_9', 'down_TUW_12', 'down_TUW_24',
               'down_TW_3', 'down_TW_6', 'down_TW_9', 'down_TW_12', 'down_TW_24',
               'down_M_3', 'down_M_6', 'down_M_9', 'down_M_12', 'down_M_24',
        'down_MUW_3', 'down_MUW_6', 'down_MUW_9', 'down_MUW_12', 'down_MUW_24',
        'down_MW_3', 'down_MW_6', 'down_MW_9', 'down_MW_12', 'down_MW_24',
        'up_W_3', 'up_W_6','up_W_9', 'up_W_12', 'up_W_24',
            'up_T_3', 'up_T_6', 'up_T_9', 'up_T_12', 'up_T_24',
            'up_TUW_3', 'up_TUW_6', 'up_TUW_9','up_TUW_12', 'up_TUW_24',
            'up_TW_3', 'up_TW_6', 'up_TW_9', 'up_TW_12','up_TW_24',
            'up_M_3', 'up_M_6', 'up_M_9', 'up_M_12', 'up_M_24',
                'up_MUW_3', 'up_MUW_6', 'up_MUW_9','up_MUW_12', 'up_MUW_24',
                'up_MW_3', 'up_MW_6', 'up_MW_9', 'up_MW_12','up_MW_24']]

finalsets = [#
        'ER to Golgi transport vesicle (GO:0030134)',
        'Golgi membrane (GO:0000139)',
        
        
        #ER
        'endoplasmic reticulum membrane (GO:0005789)',
        'mitochondrial matrix (GO:0005759)',
        'extracellular vesicular exosome (GO:0070062)',
    
        
        'lysosome (GO:0005764)',
        'cytoplasmic vesicle membrane (GO:0030659)',
        'cytoplasmic vesicle part (GO:0044433)', 
        
        'pre-autophagosomal structure (GO:0000407)',
        'pre-autophagosomal structure membrane (GO:0034045)',
        
        
        'COPII vesicle coat (GO:0030127)',
        'TFEB_TFE3_Transcriptional Activators',
        
        'nuclear pore (GO:0005643)',
        'spliceosomal complex (GO:0005681)',
        'nucleolus (GO:0005730)',
        'nucleoplasm (GO:0005654)',
        'MCM complex (GO:0042555)',
        
        'PLD_Down', 'PLD_Up', 
       ]

w = w.T[finalsets].T
#w = w.reindex(index=w.index[::-1])
#plt.figure(figsize=(8,6))


#def make_rgb_transparent(color, alpha):
#    rgb = colors.colorConverter.to_rgb(color)
#    return [alpha * c1 + (1 - alpha) * c2
#            for (c1, c2) in zip(rgb, (1,1,1))]
#    
#let's stick with alpha = dose, use size for time
c = [el for ls in [[x]*5 for x in ['gold','b','gold','g','r','gold','orange']] for el in ls]*2

#clist = [make_rgb_transparent(color,alpha) for color in c for alpha in [0.2,0.4,0.6,0.8,1] ]

lut = dict(zip(w.columns.tolist(),c))

hm=sns.clustermap(w,row_cluster=False,col_cluster=False,cmap='hot_r',vmax=6,
               yticklabels=True,cbar_kws={'label':'-log10(fdr)'},figsize=(12,6),
               col_colors=w.columns.map(lut))

fig = hm.fig
fig.subplots_adjust(left=0.01,bottom=0.5,right=0.6,top=0.99)
fig.savefig('w_cc2015.pdf')
