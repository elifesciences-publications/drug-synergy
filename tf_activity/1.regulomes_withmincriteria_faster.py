#get active transcription factors according to their regulomes from MCF7 aracne network
#doing multiple fisher's exact tests and a FDR correction
#try imposing cutoffs on contingency table inputs
#august 31 2015

import pandas as pd
from scipy import stats
from statsmodels.sandbox.stats.multicomp import multipletests as mult
import numpy as np
import itertools
from datetime import datetime
startTime = datetime.now()

tx = 'TW'

reg = pd.read_csv('../network_new/MCF7_Woo_et_al-Li_et_al_genelevel_noregulatorsastargets_20160126_JELD.txt',
                  sep='\t')

de = pd.read_csv('../gene_expression/diff_limma_december/diffexp_'+ 
tx+'_alltimes_JELD.csv')

allRNAseq = set(pd.read_csv('../gene_expression/diff_limma_december/limma_t0_M_0.txt',header=0,
                            sep='\t',index_col='1').index)
regset = set(reg.Regulator).union(set(reg.Target))

inter = regset.intersection(allRNAseq)

reg = reg[(reg.Regulator.isin(inter))&(reg.Target.isin(inter))]

de = de[de.gene.isin(inter)]

dirs = ['up','down'] #

txtypes = ['Transcriptional Activator','Transcriptional Repressor']

def Fishers_mins_w_regsize(row,minlist,minov):

    direction = row.direction
    txtype = row.regulation_type
    regulator = row.regulator
    time = row.time
    
    dx = set(de[(de.direction==direction)&(de.time==time)].gene)
    notdx = inter - dx
    r = set(reg[(reg.Regulator==regulator)&(
    reg.Type==txtype)].Target)
    notr = set(reg[(reg.Regulator!=regulator)|(
    reg.Type!=txtype)].Target)
    
    ##contingency table: 
            ## [[diffex&inregulome,diffex&not in regulome],
            ## [not diffex&in regulome, not diffex & not in regulome]]
    if len(dx)>=minlist and len(r)>=minlist and len(dx.intersection(r))>=minov:
        cont = [[len(dx.intersection(r)),len(dx.intersection(notr))], \
                [len(notdx.intersection(r)),len(notdx.intersection(notr))]]
        p = stats.fisher_exact(cont)[1]
    else:
        p = np.nan
    return (p,len(r))

##set regulome and dysregulated genes list minimum length
minlist = 3
##set minimum overlap 
minov = 2

#create dataframe for tx with fisher's p values
print('creating dataframe...')
x = {'time':de.time.unique().tolist(),'regulator':reg.Regulator.unique().tolist(),'direction':dirs,
     'regulation_type':txtypes}
df = pd.DataFrame(list(itertools.product(*x.values())), columns=x.keys())
print('getting fishers results...')
temp = df.apply(Fishers_mins_w_regsize,axis=1,minlist = minlist, minov = minov)
fisher = pd.DataFrame(temp.tolist(),columns = ['Fisher_p_value','regulome_size'])

df = df.merge(fisher,how='outer',left_index=True,right_index=True)

#fdr only for the tests that were actually done
print('getting fdr...')
nan = df[np.isnan(df.Fisher_p_value)]
notnan = df[~np.isnan(df.Fisher_p_value)]
notnan['Fisher_FDR'] = mult(notnan['Fisher_p_value'],method='fdr_bh')[1]
df = pd.concat([nan,notnan])
df.Fisher_FDR = df.Fisher_FDR.astype(float)

#interpret on/off
def ONOFF(row):
    if row.Fisher_FDR < 0.05:
        if row.direction == 'up':
            if row.regulation_type == 'Transcriptional Activator':
                return 'ON'
            elif row.regulation_type == 'Transcriptional Repressor':
                return 'OFF'
        elif row.direction =='down':
            if row.regulation_type == 'Transcriptional Activator':
                return 'OFF'
            elif row.regulation_type == 'Transcriptional Repressor':
                return 'ON'
    else:
        return 'unchanged'

print('getting on/off...')
df['TF_state'] = df.apply(ONOFF,axis=1)

#write file
df.to_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'.csv',index=False)


print(str(len(df.regulator.unique()))+' total tfs')

print('ON/OFF TFs:')
for time in de.time.unique():
    timedf = df[df.time==time]
    print(str(len(timedf[timedf["Fisher_FDR"] < 0.05].regulator.unique()))+
    ' altered tfs at time '+str(time))

print(datetime.now() - startTime)

#935 total tfs

#M
##ON/OFF TFs:
##23 altered tfs at time 3
##50 altered tfs at time 6
##50 altered tfs at time 9
##52 altered tfs at time 12
##49 altered tfs at time 24

#T
##18 altered tfs at time 3
##32 altered tfs at time 6
##26 altered tfs at time 9
##49 altered tfs at time 12
##73 altered tfs at time 24

#TM
##ON/OFF TFs:
##117 altered tfs at time 3
##210 altered tfs at time 6
##223 altered tfs at time 9
##258 altered tfs at time 12
##243 altered tfs at time 24

#W
##ON/OFF TFs:
##0 altered tfs at time 0
##99 altered tfs at time 3
##176 altered tfs at time 6
##179 altered tfs at time 9
##190 altered tfs at time 12
##171 altered tfs at time 24

#TW
##ON/OFF TFs:
##0 altered tfs at time 0
##101 altered tfs at time 3
##156 altered tfs at time 6
##169 altered tfs at time 9
##160 altered tfs at time 12
##159 altered tfs at time 24

#MW
##ON/OFF TFs:
##0 altered tfs at time 0
##96 altered tfs at time 3
##155 altered tfs at time 6
##165 altered tfs at time 9
##175 altered tfs at time 12
##175 altered tfs at time 12
##157 altered tfs at time 24
