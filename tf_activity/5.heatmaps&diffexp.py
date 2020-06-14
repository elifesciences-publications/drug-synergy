#make heatmaps and correlate activity with diffexp
#march 14 2016

import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

def signReverse(row):
    if row.TF_state=='ON':
        return row.logFDR
    elif row.TF_state=='OFF':
        return -row.logFDR

txs = ['M','T','W','TM','MW','TW']

times = ['3','6','9','12','24']

de = pd.DataFrame()
for tx in txs:
    for time in times:
        name = '../gene_expression/diff_limma_december/limma_t0_'+tx+'_'+time+'.txt'
        df = pd.read_csv(name,header=0,sep='\t')#,index_col='1'
        df.rename(columns={'1':'regulator'},inplace=True)
        df['time'] = time
        df['treatment'] = tx
        de = pd.concat([de,df])

##act = pd.DataFrame()
##shad = pd.DataFrame()
cat = pd.DataFrame()
for tx in ['T','M','TM']:
#make heatmap
    final = pd.read_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv')
    final['Unique_Identifier'] = final.regulator + '_' + final.regulation_type.str.split(' ').str[1]
    final['logFDR'] = -np.log10(final.Fisher_FDR)
    final['logFDR_toplot'] = final.apply(signReverse,axis=1)

    df = final[['Unique_Identifier','time','logFDR_toplot']]
    df = df.pivot('Unique_Identifier','time','logFDR_toplot')
    df.columns = [tx+'_'+str(int(col)) for col in df.columns]
    
    cat = pd.concat([cat,df],axis=1)
cat.fillna(0.,inplace=True)
cat.to_csv('MCF7_aracne_T-M-TM_shadowed_TFs_toplot.tsv',sep='\t')
tm = cat.copy()

cat = pd.DataFrame()
for tx in ['W','MW','TW']:
#make heatmap
    final = pd.read_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv')
    final['Unique_Identifier'] = final.regulator + '_' + final.regulation_type.str.split(' ').str[1]
    final['logFDR'] = -np.log10(final.Fisher_FDR)
    final['logFDR_toplot'] = final.apply(signReverse,axis=1)

    df = final[['Unique_Identifier','time','logFDR_toplot']]
    df = df.pivot('Unique_Identifier','time','logFDR_toplot')
    df.columns = [tx+'_'+str(int(col)) for col in df.columns]

    cat = pd.concat([cat,df],axis=1)
cat.fillna(0.,inplace=True)
cat.to_csv('MCF7_aracne_W-MW-TW_shadowed_TFs_toplot.tsv',sep='\t')

al = tm.merge(cat,left_index=True,right_index=True,how='outer').fillna(0).applymap(abs)
al['TF'] = al.index.str.split('_').str[0]

new = pd.DataFrame()
for tf in al.TF.unique():
    df = al[al.TF==tf]
    new = pd.concat([new,pd.DataFrame(df.max()).T])

new.set_index('TF',drop=True,inplace=True)
new.to_csv('MCF7_aracne_all_final_TFs.tsv',sep='\t')

##    final['treatment'] = tx
##
##    original = pd.read_csv('MCF7_aracne_'+tx+'_regulators_fishers_mins_3_2.csv')
##
##    original['logFDR'] = -np.log10(original.Fisher_FDR)
##
##    original['treatment'] = tx
##
##    act = pd.concat([act,original])
##    shad = pd.concat([shad,final])
##
##act['Unique_Identifier'] = act.regulator + '_'+ act.time.astype(int).astype(str)+'_'+act.treatment
##act = act.set_index('Unique_Identifier')[['regulator','time','treatment',
##                                          'regulation_type','TF_state','Fisher_FDR','logFDR']]
##
##shad['Unique_Identifier'] = shad.regulator + '_'+ shad.time.astype(int).astype(str)+'_'+shad.treatment
##shad = shad.set_index('Unique_Identifier')[['regulator','time','treatment',
##                                          'regulation_type','TF_state','Fisher_FDR','logFDR']]
##
##act['Full_ID'] = act.regulator + '_' + act.regulation_type.str.split(' ').str[1] + '_'+ \
##                 act.time.astype(int).astype(str)+'_'+act.treatment
##
##shad['Full_ID'] = shad.regulator + '_' + shad.regulation_type.str.split(' ').str[1] + '_'+ \
##                  shad.time.astype(int).astype(str)+'_'+shad.treatment
##
##act = act[~act.Full_ID.isin(shad.Full_ID)]
##removed = act[act.TF_state!='unchanged']
##removed.TF_state = 'unchanged'
##removed['Removed After Shadow'] = True
##act = pd.concat([act[act.TF_state=='unchanged'],removed])
##print(act[act.TF_state!='unchanged'])
##
##act = pd.concat([act,shad])

###correlate p values for all TFs
##de['logFDR_de'] = -np.log10(de['adj.P.Val'])
##de['Unique_Identifier'] = de.regulator+'_'+de.time.astype(str)+'_'+de.treatment
##de = de.set_index('Unique_Identifier')[['adj.P.Val','logFC','logFDR_de']]
##
##act = act.merge(de,how='left',left_index=True,right_index=True).fillna(0.)
##
##plt.scatter(act['logFDR_de'],act['logFDR'])
##plt.xlabel('Differential Expression FDR Corrected p-value (-log10)')
##plt.ylabel('TF Activity FDR Corrected p-value (-log10)')
##corr = stats.spearmanr(act['adj.P.Val'],act['Fisher_FDR'])
##plt.text(0.1,24,'Spearman: Correlation: '+str(corr[0]))
##plt.text(14.5,23,'p-value: '+str(corr[1]))
##plt.xlim([0,100])
##plt.ylim([0,25])
##plt.title('Activity vs. Differential Expression for  TFs')
##plt.show()
##
###correlate p values for significant TFs
##shad = shad.merge(de,how='left',left_index=True,right_index=True).fillna(0.)
##
##plt.scatter(shad['logFDR_de'],shad['logFDR'])
##plt.xlabel('Differential Expression FDR Corrected p-value (-log10)')
##plt.ylabel('TF Activity FDR Corrected p-value (-log10)')
##corr = stats.spearmanr(shad['adj.P.Val'],shad['Fisher_FDR'])
##plt.text(0.1,24,'Spearman: Correlation: '+str(corr[0]))
##plt.text(14.5,23,'p-value: '+str(corr[1]))
##plt.xlim([0,100])
##plt.ylim([0,25])
##plt.title('Activity vs. Differential Expression for Significant TFs')
##plt.show()
##
###histograms
##on = act[act.TF_state=='ON']
##off = act[act.TF_state=='OFF']
##un = act[act.TF_state=='unchanged']
##plt.hist(on['logFDR_de'].tolist(),bins=32,color='r',alpha=0.5)
##plt.hist(off['logFDR_de'].tolist(),bins=24,color='b',alpha=0.5)
##plt.hist(un['logFDR_de'].tolist(),bins=32,color='k',alpha=0.2)
##plt.legend(['enabled TFs','disabled TFs','unchanged TFs'])
###plt.ylim([0,500]) #remove this to get the top of the graph
##plt.xlabel('Differential Expression FDR Corrected p-value (-log10)')
##plt.ylabel('Frequency')
##plt.show()
##
##plt.hist(on['logFC'].tolist(),bins=15,color='r',alpha=0.5)
##plt.hist(off['logFC'].tolist(),color='b',alpha=0.5)
##plt.hist(un['logFC'].tolist(),bins = 20,color='k',alpha=0.2)
##plt.legend(['enabled TFs','disabled TFs','unchanged TFs'])
###plt.ylim([0,1100]) #remove this to get the top of the graph
##plt.xlabel('Log Fold Change of Expression')
##plt.ylabel('Frequency')
##plt.show()



