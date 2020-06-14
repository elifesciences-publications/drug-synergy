#synergistic/diff expressed gene sets, for each combo
#jul 8 2015

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn as venn
import csv

##set p value cutoff
pval = 1E-18

txs = ['M','T','W','TM','TW','MW']
times = ['3','6','9','12','24']

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

#TFs from TM not explained by diffexp or splic
tfdict = {}
neitherdict = {}
for tx in txs:
    df = tfdict[tx] = pd.read_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv')
    neitherdict[tx] = pd.DataFrame() 
    for time in times:
        timedf = df[df.time==int(time)]
        sets[tx+'_'+time+'_TF'] = set(timedf[timedf.Fisher_FDR < 0.05].regulator)
        timedf = timedf[~timedf.regulator.isin(sets[tx+'_'+time+'_splic'].union(sets[tx+'_'+time+'_diff']))]
        neitherdict[tx] = pd.concat([neitherdict[tx],timedf])
        
neitherdict['TM'].sort('Fisher_FDR').head(n=20).to_csv('top_TFs_notexplained.csv')






#diffexp vs splicing stackplots and data

for tx in txs:
    both= [0]
    splic = [0]
    diffexp = [0]
    neither = [0]
    
    for time in times: 
        diffexp.append(len(((sets[tx+'_'+time+'_diff']).difference(sets[tx+'_'+time+'_splic'])).intersection(sets[tx+'_'+time+'_TF'])))
        splic.append(len(((sets[tx+'_'+time+'_splic']).difference(sets[tx+'_'+time+'_diff'])).intersection(sets[tx+'_'+time+'_TF'])))
        both.append(len(((sets[tx+'_'+time+'_splic']).intersection(sets[tx+'_'+time+'_diff'])).intersection(sets[tx+'_'+time+'_TF'])))
        neither.append(len(((sets[tx+'_'+time+'_TF']).difference(sets[tx+'_'+time+'_splic'])).difference(sets[tx+'_'+time+'_diff'])))
        
        

    plt.figure()
    plt.stackplot([0,1,2,3,4,5],diffexp,both,splic,neither,linewidth = 0.1,colors=('r','m','b','gray'),alpha=0.5)
    plt.xticks([0,1,2,3,4,5],times)
    plt.ylabel('Genes')
    plt.xlabel('Time')
    plt.title('Differentially Active TFs in '+tx)
    plt.legend(['Differentially Expressed','Both','Differentially Spliced','Neither'],bbox_to_anchor=[0.37,0.99],fontsize='smaller')#
    plt.savefig('TFs_Genes_'+tx+'.png')


#pie chart

##set fonts
font = {'fontname':'Arial','size':19}



plt.figure(figsize=(6,6))
sections, texts = plt.pie([sum(diffexp),sum(both),sum(splic),sum(neither)],labels = ['Differentially\nExpressed',
                                      'Differentially\nExpressed\n& Spliced',
                                      '     Differentially\n     Spliced','Other'],
        colors=['b','g','y','gray'],startangle=90)

for i in sections:
    i.set_linewidth(0)
    i.set_alpha(0.7)

for i in texts:
    i.set_family('Arial')
    i.set_size(16)

#plt.title('Gene Status of Differentially Active Transcription Factors',**font)
plt.savefig('TF_pie.png',dpi=300)




