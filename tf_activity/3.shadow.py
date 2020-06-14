#shadow analysis
#for tfs at tx, time:
#compare the regulomes of each 2 tfs
#if significant overlap:
#calculate size of overlap and store
#order by p-value of overlap,
#for each pair

#sept 30 2015

import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
from statsmodels.sandbox.stats.multicomp import multipletests as mult
from datetime import datetime
startTime = datetime.now()

#####functions
def sign(row,hour):
    if row[hour] >= 0:
        return 'ON'
    if row[hour] < 0:
        return 'OFF'

def direction(row):
    if row.TF_state == 'ON':
        if row.regulation_type == 'Transcriptional Activator':
            return 'up'
        if row.regulation_type == 'Transcriptional Repressor':
            return 'down'
    if row.TF_state == 'OFF':
        if row.regulation_type == 'Transcriptional Activator':
            return 'down'
        if row.regulation_type == 'Transcriptional Repressor':
            return 'up'

def Fishers_mins_w_ovsize(row,minlist,minov):
    txtype1 = row.regulation_type1
    txtype2 = row.regulation_type2
    regulator1 = row.regulator1
    regulator2 = row.regulator2
    
    r1 = set(reg[(reg.Regulator==regulator1)&(
        reg.Type==txtype1)].Target)
    notr1 = set(reg[(reg.Regulator!=regulator1)|(
        reg.Type!=txtype1)].Target)
    r2 = set(reg[(reg.Regulator==regulator2)&(
        reg.Type==txtype2)].Target)
    notr2 = set(reg[(reg.Regulator!=regulator2)|(
        reg.Type!=txtype2)].Target)
    ##contingency table: 
            ## [[r1&r2,diffex&notr1&r2,
            ## [r1&notr2, notr1&notr2]]
    overlap = r1.intersection(r2)
    if len(r1)>=minlist and len(r2)>=minlist and len(overlap)>=minov:
        cont = [[len(overlap),len(r1.intersection(notr2))], \
                [len(notr1.intersection(r2)),len(notr1.intersection(notr2))]]
                
        return (stats.fisher_exact(cont)[1], len(overlap))
    else:
        return (np.nan, len(overlap))


def Fishers_shadow_choose(row,minlist,minov):
    
    if row.direction1 ==row.direction2:
        direction = row.direction1
    else:
        return False

    time = row.time
    
    txtype1 = row.regulation_type1
    txtype2 = row.regulation_type2
    regulator1 = row.regulator1
    regulator2 = row.regulator2
    
    r1 = set(reg[(reg.Regulator==regulator1)&(
        reg.Type==txtype1)].Target)
    r2 = set(reg[(reg.Regulator==regulator2)&(
        reg.Type==txtype2)].Target)
    
    r1r2intersection = r1.intersection(r2)
    left = r1 - r2
    right = r2 - r1
    
    notr1r2inter = inter-r1r2intersection
    notleft = inter-left
    notright = inter-right

    dx = set(de[(de.direction==direction)&(de.time==time)].gene)
    notdx = inter - dx
    
    ##contingency table: 
            ## [[diffex&inregulome,diffex&not in regulome],
            ## [not diffex&in regulome, not diffex & not in regulome]]
    
    #left regulome
    #print('left',len(dx),len(left))
    if len(dx)>=minlist and len(left)>=minlist and len(dx.intersection(left))>=minov:
        cont = [[len(dx.intersection(left)),len(dx.intersection(notleft))], \
                [len(notdx.intersection(left)),len(notdx.intersection(notleft))]]
        leftp = stats.fisher_exact(cont)[1]
    else:
        leftp = np.nan

    #right regulome
    #print('right',len(dx),len(right))
    if len(dx)>=minlist and len(right)>=minlist and len(dx.intersection(right))>=minov:
        cont = [[len(dx.intersection(right)),len(dx.intersection(notright))], \
                [len(notdx.intersection(right)),len(notdx.intersection(notright))]]
        rightp = stats.fisher_exact(cont)[1]
    else:
        rightp = np.nan

    if len(dx)>=minlist and len(r1r2intersection)>=minlist and len(dx.intersection(r1r2intersection))>=minov:
        cont = [[len(dx.intersection(r1r2intersection)),len(dx.intersection(notr1r2inter))], \
                [len(notdx.intersection(r1r2intersection)),len(notdx.intersection(notr1r2inter))]]
        interp = stats.fisher_exact(cont)[1]
    else:
        interp = np.nan

    return (leftp,rightp,interp)


#
#####code

#set treatment to analyze
tx = 'W'

##set regulome and dysregulated genes list minimum length
minlist = 3
##set minimum overlap 
minov = 2

#import differentially expressed
de = pd.read_csv('../gene_expression/diff_limma_december/diffexp_'+
                 tx+'_alltimes_JELD.csv')

#import network and limit to genes for which we have data
reg = pd.read_csv('../network_new/MCF7_Woo_et_al-Li_et_al_genelevel_noregulatorsastargets_20160126_JELD.txt',
                  sep='\t')
allRNAseq = set(pd.read_csv('../gene_expression/diff_limma_december/limma_t0_M_0.txt',header=0,
                            sep='\t',index_col='1').index)
regset = set(reg.Regulator).union(set(reg.Target))
inter = regset.intersection(allRNAseq)
reg = reg[(reg.Regulator.isin(inter))&(reg.Target.isin(inter))]

de = de[de.gene.isin(inter)]

#run the tests for each direction, timepoint
dirs = ['up','down']
times = [0,3,6,9,12,24]#

#fishers exact test data - throw out nonsensical
tfs = pd.concat([pd.read_csv(
    'MCF7_aracne_'+tx+'_regulators_fishers_mins_3_2concordant.csv'),
                 pd.read_csv(
                     'MCF7_aracne_'+tx+'_regulators_fishers_mins_3_2unique.csv'),
                 pd.read_csv('MCF7_aracne_'+tx+'_regulators_fishers_mins_3_2discordant.csv')])

overlaps = pd.DataFrame(columns = ['regulator1','direction1','regulation_type1','TF_state1',
                                           'regulator2','direction2','regulation_type2','TF_state2',
                                           'time','shadow_p_value','shadow_overlap'])

##check the overlap of each TF regulome with each other TF regulome
for time in times:
    print('TIME: ',time)
    tf = tfs[tfs.time==time]
    for di in dirs:
        erlist = []
        print('running Fishers exact tests for '+di+' ...')
        s = tf[tf.direction==di]
        error = False
        for i, row in s.iterrows():
            df = s[s.regulator==row.regulator]
            if df.shape[0] > 1:
                er = df.regulator.reset_index(drop=True)[0]
                if er not in erlist:
                    print('Too many rows for this TF',er)
                    erlist.append(er)
        #print('getting combinations...')
        x = s.T.to_dict()
        temp = pd.DataFrame(list(combinations(x.values(),2)))
        if temp.shape[0]>0:
            #print('getting left and right...')
            l = pd.DataFrame(temp[0].tolist())
            r = pd.DataFrame(temp[1].tolist())
            newcols = [col+'1' for col in l.columns if col != 'time'] +['time']
            l.columns = newcols
            newcols = [col+'2' for col in r.columns if col != 'time'] + ['time']
            r.columns = newcols
            r.drop('time',axis=1,inplace=True)
            #print('concatting...')
            new = pd.concat([l,r],axis=1)
            temp = pd.DataFrame(
                new.apply(Fishers_mins_w_ovsize,axis=1,minlist=minlist,minov=minov).tolist())
            temp.columns = ['shadow_p_value','shadow_overlap']
            new = pd.concat([new,temp],axis=1)
            overlaps = pd.concat([overlaps,new])

## of altered TF in Fishers
print('FISHERS:')
for time in times:
    timedf = tfs[tfs.time==time]
    print(str(len(timedf.regulator.unique()))+
    ' altered tfs at time '+str(time))
    
#fdr only for the tests that were actually done
print('Calculating FDR...')
nan = overlaps[np.isnan(overlaps.shadow_p_value)]
notnan = overlaps[~np.isnan(overlaps.shadow_p_value)]
notnan['shadow_FDR'] = mult(notnan['shadow_p_value'],method='fdr_bh')[1]
overlaps = pd.concat([nan,notnan])
overlaps.shadow_FDR = overlaps.shadow_FDR.astype(float)

overlaps.to_csv('TFs_shadow_'+tx+'.csv',index=False)

##for TF pairs that overlap, choose which TF is the better one
print('Calculating new Fishers exact tests with differential expression...')
df = overlaps[overlaps.shadow_FDR<0.05]
new = pd.DataFrame(df.apply(Fishers_shadow_choose,axis =1,minlist = minlist, minov= minov).tolist())
new.columns = ['p_value_1','p_value_2','p_value_intersection']
df = pd.concat([df.reset_index(drop=True),new],axis=1)

#combine the two new columns of p values, take fdr, then put the fdrs back in the dataframe
print('Calculating FDR...')
ser1= df.p_value_1
ser2 = df.p_value_2
ser3 = df.p_value_intersection
ser2.index=ser2.index+ser1.shape[0]
ser3.index=ser3.index+ser1.shape[0]+ser2.shape[0]
p = pd.DataFrame(pd.concat([ser1,ser2,ser3]).dropna())
p['FDR1'] = mult(p[0],method='fdr_bh')[1]
df= df.merge(p[['FDR1']],how='left',left_index=True,right_index=True)
p.rename(columns = {'FDR1':'FDR2'},inplace=True)
p.index=p.index-ser1.shape[0]
df= df.merge(p[['FDR2']],how='left',left_index=True,right_index=True)
p.rename(columns = {'FDR2':'FDRi'},inplace=True)
p.index=p.index-ser2.shape[0]
df = df.merge(p[['FDRi']],how='left',left_index=True,right_index=True)

#choose the cases where either left or right is significant, not both
answer = df[((df.FDR1<0.05)&((df.FDR2>=0.05)|(np.isnan(df.FDR2))))|(
             ((df.FDR1>=0.05)|(np.isnan(df.FDR1)))&(df.FDR2<0.05))]
#cases where none of left, right, or intersection are significant and to be removed
bad = df[((df.FDR1>=0.05)|(np.isnan(df.FDR1)))&((df.FDR2>=0.05)|(np.isnan(df.FDR2)))&((df.FDRi>=0.05)|(np.isnan(df.FDRi)))]
answer[['FDR1','FDR2']] = answer[['FDR1','FDR2']].fillna(1)
answer['todrop'] = answer[['FDR1','FDR2']].T.idxmax()
badcopy = bad.copy()
bad['todrop'] = 'FDR1'
badcopy['todrop'] = 'FDR2'

answer = pd.concat([answer,bad,badcopy])
answer.to_csv(tx+'_choosing_after_shadow.csv')

print('Getting final TFs to drop...')
todropdf = pd.DataFrame()
for i, row in answer.iterrows():
    ls = [name for name in row.index if row.todrop[-1] in name] +['time']
    temp = pd.DataFrame(row[ls]).T
    temp.columns = [col.strip('_12') for col in ls]
    todropdf = pd.concat([todropdf,temp])

print('Dropping...')
print(tfs.shape)
tfs = tfs[~(tfs.TF_state+'_'+tfs.regulation_type+'_'+tfs.regulator+'_'+tfs.time.astype(int).astype(
    str)).isin(
    (todropdf.TF_state+'_'+todropdf.regulation_type+'_'+todropdf.regulator+'_'+todropdf.time.astype(
        int).astype(str)).tolist())]
print(tfs.shape)
tfs.to_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv',index=False)

print('FINAL AFTER SHADOW:')
for time in times:
    timedf = tfs[tfs.time==time]
    print(str(len(timedf.regulator.unique()))+
    ' altered tfs at time '+str(time))
    
print(datetime.now() - startTime)


##MEFLOQUINE
##FISHERS:
##0 altered tfs at time 0
##23 altered tfs at time 3
##49 altered tfs at time 6
##50 altered tfs at time 9
##52 altered tfs at time 12
##49 altered tfs at time 24

##FINAL AFTER SHADOW:
##0 altered tfs at time 0
##4 altered tfs at time 3
##30 altered tfs at time 6
##28 altered tfs at time 9
##38 altered tfs at time 12
##33 altered tfs at time 24


##TAMOXIFEN
##FISHERS:
##0 altered tfs at time 0
##18 altered tfs at time 3
##32 altered tfs at time 6
##26 altered tfs at time 9
##49 altered tfs at time 12
##73 altered tfs at time 24
##
##FINAL AFTER SHADOW:
##0 altered tfs at time 0
##4 altered tfs at time 3
##9 altered tfs at time 6
##13 altered tfs at time 9
##27 altered tfs at time 12
##56 altered tfs at time 24


##TM
##FISHERS:
##0 altered tfs at time 0
##119 altered tfs at time 3
##234 altered tfs at time 6
##262 altered tfs at time 9
##338 altered tfs at time 12
##333 altered tfs at time 24
##
##FINAL AFTER SHADOW:
##0 altered tfs at time 0
##63 altered tfs at time 3
##143 altered tfs at time 6
##167 altered tfs at time 9
##190 altered tfs at time 12
##175 altered tfs at time 24


##W
##FISHERS:
##0 altered tfs at time 0
##96 altered tfs at time 3
##172 altered tfs at time 6
##173 altered tfs at time 9
##182 altered tfs at time 12
##166 altered tfs at time 24
##
##FINAL AFTER SHADOW:
##0 altered tfs at time 0
##61 altered tfs at time 3
##135 altered tfs at time 6
##143 altered tfs at time 9
##145 altered tfs at time 12
##139 altered tfs at time 24


##TW
##FISHERS:
##0 altered tfs at time 0
##95 altered tfs at time 3
##151 altered tfs at time 6
##164 altered tfs at time 9
##156 altered tfs at time 12
##157 altered tfs at time 24
##
##FINAL AFTER SHADOW:
##0 altered tfs at time 0
##65 altered tfs at time 3
##122 altered tfs at time 6
##136 altered tfs at time 9
##133 altered tfs at time 12
##134 altered tfs at time 24


##MW
##FISHERS:
##0 altered tfs at time 0
##93 altered tfs at time 3
##151 altered tfs at time 6
##160 altered tfs at time 9
##170 altered tfs at time 12
##150 altered tfs at time 24
##
##FINAL AFTER SHADOW:
##0 altered tfs at time 0
##55 altered tfs at time 3
##124 altered tfs at time 6
##137 altered tfs at time 9
##138 altered tfs at time 12
##127 altered tfs at time 24
