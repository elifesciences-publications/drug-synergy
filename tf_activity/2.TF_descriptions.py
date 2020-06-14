#make figures describing the altered state TFs
#sep 1 2015

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn as venn
import itertools
from scipy import stats

txs = ['T','M','TM','W','MW','TW']

txtypes = ['Transcriptional Activator','Transcriptional Repressor']

states = ['ON','OFF']

times = [0,3,6,9,12,24]

minlist = 3

minov = 2

setsmins = {}

sets = {}

dfs= {}

#Making figures for all treatments, all types
##with minimums
c = []
for i in range(12):
    c.append(plt.cm.jet(22*i))
colors = itertools.cycle(c)
##
##fig = plt.figure(1)
##ax = fig.add_subplot(111)
##
for tx in txs:
    df = pd.read_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'.csv')



##    df = df[df.Fisher_FDR < 0.05]
##    for txtype in txtypes:
##        typedf = df[df.regulation_type==txtype]
##        for state in states:
##            statedf = typedf[typedf.TF_state==state]
####            ls = []
##            for time in times:
##                timedf = statedf[statedf.time==time]
####                ls.append(timedf.regulator.unique().shape[0])
##                setsmins[tx+'_'+txtype+'_'+state+'_'+str(time)]=set(timedf.regulator)
####            ax.plot(ls,color = next(colors),label=state+' '+txtype+'s '+' in '+tx)



##handles, labels = ax.get_legend_handles_labels()
##lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.4,0.9))
##ax.set_xticks([0,1,2,3,4,5])
##ax.set_xticklabels(times)
##ax.set_xlabel('Time')
###plt.legend(bbox_to_anchor=(1,0.5),loc = 'center left')
##ax.set_ylabel('Number of Transcription Factors')
##ax.set_title("TFs by Fisher's Exact Test with Minimum Criteria")
###plt.legend(bbox_to_anchor=(1,0.5),loc = 'center left')
##fig.savefig(tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)
##            ,bbox_extra_artists=(lgd,), bbox_inches='tight')
##            
##
####Make venn diagrams for each timepoint
####with minimums
##for time in times:
##    for txtype in txtypes:
##        for state in states:
##            plt.figure()
##            v = venn.venn3([setsmins['M_'+txtype+'_'+state+'_'+str(time)],
##                           setsmins['T_'+txtype+'_'+state+'_'+str(time)],
##                           setsmins['TM_'+txtype+'_'+state+'_'+str(time)]],
##                           set_labels = ('Mefloquine','Tamoxifen','TM'))
##            plt.title(state+' '+txtype+'s at time '+str(time))
##            plt.savefig('venn_mins/'+state+' '+txtype+'s at time '+str(time))
    

nonsens = set()
##concordant and discordant ON/OFF TFs
fig = plt.figure(2)
ax = fig.add_subplot(111)
for tx in ['T','M','TM']:
    df = pd.read_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'.csv')
    df = df[df.TF_state!='unchanged']
    concordant = pd.DataFrame()
    discordant = pd.DataFrame(columns=['Fisher_FDR', 'Fisher_p_value', 'direction', 'regulation_type',
       'regulator', 'regulome_size', 'time', 'TF_state'])
    nonsensical = pd.DataFrame(columns = ['Fisher_FDR', 'Fisher_p_value', 'direction', 'regulation_type',
       'regulator', 'regulome_size', 'time', 'TF_state'])
    unique = pd.DataFrame()
    extras = pd.DataFrame()
    for time in times:
        timedf = df[df.time==time]
        dups = timedf[timedf.duplicated('regulator')]
        duplicated = timedf[timedf.regulator.isin(dups.regulator.tolist())]
        for regulator in duplicated.regulator.unique().tolist():
            regdf = duplicated[duplicated.regulator==regulator]
            if len(regdf.TF_state.unique()) > 1:
                if len(regdf.regulation_type.unique())==1:
                    nonsensical = pd.concat([nonsensical,regdf])
                elif len(regdf.regulation_type.unique())==2:
                    if regdf.shape[0]==2:
                        discordant = pd.concat([discordant,regdf])
                    else:
                        one = regdf[regdf.regulation_type==regdf.regulation_type.unique()[0]]
                        if one.shape[0] > 1:
                            nonsensical = pd.concat([nonsensical,one])
                        else:
                            unique = pd.concat([unique,one])
                        two = regdf[regdf.regulation_type==regdf.regulation_type.unique()[1]]
                        if two.shape[0] > 1:
                            nonsensical = pd.concat([nonsensical,two])
                        else:
                            unique = pd.concat([unique,two])
            else:
                concordant = pd.concat([concordant,regdf])
        unique = pd.concat([unique,timedf[~timedf.regulator.isin(dups.regulator.tolist())]])
    nonsens = nonsens.union(set(nonsensical.regulator))
####    concordant.to_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'concordant.csv',index=False)
####    discordant.to_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'discordant.csv',index=False)
####    nonsensical.to_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'nonsensical.csv',index=False)
####    unique.to_csv('MCF7_aracne_'+tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'unique.csv',index=False)       
####plot
##    unls = []
##    conls = []
##    disls = []
##    nonls = []
##    for time in times:
##        unls.append(len(unique[unique.time==time].regulator.unique()))
##        conls.append(len(concordant[concordant.time==time].regulator.unique()))
##        if discordant.shape[0]>0:
##            disls.append(len(discordant[discordant.time==time].regulator.unique()))
##        nonls.append(len(nonsensical[nonsensical.time==time].regulator.unique()))
##    ax.plot(unls,color = next(colors),label='unique tfs in '+tx)
##    ax.plot(conls,color = next(colors),label='concordant tfs in '+tx)
##    ax.plot(disls,color = next(colors),label='discordant tfs in '+tx)
##    ax.plot(nonls,color = next(colors),label='nonsensical tfs in '+tx)
##handles, labels = ax.get_legend_handles_labels()
##lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.4,0.9))
##ax.set_xticks([0,1,2,3,4,5],times)
##ax.set_xlabel('Time')
##ax.set_ylabel('Number of Transcription Factors')
##ax.set_title("TFs by Fisher's Exact Test with Minimum Criteria")
##fig.savefig(tx+'_regulators_fishers_'+'mins_'+str(minlist)+'_'+str(minov)+'_concordance'
##            ,bbox_extra_artists=(lgd,), bbox_inches='tight')
##fig.show()
