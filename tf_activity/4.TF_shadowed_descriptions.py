#make figures describing the altered state TFs
#sep 1 2015

import pandas as pd
from matplotlib import pyplot as plt
import matplotlib_venn as venn
import itertools
from matplotlib.colors import ColorConverter
import numpy as np
from scipy import stats


#access color dict for color mixing
cl = ColorConverter.colors

#using matplotlib-venn colormixing similar to venn.__venn3.compute_colors
def mix(col1,col2,col3=None):
    #input a color as string
    if col3==None:
        return venn._common.mix_colors(np.array(cl[col1]),np.array(cl[col2]))
    else:
        return venn._common.mix_colors(np.array(cl[col1]),np.array(cl[col2]),np.array(cl[col3]))


##set fonts
font = {'fontname':'Arial', 'size':19}


txs = ['T','M','TM','W','MW','TW']

combs = ['TM','MW','TW']

txtypes = ['Transcriptional Activator','Transcriptional Repressor']

states = ['ON','OFF']

times = [0,3,6,9,12,24]

minlist = 3

minov = 2

setsmins = {}

sets = {}

#Making figures for all treatments, all types
##with minimums
##c = []
c ={
'T' : [i/255. for i in [86, 112, 255]],
'M' : [i/255. for i in [255, 9, 60]],
'W': 'y',
#'W' : [i/255. for i in [249, 255, 60]],
'TM' : [i/255. for i in [255,0,255]],
'TW' : [i/255. for i in [0,255,79]],
'MW' : [i/255. for i in [255,170,0]]}


#find correlations
for comb in combs:
    print(comb)
    for time in times:
        print(time)
        two = pd.DataFrame()
        for tx in [comb[0],comb[1]]:
            df = pd.read_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv')
            df['idx'] = df['regulator']+'_'+df['regulation_type']
    
            timedf = df[df.time==time]
            timedf.drop_duplicates('idx',inplace=True)
            timedf.set_index('idx',inplace=True)

        
            #print(dfs[comb[0]+'_'+str(time)].Fisher_FDR.head())
            #print(dfs[comb[1]+'_'+str(time)].Fisher_FDR.head())
            two = two.merge(timedf[['Fisher_FDR']],how='outer',left_index=True,right_index=True)
        two.dropna(inplace=True)
        #print(two.head())
        print(stats.spearmanr(two.as_matrix())[0])

##        
##for i in range(12):
##    c.append(plt.cm.jet(22*i))
##colors = itertools.cycle(c)

l = ['-','--',':','-.']
lines = itertools.cycle(l)

for comb in combs:
    fig = plt.figure(combs.index(comb)+1)
    #ax = fig.add_subplot(111)

    for tx in [comb[0],comb[1],comb]:
        df = pd.read_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv')

        
        df = df[df.Fisher_FDR < 0.05]
        for txtype in txtypes:
            typedf = df[df.regulation_type==txtype]
            for state in states:
                statedf = typedf[typedf.TF_state==state]
                ls = []
                for time in times:
                    timedf = statedf[statedf.time==time]
                    ls.append(timedf.regulator.unique().shape[0])
                    setsmins[tx+'_'+txtype+'_'+state+'_'+str(time)]=set(timedf.regulator)
##                ax.plot(times,ls,color = c[tx],linestyle=next(lines),label=state+' '+txtype+'s '+' in '+tx)
##
##    handles, labels = ax.get_legend_handles_labels()
##    labels = [label.replace('ON','Activated').replace('OFF','Inactivated').replace(
##        'Transcriptional Repressors ','Negative Effectors').replace('Transcriptional Activators ','Positive Effectors') for label in labels]
##    print(labels)
##    print(type(labels))
##    lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.36,0.9))
##    for text in lgd.get_texts():
##        text.set_family('Arial')
##    ax.set_xticks(times)
##    ax.set_xticklabels(times,**font)
##    ax.set_xlim([0,24])
##    ax.set_xlabel('Time (hours)',**font)
##    ax.set_yticklabels([int(tick) for tick in ax.get_yticks()],**font)
##    print(ax.get_yticks())
##    #plt.legend(bbox_to_anchor=(1,0.5),loc = 'center left')
##    ax.set_ylabel('Number of Transcription Factors',**font)
##    #ax.set_title("TFs by Fisher's Exact Test with Minimum Criteria",**font)
##    #plt.legend(bbox_to_anchor=(1,0.5),loc = 'center left')
##    fig.savefig('shadowed_results/'+comb+'_regulators_fishers_shadowed'
##                ,bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)
##
for tx in txs:
    for time in times:
        setsmins[tx+'_'+str(time)] = setsmins[tx+'_Transcriptional Activator_ON_'+str(time)].union(
                                     setsmins[tx+'_Transcriptional Activator_OFF_'+str(time)]).union(
                                     setsmins[tx+'_Transcriptional Repressor_ON_'+str(time)]).union(
                                     setsmins[tx+'_Transcriptional Repressor_OFF_'+str(time)])
        for state in states:
            setsmins[tx+'_'+state+'_'+str(time)] = setsmins[tx+'_Transcriptional Activator_'+state+'_'+
                                                           str(time)].union(
                                     setsmins[tx+'_Transcriptional Repressor_'+state+'_'+str(time)])
        
##Make venn diagrams for each timepoint, for all differential TFs
##with minimums
##            ##figure 5
for comb in ['TM','MW','TW']:
    if comb == 'TM':
        labels = ('Tamoxifen','Mefloquine','TM')
    elif comb == 'MW':
        labels = ('Mefloquine','    Withaferin','MW      ')
    elif comb == 'TW':
        labels = ('Tamoxifen','Withaferin      ','           TW')
    for time in times:
        plt.figure()
        v = venn.venn3([setsmins[comb[0]+'_'+str(time)],
                        setsmins[comb[1]+'_'+str(time)],
                        setsmins[comb+'_'+str(time)],
                       ],
                       set_labels = labels)#['','','']) ##turn labels on here
        for i in ['001','100','010','110','011','101','111']:
            if v.get_label_by_id(i) != None:
                v.get_label_by_id(i).set_text('')
        for text in v.set_labels:
            text.set_fontsize(46)
            text.set_family('Arial')
        #plt.title('Differentially Active TFs at time '+str(time))
        #plt.show()
        plt.savefig('shadowed_results/'+comb+'_'+str(time)+'.pdf') #withnumbers
#####figure 6
##    if comb == 'TM':
##        labels = ('Tamoxifen','TM','Mefloquine')
##    elif comb == 'MW':
##        labels = ('Mefloquine','MW  ','  Withaferin')
##    elif comb == 'TW':
##        labels = ('Tamoxifen','TW     ' ,'    Withaferin')
##    for time in times:
##        plt.figure(time+100)
##        v = venn.venn3([setsmins[comb[0]+'_'+str(time)],
##                        setsmins[comb+'_'+str(time)],
##                        setsmins[comb[1]+'_'+str(time)],
##                        
##                       ],
##                       set_labels = ['','','']) ##turn labels on here labels)#
##        for i in ['001','100','010','110','011','101','111']:
##            if v.get_label_by_id(i) != None:
##                v.get_label_by_id(i).set_text('')
##        if v.get_patch_by_id('010') != None:
##            v.get_patch_by_id('010').set_color('white')
##            v.get_patch_by_id('010').set_edgecolor('k')
##        if v.get_patch_by_id('110') != None:
##            v.get_patch_by_id('110').set_color([i/255. for i in[86, 112, 255]])
##            v.get_patch_by_id('110').set_alpha(0.9)
##            #v.get_patch_by_id('110').set_edgecolor('k')
##        if v.get_patch_by_id('011') != None:
##            v.get_patch_by_id('011').set_color([i/255. for i in[255, 9, 60]])
##            v.get_patch_by_id('011').set_alpha(0.9)
##            #v.get_patch_by_id('110').set_edgecolor('k')
##        if v.get_patch_by_id('111') != None:
##            v.get_patch_by_id('111').set_color([i/255. for i in[255,0,255]])
##            v.get_patch_by_id('111').set_alpha(0.9)
##            #v.get_patch_by_id('110').set_edgecolor('k')
##        if v.get_patch_by_id('100') != None:
##            v.get_patch_by_id('100').set_color([i/255. for i in[86, 112, 255]])
##            #v.get_patch_by_id('110').set_alpha(1.)
##            #v.get_patch_by_id('110').set_edgecolor('k')
##        if v.get_patch_by_id('001') != None:
##            v.get_patch_by_id('001').set_color([i/255. for i in[255, 9, 60]])
##            #v.get_patch_by_id('011').set_alpha(1.)
##            #v.get_patch_by_id('110').set_edgecolor('k')
##        if v.get_patch_by_id('101') != None:
##            v.get_patch_by_id('101').set_color([i/255. for i in[255,0,255]])
##            #v.get_patch_by_id('111').set_alpha(1.)
##            #v.get_patch_by_id('110').set_edgecolor('k')
##        for text in v.set_labels:
##            text.set_fontsize(30)
##            text.set_family('Arial')
##        #plt.title('Differentially Active TFs at time '+str(time))
##        #plt.show()
##        plt.savefig('shadowed_results/'+comb+'_'+str(time)+'nolabels.pdf')
##        
##        for state in states:
##            plt.figure()
##            v = venn.venn3([setsmins[comb[0]+'_'+state+'_'+str(time)],
##                       setsmins[comb[1]+'_'+state+'_'+str(time)],
##                       setsmins[comb+'_'+state+'_'+str(time)]],
##                       set_labels = labels)
##            plt.title(state+' TFs at time '+str(time))
##            plt.savefig('shadowed_results/'+comb+'/'+comb+'_'+state+'_'+str(time)+'.pdf')

        
##get all possible sets for stackplot
tm = [0]
tw = [0]
mw = [0]

T1 = [0]
M1 = [0]
TM = [0]
TMnT = [0]
TMnM = [0]
TnM = [0]
TnMnTM = [0]

W2= [0]
T2 = [0]
TW = [0]
TWnW = [0]
TWnT = [0]
TnW = [0]
TWnTnW = [0]

W3= [0]
M3 = [0]
MW = [0]
MWnM = [0]
MWnW = [0]
MnW = [0]
MWnMnW = [0]

sets = setsmins

for hour in times[1:]:
    time = str(hour)
    #TM
    T1.append(len(((sets['T_'+time]).difference(sets['M_'+time])).difference(
        sets['TM_'+time])))
    M1.append(len(((sets['M_'+time]).difference(sets['T_'+time])).difference(
        sets['TM_'+time])))
    TnM.append(len(((sets['T_'+time]).intersection(sets['M_'+time])).difference(
        sets['TM_'+time])))
    TnMnTM.append(len(((sets['T_'+time]).intersection(sets['M_'+time])).intersection(
        sets['TM_'+time])))
    TMnT.append(len(((sets['T_'+time]).intersection(sets['TM_'+time])).difference(
        sets['M_'+time])))
    TMnM.append(len(((sets['M_'+time]).intersection(sets['TM_'+time])).difference(
        sets['T_'+time])))
    TM.append(len((sets['TM_'+time].difference(sets['T_'+time])).difference(
        sets['M_'+time])))
    #TW
    W2.append(len(((sets['W_'+time]).difference(sets['T_'+time])).difference(
        sets['TW_'+time])))
    T2.append(len(((sets['T_'+time]).difference(sets['W_'+time])).difference(
        sets['TW_'+time])))
    TnW.append(len(((sets['T_'+time]).intersection(sets['W_'+time])).difference(
        sets['TW_'+time])))
    TWnTnW.append(len(((sets['T_'+time]).intersection(sets['W_'+time])).intersection(
        sets['TW_'+time])))
    TWnW.append(len(((sets['TW_'+time]).intersection(sets['W_'+time])).difference(
        sets['T_'+time])))
    TWnT.append(len(((sets['TW_'+time]).intersection(sets['T_'+time])).difference(
        sets['W_'+time])))
    TW.append(len((sets['TW_'+time].difference(sets['T_'+time])).difference(
        sets['W_'+time])))
    #MW
    W3.append(len(((sets['W_'+time]).difference(sets['M_'+time])).difference(
        sets['MW_'+time])))
    M3.append(len(((sets['M_'+time]).difference(sets['W_'+time])).difference(
        sets['MW_'+time])))
    MnW.append(len(((sets['W_'+time]).intersection(sets['M_'+time])).difference(
        sets['MW_'+time])))
    MWnMnW.append(len(((sets['W_'+time]).intersection(sets['M_'+time])).intersection(
        sets['MW_'+time])))
    MWnW.append(len(((sets['MW_'+time]).intersection(sets['W_'+time])).difference(
        sets['M_'+time])))
    MWnM.append(len(((sets['MW_'+time]).intersection(sets['M_'+time])).difference(
        sets['W_'+time])))
    MW.append(len((sets['MW_'+time].difference(sets['M_'+time])).difference(
        sets['W_'+time])))

plt.figure()
plt.stackplot([0,3,6,9,12,24],T1,TnM,M1,TMnT,TnMnTM,TMnM,TM,linewidth=0,
              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
plt.xticks([0,3,6,9,12,24],times,**font)
plt.ylabel('Differentially Active Transciption Factors',**font)
plt.xlabel('Time (hours)',**font)
#plt.title('Differentially Active Transcription Factors in Combination TM',**font)
#plt.legend(['Tamoxifen','Tamoxifen-Mefloquine','Mefloquine','Tamoxifen-TM','Mefloquine-TM','Intersection','TM'])#,bbox_to_anchor=[0.5,1] 
plt.xlim([0,24])
plt.ylim([0,250])
plt.yticks(**font)
plt.savefig('TFs_TM_stackplot.pdf',dpi=300)

plt.figure()
plt.stackplot([0,3,6,9,12,24],T2,TnW,W2,TWnT,TWnTnW,TWnW,TW,linewidth=0,
              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
plt.xticks([0,3,6,9,12,24],times,**font)
plt.ylabel('Differentially Active Transciption Factors',**font)
plt.xlabel('Time(hours)',**font)
#plt.title('Differentially Active Transcription Factors in Combination TW',**font)
#plt.legend(['Withaferin','Intersection','Withaferin-TW','TW'],bbox_to_anchor=[0.5,1])
plt.xlim([0,24])
plt.ylim([0,250])
plt.yticks(**font)
plt.savefig('TFs_TW_stackplot.pdf',dpi=300)

plt.figure()
plt.stackplot([0,3,6,9,12,24],M3,MnW,W3,MWnM,MWnMnW,MWnW,MW,linewidth=0,
              colors=('r',mix('r','g'),'g',mix('r','b'),mix('r','b','g'),mix('b','g'),'b'),alpha=0.4)
plt.xticks([0,3,6,9,12,24],times,**font)
plt.ylabel('Differentially Active Transciption Factors',**font)
plt.xlabel('Time (hours)',**font)
#plt.title('Differentially Active Transcription Factors in Combination MW',**font)
#plt.legend(['Withaferin','Intersection','Withaferin-MW','MW'],bbox_to_anchor=[0.5,1])
plt.xlim([0,24])
plt.ylim([0,250])
plt.yticks(**font)
plt.savefig('TFs_MW_stackplot.pdf',dpi=300)
##
####concordant and discordant ON/OFF TFs
##l = ['-','--',':']
##lines = itertools.cycle(l)
##for comb in combs: 
##    fig = plt.figure(4+combs.index(comb))
##    ax = fig.add_subplot(111)
##    for tx in [comb[0],comb[1],comb]:
##        df = pd.read_csv('MCF7_aracne_'+tx+'_Fishers_shadowed_TFs.csv')
##        df = df[df.TF_state!='unchanged']
##        concordant = pd.DataFrame()
##        discordant = pd.DataFrame(columns=['Fisher_FDR', 'Fisher_p_value', 'direction', 'regulation_type',
##           'regulator', 'regulome_size', 'time', 'TF_state'])
##        nonsensical = pd.DataFrame(columns = ['Fisher_FDR', 'Fisher_p_value', 'direction', 'regulation_type',
##           'regulator', 'regulome_size', 'time', 'TF_state'])
##        unique = pd.DataFrame()
##        extras = pd.DataFrame()
##        for time in times:
##            timedf = df[df.time==time]
##            dups = timedf[timedf.duplicated('regulator')]
##            duplicated = timedf[timedf.regulator.isin(dups.regulator.tolist())]
##            for regulator in duplicated.regulator.unique().tolist():
##                regdf = duplicated[duplicated.regulator==regulator]
##                if len(regdf.TF_state.unique()) > 1:
##                    if len(regdf.regulation_type.unique())==1:
##                        nonsensical = pd.concat([nonsensical,regdf])
##                    elif len(regdf.regulation_type.unique())==2:
##                        if regdf.shape[0]==2:
##                            discordant = pd.concat([discordant,regdf])
##                        else:
##                            one = regdf[regdf.regulation_type==regdf.regulation_type.unique()[0]]
##                            if one.shape[0] > 1:
##                                nonsensical = pd.concat([nonsensical,one])
##                            else:
##                                unique = pd.concat([unique,one])
##                            two = regdf[regdf.regulation_type==regdf.regulation_type.unique()[1]]
##                            if two.shape[0] > 1:
##                                nonsensical = pd.concat([nonsensical,two])
##                            else:
##                                unique = pd.concat([unique,two])
##                else:
##                    concordant = pd.concat([concordant,regdf])
##            unique = pd.concat([unique,timedf[~timedf.regulator.isin(dups.regulator.tolist())]])
####        concordant.to_csv('shadowed_results/'+comb+'/'+'MCF7_aracne_'+tx+'_Fishers_shadowed_TFs_'+'concordant.csv',index=False)
####        discordant.to_csv('shadowed_results/'+comb+'/'+'MCF7_aracne_'+tx+'_Fishers_shadowed_TFs_'+'discordant.csv',index=False)
####        nonsensical.to_csv('shadowed_results/'+comb+'/'+'MCF7_aracne_'+tx+'_Fishers_shadowed_TFs_'+'nonsensical.csv',index=False)
####        unique.to_csv('shadowed_results/'+comb+'/'+'MCF7_aracne_'+tx+'_Fishers_shadowed_TFs_'+'unique.csv',index=False)
##            if nonsensical.shape[0] > 0:
##                print(nonsensical.head())
##    ##plot
##        unls = []
##        conls = []
##        disls = []
##        nonls = []
##        for time in times:
##            unls.append(len(unique[unique.time==time].regulator.unique()))
##            conls.append(len(concordant[concordant.time==time].regulator.unique()))
##            if discordant.shape[0]>0:
##                disls.append(len(discordant[discordant.time==time].regulator.unique()))
##            else:
##                disls.append(0)
##            nonls.append(len(nonsensical[nonsensical.time==time].regulator.unique()))
##        ax.plot(times,unls,color = c[tx],linestyle=next(lines),label='unique tfs in '+tx)
##        ax.plot(times,conls,color = c[tx],linestyle=next(lines),label='concordant tfs in '+tx)
##        ax.plot(times,disls,color = c[tx],linestyle=next(lines),label='discordant tfs in '+tx)
## #       ax.plot(times,nonls,color = c[tx],linestyle=next(lines),label='nonsensical tfs in '+tx)
##    handles, labels = ax.get_legend_handles_labels()
##    lgd = ax.legend(handles, labels, loc='upper center', bbox_to_anchor=(1.23,0.9))
##    for text in lgd.get_texts():
##        text.set_family('Arial')
##    ax.set_xticks(times)
##    ax.set_xticklabels(times,**font)
##    ax.set_xlim([0,24])
##    ax.set_yticklabels([int(tick) for tick in ax.get_yticks()],**font)
##    ax.set_xlabel('Time (hours)',**font)
##    ax.set_ylabel('Number of Transcription Factors',**font)
##    #ax.set_title("TFs by Fisher's Exact Test with Minimum Criteria",**font)
##    fig.savefig('shadowed_results/'+tx+'/'+tx+'_Fishers_shadowed_TFs_concordance'
##                ,bbox_extra_artists=(lgd,), bbox_inches='tight',dpi=300)
