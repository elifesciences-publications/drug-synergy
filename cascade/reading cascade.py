#reading cascade files

import pandas as pd
import csv

pd.set_option('display.max_columns', 500)

types = {'0128':'time t only','0313':'time t and t-1'}

times = {1:0,2:3,3:6,4:9,5:12,6:24}

dfs = {'time t only':pd.DataFrame(index=times.values(),columns=['top 101:','top 111:',
                                                                'top 011:','past 001:',
                                                                'mid 111:','mid 101-111-011:',
                                                                'mid 101-011','bot 001',
                                                                'past 001->001','unexplained'],
                                  dtype=object),
       'time t and t-1':pd.DataFrame(index=times.values(),columns=['top 101:','top 111:',
                                                                'top 011:','past 001:',
                                                                'mid 111:','mid 101-111-011:',
                                                                'mid 101-011','bot 001',
                                                                'past 001->001','unexplained'],
                                  dtype=object)}

for date in ['0128','0313']:
    print(types[date])
    df = dfs[types[date]]
    for time in list(range(1,7)):
        print(times[time])
        with open('20170313_mcf7_transcriptional_cascade_layers/cascade_layer_'+
                  str(time)+'_2017'+date+'.txt','r') as file:
            file = csv.reader(file,delimiter=' ')
            file= list(file)
            for line in file:
                line = [i for i in line if i != '']
                if line[0] in ['top','past','mid','bot',]:
                    if len(line)>2:
                        if line[2]=='->':
                            df.ix[times[time],line[0]+' '+line[1]+line[2]+line[3]] = line[4:]
                        else:
                            df.ix[times[time],line[0]+' '+line[1]] = line[2:]
                    else:
                        df.ix[times[time],line[0]+' '+line[1]] = line[2:]
                else:
                    df.ix[times[time],line[0]] = line[1:]   

diff = pd.DataFrame(index=times.values(),columns=['top 101:','top 111:',
                                                                'top 011:','past 001:',
                                                                'mid 111:','mid 101-111-011:',
                                                                'mid 101-011','bot 001',
                                                                'past 001->001','unexplained'],
                    dtype=object)
diffr = diff.copy()

for row in diff.index:
    for col in diff.columns:
            diff.ix[row,col] = set(dfs['time t and t-1'].ix[row,col]).difference(dfs['time t only'].ix[
                row,col])
print('****ADDED WHEN INCLUDING PREVIOUS TIME****')
print(diff)

for row in diffr.index:
    for col in diffr.columns:
            diffr.ix[row,col] = set(dfs['time t only'].ix[row,col]).difference(dfs['time t and t-1'].ix[
                row,col])
print('****ADDED WHEN EXCLUDING PREVIOUS TIME****')
print(diffr)

dfs['N time t and t-1'] = dfs['time t and t-1'].applymap(len)
dfs['N time t only'] = dfs['time t only'].applymap(len)

for frame in ['N time t and t-1','N time t only']:
    df = dfs[frame]
    df['explained']=df[['mid 111:','mid 101-111-011:','mid 101-011','bot 001','past 001->001']].sum(axis=1)
    df['total'] = df[['explained','unexplained']].sum(axis=1)
    df['percent explained'] = 100*df['explained']/df['total']
    dfs[frame] = df

print('****NUMBERS****')
print(pd.concat([dfs['N time t only'],dfs['N time t and t-1']],axis=1,keys=['N time t only','N time t and t-1']))

Ndiff = dfs['N time t and t-1']-dfs['N time t only']

print('****DIFFERENCE****')
print(Ndiff)


