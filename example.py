import pandas as pd 
import temporal_simulation as ts
import numpy.random as rdm


data='looped_primary_school'
original_df=pd.read_csv('data/'+data+'.txt', sep='\t', header=None, names=['ID1','ID2','start_time','end_time'])
# read it again with edge directions reversed (so that the disease can go in both directions)
reverse_df=pd.read_csv('data/'+data+'.txt', sep='\t', header=None, names=['ID2','ID1','start_time','end_time'])    
#put them together
df=pd.concat([original_df,reverse_df])    

ID_list=list(set(df['ID1']))     
contacts_dict={}
for node in ID_list:
    node_df=df[(df['ID1']==node)]
    names=node_df['ID2'].tolist()
    start_times=node_df['start_time'].tolist()
    end_times=node_df['end_time'].tolist()
    contacts_dict[node]=[]
    for i in range(len(names)):
        contacts_dict[node].append([names[i],start_times[i],end_times[i],0])

parameters={'beta':0.001,
            'l_mode':22, 
            'l_dispersion':1.1,
            'i_mode':2,
            'i_shape':5,
            'asymptomatic_proportion':0.0,
            }

seed=ID_list[int(len(ID_list)*rdm.random())]
time_of_infection=min([contact[1] for contact in contacts_dict[seed]])
tree=ts.get_infection_tree(seed,contacts_dict,time_of_infection,parameters)    

for t in tree:
    print(t[0],'infected by',t[2],'at time',t[1])

