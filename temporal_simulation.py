import numpy as np
import random

def get_infection_tree(seed,contacts,time,params):
    #print()
    susceptible_nodes=list(set(contacts.keys()))
    N=len(susceptible_nodes)
    
    time_of_infection={}
        
    source_of_infection={seed:None}
    #make a proportion of the nodes asymptomatic
    asymptomatic_nodes=random.sample(susceptible_nodes,int(N*params['asymptomatic_proportion']))

    #infection_tree is the tree network of infections that the function returns
    infection_tree=[]    
    #keep a list of infected nodes (at the beginning it contains only 'node')    
    infections=[[seed,time]]
    
    #get the parameters for the lognormal distribution
    sigma=np.log(params['l_dispersion'])
    mu=(sigma**2)+np.log(params['l_mode'])

    while len(infections) > 0:
        
        #get earliest infection event on the list
        next_infection=min(infections,key=lambda x: x[1])        
        #remove the chosen infection
        infections.remove(next_infection)
        
        infectious_node=next_infection[0]
        #add it to the immune list
        susceptible_nodes.remove(infectious_node)     
        #update the time
        time=next_infection[1]
        
        #add to the tree
        infection_tree.append(next_infection+[source_of_infection[infectious_node]])
        
        #randomly select a latent duration 
        if infectious_node==seed:
            latent_duration=0
        else:
            latent_duration=int(60*60*np.random.lognormal(mu,sigma))
            
        #randomly select an infectious duration for each type
        infectious_duration={}
        for i in range(2):
            if infectious_node in asymptomatic_nodes:
                infectious_duration=24*60*60
            else:
                infectious_duration=int(60*60*np.random.gamma(params['i_shape'], scale=params['i_mode']/(params['i_shape']-1)))
                
        #contact_list=contacts[infectious_node].copy()#
        contact_list=[contact for contact in contacts[infectious_node] if contact[1]<time+latent_duration+infectious_duration and contact[2]>time+latent_duration]
        #loop over all the potentially infectious interactions
        while len(contact_list) > 0:
            contact=contact_list.pop()
            
            name=contact[0]
            contact_start=contact[1]
            contact_end=contact[2]
                            
            #exposure starts either at the start of infectious period or the start of interaction, whichever is later
            exposure_start=max(time+latent_duration,contact_start)
            #exposure endes at the end of infectious period or the end of interaction, whichever is earlier
            exposure_end=min(time+latent_duration+infectious_duration,contact_end)                          
            
            #select a random time for the infection to occur
            r=random.random()
            #this is equivalent to an attempt at transmission occuring each second (different beta depending on the location)           
            l=params['beta']/(1-params['beta'])            
            #l=-np.log(1-b)                              
            infection_time=exposure_start+int(-(1/l)*np.log(1-r))             

            #check that the interaction is with a susceptible node, that it is not a sick day at work                                  
            if infection_time<exposure_end and name in susceptible_nodes:
                #if this is the case then a potential infection occurs!
                
                #if 'name' is already in the list of infections then check to see which happened first
                if name in [i[0] for i in infections]:
                    if infection_time<time_of_infection[name]:
                        #remove the later infection time
                        infections.remove([name,time_of_infection[name]])
                        #add the earlier infection time
                        infections.append([name,infection_time])
                        #keep track of where it came from
                        source_of_infection[name]=infectious_node
                        #update the list
                        time_of_infection[name]=infection_time
                #if this is the first time 'name' has been infected then add it to the list
                else:
                    #keep track of the time of the infection
                    time_of_infection[name]=infection_time
                    #keep track of who it came from
                    source_of_infection[name]=infectious_node
                    #update the list
                    infections.append([name,infection_time])

    return infection_tree
    
