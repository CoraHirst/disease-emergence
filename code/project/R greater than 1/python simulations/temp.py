# -*- coding: utf-8 -*-
"""
Spyder Editor

This function runs spills iterations of iterable function and returns 6 features calculated from all runs: the probability of emergence, the probability of extinction, the size if emergent, the size if extinct, the generation of emergence, and the generation of extinction
these features are saved in a dictionary which is the output of this file. 
"""

#packages
import numpy as np
import pandas as pd
import multiprocessing as mp

#import iterable funciont
import iterable_function

# Parms should be defined when running script 
#spills = 100 #number of spillovers 
#Rwt = 1.3 #wt r
#deltaR = 0.3 #change in r with mut  
#Rfinal = 3 #r of final type 
#mu = 10**(-3) #mutation rate
#gens = 100 #number of generations

#parm vector
def single_regime(Rwt, deltaR, Rfinal, mu, gens, spills):
    
    parm_vector = {'Rwt': Rwt, 'deltaR': deltaR, 'Rfinal': Rfinal, 'mu': mu, 'gens': gens}
    
    
    # we would like to run these processes in parallel, as each iteration can be run independent of the outcome of the others. We want to return a vector coding the outcome, size, and generation of an iteration, and we want to output finally a matrix with rows being the output vector of each iteration
    
    if __name__ == '__main__':
        num_processes = spills # Use all available cores
        with mp.Pool(processes=num_processes) as pool:
            results = pool.map(iterable_function.iterable_bp_deltaR, [parm_vector] * num_processes) # Execute my_function on each core, with the same argument
        # Process the results
    
    resultsDF = pd.DataFrame(np.stack(results), columns = ['outcome', 'total','generation'])
    
    
    info = {'pEmergent': resultsDF[resultsDF['outcome'] == 1]['outcome'].sum()/len(resultsDF),
            'pExtinct': 1 - ((resultsDF[resultsDF['outcome'] == 1])['outcome'].sum()/len(resultsDF)),
            'Emergent_size': resultsDF[resultsDF['outcome'] == 1]['total'].mean(),
            'Extinct_size': resultsDF[resultsDF['outcome'] == 0]['total'].mean(),
            'Emergent_generation': resultsDF[resultsDF['outcome'] == 1]['generation'].mean(),
            'Extinct_generation': resultsDF[resultsDF['outcome'] == 0]['generation'].mean()}
    
    
    return info