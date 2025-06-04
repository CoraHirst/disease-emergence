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


#import iterable function 
import iterable_function as it

  
#define simulation
def single_regime(Rwt, deltaR, Rfinal, mu, gens, spills):
    
    parm_vector = {'Rwt': Rwt, 'deltaR': deltaR, 'Rfinal': Rfinal, 'mu': mu, 'gens': int(gens)} #input of iterable function
    
    # we would like to run these processes in parallel, as each iteration can be run independent of the outcome of the others. We want to return a vector coding the outcome, size, and generation of an iteration, and we want to output finally a matrix with rows being the output vector of each iteration
    
    if __name__ == '__main__':
        num_processes = int(spills) # Use all available cores
        with mp.Pool(processes=num_processes) as pool:
            results = pool.map(it.iterable_bp_deltaR, [parm_vector] * num_processes) # Execute my_function on each core, with the same argument
            pool.close()
            pool.join()
            # Process the results
    
            resultsDF = pd.DataFrame(np.stack(results), columns = ['outcome', 'total','generation'])
    
            info = {'nTooBig':len(resultsDF[resultsDF['total'] == np.nan]),'pUnfinished':(len(resultsDF[resultsDF['outcome'] == 2])/len(resultsDF)), 'pEmergent': resultsDF[resultsDF['outcome'] == 1]['outcome'].sum()/len(resultsDF),'pExtinct': (len(resultsDF[resultsDF['outcome'] == 0])/len(resultsDF)),'Emergent_size': resultsDF[resultsDF['outcome'] == 1]['total'].mean(skipna=True),'Extinct_size': resultsDF[resultsDF['outcome'] == 0]['total'].mean(skipna = True),'Emergent_generation': resultsDF[resultsDF['outcome'] == 1]['generation'].mean(skipna=True),'Extinct_generation': resultsDF[resultsDF['outcome'] == 0]['generation'].mean(skipna=True)}

            return info

##### simulation
outcomes = [] #empty matrix to store results 

# import parms 
parms_list = np.loadtxt("parms.txt", dtype = float)

# iterate through all parameter regimes
for regime in range(len(parms_list)):
    outcome = single_regime(*parms_list[regime])
    outcomes.append(outcome) 
    
    print("One done, on to next")




outcomes # lets see how we did 

        
# chill for tonight! But then lets get on with the rest of this. 
#### moving on 
#### we need to concatenate all dictionaries into a data frame with keys as the column titles 
#### then we need to add two more columns by column binding the ind_parms array so that we know what our parameters were!
#### ideally we can try to speed it up even more so that we can run more iterations per regime
#### ultimately I'd like to make it modular enough that we can run it all straight from command line. But this is great! 

#### still slow - three days for 200 sims vs 3.6 days for 68 though so there's that!

# save output as dataframe
outcomesDF = pd.DataFrame(outcomes)
#save parms as dataframe
parmsDF = pd.DataFrame(parms_list, columns = ['Rwt', 'deltaR', 'Rfinal', 'mu', 'allottedGens', 'spills'])


#combine columnwise outcomes and parms for saving to analyse in r
outputDF = pd.concat([parmsDF, outcomesDF], axis = 1)

#### save output dataframe to outcomes.csv
file_path = "output3.csv"
outputDF.to_csv(file_path, index = False) #index=false prevents the csv file from writing the line indices



# code for saving list of dictionaries to a text file
#with open(file_path, 'w') as file: 
  #  for item in outcomes:
     #   file.write(str(item) + '\n')
        
######## no i did a dumb thing I have to take the mean of emergence size IGNORING NAs



