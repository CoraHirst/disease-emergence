#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 13:41:45 2025

@author: chirst
Run this file to generate a text file of desired parameters delineated by spaces
"""

#required libraries and functions
import numpy as np

#arguments 

# define matrix of arguments 
spills = [1000] #number of spillovers 
Rfinal = [3] #r of final type 
mu = [10**(-3)] #mutation rate
gens = [100] #number of generations

#only change Rwt and deltaR with each run
Rwt = (np.linspace(1.0,1.5, num = 20 , endpoint = False)).tolist()
deltaR = [.02,.05,.1,.25, 1]       



### create a grid of all possible parameter combinations
ind_parms = np.array([(a,b,c,d,e,f) for a in Rwt for b in deltaR for c in Rfinal for d in mu for e in gens for f in spills]) #first column is Rwt, second column is deltaR

#save to txt file 
np.savetxt("parms.txt", ind_parms, delimiter = ' ', fmt = '%f')








