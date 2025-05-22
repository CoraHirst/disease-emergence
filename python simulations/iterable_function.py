#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 18 12:46:14 2025

@author: chirst
this function defines a branching process model of adaptive evolution during infection chains that takes vector parm_vector as input and returns a matrix where columns are the number of infections with variants 1:m and rows are generations of infections
parm_vector takes in [Rwt, deltaR, Rfinal, mu, and gens] 
"""

#packages
import numpy as np
import math
import numpy.random as random



#iterable function 
def iterable_bp_deltaR(parm_vector): # model goes to evolution of high transmissibility or extinction after a single introduction and returns the outcome, final size, and total number of generations of the stochastic run
    """Branching process model of adaptive virus evolution to high transmissibility during transmission after spillover. Inputs initial transmissibility, final transmissibility, change in transmissibility with a mutation, mutation rate, and number of generations to run the branching process."""
    #define parameters 
    Rwt = parm_vector['Rwt']  #wt r
    deltaR = parm_vector['deltaR']#change in r with mut  
    Rfinal = parm_vector['Rfinal']#r of final type 
    mu = parm_vector['mu']#mutation rate
    gens = parm_vector['gens']#number of generations
    
    m = math.ceil((Rfinal - Rwt)/deltaR) + 1 # number of steps is the interval between final and initital divided by delta R. Number of types is number of steps + 1.
    R0 = Rwt + np.arange(0,m+1)*deltaR #+1 because np.arange is exclusive
    
    #empty data frame to store numbers of individuals of each type over gens generations
    X = np.full(shape = (m,gens+1), fill_value= 0) #m rows means 0 - (m-1) indexing because, python. #gens+1 columns means 0-(gens) indexing. 
    
    #initialize
    X[0,0] = 1 # start with only a single infection of the wild type

    ##### simulation
    ### stochastic simulation 
    for gen in range(1,gens):
      X[0,gen] = sum(random.poisson(lam = (1-mu)*R0[0], size = X[0,gen-1]))
      for i in range(1,m-2):
        X[i,gen] = sum(random.poisson(lam = (1-mu)*R0[i], size = X[i,gen-1])) + sum(random.poisson(lam = mu*R0[i-1], size = X[i-1,gen-1]))
      X[m-1,gen] = sum(random.poisson(lam = mu*R0[m-2], size = X[m-2,gen-1])) + sum(random.poisson(lam = R0[m-1], size = X[m-1, gen-1]))
     
      if gen == gens: #if we reach the end of sim before extinction or emergence, typ = 1 (we can discuss this later) 
        typ = 1
        size = np.sum(X)
        generation = gen
        break
      if all(X[:,gen] == 0): # if extinct
        typ = 0
        size = np.sum(X)
        generation = gen
        break
      elif any(X[m-1,:] > 0):
        typ = 1
        size = np.sum(X)
        generation = gen
        break
      elif np.sum(X) > 10**(8): #if too big
        typ = 1
        size = np.nan
        generation = gen
        break
      else:
          typ = 2
          size = np.sum(X)
          generation = gen
    iteration_result = [typ,size,generation]
    return iteration_result  # Optional return statement