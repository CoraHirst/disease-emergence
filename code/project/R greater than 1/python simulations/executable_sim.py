#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 19 16:25:13 2025

@author: chirst
"""

### execute single_regime

import single_regime as sr
import sys


if __name__ == "__main__":
  Rwt = sys.argv[1]
  deltaR = sys.argv[2]
  Rfinal = sys.argv[3]
  mu = sys.argv[4]
  gens = sys.argv[5]
  spills = sys.argv[6]
  
  
Rwt = 1.3
deltaR = 0.5
Rfinal = 3
mu = 0.003
gens = 100
spills = 10 
 
outcome = sr.single_regime(Rwt=1.3, deltaR=0.5, Rfinal=3, mu=0.003, gens=100, spills=10)

print(outcome)

sr.single_regime(Rwt, deltaR, Rfinal, mu, gens, spills)
