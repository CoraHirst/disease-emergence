######################### adaptive evolution function in terms of increase in R0 with mutation, sigma
######## evenly dispersed distribution described by generating functions for secondary cases

pEmergence = function(sigma, mu, Re_spill, Re_evolved = 3) {
  ntypes = ceiling((log10(Re_evolved) - log10(Re_spill))/log10(1+sigma)) + 1 #number of types in the branching process #note that we need a slight overshoot of R0 = 1 to get supercritical
  
  ### define R0 vector 
  Re = Re_spill*(1+sigma)^(0:(ntypes-1))  # add the changes from wtR0 to give R0 of each variant
  Re[ntypes] = 1000 #final R0 >> 1 to account for differences in emergence probabilities arising from final R0
  

  
  ### define probability of emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} # 1 - product of the probabilities of extinction for each starting lineage
  
  ### define fixed point equations from PGFs
  ## system of nl equations to solve
  multi_mut <- function(x) { # for poisson dist. the mean is R_0 and the generating function is exp(Re(s -1))
    y <- numeric(ntypes)
    for(i in 1:(ntypes-1)){
      y[i] = exp(-(1-mu)*(Re[i])*(1-x[i]))*exp(-mu*Re[i]*(1-x[i+1])) - x[i] # fixed point from probability generating function 
    }
    y[ntypes] = exp(-Re[ntypes]*(1-x[ntypes])) - x[ntypes] # can only make type ntypes, regardless of mutation. 
    y
  }
  
  ### define initial conditions
  init = rep(0, ntypes) #start the initial number of cases of each type to 0 
  init[1] = 1 #except for 1 case of the starting type
  
  ### define initial guess for fixed point vector
  xstart = rep(0,ntypes) #easiest to set all probs to 0
  xstart[1] = 0 #except for the prob of a linneage started by our 1 starting case type 1
  
  ### solve for fixed points
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x #solve system of nonlinear equations for fixed point
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) #solve for 1-sum(q_i^init(i)) - should be q_1^1 for a single starting case of type 1 
  
  return(prob.emergence) 
}

########################## total emergencce probability ########################## 
total_P = function(p, S, years) { 1 - (1-p)^(S*years)}
