######################### functions for multi-type branching process with even dispersion ##########################

#################### fixed-step model ######################
### Additive (fixed-step) R0 model where wtR0 determines how far from final type; final-type R0 >> 1
################### probability of emergence with fixed increase in R0 each mutation ########################

pEmergence_fixed_evenlydispersed = function(m, mu, R0_1, R0_2, Rfinal) {
  ### calculate step size and breaks
  Rstep = R0_2/(m-1) #increase in R0 with each mutation
  breaks = Rstep*(1:(m)) #intervals describing the range of R0 for types (variants) with increasing number of mutations
  
  ### define R0 vector 
  ## determine into which interval the wtR0 falls
  diffs = breaks - R0_1 #distance from wtR0 to each interval break
  type = which(diffs == min(diffs[which(diffs >= 0)]))[1] # position of the least positive distance gives the interval - type - into which the wtR0 falls
  
  ## calculate change in R0 from wt to each other variant
  delta_R0 = Rstep*c(-(type-1):0, 1:(m-type)) #vector of changes from wtR0 to the R0 of each other variant
  R0 = R0_1 + delta_R0 # add the changes from wtR0 to wtR0 to give R0 of each variant
  R0[m] = Rfinal #final R0 >> 1 to account for differences in emergence probabilities arising from final R0
  
  ### define probability of emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} # 1 - product of the probabilities of extinction for each starting lineage
  
  ### define fixed point equations from PGFs
  ## system of nl equations to solve
  multi_mut <- function(x) { # for poisson dist. the mean is R_0 and the generating function is exp(R0(s -1))
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0[i])*(1-x[i]))*exp(-mu*R0[i]*(1-x[i+1])) - x[i] # fixed point from probability generating function 
    }
    y[m] = exp(-R0[m]*(1-x[m])) - x[m] # can only make type m, regardless of mutation. 
    y
  }
  
  ### define initial conditions
  init = rep(0, m) #start the initial number of cases of each type to 0 
  init[type] = 1 #except for 1 case of the starting type
  
  ### define initial guess for fixed point vector
  xstart = rep(0,m) #easiest to set all probs to 0
  xstart[type] = 1 #except for the prob of a linneage started by our 1 starting case type
  
  ### solve for fixed points
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init)
  
  return(prob.emergence) 
}
