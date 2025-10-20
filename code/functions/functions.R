######################### adaptive evolution function in terms of increase in R0 with mutation, delta R
######## evenly dispersed distribution described by generating functions for secondary cases
pEmergence_deltaR = function(delta_R, mu, R_wt, R_adapted) {
  ### calculate number of variants
  m = ceiling(abs(log10(R_wt + 10^-8))/log10(1+delta_R)) +1  #number of types in the branching process #note that we need a slight overshoot of R0 = 1 to get supercritical
  
  ### define R0 vector 
  R0 = R_wt*(1+delta_R)^(0:(m-1))  # add the changes from wtR0 to give R0 of each variant
  R0[m] = R_adapted #final R0 >> 1 to account for differences in emergence probabilities arising from final R0
  
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
  init[1] = 1 #except for 1 case of the starting type
  
  ### define initial guess for fixed point vector
  xstart = rep(0,m) #easiest to set all probs to 0
  xstart[1] = 1 #except for the prob of a linneage started by our 1 starting case type 1
  
  ### solve for fixed points
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x #solve system of nonlinear equations for fixed point
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) #solve for 1-sum(q_i^init(i)) - should be q_1^1 for a single starting case of type 1 
  
  return(prob.emergence) 
}


#nonlinear equation solver to solve for fixed point of the generating function for the single-type branching process with a supercritical R
## define branching process generating function for offspring distribution
pEvolution_singletype = function(Rwt) {
  bp_singletype = function(x) {y = exp(Rwt*(x-1)) - x} #fixed point solution will be this function set to 0
  
  ## define initial conditions 
  init = 1 #start with single spillover infection meep meep
  
  ## start fixed point guess
  xstart = 0.5 #50 percent will go extinct so lets see
  
  ## solve for non extinction probability (1-Pextinction)
    Rwt = Rwt #set wildtype R
    return(1 - nleqslv(xstart, bp_singletype, method="Newton", global="none", control=list(trace=1,stepmax=2))$x) #solve for nonextinction probability
  }
  

#fraction to vaccinate to reduce R to some threshold below R=1
frac_vacc = function(Rt,Reff) {1-(Rt/Reff)}

#return non-extinction probability of single-type bp
non_extinction = function(Rwt) {
  #single-type BP 
  bp_singletype = function(x) {y = exp(Rwt*(x-1)) - x}
  #set start guess
  xstart = 0.5
  #solve for nonextinction prob
  P = 1 - nleqslv(xstart, bp_singletype, method="Newton", global="none", control=list(trace=1,stepmax=2))$x
  return(P) #solve for nonextinction probability
}

######## evenly dispersed distribution described by generating functions for secondary cases
pEmergence_supercrit_deltaR = function(delta_R, mu, R_wt, R_adapted, supercritical_R = 3) {
  m = ceiling((log10(supercritical_R) - log10(R_wt))/log10(1+delta_R)) + 1 #number of types in the branching process #note that we need a slight overshoot of R0 = 1 to get supercritical
  
  ### define R0 vector 
  R0 = R_wt*(1+delta_R)^(0:(m-1))  # add the changes from wtR0 to give R0 of each variant
  R0[m] = R_adapted #final R0 >> 1 to account for differences in emergence probabilities arising from final R0
  
  

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
  init[1] = 1 #except for 1 case of the starting type
  
  ### define initial guess for fixed point vector
  xstart = rep(0,m) #easiest to set all probs to 0
  xstart[1] = 0 #except for the prob of a linneage started by our 1 starting case type 1
  
  ### solve for fixed points
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x #solve system of nonlinear equations for fixed point
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) #solve for 1-sum(q_i^init(i)) - should be q_1^1 for a single starting case of type 1 
  
  return(prob.emergence) 
}

#calculates the first supercritical R for a lineage starting with wildtype R_wt and increasing by delta_R percent with each mutation

supercrit_R = function(Rwt, deltaR){
  m = ceiling(abs(log10(Rwt + 10^-12))/log10(1+deltaR)) #number of types in the branching process #note that we need a slight overshoot of R0 = 1 to get supercritical
  ### define R0 vector 
  super_R = Rwt*(1+deltaR)^(m)
  return(super_R) # add the changes from wtR0 to give R0 of first supercritical variant
}

########################## total emergencce probability ########################## 
total_P = function(p, S, years) { 1 - (1-p)^(S*years)}
