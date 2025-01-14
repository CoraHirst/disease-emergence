### single-mutation probability of extinction function 
pEmergence_single = function(mu, R0_1, R0_2, xstart, init){
  ### ## Define fixed point equations
  # pgfs - qi
  fpeq_1 = function(q1,q2) {exp(-(1-mu)*R0_1*(1-q1))*exp(-mu*R0_1*(1-q2)) - q1}
  fpeq_2 = function(q1,q2) {exp(-R0_2*(1-q2)) - q2}

  # define system of nl equations
  single_mut <- function(x) {
    y <- numeric(2)
    y[1] <- fpeq_1(x[1], x[2])
    y[2] <- fpeq_2(x[1], x[2])
    y
  }

  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} #we could vectorize this for qs longer than 

  # parameters
  mu = mu #mutation rate
  R0_1 = R0_1 #initial R0
  R0_2 = R0_2 #final R0_2

  # define initial conditions 
  init = init #start with 1 wild type and no evolved

  # initial guess
  xstart = xstart
  
  ### calculate fixed points and emergence prob
  # newton start
  qs = nleqslv(xstart, single_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) 
  
  return(prob.emergence)
}

###### Multiple mutations required
#### jackpot model
pEmergence_jackpot = function(m, mu, R0_1, R0_2, xstart, init){ #m is number of types, mu is mutation rate, R0_1 is wt R0, Ro_2 is adapted type R0, xstart is a vector of initial guesses for fixed point equations (length(m)), and init is initial number of cases of each type
  ## parms 
  R0_1 = R0_1 #R0_1 of initial wildtype
  R0_2 = R0_2 #final adapted type R0_2
  mu = mu #mutation rate (same for each type)
  m = m #number of types here, number of intermediates in paper
  init = init #vector initial numbers of infections of each type
  xstart = xstart #vector of initial guesses (between 0 and 1) for fixed point solution
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)}
  
  
  ## solve system of equations
  # define system of nl equations
  multi_mut <- function(x) {
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0_1)*(1-x[i]))*exp(-mu*R0_1*(1-x[i+1])) - x[i]
    }
    y[m] = exp(-R0_2*(1-x[m])) - x[m]
    y
  }
  
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init)
  
  return(prob.emergence)
}

###### Multiple mutations required
#### additive model
pEmergence_additive_paper = function(mu, R0_1, R0_2, Rstep, xstart, init){ #m is number of types, mu is mutation rate, R0_1 is wt R0, Ro_2 is adapted type R0, xstart is a vector of initial guesses for fixed point equations (length(m)), and init is initial number of cases of each type
  ## parms 
  R0_1 = R0_1 #R0_1 of initial wildtype
  R0_2 = R0_2 #final adapted type R0_2
  R0_int = (R0_1+R0_2)/2 #intermediate R0
  mu = mu #mutation rate (same for each type)

  init = init #vector initial numbers of infections of each type
  xstart = xstart #vector of initial guesses (between 0 and 1) for fixed point solution
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)}
  
  
  ## solve system of equations
  # define system of nl equations
  multi_mut <- function(x) {
    y <- numeric(m)
    y[1] = exp(-(1-mu)*R0_1*(1-x[1]))*exp(-mu*R0_1*(1-x[2])) - x[1]
    y[2] = exp(-(1-mu)*R0_int*(1-x[2]))*exp(-mu*R0_int*(1-x[3])) - x[2]
    y[3] = exp(-R0_2*(1-x[3])) - x[3]
    y
  }
  
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init)
  
  return(prob.emergence)
}


###### Multiple mutations required
#### additive model
pEmergence_additive = function(m, mu, R0_1, R0_2, xstart, init){ #m is number of types, mu is mutation rate, R0_1 is wt R0, Ro_2 is adapted type R0, xstart is a vector of initial guesses for fixed point equations (length(m)), and init is initial number of cases of each type
  ## parms 
  R0_1 = R0_1 #R0_1 of initial wildtype
  R0_2 = R0_2 #final adapted type R0_2
  mu = mu #mutation rate (same for each type)
  m = m #number of types here, number of intermediates in paper
  Rstep = (R0_2-R0_1)/(m-1) # increase in R0 with each step
  init = init #vector initial numbers of infections of each type
  xstart = xstart #vector of initial guesses (between 0 and 1) for fixed point solution
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)}
  
  ## solve system of equations
  # define system of nl equations
  multi_mut <- function(x) {
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0_1+(Rstep*(i-1)))*(1-x[i]))*exp(-mu*(R0_1+(Rstep*(i-1)))*(1-x[i+1])) - x[i]
    }
    y[m] = exp(-R0_2*(1-x[m])) - x[m]
    y
  }
  
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init)
  
  return(prob.emergence)
}


#################### additive model ######################
### Additive (fixed-step) R0 model where wtR0 determines how far from final type
################### probability of emergence with fixed increase in R0 each mutation ########################
pEmergence_fixed = function(m, mu, R0_1, R0_2) {
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
  
  ### define probability of emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} # 1 - product of the probabilities of extinction for each starting linneage
  
  ### define fixed point equations from PGFs
  ## system of nl equations to solve
  multi_mut <- function(x) {
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0[i])*(1-x[i]))*exp(-mu*R0[i]*(1-x[i+1])) - x[i]
    }
    y[m] = exp(-R0[m]*(1-x[m])) - x[m]
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


