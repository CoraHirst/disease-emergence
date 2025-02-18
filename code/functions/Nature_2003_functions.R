####### Figure 2: Single mutation required for emergence #######
### single-mutation probability of extinction function 
pEmergence_single = function(mu, R0_1, Rfinal){
  ### ## Define fixed point equations
  # pgfs - qi
  fpeq_1 = function(q1,q2) {exp(-(1-mu)*R0_1*(1-q1))*exp(-mu*R0_1*(1-q2)) - q1} # pgf of secondary cases from type 1 set to and evaluated at q1, q2 (extinction probs for lineages starting with type 1 or type 2 infection)
  fpeq_2 = function(q1,q2) {exp(-Rfinal*(1-q2)) - q2} # pgf of secondary cases from type 2 set to and evaluated at q1, q2
  
  # define system of nl equations
  single_mut <- function(x) { # input function for nleqslv equation solver
    y <- numeric(2)
    y[1] <- fpeq_1(x[1], x[2]) #first element is extinction probability for starting with infection of type 1
    y[2] <- fpeq_2(x[1], x[2]) #second element is extinction probability for starting with infection of type 2
    y
  }
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} # defined as 1- extinction probability; for m-type branching processes, this is the product of the probabilities of extinction of a lineage starting with each type raised to the number of those types the process starts with. In other words, its the total probability that at least one lineage the process starts with does not go extinct 
  
  # parameters
  mu = mu #mutation rate
  R0_1 = R0_1 # R_0 of type 1 infections, aka R_0 of introduced pathogen
  Rfinal = Rfinal # R_0 of type 2 infections
  
  # define initial conditions 
  init = c(1,0) #start with 1 wild type and no evolved
  
  # initial guess for extinction probabilities for lineages starting with type 1 and type 2 
  xstart = c(0.99, 0.01) #start with a guess of 0.99 extinction starting with lineage 1 and 0 for type 2
  
  ### calculate fixed points and emergence prob
  # newton start
  qs = nleqslv(xstart, single_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x # solve for extinction probabilities from nleqslv package - see documentation
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) # because our model starts with 1 infection of type 1, this simplifies to q1^1. But its good to be explicit
  
  return(prob.emergence)
}

###### Figure 3: Multiple mutations required for emergence ####
#### jackpot model
pEmergence_jackpot = function(n, mu, R0_1, Rfinal){ #n is number of intermediate types, m-2, mu is mutation rate, R0_1 is wt R0, Rfinal is adapted type R0
 
  # parms 
  R0_1 = R0_1 #R0_1 of initial wildtype
  Rfinal = Rfinal #final adapted type R0_m
  mu = mu #mutation rate (same for each type)
  m = n+2 #number of types here, number of intermediates in paper
  
  # initial conditions
  init = c(1, rep(0, m-1)) #vector initial numbers of infections of each type, length(m); start with only 1 infection of wildtype, type 1
  xstart = c(1, rep(0, m-1)) #vector of initial guesses (between 0 and 1) for fixed point solution, (length(m))
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} #defined as 1- probability of extinction, which is defined as the product of the probabilities that the lineage from each starting type goes extinct
  
  ## solve system of equations
  # define system of nl equations
  multi_mut <- function(x) {
    y <- numeric(m)
    for(i in 1:(m-1)){ #vectorize, fixed point equations i:(m-1) will have the same form
      y[i] = exp(-(1-mu)*(R0_1)*(1-x[i]))*exp(-mu*R0_1*(1-x[i+1])) - x[i] #type i only gives rise to type i or type i+1
    }
    y[m] = exp(-Rfinal*(1-x[m])) - x[m] #type m only gives rise to type m
    y
  }
  
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x #extinction probabilities
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) #calculate emergence prob from extinction prob
  
  return(prob.emergence)
}

#### additive model
pEmergence_additive = function(n, mu, R0_1, Rfinal){ #m is number of types, mu is mutation rate, R0_1 is wt R0, Rfinal is adapted type R0
  
  # parms 
  R0_1 = R0_1 #R0_1 of initial wildtype
  Rfinal = Rfinal #final adapted type R0_m
  mu = mu #mutation rate (same for each type)
  m = n+2 #number of types here, number of intermediates in paper
  Rstep = (1-R0_1)/(m-1) # increase in R0 with each step
  
  # initial conditions
  init = c(1, rep(0, m-1)) #vector initial numbers of infections of each type, length(m); start with only 1 infection of wildtype, type 1
  xstart = c(1, rep(0, m-1)) #vector of initial guesses (between 0 and 1) for fixed point solution, (length(m))
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)}
  
  ## solve system of equations
  # define system of nl equations
  multi_mut <- function(x) {
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0_1+(Rstep*(i-1)))*(1-x[i]))*exp(-mu*(R0_1+(Rstep*(i-1)))*(1-x[i+1])) - x[i]
    }
    y[m] = exp(-Rfinal*(1-x[m])) - x[m]
    y
  }
  
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x #extinction probabilities
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) #calculate emergence prob from extinction prob
  
  return(prob.emergence)
}

#### additive model
pEmergence_additive = function(n, mu, R0_1, Rfinal){ #m is number of types, mu is mutation rate, R0_1 is wt R0, Rfinal is adapted type R0
  
  # parms 
  R0_1 = R0_1 #R0_1 of initial wildtype
  Rfinal = Rfinal #final adapted type R0_m
  mu = mu #mutation rate (same for each type)
  m = n+2 #number of types here, number of intermediates in paper
  Rstep = (1-R0_1)/(m-1) # increase in R0 with each step
  
  # initial conditions
  init = c(1, rep(0, m-1)) #vector initial numbers of infections of each type, length(m); start with only 1 infection of wildtype, type 1
  xstart = c(1, rep(0, m-1)) #vector of initial guesses (between 0 and 1) for fixed point solution, (length(m))
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)}
  
  ## solve system of equations
  # define system of nl equations
  multi_mut <- function(x) {
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0_1+(Rstep*(i-1)))*(1-x[i]))*exp(-mu*(R0_1+(Rstep*(i-1)))*(1-x[i+1])) - x[i]
    }
    y[m] = exp(-Rfinal*(1-x[m])) - x[m]
    y
  }
  
  # newton start
  qs = nleqslv(xstart, multi_mut, method="Newton", global="none", control=list(trace=1,stepmax=2))$x #extinction probabilities
  
  # solve for prob.emergence
  prob.emergence = prob_emergence(qs = qs, init = init) #calculate emergence prob from extinction prob
  
  return(prob.emergence)
}


