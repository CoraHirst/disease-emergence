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



#################### additive model ######################
### multiple-mutation probability of extinction function with additive R0

pEmergence_additive = function(m, mu, R0_1, R0_2, xstart, init){ #m is number of types, mu is mutation rate, R0_1 is wt R0, Ro_2 is adapted type R0, xstart is a vector of initial guesses for fixed point equations (length(m)), and init is initial number of cases of each type
   
  ### define system of nl equations
  single_mut <- function(x) {
    y <- numeric(m) #vector to store fixed points
    ### define fixed point equations
    for(i in 1:(m-1)){
      y[i] = exp(-(1-mu)*(R0_1+Rstep*(1-i))*(1-x[i]))*exp(-mu*R0_1*(1-x[i+1])) - x[i] #all variants from wt to just before adapted variant
    }
    y[ms[j]] = exp(-R0_2*(1-x[m])) - x[m] #adapted variant
    y
  }
  
  # function for prob.emergence
  prob_emergence = function(qs,init) {1 - prod(qs^init)} #qs are solutions to fixed point equation, init is vector of initial number of infections o each type  
  
  # parameters
  m = m #number of variant types (number of mutations required for emergence +1)
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
