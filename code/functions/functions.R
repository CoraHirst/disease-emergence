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
