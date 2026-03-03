################### probability of emergence with fixed increase in R0 each mutation ########################
pEmergence_fixed = function(m, mu, R0_1, R0_2, Rfinal) {
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

#### lets see if that works

########################### run sim ###################################
## parameters
R0_1s = seq(0, 1.2, by = 0.01) #wt R0s to test
R0_2 = 1.2 #R0 at which human-human transmission is achieved/end of interval
Rfinal = 1000
m = 3 #how many types? m-1 mutations required, meaning m-2 intermediate types
mu = 10^-3 #mutation rate

probs.emergence = NaN # vector to store probabilities of emergence

# solve for probabilities given each starting R_0
for(i in 1:length(R0_1s)){
  probs.emergence[i] = pEmergence_fixed(m = m, mu = mu, R0_1 = R0_1s[i], R0_2 = R0_2, Rfinal = Rfinal)
}

## format df for plotting
probs.emergence = cbind("R0" = R0_1s, "PE" = probs.emergence) %>% as.data.frame()

# lets compare this to the additive model where m = 3
#### additive model
probs.emergence.additive = NaN
# initial guess
xstart = c(1, rep(0, (m-1)))
#set init
init = c(1, rep(0, m-1))

# choose R0_1
for(i in 1:length(R0_1s)){
  # add solution to vector
  probs.emergence.additive[i] = pEmergence_additive(m=m, mu=mu, R0_1=R0_1s[i], R0_2=R0_2, xstart=xstart, init=init)
}

### prepare probs.emergence.additive for plotting
# format wide df
probs.emergence.additive = cbind("R0" = R0_1s, "PE" = probs.emergence.additive) %>% as.data.frame() #add R01_s column and format as df


### plot oh please work 
ggplot() + 
  geom_line(data = probs.emergence, aes(x = R0, y = PE, col = "stepwise_mutations")) +
  geom_line(data = probs.emergence.additive, aes(x = R0, y = PE, col = "additive")) +
  theme_bw() +
  labs(x = "R0 of introduced pathogen", y = "Probability of emergence", title = "3 types, 2 mutations required") +
  scale_y_continuous(trans = "log10") + 
  scale_color_manual(values = c("black", "forestgreen"), labels = c("additive","stepwise_mutations"))

