######################### functions for multi-type branching process with over dispersion ##########################


#################### fixed-step model ######################
### Additive (fixed-step) R0 model where wtR0 determines how far from final type; final-type R0 >> 1
################### probability of emergence with fixed increase in R0 each mutation ########################

pEmergence_fixed_overdispersed = function(m, mu, R0_1, R0_2, Rfinal, r) { #r is overdispersion param for neg binom, limit r-> of neg binom is poisson dist. 
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
  
  ### calculate probabilities of success for neg binoms
  p1 = r/((1-mu)*R0 + r) #given mean of neg binom is (1-mu)R0
  p2 = r/((mu)*R0 + r) #given mean of neg binom is (mu)R0
  
  pfinal = r/(R0 + r) # given mean of neg binom is R0
  
  ### define fixed point equations from PGFs
  ## system of nl equations to solve
  multi_mut <- function(x) { # for poisson dist. the mean is R_0 and the generating function is exp(R0(s -1))
    y <- numeric(m)
    for(i in 1:(m-1)){
      y[i] = ((p1[i]/(1-(1-p1[i])*x[i]))^r)*((p2[i]/(1-(1-p2[i])*x[i+1]))^r) # fixed point from probability generating function 
    }
    y[m] = ((p1[m]/(1-(1-p1[m])*x[m]))^r)# can only make type m, regardless of mutation. 
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


#### lets test this out 


## parameters
R0_1s = seq(0, 1, by = 0.01) #wt R0s to test
R0_2 = 1 #R0 at which human-human transmission is achieved/end of interval
Rfinal = 1000 #R0 of final phenotype
ms = c(2,5,10) #how many types? m-1 mutations required, meaning m-2 intermediate types
mu = 10^-3 #mutation rate

## matrix to store emergence probabilities 
probs.emergence = matrix(nrow = length(R0_1s), ncol = length(ms)) # vector to store probabilities of emergence

## simulate
for(j in 1:length(ms)){
  for(i in 1:length(R0_1s)){
    probs.emergence[i,j] = pEmergence_fixed(m = ms[j], mu = mu, R0_1 = R0_1s[i], R0_2 = R0_2, Rfinal = Rfinal) #run simulation for each m,R0_1 combination
  }
}

## format matrix for plotting
probs.emergence = cbind(R0_1s, probs.emergence) %>% as.data.frame() #format as dataframe with column denoting the R0_1 used to generate the probabilities row-wise and column for each m used to generate the probabilities column-wise
colnames(probs.emergence) = c("R0", ms) #name the columns for easy referencing
probs.emergence = probs.emergence %>% #pivot to long-form
  pivot_longer(cols = !R0, #keep R0_1 column
               names_to = "m", #add feature denoting "m" used row-wise
               values_to = "P.Emergence") #add feature denoting the pemergence calculated from m in m column and R0_1 in R0_1 column

## plot

# colors for plotting
#my.colors = c("black", "red", "forestgreen", "blue", "orange", "purple") #because pretty
# plot
ggplot() + 
  geom_line(data = probs.emergence, aes(x = R0, y = P.Emergence, group = m, col = m)) + #plot p emergence as function of R0_1 for each m
  theme_bw() + #white background with grey gridlines
  labs(x = "R0 of introduced pathogen", y = "Probability of Emergence", title = "probability of emergence with m phenotypes") + #axis and title labels
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(trans = "log10") + #log10 y scale
  scale_color_manual(values = my.colors) #because pretty


### lets try the neg binom with really large r

## parameters
R0_1s = seq(0, 1, by = 0.01) #wt R0s to test
R0_2 = 1 #R0 at which human-human transmission is achieved/end of interval
Rfinal = 1000 #R0 of final phenotype
ms = c(2,5,10) #how many types? m-1 mutations required, meaning m-2 intermediate types
mu = 10^-3 #mutation rate
r = 5 #overdispersion parameter

## matrix to store emergence probabilities 
probs.emergence = matrix(nrow = length(R0_1s), ncol = length(ms)) # vector to store probabilities of emergence

## simulate
for(j in 1:length(ms)){
  for(i in 1:length(R0_1s)){
    probs.emergence[i,j] = pEmergence_fixed_overdispersed(m = ms[j], mu = mu, R0_1 = R0_1s[i], R0_2 = R0_2, Rfinal = Rfinal, r = r) #run simulation for each m,R0_1 combination
  }
}

## format matrix for plotting
probs.emergence = cbind(R0_1s, probs.emergence) %>% as.data.frame() #format as dataframe with column denoting the R0_1 used to generate the probabilities row-wise and column for each m used to generate the probabilities column-wise
colnames(probs.emergence) = c("R0", ms) #name the columns for easy referencing
probs.emergence = probs.emergence %>% #pivot to long-form
  pivot_longer(cols = !R0, #keep R0_1 column
               names_to = "m", #add feature denoting "m" used row-wise
               values_to = "P.Emergence") #add feature denoting the pemergence calculated from m in m column and R0_1 in R0_1 column

## plot

# colors for plotting
#my.colors = c("black", "red", "forestgreen", "blue", "orange", "purple") #because pretty
# plot
ggplot() + 
  geom_line(data = probs.emergence, aes(x = R0, y = P.Emergence, group = m, col = m)) + #plot p emergence as function of R0_1 for each m
  theme_bw() + #white background with grey gridlines
  labs(x = "R0 of introduced pathogen", y = "Probability of Emergence", title = "probability of emergence with m phenotypes") + #axis and title labels
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(trans = "log10") + #log10 y scale
  scale_color_manual(values = my.colors) #because pretty



