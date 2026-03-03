matrix_model = function(births, deaths, init_pop, gens){#initialize leslie matrix with parameters from 1980
    L = diag(x = 1 - deaths[1, 1:100]) #firts year survivorship 1-death rate on diagonal
    L = rbind(births[1, 1:100], L) #add first year fecundities to top row
    L = cbind(L, c(births[1,101], rep(0,99), 1-deaths[1,101])) # add column ensures last row includes final death rate 
    tail(L) #looks good!
    
    # matrix to store population densities by age each year. rows will be age, columns will be year. 
    m = matrix(nrow = length(init_pop), ncol = gens+1) #44 columns for 44 years
    m[,1] = init_pop #initialize first column with populations structure in 1980
    
    #simulate population structure since 1980
    for(year in 1:gens){
      m[,year+1] = L %*% m[,year] #update population distribution
      
      #adjust leslie matrix with parameters from year
      L = diag(x = 1 - deaths[year+1, 1:100]) #firts year survivorship 1-death rate on diagonal
      L = rbind(births[year+1, 1:100], L) #add first year fecundities to top row
      L = cbind(L, c(births[year+1,101], rep(0,99), 1-deaths[year+1,101])) # add column ensures last row includes final survivorship
      
    }
    return(m)
    
}