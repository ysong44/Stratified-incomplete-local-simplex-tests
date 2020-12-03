rm(list = ls())

nSim_gap = 400 #
nSim_rep = 1:4 #nSim_rep*nSim_gap: the total number of simulation repetitions
ds = c(2) #dimension of X
ns = c(1000)#sample size
kers = c(1, 2) # Simplex statistic kernel: ker=1 sign function; ker=2 identiy function
localKers = c(1) #Local kernels
b_ns = c(0.5)#bandwidth
regExponents = c(1)#parameter for polynomial regression function
reg_opt = 1 #parameter for generating samples from different regerssion functions
sigma_es = c(0.2)#noise level

#-------------------------
cnt = 1 
for(sigma_e in sigma_es){
 for(simplexKer in kers){
  for(localKerInd in localKers){
    for(regExponent in regExponents){
      for(d in ds){
        for(n in ns){
          for(b_n in b_ns){
            for(sim_ind in nSim_rep){
 	            filename = paste('inputConfig/inputConfig_', cnt ,'.RData', sep = '')
	            nSim = (nSim_gap*(nSim_rep[sim_ind]-1)+1) :(nSim_gap *nSim_rep[sim_ind])
	            r = d+2
	            p = ifelse(d==1, 10,ifelse(d == 2, 5, ifelse(d == 3, 3, 3)))# Discretization of U-process to p-dimensional U-statistic in each dimension
	            N =  min(10*p^d *n*b_n^(-d*r), choose(n,r)) #bernoulli sampling budget
	            JK_size = ifelse(n = 500, 20*n*b_n^(-d*r), 10*n*b_n^(-d*r)) #in DC step, how many combinations are chosen
	            NumCores = 32
              save('NumCores','nSim','n','d','r','p', 'b_n','N','JK_size','sigma_e','reg_opt','simplexKer', 'localKerInd','regExponent', file = filename)
              cnt = cnt + 1
            }
          }
        }  
      }
    }
  }
 }
}
