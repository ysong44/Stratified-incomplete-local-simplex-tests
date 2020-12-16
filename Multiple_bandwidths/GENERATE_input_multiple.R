source('util_funs.R') 

NumCores = 40 #specify the number of computing cores
MultiBns = c(0.6, 0.8, 1) #specify multiple bandwidths
simplexKer = 2 #2: use the sign of simplex stat; 1: use directly the simplex stat

#Specify the sample size and dimension of the covariates
n = 1000
d = 2
r = d+2

#Generate data. X: n*d, Y: vector of length d.  
X = matrix(runif(n*d,0,1), n, d)
y = reg_fun_meta(X, 1, 1.5) + 0.2 * rnorm(n)#Ex: Polynomial function, Gaussian noise

#Generate the query points: queryPoints a matrix of size d*(Num of Query Points)
#in the example, eachDim is the number of query points per dimension
eachDim = ifelse(d==1, 10,ifelse(d == 2, 5, ifelse(d == 3, 3, 3)))
queryPoints = generate_Query(0.3,0.7,eachDim,d)


#Computation parameters N and N_2 for each bandwidth
#The computational time is mainly determined by MultiJKs
MultiNs =  10*ncol(queryPoints) *n*MultiBns^(-d*r) #multiple Ns
MultiJKs = ifelse(n <= 500, 20*n*MultiBns^(-d*r), 10*n*MultiBns^(-d*r)) #Multiple N_2s



#Other parameters
localKerInd = 1 #kernel for localization; might just use default the uniform kernel


#Save the output to a file SILS_input.RData
#Then run Run_SILS_test.R to obtain a p value
save('NumCores','MultiBns','n','d','r', 'X','y', 'queryPoints', 'MultiNs','MultiJKs', 
     'simplexKer', 'localKerInd', file = 'INPUT_TO_SILS_Multiple.RData')

