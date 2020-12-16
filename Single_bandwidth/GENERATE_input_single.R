source('util_funs.R') 

NumCores = 40 #specify the number of computing cores
b_n = 0.6 #the bandwidth parameter
simplexKer = 2 #2: use the sign of simplex stat; 1: use directly the simplex stat

#Specify the sample size and dimension of the covariates
n = 1000
d = 2
r = d+2

#Generate data. X: n*d, Y: vector of length d.  
X = matrix(runif(n*d,0,1), n, d)
y = reg_fun_meta(X, 1, 2) + 0.2 * rnorm(n)#Ex: Polynomial function, Gaussian noise

#Generate the query points: queryPoints a matrix of size d*(Num of Query Points)
#in the example, eachDim is the number of query points per dimension
eachDim = ifelse(d==1, 10,ifelse(d == 2, 5, ifelse(d == 3, 3, 3)))
queryPoints = generate_Query(0.3,0.7,eachDim,d)


#Computation parameters N and N_2 in the paper
#The computational time is mainly determined by JK_size
N =  10*ncol(queryPoints) *n*b_n^(-d*r) #N
JK_size =  ifelse(n <= 500, 20*n*b_n^(-d*r), 10*n*b_n^(-d*r)) #N_2


#Other parameters
localKerInd = 1 #kernel for localization; might just use default the uniform kernel

       
#Save the output to a file SILS_input.RData
#Then run Run_SILS_test.R to obtain a p value
save('NumCores','b_n','n','d','r', 'X','y', 'queryPoints', 'N','JK_size', 
     'simplexKer', 'localKerInd', file = 'INPUT_TO_SILS_Single.RData')
