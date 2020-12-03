# Simplex test for concavity of regression

rm(list=ls())
library(doParallel)
library(MASS)
library(mvtnorm)
library(EL)
library(Matrix)
library(Rfast)

source('bernoulliSampling.R') #suppose we only consider Bernoulli sampling
source('util_funs.R') 

args = commandArgs(trailingOnly=TRUE)
cat(paste('inputConfig/inputConfig_',args[1],'.RData\n',sep=''))
load(paste('inputConfig/inputConfig_',args[1],'.RData',sep=''))

cat('d:',d,' n:', n, ' b_n:', b_n,' regExponent:', regExponent, ' sigma_e:', sigma_e,'\n')

stopImplicitCluster()
registerDoParallel(NumCores)

Verbose = TRUE

# Simulation setups
#nSim = 1:300		# Number of simulations, just for testing!!!!
#d = 3		# Number of covariates
#r = d+2			# Kernel order of simplex test statistic
#n = 1000			# Sample size
#b_n = 0.65 #n^(-1/(d+5))
#p = ifelse(d==1, 10,ifelse(d == 2, 5, ifelse(d == 3, 3, 3)))# Discretization of U-process to p-dimensional U-statistic in each dimension
#N =  min(10*p^d *n*b_n^(-d*r), choose(n,r)) #bernoulli sampling budget
#JK_size = p^d*n*b_n^(-d*r) #in DC step, how many combinations are chosen
n_S1 = n	# Size of S_1
B = 1500			# Number of bootstraps
#simplexKer = 1			# Simplex statistic kernel: simplexKer=1 sign function; simplexKer=2 identiy function
#localKerInd = 1#local kernel index: ker =1 uniform; ker =2 Epanechnikov; ker = 3 Quartic ; ker = 4 Cosine 
#reg_opt = 1		# Regression function options (see details below)
#regExponent = 1 #if reg_opt == 1
# Data generation parameters
#sigma_e = 0.1		# Standard deviation of the mean-zero error term


#Query points
queryPoints = generate_Query(0.3,0.7,p,d)
p = ncol(queryPoints)

# Bernoulli sampling parameters
p_n = min(N /choose(n,r),1)			# Bernoulli sampling parameter
q_n = min(JK_size/choose(n-1,r-1), 1)
alpha_n = n / N			


# Output parameters
alpha = seq(0.01,0.99,by=0.01)
L_alpha = length(alpha)
alpha_prob0 = matrix(0,length(nSim),L_alpha)
alpha_prob0_without_normalization_with_B = matrix(0,length(nSim),L_alpha)
alpha_prob0_without_normalization_without_B = matrix(0,length(nSim),L_alpha)

U_max = c()
U_max_without_normalization = c()

Uxi_max = matrix(0,length(nSim),B)
Uxi_max_without_normalization_with_B = matrix(0,length(nSim),B)
Uxi_max_without_normalization_without_B = matrix(0,length(nSim),B)

# Simulation loop
dur = c()
#SimRes = foreach(ii = 1:length(nSim), .combine=rbind, .packages=c('Matrix','Rfast')) %dopar% {
for(ii in  1:length(nSim)){
	cat('Simulation', ii, '\n')
	# Make data
	X = matrix(runif(n*d,0,1), n, d)

	bin_sel = rbinom(n,1,0.5)
	epsilon = bin_sel*rnorm(n, -0.1, 0.06) + (1-bin_sel)* rnorm(n, 0.1, 0.24)

	y = reg_fun_meta(X, reg_opt, regExponent) + epsilon
	V = t(cbind(X, rep(1,n), y))		# V: (d+2)*n, the all one column is for simplex stat

	start_time = proc.time()
	
	#pre-processing
	localWeightMat = Matrix(0, n, p, sparse = TRUE) #data by query points
	B_q_d = vector("list", length(p)) #Find for each query neighbor data points
	for(id_q in 1:p){
	  res = localKerVal(V[1:d, ,drop=FALSE], queryPoints[, id_q, drop=FALSE], b_n, localKerInd) 
	  res_nzero = which(abs(res) > 10^(-20))
	  localWeightMat[res_nzero, id_q] = res[res_nzero]
	  B_q_d[[id_q]] = res_nzero
	}

	B_data_query = vector("list", length(n)) #Find for each query neighbor data points
	for(idx in 1:n){
	  B_data_query[[idx]] = which( abs(localWeightMat[idx, ]) > 10^(-20) )
	}
	
	# Bernoulli sampling to compute the incomplete U-statistic
	#U = rep(0,p)
	#gamma_hat_B = rep(0,p)
	#N_hat =  rep(0,p)
	
#	for(id_q in 1:p){ #partitioned sampling
	U_gamma_N = foreach(id_q = 1:p, .combine=cbind, .packages=c('Matrix','Rfast')) %dopar% {
	  
	  tmpUGN = matrix(0,nrow = 3, ncol = 1)
	  localDP = B_q_d[[id_q]]
	  tmp_len = length(localDP)
	  if(tmp_len < r-1){ next }
	  subSamples = bernoulliSampling(tmp_len, r, p_n)
	  T_1 = ifelse(is.null(subSamples), 0, nrow(subSamples))
	  T_2 = rbinom(1, choose(n,r) - choose(tmp_len,r), p_n)
	  #N_hat[id_q] = T_1 + T_2
	  tmpUGN[3] = T_1 + T_2
	  if(T_1 == 0){ next }
	  
	  tmplocalWeightMat = as.matrix(localWeightMat[, id_q])
	  tmpW = matrix(1, T_1, 1)
	  conv_data_index = matrix(1, T_1, r)
	  tmp_max = matrix(0, d, T_1)
	  tmp_min = matrix(1, d, T_1)
	  simplexLog = matrix(FALSE,T_1,r)
	  for(tmpi in 1:r){
	    tmpInd = localDP[subSamples[,tmpi]]
	    conv_data_index[,tmpi] = tmpInd #Find correct index, and weight
	    tmpW = tmpW * tmplocalWeightMat[tmpInd]
	    tmp_V = V[1:d, tmpInd]
	    tmp_max = pmax(tmp_max, tmp_V)
	    tmp_min = pmin(tmp_min, tmp_V)
	  }
	  for(tmpi in 1:r){
	    tmp_V = V[1:d, conv_data_index[,tmpi] ]
	    simplexLog[ , tmpi] = colAll(tmp_V < tmp_max) & colAll(tmp_V > tmp_min)
	    #(apply(tmp_V < tmp_max, 2, all)) & (apply(tmp_V > tmp_min,2,all))
	  } 
	  
	  for(ell in 1:T_1){
	     if(tmpW[ell] == 0 || !any(simplexLog[ell, ])) { next }
	     temp_res = compute_hvals_with_ind(V[, conv_data_index[ell, ]], tmpW[ell], simplexKer, simplexLog[ell,])
	    if(temp_res$FLAG){
	      tmpUGN[1] = tmpUGN[1] + temp_res$h_val
	      tmpUGN[2] = tmpUGN[2] + temp_res$h_val^2
	      #U[id_q] = U[id_q] + temp_res$h_val
	      #gamma_hat_B[id_q] = gamma_hat_B[id_q] + temp_res$h_val^2
	    }
	  }

	  if(Verbose){
	    print(cat(id_q, T_1))
	  }
	  tmpUGN
	}
#	U = U/N_hat
#	gamma_hat_B = gamma_hat_B/N_hat - U^2
	U = U_gamma_N[1,]/U_gamma_N[3,]
	gamma_hat_B =  U_gamma_N[2,]/U_gamma_N[3,] - U^2
	###############################################################################
	
	# Estimation of the Hajek projection

	#g_hat = matrix(0,p,n_S1)
	#g_breve = rep(0,p)
	#gamma_hat_A = rep(0,p)
  #for(idx in 1:n_S1){#algorithm_two modified
	g_hat = foreach(idx = 1:n_S1, .combine=cbind, .packages=c('Matrix','Rfast')) %dopar% {
  	  oneColRes = matrix(0,nrow = p, ncol = 1)
      query_index = B_data_query[[idx]]
      for(qi in query_index){
        nn_dp = setdiff(B_q_d[[qi]], idx) #index for nearby query points, excluding idx
        nn = length(nn_dp)
        if(nn < r-1){ next }
        JK_incomp_ustat_idx = bernoulliSampling(nn, r-1, q_n)
        if(is.null(JK_incomp_ustat_idx)){ next }
        JK_T_1 = nrow(JK_incomp_ustat_idx)
        JK_T_2 = rbinom(1, choose(n-1,r-1) - choose(nn,r-1), q_n)
        tmplocalWeightMat =  as.matrix(localWeightMat[, qi])

        ##########Do some preparation###########
        conv_data_index = matrix(idx,JK_T_1,r)
        vector_idx = V[1:d, idx] %*% matrix(1,1,JK_T_1)
        tmp_max = vector_idx
        tmp_min = vector_idx
        tmp_weight = rep(tmplocalWeightMat[idx], JK_T_1)

        for(tmpi in 2:r){
          tmp_ind = nn_dp[JK_incomp_ustat_idx[,tmpi-1]] #converted back to the correct indexing
          conv_data_index[,tmpi] = tmp_ind
          tmp_weight = tmp_weight * tmplocalWeightMat[tmp_ind]
          tmp_V = V[1:d, tmp_ind]
          tmp_max = pmax(tmp_max, tmp_V)
          tmp_min = pmin(tmp_min, tmp_V)
        }

        simplexLog =  matrix(FALSE, JK_T_1, r) #some points are clearly not an inner point
        simplexLog[ , 1] = colAll(vector_idx < tmp_max) & colAll(vector_idx > tmp_min)
        #(apply(vector_idx < tmp_max, 2, all)) & (apply(vector_idx > tmp_min,2,all))

        for(tmpi in 2:r){
          tmp_V = V[1:d, conv_data_index[,tmpi] ]
          simplexLog[ , tmpi] = colAll(tmp_V < tmp_max) & colAll(tmp_V > tmp_min) 
        #(apply(tmp_V < tmp_max, 2, all)) & (apply(tmp_V > tmp_min,2,all))
        }
      
        #####For each selected tuple
        for(ell in 1:JK_T_1){
          if(tmp_weight[ell] == 0 || !any(simplexLog[ell, ]) ){ next } 
          tmp_res = compute_hvals_with_ind(V[, conv_data_index[ell, ]], tmp_weight[ell], simplexKer, simplexLog[ell, ])
          
          if(tmp_res$FLAG){
            oneColRes[qi] = oneColRes[qi] + tmp_res$h_val
          }
        }
       oneColRes[qi] = oneColRes[qi]/(JK_T_1+JK_T_2)
    }
    if(Verbose){
     cat('Sim_', nSim[ii] ,'_U_A', idx,' ',length(query_index),'\n')
    }
    oneColRes
	}
	g_breve = rowMeans(g_hat)
	gamma_hat_A = rowMeans(g_hat^2) - g_breve^2
  #g_breve = g_breve/n_S1
  #gamma_hat_A = gamma_hat_A/n_S1 - g_breve^2

	
  ### compute test statistics
	const_hat = sqrt(r^2*gamma_hat_A + alpha_n*gamma_hat_B)
	#const_hat = r*sqrt(gamma_hat_A)
  U_max[ii] = sqrt(n)*max(U/const_hat)
	U_max_without_normalization[ii] = sqrt(n)*max(U)
	# Bootstrap loop
	for(b in 1:B){
	  #bootstrapRes = foreach(i = 1:B, .combine = cbind) %dopar% {

		# Bootstrap \gamma_B
	  #xi = rnorm(M_B_len)
		#U_B = #M_B %*% matrix(xi, ncol = 1) - (sum(xi) + rnorm(1,0,sqrt(N_hat - M_B_len))) * U
		#U_B = as.matrix(U_B / sqrt(N_hat)) 
		U_B = sqrt(gamma_hat_B) * rnorm(p)
		
		# Bootstrap \gamma_A
		xi_prime = rnorm(n_S1)
		U_A =  g_hat %*% matrix(xi_prime, ncol = 1) -sum(xi_prime)*g_breve
		U_A = as.matrix(U_A / sqrt(n_S1))
		
		# Combining U_A and U_B to get the bootstrapped U-statistic
		Uxi = r* U_A + sqrt(alpha_n) * U_B
		#Uxi = r* U_A
		Uxi_max[ii,b] = max(Uxi/const_hat)
		Uxi_max_without_normalization_with_B[ii,b] = max(Uxi)
  	Uxi_max_without_normalization_without_B[ii,b] = max(r* U_A)
		
		#max(Uxi/const_hat)
	
	}
	#Uxi_max[ii, ] = bootstrapRes

	end_time = proc.time()
	dur[ii] = end_time[1] - start_time[1]

	# Compare the quantiles
	 quant0 = quantile(Uxi_max[ii,], alpha)
	 quant0_without_normalization_with_B = quantile(Uxi_max_without_normalization_with_B[ii,], alpha)
	 quant0_without_normalization_without_B = quantile(Uxi_max_without_normalization_without_B[ii,], alpha)
	 
	 
	 alpha_prob0[ii,] = (U_max[ii] <= quant0)
	 alpha_prob0_without_normalization_with_B[ii,] = (U_max_without_normalization[ii] <= quant0_without_normalization_with_B)
	 alpha_prob0_without_normalization_without_B[ii,] = (U_max_without_normalization[ii] <= quant0_without_normalization_without_B)
	 
	 data2Save = list(U_max = U_max[ii], Uxi_max = Uxi_max[ii,], alpha_prob0 = alpha_prob0[ii,], dur = end_time - start_time,
	                  alpha_prob0_without_normalization_with_B =  alpha_prob0_without_normalization_with_B[ii,],
	                  alpha_prob0_without_normalization_without_B = alpha_prob0_without_normalization_without_B[ii,],
	                  r_gamma_hat_A =  r^2* gamma_hat_A, alpha_n_gamma_hat_A = alpha_n * gamma_hat_B, bn = b_n,
	                  ind_max = which.max(U/const_hat), ind_max_without_norm = which.max(U))
	 #save('data2Save', file = paste('tmpFiles/n_', n, '_d_',d, '_nsim_', nSim[ii], '.RData', sep = ''))
	 save('data2Save', file = paste('tmpFiles/conf_', args[1], '_n_', n, '_d_',d, '_nsim_', nSim[ii], '.RData', sep = ''))
	 
	 alpha_prob0[ii,]
}

alpha_prob0 = SimRes
alpha_hat0 = colMeans(alpha_prob0)

i_alpha = seq(1, L_alpha, by=1)
plot(alpha[i_alpha], alpha_hat0[i_alpha],xlab='alpha',ylab='Bootstrap approximation', main='Gaussian multiplier bootstrap', pch='+', type='o')
abline(0,1,lty=2)

#save.image('concavity_test.RData')
save('alpha', 'alpha_hat0','dur', 'nSim', 'B','d', 'n','p','simplexKer', 'b_n', file = 'Desktop_concavity_test.RData')
stopImplicitCluster()
