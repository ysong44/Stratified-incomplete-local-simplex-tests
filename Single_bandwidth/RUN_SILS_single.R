# Simplex test for concavity of regression
rm(list=ls())
library(doParallel)
library(MASS)
library(mvtnorm)
library(EL)
library(Matrix)
library(Rfast)

source('bernoulliSampling.R') 
source('util_funs.R') 

load('INPUT_TO_SILS_Single.RData')

cat('NumCores: ', NumCores,' d:',d,' n:', n, ' b_n:', b_n,'\n')

stopImplicitCluster()
registerDoParallel(NumCores)

Verbose = TRUE

# Simulation setups
n_S1 = n	# Size of S_1
B = 1500			# Number of bootstraps
p = ncol(queryPoints)

# Bernoulli sampling parameters
p_n = min(N /choose(n,r),1)			# Bernoulli sampling parameter
q_n = min(JK_size/choose(n-1,r-1), 1)
alpha_n = n / N			
V = t(cbind(X, rep(1,n), y))		# V: (d+2)*n, the all one column is for simplex stat


#Algorithm starts
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
U_gamma_N = foreach(id_q = 1:p, .combine=cbind, .packages=c('Matrix','Rfast')) %dopar% {
	  
	  tmpUGN = matrix(0,nrow = 3, ncol = 1)
	  localDP = B_q_d[[id_q]]
	  tmp_len = length(localDP)
	  if(tmp_len < r-1){ next }
	  subSamples = bernoulliSampling(tmp_len, r, p_n)
	  T_1 = ifelse(is.null(subSamples), 0, nrow(subSamples))
	  T_2 = my_rbinom(choose(n,r) - choose(tmp_len,r), p_n)
	  #T_2 = rbinom(1, choose(n,r) - choose(tmp_len,r), p_n)
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
	  } 
	  
	  for(ell in 1:T_1){
	     if(tmpW[ell] == 0 || !any(simplexLog[ell, ])) { next }
	     temp_res = compute_hvals_with_ind(V[, conv_data_index[ell, ]], tmpW[ell], simplexKer, simplexLog[ell,])
	    if(temp_res$FLAG){
	      tmpUGN[1] = tmpUGN[1] + temp_res$h_val
	      tmpUGN[2] = tmpUGN[2] + temp_res$h_val^2
	    }
	  }

	  if(Verbose){
	    print(cat(id_q, T_1))
	  }
	  tmpUGN
}

U = U_gamma_N[1,]/U_gamma_N[3,]
gamma_hat_B =  U_gamma_N[2,]/U_gamma_N[3,] - U^2
###############################################################################
	
# Estimation of the Hajek projection
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
    JK_T_2 = my_rbinom(choose(n-1,r-1) - choose(nn,r-1), q_n)
    #JK_T_2 = rbinom(1, choose(n-1,r-1) - choose(nn,r-1), q_n)
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
    cat(idx,' ',length(query_index),'\n')
  }
  oneColRes
}
g_breve = rowMeans(g_hat)
gamma_hat_A = rowMeans(g_hat^2) - g_breve^2

	
### Bootstrap
U_max_without_normalization = sqrt(n)*max(U)
cnt_gt_times = 0
	for(b in 1:B){
		U_B = sqrt(gamma_hat_B) * rnorm(p)
		
		# Bootstrap \gamma_A
		xi_prime = rnorm(n_S1)
		U_A =  g_hat %*% matrix(xi_prime, ncol = 1) -sum(xi_prime)*g_breve
		U_A = as.matrix(U_A / sqrt(n_S1))
		
		# Combining U_A and U_B to get the bootstrapped U-statistic
		Uxi = r* U_A + sqrt(alpha_n) * U_B
		if(U_max_without_normalization <= max(r* U_A)){
		  cnt_gt_times = cnt_gt_times + 1
		}
	}
end_time = proc.time()
dur  = end_time - start_time
###########################################










#####OUTPUT the p-value
pval = cnt_gt_times/B
cat('p-value is', pval)



stopImplicitCluster()
