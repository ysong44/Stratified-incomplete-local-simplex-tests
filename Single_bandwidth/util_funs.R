

# Check simplex qualification for (d+2) points in R^{d}
# Test code:
# d = 2
# V=matrix(rnorm((d+1)*(d+2)),d+1,d+2)
# V[d+1,] = V[d+1,] + colSums(V[1:d,])
# simplexStatValue(V,1)
#cnt = 1
#s = 0
#nsim = 10^5
#for(tt in 1:nsim){
#  d = 3
#  V=matrix(runif((d+1)*(d+2)),d+1,d+2)
#  V[d+1,] = 0.1*rnorm(d+2) + colSums(V[1:d,]^2)
#  tmp = simplexStatValue(V,1)
#  if(tmp != 0){
#    cnt = cnt + 1
#  }
#  s = s + tmp
#}
#c(cnt/nsim, s/nsim)
####
# To do: add a singularity check
#calculate the simplex statistics, without local kernel weighting
simplexStatValue = function(V, kerIndex){
  val = 0
  d_plus_one = nrow(V) - 1

  for(j in 1:(d_plus_one+1))
  {
    a = solve.default(V[1:d_plus_one,-j],V[1:d_plus_one,j])
    if(min(a)>0)
    {
      if(1==kerIndex)		# sign kernel
      {
        val = sign(sum(a*V[d_plus_one+1,-j])-V[d_plus_one+1,j]) 
      }
      if(2==kerIndex)		# identity kernel
      {
        val= (sum(a*V[d_plus_one+1,-j])-V[d_plus_one+1,j])
      }
      break
    }
  }
  
  return(val)
}

simplexStatValue_with_ind = function(V, kerIndex, localLogic){
  val = 0
  d_plus_one = length(localLogic) - 1
  
  for(j in 1:length(localLogic) )
  {
    if(!localLogic[j]) { next }
    a = solve.default(V[1:d_plus_one,-j],V[1:d_plus_one,j])
    if(min(a)>0)
    {
      if(1==kerIndex)		# sign kernel
      {
        val = sign(sum(a*V[d_plus_one+1,-j])-V[d_plus_one+1,j]) 
      }
      if(2==kerIndex)		# identity kernel
      {
        val= (sum(a*V[d_plus_one+1,-j])-V[d_plus_one+1,j])
      }
      break
    }
  }
  
  return(val)
}

####calculate the local kernel value for single observation at different discretization points
#res = localKerVal(matrix(runif(3*10^5),3,10^5), matrix(runif(3),3,1), 0.3, 1)
localKerVal = function(DiscrPts, dataX, b, localKerInd){#b is the bandwidth
  #first center and scale
  p = ncol(DiscrPts)
  d = nrow(DiscrPts)
  csVals = (DiscrPts - dataX %*% matrix(1,nrow =1 , ncol = p))/(b/2)
  #val = matrix(0, nrow = p, ncol = 1)  
  if(localKerInd == 1){ #uniform kernel
    val = (abs(csVals) <= 1)
  }
  if(localKerInd == 2){#Epanechnikov kernel
    val = ((abs(csVals) <= 1)* 3/4*(1-csVals^2))
  }
  if(localKerInd == 3){#Quartic  kernel
    val = ((abs(csVals) <= 1)* 15/16*((1-csVals^2)^2))
  } 
  if(localKerInd == 4){ #consine kernel
    val = ((abs(csVals) <= 1)* pi/4*cos(pi/2*csVals))
  }
  #  val = as.matrix(val)
  #resVal = matrix(1, nrow = p, ncol = 1)
  #for(id in 1:d){
  #  resVal = resVal *  val[id,]/b
  #}
  return(apply(val/b,2,prod))
}


# Regression function 1 (d=1): f(x) = x^{alpha}
# Concave function (alpha>=1): linear when alpha=1
reg_fun_1 = function(X, alpha)
{
  f = rowSums(X^(alpha))
  return(f)
}

# Regression function 2 (d=1): f(x) = 4*(x-0.5)^{3}
# Monotone function (contains both concave and convex parts): still quasi-concave
reg_fun_2 = function(X)
{
  f = 8*(X-0.5)^3
}

# Regression function 3 (d=1): GSV m_{4}
# Not quasi-concave: mostly concave, but there is a sharp dip at x=0.25
reg_fun_3 = function(X)
{
  X_l = X[X<0.5]
  X_u = X[X>=0.5]
  f_l = 10*(X_l-0.5)^3 - exp(-100*(X_l-0.25)^2)
  f_u = 0.1*(X_u-0.5) - exp(-100*(X_u-0.25)^2)
  f = c(f_l, f_u)
}

# Regression function 4 (d=1): 
# Quasi-concave but not concave: discontinuous
reg_fun_4 = function(X)
{
}


# Regression functions
# X: n*d design matrix (entries of X are bounded between 0 and 1)
reg_fun_meta = function(X, reg_opt, alpha)
{
  n = nrow(X)
  mf = rep(0,n)
  if(0==reg_opt) {
    mf = rep(0,n)		# zero function
  }
  if(1==reg_opt){
    mf = reg_fun_1(X,alpha) # Regression function
  }
  if(2==reg_opt){ mf = reg_fun_2(X)
  }
  if(3==reg_opt){ mf = reg_fun_3(X)
  }
  return(mf)
}


# Local simplex test statistic kernel for concavity
# Input parameters:
# V: (d+1)*(d+2), containing r=(d+2) vectors of v such that v = (x,y), a (d+1) vector
# ker=1 sign function; ker=2 identiy function
# localWeight[i, j]: L_(b_n)(query_i - DataPoint_j)
# Output: a flag indicating whether values are all zero
# (p*1) vector containing the values of the U-process evaluated at the p discretized points
####
compute_hvals = function(V, localW, ker_ind)
{
  res = list(FLAG = FALSE)
  if(any(abs(localW) > 10^(-16))){
    simpStatVal = simplexStatValue(V, ker_ind) #compute the simplex statistic
    if(abs(simpStatVal) > 10^(-16)){
      res$FLAG = TRUE
      res$h_val = localW*simpStatVal
    }
  }
  return(res) 
}	


compute_hvals_with_ind = function(V, localW, ker_ind, localLogic)
{
  res = list(FLAG = FALSE)
  if(any(abs(localW) > 10^(-16))){
    simpStatVal = simplexStatValue_with_ind(V, ker_ind, localLogic) #compute the simplex statistic
    if(abs(simpStatVal) > 10^(-16)){
      res$FLAG = TRUE
      res$h_val = localW*simpStatVal
    }
  }
  return(res) 
}	



#Generate query point
#equal space points in [l,u]^d, each dimension has p points
generate_Query = function(l, u, p, d){
  oneD = matrix(seq(l,u,length=p),nrow = 1)
  if(d == 1){
    queryPoints = oneD
  }else if(d == 2){
    tmpM = t(oneD) %*% rep(1,p)
    queryPoints =  rbind(matrix(tmpM, nrow = 1), matrix(t(tmpM), nrow = 1))
  }else if(d == 3){
    queryPoints = matrix(0, 3, p^3)
    queryPoints[1,] = rep(oneD,p^2)
    queryPoints[2,] = rep(matrix(matrix(1,p,1) %*% oneD, ncol = 1),p)
    queryPoints[3,] = matrix(matrix(1,p^2,1) %*% oneD, nrow = 1)
  }else if(d == 4){
    queryPoints = matrix(0, 4, p^4)
    queryPoints[1,] = rep(oneD,p^3)
    queryPoints[2,] = rep(matrix(matrix(1,p,1) %*% oneD, ncol = 1),p^2)
    queryPoints[3,] = rep(matrix(matrix(1,p^2,1) %*% oneD, nrow = 1),p)
    queryPoints[4,] = matrix(matrix(1,p^3,1) %*% oneD, nrow = 1)
  }else{
    print("d needs to be <= 4!")
  }
  return(queryPoints)
}

my_rbinom = function(ntrial, np){
  tmp_d = 3000
  tmp_q = ntrial %/% tmp_d
  tmp_r = ntrial %% tmp_d
  return(sum(rbinom(tmp_d, tmp_q, np)) + rbinom(1,tmp_r,np))
}
