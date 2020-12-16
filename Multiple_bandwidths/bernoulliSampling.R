#test 1: library(tictoc); tic(); res = bernoulliSampling(500,5,500000); toc()
#test 2: library(tictoc); tic(); res = bernoulliSampling(1000,4,333*1000); toc()
#I_{n,r} = {(i_1,..., i_r): 1 <= i_1 < i_2 < ... < i_r <= n} ordered in the obvious way; 
#each subsample is included with prob p_n = N/nchoosek(n,r)
bernoulliSampling = function(n,r,p_n){
  tt= choose(n,r)
  subSamples = selSubsets(n,r,matrix(sample(tt,rbinom(1, tt, p_n)), ncol = 1))
  return(subSamples)
}



#test 1:  res = selSubsets( 10, 5, matrix(c(2, 5, 10, 75, 100, 125,  150, 175, 200,  250,  252), ncol = 1) )
#test 2:  res = selSubsets( 500, 4, matrix(seq(1, by = 10^4, length.out = 250000), ncol = 1) )
#I_{n,r} = {(i_1,..., i_r): 1 <= i_1 < i_2 < ... < i_r <= n} ordered in the obvious way; 
#indexM is a collection of linear index (!!!!!a column vector!!!!)
#selection subsets from I_{n,r} indexed by indexM
selSubsets = function ( n, r, indexM )
{
  if(length(indexM) == 0){
    return(c())
  }
  subsamples = matrix(0, nrow(indexM),r)
  nM = matrix(n, nrow(indexM),1)
  for(id in 1:r) #recursively find the indexes
    {#tmpSum is the total number with the first item < subsamples(:,id)-1
      tmpRes = findFirstItem(nM, r + 1 - id, indexM)
      subsamples[,id] = tmpRes[,1]
      indexM = indexM - tmpRes[,2]
      nM = nM - tmpRes[,1]
  }
  subsamples = t(apply(subsamples,1, cumsum));
  return(subsamples)
}



#find the first item given linear index
#firstItem is the first item with indexM
#tmpNum is the total combination with first item < firstItem -1
findFirstItem = function(nM, selR, indexM)
{
  tmpNum  = matrix(0, nrow(indexM), 1);
  if(selR == 1){
    firstItem = indexM
  }
  else{
    firstItem = matrix(0, nrow(indexM), 1)
    cumNum = matrix(0, nrow(indexM), 1)
    activeInd = matrix(1:nrow(indexM), nrow(indexM), 1)  
    for(id in 1:max(nM)){
      #nchoosek(nM - id, selR-1)
      tmpCumNum = cumNum + choose(nM-id,selR-1); #the number of combinations with id being the first
    
      tmpInd = (cumNum < indexM[activeInd]) & (tmpCumNum >= indexM[activeInd])
      firstItem[activeInd[tmpInd]] = id
      tmpNum[activeInd[tmpInd]] = cumNum[tmpInd]
      
      if(all(tmpInd)){
        break
      }
      else{
        activeInd = activeInd[!tmpInd]
        cumNum = tmpCumNum[!tmpInd]
        nM = nM[!tmpInd]
      }    
    }
  }
  return(cbind(firstItem, tmpNum))
}


    

