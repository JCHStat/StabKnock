
####### functions for DS and MDS.
library(glmnet)
# library(DSfdr)

#' @param X Design matrix
#' @param y Response vector
#' @param num_split Repeated number of DS procedure for MDS
#' @param q FDR control level
#' @return A list containing the selected variables using DS or MDS.
####### auxiliary function
analys <- function(mm, ww, q){
  ### mm: mirror statistics
  ### ww: absolute value of mirror statistics
  ### q:  FDR control level
  cutoff_set <- max(ww)
  for(t in ww){
    ps <- length(mm[mm > t])
    ng <- length(na.omit(mm[mm < -t]))
    rto <- (ng + 1)/max(ps, 1)
    if(rto <= q){
      cutoff_set <- c(cutoff_set, t)
    }
  }
  cutoff <- min(cutoff_set)
  selected_index <- which(mm > cutoff)

  return(selected_index)
}                                                      

DS <- function(X, y, num_split, q){
  n <- dim(X)[1]; p <- dim(X)[2]
  inclusion_rate <- matrix(0, nrow = num_split, ncol = p)
  num_select <- rep(0, num_split)
  e=matrix(rep(0,num_split*p),nrow=num_split)
  index_e_na=c()
  for(iter in 1:num_split){
    M <- rep(0, p)
    cutoff <- Inf
    ### randomly split the data
    sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
    sample_index2 <- setdiff(c(1:n), sample_index1)
    
    ### get the penalty lambda for Lasso
    cvfit <- cv.glmnet(X[sample_index1, ], y[sample_index1], type.measure = "mse", nfolds = 10)
    lambda <- cvfit$lambda.min
    
    ### run Lasso on the first half of the data
    beta1 <- as.vector(glmnet(X[sample_index1, ], y[sample_index1], family = "gaussian", alpha = 1, lambda = lambda)$beta)
    nonzero_index <- which(beta1 != 0)
    if(length(nonzero_index)!=0){
      ### run OLS on the second half of the data, restricted on the selected features
      beta2 <- rep(0, p)
      beta2[nonzero_index] <- as.vector(lm(y[sample_index2] ~ X[sample_index2, nonzero_index] - 1)$coeff)
      
      ### calculate the mirror statistics
      M <- sign(beta1 * beta2) * (abs(beta1) + abs(beta2))
      # M <- abs(beta1 + beta2) - abs(beta1 - beta2)
      
      # selected_index = analys(M, abs(M), q)
      # DS_selected_index = selected_index
      mm=M
      ww=abs(M)
      cutoff_set <- max(ww)
      for(t in ww){
        ps <- length(mm[mm > t])
        ng <- length(na.omit(mm[mm < -t]))
        rto <- (ng + 1)/max(ps, 1)
        if(rto <= q){
          cutoff_set <- c(cutoff_set, t)
        }
      }
      cutoff <- min(cutoff_set)
      selected_index <- which(mm > cutoff)
      DS_selected_index = selected_index
      ### number of selected variables
      if(length(selected_index)!=0){
        num_select[iter] <- length(selected_index)
        inclusion_rate[iter, selected_index] <- 1/num_select[iter]
      }
    }else{
      DS_selected_index = NULL
    }
    ###compute e-value
    e[iter,]=p*(M>cutoff)/(1+sum(M<(-cutoff)))
    # e[iter,]=p*(M>cutoff)/(q*sum(M>=cutoff))
    if(sum(is.na(e[iter,]))>0){index_e_na=c(index_e_na,iter)}
  }
  ##由于Lasso可能出现初筛变量数大于n/2,导致部分e值为NA,下面先将NA的e值删掉
  if(length(index_e_na!=0)){e=e[-index_e_na,,drop=FALSE]}
  if(nrow(e)>0){
    e_value=colSums(e)/nrow(e)
  }else{
    e_value=rep(0,p)
  }
  ### multiple data-splitting (MDS) result
  inclusion_rate <- apply(inclusion_rate, 2, mean)
  
  ### rank the features by the empirical inclusion rate
  feature_rank <- order(inclusion_rate)
  feature_rank <- setdiff(feature_rank, which(inclusion_rate == 0))
  if(length(feature_rank)!=0){
    null_feature <- numeric()
    
    ### backtracking
    for(feature_index in 1:length(feature_rank)){
      if(sum(inclusion_rate[feature_rank[1:feature_index]]) > q){
        break
      }else{
        null_feature <- c(null_feature, feature_rank[feature_index])
      }
    }
    MDS_selected_index <- setdiff(feature_rank, null_feature)
  }else{
    MDS_selected_index = NULL
  }
  return(list(DS_feature = DS_selected_index, MDS_feature = MDS_selected_index, e_value=e_value))
}

