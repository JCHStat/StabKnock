

### main function
Knock_SPD <- function(X, y, L = 100, alpha = 1, q = 0.1, cv=1, lambda=NULL){
      n = dim(X)[1]
      p = dim(X)[2]
      knock <- create.fixed(X, method= "equi")
      X <- knock$X
      Xk <- knock$Xk
      X_aug <- cbind(X,Xk)
      S_hat <- rep(0,2*p)
      if(cv==1){
            cv.fit = cv.glmnet(x = X_aug,y = y,family = "gaussian", intercept=F,alpha = 1,nfolds = 10, standardize=F)
            lambda = alpha * cv.fit$lambda.min
      }
      if(cv==0 & is.null(lambda)) lambda = alpha * sqrt(log(2*p)/n)
      if(cv==0 & is.null(lambda)==F) lambda=lambda
      for(i in 1:L){
        index <- sample((1:n), floor(n/2),replace = F)

          my.fit = glmnet(x = X_aug[index,],y = y[index],family = "gaussian", intercept=F,alpha = 1,lambda = lambda, standardize=F)
          b.hat <- as.vector(coef(my.fit)[-1])
          S_index1 <- which(b.hat!=0)
          
          my.fit = glmnet(x = X_aug[-index,],y = y[-index],family = "gaussian", intercept=F,alpha = 1,lambda = lambda, standardize=F)
          b.hat <- as.vector(coef(my.fit)[-1])
          S_index2 <- which(b.hat!=0)
        
        
        #S_index <- union(S_index1, S_index2)
        S_index <- intersect(S_index1, S_index2)
        
        S = rep(0,2*p)
        S[S_index]=1
        
        S_hat = S_hat + S
      }
      S_hat <- S_hat/L
    
      #compute Knockoff statistics
      W <- S_hat[1:p]-S_hat[-(1:p)]
      
      # find the knockoff threshold T
      t = sort(c(0, abs(W)))
      ratio = c(rep(0, p))
      
      for (j in 1:p) {
            ratio[j] = (sum(W <= -t[j]))/max(1, sum(W >= t[j]))
      }
      id = which(ratio <= q)[1]
      if(length(id) == 0){
            T = Inf
      } else {
            T = t[id]
      }
      
      # find the knockoff plus threshold T plus
      t_plus = sort(c(0, abs(W)))
      ratio_plus = c(rep(0, p))
      
      for (j in 1:p) {
            ratio_plus[j] = (1+sum(W <= -t_plus[j]))/max(1, sum(W >= t_plus[j]))
      }
      id_plus = which(ratio_plus <= q)[1]
      if(length(id_plus) == 0){
            T_plus = Inf
      } else {
            T_plus = t_plus[id_plus]
      }
      
      # set of discovered variables
      S <- which(W >= T)
      S_plus <- which(W >= T_plus)
      
      return(list(S=S,S_plus=S_plus,S_hat=S_hat,lambda=lambda, W=W, T=T))
}




### anciliary functions
fdp_power <- function(selected_index, signal_index){
  num_selected <- length(selected_index)
  tp <- length(intersect(selected_index, signal_index))
  fp <- num_selected - tp
  fdp <- fp / max(num_selected, 1)
  power <- tp / length(signal_index)
  return(list(fdp = fdp, power = power))
}


AIC=function(y,y_hat,n,d){
  RSS=sum((y-y_hat)^2)
  sigma2_hat=var(y-y_hat)
  AIC=(RSS+2*d*sigma2_hat)/(n*sigma2_hat)
  return(AIC)
}



