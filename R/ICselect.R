# Multi-cores implements for features selection.

ebicc<-function(Y,X_s,family,tune,codingtype,keyset,
                
                k_min,k_max,n,pp,gamma,para,num_cores){
  
  fss<-function(v,family,codingtype,keyset){
    
    # Fit smle to retain v features of the model Y ~ X_s.
    # Return likelihood and degree of the freedom of the model.
    
    ff<- SMLE(Y=Y, X=X_s, k=v, family=family,codingtype = codingtype, group=TRUE,keyset=keyset)
    
    return(list(d_f= ff$df, likelihood = logLik(ff)))
    
  }
  
  
  if(para){
    
    select<-mclapply(k_min:k_max,fss,family=family,codingtype=codingtype,keyset=keyset,mc.cores = num_cores)
    
  }else {
    
    select<-lapply(k_min:k_max,fss,family=family,codingtype=codingtype,keyset=keyset)
    
  }
  
  
  select_crit<-switch(tune,
                      'ebic'=function(select,v){return(-2 *  select[[v]][[2]]  +   select[[v]][[1]]* log(n) + 2 *  select[[v]][[1]] * gamma * log(pp))},
                      'bic' =function(select,v){return(-2 *  select[[v]][[2]]  +   select[[v]][[1]]* log(n))},
                      'aic' =function(select,v){return(-2 *  select[[v]][[2]]  +   select[[v]][[1]]*2)}
  )
  
  
  return(unlist(lapply(1:(k_max-k_min+1),select_crit,select=select)))
   
}

