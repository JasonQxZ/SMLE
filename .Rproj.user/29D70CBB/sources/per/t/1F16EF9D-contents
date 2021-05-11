#tools function for selection.
#Multicores implenments for features selection.
ebicc<-function(Y,X_s,family,tune,k_min,k_max,

                n,pp,gamma,para,num_cores){
  #fast_smle_by_sparisty function
  #return retained feature ids for given sparisity
  fss<-function(v){

    ff<- SMLE(Y=Y, X=X_s, k=v, family=family,fast=TRUE)


    return(lh(Y=Y, X=as.matrix(X_s[,ff$ID_Retained]),
              beta=as.vector(glm(Y ~ as.matrix(X_s[,ff$ID_Retained])-1, family=family)$coefficients),
              family=family))
  }

  ll<-unlist(lapply(k_min:k_max,fss))

  select_crit<-switch(tune,
                      'ebic'=function(ll,v){return(-2 * ll  + v* log(n) + 2 * v * gamma * log(pp))},
                      'bic' =function(ll,v){return(-2 * ll +  v* log(n))},
                      'aic' =function(ll,v){return(-2 * ll +  v*2)}
                      )
  return(mapply(select_crit,ll,k_min:k_max))
}


ctg_ebicc<-function(Y,X_s,family,tune,codingtype,

                   k_min,k_max,n,pp,gamma,para,num_cores){

  #fast_smle_by_sparisty function

  #X_s contains categorical features.

  fss<-function(v,codingtype){

    ff<- SMLE(Y=Y, X=X_s, k=v, family=family,codingtype = codingtype,categorical = T,group=T)

    X_v <- X_s[,ff$ID_Retained]

    Ci <- sapply(X_v,is.factor)

    if( any(sapply(X_v,is.factor)) ){
      # if X_v contains categorical  features

      X_dummy <- as.matrix(suppressWarnings(dummy.data.frame(X_v ,sep="_",codingtype = codingtype)))

      if(codingtype =="all"){

        dummy_f <- sum(sapply(list(X_v[,Ci]),nlevels)-1)

      }else{

        dummy_f <- sum(sapply(list(X_v[,Ci]),nlevels)-2)

      }

      ll <- lh(Y=Y, X=X_dummy,   beta=as.vector(glm(Y ~ X_dummy-1, family=family)$coefficients), family=family)
      return(list(d_f= dummy_f + v, likelihood = ll))

    }else{
      #X_v is a matrix
      X_v =  as.matrix(X_v)

      ll <- lh(Y=Y, X=X_v,
                beta=as.vector(glm(Y ~ X_v-1, family=family)$coefficients),
                family=family)

      return(list(d_f= v, likelihood = ll))

    }
  }

  select<-lapply(k_min:k_max,fss,codingtype=codingtype)
  select_crit<-switch(tune,
                      'ebic'=function(select,v){return(-2 *  select[[v]][[2]]  +   select[[v]][[1]]* log(n) + 2 *  select[[v]][[1]] * gamma * log(pp))},
                      'bic' =function(select,v){return(-2 *  select[[v]][[2]]  +   select[[v]][[1]]* log(n))},
                      'aic' =function(select,v){return(-2 *  select[[v]][[2]]  +   select[[v]][[1]]*2)}
  )


  return(unlist(lapply(k_min:k_max,select_crit,select=select)))
}
