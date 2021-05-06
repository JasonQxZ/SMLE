#' Extract log-likelihood
#' @export
#' @param object Object of class \code{'smle'} or \code{'sdata'}. 
#' @param ... Forwarded arguments.
#' @return Value of log-likelihood
#' 
#' @method logLik smle
#' @export 
logLik.smle<-function(object,...){
  codingtype = object$codingtype
  Y<-object$Y
  X_s <- object$X
  X_v <- X_s[,object$ID_Retained]
  Ci <- sapply(X_v,is.factor)
  family <- object$family
  if( any(sapply(X_v,is.factor)) ){
    # if X_v contains categorical  features
    
    X_dummy <- as.matrix(suppressWarnings(dummy.data.frame(X_v ,sep="_",codingtype =codingtype)))
    
    if(codingtype =="all"){
      
      dummy_f <- sum(sapply(list(X_v[,Ci]),nlevels)-1)
      
    }else{
      
      dummy_f <- sum(sapply(list(X_v[,Ci]),nlevels)-2)
      
    }
    
    ll <- lh(Y=Y, X=X_dummy,   beta=as.vector(glm(Y ~ X_dummy-1, family=family)$coefficients), family=family)
    return(ll)
    
  }else{
    #X_v is a matrix
    X_v =  as.matrix(X_v)
    ll <- lh(Y=Y, X=X_v,
             beta=as.vector(glm(Y ~ X_v-1, family=family)$coefficients),
             family=family)
    
    return(ll)
    
  }
}