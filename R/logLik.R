#' Extract log-likelihood
#' 
#' @import stats
#' @param object Object of class \code{'smle'} or \code{'sdata'}. 
#' @param ... Forwarded arguments.
#' @return Returns an object of class logLik. This is a number with at least one attribute, "df" (degrees of freedom), giving the number of (estimated) parameters in the model.
#' @rdname logLik
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
    
    refit <-glm(Y ~ X_dummy, family=family)
    ll <- stats::logLik(refit,...)
    return(ll)
    
  }else{
    #X_v is a matrix
    X_v =  as.matrix(X_v)
    refit <-glm(Y ~ X_v, family=family)
    ll <- stats::logLik(refit,...)
    return(ll)
    
  }
}

#' @import stats
#' @rdname logLik
#' @method logLik selection
#' @export 
logLik.selection<-function(object,...){
  codingtype = object$codingtype
  Y<-object$Y
  X_s <- object$X
  X_v <- X_s[,object$ID_Selected]
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
    
    refit <-glm(Y ~ X_dummy, family=family)
    ll <- stats::logLik(refit,...)
    return(ll)
    
  }else{
    #X_v is a matrix
    X_v =  as.matrix(X_v)
    refit <-glm(Y ~ X_v, family=family)
    ll <- stats::logLik(refit,...)
    return(ll)
    
  }
}