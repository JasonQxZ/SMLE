#' Extract log-likelihood
#' 
#' This is a method function written to extract log-liklihood from \code{'smle'} and \code{'selection'} objects. 
#' It refits the model by \code{\link[stats]{glm}} based on the response and the selected features after screening(selection), 
#' and returns an object of \code{\link[stats]{logLik}} from the generic.
#' @import stats
#' @param object Object of class \code{'smle'} or \code{'sdata'}. 
#' @param ... Forwarded arguments.
#' @return Returns an object of class \code{'logLik'}. This is a number with at least one attribute,
#'  \code{"df"} (degrees of freedom), giving the number of (estimated) parameters in the model. For more details, see the generic \code{'logLik'} in \pkg{stats}.
#' @rdname logLik
#' @method logLik smle
#' @examples
#' set.seed(1)
#' Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
#' fit<-SMLE(Y=Data$Y, X=Data$X, k=9, family = "gaussian")
#' logLik(fit)
#' 
#' @export 
logLik.smle<-function(object,...){
  codingtype = object$codingtype
  Y<-object$Y
  X_s <- object$X
  X_v <- X_s[,object$ID_retained]
  Ci <- sapply(X_v,is.factor)
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
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
  X_v <- X_s[,object$ID_selected]
  Ci <- sapply(X_v,is.factor)
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
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