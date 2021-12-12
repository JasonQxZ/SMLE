#' Extract log-likelihood
#' 
#' This is a method written to extract the log-likelihood from \code{'smle'} and \code{'selection'} objects. 
#' It refits the model by \code{\link[stats]{glm}()} based on the response and the selected features after screening (selection), 
#' and returns an object of \code{'logLik'} from the generic.
#' @import stats
#' @param object An object of class \code{'smle'} or \code{'sdata'}. 
#' @param ... Forwarded arguments.
#' @return Returns an object of class \code{'logLik'}. This is a number with at least one attribute,
#'  \code{"df"} (degrees of freedom), giving the number of (estimated) parameters in the model. For more details, see the generic \code{\link[stats]{logLik}()} in \pkg{stats}.
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
  
  if(is.null(object$data)){
    Y<-object$Y
    X_s <- object$X
    X_v <- X_s[,object$ID_retained]
    
  }else{
    Y<-object$Y
    X_v <- object$data[,object$ID_retained]
  }
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())

  data = data.frame(Y = object$Y, X_v)
  if( !is.null(object$iteration_data$feature_name)){  names(data) <- c("Y",object$iteration_data$feature_name)}
 
  fit<-glm(Y~.,data = data ,family = family)
  

  ll <- stats::logLik(fit,...)
    return(ll)
    
  }


#' @import stats
#' @rdname logLik
#' @method logLik selection
#' @export 
logLik.selection<-function(object,...){
  if(is.null(object$data)){
    Y<-object$Y
    X_s <- object$X
    X_v <- X_s[,object$ID_selected]
    feature_name <- colnames(X_s)[object$ID_selected]
  }else{
    Y<-object$Y
    X_v <- object$data[,object$ID_selected]
    feature_name <- names(data)[object$ID_selected]
  }
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X_v)
  
  if( !is.null(feature_name)){  names(data) <- c("Y",feature_name)}
  
  fit<-glm(Y~.,data = data ,family = family)
    
  ll <- stats::logLik(fit,...)
    return(ll)

}