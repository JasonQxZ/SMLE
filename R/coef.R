#' Extract coefficients from fitted model
#' 
#' Extract coefficients from fitted model for either a \code{'smle'} or \code{'selection'} object.
#' 
#' @param object Returned object from either the function \code{\link{SMLE}} or \code{\link{smle_select}}.
#' @param ... This argument is not used and listed for method consistency.
#' @return Fitted coefficients based on the screened or selected model specified in the object.  
#' 
#' @rdname coef
#' @method coef smle
#' @examples
#' 
#' set.seed(1)
#' Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
#' fit<-SMLE(Y=Data$Y, X=Data$X, k=9, family = "gaussian")
#' coef(fit)
#' fits<-smle_select(fit)
#' coef(fits)

#' @export
coef.smle<-function(object,...)
{
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_retained])
    
  fit<-glm(Y~.,data = data ,family = family)
  
  coef(fit)
  
  }
  
  
         
        
#' @rdname coef
#' @method coef selection
#' @export
coef.selection<-function(object,...)
{  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_selected])
  
  fit<-glm(Y~.,data = data ,family = family)
  
  coef(fit)
}
