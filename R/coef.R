#' Extract coefficients from fitted model
#' 
#' Extract coefficients from fitted model for either a \code{'smle'} or \code{'selection'} object.
#' 
#' @param object Returned object from either the function \code{\link{SMLE}()} or \code{\link{smle_select}()}.
#' @param ... This argument is not used and listed for method consistency.
#' @return Fitted coefficients based on the screened or selected model specified in the object.  
#' 
#' @rdname coef
#' @method coef smle
#' @examples
#' 
#' set.seed(1)
#' Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
#' fit<-SMLE(Y=Data$Y, X=Data$X, k=15, family = "gaussian")
#' coef(fit)
#' fit_s<-smle_select(fit)
#' coef(fit_s)

#' @export
coef.smle<-function(object,...)
{
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  model <- object$X[,object$ID_retained]
  
  if(is.null(colnames(model))){ colnames(model) <- object$ID_retained }
  
  data = data.frame(Y = object$Y, X= model)
    
  fit<-glm(Y~.,data = data ,family = family)
  
  coef(fit)
  
  }
  
  
         
        
#' @rdname coef
#' @method coef selection
#' @export
coef.selection<-function(object,...)
{  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  model <- object$X[,object$ID_selected]
  
  if(is.null(colnames(model))){ colnames(model) <- object$ID_selected }
  
  data = data.frame(Y = object$Y, X= model)
  
  fit<-glm(Y~.,data = data ,family = family)
  
  coef(fit)
}
