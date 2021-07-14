#' Extract coefficients from fitted model
#' 
#' Extract coefficients from fitted model for either a \code{'smle'} or \code{'selection'} object.
#' 
#' @param object Returned fitted model object from either the function \code{SMLE} or \code{smle_select}.
#' @param ... This argument is not used and listed for method consistency.
#' @return Coefficients extracted from the fitted model object of class \code{'smle'} or class \code{'selection'}. 
#' 
#' @rdname coef
#' @method coef smle
#' @examples
#' 
#' Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
#' fit<-SMLE(Y=Data$Y, X=Data$X, k=9, family = "gaussian")
#' coef(fit)

#' @export
coef.smle<-function(object,...)
{return(object$coef_retained)}
#' @rdname coef
#' @method coef selection
#' @seealso SMLE
#' @export
coef.selection<-function(object,...)
{return(object$coef_selected)}