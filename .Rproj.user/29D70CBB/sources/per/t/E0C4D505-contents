#' Extract coefficient from fitted model
#' 
#' @param object Fitted model
#' @param ... Forwarded arguments
#' @return Coefficients after screening or selection.
#' 
#' @rdname coef
#' @method coef smle
#' @export
coef.smle<-function(object,...)
{return(object$Coef_Retained)}
#' @rdname coef
#' @method coef selection
#' @export
coef.selection<-function(object,...)
{return(object$Coef_Selected)}
