#' Summarize SMLE-screening and selection
#'
#' @rdname summary
#' @description This function prints a summary of a \code{'smle'} (or a \code{'selection'}) object.
#' In particular, it shows the features retained after SMLE-screening (or selection) with the related convergence information.
#' @import stats
#' @param object A \code{'smle'} or \code{'selection'} object.
#'
#' @param ... This argument is not used and listed for method consistency.
#'
#' @return
#' No return value.
#' @export
#' @method summary smle
#' @examples
#' set.seed(1)
#' Data <- Gen_Data(correlation = "MA", family = "gaussian")
#' fit <- SMLE(Y = Data$Y, X = Data$X, k = 20, family = "gaussian")
#' summary(fit)
#' fit_s <- smle_select(fit)
#' summary(fit_s)
#'
summary.smle <- function(object, ...){
  
  lg<-logLik(object)
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial" = binomial(), "poisson"=  poisson())
  
  ans=list()
  ans$family = object$family
  ans$steps = object$steps
  ans$loglikelihood = lg
  ans$DimY= length(object$Y)
  ans$DimX= dim(object$X)
  ans$size = object$num_retained
  ans$ID_retained  =object$ID_retained
  ans$coef_retained = object$coef_retained
  ans$X = object$X
  
  
  if( !is.null(object$Intercept) ){
    ans$intercept= object$intercept
    
  }else{ans$intercept =NULL}
  
  
  if(object$ctg==TRUE){
    
    ans$ctg = TRUE
    ans$CI = object$I$CI
    ans$levels = nlevels(object$X[,object$I$CI])
  
  }
  
  ans$call <- object$call
  class(ans)= "summary.smle"
  ans
}

#' @import stats
#' @method summary selection
#' @export
#' @rdname summary
summary.selection <- function(object, ...){
  
  lg<-logLik(object)
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  ans=list()
  ans$family = object$family
  ans$steps = object$steps
  ans$loglikelihood = lg
  ans$DimY= length(object$Y)
  ans$DimX= dim(object$X)
  ans$size = object$num_selected
  ans$coef = object$coef_selected
  ans$vote = object$vote
  ans$ID_voted = object$ID_voted
  ans$criterion = object$criterion
  ans$ID_selected  =object$ID_selected
  ans$coef_selected = object$coef_selected
  ans$gamma_ebic = object$gamma_ebic 
  ans$X = object$X
  
  if( !is.null(object$Intercept) ){
    ans$intercept= object$intercept
    
  }else{ans$intercept =NULL}
  
  if(object$ctg==TRUE){
    ans$ctg = TRUE
    ans$CI = object$I$CI
    ans$levels = nlevels(object$X[,object$I$CI])
    
  }
  
  ans$call <- object$call
  class(ans)= "summary.selection"
  ans
}



