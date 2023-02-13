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
  
  ans=list( family = object$family, steps = object$steps , loglikelihood = lg,
            DimY= length(object$Y), DimX= dim(object$X) , size = object$num_retained,
            ID_retained  =object$ID_retained, coef_retained = object$coef_retained,
            X = object$X , data = object$data , iteration_data = object$iteration_data)
  
  if( !is.null(object$intercept) ){
     ans$intercept= object$intercept
    
  }else{ans$intercept =NULL}
  
  
  if(object$categorical==TRUE){
    
    ans$categorical = TRUE
    ans$CI = object$modified_data$CI
    ans$levels = nlevels(object$X[,object$modified_data$CI])
  
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
  
  ans=list( family = object$family, steps = object$steps , loglikelihood = lg,
            DimY= length(object$Y), DimX= dim(object$X) , size = object$num_selected,
            ID_selected  =object$ID_selected , coef_selected = object$coef_selected,
            vote = object$vote ,  ID_voted = object$ID_voted , criterion = object$criterion,
            criterion = object$criterion , gamma_ebic = object$gamma_ebic, 
            X = object$X , data = object$data , iteration_data = object$iteration_data)

  if( !is.null(object$intercept) ){
    ans$intercept= object$intercept
    
  }else{ans$intercept =NULL}
  
  if(object$categorical==TRUE){
    ans$categorical = TRUE
    ans$CI = object$modified_data$CI
    ans$levels = nlevels(object$X[,object$modified_data$CI])
    
  }
  
  ans$call <- object$call
  class(ans)= "summary.selection"
  ans
}



