#' summary a SMLE object from SMLE
#'
#' @rdname summary
#' @description This functions prints a summary of a SMLE object.
#' In particular, it shows the features retained after SMLE-screening
#' and the related convergence information.
#' @import stats
#' @param object Fitted '\code{smle}' object.
#'
#' @param ... This argument is not used and listed for method consistency.
#'
#' @return
#' No return value.
#' @export
#' @method summary smle
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' summary(fit)
#'
summary.smle <- function(object, ...){
  
  lg<-logLik(object)
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  ans=list()

  
  if(object$ctg==TRUE){
    
    ans$ctg = TRUE
    ans$CI = object$I$CI
    ans$levels = nlevels(object$X[,object$I$CI])
    
    
    
  }else{
    

    
    
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
  
  
  ans=list()
  
  
  if(object$ctg==TRUE){
    
    data = data.frame(Y = object$Y, X= object$X[,object$ID_Selected])
    
    data<- suppressWarnings(dummy.data.frame(data ,sep="."))
    

    ans$ctg = TRUE
    ans$CI = object$I$CI
    ans$levels = nlevels(object$X[,object$I$CI])
    
    
    
  }else{
    
    if(is.null(colnames(object$X))){
      
      X = object$X
      
      colnames(X) <- paste0("X.",seq(1,ncol(X)))
    }
    

    
    
  }
  
  ans$call <- object$call
  class(ans)= "summary.smle"
  ans
}



