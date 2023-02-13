#' Extract coefficients from fitted model
#' 
#' Extract coefficients from fitted model for either a \code{'smle'} or \code{'selection'} object.
#' 
#' @param object Returned object from either the function \code{\link{SMLE}()} or \code{\link{smle_select}()}.
#' @param refit A logical flag that controls what coefficients are being return. Default is \code{TRUE}.
#' @param ... This argument is not used and listed for method consistency.
#' @return Fitted coefficients based on the screened or selected model specified in the object. 
#' If \code{refit = TRUE}, the coefficients are estimated by re-fitting the final 
#' screened/selected model with \code{\link{glm}()}. If \code{refit = FALSE} the coefficients estimated by the IHT algorithm are returned.  
#' 
#' @rdname coef
#' @method coef smle
#' @examples
#' 
#' set.seed(1)
#' Data<-Gen_Data(n=100, p=5000, family = "gaussian", correlation="ID")
#' fit<-SMLE(Y = Data$Y, X = Data$X, k=15, family = "gaussian")
#' coef(fit)
#' fit_s<-smle_select(fit)
#' coef(fit_s)

#' @export
coef.smle<-function(object,refit = TRUE, ...){
  
  intercept <- object$intercept
  
  coefficients <- object$coef_retained
  
  data <- object$data
  
  ID <- object$ID_retained
  
  if(refit == FALSE){    
    
    COEF<- c(intercept, coefficients)
    
    if(is.null(names(coefficients))){names(coefficients) <- paste0("`",ID,"`")}
    
    names(COEF)<- c('(intercept)',names(coefficients))
    
    return(COEF)
    
  }else{
    
    if(is.null(data)){  
      
      # This case for matrix/data.frame input. ID are indices in data.
      
      model = object$X[,ID]
    
      if(is.null(colnames(model))){ colnames(model) <- ID }
    
      data = data.frame( Y = object$Y, model)
  
      names(data) <- c("Y",colnames(model))

      }else{
        
        # This case for formula input. ID are indices in data.
      
        model = object$data[,ID]

        data = data.frame( Y = object$Y, model)
      
        names(data) <- c("Y",names(model))
      
      }
    
    fit<-glm(Y~.,data = data ,family = object$family)
    
    return(coef(fit))
    
    }
   }

  
  
         
        
#' @rdname coef
#' @method coef selection
#' @export
coef.selection<-function(object, refit = TRUE, ...)
{  

  if(object$vote == TRUE){ 
    
    object$ID_retained <- object$ID_voted
    
  }else{
      
    object$ID_retained <- object$ID_selected
    
  }
  
  object$coef_retained <- object$coef_selected
  
  class(object) <- 'smle' 
 
  return(coef(object, refit = refit ,...))
}
