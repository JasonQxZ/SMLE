#' Extract coefficients from fitted model
#' 
#' Extract coefficients from fitted model for either a \code{'smle'} or \code{'selection'} object.
#' 
#' @param object Returned object from either the function \code{\link{SMLE}()} or \code{\link{smle_select}()}.
#' @param refit A logical flag that controls what coefficients are being return. Default is \code{TRUE}.
#' @param ... This argument is not used and listed for method consistency.
#' @return Fitted coefficients based on the screened or selected model specified in the object. 
#' If \code{refit=TRUE}, the coefficients are estimated by re-fitting the final 
#' screened/selected model with \code{\link{glm}()}. If \code{refit=FALSE} the coefficients estimated by the IHT algorithm are returned.  
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
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"= binomial(), "poisson" = poisson())
  
  if(refit == FALSE){    
    
    COEF<- c(object$intercept, object$coef_retained)
    
    names(COEF)[1]<- '(intercept)'
    
    return(COEF)
    
  }else{
    
    if(is.null(object$data)){  
    
    model = object$X[,object$ID_retained]
    
    if(is.null(colnames(model))){ colnames(model) <- object$ID_retained }
    
    data = data.frame( Y = object$Y, model)
    
    names(data) <- c("Y",colnames(model))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    return(coef(fit))
    
    }else{
      
      model = object$data[,object$ID_retained]

      data = data.frame( Y = object$Y, model)
      
      names(data) <- c("Y",names(object$data)[object$ID_retained])
      
      fit<-glm(Y~.,data = data ,family = family)
      
      return(coef(fit))
      }
    }
   }

  
  
         
        
#' @rdname coef
#' @method coef selection
#' @export
coef.selection<-function(object,refit = TRUE,...)
{  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  if(object$vote == TRUE){ ID <- object$ID_voted}else{ID<- object$ID_selected}
  
  if(refit == FALSE){
    
    COEF<- c(object$intercept, object$coef_selected)
    
    names(COEF)[1]<- '(intercept)'
    
    return(COEF)
    
    
    }else{
      
      
        if(is.null(object$data)){  
          
          model = object$X[,ID]
          
          if(is.null(colnames(model))){ colnames(model) <- ID }
          
          data = data.frame(Y = object$Y, model)
          
          names(data) <- c("Y",colnames(model))
          
          fit<-glm(Y~.,data = data ,family = family)
          
          return(coef(fit))
          
        }else{
          
          model = object$data[,ID]
          
          data = data.frame( Y = object$Y, model)
          
          names(data) <- c("Y", names(object$data)[ID])
          
          fit<-glm(Y~.,data = data ,family = family)
          
          return(coef(fit))
        }
      
  }
}
