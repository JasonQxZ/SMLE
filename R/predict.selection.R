#' Prediction based on SMLE selection
#' 
#' @description
#' This function makes prediction for the response in new data based on the 
#' features selected by \code{'selection'} object.
#' 
#' @importFrom stats binomial gaussian glm.fit poisson predict.glm
#'
#'
#' @param object A fitted object of class \code{'selection'} as the output 
#' from \code{smle_select()}.
#'
#'
#' @param newdata Matrix of new values for x at which predictions are to be made,
#' without the intercept term. If omitted, the function returns the fitted 
#' response values based on the training data and the retained features in \code{'selection'}.
#'
#'
#' @param type  The type of prediction required. Type "link" leads to a prediction 
#' on the scale of linear predictors (i.e. linear combination of retained features); 
#' type \code{'response'} leads to a prediction on the mean value of the response.
#' When \code{'family=gaussian'} in \code{'selection'}, the two predication types are equivalent.
#'
#' @param ... 	Further arguments pass to \code{predict.glm()}.
#'
#' @return
#' 
#' The object returned depends on type.
#'
#' @export
#' @method predict selection
#' @examples
#'
#' set.seed(123.456)
#'
#' Data_sim<-Gen_Data(n= 200, p =1000, correlation="AR",family = "gaussian")
#'
#' fit<-SMLE(Data_sim$Y,Data_sim$X, family = "gaussian")
#
#' E<-smle_select(fit, tune="ebic")
#'
#' predict(E , type ="link")
#'
#' @export
predict.selection<-function(object,newdata = NULL,type = c("link", "response"),...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_Selected])
  
  if( object$Ctg ){
    # Categorical prediction
    
    data<- suppressWarnings(dummy.data.frame(data ,sep="."))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = type,...))
      
    }else{
      
      new_data = suppressWarnings(data.frame(Y = object$Y, X= newdata[,object$ID_Selected]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep="."))
      
      return(predict.glm(fit, newdata=new_data ,type = type,...))
      
    }
    
    
    
  }else{
    
    # Numerical prediction
    data = data.frame(Y = object$Y, X= object$X[,object$ID_Selected])
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = type,...))
      
    }else{
      
      new_data = suppressWarnings(data.frame(Y = object$Y, X= newdata[,object$ID_Selected]))
      
      return(predict.glm(fit, newdata=new_data ,type = type,...))
      
    }
    
    
  }
}

  

