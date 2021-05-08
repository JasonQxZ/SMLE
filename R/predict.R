#' Prediction based on SMLE screening and selection
#' 
#' @description
#' This function makes prediction for the response in new data based on the 
#' features retained in \code{'smle'} object.
#' 
#' @importFrom stats binomial gaussian glm.fit poisson
#'
#' @param object A fitted object of class \code{'smle'} as the output from
#' SMLE.
#'
#' @param newdata Matrix of new values for x at which predictions are to be made,
#' without the intercept term. If omitted, the function returns the fitted 
#' response values based on the training data and the retained features in \code{'smle'}.
#'
#' @param ... 	Further arguments pass to \code{predict.glm()}.
#'
#' @return The object returned depends on type.
#' @rdname predict
#' 
#'
#' @export
#' @method predict smle
#' @examples
#'
#' set.seed(123.456)
#'
#' Data_sim<-Gen_Data(n= 200, p =1000, correlation="AR",family = "gaussian", num_ctgidx =5)
#'
#' fit<-SMLE(Data_sim$Y,Data_sim$X, family = "gaussian" , k=10)
#'
#' predict(fit,newdata= Data_sim$X[10:20,])
#'
predict.smle<-function(object,newdata = NULL,...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_Retained])
  
  if(!exists("type")){type = "response"}
  
  if( object$ctg ){
    # Categorical prediction
    
    data<- suppressWarnings(dummy.data.frame(data ,sep="."))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = type,...))
      
    }else{
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_Retained]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep="."))
      
      return(predict.glm(fit, newdata=new_data ,type = type,...))
      
    }
    

    
  }else{
    
    # Numerical prediction
    data = data.frame(Y = object$Y, X= object$X[,object$ID_Retained])
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
        
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_Retained]))
      
      return(predict.glm(fit, newdata=new_data ,type = "response",...))
      
      }
       

  }
}
#' @rdname predict
#'
#' @export
#' @method predict selection
#' @examples
#'
#'
predict.selection<-function(object,newdata = NULL,...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_Selected])
  
  if(!exists("type")){type = "response"}
  
  if( object$ctg ){
    # Categorical prediction
    
    data<- suppressWarnings(dummy.data.frame(data ,sep="."))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = type,...))
      
    }else{
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_Selected]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep="."))
      
      return(predict.glm(fit, newdata=new_data ,type = type,...))
      
    }
    
    
    
  }else{
    
    # Numerical prediction
    data = data.frame(Y = object$Y, X= object$X[,object$ID_Selected])
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_Selected]))
      
      return(predict.glm(fit, newdata=new_data ,type = "response",...))
      
    }
    
    
  }
}





  
  

