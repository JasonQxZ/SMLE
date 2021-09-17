#' Prediction based on SMLE screening and selection
#' 
#' @description
#' This function returns the predicted response values for a fitted model of class \code{'smle'} or \code{'selection'}.
#' 
#' @importFrom stats binomial gaussian glm.fit poisson
#'
#' @param object A \code{'smle'} or \code{'selection'} object.
#'
#' @param newdata Matrix of new values for the features at which predictions are to be made using the final model 
#' from \code{\link{SMLE}} or \code{\link{smle_select}}. If omitted, the function returns the fitted 
#' response values based on the training data.
#'
#' @param ... 	Further arguments passed to \code{\link[stats]{predict.glm}}.
#'
#' @return A prediction vector. The length of the vector equals to the number of observations of the data fitted in.
#' @rdname predict
#' 
#'
#' @export
#' @method predict smle
#' @examples
#' set.seed(1)
#' Data_sim1<-Gen_Data(n= 200, p =1000, correlation="AR",family = "gaussian", num_ctgidx =5)
#' fit1<-SMLE(Data_sim1$Y,Data_sim1$X, family = "gaussian" , k=10)
#' predict(fit1,newdata= Data_sim1$X[10:20,])
#' fit1_s <- smle_select(fit1)
#' predict(fit1_s)
#' 
#' 
#' #Prediction example with influential categorical feature. 
#' set.seed(2)
#' Data_sim2<-Gen_Data(n= 220, p =1000, correlation="AR",family = "gaussian", num_ctgidx =5)
#' fit2<-SMLE(Data_sim2$Y[1:200],Data_sim2$X[1:200,], family = "gaussian" , k=10)
#' predict(fit2)
#' predict(fit2,newdata= Data_sim2$X[201:220,])
#' 
#' 
#' 
predict.smle<-function(object,newdata = NULL,...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_retained])
  
  if(!exists("type")){type = "response"}
  
  if( object$ctg ){
    # Categorical prediction
    
    data<- suppressWarnings(dummy.data.frame(data ,sep=".",codingtype = object$codingtype))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_retained]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep=".",codingtype = object$codingtype))
      
      return(predict.glm(fit, newdata=newdata_dummy ,type = "response",...))
      
    }
    
    
    
  }else{
    
    # Numerical prediction
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      new_data = data.frame( X= newdata[,object$ID_retained])
      
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
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_selected])
  
  if(!exists("type")){type = "response"}
  
  if( object$ctg ){
    # Categorical prediction
    
    data<- suppressWarnings(dummy.data.frame(data ,sep=".",codingtype = object$codingtype))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_selected]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep=".",codingtype = object$codingtype))
      
      return(predict.glm(fit, newdata=newdata_dummy ,type = "response",...))
      
    }
    
    
    
  }else{
    
    # Numerical prediction
    data = data.frame(Y = object$Y, X= object$X[,object$ID_selected])
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_selected]))
      
      return(predict.glm(fit, newdata=new_data ,type = "response",...))
      
    }
  }
}
