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
#' from \code{\link{SMLE}()} or \code{\link{smle_select}()}. If omitted, the function returns the fitted 
#' response values based on the training data.
#'
#' @param ... 	Further arguments passed to \code{\link[stats]{predict.glm}()}.
#'
#' @return A prediction vector. The length of the vector equals to the number of observations of the data fitted in.
#' @rdname predict
#' 
#'
#' @export
#' @method predict smle
#' @examples
#' 
#' set.seed(1)
#' Data_sim <- Gen_Data(n = 420, p = 1000, sigma = 0.5, family = "gaussian")
#' train_X <- Data_sim$X[1:400,]; test_X <- Data_sim$X[401:420,]
#' train_Y <- Data_sim$Y[1:400]; test_Y <- Data_sim$Y[401:420]
#' fit1 <- SMLE(train_Y, train_X, family = "gaussian", k = 10)
#' 
#' #Fitted responses vs true responses in training data
#' predict(fit1)[1:10]
#' train_Y[1:10]
#' 
#' #Predicted responses vs true responses in testing data
#' predict(fit1, newdata = test_X)
#' test_Y
#' 
#' 
#' 
predict.smle<-function(object,newdata = NULL,...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_retained])
  
  if(!exists("type")){type = "response"}
  
  if( object$ctg ){
    # Categorical prediction
    
    CT <- object$codingtype
    
    if(CT == "all"){CT <- "standard"}
    
    data<- suppressWarnings(dummy.data.frame(data ,sep=".",codingtype = CT))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      
      
      # check ctg order for newdata
      
      CI<-(1:dim(object$X[,object$ID_retained])[2])[sapply(object$X[,object$ID_retained],is.factor)]
      
      CI_new<-(1:dim(newdata[,object$ID_retained])[2])[sapply(newdata[,object$ID_retained],is.factor)]
      
      if(!identical(CI, CI_new) ) stop("Testing features should have the same order with training features")
      
      # check ctg levels for newdata
      
      if(length(CI) == 1){   
        
        dum_col <- nlevels(object$X[,object$ID_retained][,CI])
      
        dum_col_new <- nlevels(newdata[,object$ID_retained][,CI])
        
      }else{
        
        dum_col <- sapply(object$X[,object$ID_retained][,CI],nlevels)
        
        dum_col_new <- sapply(newdata[,object$ID_retained][,CI],nlevels)
        
      }
      
      
      
      if(!identical(dum_col, dum_col_new) ) stop("Testing categorical features should have the same levels with training categroical features")
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_retained]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep=".",codingtype = CT))
      
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
predict.selection<-function(object,newdata = NULL, ...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  data = data.frame(Y = object$Y, X= object$X[,object$ID_selected])
  
  if(!exists("type")){type = "response"}
  
  if( object$ctg ){
    
    # Categorical prediction
    
    data<- suppressWarnings(dummy.data.frame(data ,sep=".",codingtype = CT))
    
    fit<-glm(Y~.,data = data ,family = family)
    
    if(is.null(newdata)){
      
      return(predict.glm(fit, newdata=NULL ,type = "response",...))
      
    }else{
      
      CT <- object$codingtype
      
      if(CT == "all"){CT <- "standard"}
      
      # check ctg order for newdata
      
      CI<-(1:dim(object$X[,object$ID_selected])[2])[sapply(object$X[,object$ID_selected],is.factor)]
      
      CI_new<-(1:dim(newdata[,object$ID_selected])[2])[sapply(newdata[,object$ID_selected],is.factor)]
      
      if(!identical(CI, CI_new) ) stop("Testing features should have the same order with training features")
      
      # check ctg levels for newdata
      
      if(length(CI) == 1){   
        
        dum_col <- nlevels(object$X[,object$ID_selected][,CI])
        
        dum_col_new <- nlevels(newdata[,object$ID_selected][,CI])
        
      }else{
        
        dum_col <- sapply(object$X[,object$ID_selected][,CI],nlevels)
        
        dum_col_new <- sapply(newdata[,object$ID_selected][,CI],nlevels)
        
      }
      
      if(!identical(dum_col, dum_col_new) ) stop("Testing categorical features should have the same levels with training categroical features")
      
      new_data = suppressWarnings(data.frame( X= newdata[,object$ID_selected]))
      
      newdata_dummy  <- suppressWarnings(dummy.data.frame(new_data ,sep=".",codingtype = CT))
      
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
