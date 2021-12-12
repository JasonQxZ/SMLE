#' Prediction based on SMLE screening and selection
#' 
#' @description
#' For a model object of class \code{'smle'} or \code{'selection'}, this function returns the predicted response values after re-fitting the model using \code{\link{glm}}.
#' 
#' @importFrom stats binomial gaussian glm.fit poisson
#'
#' @param object A \code{'smle'} or \code{'selection'} object.
#'
#' @param newdata Matrix of new values for the features at which predictions are to be made. If omitted, the fitted linear predictors are used.
#' 
#' @param type The type of prediction required by \code{\link[stats]{predict.glm}()}.
#'
#' @param ... 	Further arguments passed to \code{\link[stats]{predict.glm}()}.
#'
#' @return A prediction vector with length equal to the number of rows of \code{newdata}.
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
#' fit1 <- SMLE(Y = train_Y, X = train_X, family = "gaussian", k = 10)
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
predict.smle<-function(object,newdata = NULL, type = c("link", "response", "terms"), ...){
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson" = poisson())
  
  if(!is.null(object$data)){
    
    model = object$data[,object$ID_retained]
    
    data = data.frame( Y = object$Y, model)
    
    names(data) <- c("Y",names(object$data)[object$ID_retained])
    
    if(!is.null(newdata)){
      
      newdata<- newdata[(1:dim(newdata)[2])[names(newdata) %in% object$iteration_data$feature_name]]
      
      }
     
    fit<-glm(Y~.,data = data ,family = family)
    
    return(predict.glm(fit, newdata=newdata ,type = type,...))
    
  }else{
    
    model = object$X[,object$ID_retained]
    
    if(is.null(colnames(model))){ colnames(model) <- object$ID_retained }
    
    data = data.frame( Y = object$Y, model)
    
    names(data) <- c("Y",colnames(model))
    
    if(!is.null(newdata)){
      
      if(!is.null(colnames(newdata))){
        
        newdata<- newdata[(1:dim(newdata)[2])[names(newdata) %in% colnames(object$X[,object$ID_retained])]]
        
      }else{
        
        newdata <- data.frame(newdata[,object$ID_retained])
        
        names(newdata)<-names(data)[-1]
        
      }
     }
    
    fit<-glm(Y~.,data = data ,family = family)
    
    return(predict.glm(fit, newdata=newdata ,type = type,...))
    }
  }
#' @rdname predict
#'
#' @export
#' @method predict selection
#' @examples
#'
#'
predict.selection<-function(object,newdata = NULL,type = c("link", "response", "terms"), ...){
  
  if(object$vote == TRUE){ ID <- object$ID_voted}else{ID<- object$ID_selected}
  
  family<-switch(object$family, "gaussian" = gaussian(),  "binomial"=binomial(), "poisson"=poisson())
  
  if(!is.null(object$data)){
    
    model = object$data[,ID]
    
    data = data.frame( Y = object$Y, model)
    
    names(data) <- c("Y",names(object$data)[ID])
    
    if(!is.null(newdata)){
      
      newdata<- newdata[(1:dim(newdata)[2])[names(newdata) %in% names(object$data[ID])]]

    }
    
    fit<-glm(Y~.,data = data ,family = family)
    
    return(predict.glm(fit, newdata=newdata ,type = type,...))
    
  }else{
    
    model = object$X[,ID]
    
    if(is.null(colnames(model))){ colnames(model) <- ID }
    
    data = data.frame( Y = object$Y, model)
    
    names(data) <- c("Y",colnames(model))
    
    if(!is.null(newdata)){
      
      if(!is.null(colnames(newdata))){
        
        newdata<- newdata[(1:dim(newdata)[2])[names(newdata) %in% colnames(object$X[,ID])]]
        
      }else{     
        
        newdata <- data.frame(newdata[,object$ID])
      
        names(newdata)<-names(data)[-1]
        
        
      }
      


    }
    
    fit<-glm(Y~.,data = data ,family = family)
    
    return(predict.glm(fit, newdata=newdata ,type = type,...))
  }
}
    

