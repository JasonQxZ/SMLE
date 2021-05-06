#' summary a SMLE object from SMLE
#'
#' @description This functions prints a summary of a SMLE object.
#' In particular, it shows the features retained after SMLE-screening
#' and the related convergence information.
#'
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

  cat("\nCall: ", deparse(object$call), "\n\n")

  Description<-data.frame("Dim_of_Y :" = paste(c(dim(object$X)[1],1),collapse = ' x '),

                          "Dim_of_X :" = paste(dim(object$X),collapse = ' x '),

                          "Model_Type :"= .simpleCap(object$family),

                          "Retained_model_size :"= as.character(object$Num_Retained),

                          "Retained_features :"=paste("V",object$ID_Retained,sep='',collapse=' '),

                          "Coefficients :" = paste(format(object$Coef_Retained, digits = 3),collapse = ' '),

                          "Number_of_steps :" = as.character(object$steps))

  if( !is.null(object$Intercept) ){

    Description <- cbind ( Description, "Intercept :" = as.character(format(object$Intercept, digits = 3)) )

    }

  if("ctg" %in% class(object)){
    
    if(length(object$I$CI)==1){
      
      Description <- cbind(Description,
                                               "Categorical_features :" =paste("C", object$I$CI, sep='',collapse = ' '      ),
                                               "Level_of_categories  :" =nlevels(object$I$CM[,object$I$CI])) 
    }else{    
      
      Description <- cbind(Description,
                                   "Categorical_features :" =paste("C", object$I$CI, sep='',collapse = ' '      ),
                                   "Level_of_categories  :" =paste(sapply(object$I$CM[,object$I$CI],nlevels) , collapse = ', '))
      }
    

    }

  message(paste0(capture.output(t(Description)), collapse = "\n"))

}

