#' Summary a selection object from smle_select
#' @description This function prints a summary of a \code{'selection'} object.
#' In particular, it gives the selected features along with their re-fitted model coefficients.
#' For reference, it also shows the values of the selection criterion used in selection for all candidate models.
#'
#' @importFrom utils capture.output
#'
#' @param object Fitted \code{'selection'} object.
#'
#' @param ... This argument is not used and listed for method consistency.
#'
#' @return
#' No return value, called for side effects.
#'
#' @export
#'
#'
#' @method summary selection
#'
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' E<-smle_select(fit)
#' summary(E)
#'
summary.selection<-function( object,...){

  if(object$vote != T){

    Description<-data.frame("Dim_of_Y :" = paste(c(dim(object$Y)[1],1),collapse = ' x '),
                            
                            "Dim_of_X :" = paste(dim(object$X),collapse = ' x '),
                            
                            "Model_Type :"= .simpleCap(object$family),
                            
                            "Retained_model_size :"= as.character(object$Num_Selected),
      
                            "Selection_criterion" = as.character(object$criterion),

                            "Selected_features" = paste("V",object$ID_Selected,sep='',collapse=' '))

    if(object$criterion=='ebic') {

      Description<-cbind(Description,"Gamma_for_ebic"=as.character(object$gamma_ebic))

    }

  }else{

    Description<-data.frame("Dim_of_Y :" = paste(c(dim(object$Y)[1],1),collapse = ' x '),
                            
                            "Dim_of_X :" = paste(dim(object$X),collapse = ' x '),
                            
                            "Model_Type :"= .simpleCap(object$family),
                            
                            "Gamma_candidate :"=  paste(object$gamma_seq,sep='',collapse=' '),
                            
                            "features_selected_by_voting"=paste("V",object$ID_Voted,sep='',collapse=' '))
                        
  }


  message(paste0(capture.output(t(Description)), collapse = "\n"))
}
