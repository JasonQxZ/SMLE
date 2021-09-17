#' Extract and adjust voting from a \code{'selection'} object
#' 
#' 
#' When \code{'selection'} is used with \code{criterion="ebic"} and \code{vote=TRUE}, 
#' users can use \code{vote_update} to adjust the voting threshold without a need 
#' of rerun \code{smle_select}.
#' 
#' @param object A \code{'selection'} object as the output from \code{\link{smle_select}}.
#' @param ... This argument is not used and listed for method consistency
#' @return The function returns a vector indicating the features selected by
#'  EBIC voting with the specified \code{vote_threhold}.
#' @export
#' @examples 
#' set.seed(1)
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' fit_s<-smle_select(fit,vote=TRUE)
#' plot(fit_s)
#' vote_update(fit_s,vote_threshold = 0.3)

vote_update<-function(object, ...){
  UseMethod("vote_update")
}
#' @method vote_update selection
#' @rdname vote_update
#' @param vote_threshold A voting threshold in percentage. A feature is
#'  considered to be important when it receives votes passing the threshold.
#'  Default is 0.6.
#' @export
vote_update.selection<-function(object,vote_threshold=NULL,...){
  
  if(is.null(vote_threshold)){
    
    return(object$ID_voted)
    
  }else{
    IP<-object$ID_pool
    
    IP_f<-summary(IP)[order(summary(IP),decreasing= T)]/max(summary(IP))
    
    ID_voted<-as.numeric(names(IP_f[IP_f>=vote_threshold]))  
    
    return(object$subset[ID_voted])
    
    }
}