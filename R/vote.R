#' Extract and adjust voting from selection model
#' 
#' This function is a method for \code{'selection'} object from \code{smle_select()}.
#' It extracts voting result for \code{'selection'} and can be used to change the voting percentage threshold.
#' @param object Object of class \code{'selection'}.
#' @param ... This argument is not used and listed for method consistency
#' @return Vector containing the features selected.
#' @export
#' @examples 
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' E<-smle_select(fit,vote=T)
#' plot(E)
#' vote(E,vote_threshold = 0.8)

vote<-function(object, ...){
  UseMethod("vote")
}
#' @method vote selection
#' @rdname vote
#' @param vote_threshold A relative voting threshold in percentage. A feature is
#'  considered to be important when it receives votes passing the threshold. Default is 0.8.
#' @export
vote.selection<-function(object,vote_threshold=NULL,...){
  
  if(is.null(vote_threshold)){
    
    return(object$ID_voted)
    
  }else{
    IP<-object$ID_pool
    
    IP_f<-summary(IP)[order(summary(IP),decreasing= T)]/max(summary(IP))
    
    ID_voted<-as.numeric(names(IP_f[IP_f>=vote_threshold]))  
    
    return(object$subset[ID_voted])
    
    }
}