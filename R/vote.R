#' Extract and adjust voting from selection model
#' 
#' This function is a method for \code{'selection'} object from \code{smle_select()}.
#' It extract voting result for \code{'selection'} and is able to change the threshold.
#' @param object Object of class \code{'selection'}.
#' @param ... keep for consistency.
#' @return Vector containing the features selected
#' @export
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