#' Extract and adjust voting from SMLE selection
#' 
#' 
#' When \code{\link{smle_select}()} is used with \code{criterion = "ebic"} and \code{vote = TRUE}, users 
#' can use \code{\link{vote_update}()} to adjust the voting threshold without a need 
#' of rerun \code{\link{smle_select}()}.
#' 
#' @param object A \code{'selection'} object as the output from \code{\link{smle_select}()}.
#' @param ... This argument is not used and listed for method consistency.
#' @return The function returns a vector indicating the features selected by
#'  EBIC voting with the specified \code{vote_threhold}.
#' @export
#' @examples 
#' set.seed(1)
#' Data <- Gen_Data(n = 100, p = 3000, correlation = "MA", rho = 0.7, family = "gaussian")
#' colnames(Data$X)<- paste("X.",seq(3000) , sep = "")
#' fit <- SMLE(Y = Data$Y, X = Data$X, k = 20, family = "gaussian")
#' fit_s <- smle_select(fit, criterion = "ebic", vote = TRUE)
#' plot(fit_s)
#' fit_s
#' vote_update(fit_s, vote_threshold = 0.4)

vote_update<-function(object, ...){
  UseMethod("vote_update")
}
#' @method vote_update selection
#' @rdname vote_update
#' @param vote_threshold A voting threshold in percentage. A feature is
#'  considered to be important when it receives votes passing the threshold.
#'  Default is 0.6.
#' @export
vote_update.selection<-function(object, vote_threshold = 0.6, ...){
  
    IP<-object$ID_pool
    
    IP_f<-summary(IP)[order(summary(IP),decreasing= T)]/max(summary(IP))
    
    ID_names<-object$subset[sort(as.numeric(names(IP_f[IP_f>=vote_threshold])))] 
    
    if(!is.null(object$data)){
      
      ID_names <- colnames(object$X)[ID_names]
      
    }
    
    return(ID_names)

}