#' Plots to visualize the post-screening selection
#'
#' @description
#' This function constructs a sparsity vs. selection criterion curve for a \code{'selection'} object.
#'  When EBIC is used with voting, it also constructs a histogram showing the voting result.
#' @param x A \code{'selection'} object as the output from \code{\link{smle_select}()}.
#' @param ... Additional arguments to the \code{\link[base]{plot}()} function.
#' @method plot selection
#' @return
#' No return value.
#' 
#' @importFrom graphics plot.new
#' @examples
#' set.seed(1)
#' Data <- Gen_Data(correlation = "MA", family = "gaussian")
#' fit <- SMLE(Y = Data$Y, X = Data$X, k = 20, family = "gaussian")
#' fit_s <- smle_select(fit, vote = TRUE)
#' plot(fit_s)
#'
#' @export
#'
#'
plot.selection<-function(x,...){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  plot(x$criterion_value,xlab="Model sparsity", ylab= paste(x$criterion,"value"), cex.axis = 1.5, cex.lab = 1.5,...)
  if(x$vote ==TRUE){
    dev.new()
    y<-data.frame("Proportion"= sort(summary(x$ID_pool),decreasing = T)/length(x$gamma_seq))
    ID_names<- x$subset[as.numeric(names(summary(x$ID_pool)[order(summary(x$ID_pool),decreasing= T)]))]
    if(!is.null(x$data)){
      ID_names <- colnames(x$X)[ID_names]
    }
    barplot(y$Proportion,names.arg =ID_names ,
            xlab = "Candidate Features IDs",ylab="Features Voting Proportion",main="Voting results", 
            cex.axis = 1.5, cex.lab = 1.5)
  }
}
