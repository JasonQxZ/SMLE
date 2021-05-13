#' Plots to visualize the selection
#'
#' @description
#' This function constructs a sparsity vs. selection criterion curve for a selection object.
#'  When EBIC is used with voting, it also constructs a histogram showing the voting result.
#' @param x Fitted \code{'selection'} object from \code{smle_select()}.
#' @param ... Additional arguments to the plot function.
#' @method plot selection
#' @return
#' No return value.
#' 
#' @importFrom graphics plot.new
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' E<-smle_select(fit)
#' #Then E is a object of "selection"
#' plot(E)
#'
#' @export
#'
#'
plot.selection<-function(x,...){
  ## new plot
  plot.new()
  plot(x$criterion_value,xlab="Model sparsity", ylab= paste(x$criterion,"value"),...)
  if(x$vote ==TRUE){

    dev.new()
    percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits), "%")
    }
    y<-data.frame("Proportion"= sort(summary(x$ID_pool),decreasing = T)/length(x$gamma_seq))
    ID_names<- x$sub_model[as.numeric(names(summary(x$ID_pool)[order(summary(x$ID_pool),decreasing= T)]))]
    barplot(y$proportion,names.arg =ID_names ,
            xlab = "Candidate Features IDs",ylab="Featrues Voting Proportion",main="Voting results"
            )
  }
}
