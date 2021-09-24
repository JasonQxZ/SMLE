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
#' fit <- SMLE(Data$Y, Data$X, k = 20, family = "gaussian")
#' fit_s <- smle_select(fit, vote = TRUE)
#' plot(fit_s)
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
    ID_names<- x$subset[as.numeric(names(summary(x$ID_pool)[order(summary(x$ID_pool),decreasing= T)]))]
    barplot(y$Proportion,names.arg =ID_names ,
            xlab = "Candidate Features IDs",ylab="Features Voting Proportion",main="Voting results"
            )
  }
}
