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
  
  plt1<-function(){
    
    invisible(plot( x= seq(from =x$k_min, to =x$k_max),y = x$criterion_value,xlab="number of features in model", ylab= paste(x$criterion,"value"),...))
    
  }
  
  plt2<-function(){
    
    y<-data.frame("Proportion"= sort(summary(x$ID_pool),decreasing = T)/length(x$gamma_seq))
    ID_names<- x$subset[as.numeric(names(summary(x$ID_pool)[order(summary(x$ID_pool),decreasing= T)]))]
    if(!is.null(x$data)){
      ID_names <- colnames(x$X)[ID_names]
    }
    invisible(barplot(y$Proportion,names.arg =ID_names ,
            xlab = "candidate features IDs",ylab="features voting proportion",main="voting results", 
            ...))
  }
  

  if(x$vote ==TRUE){
    plt1()
    plt2()
  }else{

    plt1()
    
  }
  
  
 
}
