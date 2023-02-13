#' Plots to visualize SMLE screening
#'
#' This function returns two plot windows. By default, the first shows 
#' 1) the solution path (estimated coefficient by iteration step) for 
#' the retained features.
#' By default, the second plot contains 4 plots to assess convergence:
#' 2) log-likelihood,  
#' 3) Euclidean distance between the current and the previous coefficient estimates,   
#' 4) the number of tries in u-search (see details of \code{\link{SMLE}()}),  
#' and 5) the number of features changed in the current active set.
#'
#' @param x A \code{'smle'} object as the output from \code{\link{SMLE}()}.
#' 
#' @param num_path The number of top coefficients to be shown.
#' Default is equal to the number of features retained in the model.
#' @param label Logical flag for whether to label each curve with the feature index. Default is \code{TRUE}.
#' @param which_path A vector to control which features are shown in addition to the paths for the most significant coefficients.
#' @param out_plot A number from 1 to 5 indicating which plot is to be shown in the separate window; the default for solution path plot is "1".
#' See Description for plot labels 2-5.
#' 
#' @param ... Additional arguments passed to the second plot.
#' @return
#' No return value.
#' @export
#' @importFrom graphics plot.new
#' @importFrom graphics text
#' @method plot smle
#' @examples
#' \donttest{
#' set.seed(1)
#' Data <- Gen_Data(correlation = "CS")
#' fit <- SMLE(Y = Data$Y,X = Data$X, k = 20, family = "gaussian")
#' plot(fit)
#'}

plot.smle<-function(x,num_path=NULL,label = TRUE,which_path=NULL,out_plot=1,...){
  
  ## Save old plot settings before changing the settings
  
  old.par <- par(no.readonly = TRUE) # all par settings which could be changed.

  nsteps<-x$steps

  Feature_path<-x$path_retained
  
  if(!is.null(x$data)){
    
    # If formula is provided, extract the actual indices by name.
  
    x$ID_retained<- sort((1:length(colnames(x$X)))[colnames(x$X) %in% x$iteration_data$feature_name])
  
  }
  
  plt1<- function(){
    
    if(is.null(num_path)){num_path=length(x$ID_retained)}
    
    TOP_index<- x$ID_retained[sort(abs(x$coef_retained),decreasing = T,index.return=T)$ix][1:num_path]
    
    TOP_index<-unique(c(TOP_index,which_path))
    
    Feature_path<-Feature_path[TOP_index,]
    
    invisible(plot(NULL, xlim=c(1,nsteps+1 + 0.1*nsteps), ylim=c(floor(min(Feature_path[,nsteps])),ceiling(max(Feature_path[,nsteps]))),xlab="steps",ylab="coefficients",main = "Solution path",...))
    
    lines(rep(0,nsteps+1),lty=1,lwd=1,col="black")
    
    A<-lapply(1:length(TOP_index),function(i){lines(Feature_path[i,],lty=i,col=rainbow(length(TOP_index))[i])})
    
    if(label){
      nnz=length(TOP_index)
      xpos=nsteps+1
      pos=4
      xpos=rep(xpos+0.02*xpos,nnz)
      ypos=Feature_path[,nsteps+1]
      text(xpos,ypos,paste(TOP_index),pos=pos,...)
    }
  }
  
  plt2 <- function(){
    
    plot(y=x$likelihood_iter,x=1:nsteps,xlab="steps",ylab="log-likelihood",type="b", main= "Likelihood convergence",...)
    
  }
  
  plt3 <- function(){
    
    plot(x$coef_dist, x=1:length(x$coef_dist),xlab="steps",ylab="L2-dist b/w beta updates",
         type="b" ,main ="Coefficient convergence",...)

  }
  
  plt4<-function(){
    
    plot(y=x$Usearch,x=1:nsteps,xlab="steps",ylab="No. of tries",main ="U-search",...)
    
  }
  
  plt5 <- function(){
    
    invisible(plot(y=x$iteration_data$FD,x=1:nsteps,xlab="steps",ylab="No. of changes",main = "Retained feature change",...))

  }
  
  plots <- paste0('plt',seq(1,5),'()')
  
  out_plot <- paste0('plt',out_plot,'()')
  
  #--------------------plot-------------
  
  eval(parse(text = out_plot))

  par(mfrow=c(2,2))

  lapply(plots[plots != out_plot],function(p){eval(parse(text = p))})
  
  par(old.par)

 
}
  