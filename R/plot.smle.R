#' Plots to visualize the SMLE screening step
#'
#' This function returns two plot windows. By default, the first contains 4 plots to assess:
#' 1) log-likelihood,  2) Euclidean distance between the current
#' and the previous coefficient estimates,   3)  the number of tries in tuning parameter "u" in IHT
#' algorithm (see details of \code{SMLE()}),  and 4) the number of features changed in the current active set.
#' By default, the second plot shows the solution path (estimated coefficient by iteration step) for selected features.
#'
#' @param x Fitted  \code{"smle"} object from SMLE.
#' 
#' @param num_path Number of top coefficients to be shown in solution path plot.
#' Default in solution path plot is equal to the number of features retained in the model.
#' @param which_path A vector to control which features are shown in addition to the paths for the most significant coefficients.
#' @param out_plot A number from 1 to 5 indicating which plot is to be shown in the separate window; the default for solution path plot is "5".
#' See Description for plot labels 1-4.
#' @param ... Additional arguments to the second plot.
#' @return
#' The function returns a vector of the predicted or fitted response values depending on the specified arguments.
#' @export
#' @importFrom graphics plot.new
#' @method plot smle
#' @examples
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' plot(fit)
#'

plot.smle<-function(x,num_path=NULL,which_path=NULL,out_plot=5,...){
  
  ## new plot
  plot.new()
  oldpar <- par(no.readonly = TRUE)

  on.exit(par(oldpar))

  #-------------------check
  nsteps<-x$steps

  Feature_path<-x$path_retained

  #--------------------plot-------------
  if(out_plot==5){

  par(mfrow=c(2,2))

  plot(y=x$likelihood_iter,x=1:nsteps,xlab="steps",ylab="log-likelihood",type="b")

  title("Likelihood convergence")

  plot(x$coef_dist, x=1:length(x$coef_dist),xlab="steps",ylab="L2-dist b/w beta updates",
       main=" Coefficient convergence",type="b" )

  plot(y=x$usearch,x=1:nsteps,xlab="steps",ylab="No. of tries")

  title("U-search")

  plot(y=x$FD,x=1:nsteps,xlab="steps",ylab="No. of changes")

  title("Retained feature change")

  dev.new()

  if(is.null(num_path)){num_path=x$k}

  TOP_index<- x$ID_retained[sort(abs(x$coef_retained),decreasing = T,index.return=T)$ix][1:num_path]

  TOP_index<-unique(c(TOP_index,which_path))

  Feature_path<-Feature_path[TOP_index,]

  plot(NULL, xlim=c(1,nsteps+1), ylim=c(floor(min(Feature_path[,nsteps])),ceiling(max(Feature_path[,nsteps]))),xlab="steps",ylab="coefficients",...)

  title("Solution path")

  lines(rep(0,nsteps),lty=1,lwd=1,col="black")
  
  A<-lapply(1:length(TOP_index),function(i){lines(Feature_path[i,],lty=i,col=rainbow(length(TOP_index))[i])})

  legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
           legend=TOP_index,bty="n")
    
  }else if(out_plot==1){

    par(mfrow=c(2,2))
    
    if(is.null(num_path)){num_path=5}
    
    TOP_value<-rep(0,num_path)
    
    TOP_value<-sort(abs(Feature_path[,ncol(Feature_path)]),decreasing = T)[1:num_path]
    
    TOP_index<-(1:dim(x$X)[2])[Feature_path[,ncol(Feature_path)]%in% unique(c(TOP_value,-TOP_value))]
    
    TOP_index<-unique(c(TOP_index,which_path))
    
    Feature_path<-Feature_path[TOP_index,]
    
    plot(NULL, xlim=c(1,nsteps+1), ylim=c(floor(min(Feature_path[,nsteps])),ceiling(max(Feature_path[,nsteps]))),xlab="steps",ylab="coefficients",...)
    
    title("Solution path")
    
    lines(rep(0,nsteps),lty=1,lwd=1,col="black")
    
    A<-lapply(1:length(TOP_index),function(i){lines(Feature_path[i,],lty=i,col=rainbow(length(TOP_index))[i])})
    
    legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
           legend=TOP_index,bty="n")
    
    plot(x$coef_dist, x=1:length(x$coef_dist),xlab="steps",ylab="L2-dist b/w beta updates",
         main=" Coefficient convergence",type="b" )
    title(' Coefficient convergence')
    plot(y=x$usearch,x=1:nsteps,xlab="steps",ylab="No. of tries")
    title("U-search")
    #
    plot(y=x$FD,x=1:nsteps,xlab="steps",ylab="No. of changes")
    title("Retained feature change")

  dev.new()
  plot(y=x$likelihood_iter,x=1:nsteps,xlab="steps",ylab="Soulution path",type="b",...)
  title("Likelihood convergence")
  }else if(out_plot==2){
  
    par(mfrow=c(2,2))
    
    if(is.null(num_path)){num_path=5}
    
    TOP_value<-rep(0,num_path)
    
    TOP_value<-sort(abs(Feature_path[,ncol(Feature_path)]),decreasing = T)[1:num_path]
    
    TOP_index<-(1:dim(x$X)[2])[Feature_path[,ncol(Feature_path)]%in% unique(c(TOP_value,-TOP_value))]
    
    TOP_index<-unique(c(TOP_index,which_path))
    
    Feature_path<-Feature_path[TOP_index,]
    
    plot(NULL, xlim=c(1,nsteps+1), ylim=c(floor(min(Feature_path[,nsteps])),ceiling(max(Feature_path[,nsteps]))),xlab="steps",ylab="coefficients",...)
    
    title("Solution path")
    
    lines(rep(0,nsteps),lty=1,lwd=1,col="black")
    
    A<-lapply(1:length(TOP_index),function(i){lines(Feature_path[i,],lty=i,col=rainbow(length(TOP_index))[i])})
    
    legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
           legend=TOP_index,bty="n")
    
    plot(y=x$likelihood_iter,x=1:nsteps,xlab="steps",ylab="log-likelihood",type="b")
    title("Likelihood convergence")
    plot(y=x$usearch,x=1:nsteps,xlab="steps",ylab="No. of tries")
    title("U-search")
    plot(y=x$FD,x=1:nsteps,xlab="steps",ylab="No. of changes")
    title("Retained feature change")

    dev.new()
    plot(x$coef_dist, x=1:length(x$coef_dist),xlab="steps",ylab="L2-dist b/w beta updates",
         main=" Coefficient convergence",type="b" )
    title(' Coefficient convergence')
  
    }else if(out_plot==3){
    
    par(mfrow=c(2,2))
      
      if(is.null(num_path)){num_path=5}
      
      TOP_value<-rep(0,num_path)
      
      TOP_value<-sort(abs(Feature_path[,ncol(Feature_path)]),decreasing = T)[1:num_path]
      
      TOP_index<-(1:dim(x$X)[2])[Feature_path[,ncol(Feature_path)]%in% unique(c(TOP_value,-TOP_value))]
      
      TOP_index<-unique(c(TOP_index,which_path))
      
      Feature_path<-Feature_path[TOP_index,]
      
      plot(NULL, xlim=c(1,nsteps+1), ylim=c(floor(min(Feature_path[,nsteps])),ceiling(max(Feature_path[,nsteps]))),xlab="steps",ylab="coefficients",...)
      
      title("Solution path")
      
      lines(rep(0,nsteps),lty=1,lwd=1,col="black")
      
      A<-lapply(1:length(TOP_index),function(i){lines(Feature_path[i,],lty=i,col=rainbow(length(TOP_index))[i])})
      
      legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
             legend=TOP_index,bty="n")
      title('Coefficient convergence')
    plot(y=x$FD,x=1:nsteps,xlab="steps",ylab="No. of changes")
    title("Retained feature change")
    dev.new()
    plot(y=x$usearch,x=1:nsteps,xlab="steps",ylab="Soulution path",...)
    title("U-search")
  }else if(out_plot==4){
    par(mfrow=c(2,2))
    
    if(is.null(num_path)){num_path=5}
    
    TOP_value<-rep(0,num_path)
    
    TOP_value<-sort(abs(Feature_path[,ncol(Feature_path)]),decreasing = T)[1:num_path]
    
    TOP_index<-(1:dim(x$X)[2])[Feature_path[,ncol(Feature_path)]%in% unique(c(TOP_value,-TOP_value))]
    
    TOP_index<-unique(c(TOP_index,which_path))
    
    Feature_path<-Feature_path[TOP_index,]
    
    plot(NULL, xlim=c(1,nsteps+1), ylim=c(floor(min(Feature_path[,nsteps])),ceiling(max(Feature_path[,nsteps]))),xlab="steps",ylab="coefficients",...)
    
    title("Solution path")
    
    lines(rep(0,nsteps),lty=1,lwd=1,col="black")
    
    A<-lapply(1:length(TOP_index),function(i){lines(Feature_path[i,],lty=i,col=rainbow(length(TOP_index))[i])})
    
    legend("topright",lty=1:length(TOP_index),cex=0.5,col=rainbow(length(TOP_index)),
           legend=TOP_index,bty="n")
    

    plot(y=x$likelihood_iter,x=1:nsteps,xlab="steps",ylab="log-likelihood",type="b")
    title("Likelihood convergence")
    plot(x$coef_dist, x=1:length(x$coef_dist),xlab="steps",ylab="L2-dist b/w beta updates",
         main=" Coefficient convergence",type="b" )
    title(' Coefficient convergence')
    plot(y=x$usearch,x=1:nsteps,xlab="steps",ylab="No. of tries")
    title("U-search")
    dev.new()
    plot(y=x$FD,x=1:nsteps,xlab="steps",ylab="No. of changes")
    title("Retained feature change")
  }
}

