#' Print an object
#' @description This function prints information about the fitted model from a call to \code{\link{SMLE}()} or \code{\link{smle_select}()},
#'  or about the simulated data from a call to \code{\link{Gen_Data}()}. The object passed as an argument to print is returned invisibly. 
#' @rdname print
#' @importFrom utils capture.output
#'
#' @param x Fitted object.
#'
#' @param ... This argument is not used and listed for method consistency.
#' 
#' @return Return argument invisibly.
#' 
#' @examples
#' set.seed(1)
#' Data<-Gen_Data(correlation = "MA", family = "gaussian")
#' Data
#' fit<-SMLE(Y = Data$Y, X = Data$X, k = 20, family = "gaussian")
#' print(fit)
#' summary(fit)
#' @export
#' @method print smle



print.smle<-function(x,...){
  
  catln <- function (...)  base::cat(..., "\n", sep = "")
  
  ## call
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Subset:")
  catln("  Model size: ", as.character(x$num_retained))
  if( !is.null(x$iteration_data$feature_name)){ catln("  Feature name: ", paste(x$iteration_data$feature_name, sep = "", collapse = ", "))}
  catln("  Feature index: ", paste(x$ID_retained, sep = "\n", collapse = ", "))
  

  ## done
  invisible(x)
  
}


#' @export
#' @method print selection

#' @rdname print
print.selection<-function(x,...){
  catln <- function (...)  base::cat(..., "\n", sep = "")
  
  ## call
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Subset:")
  if(is.null(x$data)){  
    if(x$vote==TRUE){
      if( !is.null(colnames(x$X))){ catln("  Feature name: ", paste(colnames(x$X)[x$ID_voted], sep = "", collapse = ", "))}
      catln("  Features selected by voting: ",  paste(x$ID_voted, sep = "\n", collapse = ", "))
    }else{
      if( !is.null(colnames(x$X))){ catln("  Feature name: ", paste(colnames(x$X)[x$ID_selected], sep = "", collapse = ", "))}  
      catln("  Selected feature index: ", paste(x$ID_selected, sep = "\n", collapse = ", "))}
    ## done
    return(invisible(x))
    
  }else{
    
    if(x$vote==TRUE){
      catln("  Feature name: ", paste(names(x$data)[x$ID_voted], sep = "", collapse = ", "))
      catln("  Features selected by voting: ",  paste(x$ID_voted, sep = "\n", collapse = ", "))
    }else{
      catln("  Feature name: ", paste(names(x$data)[x$ID_selected], sep = "", collapse = ", "))  
      catln("  Selected feature index: ", paste(x$ID_selected, sep = "\n", collapse = ", "))}
    ## done
    return(invisible(x))
    
    
    }

  
}

#' 
#' @export
#' @method print summary.smle
#' @rdname print

print.summary.smle <- function(x, ...){
  
  catln <- function (...)  base::cat(..., "\n", sep = "")
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Summary:")
  cat("\n")
  catln("  Length of response: " , x$DimY)
  catln("  Dim of features: " , paste(x$DimX,collapse = ' x '))
  catln("  Model type: " , x$family$family)
  catln("  Model size: ", x$size)
  
  if( !is.null(x$iteration_data$feature_name)){ catln("  Feature name: ", paste(x$iteration_data$feature_name, sep = "", collapse = ", "))}
  catln("  Feature index: ", paste(x$ID_retained, sep = "\n", collapse = ", "))
  if( !is.null(x$intercept)){ catln("  Intercept: ", x$intercept)}
  catln("  Coefficients estimated by IHT: ", paste(format(x$coef_retained, digits = 3,trim = TRUE),collapse = ', '))
  catln("  Number of IHT iteration steps: ",as.character(x$steps))
  
  
  invisible(x)
}

#' @export
#' @method print summary.selection
#' @rdname print
print.summary.selection <- function(x, ...){
  
  catln <- function (...)  base::cat(..., "\n", sep = "")
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Summary:")
  catln(" ")
  catln("  Length of response: " , x$DimY)
  catln("  Dim of features: " , paste(x$DimX,collapse = ' x '))
  catln("  Model type: " , x$family$family)
  catln("  Selected model size: ", x$size)
  if( !is.null(x$data)){ catln("  Feature name: ", paste(names(x$data)[x$ID_selected], sep = "", collapse = ", "))
    }else if( !is.null(colnames(x$X))){ catln("  Feature name: ", paste(colnames(x$X)[x$ID_selected], sep = "", collapse = ", "))}
  catln("  Selected feature index: ", paste(x$ID_selected, sep = "\n", collapse = ", "))
  catln("  Selection criterion: ", x$criterion)
  if(x$criterion=='ebic'){  catln("  Gamma for ebic: ", as.character(x$gamma_ebic))}
  if(x$vote==TRUE){  catln("  Features selected by voting: ",  paste(x$ID_voted, sep = "\n", collapse = ", "))}
  
  invisible(x)
}



#' @export
#' @method print sdata
#' @rdname print
print.sdata<-function(x,...){
  catln <- function (...)  base::cat(..., "\n", sep = "")

  ## call
  catln("Call:")
  catln(" ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Simulated Dataset Properties:")
  catln(" Length of response: " , length(x$Y))
  catln(" Dim of features: " , paste(dim(x$X),collapse = ' x '))
  catln(" Correlation: " , x$correlation)
  if(x$correlation != "independent"){catln(" Rho: ", x$rho)}
  catln(" Index of causal features: " , paste(x$subset_true,collapse = ', '))
  if(x$categorical ==T){catln(" Contain categorical features: TRUE" )}
  catln(" Model type: ", x$family)

  ## done
  invisible(x)
  
}
