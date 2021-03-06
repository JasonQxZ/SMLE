#' Print an object
#' @description This functions prints information about the fitted model from a call to SMLE() or smle_select(),
#'  or about the simulated data from a call to Gen_data(). The object passed as an argument to print is returned invisibly. 
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
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' Data
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' fit
#' summary(E)
#' E<-smle_select(fit,vote=T)
#' E
#' summary(E)
#' @export
#' @method print smle



print.smle<-function(x,...){
  catln <- function (...)  base::cat(..., "\n", sep = "")
  print <- function (..., skip = 0, indent = 0){
    output <- capture.output(base::print(...))
    if (skip > 0) output <- output[-seq_len(skip)]
    indent <- paste0(rep(" ", indent), collapse = "")
    cat(paste0(indent, output, "\n"), sep = "")
  }
  
  ## call
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Subset:")
  catln("  Model size: ", as.character(x$num_retained))
  if( !is.null(colnames(x$X))){ catln("  Feature Name: ", paste(colnames(x$X)[x$ID_retained], sep = " ", collapse = ","))}
  catln("  Feature Index: ", paste(x$ID_retained, sep = "\n", collapse = ","))
  

  ## done
  invisible(x)
  
}


#' @export
#' @method print selection

#' @rdname print
print.selection<-function(x,...){
  catln <- function (...)  base::cat(..., "\n", sep = "")
  print <- function (..., skip = 0, indent = 0) {
    output <- capture.output(base::print(...))
    if (skip > 0) output <- output[-seq_len(skip)]
    indent <- paste0(rep(" ", indent), collapse = "")
    cat(paste0(indent, output, "\n"), sep = "")
  }
  
  ## call
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Subset:")
  if( !is.null(colnames(x$X))){ catln("  Feature Name: ", paste(colnames(x$X)[x$ID_selected], sep = " ", collapse = ","))}
  catln("  Feature Index: ", paste(x$ID_selected, sep = "\n", collapse = ","))
  
  ## done
  invisible(x)
  
}

#' 
#' @export
#' @method print summary.smle
#' @rdname print

print.summary.smle <- function(x, ...){
  
  catln <- function (...)  base::cat(..., "\n", sep = "")
  print <- function (..., skip = 0, indent = 0){
    output <- capture.output(base::print(...))
    if (skip > 0) output <- output[-seq_len(skip)]
    indent <- paste0(rep(" ", indent), collapse = "")
    cat(paste0(indent, output, "\n"), sep = "")
  }
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Summary:")
  cat("\n")
  catln("  Dim of Y: " , paste(x$DimY,collapse = ' x '))
  catln("  Dim of X: " , paste(x$DimX,collapse = ' x '))
  catln("  Model type: " , x$family)
  catln("  Retained model size : ", x$size)
  
  if( !is.null(colnames(x$X))){ catln("  Feature Name: ", paste(colnames(x$X)[x$ID_retained], sep = " ", collapse = ","))}
  catln("  Feature Index: ", paste(x$ID_retained, sep = "\n", collapse = ","))
  catln("  Coefficients estimated by IHT : ", paste(format(x$coef_retained, digits = 3),collapse = ' '))
  catln("  Number of IHT iteration steps : ",as.character(x$steps))

  
  invisible(x)
}

#' @export
#' @method print summary.selection
#' @rdname print
print.summary.selection <- function(x, ...){
  
  catln <- function (...)  base::cat(..., "\n", sep = "")
  print <- function (..., skip = 0, indent = 0){
    output <- capture.output(base::print(...))
    if (skip > 0) output <- output[-seq_len(skip)]
    indent <- paste0(rep(" ", indent), collapse = "")
    cat(paste0(indent, output, "\n"), sep = "")
  }
  catln("Call:")
  catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Summary:")
  catln(" ")
  catln("  Dim of Y: " , paste(x$DimY,collapse = ' x '))
  catln("  Dim of X: " , paste(x$DimX,collapse = ' x '))
  catln("  Model type: " , x$family)
  catln("  K selected : ", x$size)
  if( !is.null(colnames(x$X))){ catln("  Feature Name: ", paste(colnames(x$X)[x$ID_selected], sep = " ", collapse = ","))}
  catln("  K features index: ", paste(x$ID_selected, sep = "\n", collapse = ","))
  catln("  Selection criterion : ", x$criterion)
  if(x$criterion=='ebic'){  catln("  Gamma for ebic : ", as.character(x$gamma_ebic))}
  if(x$vote=='T'){  catln("  Features selected by voting : ", x$ID_voted)}
  
  invisible(x)
}



#' @export
#' @method print sdata
#' @rdname print
print.sdata<-function(x,...){
  catln <- function (...)  base::cat(..., "\n", sep = "")
  print <- function (..., skip = 0, indent = 0) {
    output <- capture.output(base::print(...))
    if (skip > 0) output <- output[-seq_len(skip)]
    indent <- paste0(rep(" ", indent), collapse = "")
    cat(paste0(indent, output, "\n"), sep = "")
  }
  
  ## call
  catln("Call:")
  catln(" ", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catln(" ")
  catln(paste("An object of class",class(x)))
  catln(" ")
  catln("Simulated Dataset Properties:")
  catln(" Dim of Y: " , paste(c(length(x$Y),1),collapse = ' x '))
  catln(" Dim of X: " , paste(dim(x$X),collapse = ' x '))
  catln(" Correlation: " , x$correlation)
  if(x$correlation != "independent"){catln(" Rho: ", x$rho)}
  catln(" Index of Causal Features: " , paste(x$subset_true,collapse = ','))
  if(x$ctg ==T){catln("Design matrix concludes categorical features" )}
  catln(" Model Type: ", x$family)

  ## done
  invisible(x)
  
}
