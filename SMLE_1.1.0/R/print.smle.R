#' Print a SMLE object
#'
#' @description This functions prints the retained subset of a SMLE object.
#'
#' @param x Fitted '\code{smle}' object.
#'
#' @param ... This argument is not used and listed for method consistency.
#'
#' @return
#' No return value.
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
  
  catln("Subset:")
  if( !is.null(colnames(x$X))){ catln("  Feature Name： ", paste(colnames(x$I$CM)[x$ID_Retained], sep = " ", collapse = ","))}
  catln("  Feature Index： ", paste(x$ID_Retained, sep = "\n", collapse = ","))
  

  ## done
  invisible(x)
  
}