#' Print a SMLE object
#' 
#' @rdname print
#' @importFrom utils capture.output
#' @description This functions prints the retained subset of a SMLE object.
#'
#' @param x Fitted object.
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
  if( !is.null(colnames(x$X))){ catln("  Feature Name: ", paste(colnames(x$I$CM)[x$ID_Retained], sep = " ", collapse = ","))}
  catln("  Feature Index: ", paste(x$ID_Retained, sep = "\n", collapse = ","))
  

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
  
  catln("Subset:")
  
  catln("  Feature Name: ", paste(colnames(x$X)[x$ID_Selected], sep = " ", collapse = ","))
  catln("  Feature Index: ", paste(x$ID_Selected, sep = "\n", collapse = " "))
  
  
  ## done
  invisible(x)
  
}

#' 
#' @export
#' @method print summary.smle
#' @rdname print

print.summary.smle <- function(x, ...){
  
  if(!exists("digits")){digits = max(3L, getOption("digits") - 3L)}
  if(!exists("symbolic.cor")){symbolic.cor = x$symbolic.cor}
  if(!exists("signif.stars")){signif.stars = getOption("show.signif.stars")}
  
  
  
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Deviance Residuals: \n")
  if(x$df.residual > 5) {
    x$deviance.resid <- setNames(quantile(x$deviance.resid, na.rm = TRUE),
                                 c("Min", "1Q", "Median", "3Q", "Max"))
  }
  xx <- zapsmall(x$deviance.resid, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  
  if(length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  } else {
    ## df component added in 1.8.0
    ## partial matching problem here.
    df <- if ("df" %in% names(x)) x[["df"]] else NULL
    if (!is.null(df) && (nsingular <- df[3L] - df[1L]))
      cat("\nCoefficients: (", nsingular,
          " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if(!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L,
                      dimnames=list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
  }
  ##
  cat("\n(Dispersion parameter for ", x$family$family,
      " family taken to be ", format(x$dispersion), ")\n\n",
      apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                        "deviance:"),
                  format(unlist(x[c("null.deviance","deviance")]),
                         digits = max(5L, digits + 1L)), " on",
                  format(unlist(x[c("df.null","df.residual")])),
                  " degrees of freedom\n"),
            1L, paste, collapse = " "), sep = "")
  if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n\n",
      "Number of Fisher Scoring iterations: ", x$iter,
      "\n", sep = "")
  
  correl <- x$correlation
  if(!is.null(correl)) {
    # looks most sensible not to give NAs for undefined coefficients
    #         if(!is.null(aliased) && any(aliased)) {
    #             nc <- length(aliased)
    #             correl <- matrix(NA, nc, nc, dimnames = list(cn, cn))
    #             correl[!aliased, !aliased] <- x$correl
    #         }
    p <- NCOL(correl)
    if(p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if(is.logical(symbolic.cor) && symbolic.cor) {# NULL < 1.7.0 objects
        print(symnum(correl, abbr.colnames = NULL))
      } else {
        correl <- format(round(correl, 2L), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop=FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
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
  
  catln("Subset:")
  catln(" Dim_of_Y: " , paste(c(length(x$Y),1),collapse = ' x '))
  catln(" Dim_of_X: " , paste(dim(x$X),collapse = ' x '))
  catln(" Correlation: " , x$correlation)
  catln(" Index_of_Casual_Featrues: " , paste(x$subset_true,collapse = ','))
  if(x$Ctg ==T){catln("Design matrix concludes categorical features" )}
  catln(" Model_Type: ", x$family)
  ## done
  invisible(x)
  
}
