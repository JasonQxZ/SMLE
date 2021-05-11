## Copyright  2009-2020  Marc Hofmann and Achim Zeileis
##
## This file is part of the 'lmSubsets' R extension.
##
## 'lmSubsets' is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## 'lmSubsets' is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with 'lmSubsets'.  If not, see <http://www.gnu.org/licenses/>.



#################
##  LMSUBSETS  ##
#################


## summary for 'lmSubsets' objects
##
## Arguments:
##   object - (lmSubsets)
##   ...    - ignored
##   na.rm  - (logical)
##
## Result: (summary.lmSubsets)
##
summary.lmSubsets <- function (object, ..., na.rm = TRUE) {
    ans <- object[c("call", "terms", "nvar", "nbest", "size",
                    if (!is.null(object$weights)) "weights")]

    ans$stats <- local({
        s <- stats_lmSubsets(object)

        r2 <- r2_stats(s)
        r2adj <- r2adj_stats(s, r2 = r2)

        data.frame(SIZE = object$submodel$SIZE, BEST = object$submodel$BEST,
                   sigma = sigma_stats(s), R2 = r2, R2adj = r2adj,
                   pval = pval_stats(s), Cp = cp_stats(s), AIC = aic_stats(s),
                   BIC = bic_stats(s))
    })

    if (na.rm) {
        ans$stats <- na.omit(ans$stats)
    }

    class(ans) <- "summary.lmSubsets"

    ans
}


## print 'lmSubsets' summary
##
## Arguments:
##   x   - (summary.lmSubsets)
##   ... - forwarded
##
## Result: (summary.lmSubsets) invisible
##
print.summary.lmSubsets <- function (x, ...) {
    catln <- function (...)  base::cat(..., "\n", sep = "")
    print <- function (..., skip = 0, indent = 0) {
        output <- capture.output(base::print(...))
        if (skip > 0)  output <- output[-seq_len(skip)]
        indent <- paste0(rep(" ", indent), collapse = "")
        cat(paste0(indent, output, "\n"), sep = "")
    }

    ## call
    catln("Call:")
    catln("  ", paste(deparse(x$call), sep = "\n", collapse = "\n"))

    catln()

    ## stats
    stats <- within(x$stats, {
        SIZE <- ifelse(SIZE != head(c(-1, SIZE), -1), SIZE, NA)
    })
    stats <- sapply(names(stats), function (nm) {
        x <- stats[[nm]]

        if (nm == "pval") format_pval(x, na.encode = FALSE, ...)
        else format_default(x, na.encode = FALSE, ...)
    })
    stats <- ifelse(is.na(stats), "", stats)
    rownames(stats) <- rep("", nrow(stats))

    catln("Statistics:")
    print(stats, quote = FALSE, indent = 2)

    ## done
    invisible(x)
}



################
##  LMSELECT  ##
################


## summary for 'lmSelect' objects
##
## Arguments:
##   object - (lmSelect)
##   ...    - ignored
##   na.rm  - (logical)
##
## Result: (summary.lmSelect)
##
summary.lmSelect <- function (object, ..., na.rm = TRUE) {
    ans <- object[c("call", "terms", "nvar", "nbest", "size",
                    if (!is.null(object$weights)) "weights")]

    ans$stats <- local({
        s <- stats_lmSelect(object)

        r2 <- r2_stats(s)
        r2adj <- r2adj_stats(s, r2 = r2)

        data.frame(SIZE = object$submodel$SIZE, BEST = object$submodel$BEST,
                   sigma = sigma_stats(s), R2 = r2, R2adj = r2adj,
                   pval = pval_stats(s), Cp = cp_stats(s), AIC = aic_stats(s),
                   BIC = bic_stats(s))
    })

    if (na.rm) {
        ans$stats <- na.omit(ans$stats)
    }

    class(ans) <- "summary.lmSelect"

    ans
}


## print 'lmSelect' summary
##
## Arguments:
##   x   - (summary.lmSelect)
##   ... - forwarded
##
## Result: (summary.lmSelect) invisible
##
print.summary.lmSelect <- function (x, ...) {
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

    catln()

    ## stats
    stats <- sapply(names(x$stats), function (nm) {
        x <- x$stats[[nm]]

        if (nm == "pval") format_pval(x, na.encode = FALSE, ...)
        else format_default(x, na.encode = FALSE, ...)
    })
    stats <- ifelse(is.na(stats), "", stats)
    rownames(stats) <- rep("", nrow(stats))

    catln("Statistics:")
    print(stats, quote = FALSE, indent = 2)

    ## done
    invisible(x)
}
