\name{smle-package}
\alias{smle-package}
\docType{package}
\title{
Joint SMLE-screening for generalized linear models
}
\description{
Feature screening is a powerful tool in processing ultrahigh dimensional data. It attempts to screen
out most irrelevant features in preparation for a more elaborate analysis. This package provides an efficient implementation of SMLE-screening for linear, logistic, and Poisson models, where the joint effects among features are naturally incorporated in the screening process. The package also provides a function for conducting accurate post-screening feature selection based on an iterative hard-thresholding procedure and a user-specified selection criterion.
}
\details{
\tabular{ll}{
Package: \tab smle\cr
Type: \tab Package\cr
Version: \tab 2.1-0\cr
Date: \tab 2023-01-04\cr
License: \tab GPL-3\cr
}
Input a \eqn{n \times 1} response vector Y and a \eqn{n \times p} predictor (feature) matrix X. The package outputs a set of \eqn{k < n} features that seem to be most relevant for joint regression. Moreover, the package provides a data simulator that generates synthetic datasets from high-dimensional GLMs, which accommodate both numerical and categorical features with commonly used correlation structures. \cr

Key functions:\cr
\code{Gen_Data}\cr
\code{SMLE}\cr
\code{smle_select}\cr
}
\author{Qianxiang Zang, Chen Xu, Kelly Burkett\cr
Maintainer: Qianxiang Zang  <qzang023@uottawa.ca>
}
\references{Xu, C. and Chen, J. (2014)
  The Sparse MLE for Ultrahigh-Dimensional Feature Screening
  \emph{Journal of the American Statistical Association}, \bold{109}(507), 1257–1269.\cr

  Friedman, J., Hastie, T. and Tibshirani, R. (2010)
  Regularization Paths for Generalized Linear Models via Coordinate
    Descent
  \emph{Journal of Statistical Software}, \bold{33}(1), 1-22.\cr
    }

\examples{
set.seed(1)
#Generate correlated data
Data <- Gen_Data(n = 200, p = 5000, correlation = "MA",family = "gaussian")
print(Data)

# joint feature screening via SMLE
fit <- SMLE(Y = Data$Y, X = Data$X, k = 10, family = "gaussian")
print(fit)
summary(fit)
plot(fit)

#Are there any features missed after screening?
setdiff(Data$subset_true, fit$ID_retained)

# Elaborative selection after screening
fit_s <- smle_select(fit, gamma_ebic = 0.5, vote = FALSE)

#Are there any features missed after selection? 
setdiff(Data$subset_true, fit_s$ID_selected)
print(fit_s)
summary(fit_s)
plot(fit_s)
}
