#' @title
#' Data simulator for high-dimensional GLMs
#'
#' @description
#' This function generates synthetic datasets from GLMs with a user-specified correlation structure.
#' It permits both numerical and categorical features, whose quantity can be larger than the sample size.
#'
#'
#' @import mvnfast mvnfast
#' @import matrixcalc
#' @details
#'
#' Simulated data \eqn{(y_i , x_i)} where \eqn{ x_i = (x_{i1},x_{i2} , . . . , x_{ip})} are generated as follows:
#' First, we generate a \eqn{p} by \eqn{1} model coefficient vector beta with all 
#' entries being zero, except for the positions specified in \code{pos_truecoef},
#' on which \code{effect_truecoef} is used. When \code{pos_truecoef} is not specified, 
#' we randomly choose \code{num_truecoef} positions from the coefficient
#' vector. When \code{effect_truecoef} is not specified, we randomly set the strength 
#' of the true model coefficients as follow:
#' \deqn{(0.5+U) Z,}
#' where \eqn{U} is sampled from a uniform distribution from 0 to 1,  
#' and \eqn{Z} is sampled from a binomial distribution \eqn{P(Z=1)=1/2,P(Z=-1)=1/2}.
#'
#' Next, we generate a \eqn{n} by \eqn{p} feature matrix \eqn{X} according to the model selected with
#' \code{correlation} and specified as follows.
#'
#' Independent (ID):  all features are independently generated from \eqn{N( 0, 1)}.
#'
#' Moving average (MA): candidate features \eqn{x_1,..., x_p} are joint normal,
#' marginally \eqn{N( 0, 1)}, with
#'
#' \eqn{cov(x_j, x_{j-1}) = \rho}, \eqn{cov(x_j, x_{j-2}) = \rho/2} and \eqn{cov(x_j, x_h) = 0} for \eqn{|j-h|>3}.
#'
#' Compound symmetry (CS): candidate features \eqn{x_1,..., x_p} are joint normal,
#'  marginally \eqn{N( 0, 1)}, with \eqn{cov(x_j, x_h) =\rho/2} if \eqn{j}, \eqn{h}
#' are both in the set of important features and \eqn{cov(x_j, x_h)=\rho} when only
#' one of \eqn{j} or \eqn{h} are in the set of important features.
#' 
#' Auto-regressive (AR): candidate features \eqn{x_1,..., x_p} are joint normal, marginally \eqn{N( 0, 1)}, with
#'
#' \eqn{cov(x_j, x_h) = \rho^{|j-h|}} for all \eqn{j} and \eqn{h}. The correlation strength \eqn{\rho} is controlled by the argument \code{rho}. 
#' 
#' 
#' Then, we generate the response variable \eqn{Y} according to its response type, which is controlled by the argument \code{family}
#' For the Gaussian model, \eqn{y_i =x_i\beta + \epsilon_i} where \eqn{\epsilon_i} is \eqn{N( 0, 1)} for \eqn{i} from \eqn{1} to \eqn{n}. 
#' For the binary model let \eqn{\pi_i = P(Y = 1|x_i)}. We sample \eqn{y_i} from Bernoulli(\eqn{\pi_i}) where logit\eqn{(\pi_i) = x_i \beta}.
#' Finally, for the Poisson model, \eqn{y_i} is generated from the Poisson distribution with the link \eqn{\pi_i} = exp\eqn{(x_i\beta )}.
#' For more details see the reference below.
#'
#' @param n Sample size, number of rows for the feature matrix to be generated.
#'
#' @param p Number of columns for the feature matrix to be generated.
#'
#' @param num_ctgidx The number of features that are categorical. Set to \code{FALSE} for only numerical features. Default is \code{FALSE}.
#'
#' @param num_truecoef The number of features (columns) that affect response. Default is 5.
#'
#' @param level_ctgidx  Vector to indicate the number of levels for the categorical features in \code{pos_ctgidx}. Default is 2.
#'
#' @param effect_truecoef  Effect size corresponding to the features in \code{pos_truecoef}. 
#' If not specified, effect size is sampled based on a uniform distribution and direction is randomly sampled. See Details.
#'
#' @param pos_ctgidx Vector of indices denoting which columns are categorical.
#'
#' @param pos_truecoef Vector of indices denoting which features (columns) affect the response variable. 
#' If not specified, positions are randomly sampled. See Details for more information.
#'
#' @param family Model type for the response variable.
#' \code{"gaussian"} for normally distributed data, \code{poisson} for non-negative counts,
#' \code{"binomial"} for binary (0-1).
#'
#' @param correlation Correlation structure among features. \code{correlation = "ID"} for independent,
#' \code{correlation = 'MA'} for moving average, \code{correlation = "CS"} for compound symmetry, \code{correlation = "AR"} 
#' for auto regressive. Default is \code{"ID"}. For more information see Details.
#'
#' @param rho Parameter controlling the correlation strength, default is \code{0.2}. See Details.
#'
#' @param sigma Parameter for noise level.
#'
#'
#'
#'
#' @references
#' Xu, C. and Chen, J. (2014). The Sparse MLE for Ultrahigh-Dimensional Feature
#' Screening, \emph{Journal of the American Statistical Association}, \bold{109}(507), 1257-1269
#'
#'
#'
#'
#' @return
#' \item{call}{The call that produced this object.}
#' \item{Y}{Response variable vector of length \eqn{n}.}
#'
#' \item{X}{Feature matrix or data.frame (matrix if \code{num_ctgidx =FALSE} and data.frame otherwise).}
#'
#' \item{subset_true}{Vector of column indices of X for the features that affect the response variables (relevant features).}
#'
#' \item{coef_true}{Vector of effects for the features that affect the response variables.}
#' 
#' \item{categorical}{Logical flag whether the model contains categorical features.}
#' 
#' \item{CI}{Indices of categorical features when \code{categorical = TRUE}.}
#'
#' rho,family,correlation are return of arguments passed in the function call.
#' 
#'
#' @export
#'
#' @examples
#' \donttest{
#' #Simulating data with binomial response and auto-regressive structure.
#' set.seed(1)
#' Data <- Gen_Data(n = 500, p = 2000, family = "binomial", correlation = "AR")
#' cor(Data$X[,1:5])
#' print(Data)
#' }
#'
Gen_Data<-function(n = 200,p = 1000,sigma = 1,
                   num_ctgidx = NULL,  pos_ctgidx = NULL,
                   num_truecoef = NULL,  pos_truecoef = NULL,
                   level_ctgidx = NULL,  effect_truecoef= NULL,
                   correlation = c("ID","AR","MA","CS"),
                   rho = 0.2, family = c("gaussian","binomial","poisson")){
  
  
  
  ##--------------Argument Check and Initialize Parameters----------
  
  correlation=match.arg(correlation)
  
  if(is.null(correlation)){correlation <- 'ID'}
  
  family=match.arg(family)
  
  cl<-match.call()
  
  cl[[1]] <- as.name("Gen_Data")
  
  #Argument check for coefficients
  
  if (n < 0 || p < 0 || n%%1!=0 || p%%1!=0){
    
    stop("Both n and p must be positive integers")
    
  }
  
  
  
  if (is.null(num_truecoef) ) {
    
    if (is.null(pos_truecoef) ) {
      
      ## If no input value for number and position of coefficient, 5 is assigned as numbers and position is randomly sampled from 1:p.
      
      num_truecoef <- 5
      
      pos_truecoef <- sample(1:p , size = num_truecoef , replace = FALSE)
      
    }else{
      
      num_truecoef = length(pos_truecoef)
      
    }
  }else{
    
    if (!(class(num_truecoef) %in% c('numeric','integer'))| length(num_truecoef) != 1 | num_truecoef%%1 != 0 | num_truecoef < 0){
      
      #  Catch the error of negative non-integer input.
      
      stop("The number of features affecting the response variable should be a positive integer.")
      
    }
    if (is.null(pos_truecoef)) {
      
      # Randomly sample from 1:p.
      
      pos_truecoef <- sample(1:p , size = num_truecoef , replace = F)
      
    }else{
      
      if (!(class(pos_truecoef) %in% c('numeric','integer'))  | any(pos_truecoef%%1 != 0) |  any(pos_truecoef <= 0) ){
        
        stop("Position of features affecting the response variable should be a list of positive integers.")
      }
      
      if ( length( pos_truecoef ) != length( unique( pos_truecoef ) ) ) {
        
        # Duplicated elements will be ignored. 
        
        pos_truecoef <- unique( pos_truecoef )
        
        warning("Position of features affecting the response variable should be unique.")
        
      }
      if ( length(pos_truecoef) != num_truecoef ) {
        
        stop("Value of 'num_truecoef' must match the length of 'pos_truecoef'.")
        
      }
    }
  }
  
  
  # Initial effects check
  
  if ( is.null(effect_truecoef )){
    
    # A good sample setting for SMLE screening.
    
    effect_truecoef <- ( 0.5 + abs( rnorm( num_truecoef ))) * sample( c(1,-1),
                                                    num_truecoef , replace = TRUE)
    
  }else{
    
    if (length(effect_truecoef) != num_truecoef){
      
      stop("Length of 'effect_truecoef' does not match the specified number of causal features ('num_truecoef').")
      
    }
    
  }
  
  # Numerical data and Categorical parameter check
  
  output_ctg = TRUE
  
  if ( is.null(num_ctgidx) ) {
    
    if ( is.null(pos_ctgidx) ){
      
      # Output numerical data if nothing specified. 
      
      num_ctgidx <- 0
      
      output_ctg <- FALSE
      
    }else{
      
      num_ctgidx <- length( pos_ctgidx )
      
    }
    
  }else{
    
    if ( !num_ctgidx ) {
      
      output_ctg <- FALSE
      
    }
    
    if (!(class(num_ctgidx) %in% c('numeric','integer'))| length(num_ctgidx) != 1 | num_ctgidx %%1 != 0 | num_ctgidx < 0){
      
      stop("The number of categorical features must be a positive integer.")
    }
    
    if( is.null(pos_ctgidx)  ){
      
      pos_ctgidx <- sample( 1:p, num_ctgidx, replace = FALSE )
      
    }else{
      
      if ( !(class(pos_ctgidx) %in% c('numeric','integer')) | any(pos_ctgidx%%1 != 0) |  any(pos_ctgidx <= 0) ){
        
        stop("The vector 'pos_ctgidx' must be a list of positive integers.")
      }
      
      if ( length( pos_ctgidx ) != length( unique( pos_ctgidx ) ) ) {
        
        pos_ctgidx <- unique( pos_ctgidx )
        
        warning("Positions provided for the categorcial features should be unique.")
        
      }
      if ( length(pos_ctgidx) != num_ctgidx ) {
        
        stop("The length of 'pos_ctgidx' does not match the specified number of categorical features (num_ctgidx).")
        
      }
    }
  }
  
  
  #Default level of categorical data
  
  if ( is.null(level_ctgidx )){
    
    level_ctgidx <- 2
  }
  
  ##------------------------------------------Create Data-----------------------------------------------
  
  #Numerical data
  
  D<-Numeric_Gen(n, p, pos_truecoef, effect_truecoef, correlation, rho, family, sigma,cl )
  
  correlation = switch(correlation,
                     
                     'ID' = 'independent',
                     
                     'MA' = 'moving average',
                     
                     'AR' = 'auto regressive',
                     
                     'CS' = 'compound symmetry')
  
  if( output_ctg ){
    
    if (length( level_ctgidx ) == 1 ){
      
      # Same level for all categorical data
      
      level_ctgidx <- rep(level_ctgidx, times = length( pos_ctgidx ) )
      
    }
    
    if(length( level_ctgidx ) != length( pos_ctgidx )){
      
      stop("The length of the vector 'level_ctgidx' should be the same as the length of the vector 'pos_ctgidx'.")
    
      }
      
    Z <- data.frame(D$X)
      
    for(i in 1:length(pos_ctgidx)){
        
        Z[,pos_ctgidx[i]] <- cut(D$X[,pos_ctgidx[i]],breaks = level_ctgidx[i],labels = LETTERS[1:level_ctgidx[i]])
        
      }
      
    colnames(Z)[pos_ctgidx] <- paste0('C',pos_ctgidx)

    # Set name
    
    Data_out <- list(call = cl, Y = D$Y , X = Z , subset_true = sort(pos_truecoef),
                     
                     coef_true = effect_truecoef , family = family, categorical = TRUE , 
                     
                     correlation = correlation , rho = rho , CI = sort(pos_ctgidx))
    
  }else{ 
    
    Data_out <- D 
    
    Data_out$correlation = correlation
  }
  
  class(Data_out) <- 'sdata'
  
  return(Data_out)
  
}


Gen_DesignMatrix<-function(correlation=c('ID','AR','MA','CS'),
                           N,P,pos_truecoef,rho){
  
  correlation<-match.arg(correlation)
  
  #In all cases, a matrix consisting of the numeric features is returned.
  if(correlation=='ID'){
    
    return(matrix(rnorm(N*P), nrow=N, ncol=P))
    
  }else if(correlation=='MA'){
    
    Z<-matrix(rnorm(N*(P+2)), nrow=N, ncol=P+2)
    
    return(sapply(seq(1:P),function(i){
      
      sqrt(1-rho)*Z[,i]+sqrt(rho/2)*Z[,i+1]+sqrt(rho/2)*Z[,i+2]
      
    }))
    
  }else if( correlation=='AR'){
    
    ar1_cor <- function(P, rho) {
      
      exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))
      
      rho^exponent
      
    }
    
    return(mvnfast::rmvn(n = N, mu = rep(0,P), ar1_cor(P,rho)))
    
  }else{
    
    CC<-matrix(rho, ncol=P, nrow=P)
    
    CC[pos_truecoef, pos_truecoef]<- rho/2
    
    diag(x=CC)<-1
    
    if(!is.positive.definite(CC)){ stop("The given correlation matrix is not 
        positive difinite. Please try a different rho to get a valid matrix.") }
    
    return(mvnfast::rmvn(n = N, mu = rep(0,P), CC))
    
  }
  
}


Numeric_Gen<-function(n , p , pos_truecoef , effect_truecoef, correlation, 
                      rho, family, sigma, cl ){
  
  numeric_data<-Gen_DesignMatrix(correlation,n,p,pos_truecoef=pos_truecoef,rho)
  
  Beta <- rep(0,p)
  
  Beta[pos_truecoef] <- effect_truecoef
  
  theta <- numeric_data %*% Beta
  
  if(family=='gaussian'){
    
    Y<-  rnorm(n, mean = theta, sd = sigma)
    
  }else if(family=='binomial'){
    
    pi <- exp(theta) / (1 + exp(theta))
    
    Y  <- rbinom(n, size = 1, prob = pi)
    
  }else{
    
    mmu <- exp(theta)
    
    Y  <- rpois(n, lambda = mmu)
    
  }
  
  return(list(call = cl, Y = Y,X = numeric_data , subset_true = sort(pos_truecoef), 
              coef_true= effect_truecoef[order(pos_truecoef)],
              family = family, categorical = FALSE, correlation = correlation,
              rho=rho))
  
}
