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
#' First, we generate a \eqn{p} by \eqn{1} model coefficient vector beta with all entries being zero, except for the positions specified in \code{pos_truecoef},
#' on which \code{effect_truecoef} is used. When \code{pos_truecoef} is not specified, we randomly choose \code{num_truecoef} positions from the coefficient
#' vector. When \code{effect_truecoef} is not specified, we randomly set the strength of the true model coefficients as follow:
#' \deqn{(0.5+U) Z,}
#' where \eqn{U} is sampled from a uniform distribution from 0 to 1,  and \eqn{Z} is sampled from a binomial distribution \eqn{P(Z=1)=1/2,P(Z=-1)=1/2}.
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
#' Compound symmetry (CS): candidate features \eqn{x_1,..., x_p} are joint normal, marginally \eqn{N( 0, 1)}, with \eqn{cov(x_j, x_h) =\rho/2} if \eqn{j}, \eqn{h}
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
#' @param effect_truecoef  Effect size corresponding to the features in \code{pos_truecoef}. If not specified, effect size is sampled based on a uniform distribution and direction is randomly sampled.  See Details.
#'
#' @param pos_ctgidx Vector of indices denoting which columns are categorical.
#'
#' @param pos_truecoef Vector of indices denoting which features (columns) affect the response variable. If not specified, positions are randomly sampled. See Details for more information.
#'
#' @param family Model type for the response variable.
#' \code{"gaussian"} for normally distributed data, \code{poisson} for non-negative counts,
#' \code{"binomial"} for binary (0-1).
#'
#' @param correlation Correlation structure among features. \code{correlation = "ID"} for independent,
#' \code{correlation = 'MA'} for moving average, \code{correlation = "CS"} for compound symmetry, \code{correlation = "AR"} 
#' for auto regressive Default is \code{"ID"}. For more information see Details.
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
#' \item{X}{Feature matrix or dataframe (matrix if \code{num_ctgidx =FALSE} and dataframe otherwise).}
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



  ##-----------------------------------------------------Argument Check and Initialize Parameters--------------------------------------------------------
  correlation=match.arg(correlation)
  if(is.null(correlation)){correlation <- 'ID'}
  family=match.arg(family)
  cl<-match.call()
  cl[[1]] <- as.name("Gen_Data")
  
  #Argument check for coefficients
  if (is.null(num_truecoef) ) {

    if (is.null(pos_truecoef) ) {

      ## No value for num_,pos_ truecoef, num_truecoef =5 for default and pos is randomly sampled from 1:p

      num_truecoef<-5; pos_truecoef <- sample(1:p , size = num_truecoef , replace = F)

      }else{

      num_truecoef = length(pos_truecoef)

      }
  }else{

    if (!(class(num_truecoef) %in% c("numeric",'integer'))| length(num_truecoef) != 1 | num_truecoef%%1 != 0 | num_truecoef < 0){

      stop("The number of features affect the response variable should be an postive integer.")

      }
    if (is.null(pos_truecoef)) {

      # Randomly sample from 1:p.

      pos_truecoef <- sample(1:p , size = num_truecoef , replace = F)

      }else{

        if (!(class(pos_truecoef) %in% c("numeric",'integer'))  | any(pos_truecoef%%1 != 0) |  any(pos_truecoef <= 0) ){

          stop("Position for features affects the response variable should be a list of postive integer.")
        }
        if ( length(pos_truecoef) != num_truecoef ) {

          stop("Coeffiention length does not match")

        }
        if ( length( pos_truecoef ) != length( unique( pos_truecoef ) ) ) {

          pos_truecoef <- unique( pos_truecoef )

          warning("Position for features affects the response variable should be unique.")

        }
      }
  }


  #effects

  if ( is.null(effect_truecoef )){

      #effect_truecoef <- runif( num_truecoef , min = -1 , max = 1 )

      effect_truecoef <- (0.5+abs(rnorm(num_truecoef))) * sample(c(1,-1),num_truecoef,replace = T)

  }else{

    if (length(effect_truecoef) != num_truecoef){

      stop("Effects should match number of non-zero features.")

    }

  }

  ##------------------------------------------Create Data-----------------------------------------------

  #Numerical data

  D<-Numeric_Gen(n, p, pos_truecoef, effect_truecoef, correlation, rho, family, sigma,cl )

  # Numerical data and Categorical parameter check


  if ( is.null(num_ctgidx) ) {

    if ( is.null(pos_ctgidx) ){

      # Numerical data as default

      num_ctgidx<- 0

      return(  D  )

      }else{

        num_ctgidx <- length( pos_ctgidx )

      }

  }else{
    if ( num_ctgidx== 0 ) { return( D ) }

    if (!(class(num_ctgidx) %in% c("numeric",'integer'))| length(num_ctgidx) != 1 | num_ctgidx %%1 != 0 | num_ctgidx < 0){

      stop("The number of categorical features should be an postive integer.")
    }

    if( is.null(pos_ctgidx)  ){

    pos_ctgidx <- sample( 1:p, num_ctgidx, replace = FALSE )

    }else{

      if ( !(class(pos_ctgidx) %in% c("numeric",'integer')) | any(pos_ctgidx%%1 != 0) |  any(pos_ctgidx <= 0) ){

        stop("Position for categorical features should be a list of postive integer.")
      }
      if ( length(pos_ctgidx) != num_ctgidx ) {

        stop("Position length does not match the number of categorical features")

      }
      if ( length( pos_ctgidx ) != length( unique( pos_ctgidx ) ) ) {

        pos_ctgidx <- unique( pos_ctgidx )

        warning("Position for categorcial features should be unique.")

      }
      }

    }


    #level of categorical data

    if ( is.null(level_ctgidx )){

      level_ctgidx <- 2

    }


    if (length( level_ctgidx ) == 1 & all(level_ctgidx >= 0) ){

      # Same level for all categorical data


      Z<-lapply(seq(1:p),function(i){if (i %in% pos_ctgidx){

        cut(D$X[,i],breaks = level_ctgidx,labels = LETTERS[1:level_ctgidx])

      }else{D$X[,i]}


        }
        )

    }else if( length( level_ctgidx ) == length( pos_ctgidx )){


      Z<-lapply(seq(1:p),function(i){if (i %in% pos_ctgidx){

        j<-match(i,pos_ctgidx); cut(D$X[,i],breaks = level_ctgidx[j],labels = LETTERS[1:level_ctgidx[j]])

        }else{D$X[,i]}
        }
        )

    }else{

       stop("level of categorical features does not match position of categorical feature" )

    }

    cn<-unlist(lapply(seq(1:p),function(i){if (i %in% pos_ctgidx){paste0("C",i)}else(paste0("N",i))}))

    Z<-as.data.frame(Z,stringsAsFactors = FALSE)

    colnames(Z)<-cn
    
    correlation=switch(correlation,
                       
                       'ID'='independent',
                       
                       'MA'='moving average',
                       
                       "AR"="auto regressive",
                       
                       'CS'='compound symmetry')

    D<-list(call = cl, Y = D$Y,X = Z , subset_true = sort(pos_truecoef), coef_true= effect_truecoef,

            family = family, categorical = TRUE, correlation = correlation, rho = rho, CI= sort(pos_ctgidx))

    class(D) <- "sdata"

    return(D)

  }


  Gen_DesignMatrix<-function(correalation=c('ID','AR','MA','CS'),
                             N,P,pos_truecoef,rho){

    correalation<-match.arg(correalation)

    if (N < 0 || P < 0 || N%%1!=0 || P%%1!=0){

      stop("N,P should be postive integers")

    }

    if(correalation=='ID'){

      return(matrix(rnorm(N*P), nrow=N, ncol=P))

    }else if(correalation=='MA'){

      Z<-matrix(rnorm(N*(P+2)), nrow=N, ncol=P+2)

      return(sapply(seq(1:P),function(i){

        sqrt(1-rho)*Z[,i]+sqrt(rho/2)*Z[,i+1]+sqrt(rho/2)*Z[,i+2]

      }))

    }else if( correalation=='AR'){

      ar1_cor <- function(P, rho) {

        exponent <- abs(matrix(1:P - 1, nrow = P, ncol = P, byrow = TRUE) - (1:P - 1))

        rho^exponent

      }

      return(mvnfast::rmvn(n = N, mu = rep(0,P), ar1_cor(P,rho)))

    }else{

      CC<-matrix(rho, ncol=P, nrow=P)

      CC[pos_truecoef, pos_truecoef]<- rho/2

      diag(x=CC)<-1
      if(!is.positive.definite(CC)){ stop("The Given Correlation matrix is not positive difinite") }
      return(mvnfast::rmvn(n = N, mu = rep(0,P), CC))

    }

  }


  Numeric_Gen<-function(n , p , pos_truecoef , effect_truecoef, correlation, rho, family, sigma, cl ){

    numeric_data<-Gen_DesignMatrix(correlation,n,p,pos_truecoef=pos_truecoef,rho)

    BETA<-matrix(0, nrow=p, ncol=1)

    BETA[pos_truecoef, 1] <- effect_truecoef

    theta<-tcrossprod(as.matrix(numeric_data),t(BETA))

    if(family=="gaussian"){

      Y<- theta + rnorm(n, mean = 0, sd = sigma)

    }else if(family=="binomial"){

      pi <- exp(theta) / (1 + exp(theta))

      Y  <- rbinom(n, size = 1, prob = pi)

    }else{

      mmu <- exp(theta)

      Y  <- rpois(n, lambda = mmu)

    }
    correlation=switch(correlation,

                       'ID'='independent',

                       'MA'='moving average',

                       "AR"="auto regressive",

                       'CS'='compound symmetry')

    D<-list(call = cl, Y = Y,X = numeric_data , subset_true = sort(pos_truecoef), coef_true= effect_truecoef[order(pos_truecoef)],
            family = family, categorical = FALSE, correlation = correlation,rho=rho,CI = NULL)

    class(D)<-"sdata"

    return(D)

  }










