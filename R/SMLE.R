#' Joint feature screening via sparse maximum likelihood estimation for GLMs
#' 
#' 
#' Input a \eqn{n} by \eqn{1} response \eqn{Y} and a \eqn{n} by \eqn{p} feature matrix \eqn{X};
#' the function uses SMLE to retain only a set of \eqn{k<n} features that seem
#' to be most related to the response variable. It thus serves as a pre-processing step for an
#' elaborative analysis. In SMLE, the joint effects between features are naturally
#' accounted for; this makes the screening more reliable. The function uses the
#' efficient iterative hard thresholding (IHT) algorithm with step parameter
#' adaptively tuned for fast convergence. Users can choose to further conduct
#' an elaborative selection after SMLE-screening. See \code{\link{smle_select}()} for more details.
#'
#' @importFrom glmnet glmnet
#' @importFrom grDevices dev.new dev.off rainbow
#' @importFrom graphics barplot legend lines par plot title
#' @importFrom stats dpois glm lm model.frame model.matrix quantile rbinom rnorm rpois sd .getXlevels model.response
#' @importFrom utils tail
#' 
#' @details
#' With the input \eqn{Y} and \eqn{X}, \code{\link{SMLE}()} conducts joint feature screening by running
#' iterative hard thresholding algorithm (IHT), where the default initial value is set to
#' be the Lasso estimate with the sparsity closest to the sample size minus one.
#' 
#' In \code{\link{SMLE}()}, the initial value for step size parameter \eqn{1/u} is 
#' determined as follows. When \code{coef_initial = 0}, we set \eqn{1/u = U / \sqrt{p}}.
#' When \code{coef_initial != 0}, we generate a sub-matrix \eqn{X_0} using the columns of \eqn{X} 
#' corresponding to the non-zero positions of \code{coef_initial} and set
#' \eqn{1/u = U/\sqrt{p}||X||^2_{\infty}} and recursively decrease the value of step size by 
#' \code{U_rate} to guarantee the likelihood increment. This strategy is called \eqn{u}-search.
#'
#' \code{\link{SMLE}()} terminates IHT iterations when either \code{tol} or \code{max_iter} is
#' satisfied. When \code{fast = TRUE}, the algorithm also stops when the non-zero
#' members of the coefficient estimates remain the same for 10 successive
#' iterations or the log-likelihood difference between coefficient estimates is less
#' than \eqn{0.01} times the log-likelihood increase of the first step, or 
#' \code{tol}\eqn{\sqrt k} is satisfied.
#' 
#' In \code{\link{SMLE}()}, categorical features are coded by dummy covariates with the
#' method specified in \code{codingtype}. Users can use \code{group} to specify
#' whether to treat those dummy covariates as a single group feature or as
#' individual features.
#' When \code{group = TRUE} with \code{penalize_mod = TRUE}, the effect for a group
#' of \eqn{J} dummy covariates is computed by
#' 
#' \deqn{ \beta_i = \sqrt{(\beta_1)^2+...+(\beta_J)^2}/\sqrt J,}
#' 
#' which will be treated as a single feature in IHT iterations. When \code{group = FALSE}, 
#' a group of \eqn{J} dummy covariates will be treated as \eqn{J} individual features in the IHT iterations; in this case, 
#' a categorical feature is retained after screening when at least one of the corresponding dummy covariates is retained.
#'
#' Since feature screening is usually a preprocessing step, users may wish to
#' further conduct an elaborative feature selection after screening. This can
#' be done by setting \code{selection = TRUE} in \code{\link{SMLE}()} or applying any existing
#' selection method on the output of \code{\link{SMLE}()}.
#'
#'
#' @param Y The response vector \eqn{Y} of dimension \eqn{n} by \eqn{1}. Quantitative for
#' \code{family = "gaussian"}, non-negative counts for \code{family = "poisson"},
#' binary (0-1) for \code{family = "binomial"}. Input \eqn{Y} should be \code{'numeric'}.
#'
#' @param X The \eqn{n} by \eqn{p} feature matrix \eqn{X} with each column denoting a feature
#' (covariate) and each row denoting an observation vector. The input should be
#' a \code{'matrix'} object for numerical data, and \code{'data.frame'} for categorical
#' data (or a mixture of numerical and categorical data). The algorithm will
#' treat covariates having class \code{'factor'} as categorical data and extend the data
#' frame dimension by the dummy columns needed for coding the categorical features.
#'
#' @param k Total number of features (including \code{keyset}) to be retained 
#' after screening. Default is the largest integer not 
#' exceeding \eqn{0.5}log\eqn{(n) n^{1/3}}.
#'
#' @param family Model assumption between \eqn{Y} and \eqn{X}; 
#' either a character string representing one of the built-in families, or else a glm() family object. 
#' The default model is Gaussian linear.
#' 
#' @param U A numerical multiplier of initial tuning step parameter in 
#' IHT algorithm. Default is 1. For binomial model, a larger initial value is recommended; a
#' smaller one is recommended for poisson model.
#' 
#' @param U_rate Decreasing rate in tuning step parameter \eqn{1/u} in IHT
#' algorithm. See Details.
#'
#' @param keyset A numeric vector with column indices for the key features that 
#' do not participate in feature screening and are forced to remain in the model.
#' The column indices for the key features should be from \code{data} if \code{'formula'} is used
#' or in \code{X} if \code{X} and \code{Y} are provided. The class of \code{keyset} 
#' can be \code{'numeric'},\code{'integer'} or \code{'character'}. Default is \code{NULL}.
#'
#' @param intercept A logical flag to indicate whether to an intercept be used in
#' the model. An intercept will not participate in screening.
#'
#' @param categorical A logical flag for whether the input feature matrix includes 
#' categorical features( either \code{'factor'} or \code{'character'}). \code{FALSE} 
#' treats all features as numerical and not check for whether there are categorical features; 
#' \code{TRUE} treats the data as having some categorical features and the algorithm 
#' determines which columns contain the categorical features. If all features are 
#' known to be numerical, it will be faster to run SMLE with this argument set to 
#' \code{FALSE}. we will need to find which columns are the categorical features.
#' Default is \code{TRUE}.
#' 
#' @param group Logical flag for whether to treat the dummy covariates of a
#' categorical feature as a group. (Only for categorical data, see Details).
#' Default is \code{TRUE}.
#'
#' @param codingtype Coding types for categorical features; default is \code{"DV"}.
#' \code{codingtype = "all"} convert each level to a 0-1 vector.
#' \code{codingtype = "DV"} conducts deviation coding for each level in
#' comparison with the grand mean.
#' \code{codingtype = "standard"} conducts standard dummy coding for each level
#' in comparison with the reference level (first level).
#'
#' @param coef_initial A \eqn{p}-dimensional vector for the initial coefficient 
#' value of the IHT algorithm.  The default is to use Lasso with the sparsity 
#' closest to \eqn{n-1}. 
#' 
#' @param penalize_mod A logical flag to indicate whether adjustment is used in
#' ranking groups of features. This argument is applicable only when
#' \code{categorical = TRUE} with \code{group = TRUE}. When \code{penalize_mod = TRUE}, 
#' a factor of \eqn{\sqrt J} is divided from the \eqn{L_2} effect of a group with \eqn{J} members. 
#' Default is \code{TRUE}.
#' @param standardize A logical flag for feature standardization, prior to
#' performing feature screening. The resulting coefficients are
#' always returned on the original scale. 
#' If features are in the same units already, you might not wish to
#' standardize. Default is \code{standardize = TRUE}.
#'
#' @param fast Set to \code{TRUE} to enable early stop for SMLE-screening. It may help
#' to boost the screening efficiency with a little sacrifice of accuracy.
#' Default is \code{FALSE}, see Details.
#'
#' @param max_iter Maximum number of iteration steps. Default is 500. 
#'
#' @param tol A tolerance level to stop the iterations, when the squared sum of
#' differences between two successive coefficient updates is below it.
#' Default is \eqn{10^{-3}}.
#' 
#' @param selection A logical flag to indicate whether an elaborate selection
#' is to be conducted by \code{\link{smle_select}()} after screening.
#'  If \code{TRUE}, the function will return a \code{'selection'} object, 
#'  see \code{\link{smle_select}()} documentation. Default is \code{FALSE}.
#'
#' @param ... Additional arguments to be passed to \code{\link{smle_select}()} if \code{selection = TRUE}. 
#' See \code{\link{smle_select}()} documentation for more details. 
#' 
#' @references
#' UCLA Statistical Consulting Group. \emph{coding systems for categorical
#' variables in regression analysis}. \url{https://stats.oarc.ucla.edu/r/library
#' /r-library-contrast-coding-systems-for-categorical-variables/}.
#' Retrieved May 28, 2020.
#'
#' Xu, C. and Chen, J. (2014). The Sparse MLE for Ultrahigh-Dimensional Feature
#' Screening, \emph{Journal of the American Statistical Association}, \bold{109}(507), 1257-1269.
#'
#'
#'
#' @return
#' \item{call}{The call that produced this object.}
#' 
#' \item{ID_retained}{A vector indicating the features retained after SMLE-screening.
#' The output includes both features retained by \code{\link{SMLE}()} and the features specified in \code{keyset}.}
#'
#' \item{coef_retained}{The vector of coefficients estimated by IHT for the retained features. When the 
#' retained set contains a categorical feature, the value returns a group effect if
#' \code{group = TRUE}, or returns the strongest dummy covariate effect if \code{group = FALSE}.}
#'
#' \item{path_retained}{IHT iteration path with columns recording the coefficient updates.}
#'
#' \item{num_retained}{Number of retained features after screening.}
#'
#' \item{intercept}{The estimated intercept value by IHT, if \code{intercept = TRUE}.}
#'
#' \item{steps}{Number of IHT iterations.}
#'
#' \item{likelihood_iter}{A list of log-likelihood updates over the IHT iterations. }
#'
#' \item{Usearch}{A vector giving the number of attempts to find a proper \eqn{1/u} at each iteration step.}
#'
#' \item{modified_data}{A list containing data objects generated by SMLE. 
#'
#' \code{CM}: Design matrix of class \code{'matrix'} for numeric features (or \code{'data.frame'} with categorical features).
#'
#' \code{DM}: A matrix with dummy variable features added. (only if there are categorical features).
#'
#' \code{dum_col}: Number of levels for all categorical features.
#'
#' \code{CI}: Indices of categorical features in \code{CM}.
#'
#' \code{DFI}: Indices of categorical features in \code{IM}.
#' 
#'  }
#' \item{iteration_data}{A list containing data objects that track the coefficients over iterations.
#' 
#' \code{IM}: Iteration path matrix with columns recording IHT coefficient updates.
#' 
#' \code{beta0}: Inital value of regression coefficient for IHT.
#' 
#' \code{feature_name}: A list contains the names of selected features.
#' 
#' \code{FD}: A matrix that contains feature indices retained at each iteration step.
#' } 
#' \code{X}, \code{Y}, \code{data}, \code{family}, \code{categorical} and \code{codingtype} are return of arguments passed in the function call.
#' @export
#' 
#' @examples
#' 
#' # Example 1:
#' set.seed(1)
#' Data <- Gen_Data( n= 200, p = 5000, family = "gaussian", correlation = "ID")
#' fit <- SMLE( Y = Data$Y , X = Data$X, k = 9,family = "gaussian")
#' summary(fit)
#' Data$subset_true %in% fit$ID_retained # Sure screening check.
#' plot(fit)
#' \donttest{
#' # Example 2:
#' set.seed(1)
#' Data_sim2 <- Gen_Data(n = 420, p = 1000, family = "gaussian", num_ctgidx = 5, 
#'                       pos_ctgidx = c(1,3,5,7,9), effect_truecoef= c(1,2,3,-4,-5),
#'                       pos_truecoef = c(1,3,5,7,8), level_ctgidx = c(3,3,3,4,5))
#' train_X <- Data_sim2$X[1:400,]; test_X <- Data_sim2$X[401:420,]
#' train_Y <- Data_sim2$Y[1:400]; test_Y <- Data_sim2$Y[401:420]
#' fit <- SMLE(Y = train_Y, X = train_X, family = "gaussian", group = TRUE, k = 15)
#' predict(fit, newdata = test_X)
#' test_Y
#' 
#' # Example 3:
#' library(datasets)
#' data("attitude")
#' set.seed(1)
#' noise <- matrix(rnorm(30*100, mean = mean(attitude$rating) , sd = 1), ncol = 100)
#' colnames(noise) <- paste("Noise", seq(100), sep = ".")
#' df <- data.frame(cbind(attitude, noise))
#' fit <- SMLE(rating ~., data = df)
#' fit
#' }
#' 
#' 
SMLE <- function(formula=NULL, ...)
  UseMethod("SMLE")

#' @rdname SMLE
#' @export
SMLE.default<-function(formula=NULL, X=NULL,Y=NULL, data=NULL, k=NULL, 
                       family=c("gaussian","binomial","poisson"),
                       keyset = NULL, intercept = TRUE, categorical = TRUE,
                       group = TRUE , codingtype = NULL , coef_initial=NULL,
                       max_iter = 500 , tol = 10^(-3) ,selection = F ,
                       standardize = TRUE , fast = FALSE ,U = 1, U_rate=0.5 ,
                       penalize_mod = TRUE,...){
  
  
  #----------------------------------------------------------------------------#
  #----------------------------Input preprocess--------------------------------#
  #----------------------------------------------------------------------------#
  
  cl <- match.call()
  
  cl[[1]] <- as.name( "SMLE" )
  
  ##---------------------------Data checking----------------------------------##
  
  ### Input object is a design matrix
  
  if( is.null( Y ) & is.null( formula )){
    
    stop( "One of the response variable ('Y') or the formula must be provided. ",call. = FALSE)
    
  }
  
  ### User has provided minimal required arguments
  ### Formula is provided but is not of class(formula)
  
  if( (is.null(X) & is.null(data)) | (is.null(Y) & !is.null( formula)) ){
    
    stop(" The assignment is not clear, please use either 'Y = , X = ' or 'formula = , data = ' to specify the response variable and feature matrix.",call. = FALSE)
    
  }
  
  ##------------------------Argument Checking---------------------------------##
  
  ### Check if input model assumption valid.
  
  if(!inherits(family,"family")){
    
    # Convert it to the family object with default link.
    # 'identity' for Gaussian,'logit' for Binomial and 'log' for Poisson.
    
    if( is.character(family) ) {
      
      # Case 1 that input argument is character.
      
      family <- match.arg( family )
      
      family <- eval(parse(text = paste0(family,"()")))
      
    }else if( inherits( family,'function')){
      
      # Case 2 that input argument is function.
      
      family <- family()
      
    }else{
      
      stop(" Unexpected family input.",call. = FALSE)
      
    }
    
  }
  
  ### Check input model link valid.
  
  if(!((family$family == "gaussian" & family$link=='identity') |
       (family$family == "binomial" & family$link=='logit') |
       (family$family == "poisson" & family$link=='log'))){
    
    stop("Only canonical links are accepted at this time.",call. = FALSE) 
    
  }
  
  ### Check if keyset valid
  
  if( !is.null( keyset ) & !inherits( keyset , c('numeric','integer','character'))){ 
    
    stop( 'Keyset should be numeric, integer or character',call. = FALSE )
    
  }
  
  if(inherits(keyset,'character' )){
    
    keyset <- which( names( X ) %in% keyset )
    
  }
  
  # Check if initial coefficient valid.
  
  if( !is.null( coef_initial ) & length( coef_initial ) != dim( X )[2]){
    
    stop( "Initial coefficient should have the same number of features as feature matrix.",call. = FALSE )
    
    
  }
  
  ##--------------------------Data Preprocess---------------------------------##
  
  ### Backup original data matrix
  
  Data_X <- X
  
  ### Assign defult k
  
  if( is.null( k )){
    
    k <- floor( 1/2 * log( dim( X )[1]) * dim( X )[1] ^ (1/3))
    
  }
  
  ### Check if the data contain categorical features.
  
  CI <- NULL
  
  if(categorical){
    
    ### Determine which features are categorical.
    
    CI <-  which( sapply( X, is.factor ))
    
    if(length(CI) == 0){
      
      CI<-NULL
      
      categorical = FALSE
      
    }
    
  }
  
  if( categorical ){
    
    X_N <- as.matrix(X[,-CI])
    
    colnames(X_N) <- colnames(X)[-CI]
    
  }else{
    
    X_N <-X
    
  }
  
  ### Standardize numeric portion of data matrix X_N.
  if(standardize){
    
    ### Remove categorical features from data matrix before standardizing
    
    Xs <-scale(X_N)
    
    X_mean <- attributes(Xs)$`scaled:center`
    
    X_sd<- attributes(Xs)$`scaled:scale`
    
  }else{
    
    Xs <- X_N
    
    X_mean<- NULL
    
    X_sd <- NULL}
  
  ### If there are categorical features, check coding type and group argument,
  ### and add categorical features them back to the data matrix.
  
  if( categorical ){
    
    if( is.null( codingtype )){
      
      if( group ){ codingtype <- "DV"}else{ codingtype <- "all" }}
    
    if( !group & codingtype != "all" ){
      
      stop("codingtype should be 'all' when group = FALSE",call. = FALSE)
      
    }
    
    for( i in 1:length( CI )){
      
      if( CI[i] == 1 ){ ### First column of X is a factor
        
        name <- colnames( Xs )
        
        Xs <- data.frame( X[,1] , Xs )
        
        colnames( Xs ) <- c( names( X )[1], name)
        
      }else if( CI[i] == dim( Xs )[2] + 1){ ### Last column of X is a factor
        
        name <- colnames( Xs )
        
        Xs <- data.frame( Xs, X[ , CI[i]] )
        
        colnames( Xs )<-c( name, names( X )[CI[i]] )
        
      }else{ ### Putting back the categorical columns that are not the first or last columns.
        
        name <- colnames( Xs )
        
        Xs <- data.frame( Xs[, 1:CI[i]-1 ], X[, CI[i]] ,
                          
                          Xs[, CI[i]:dim( Xs )[2]])
        
        colnames(Xs) <- c( name[1:CI[i]-1], names( X )[ CI[i] ],
                           
                           name[ CI[i]:(length( name ))])
      }
    }
    
  }
  
  #-------------------- Fit model using IHT algorithm -----------------------#

  fit<-SMLE_fit(Y =Y , X = Xs, k = k, family = family, keyset = keyset,
                  
                categorical = categorical, codingtype = codingtype,
                  
                intercept = intercept, max_iter = max_iter, tol = tol, 
                  
                fast = fast,group = group, CI = CI, penalize_mod = penalize_mod,
                  
                U = U, U_rate = U_rate, X_mean = X_mean, X_sd = X_sd, 
                  
                coef_initial = coef_initial, Data_X = Data_X)
  
  
  #------Adding method-----
  
  fit$call = cl
  
  class(fit) = "smle"
  
  #-----Stage II ---------
  if(!selection ){ ### Feature screening only
    
    return(fit)
    
  }else{ ### Feature selection after feature screening
    
    return(smle_select(fit, ...))
  }
}

#' @rdname SMLE
#' @param formula An object of class \code{'formula'} (or one that can be coerced to that class): a symbolic description of the model to be fitted. It should be \code{NULL} when \code{X} and \code{Y} are used.
#' @param data An optional data frame, list or environment (or object coercible by \code{\link[base]{as.data.frame}()} to a \code{'data.frame'}) containing the features in the model. It is required if \code{'formula'} is used.
#' @export
SMLE.formula<- function(formula, data, k=NULL,keyset = NULL, categorical = NULL,...) {
  
  
  # Backup of the vector keyset as the type of the vector will be changed
  
  if( !is.null( keyset ) & !inherits( keyset , c('numeric','integer','character'))){ 
    
    stop( 'Keyset should be numeric, integer or character',call. = FALSE )
    
  }
  
  KEYSET = keyset
  
  if( inherits( keyset , c('integer','numeric')) ){
    
    keyset <- as.character(names(data)[keyset])
    
  }
  
  #--Extract the features from the data matrix that are in the model formula--#
  
  ## keep call (of generic)
  
  cl <- match.call()
  
  cl[[1]] <- as.name("SMLE")
  
  ## Extract model frame
  
  mf <- match.call( expand.dots = FALSE )
  
  m <- c("formula", "data")
  
  m <- match( m, names( mf ), 0L )
  
  mf <- mf[c( 1L, m )]
  mf[[1L]] <- quote( model.frame )
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf, parent.frame())
  
  ## Extract model terms
  
  mt <- attr(mf, "terms")
  
  ## Assign model matrix and response to variables x and y
  
  x <- model.frame(mt, mf, contrasts= NULL)[,-1]
  
  y <- model.response(mf, "numeric")
  
  #---------------------------------------------------------------------------#
  
  key_ind<- as.numeric((1:length(names(x)))[names(x) %in% keyset])
  
  if(is.null(k)){
    
    k<-floor(1/2*log(dim(x)[1])*dim(x)[1]^(1/3))
    
  }
  
  ## Call SMLE on data matrix constructed using the formula object
  
  ans <- SMLE(Y = y, X = x, k = k, keyset = key_ind,...)
  
  ## Adjust the order of formula input.
  
  ans$keyset <- KEYSET
  ans$ID_retained <- as.numeric((1:length(names(data)))[names(data) %in% ans$iteration_data$feature_name])
  ans$iteration_data$feature_name <- names(data)[ans$ID_retained]
  ans$coef_retained<-ans$coef_retained[order(match(names(ans$coef_retained),names(data)))]
  ans$data <- data
  ans$call <- cl
  
  ## Done
  
  return(ans)
}

