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
#' In \code{\link{SMLE}()}, the initial value for parameter \eqn{1/u} is set to be
#' \eqn{1/||X||^2_{\infty}} and recursively decrease the value of 
#' \eqn{1/u} by \code{U_rate} to guarantee the likelihood increment.This strategy is called \eqn{u}-search.
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
#' When \code{group=TRUE} with \code{penalize_mod = TRUE}, the effect for a group
#' of \eqn{J} dummy covariates is computed by
#'
#' \deqn{ \beta_i = \sqrt{(\beta_1)^2+...+(\beta_J)^2}/\sqrt J,}
#'
#' which will be treated as a single feature in IHT iterations. When \code{group = FALSE}, 
#' a group of \eqn{J} dummy covariates will be treated as \eqn{J} individual features in the IHT iterations;  in this case, 
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
#' @param family Model assumption between \eqn{Y} and \eqn{X}; the default model is Gaussian
#' linear.
#'
#' @param categorical A logical flag whether the input feature matrix includes
#' categorical features. If \code{categorical = TRUE}, a model intercept will
#' be used in the screening process. Default is \code{NULL}.
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
#' @param group Logical flag for whether to treat the dummy covariates of a
#' categorical feature as a group. (Only for categorical data, see Details).
#' Default is \code{TRUE}.
#'
#' @param codingtype Coding types for categorical features; default is \code{"DV"}.
#' \code{Codingtype = "all"} Convert each level to a 0-1 vector.
#' \code{Codingtype = "DV"} conducts deviation coding for each level in
#' comparison with the grand mean.
#' \code{Codingtype = "standard"} conducts standard dummy coding for each level
#' in comparison with the reference level (first level).
#'
#' @param coef_initial A \eqn{p}-dimensional vector for the initial coefficient 
#' value of the IHT algorithm.  The default is to use Lasso with the sparsity 
#' closest to \eqn{n-1}. 
#' 
#' @param penalize_mod A logical flag to indicate whether adjustment is used in
#' ranking groups of features. This argument is applicable only when
#' \code{categorical= TRUE} with \code{group=TRUE}. When \code{penalize_mod=TRUE}, 
#' a factor of \eqn{\sqrt J} is divided from the \eqn{L_2} effect of a group with \eqn{J} members. 
#' Default is \code{TRUE}.
#' @param standardize A logical flag for feature standardization, prior to
#' performing feature screening. The resulting coefficients are
#' always returned on the original scale. 
#' If features are in the same units already, you might not wish to
#' standardize. Default is \code{standardize=TRUE}.
#'
#' @param fast Set to \code{TRUE} to enable early stop for SMLE-screening. It may help
#' to boost the screening efficiency with a little sacrifice of accuracy.
#' Default is \code{FALSE}, see Details.
#'
#' @param max_iter Maximum number of iteration steps. Default is 500. 
#'
#' @param tol A tolerance level to stop the iterations, when the squared sum of
#' differences between two successive coefficient updates is below it.
#' Default is \eqn{10^{-2}}.
#' 
#'
#' @param selection A logical flag to indicate whether an elaborate selection
#' is to be conducted by \code{\link{smle_select}()} after screening.
#'  If \code{TRUE}, the function will return a \code{'selection'} object, 
#'  see \code{\link{smle_select}()} documentation. Default is \code{FALSE}.
#'
#' @param ... Additional arguments to be passed to \code{\link{smle_select}()} if \code{selection=TRUE}. 
#' See \code{\link{smle_select}()} documentation for more details. 
#' 
#' @references
#' UCLA Statistical Consulting Group. \emph{coding systems for categorical
#' variables in regression analysis}. \url{https://stats.idre.ucla.edu/spss
#' /faq/coding-systems-for-categorical-variables-in-regression-analysis-2/}.
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
#' \code{group=TRUE}, or returns the strongest dummy covariate effect if \code{group=FALSE}.}
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
#' Data <- Gen_Data(n= 200, p = 5000, family = "gaussian", correlation = "ID")
#' fit <- SMLE(Y = Data$Y, X = Data$X, k = 9, family = "gaussian")
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
SMLE.default<-function(formula=NULL, X=NULL,Y=NULL,data=NULL, k=NULL, 
                       family=c("gaussian","binomial","poisson"),
                       categorical = NULL , keyset = NULL, intercept = TRUE ,
                       group = TRUE , codingtype = NULL , coef_initial=NULL,
                       max_iter = 500 , tol = 0.01 ,selection = F ,
                       standardize = TRUE , fast = FALSE , U_rate=0.5 ,
                       penalize_mod = TRUE,...){
  
  #-------Input preprocess-------
  
  family<-match.arg(family)
  
  cl<-match.call()
  cl[[1]] <- as.name("SMLE")
  
  ####Input object is a design matrix
  if(is.null(Y)&is.null(formula)){stop("Response required")}
  if(is.null(X)&is.null(data)){stop("Feature Matrix required")}
  if(is.null(Y)&!is.null(formula)){Y<-formula}
  if(is.null(X)&!is.null(data)){X<-data}
  if(!is.null(keyset)& class(keyset)!='numeric'& class(keyset)!='integer' & class(keyset)!='character'){stop('Keyset should be numeric, integer or character')}
  if(class(keyset)=='character'){keyset <-(1:length(names(X)))[names(X) %in% keyset]}
  Data_X<-X
  ####
  if(is.null(k)){
    
    k<-floor(1/2*log(dim(X)[1])*dim(X)[1]^(1/3))
  }
  
  if(is.null(categorical)){
    
    CI<-(1:dim(X)[2])[sapply(X,is.factor)]
    
    if( any(CI) ){
      
      categorical = TRUE
      
    }else{
      
      categorical= FALSE
      
      X<- as.matrix(X, ncol=dim(X)[2])
      
    }
    
  }
  
  
  #-------Run Algoriathm------
  
  if(!is.null(coef_initial)){
    
    stopifnot(length(coef_initial)==dim(X)[2])
  }
  
  
  if(categorical == TRUE){
    
    CI<-(1:dim(X)[2])[sapply(X,is.factor)]
    
    if(is.null( codingtype )){if(group==TRUE){codingtype<-"DV"}else{codingtype<-"all"}}
    
    if(group==FALSE & codingtype!="all"){stop("codingtype should be 'all' when group = FALSE")}
    
    X_N <-as.matrix(X[,-CI])
    
    center = colMeans(X_N)
    
    X.c = sweep(X_N, 2, center)
    
    unit.var = apply(X.c, 2, sd)
    
    if(sum(unit.var==0)!=0){
      zero_var = (1:dim(X)[2])[unit.var==0]
      X.c<-X.c[,-zero_var]
      unit.var<-unit.var[-zero_var]
    }
    
    Xs = as.matrix(sweep(X.c, 2, unit.var, "/"))
    
    colnames(Xs) <- names(X[-CI])
    
    X_mean<-center
    
    X_sd<- unit.var
      
      for(i in 1:length(CI)){
      
        if( CI[i] == 1 ){
        
          name<-colnames(Xs)
        
          Xs <- data.frame( X[,1] , Xs )
        
          colnames(Xs)<-c(names(X)[1],name)
        
       }else if( CI[i]>dim(Xs)[2]){
        
          name<-colnames(Xs)
        
          Xs<-data.frame(Xs,X[,CI[i]])
        
          colnames(Xs)<-c(name,names(X)[CI[i]])
        
        }else{ 
        
          name<- colnames(Xs)
        
          Xs<-data.frame(Xs[,1:CI[i]-1],X[,CI[i]],Xs[,CI[i]:dim(Xs)[2]])
        
         colnames(Xs)<-c(name[1:CI[i]-1],names(X)[CI[i]],name[CI[i]:(length(name))])
        }
      }
  
    
    fit<- ctg_fit(Y, Xs, k, family, keyset, CI, max_iter, tol, fast,
                  
                  intercept, group, codingtype, penalize_mod,
                  
                  U_rate, X_mean , X_sd, coef_initial,Data_X)
    
    fit$call<-cl
    
    if(selection==F){
      
      return(fit)
      
    }else{
      
      return(smle_select(fit, ...))
    }
        return(fit)
    
  }else{
    
    X_mean = NULL;X_sd=NULL;Xs<-X
    
    if(standardize== TRUE){
      
      name<-dimnames(X)
      
      center = colMeans(X)
      
      X.c = sweep(X, 2, center)
      
      unit.var = apply(X.c, 2, sd)
      
      if(sum(unit.var==0)!=0){
        zero_var = (1:dim(X)[2])[unit.var==0]
        X.c<-X.c[,-zero_var]
        unit.var<-unit.var[-zero_var]
      }
      
      val <- sweep(X.c, 2, unit.var, "/")
      
      Xs <- as.matrix(val,dimnames=name)

      X_mean<-center
      
      X_sd<- unit.var
      
    }
    
    fit<-SMLE_fit(Y=Y, X= Xs, k= k, family=family, keyset=keyset,
                  
                  intercept=intercept, max_iter =max_iter, tol =tol, fast =fast,
                  
                  U_rate=U_rate, X_mean = X_mean, X_sd=X_sd, coef_initial=coef_initial, Data_X=Data_X)
  }
  
  #------Adding method-----
  fit$call=cl
  
  class(fit)="smle"
  
  #-----Stage II ---------
  if(selection==F){
    
    return(fit)
    
  }else{
    
    return(smle_select(fit, ...))
  }
}

#' @rdname SMLE
#' @param formula An object of class \code{'formula'} (or one that can be coerced to that class): a symbolic description of the model to be fitted. It should be \code{NULL} when \code{X} and \code{Y} are used.
#' @param data An optional data frame, list or environment (or object coercible by \code{\link[base]{as.data.frame}()} to a \code{'data.frame'}) containing the features in the model. It is required if \code{'formula'} is used.
#' @export
SMLE.formula<- function(formula, data, categorical=NULL, k=NULL,keyset = NULL, ...) {
  
  if(!is.null(keyset)& class(keyset)!='numeric'& class(keyset)!='integer' & class(keyset)!='character'){stop('Keyset should be numeric, integer or character')}
  KEYSET = keyset
  if(class(keyset)=='integer'|class(keyset)=='numeric'){keyset <-as.character(names(data)[keyset])}
  
  ## keep call (of generic)
  cl <- match.call()
  cl[[1]] <- as.name("SMLE")
  ## model frame
  mf <- match.call(expand.dots = FALSE)
  m <- c("formula", "data")
  m <- match(m, names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(model.frame)
  mf$drop.unused.levels <- TRUE
  mf <- eval(mf, parent.frame())
  
  ## model terms
  mt <- attr(mf, "terms")
  
  ## model matrix and response
  x <- model.frame(mt, mf, contrasts= NULL)[,-1]
  y <- model.response(mf, "numeric")
  
  keyset<- as.numeric((1:length(names(x)))[names(x) %in% keyset])
  
  
  if(is.null(k)){
    
    k<-floor(1/2*log(dim(x)[1])*dim(x)[1]^(1/3))
  }
  
  
  if(is.null(categorical)){
    
    CI<-(1:dim(x)[2])[sapply(x,is.factor)]
    
    if( any(CI) ){
      
      ans <- SMLE(Y = y, X = x, categorical = TRUE, k = k, keyset = keyset,...)
      
    }else{
      
      ## fit subsets
      ans <- SMLE(Y = y, X = as.matrix(x, ncol=dim(x)[2]),categorical = FALSE, k = k, keyset = keyset, ...)
    }
  }else{
    
    ans <- SMLE(Y = y, X = x,categorical = categorical,  k = k, keyset = keyset, ...)  
    
    }
  ## Adjust the order of formula input.
  ans$keyset <- KEYSET
  ans$ID_retained <- as.numeric((1:length(names(data)))[names(data) %in% ans$iteration_data$feature_name])
  ans$iteration_data$feature_name <- names(data)[ans$ID_retained]
  ans$coef_retained<-ans$coef_retained[order(match(names(ans$coef_retained),names(data)))]
  ans$data <- data
  ## done
  ans
}

SMLE_fit<-function(Y,X, k, family="gaussian", keyset=NULL,
                   intercept=TRUE, max_iter=500, tol=0.01, fast=FALSE,
                   U_rate=0.5,X_mean=NULL,X_sd=NULL,coef_initial=NULL, Data_X){
  
  LH<-rep(0,max_iter)
  
  number_of_Ucheck<-rep(0,max_iter)
  
  n<-dim(X)[1];p<-dim(X)[2]
  
  if( is.null( coef_initial ) ){
    
    fit_pre<-glmnet(x=X,y=Y,family=family)
    
    Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
    
  }else{ 
    
    Beta0 <- coef_initial }
  
  FD<-NULL
  
  if(is.null(keyset)){
    
    number_of_ID_retained <-k
    
  }
  else{
    
    number_of_ID_retained<-k-(length(keyset)) 
    
  }
  if(is.null(Data_X)){Data_X =X}
  if(intercept==T){
    
    Beta0<-c(mean(Y-tcrossprod(as.matrix(X),t(Beta0))),Beta0)
    
    names(Beta0)[1]="(intercept)"
    
    X_iter<-cbind(matrix(1,nrow  = n, ncol = 1),X)
    
    ID_None0<-(1:(p+1))[Beta0!=0]
    
    keyset<-c(1,keyset+1)
    
  }else{
    
    X_iter<-X
    
    ID_None0<-(1:p)[Beta0!=0]
    
  }
  
  pp<-p+intercept
  
  Screening_index<-sub_off(1:pp,keyset)
  
  beta_path<-as.matrix(Beta0,nrow=pp,ncol=1)
  
  I<-list(CM=X,Y=Y,IM=X_iter,beta0=Beta0,family=family)
  
  pp<-p+intercept
  
  coef_None0 <- as.matrix(Beta0[ID_None0],ncol=1)
  
  Xs_0 <- X_iter[, ID_None0]
  
  if(!is.null(coef_initial) & sum(coef_initial!=0)==0){ 
    
    R_0<-matrix(0, ncol=1, nrow=n) 
    
  }else{
    R_0  <- Xs_0 %*% coef_None0
  }
  
  
  R_0<-switch(family,
              
              "gaussian"=R_0,
              
              "poisson"=exp(R_0),
              
              'binomial'=exp(R_0)/(1+exp(R_0)))
  
  V_0<-crossprod(X_iter,Y - R_0)
  
  
  if(!is.null(coef_initial) & sum(coef_initial!=0)==0){
    
    if(intercept==T){
      
      #When starts from zero with uu is 1/||X||\infit = 1/sqrt(n) 
      uu<-1/(sqrt(n))
      
      
    }else{
      
      uu<- 1/max(colSums(Xs_0^2))
      
    }
    
    
  }else{    
    
    uu<-1/max(colSums(Xs_0^2))
    
  }
  
  
  
  U_i <- uu
  
  #iteration start form 1
  i<-1
  
  Beta_s<-Beta0
  
  
  #######################################################################
  # iteration start---------------------------
  
  
  repeat{
    
    count<-0
    
    repeat{
      
      Beta_t<- Beta_s + uu * V_0
      
      Beta_t[Screening_index]<-Hard(t=Beta_t[Screening_index],k= number_of_ID_retained)
      
      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)
      

      if (ucheck >= 1){
        
        break
        
      }else{
        
        uu <- U_rate * uu

        count<-count+1
        
      }
      
    }
    
    sindex<-(1:pp)[Beta_s!= 0]
    
    tindex<-(1:pp)[Beta_t!= 0]
    
    fs<-sum(!(tindex %in% sindex))
    
    FD<-c(FD,fs)
    
    beta_path<-cbind(beta_path,as.matrix(Beta_t,ncol=1,nrow=pp))
    
    LH[i]<-lh(Y, X_iter, Beta_t,family=family)
    
    valid_LH_diff<- 0.01 *(LH[2] - LH[1])
    
    number_of_Ucheck[i]<-count
    
    
    ######## convergence check #######################################
    if(i>1){
      
      MSE<- sqrt(sum(( Beta_s-Beta_t )^2))
      
      if(fast == TRUE){
        
        
        if( MSE/number_of_ID_retained < tol){break}
        else if( (LH[i]-LH[i-1])< valid_LH_diff){break}
        else if(i>10){if(sum(tail(FD,10))==0){break}}
        
      }else{
        
        if(MSE < tol || i >= max_iter){
          break}
      }
      
      
    }
    
    ##################################################################
    
    
    Beta_s<-Beta_t
    
    ID_None0<- sort((1:pp)[Beta_s!= 0])
    
    coef_None0<-as.matrix(Beta_s[ID_None0],ncol=1)
    
    Xs_0<-X_iter[, ID_None0]
    
    R_0<-crossprod(t(Xs_0),coef_None0)
    
    R_0<-switch(family,
                
                "gaussian"=R_0,
                
                "poisson"=exp(R_0),
                
                'binomial'=exp(R_0)/(1+exp(R_0))
    )
    
    V_0<-crossprod(X_iter,Y - R_0)
    
    uu<-U_i
    
    i<-i+1
    
  }
  
  
  Coef_dist<-apply((beta_path[,-1]-beta_path[,-ncol(beta_path)]),2,function(x){sqrt(sum(x^2))})
  
  # Rescale all output
  
  if( intercept == T ){
    
    #resale beta_path
    
    beta_path<-beta_path[-1,]             
    
    Intercept_value<- as.vector(coef_None0)[1]
    
    ID_None0<-ID_None0[-1]-1
    
    coef_None0 <- as.vector(coef_None0)[-1]
    
    if(!is.null(X_mean)){
      
      coef_None0<- coef_None0/X_sd[ID_None0]
      
      Intercept_value <- Intercept_value-X_mean[ID_None0]%*%coef_None0
      
    }
  }else{
    
    if(!is.null(X_mean)){
      
      coef_None0<- coef_None0/X_sd[ID_None0]
      
      beta_path<- beta_path/X_sd
    }
    
    Intercept_value = NULL
    
  }
  
  feature_name =NULL
  
  if( !is.null(colnames(Data_X))){ feature_name=colnames(Data_X)[ID_None0]}
  
  fit<-list(iteration_data = list(IM = I$IM, beta0 = Beta0, FD = FD,            
                                 
                                   feature_name =feature_name ),
            X=Data_X, Y=Y, 
            
            keyset = keyset, family = family, k = k,
            
            intercept=Intercept_value,
            
            steps = i,
            
            likelihood_iter=LH[1:i],
            
            Usearch=number_of_Ucheck[1:i],
            
            path_retained  = beta_path,
            
            num_retained   = length(ID_None0),
            
            ID_retained    = ID_None0,
            
            coef_retained  = coef_None0,
            
            categorical = FALSE, coef_dist= Coef_dist,
            
            fast=fast
  )
  fit
}

ctg_fit<-function(Y , X , k ,
                  
                  family = c("gaussian","binomial","poisson"), keyset ,CI, 
                  
                  max_iter, tol ,fast,
                  
                  intercept , group ,
                  
                  codingtype , penalize_mod,
                  
                  U_rate , X_mean , X_sd, coef_initial=NULL,Data_X=NULL){
  
  family<-match.arg(family)
  
  if(is.null(Data_X)){Data_X =X}
  n<-dim(X)[1];p<-dim(X)[2]
  
  #--------------------------------------------------------------#
  
  if(group == FALSE){if(sum(CI %in% keyset)){stop("Group shoule be TRUE if keyset contains categorical features.")}}
  
  
  if(length(CI)== 1){
    
    if(codingtype=="all"){
      
      dum_col = nlevels(X[,CI])
      
    }else{ 
        dum_col = nlevels(X[,CI])-1
    }
  }else{
      if(codingtype=="all"){
        
        dum_col = sapply(X[,CI],nlevels)
        
      }else{ 
        dum_col = sapply(X[,CI],nlevels)-1
      }
    }
                               

  X_dummy<-suppressWarnings(dummy.data.frame(X,sep="_",codingtype = codingtype))
  
  
  
  #--------------------------------------------------------------#
  
  #Tracking index in Dummy working matrix
  
  Dummy_index<-c()
  
  Dummy_sum<-0
  
  for(i in 1:length(CI)){
    
    Dummy_index<-c(Dummy_index,list(CI[i]-1+seq(dum_col[i])+Dummy_sum))
    
    Dummy_sum<-Dummy_sum+dum_col[i]-1
    
  }
  
  DFI<-Dummy_index
  
  DI<-unlist(lapply(DFI, function(l) l[[1]]))
  
  #--------------------------------------------------------------#
  
  
  if( is.null( coef_initial ) ){
    
    fit_pre<-glmnet(x=as.matrix(X_dummy,dimnames = dimnames(X_dummy)),y=Y,family=family)
    
    Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
    
  }else{ 
    
    Beta0<- rep(0,dim(X_dummy)[2])
    
    Beta0[-unlist(DFI)] <- coef_initial[-CI]
    
    for(i in 1:length(CI)){
      
      Beta0[Dummy_index[[i]]]<-rep(coef_initial[CI[i]],dum_col[i])
      
    }
    
  }
  
  
  Beta0<-c(mean(Y-tcrossprod(as.matrix(X_dummy),t(Beta0))),Beta0)
  
  names(Beta0)[1]="(intercept)"
  
  
  #--------------------------------------------------------------#
  
  X_iter<-as.matrix(cbind(matrix(1,nrow  = n, ncol = 1),X_dummy))
  
  colnames(X_iter)[1]<-'(intercept)'
  
  pp<-dim(X_iter)[2]
  
  I<-list(Y=Y,CM=X,CI=sort(CI),dum_col=dum_col,IM=X_iter,
          DFI= lapply(DFI,function(x) x+1),DI= DI+1,family=family,codingtype=codingtype)
  
  
  #--------------------------------------------------------------#
  
  
  
  
  #--------------------------------------------------------------#
  
  ID_None0<-sort((1:pp)[Beta0!=0])  #length : = pp
  
  coef_None0 <- as.matrix(Beta0[ID_None0],ncol=1)
  
  Xs_0 <- X_iter[, ID_None0]
  
  
  if(!is.null(coef_initial) & sum(coef_initial!=0)==0){ 
    
    R_0<-matrix(0, ncol=1, nrow=n) 
    
  }else{
    R_0  <- Xs_0 %*% coef_None0
  }
  
  R_0<-switch(family,
              "gaussian"=R_0,
              "poisson"=exp(R_0),
              'binomial'=exp(R_0)/(1+exp(R_0))
  )
  
  
  V_0<-crossprod(X_iter,Y - R_0)
  
  if(!is.null(coef_initial) & sum(coef_initial!=0)==0){
    
    if(intercept==T){
      
      #When starts from zero with uu is 1/||X||\infit = 1/sqrt(n) 
      
      uu<- 1/max(colSums(as.matrix(Xs_0^2,nrow=n)))
      
      
    }else{
      
      uu<-1/(sqrt(n))
      
    }
    
    
  }else{    
    
    uu<-1/max(colSums(Xs_0^2))
    
  }
  
  ###########################################################
  
  # iteration start---------------------------
  U_i<- uu
  
  i<-1
  
  Beta_s<-Beta0# starting iteration.
  
  LH<-rep(0,max_iter)
  
  number_of_Ucheck<-rep(0,max_iter)
  
  FD<-NULL
  
  Screening_index<-sub_off(1:p,keyset)
  
  beta_path<-as.matrix(Group_Beta(Beta0,I,penalize_mod),nrow=pp,ncol=1)
  
  Screening_Dindex<-sub_off(2:pp,CI2DI(I,keyset))
  
  if(is.null(keyset)){ number_of_ID_retained <- k }else{ 
    
    
    number_of_ID_retained <- k-(length(keyset))
    
    if(number_of_ID_retained <= 0){stop("The screening size k should be larger than the number of features in keyset.")}
    
    }
  
  repeat{
    
    count<-0
    
    repeat
    {
      
      
      Beta_t<-Beta_s + uu * V_0         # length(Beta_s)  = pp
      
      if(group==T){
        
        Beta_t<-GroupHard(Beta_t,I,k = k-(length(keyset)),Screening_index,penalize_mod)
        # length(Beta_t) = pp
      }else{
        
        Beta_t[Screening_Dindex]<- Hard(t=Beta_t[Screening_Dindex], k=number_of_ID_retained)
        # length(Beta_t) = pp
      }
      
      
      ########## u-check #################
      
      ucheck<- Uh(A=X_iter, uh=uu, b0=Beta_s, b1=Beta_t, family=family)

      if(ucheck >= 1)
      {
        break
        
      }else{
        
        uu <- U_rate * uu
        
        count<-count+1
        
      }
    }
    
    sindex<-(1:pp)[Beta_s!= 0]
    
    tindex<-(1:pp)[Beta_t!= 0]
    
    fs<-sum(!(tindex %in% sindex))
    
    FD<-cbind(FD,fs)
    
    beta_path<-cbind(beta_path,as.matrix(Group_Beta(Beta_t,I,penalize_mod),ncol=1,nrow=pp))
    
    LH[i]<-lh(Y, X_iter, Beta_t,family=family)
    
    number_of_Ucheck[i]<-count
    
    ######## convergence check ##############
    
    
    if(i>1){
      
      if(fast == TRUE){
        
        MSE<- sqrt(sum((Beta_s-Beta_t)^2))/k
        
        if(  (length(FD)>10 & sum(tail(FD,10)) )|| ((LH[i]-LH[i-1])< 0.05*LH[2]-LH[1])  ){break}
        
      }else{
        
        MSE<- sqrt(sum((Beta_s-Beta_t)^2))
      }
      
      if((i>=max_iter)||(MSE < tol)) {break}
      
    }
    #########################################
    
    Beta_s<-Beta_t
    
    ID_None0<- sort((1:pp)[Beta_s!= 0])
    
    coef_None0 <- as.matrix(Beta_s[ID_None0],ncol=1)
    
    Xs_0<-X_iter[, ID_None0]
    
    R_0<-crossprod(t(Xs_0), coef_None0)
    
    R_0<-switch(family,
                "gaussian"=R_0,
                "poisson"=exp(R_0),
                'binomial'=exp(R_0)/(1+exp(R_0))
    )
    
    V_0<-crossprod(X_iter,Y - R_0)
    
    uu<-U_i
    
    i<-i+1
  }
  
  
  Coef_dist<-apply((beta_path[,-1]-beta_path[,-ncol(beta_path)]),2,function(x){sqrt(sum(x^2))})
  
  Intercept_value<-coef_None0[1]
  
  ID_None0 <- DI2CI(I,ID_None0)     
  # length pp -> p
  
  if(group == TRUE){  coef <- Group_Beta(Beta_s,I,penalize_mod)}else{ coef <- Group_Beta(Beta_s,I,max=T) }

  coef_None0 <- coef[coef!=0][-1]
    
  feature_name =NULL
    
  if( !is.null(colnames(Data_X))){ 
      
    feature_name=colnames(Data_X)[ID_None0]
      
  }else{
      
    feature_name = as.character(ID_None0)}
    
    names(coef_None0)<-feature_name
    
    
  if(length(ID_None0) < k){message("Warning: Retained model size is less than k as multiple dummy covariates of a categorical feature may be retained when group=FALSE.")}
  
  fit<-list(modified_data = list(CM = I$CM, CI = sort(I$CI), dum_col = I$dum_col, DFI = I$DFI),
            
            iteration_data = list(IM = I$IM, beta0 = Beta0, FD = FD,            
                                  feature_name =feature_name ),
            X=Data_X,Y=I$Y,
            
            ID_retained=ID_None0,
            
            coef_retained = coef_None0, keyset = keyset,
            
            Usearch=number_of_Ucheck[1:i],
            
            family = family, k = k,
            
            intercept = Intercept_value, steps = i,
            
            likelihood_iter=LH[1:i],
            
            path_retained=beta_path[-1,],
            
            num_retained = length(ID_None0),
            
            group=group, fast =FALSE,
            
            categorical =TRUE, coef_dist=Coef_dist,
            
            codingtype = codingtype,
            
            fast=fast
  )
  
  class(fit)="smle"
  fit
  
}

