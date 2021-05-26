#' Elaborative feature selection with SMLE
#'
#' Elaborative feature selection with SMLE
#'
#' @description Given a response and a set of \code{K} features, this function
#' first runs \code{SMLE (fast=TRUE)} to generate a series of sub-models with
#' sparsity \code{k} varying from \code{k_min} to \code{k_max}.
#' It then selects the best model from the series based on a selection criterion.
#' When criterion EBIC is used, users can choose to repeat the selection with
#' different values of the tuning parameter, \eqn{\gamma}, and
#' conduct importance voting for each feature.
#' 
#' When \code{vote= T}, this function fits all the models with gammas in 
#' \code{gamma_seq} and features whose frequency higher than \code{vote_threshold} 
#' will be selected in \code{ID_voted}.
#'
#' @details
#' This functions accepts three types of input for GLMdata;
#' 1. \code{'smle'} object, as the output from \code{SMLE()};
#' 2. \code{'sdata'} object, as the output from \code{Gen_Data()};
#' 3. Other response and feature matrix input by users.
#'
#' Note that this function is mainly design to conduct an elaborative selection
#' after feature screening. We do not recommend using it directly for
#' ultra-high-dimensional data without screening.
#' @references
#'
#' Chen. J. and Chen. Z. (2012). "Extended BIC for small-n-large-p sparse GLM."
#' \emph{Statistica Sinica}, \bold{22}(2), 555-574.
#'
#' @importFrom parallel detectCores mclapply
#'
#' @param object Object of class \code{'smle'} or \code{'sdata'}. Users can also
#' input a response vector and a feature matrix. See examples
#'
#'
#'
#' @return
#' Returns a \code{'selection'} object with
#' \item{ID_selected}{A list of selected features.}
#' \item{coef_selected}{Fitted model coefficients based on the selected
#' features.}
#' \item{intercept}{Fitted model intercept based on the selected features.}
#' \item{criterion_value}{Values of selection criterion for the candidate models
#' with various sparsity.}
#' \item{ID_voted}{Vector containing the features selected when \code{vote=T}.}
#'
#' @examples
#'
#' # This a simple example for Gaussian assumption.
#' Data<-Gen_Data(correlation="MA",family = "gaussian")
#' fit<-SMLE(Data$Y,Data$X,k=20,family = "gaussian")
#' E<-smle_select(fit,vote=T)
#' summary(E)
#' E1<-smle_select(Data$Y,Data$X[,fit$ID_retained])
#' plot(E1)
#' @export
#'
smle_select<-function(object, ...){
  UseMethod("smle_select")
}
#' @rdname  smle_select
#'
#' @method smle_select smle
#'
#' @export
#'
smle_select.smle<-function(object,...){
  cl<-match.call()
  Data<-structure(list( Y=object$Y, X=object$X, family=object$family),class = "sdata")
  S<-smle_select(Data,subset =object$ID_retained,...)
  S$call <- cl
  return(S)
}

#'
#' @rdname smle_select
#'
#' @method smle_select sdata
#'
#' @importFrom parallel makeCluster stopCluster
#'
#' @param k_min The lower bound of candidate model sparsity. Default is 1.
#'
#' @param k_max The upper bound of candidate model sparsity. Default is as same
#' as the number of columns in input.
#'
#' @param subset A index vector indicating which features (columns of the
#' feature matrix) are to be selected.  Not applicable if a \code{'smle'}
#' object is the input.
#'
#' @param gamma_ebic The EBIC parameter in \eqn{[0 , 1]}. Default is 0.5.
#'
#' @param vote The logical flag for whether to perform the voting procedure.
#' Only available when \code{penalty ='ebic'}.
#'fit
#' @param penalty Selection criterion.One of \code{'ebic'},\code{'bic'},\code{'aic'}. Default is \code{'ebic'}.
#'
#' @param codingtype Coding types for categorical features; details see SMLE.
#'
#' @param gamma_seq The sequence of values for gamma_ebic when \code{vote =TRUE}.
#'
#' @param vote_threshold A relative voting threshold in percentage. A feature is
#'  considered to be important when it receives votes passing the threshold. Default is 0.8.
#'
#' @param parallel Logical flag to use parallel computing to do voting selection.
#' Default is FALSE. See Details.
#'
#' @param num_clusters The number of compute clusters to use. The default will be 2 times cores
#' detected.
#'
#' @export
#'
smle_select.sdata<-function(object, k_min=1, k_max=NULL, subset=NULL,
                            gamma_ebic=0.5, vote= FALSE,
                            penalty="ebic", codingtype = NULL,
                            gamma_seq=c(seq(0,1,0.2)), vote_threshold=NULL,
                            parallel = FALSE, num_clusters=NULL,...){
  cl<-match.call()

  #input check
  X<-object$X

  Y<-object$Y

  family<-object$family

  n<-dim(X)[1]

  pp<- dim(X)[2]
  

  #check input data dimension, warning for High dimensional input.
  if(is.null(subset)){

    if(pp>1000 & pp>n ){
      message("For high dimensional data input, we recommand to do screening before selection")
    }
    X_s<-X

    }else{

      X_s<-X[,subset]


      }

  k_max=dim(X_s)[2]
  
  if ( is.null(vote_threshold) ){

    vote_threshold=0.8

  }
  if(parallel==TRUE && is.null(num_clusters)){ 
    
    num_cores    <- parallel::detectCores() 
    num_clusters <- 2 * num_cores
    }

  IP<-NULL;ID_Voted<-NULL;vs<-NULL

  if(k_min<0 || k_max<k_min ||k_min%%1!=0 ||k_max%%1 !=0||k_max>pp){

    stop("Retained model size setting error(k_min and k_max).")

    }

  #Feature selection by criterion.

  if(is.null( codingtype )){codingtype<-"DV"}

  if(!codingtype%in% c("all","standard","DV")){stop("Codingtype only in
                                                    all,standard,DV")}

  #-----------------------------Algorithm-start---------------------------------
  ctg =FALSE
  
  if( any(  (1:dim(X_s)[2])[sapply(X_s,is.factor)] )){
    #if sub matrix X_s has any categorical features
    ctg= TRUE
    
    criter_value<-ctg_ebicc(Y,X_s,family,penalty,codingtype,
                            k_min,k_max,n,pp,gamma_ebic,parallel,num_clusters)

    v_s <- which.min(criter_value)+k_min-1

    f_s<- SMLE(Y=Y, X=X_s, k=v_s, family=family,categorical = T)

    #Feature selection by ebic voting.

    ebic_selection_by_gamma<-function(gamma){

      ebic_value<-ctg_ebicc(Y,X_s,family,penalty,codingtype,
                            k_min,k_max,n,pp,gamma,parallel,num_clusters)

      v_s <- which.min(ebic_value)+k_min-1

      return(SMLE(Y=Y, X=X_s, k=v_s, family=family,categorical = T)$ID_retained)

    }

  }else{


    criter_value<-ebicc(Y,X_s,family,penalty,k_min,k_max,n,pp,gamma_ebic,parallel,num_clusters)

    v_s <- which.min(criter_value)+k_min-1

    f_s<- SMLE(Y=Y, X=X_s, k=v_s, family=family)

    #Feature selection by ebic voting.

    ebic_selection_by_gamma<-function(gamma){
      library(SMLE)

      ebic_value<-ebicc(Y,X_s,family,penalty,k_min,k_max,n,pp,gamma,parallel,num_clusters)

      v_s <- which.min(ebic_value)+k_min-1

      return(SMLE(Y=Y, X=X_s, k=v_s, family=family)$ID_retained)

    }

  }


  if(vote==T){

    stopifnot(penalty == 'ebic')

    vs<-c()

    if(parallel==TRUE){
      if(.Platform$OS.type=="windows"){
        
        cl<-parallel::makeCluster(num_clusters)
        
      }
      
        vs<-unlist(mclapply(gamma_seq,ebic_selection_by_gamma,mc.cores = 2*num_clusters))
        
        }else{
      vs<-unlist(lapply(gamma_seq,ebic_selection_by_gamma))

      }
    IP<-as.factor(vs)
    
    IP_f<-summary(IP)[order(summary(IP),decreasing= T)]/max(summary(IP))
    
    ID_Voted<-as.numeric(names(IP_f[IP_f>=vote_threshold]))
    #ID_Voted<-as.numeric(names(summary(IP)[order(summary(IP),decreasing= T)[1:min(length(summary(IP)),vote_threshold)]]))

    }
  if(is.null(subset)){

    ID_Selected=f_s$ID_retained

  }else{

    ID_Selected<-subset[f_s$ID_retained]

    ID_Voted<-subset[ID_Voted]

  }

  S<-list(
          call = cl,
    
          ID_selected=ID_Selected, family = family,

          coef_selected=f_s$coef_retained,
          
          intercept = f_s$intercept,
          
          num_selected = length(f_s$coef_retained),

          vote=vote,criterion=penalty,
          
          ID_pool= IP,
          
          subset = subset,

          criterion_value=criter_value,

          ID_voted=ID_Voted,
          
          ctg = ctg,vote_threshold=vote_threshold,

          gamma_ebic=gamma_ebic,
          
          codingtype = codingtype,

          gamma_seq=gamma_seq, X = X , Y = Y)
  class(S)<-"selection"
  return(S)
  }

#' @rdname smle_select
#'
#' @method smle_select default
#'
#' @param X Input features matrix. When feature matrix input by users.
#'
#' @param family Model assumption; see SMLE. Default is Gaussian linear.
#'
#' When input is \code{'smle'} or \code{'sdata'}, the same
#' model will be used in the selection.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#'
smle_select.default<-function(object, X=NULL, family='gaussian',...){
  cl<-match.call()
  if ( any(sapply(X,as.factor) ==TRUE) ){
   ctg = TRUE
  }else{
   ctg =FALSE
    }

  Data<-structure(list(Y=object,X=X,family=family,ctg = ctg),class = "sdata")


  S<-smle_select(Data,...)
  S$call <- cl
  return(S)
}







