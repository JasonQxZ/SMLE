#' @title
#' Elaborative post-screening selection with SMLE
#'
#' @description
#' The features retained after screening are still likely to contain some that 
#' are not related to the response. The function \code{\link{smle_select}()} is designed to 
#' further identify the relevant features using \code{\link{SMLE}()}.
#' Given a response and a set of \eqn{K} features, this function
#' first runs \code{\link{SMLE}(fast = TRUE)} to generate a series of sub-models with
#' sparsity k varying from \code{k_min} to \code{k_max}.
#' It then selects the best model from the series based on a selection criterion.
#' 
#' When criterion EBIC is used, users can choose to repeat the selection with
#' different values of the tuning parameter \eqn{\gamma}, and
#' conduct importance voting for each feature. When \code{vote = T}, this function 
#' fits all the models with \eqn{\gamma} specified in \code{gamma_seq} and features 
#' with frequency higher than \code{vote_threshold} will be selected in \code{ID_voted}.
#'
#' @details
#' This function accepts three types of input objects; 
#' 1) \code{'smle'} object, as the output from \code{\link{SMLE}()}; 
#' 2) \code{'sdata'} object, as the output from \code{\link{Gen_Data}()}; 
#' 3) other response and feature matrix input by users.
#'
#' Note that this function is mainly designed to conduct an elaborative selection
#' after feature screening. We do not recommend using it directly for
#' ultra-high-dimensional data without screening.
#' 
#' @references
#' Chen. J. and Chen. Z. (2012). "Extended BIC for small-n-large-p sparse GLM."
#' \emph{Statistica Sinica}, \bold{22}(2), 555-574.
#'
#' @importFrom parallel detectCores mclapply
#'
#' @param object Object of class \code{'smle'} or \code{'sdata'}. Users can also
#' input a response vector and a feature matrix. 
#'
#'
#'
#' @return
#' \item{call}{The call that produced this object.}
#' \item{ID_selected}{A list of selected features.}
#' \item{coef_selected}{Fitted model coefficients.}
#' \item{intercept}{Fitted model intercept.}
#' \item{criterion_value}{Values of selection criterion for the candidate models with various sparsity.}
#' \item{categorical}{A logical flag whether the input feature matrix includes categorical features}
#' \item{ID_pool}{A vector containing all features selected during voting. }
#' \item{ID_voted}{A vector containing the features selected when \code{vote = T}.}
#' \item{CI}{Indices of categorical features when \code{categorical = TRUE}.}
#' \code{X}, \code{Y}, \code{family}, \code{gamma_ebic}, \code{gamma_seq}, \code{criterion}, \code{vote},
#'  \code{codyingtype}, \code{vote_threshold} are return of arguments passed in the function call.
#' @examples
#'
#' set.seed(1)
#' Data<-Gen_Data(correlation = "MA", family = "gaussian")
#' fit<-SMLE(Y = Data$Y, X = Data$X, k = 20, family = "gaussian")
#' 
#' fit_bic<-smle_select(fit, criterion = "bic")
#' summary(fit_bic)
#' 
#' fit_ebic<-smle_select(fit, criterion = "ebic", vote = TRUE)
#' summary(fit_ebic)
#' plot(fit_ebic)
#' 
#' 
#' @export
#'
smle_select<-function(object, ...){
  UseMethod("smle_select")
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
#' @param k_max The upper bound of candidate model sparsity. Default is 
#' the number of columns in feature matrix.
#'
#' @param subset An index vector indicating which features (columns of the
#' feature matrix) are to be selected.  Not applicable if a \code{'smle'}
#' object is the input.
#'
#' @param gamma_ebic The EBIC tuning parameter, in \eqn{[0 , 1]}. Default is 0.5.
#'
#' @param vote The logical flag for whether to perform the voting procedure. Only available when \code{criterion = "ebic"}.
#' 
#' @param criterion Selection criterion. One of "\code{ebic}","\code{bic}","\code{aic}". Default is "\code{ebic}".
#'
#' @param codingtype Coding types for categorical features; for more details see \code{\link{SMLE}()} documentation.
#'
#' @param gamma_seq The sequence of values for \code{gamma_ebic} when \code{vote = TRUE}.
#'
#' @param vote_threshold A relative voting threshold in percentage. A feature is
#'  considered to be important when it receives votes passing the threshold. Default is 0.6.
#'
#' @param parallel A logical flag to use parallel computing to do voting selection.
#' Default is \code{FALSE}. See Details.
#'
#' @param keyset A numeric vector with column indices for the key features that 
#' do not participate in feature screening and are forced to remain in the model. See SMLE for details.
#' @param num_clusters The number of compute clusters to use when 
#' \code{parallel = TRUE}. The default will be 2 times cores detected.
#'
#' @export
#'
smle_select.sdata<-function(object, k_min=1, k_max=NULL, subset=NULL,
                            gamma_ebic=0.5, vote= FALSE, keyset=NULL,
                            criterion="ebic", codingtype = c("DV","standard","all"),
                            gamma_seq=c(seq(0,1,0.2)), vote_threshold=0.6,
                            parallel = FALSE, num_clusters=NULL,...){
  cl<-match.call()
  
  codingtype <- match.arg( codingtype )
  
  cl[[1]] <- as.name("smle_select")
  
  #----------------------------------------------------------------------------#
  #------------------------- Argument pre-processing --------------------------#
  #----------------------------------------------------------------------------#
  
  X<-object$X

  Y<-object$Y

  family<-object$family

  n<-dim(X)[1]

  pp<- dim(X)[2]
  
  # Check input data dimension, warning for High dimensional input.
  if(is.null(subset)){

    if(pp>999){
      
      message("For high dimensional data input, we recommand to do screening before selection. ")
      Continue <- readline(prompt="Do you want to proceed? Y or N \n ")
      
      while(T){
        
        if(Continue == "N"){
          
          tryCatch(stop(), error = function(e) message('Algorithm aborted'))
          
        }else if(Continue == "Y"){
          
          message('Algorithm processing')
          
          break
          
        }else{
          
          Continue <- readline(prompt='Please input Y or N. \n ')
        }
        
        
      }
      
      
    }
    
    X_s<-X

    }else{

      X_s<-X[,subset]

    }
  
  if(!is.null(keyset)){
    
    keyset <- which(subset %in% keyset)
    
    k_min <- max(k_min, length(keyset)+1)
    
    }
  
  
  categorical <- any(sapply(X_s,is.factor))
  
  
  if(is.null(k_max)){
    
    k_max = dim(X_s)[2]
    
  }

  if(parallel && is.null(num_clusters)){ 
    
    num_cores    <- parallel::detectCores() 
    
    num_clusters <- 2 * num_cores
    
    }

  IP <- NULL; ID_Voted <- NULL; vs <- NULL

  if(k_min<0 || k_max < k_min ||k_min%%1!=0 ||k_max%%1 !=0||k_max>pp){

    tryCatch(stop(), error = function(e) message('Retained model size setting error(k_min and k_max).'))
  }
  
  #----------------------------------------------------------------------------#
  #------------------------------ Selection -----------------------------------#
  #----------------------------------------------------------------------------#

  # Run selection algorithm with information criterion specified (EBIC default)
    
  criter_value<- ebicc(Y,X_s,family,criterion,codingtype,keyset,
                            k_min,k_max,n,pp,gamma_ebic,parallel,num_clusters)
    
  v_s <- which.min(criter_value)+k_min-1
  
  f_s<- SMLE(Y=Y, X=X_s, k=v_s, family=family,keyset=keyset)


    #------------------------------------------------------------------------#
    #------------------------------- Voting ---------------------------------#
    #------------------------------------------------------------------------#

  if( vote ){

    stopifnot("Voting available when criterion = ebic"=criterion == 'ebic')

    vs<-c()
    
    ebic_selection_by_gamma<-function(gamma){
      
      ebic_value<-ebicc(Y,X_s,family,criterion,codingtype,keyset,
                        k_min,k_max,n,pp,gamma,parallel,num_clusters)
      
      v_s <- which.min(ebic_value)+k_min-1
      
      ID <- SMLE(Y=Y, X=X_s, k=v_s, family=family,keyset=keyset)$ID_retained
      
      ID
      
    }
    
    ## Run selection algorithm with different gammas 

    if(parallel ){
      
      if(.Platform$OS.type=="windows"){
        
        cl<-parallel::makeCluster(num_clusters)
        
      }
      
        vs<-unlist(mclapply(gamma_seq,ebic_selection_by_gamma,mc.cores = 2*num_clusters))
        
        }else{
          
      vs<-unlist(lapply(gamma_seq,ebic_selection_by_gamma))

        }
    
    #----------------------------------------------------------------#
    
    IP<-as.factor(vs)
    
    # Count the frequency of features selected with different model of sparsity.
    
    IP_f<-summary(IP)[order(summary(IP),decreasing= T)]/max(summary(IP))
    
    # features with frequency higher than vote_threshold will be selected in ID_voted.
    
    ID_Voted<-as.numeric(names(IP_f[IP_f>=vote_threshold]))
   

    }
  if(is.null(subset)){

    ID_Selected=f_s$ID_retained

  }else{

    ID_Selected<-subset[f_s$ID_retained]

    ID_Voted<-subset[ID_Voted]

  }

  Out<-list(
    
          call = cl, ID_selected=ID_Selected, family = family,

          coef_selected=f_s$coef_retained,
          
          intercept = f_s$intercept,
          
          num_selected = length(f_s$coef_retained),

          vote = vote,criterion = criterion,
          
          ID_pool= IP, keyset = keyset, k_min = k_min, k_max= k_max,
          
          subset = subset,

          criterion_value=criter_value,

          ID_voted = ID_Voted,
          
          categorical = categorical,
          
          vote_threshold=vote_threshold,

          gamma_ebic=gamma_ebic,
          
          codingtype = codingtype,

          gamma_seq=gamma_seq, X = X , Y = Y)
  
  class(Out)<-"selection"
  
  return(Out)
  }

#' @rdname smle_select
#'
#' @method smle_select default
#' 
#' @param Y Input response vector (when \code{object = NULL}).
#'
#' @param X Input features matrix (when \code{object = NULL}).
#'
#' @param family Model assumption; see \code{\link{SMLE}()} documentation. Default is Gaussian linear.
#'
#' When input is a \code{'smle'} or \code{'sdata'} object, the same
#' model will be used in the selection.
#'
#' @param ... Further arguments passed to or from other methods.
#'
#' @export
#'
smle_select.default<-function(object = NULL,Y =NULL, X=NULL, family='gaussian',keyset = NULL,...){
  
  # If the data is passed, this function create the sdata object, and then call 
  # smle_select with the sdata function call.
  
  cl<-match.call()
  
  cl[[1]] <- as.name("smle_select")
  
  if ( any(sapply(X,is.factor) ) ){
  
     categorical = TRUE
  
     }else{
   
       categorical =FALSE
       
  }
  if(is.null(object)){
    
    object<-Y
  
  }
  
  if(!is.null(keyset)){
    
    if( inherits( keyset , c('integer','numeric'))  ){
      
      keyset_ind <- as.character(colnames(X)[keyset])
      
    }
    
    keyset_ind <- which(colnames(X) %in% keyset)
    
  }else{keyset_ind <- keyset}
    
  Data<-structure(list(Y=object,X=X,family=family,categorical = categorical),class = "sdata")

  Out<-smle_select(Data, keyset = keyset_ind,...)
  
  Out$call <- cl
  
  return(Out)
}
#' @rdname  smle_select
#'
#' @method smle_select smle
#'
#' @export
#'

smle_select.smle<-function(object,...){
  
  # When a smle object is called, this function create the sdata object and set 
  # the retained features ID in smle object as the argument subset in smle_select,
  # and then call smle_select with the sdata function call.
  
  cl<-match.call()
  
  cl[[1]] <- as.name("smle_select")
  
  Data<-structure(list( Y=object$Y, X=object$X, family=object$family),class = "sdata")
  
  if(!is.null(object$data)){
    
    data<-object$data
    
    X <- object$X
    
    keyset <- object$keyset
    
    if( inherits( keyset , c('integer','numeric'))  ){
      
      keyset <- as.character(names(data)[keyset])
      
    }
    
    index <- which(colnames(X) %in% object$iteration_data$feature_name)
    
    keyset_ind <- which(colnames(X) %in% keyset)
    
    Out<-smle_select(Data,subset =index, keyset = keyset_ind,...)
    
    feature_name <- colnames(X)[Out$ID_selected]
    
    Out$ID_selected <-  sort(which(names(data) %in% feature_name))
    
    if(Out$vote == TRUE){
        
      vote_name<-colnames(X)[Out$ID_voted]
      
      Out$ID_voted <-sort(which(names(data) %in% vote_name))

    }
    
    Out$data<- data
    
    Out$keyset <- keyset
    
  }else{ 
    
    Out<-smle_select(Data,subset =object$ID_retained,keyset = object$keyset,...)
    
    }
  
  Out$call <- cl
  
  Out$keyset <- object$keyset
  return(Out)
}





