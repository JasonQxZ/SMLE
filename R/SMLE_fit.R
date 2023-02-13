SMLE_fit<-function(Y,X, k, family, keyset=NULL,categorical= FALSE,CI=NULL,
                   intercept, max_iter=500, tol=10^(-3), fast=FALSE, group= group,
                   U = 1, U_rate=0.5, X_mean=NULL,X_sd=NULL,coef_initial=NULL, 
                   penalize_mod, codingtype = NULL, Data_X, standardize){
  #----------------------------------------------------------------------------#
  #----------------------------Input preprocess--------------------------------#
  #----------------------------------------------------------------------------#
  
  LH <- rep(0,max_iter)     # Likelihood
  
  FD <- NULL  # Number of differences in retained features between adjacent steps
  
  number_of_Ucheck <- rep(0,max_iter) # Vector of number of Ucheck steps per iteration.
  
  n<-dim(X)[1] # Number of observations
  
  p<-dim(X)[2] # Number of features
  
  if( !group & any(CI %in% keyset)){
    
      stop("The variable 'group' must be set to TRUE if keyset contains categorical features.")
    
  }
  if(is.null(keyset)){
    
    number_of_ID_retained <-k
    
  }else{
    
    number_of_ID_retained<-k-(length(keyset)) 
    
    if(number_of_ID_retained <= 0){
      
      stop("The screening size k should be larger than the number of features in keyset.")
      
    }
    
  }

  #----------------------------------------------------------------------------#
  #-------------------------Algorithm Initialization --------------------------#
  #----------------------------------------------------------------------------#

  if(categorical){
    
    ## Getting categorical dummy variable and their indices.
    
    X_iter<-dummy.data.frame(X,codingtype = codingtype)
    
    ## Vector tracking dummy variables indices in working matrix X_iter.
    
    DFI <- attr(X_iter,"dummy variables")
    
    DI <- attr(X_iter,"DI")
    
    if( is.null( coef_initial ) ){ ## Default start
      
      fit_pre<-glmnet(x=X_iter,y=Y,family=family)
      
      Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
      
    }else{ ## Start from user supplied initial coefficient values.
      
      Beta0<- rep(0,dim(X_iter)[2])
      
      Beta0[-unlist(DFI)] <- coef_initial[-CI]
      
      for(i in 1:length(CI)){
        
        # User supplied initial coefficients to give all categorical dummy variables the supplied value
        
        Beta0[DFI[[i]]]<-rep(coef_initial[CI[i]],length(DFI[[i]]))
        
      }
      
    }
    
    ## Intercept in the categorical case if codingtype is not 'all'.
    if(codingtype == 'all'){ 
      
      intercept = FALSE 
      
      }else{
        
        intercept = TRUE 
        
        Beta0 <- c(mean(Y-X_iter%*%Beta0),Beta0)
        
        X_iter<-cbind(matrix(1,nrow  = n, ncol = 1),X_iter)
        
        colnames(X_iter)[1]<-'(intercept)'
        
        names(Beta0)[1]="(intercept)"
        
        }
    
    pp<-dim(X_iter)[2] # Number of features in the feature matrix with dummy variables (X_iter)
    
    ID_None0<-(1:pp)[Beta0!=0] #Initialization of indices of non-zero features
    
    I<-list(Y=Y,CM=X,CI=sort(CI),dum_col=lapply(DFI,function(x) length(x)),IM=X_iter,
            DFI= lapply(DFI,function(x) x+1*(intercept)),DI= DI+1*(intercept),
            family=family,codingtype=codingtype,intercept = intercept)
    
    beta_path<-as.matrix(Group_Beta(Beta0,I,penalize_mod),nrow=pp,ncol=1)
    
    Screening_Dindex<-sub_off(1:pp,CI2DI(I,keyset))
    
  }else{# No categorical features; only numerical features
    
    # numerical feature matrix pre-processing
    
    if( is.null( coef_initial ) ){ 
      
      fit_pre <- glmnet(x=X,y=Y,family=family$family)
      
      Beta0   <- c(fit_pre$beta[,dim(fit_pre$beta)[2]])
      
    }else{ 
      
      Beta0 <- coef_initial 
      
    }
    
    if(intercept){
      
      Beta0<-c(mean(Y-X%*%Beta0),Beta0)
      
      names(Beta0)[1] = "(intercept)"
      
      X_iter<-cbind(matrix(1,nrow  = n, ncol = 1),X)
      
      ID_None0<-(1:(p+1))[Beta0!=0] # Initialization of indices of non-zero features
      
      keyset<-c(1,keyset+1)
      
    }else{
      
      X_iter<-X
      
      ID_None0<-(1:p)[Beta0!=0]
      
    }
    
    beta_path<-as.matrix(Beta0 , nrow=pp , ncol=1)
    
    pp <- p+1*(intercept==TRUE)
    
    Screening_index <- sub_off(1:pp,keyset)
    
    I<-list(CM=X,Y=Y,IM=X_iter,beta0=Beta0,family=family)
  }
  
  ID_None0 <- which(Beta0 != 0) 
  
  coef_None0 <- as.matrix( Beta0[ID_None0] , ncol=1) #Initialization of retained coefficients
  
  Xs_0 <- as.matrix(X_iter[, ID_None0])
  
  if(categorical){
    
    if(group==T){
      
      Beta0 <- GroupHard(Beta0,I,k = k-(length(keyset)),DI2CI(I,Screening_Dindex),penalize_mod)
      # length(Beta_t) = pp
    }else{
      
      Beta0[Screening_Dindex] <- Hard(t=Beta0[Screening_Dindex], k=number_of_ID_retained)
      # length(Beta_t) = pp
    }
    
  }else{
    
    Beta0[Screening_index] <- Hard(t=Beta0[Screening_index],k= number_of_ID_retained)
    
  }


  if(!is.null(coef_initial) & sum(coef_initial)==0 ){ 
    
    theta_0 <- matrix(0, ncol=1, nrow=n) 
    
    ## When starts from zero with u_s is 1/||X||\infit = 1/sqrt(n) 
    
    U_0 <-U/(sqrt(p))
    
    }else{
      theta_0  <- Xs_0 %*% coef_None0

      U_0 <- U/(sqrt(p)*norm(Xs_0, "i")^2)
      
    }
  
  theta_0 <- make.link(family$link)$linkinv(theta_0)
  
  V_0 <- crossprod(X_iter , Y - theta_0)
  #iteration start form 1
  i <- 1
  
  Beta_s<-Beta0
  
  V_s <- V_0
  
  u_s <- U_0

  #----------------------------------------------------------------------------#
  #------------------------- Algorithm Iteration ------------------------------#
  #----------------------------------------------------------------------------#
  
  repeat{
    
    # Main iteration
    
    count<-0
    
    repeat{ 
      
      # Iteration for ucheck
      
      Beta_t<- Beta_s + u_s * V_s
      
      # Hard threshold
      
      if(categorical){
        
        if(group==T){
          
          Beta_t <- GroupHard(Beta_t,I,k = k-(length(keyset)),DI2CI(I,Screening_Dindex),penalize_mod)
          # length(Beta_t) = pp
        }else{
          
          Beta_t[Screening_Dindex] <- Hard(t=Beta_t[Screening_Dindex], k=number_of_ID_retained)
          # length(Beta_t) = pp
        }
        
        
      }else{
        
        Beta_t[Screening_index] <- Hard(t=Beta_t[Screening_index],k= number_of_ID_retained)
        
      }

      ucheck<- Ucheck(Y = Y, X=X_iter, beta1=Beta_s, beta2=Beta_t, family=family)
      
      if (ucheck >= 0){
        
        # Break the U check when likelihood decreasing.
        
        break
        
      }else{
        
        u_s <- U_rate * u_s
        
        count<-count+1
        
      }
      
    }
    
    #  Retained features at step s
    
    sindex <- which(Beta_s!= 0)
    
    #  Retained features at step t
    
    tindex <- which(Beta_t!= 0)
    
    # Difference of retained features between two steps
    
    fs<-sum(!(tindex %in% sindex))
    
    FD<-c(FD,fs)
    
    if(categorical){
      
      beta_path<-cbind(beta_path,as.matrix(Group_Beta(Beta_t,I,penalize_mod),
                                           ncol=1,nrow=pp))
      
    }else{
      
      beta_path<-cbind(beta_path,as.matrix(Beta_t,ncol=1,nrow=pp))
      
    }
    
    LH[i]<-lh(Y=Y, X=X_iter, beta= Beta_t,family=family)
    
    valid_LH_diff<- 0.01 *(LH[2] - LH[1])
    
    number_of_Ucheck[i]<-count
    
    
    #------------------------ convergence check -------------------------------#
    if(i>1){
      
      MSE<- sqrt(sum(( Beta_s-Beta_t )^2))
      
      if(fast){
        
        ## Break the iteration When stopping criteria satisfied.
        
        if( MSE/number_of_ID_retained < tol){break}
        
        else if( (LH[i]-LH[i-1])< valid_LH_diff){break}
        
        else if(i>10){if(sum(tail(FD,10))==0){break}}
        
      }else{
        
        if( MSE < tol || i >= max_iter){
          break}
      }
      
      
    }
    ## Update all information if the iteration continued.
    
    Beta_s<-Beta_t
    
    ID_None_s <- sort((1:pp)[Beta_s!= 0])
    
    coef_None_s <- as.matrix(Beta_s[ID_None_s],ncol=1)
    
    Xs_s <- X_iter[, ID_None_s]
    
    theta_s  <- make.link(family$link)$linkinv(Xs_s %*% coef_None_s)
    
    V_s <- crossprod(X_iter, Y - theta_s)
    
    u_s  <- U_0
    
    i<-i+1
    
  }
  
  #----------------------------------------------------------------------------#
  #----------------------------- Post Iteration -------------------------------#
  #----------------------------------------------------------------------------#
  
  
  Coef_dist<-apply((beta_path[,-1]-beta_path[,-ncol(beta_path)]),2,function(x){sqrt(sum(x^2))})
  
  # Degree of freedom
  
  df <- length(coef_None_s)
  
  # Re-scale all output
  
  if(categorical){
    
    if(codingtype != 'all'){
      
      Intercept_value<-coef_None_s[1]
      
      ID_None_s <- DI2CI(I,ID_None_s)   # length pp -> p
      
      if(group ){ coef <- Group_Beta(Beta_s,I,penalize_mod)
      
      }else{ coef <- Group_Beta(Beta_s,I,max=T) }
      
      coef_None_s <- coef[coef!=0][-1]
      
      beta_path <- beta_path[-1,]
      
    }else{
      
      Intercept_value  <- NULL
      
      ID_None_s <- DI2CI(I,ID_None_s)  # length pp -> p
      
      if(group){coef <- Group_Beta(Beta_s,I,penalize_mod)
      
      }else{ coef <- Group_Beta(Beta_s,I,max=T) }
      
      coef_None_s <- coef[coef!=0]
    
    }
    
    # Check if retained model has numerical features
    
    Numeric_ind <- ID_None_s[- which(ID_None_s %in% CI)]
    
    if(!identical(Numeric_ind,numeric(0))){
      
      # Re-scale numerical coefficients
      
      coef_None_s[ID_None_s %in% Numeric_ind] <- coef_None_s[ID_None_s %in% Numeric_ind] * X_sd[Numeric_ind] + X_mean[Numeric_ind]
      
    } 
    # else output coef_None_s directly
    
    # re-scale beta_path
    
    beta_path[-CI,]<- beta_path[-CI,] * X_sd + X_mean
    
  }else{
    
    if( intercept ){
      
      # re-sale beta_path
      
      beta_path<-beta_path[-1,]           
      
      Intercept_value<- as.vector(coef_None_s)[1]
      
      coef_None_s <- coef_None_s[-1]
      
      ID_None_s<-ID_None_s[-1]-1
      
      if(!is.null(X_mean)){
        
        coef_None_s <- coef_None_s*X_sd[ID_None_s] + X_mean[ID_None_s]
        
        beta_path- beta_path*X_sd+X_mean
        
      }
        
    }else{
      
      if(!is.null(X_mean)){

      coef_None_s<- coef_None_s*X_sd[ID_None_s] + X_mean[ID_None_s]
        
      beta_path<- beta_path*X_sd+X_mean
      
      }
      
      Intercept_value <- NULL
      
    }
  }
  
  feature_name <- NULL
  
  if( !is.null(colnames(Data_X))){ 
    
    feature_name <- colnames(Data_X)[ID_None_s]
    
  }else{
    
    feature_name <- as.character(ID_None_s)}
  
  
  fit<-list(iteration_data = list(IM = I$IM, beta0 = Beta0, FD = FD,            
                                  
                                  feature_name =feature_name ),
            
            X=Data_X, Y=Y, keyset = keyset, family = family, k = k,
            
            intercept = Intercept_value, steps = i, df = df,
            
            likelihood_iter = LH[1:i],
            
            Usearch = number_of_Ucheck[1:i],
            
            path_retained  = beta_path,
            
            num_retained   = length(ID_None_s),
            
            ID_retained    = ID_None_s,
            
            coef_retained  = coef_None_s,
            
            categorical = categorical, coef_dist= Coef_dist,
            
            fast=fast
  )
  
  if(categorical){
    
    fit$modified_data = list(CM = I$CM, CI = sort(I$CI), dum_col = I$dum_col, DFI = I$DFI)

    
  }
  
  class(fit)="smle"
  
  fit
}
