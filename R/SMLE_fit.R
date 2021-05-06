SMLE_fit<-function(Y , X, k, family="gaussian", keyset=NULL,
                 intercept=TRUE, max_iter=500, tol=0.01, fast=FALSE,
                 U_rate=0.5,X_mean=NULL,X_sd=NULL,Coef_initial=NULL){
  
  LH<-rep(0,max_iter)
  
  number_of_Ucheck<-rep(0,max_iter)

  n<-dim(X)[1];p<-dim(X)[2]
  
  if( is.null( Coef_initial ) ){
    
    fit_pre<-glmnet(x=X,y=Y,family=family)
    
    Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
    
  }else{ 
    
    Beta0 <- Coef_initial }

  FD<-NULL

  if(is.null(keyset)){
    
    number_of_ID_retained <-k
    
    }
  else{
    
    number_of_ID_retained<-k-(length(keyset)) 
    
  }
    
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
    
    I<-list(CM=X,Y=Y,IM=X_iter,Beta0=Beta0,family=family)
    
    pp<-p+intercept
    
    coef_None0 <- as.matrix(Beta0[ID_None0],ncol=1)
    
    Xs_0 <- X_iter[, ID_None0]
    
    if(!is.null(Coef_initial) & sum(Coef_initial!=0)==0){ 
      
      R_0<-matrix(0, ncol=1, nrow=n) 
      
    }else{
      R_0  <- Xs_0 %*% coef_None0
    }

    
    R_0<-switch(family,
                
                "gaussian"=R_0,
                
                "poisson"=exp(R_0),
                
                'binomial'=exp(R_0)/(1+exp(R_0)))
    
    V_0<-crossprod(X_iter,Y - R_0)
    
  
    if(!is.null(Coef_initial) & sum(Coef_initial!=0)==0){
      
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

        if ( ucheck >= 1 ){

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

      ID_None0<- (1:pp)[Beta_s!= 0]

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
    
    #rescale beta_path
    
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
  
  
  

  fit<-list(I=I,X=I$CM, Y=I$Y,

            keyset = keyset, family = family, k = k,

            Intercept=Intercept_value,

            steps = i,

            LH=LH[1:i],

            Usearch=number_of_Ucheck[1:i],

            Path_Retained  = beta_path,

            Num_Retained   = length(ID_None0),

            ID_Retained    = ID_None0,

            Coef_Retained  = coef_None0,

            FD=FD, Ctg = FALSE, Coef_dist=Coef_dist,

            fast=fast
            )
}

