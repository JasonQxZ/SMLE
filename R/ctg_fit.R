ctg_fit<-function(Y , X , k ,
                  
                  family = c("gaussian","binomial","poisson"),
                  
                  categorical  , keyset ,CI, 
                  
                  max_iter, tol ,fast,
                  
                  intercept , group ,
                  
                  codingtype , penalize_mod,
                  
                  U_rate , X_mean , X_sd, Coef_initial,cl,Data_X){
  
  family<-match.arg(family)
  
  n<-dim(X)[1];p<-dim(X)[2]
  
  #--------------------------------------------------------------#
  
  
  if(codingtype=="all"|| group==FALSE){
    
    dum_col<-sapply(X[,CI],nlevels)
    
  }else{
    
    dum_col<-sapply(X[,CI],nlevels)-1
    
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
  
  
  if( is.null( Coef_initial ) ){
    
    fit_pre<-glmnet(x=as.matrix(X_dummy,dimnames = dimnames(X_dummy)),y=Y,family=family)
    
    Beta0<-c(fit_pre$beta[,dim(fit_pre$beta)[2]])
    
  }else{ 
    
    Beta0<- rep(0,dim(X_dummy)[2])
    
    Beta0[-unlist(DFI)] <- Coef_initial[-CI]
    
    for(i in 1:length(CI)){
      
      Beta0[Dummy_index[[i]]]<-rep(Coef_initial[CI[i]],dum_col[i])
      
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
  
  ID_None0<-(1:pp)[Beta0!=0]  #length : = pp
  
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
              'binomial'=exp(R_0)/(1+exp(R_0))
  )
  
  
  V_0<-crossprod(X_iter,Y - R_0)
  
  if(!is.null(Coef_initial) & sum(Coef_initial!=0)==0){
    
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
  
  Screening_Dindex<-sub_off(1:pp,CI2DI(I,keyset))
  
  if(is.null(keyset)){number_of_ID_retained <-k
  }else{number_of_ID_retained<-k-(length(keyset)) }
  
  repeat{
    
    count<-0
    
    repeat
    {
      
      
      Beta_t<-Beta_s + uu * V_0         # length(Beta_s)  = pp
      
      if(group==T){
        
        Beta_t<-GroupHard(Beta_t,I,k=k-(length(keyset)),Screening_index,penalize_mod)
        
        
        # length(Beta_t) = pp
        
      }else{
        
        Beta_t[Screening_Dindex]<- Hard(t=Beta_t[Screening_Dindex], k=number_of_ID_retained)
        
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
    
    ID_None0<- (1:pp)[Beta_s!= 0]
    
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
  
  if(group == TRUE){
    
    ID_None0<-DI2CI(I,ID_None0)   # issue here
    
    coef <- Group_Beta(Beta_s,I,penalize_mod)
    
    coef_None0 <- coef[coef!=0][-1]
    
  }else{
    
    ID_None0 <- DI2CI(I,ID_None0)
    
    coef_None0<-coef_None0[-1]
    
  }
  
  
  fit<-list(I=I,X=Data_X,Y=I$Y,
            ID_retained=ID_None0,
            coef_retained= coef_None0,
            keyset=keyset,
            Usearch=number_of_Ucheck[1:i],
            family=family,
            k=k,
            intercept=Intercept_value,
            steps = i,
            likelihood_iter=LH[1:i],
            path_retained=beta_path[-1,],
            num_retained   = length(ID_None0),
            group=group,
            fast =FALSE,
            FD= FD, 
            ctg =TRUE,coef_dist=Coef_dist,
            codingtype = codingtype,
            fast=fast
  )
  
  fit$call=cl
  class(fit)="smle"
  fit
  
}