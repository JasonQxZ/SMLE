
Ucheck <-function(Y,X, beta1, beta2, family){
  
  #----------------------------------------------------------------------------#
  # This is the ucheck function                                                #
  #----------------------------------------------------------------------------#
  # This function compare the likelihood between beta1 and beta2 (the          #
  # coefficient at the current step and the coefficient for the next step.     #
  #----------------------------------------------------------------------------#
  # Parameters:  
  # Y : The response vector.
  # X : Feature matrix.      
  # beta1 : Coefficient of the current step. 
  # beta2 : Coefficient of the next step. 
  # family: Model assumption.                                                  #
  #----------------------------------------------------------------------------#
  
  ll_1 <- lh(Y=Y,X=X,beta=beta1,family=family)
  
  ll_2 <- lh(Y=Y,X=X,beta=beta2,family=family)
  
  return(ll_2 - ll_1)
  
}

lh<-function(Y, X, beta, family){
  
  #----------------------------------------------------------------------------#
  # This function calculate the log-likelihood                                 #
  #----------------------------------------------------------------------------#
  # This function call the aic() from family object to calculate the likelihood#
  # given the response Y, the data matrix X the coefficient beta and the model #
  # assumption family.                                                         #
  #----------------------------------------------------------------------------#
  # Parameters:  
  # Y : The response vector.
  # X : Feature matrix.      
  # beta : Coefficient input.
  # family: Model assumption.                                                  #
  #----------------------------------------------------------------------------#
  
  #Calculating log likelihood from family object in stats package
  
  linkinv <- family$linkinv
  aic <- family$aic
  dev.resids <- family$dev.resids
  
  y <- Y
  n = length(y)
  ind <- which(beta != 0)
  eta <- X[,ind]%*%as.matrix(beta[ind],nrow = length(ind))
  mu <- linkinv(eta)
  dev <- sum(dev.resids(y, mu, 1))
  return((-1/2)*aic(y, n, mu, 1, dev))
}


CI2DI<-function(I,CI_index){
  
  #----------------------------------------------------------------------------#
  # CI convert to DI                                                           #
  #----------------------------------------------------------------------------#
  # This function converts the index of categorical features from the original # 
  # matrix to the working dummy matrix.                                        #
  #----------------------------------------------------------------------------#
  # Parameters:                                                                #
  # I : A list contain data information in a smle object.                      #
  # CI_index: Index of categorical features in the original input matrix.      #
  #----------------------------------------------------------------------------#
  
  DI_index<-c();j=1;L=0

  for(i in 1:length(I$CM)){

    if(i %in% CI_index && i %in% I$CI){DI_index<-append(DI_index,i+L+seq(I$dum_col[[j]])-1)}

    if(i %in% CI_index && !(i %in% I$CI)){DI_index<-append(DI_index,i+L)}

    if(i %in% I$CI){L<-L+I$dum_col[[j]]-1;j<-j+1}
  }
  as.vector(DI_index+1*(I$intercept ))
  #If codingtype is not 'all', intercept is TRUE (1).
  
}


DI2CI<-function(I,DI_index){
  
  #----------------------------------------------------------------------------#
  # DI convert to CI                                                           #
  #----------------------------------------------------------------------------#
  # This function converts the index of categorical features from the working  # 
  # dummy matrix to the original matrix.                                       #
  #----------------------------------------------------------------------------#
  # Parameters:                                                                #
  # I : A list contain data information in a smle object.                      #
  # CI_index: Index of categorical features in the original input matrix.      #
  #----------------------------------------------------------------------------#

  CI_index<-c();j=0;L=0
  
  #If codingtype is not 'all', intercept is TRUE (1).
  
  intercept = 1*(I$intercept)

  for(i in (1+intercept):dim(I$IM)[2]){#looping over columns in dummy variable matrix

    if(i %in% DI_index && !(i %in% unlist(I$DFI))){

      #ind not Dummy

      CI_index<-append(CI_index,i-L)
      
      }

    if(i %in% DI_index && i %in% I$DI){

      CI_index<-append(CI_index,i-L)

      }

    if(i %in% DI_index && i %in% unlist(I$DFI) && !(i %in% I$DI)){CI_index<-append(CI_index,I$DI[j]-(L-I$dum_col[[j]]+1))}

    if(i %in% I$DI){L<-L+I$dum_col[[j+1]]-1;j<-j+1}
  }
  CI_index<-unique(CI_index)

  CI_index<-CI_index-intercept

  as.vector(CI_index)
}



Hard<-function(t,k,lam =1)
{
  #----------------------------------------------------------------------------#
  # Hard threshold function                                                    #
  #----------------------------------------------------------------------------#
  # This function keep the k elements with the larges absolute value and       #
  # truncated the rest to zero.                                                #
  #----------------------------------------------------------------------------#
  # Parameters:                                                                #
  # t : A list of values.                                                      #
  # k:  number of values retained.                                             #
  #----------------------------------------------------------------------------#
  
  y<-t
  t<-abs(t)
  
  if(k > 0)
  { lam<-sort(t, decreasing=TRUE)[k] }
  
  y[t<lam]<-0
  return(y)
}



GroupHard<-function(Beta_t,I,k,Screening_index,penalize_mod){
  
  #----------------------------------------------------------------------------#
  # Hard threshold function for dummy variables as a group.                    #
  #----------------------------------------------------------------------------#
  # This function calculate the importance of dummy variable as a group        #
  #----------------------------------------------------------------------------#
  # Parameters:                                                                #
  # I : A list contain data information in a smle object.                      #
  # Beta_t: Current coefficient.                                               #
  # k : Number of features retained.                                           #
  # Screening_index : Index of features can be screened off.                   #
  # penalize_mod: Logical flag indicates whether penalize group features.      #
  #----------------------------------------------------------------------------#

  x<-Group_Beta(Beta_t,I,penalize_mod)[-1]  # length p (intercept removed)

  x[Screening_index]<-Hard(t=x[Screening_index], k=k)

  CSI<-which(x==0) # the position of coefficient zeros in CM.

  DSI<-CI2DI(I,CSI) # the position of coefficient zeros with intercept in DM.

  Beta_t[DSI]<-0

  Beta_t
}


Group_Beta<-function(Beta_t,I,penalize_mod = F,max=FALSE){
  
  #----------------------------------------------------------------------------#
  # Calculate the coefficient of a categorical feature from its dummy columns. #
  #----------------------------------------------------------------------------#
  # Parameters:                                                                #
  # I : A list contain data information in a smle object.                      #
  # Beta_t: Current coefficient.                                               #
  # max : Logical flag indicates whether chose the largest effect of the dummy #
  # variables to as the effect of the categorical feature. If FALSE, take      #
  # square sum of each dummy columns and then take a square root.              #  
  # penalize_mod: Logical flag indicates whether penalize group features.      #
  #----------------------------------------------------------------------------#

  num<-as.vector(Beta_t)
  
  num<-num[-unlist(I$DFI)]

  if(I$intercept){
    
    Intercept_value = num[1]
    
    num<-num[-1]
    
    }

  x<-c()

  ci<-I$CI

  for(i in 1:length(I$CM))
  {
    if(i %in% ci){
      j<-match(i,ci)
      if(penalize_mod){x<-append(x,sqrt(crossprod(Beta_t[I$DFI[[j]]]))/sqrt(I$dum_col[[j]]))}
      else if(max == T){x<-append(x,max(abs(Beta_t[I$DFI[[j]]])))}
      else{x<-append(x,sqrt(crossprod(Beta_t[I$DFI[[j]]])))}
    }else{
      x<-append(x,num[1])
      num<-num[-1]
    }
  }
  names(x)<-names(I$CM)
  
  if(I$intercept){x <- c(Intercept_value,x)}
  
  return(x)
  
}

sub_off<-function(all,sub){
  
  #----------------------------------------------------------------------------#
  # Take elements off a list                                                   #
  #----------------------------------------------------------------------------#
  
  return(all[!(all %in% sub)])
}

