# dummy function are inspired from R dummy package.

# Michel Ballings and Dirk Van den Poel (2015). dummy: Automatic Creation of Dummies with
# Support for Predictive Modeling. R package version 0.1.3.
# https://CRAN.R-project.org/package=dummy


dummy <- function( x, data=NULL, sep="_",
                   
                   drop=TRUE, fun=as.integer, verbose = FALSE,
                   
                   codingtype=c("standard","all","DV") ) {
  
  
  if( length(x) > 1 ) stop( "More than one variable provided to produce dummy variable." )
  
  name <- x
  
  x    <- data[ , name]
  
  #   model.matrix does not work on factor w/ one level. Here we trap for the spacial case.
  
  if( length(levels(x))<2 ) {

    if( verbose ) warning( name, " has only 1 level. Dummy variable removed" )

    return(NULL)

  }

  # Get the model matrix
  
  if(codingtype=="all"){

  mm <- model.matrix( ~ x  - 1 , model.frame( ~ x - 1 ) )

  colnames.mm <- colnames(mm)

  mm <- matrix( mm, nrow=nrow(mm), ncol=ncol(mm), dimnames=list(NULL, colnames.mm) )

  # Replace the column names 'x'... with the true variable name and a seperator
  
  colnames(mm) <- sub( "^x", paste( name, sep, sep="_" ), colnames(mm) )
  
  if(! is.null(row.names(data)) ) rownames(mm) <- rownames(data)

  return(mm)


  }

  else if(codingtype=="standard"){

    mm <- model.matrix( ~ x , model.frame( ~ x ) )

    MM <- mm[,-1]

  }else if(codingtype=="DV"){

    mm <- model.matrix( ~ x - 1 , model.frame( ~ x - 1) )

    levels<- ncol(mm)

    MM<-mm[,-1]- (1/levels)

  }else{
      
    stop(" Codingtype not found")
    
    }

  colnames.MM <- colnames(mm)[-1]

  MM <- matrix( MM, nrow=nrow(mm), dimnames=list(NULL, colnames.MM) )

  # Replace the column names 'x'... with the true variable name and a seperator
  
  colnames(MM) <- sub( "^x", paste( name, sep, sep="_" ), colnames(mm)[-1] )
  
  if(! is.null(row.names(data)) ) rownames(MM) <- rownames(data)

  return(MM)

}


dummy.data.frame <- function( data, codingtype = "standard", 
                              
                              dummy.classes=getOption("dummy.classes"),  ... ) {


  #--------------------- Initialize the data.frame ----------------------------#
  
  df <- data.frame( row.names = row.names( data ))
  
  new.attr <- list()  ### Track location of dummy variables
  
  
  if( is.null( getOption( "dummy.classes" ))){ 
    
    ### Only factor and character columns are treated as categorical features.
    
    options( "dummy.classes" = c("factor","character"))
    
  } 
  
  for( nm in names(data) ) {

        old.attr <- attr(df,'dummy variables')
        
    if( any(dummy.classes == "ALL") || class(data[,nm]) %in% dummy.classes ){
      
      dummy_variables <- dummy( nm, data, codingtype= codingtype, ... )

      if( ncol(dummy_variables)>0 ) { ## Assign dummy variables indices to attribute.
        
        new.attr[[nm]] <- (ncol(df)+1):( ncol(df) + ncol(dummy_variables) )
        
        }

    } else {
      
      dummy_variables <- data[,nm, drop=FALSE ]
      
    }

    df <- as.matrix(cbind(df, dummy_variables))

  }
  
  ## Dummy indices
  
  attr( df, 'dummy variables' ) <- new.attr
  
  ## Assign the starting indice of all dummy variables to DI.
  
  attr( df, 'DI' ) <- unlist(lapply(new.attr,function(l) l[[1]]))
  
  
  
  return(df)
}

