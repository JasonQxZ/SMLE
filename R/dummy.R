dummy <- function( x, data=NULL, sep="_", drop=TRUE, fun=as.integer, verbose = FALSE,codingtype=c("standard","all","DV") ) {


  # HANDLE IF DATA IS MISSING.
  if( is.null(data) ) {


    name <- as.character( sys.call(1) )[2]
    name <- sub( "^(.*\\$)", "", name )    # REMOVE prefix e.f
    name <- sub( "\\[.*\\]$", "", name )   # REMOVE suffix


  } else {


    if( length(x) > 1 ) stop( "More than one variable provided to produce dummy variable." )
    name <- x
    x    <- data[ , name]
  }


  # CHANGE TO FACTOR: KEEP LEVELS?
  if( drop == FALSE && class(x) == "factor" ) {
    x <- factor( x, levels=levels(x), exclude=NULL )
  } else {
    x<-factor( x, exclude=NULL )
  }


  # TRAP FOR ONE LEVEL :
  #   model.matrix does not work on factor w/ one level.  Here we trap for the spacial case.
  if( length(levels(x))<2 ) {

    if( verbose ) warning( name, " has only 1 level. Producing dummy variable anyway." )

    return(
      matrix(
        rep(1,length(x)),
        ncol=1,
        dimnames=list( rownames(x), c( paste( name, sep, x[[1]], sep="" ) ) )
      )
    )

  }

  # GET THE MODEL MATRIX
  if(codingtype=="all"){

  mm <- model.matrix( ~ x  - 1 , model.frame( ~ x - 1 ) )

  colnames.mm <- colnames(mm)

  message( " ", name, ":", ncol(mm), "dummy varibles created\n" )

  mm <- matrix( mm, nrow=nrow(mm), ncol=ncol(mm), dimnames=list(NULL, colnames.mm) )

  # Replace the column names 'x'... with the true variable name and a seperator
  colnames(mm) <- sub( "^x", paste( name, sep, sep="_" ), colnames(mm) )
  if(! is.null(row.names(data)) ) rownames(mm) <- rownames(data)

  return(mm)


  }

  else if(codingtype=="standard"){

    mm <- model.matrix( ~ x , model.frame( ~ x) )

    MM <- mm[,-1]

  }else{

    mm <- model.matrix( ~ x - 1 , model.frame( ~ x -1) )

    levels<- ncol(mm)

    MM<-mm[,-1]- (1/levels)

    }

  colnames.MM <- colnames(mm)[-1]

  MM <- matrix( MM, nrow=nrow(mm), dimnames=list(NULL, colnames.MM) )

  # Replace the column names 'x'... with the true variable name and a seperator
  colnames(MM) <- sub( "^x", paste( name, sep, sep="_" ), colnames(mm)[-1] )
  if(! is.null(row.names(data)) ) rownames(MM) <- rownames(data)

  return(MM)

}


dummy.data.frame <- function( data,codingtype, 
                              dummy.classes=getOption("dummy.classes"), all=TRUE, ... ) {


  # Initialize the data.frame
  df<-data.frame( row.names=row.names(data) )
  new.attr <- list()  # Track location of dummy variables
  if( is.null( getOption("dummy.classes") ) ) options( "dummy.classes" = c("factor","character") )
  for( nm in names(data) ) {

    # cat( nm )
    old.attr <- attr(df,'dummies')

    if((class(data[,nm])[1] %in% dummy.classes ))
     {
      dummies <- dummy( nm, data, codingtype= codingtype )

      # OMIT CONSTANT COLUMNS:
      #  Variables that are constant will return a matrix with one column

      if( ncol(dummies)>0 ) new.attr[[nm]] <- (ncol(df)+1):( ncol(df)+ncol(dummies) )

    } else {
      if( ! all ) next()
      dummies <- data[,nm, drop=FALSE ]
    }

    df <- cbind(df, dummies)

  }

  attr( df, 'dummies' ) <- new.attr
  return(df)
}
