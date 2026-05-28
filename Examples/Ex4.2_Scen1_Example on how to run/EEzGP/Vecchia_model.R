
#'
#' @export
fit_model_HEzGP2 <- function(y, locs, X = NULL, covfun_name = "HEzGP_cov",
                           NNarray = NULL, start_parms = NULL, reorder = TRUE, group = TRUE,
                           m_seq = c(10,30), max_iter = 40, fixed_parms = NULL,
                           silent = FALSE, st_scale = NULL, convtol = 1e-4,p_e,q_e,m_e){
  
  n <- length(y)
  
  # check that length of observation vector same as
  # number of locations
  locs <- as.matrix(locs)
  if( nrow(locs) != n ){
    stop("length of observation vector y not equal
              to the number of locations (rows in locs)")
  }
  
  # check if design matrix is specified
  if( is.null(X) ){
    if(!silent) cat("Design matrix not specified, using constant mean \n")
    X <- rep(1,n)
  }
  X <- as.matrix(X)
  
  if( ! covfun_name %in%
      c("HEzGP_cov")){
    stop("unrecognized covariance function name `covfun_name'.")
  }
  
  
  # detect and remove missing values
  not_missing <- apply( cbind(y,locs,X), 1,
                        function(x){
                          if( sum(is.na(x) | is.infinite(x)) > 0 ){
                            return(FALSE)
                          } else { return(TRUE) }
                        }
  )
  if( sum(not_missing) < n ){
    y <- y[not_missing]
    locs <- locs[not_missing,,drop=FALSE]
    X <- X[not_missing,,drop=FALSE]
    cat(paste0( n - sum(not_missing),
                " observations removed due to missingness or Inf\n"))
  }
  
  # redefine n
  n <- length(y)
  
  # check that start_parms is specified when fixed_parms is
  if( is.null(fixed_parms) ){
    if( is.null(start_parms) ){
      start <- get_start_parms_modif(y,X,locs,covfun_name,p_e,q_e,m_e)
      start_parms <- start$start_parms
    } else {
      # check if start_parms has the right length
      start <-get_start_parms_modif(y,X,locs,covfun_name,p_e,q_e,m_e)
      if(length(start_parms) != length(start$start_parms) ){
        stop(paste0("start_parms not correct length for ",covfun_name))
      }
    }
    # define the parameters we are not fixing
    active <- rep(TRUE, length(start_parms) )
  } else {
    if( is.null(start_parms) ){
      stop("start_parms must be specified whenever fixed_parms is")
    }
    # check if start_parms has the right length
    start <- get_start_parms_modif(y,X,locs,covfun_name,p_e,q_e,m_e)
    if(length(start_parms) != length(start$start_parms) ){
      stop(paste0("start parms not correct length for ",covfun_name))
    }
    # check whether fixed_parms has appropriate values
    if( max( fixed_parms - floor(fixed_parms) ) > 0 ){
      stop("fixed_parms must contain indices of parms you want to fix")
    }
    if( min( fixed_parms < 1 ) || max(fixed_parms) > length(start_parms) ){
      stop("fixed_parms must be between 1 and number of parameters")
    }
    # define the parameters we are not fixing
    active <- rep(TRUE, length(start_parms) )
    active[fixed_parms] <- FALSE
  } 
  
  # get link functions
  linkfuns <- get_linkfun(covfun_name)
  link <- linkfuns$link
  dlink <- linkfuns$dlink
  invlink <- linkfuns$invlink
  invlink_startparms <- invlink(start_parms)
  lonlat <- linkfuns$lonlat
  if(lonlat){
    cat("Assuming columns 1 and 2 of locs are (longitude,latidue) in degrees\n")
  }
  space_time <- linkfuns$space_time
  
  penalty <- get_penalty(y,X,locs,covfun_name) 
  pen <- penalty$pen
  dpen <- penalty$dpen
  ddpen <- penalty$ddpen
  
  # get an ordering and reorder everything
  if(reorder){
    if(!silent) cat("Reordering...")
    if( n < 1e5 ){  # maxmin ordering if n < 100000
      ord <- order_maxmin(locs, lonlat = lonlat, space_time = space_time)
    } else {        # otherwise random order
      ord <- sample(n)
    }
    if(!silent) cat("Done \n")
  } else {
    ord <- 1:n
  }
  yord <- y[ord]
  locsord <- locs[ord,,drop=FALSE]
  Xord <- as.matrix( X[ord,,drop=FALSE] )
  
  # get neighbor array if not provided
  if( is.null(NNarray) ){
    if(!silent) cat("Finding nearest neighbors...")
    
    # need to ignore category variable for categorical covfuns
    if( covfun_name %in% c("matern_categorical","matern_spacetime_categorical","matern_spacetime_categorical_local") ){
      locs_for_NN <- locsord[ , 1:(ncol(locsord)-1) ]
    } else {
      locs_for_NN <- locsord
    }
    
    NNarray <- find_ordered_nn(locs_for_NN, m=max(m_seq), lonlat = lonlat,
                               st_scale = st_scale)
    
    if(!silent) cat("Done \n")
  }
  
  # refine the estimates using m in m_seq
  for(i in 1:length(m_seq)){
    m <- m_seq[i]

    
    likfun <- function(logparms){
      
      lp <- rep(NA,length(start_parms))
      lp<- logparms
      likobj <- vecchia_profbeta_loglik_grad_info_modif(
        link(lp),covfun_name,yord,Xord,locsord,NNarray[,1:(m+1)],p_e=p_e,q_e=q_e,m_e=m_h)
      likobj$loglik <- -likobj$loglik - pen(link(lp))
      likobj$grad <- -c(likobj$grad)*dlink(lp) -
        dpen(link(lp))*dlink(lp)
      likobj$info <- likobj$info*outer(dlink(lp),dlink(lp)) -
        ddpen(link(lp))*outer(dlink(lp),dlink(lp))
      likobj$grad <- likobj$grad
      likobj$info <- likobj$info
      return(likobj)
    }
    fit <- fisher_scoring( likfun,invlink(start_parms)[active],
                           link,silent=silent, convtol = convtol, max_iter = max_iter )
    invlink_startparms[active] <- fit$logparms
    #start_parms[active] <- fit$covparms
    start_parms <- link(invlink_startparms)
    fit$loglik <- -fit$loglik - pen(start_parms)
    invlink_startparms <- invlink(start_parms)
  }
  
  # return fit and information used for predictions
  fit$covfun_name <- covfun_name
  #fit$covparms <- start_parms
  lp <- rep(NA,length(start_parms))
  lp[active] <- fit$logparms
  lp[!active] <- invlink_startparms[!active]
  fit$covparms <- link(lp)
  fit$y <- y
  fit$locs <- locs
  fit$X <- X
  class(fit) <- "GpGp_fit"
  return(fit)
}



#' Print summary of GpGp fit
#'
#' @param object Object of class "GpGp_fit", usually the return value from
#' \code{\link{fit_model}}
#' @param ... additional arguments, for compatability with S3 generic 'summary'
#' @export
summary.GpGp_fit <- function(object, ...){
  cat(paste0("Covariance Function: ",object$covfun_name,"\n\n"))
  cat(paste0("Covariance Parameters: \n"))
  cat(paste0(round(object$covparms,4)),"\n\n")
  cat(paste0("Loglikelihood: ",round(object$loglik,4),"\n\n"))
  X <- as.data.frame(object$X)
  df <- data.frame(
    variable = colnames(X),
    estimate = round(object$betahat,4),
    std_error = round(object$sebeta,4),
    t_stat = round(object$tbeta,4)
  )
  rownames(df) <- c()
  cat("Linear Mean Parameters: \n")
  print(df)
  cat("\n")
  
}

get_start_parms_modif <- function(y,X,locs,covfun_name,p_e,q_e,m_e){
  
  fitlm <- stats::lm(y ~ X - 1 )
  start_var <- summary(fitlm)$sigma^2
  n <- length(y)
  
  randinds <- sample(1:n, min(n,200))
  dmat <- fields::rdist(locs[randinds,])
  
  if(covfun_name ==("HEzGP_cov")){
    npar=1+p_e+q_e+p_e
    
    start_parms <- c(rep(var(y),q_e+1))
    dmat <- fields::rdist(locs[randinds,])
    start_parms <- c(start_parms, rep(stats::median(dmat)/4,2*p_e) )
    #start_parms <- c(start_parms, rep(stats::median(dmat),2*p_e) )
    
     start_range <- median( dmat )/4
   # start_range <- median( dmat )
    start_parms <- c(rep(start_var,q_e+1), rep(start_range, 2*p_e))
  }
  
  return( list( start_parms = start_parms ))
}






#' get penalty function
#'
#' @inheritParams get_start_parms
get_penalty <- function(y,X,locs,covfun_name){
  
  fitlm <- stats::lm(y ~ X - 1 )
  vv <- summary(fitlm)$sigma^2
  # by default, no penalty
  pen <- function(x) 0.0
  dpen <- function(x) rep(0,length(x))
  ddpen <- function(x) matrix(0,length(x),length(x))
  # nugget penalty
  pen_nug <- function(x,j){ pen_loglo(x[j],.1,log(0.001)) }
  dpen_nug <- function(x,j){
    dpen <- rep(0,length(x))
    dpen[j] <- dpen_loglo(x[j],.1,log(0.001))
    return(dpen)
  }
  ddpen_nug <- function(x,j){
    ddpen <- matrix(0,length(x),length(x))
    ddpen[j,j] <- ddpen_loglo(x[j],.1,log(0.001))
    return(ddpen)
  }
  # smoothness penalty
  pen_sm <- function(x,j){ pen_loglo(x[j],.1,log(0.2)) }
  dpen_sm <- function(x,j){
    dpen <- rep(0,length(x))
    dpen[j] <- dpen_loglo(x[j],.1,log(0.2))
    return(dpen)
  }
  ddpen_sm <- function(x,j){
    ddpen <- matrix(0,length(x),length(x))
    ddpen[j,j] <- ddpen_loglo(x[j],.1,log(0.2))
    return(ddpen)
  }
  # variance penalty
  # dangerous because vv could get redefined
  pen_var <- function(x,j){ pen_hi(x[j]/vv,1,6) }
  dpen_var <- function(x,j){
    dpen <- rep(0,length(x))
    dpen[j] <- 1/vv*dpen_hi(x[j]/vv,1,6)
    return(dpen)
  }
  ddpen_var <- function(x,j){
    ddpen <- matrix(0,length(x),length(x))
    ddpen[j,j] <- 1/vv^2*ddpen_hi(x[j]/vv,1,6)
    return(ddpen)
  }
  
  # penalty on large smoothness parameters
  pen_sm_hi <- function(x,j){sm <-8.0; bb<-0.5; tt<-3.0; pen_hi(x[j]/bb,tt,sm) }
  dpen_sm_hi <- function(x,j){
    sm <- 8.0
    bb <- 0.5
    tt <- 3.0
    dpen <- rep(0,length(x))
    dpen[j] <- 1/bb*dpen_hi(x[j]/bb,tt,sm)
    return(dpen)
  }
  ddpen_sm_hi <- function(x,j){
    sm <- 8.0
    bb <- 0.5
    tt <- 3.0
    ddpen <- matrix(0,length(x),length(x))
    ddpen[j,j] <- 1/bb^2*ddpen_hi(x[j]/bb,tt,sm)
    return(ddpen)
  }
  
  return( list( pen = pen, dpen = dpen, ddpen = ddpen ) )
}

