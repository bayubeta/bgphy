# ========================== helper functions ==========================

# logit functions
logit <- function(u){
  log(u/(1 - u))
}

# inverse logit
invlogit <- function(v){
  1/(1 + exp(-v))
}


# load parameters into a PCM model
loadParams <- function(o, p){
  # load parameters into a PCM model
  # similar functionality to PCMBase::PCMParamLoadOrStore, load = TRUE

  # o: PCM object
  # p: vector of parameters

  r <- length(attr(o, "regimes")) # regimes
  pos <- 1 # index position

  if (r == 1){
    # if only one regime
    for (name in names(o)){
      o[[name]][] <- p[pos]
      pos = pos + 1
    }
  }
  else{
    # if there are multiple regimes
    for (name in names(o)){
      plength <- length(o[[name]])

      if (pos >= attr(o, "p")){
        break
      }

      if (plength == 1){
        o[[name]][[1]][] <- p[pos]
        pos = pos + 1
      }
      else{
        for (i in 1:plength){
          o[[name]][[i]][] <- p[pos]
          pos = pos + 1
        }
      }
    }
  }
  return(o)
}

# get parameter names from model
getParamNames <- function(model){
  rnames <- attr(model$model, "regimes") # get regime names
  r <- length(rnames) # get number of regimes
  modeltypes <- attr(model, "modeltypes") # model types

  # list of parameter names
  parnames <- c()

  # check if X0 parameter is present
  if(names(model$mode)[1] == "X0"){parnames[1] <- "X0"}

  # iterate over regimes
  if (r == 1){
    parnames <- c(parnames, paramNames(modeltypes))
  }
  else{
    for (i in 1:r){
      parnames <- c(parnames, paramNames(modeltypes[i], i))
    }
  }
  return(parnames)
}

# print parameter names given a model type (optionally with index)
paramNames <- function(modeltype, i = NULL){
  stopifnot("Model name not found! Choose between BM or OU." = modeltype %in% c("BM", "OU"))
  if (modeltype == "BM"){
    if (is.null(i)){
      return("sigma")
    }else{
      return(paste0("sigma_", i))
    }
  }
  if (modeltype == "OU"){
    if (is.null(i)){
      return(c("alpha", "theta", "sigma"))
    }else{
      return(c(paste0("alpha_", i), paste0("theta_", i),  paste0("sigma_", i)))
    }
  }
}


# print model with its parameters (optionally with index)
modelprint <- function(modeltype, i = NULL){
  stopifnot("Model name not found! Choose between BM or OU." = modeltype %in% c("BM", "OU"))
  if (modeltype == "BM"){
    if (is.null(i)){
      cat(paste0("BM(sigma)"))
    }else{
      cat(paste0("BM(sigma_", i, ")"))
    }
  }
  if (modeltype == "OU"){
    if (is.null(i)){
      cat(paste0("OU(alpha, theta, sigma)"))
    }else{
      cat(paste0("OU(alpha_", i, ", theta_", i, ", sigma_", i, ")"))
    }
  }
}


# log-sum-exp function
logsumexp <- function(logW, log = TRUE){
  # calculate the (log) of sum of W from logW data
  # log(sum(W)) = log(W[1] + ... + W[N])
  # = log(exp(logW[1]) + ... + exp(logW[N]))
  max_logW <- max(logW)
  logW_shifted <- logW - max_logW
  logsumW <- log(sum(exp(logW_shifted))) + max_logW
  if (log == TRUE){
    return(logsumW)
  }else{
    return(exp(logsumW))
  }
}
