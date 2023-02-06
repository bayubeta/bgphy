######################


logit <- function(u){
  log(u/(1 - u))
}

invlogit <- function(v){
  1/(1 + exp(-v))
}



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



PCMGetParamNames <- function(o){
  r <- length(attr(o, "regimes"))

  if (r == 1){
    return(names(o)[1:4])
  }
  else{
    parnames <- numeric(attr(o, "p"))
    parnames <- parnames[-length(parnames)]

    parnames[1] <- names(o)[1]
    filled <- 1

    for (i in 1:r){
      pnames <- sapply(names(o[[i+1]]), FUN = function(x){paste0(x, "_", i)})
      parnames[(filled + 1):(filled + length(pnames))] <- pnames
      filled <- filled + length(pnames)
    }
  }

  return(parnames)
}


logsumexp <- function(logW){
  # calculate the log of sum of W from logW data
  # log(sum(W)) = log(W[1] + ... + W[N])
  # = log(exp(logW[1]) + ... + exp(logW[N]))
  max_logW <- max(logW)
  logW_shifted <- logW - max_logW
  log(sum(exp(logW_shifted))) + max_logW
}




