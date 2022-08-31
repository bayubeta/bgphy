######################

logit <- function(u){
  log(u/(1 - u))
}

invlogit <- function(v){
  1/(1 + exp(-v))
}


