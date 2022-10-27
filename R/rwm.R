
#' @export
rwm <- function(model, X, tree, priors, initial, iter, scale = 0.01, progress = TRUE, ncheck = FALSE){

  # number of parameters
  d <- length(initial)

  # a vector of rejection rates, one for each step
  R <- log(stats::runif(iter))

  # a matrix of parameters, each row is a set of parameter for each step
  P <- matrix(NA, ncol = d, nrow = iter)
  colnames(P) <- PCMGetParamNames(model)


  ############
  # priors_tr: f_X to f_Y
  # tr: X to Y and vice versa

  # apply change of variables to priors into unbounded space
  priors_tr <- lapply(priors, varchange)
  ############


  ############
  # define tr, functions to transform parameters
  # tr$f: to unbounded
  # tr$g: from unbounded
  tr <- trfunc(priors_tr)
  ############

  # convert initial values into unbounded space
  initial <- tr$f(initial)

  # assign the first row as the initial parameter
  P[1,] <- initial


  # simplify lupost as a function of parameters
  lu_post <- function(p){lupost(p, model, X, tree, priors_tr, tr)}


  # ==================== for progress bar ====================
  if (progress){
    total <- iter
    pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = total)
    pb$tick(0)
    pb$tick(1)
    Sys.sleep(1/1000)
  }
  # ==========================================================


  # =================== for checkpoints ======================
  if (ncheck){
    icheck <- round(iter / ncheck)
  }
  # ==========================================================


  # record acceptance rate
  ar <- c(1)

  # ================= start RWMH algorithm ==================

  # save values for current state
  pars0 <- initial
  logpost0 <- lu_post(pars0)


  for (i in 2:iter){

    # random walk move
    pars1 <- mcmcmove(pars0, scale)

    # calculate log posterior on the new set of parameters
    logpost1 <- lu_post(pars1)

    # accept
    if (logpost1 - logpost0 > R[i]){
      P[i,] <- pars1
      pars0 <- pars1
      logpost0 <- logpost1

      ar[i] <- 1
    }
    else{
      # or reject
      P[i,] <- pars0

      ar[i] <- 0
    }

    # progress bar
    if (progress){
      pb$tick(1)
      Sys.sleep(1 / 1000)
    }

    # ================= check point =================
    if (ncheck && (i %% icheck == 0)){
      save(P, file = sprintf("check%02d.rda", as.integer(i/icheck)))
    }
  }


  class(P) <- "posterior"


  # ================= transform samples back =================
  P <- t(apply(P, 1, tr$g))


  attr(P, "accept") <- mean(ar)

  return(P)

}
