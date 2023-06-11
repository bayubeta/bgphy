
#' @export
rwm <- function(model, X, tree, priors, initial, iter, burn, scale = 0.1, progress = TRUE, ncheck = FALSE){

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
  priors_tr <- varchange(priors)
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

  # set a starting point for the adaptive phase
  adapt_start = 20*d + 1


  for (i in 2:iter){

    # random walk move
    if (i <= 20*d){
      # fixed MCMC move for i in 1:2d
      pars1 <- fixedMCMC(pars0, scale)
    }
    else if (i == adapt_start){
      # start the adaptive phase
      # compute the sample covariance for the current time (n = i-1)
      Sn <- cov(P[1:(i-1),])

      # compute the sample mean for the current time (n = i-1)
      mu_n <- colMeans(P[1:(i-1),])

      # adaptive MCMC move, sample proposal values
      pars1 <- adaptMCMC(pars0, Sn, 0.05, scale)

    }
    else{
      # adaptive MCMC move for i in adapt_start+1:iter
      n = i - 1
      # update the sample mean for the current time (n = i-1)
      mu_n <- (n*mu_n + pars0)/(n+1) # i is the next iteration, n+1

      # update the sample covariance for the current time
      Sn = covUpdate(pars0, Sn, mu_n, n)

      # adaptive MCMC move, sample proposal values
      pars1 <- adaptMCMC(pars0, Sn, 0.05, scale)
    }


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


  # ================= transform samples back =================
  P <- t(apply(P, 1, tr$g))


  # acceptance rate
  attr(P, "accept") <- mean(ar)


  return(P)

}
