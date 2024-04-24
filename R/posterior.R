# ------------------------- Inference functions -------------------------
# function to run the Importance Sampling method
IS <- function(model, X, nsample, scale, parallel = TRUE){
  #================ Laplace approximation + Importance sampling ================

  # number of parameters
  d <- length(model$priors)

  # ================== change of variables routine (to log space)
  ############
  # priors_tr: f_X to f_Y
  # tr: X to Y and vice versa
  # apply change of variables to priors into unbounded space
  priors_tr <- lapply(model$priors, prior_transform)

  ############
  # define tr, functions to transform parameters
  # tr$f: to unbounded
  # tr$g: from unbounded
  tr <- trfunc(priors_tr)


  repeat{ # in case of singular Hessian, repeat from different starting point
    # find initial position for optim(), by the mean of the prior distribution
    init <- apply(prior_sampler(model$priors)(1000), 2, stats::median) # set initial value as median of samples drawn from priors

    # transform to unbounded space
    initial <- tr$f(init)

    # simplify lupost as a function of parameters
    lu_post <- function(p){lupost(p, model$model, X, model$tree, priors_tr, tr)[[1]]}

    # ================== begin the Laplace approximation routine
    # search posterior mode
    optRes <- stats::optim(par = initial, fn = lu_post, method = "BFGS",
                           control = list(fnscale=-1), hessian = TRUE)
    # posterior mode (log space)
    post_mode <- optRes$par
    # approximated covariance, scale to focus on the important area around mode
    appr_cov <- try(solve(-optRes$hessian)*scale + diag(d)*1e-8, silent = TRUE) # add diagonal jitters for numerical stability

    # check if appr_cov is already positive definite
    if (all(eigen(appr_cov)$values >= 1e-13) & all(class(appr_cov) != "try-error")){
      break()
    }
  }


  # ================== begin the (self-normalized) Importance sampling routine
  # sample from the mixture of normal distribution
  q <- mvtnorm::rmvnorm(nsample, mean = post_mode, sigma = appr_cov)

  if (parallel){
    cl <- parallel::makeCluster(parallel::detectCores(),"PSOCK")
    parallel::clusterExport(cl, varlist = c("PCMloglik", "setParams", "loadParams"), envir = environment())
    # log-unnormalized posterior
    logp <- t(parallel::parApply(cl, q, 1, lupost,
                                 model$model, X, model$tree, priors_tr, tr)) # log-unnormalized posterior & likelihood
    logp <- as.vector(logp) # convert to vector

    # log proposal
    logq <- parallel::parApply(cl, q, 1, mvtnorm::dmvnorm,
                               mean = post_mode, sigma = appr_cov, log = TRUE) # log density of normal
    on.exit(parallel::stopCluster(cl))

  }else{
    # log-unnormalized posterior
    logp <- t(apply(q, 1, lupost, model$model, X, model$tree, priors_tr, tr)) # log-unnormalized posterior & likelihood
    logp <- as.vector(logp) # convert to vector
    logq <- apply(q, 1, mvtnorm::dmvnorm,
                  mean = post_mode, sigma = appr_cov, log = TRUE) # log density of normal
  }

  logW <- logp - logq # log of weights
  # W <- exp(logW) # weights

  logWst <- logW - logsumexp(logW) # log of normalized weights
  Wst <- exp(logWst) # normalized weights



  # ------------------------ second pass ------------------------



  # resample (with replacement) using SIR, use normalized weights as the prob.
  new_idx <- sample(1:nsample, nsample, prob = Wst, replace = TRUE)

  # newly resampled particles
  q_new <- q[new_idx,]

  # ======== Sample from a MVN, where the marginal means are q_new
  # covariance of the perturbation
  S_p <- diag(apply(q_new, 2, var))

  # perturbations
  e <- mvtnorm::rmvnorm(nsample, mean = rep(0,d), sigma = S_p/nsample)
  q <- q_new + e



  if (parallel){
    cl <- parallel::makeCluster(parallel::detectCores(),"PSOCK")
    parallel::clusterExport(cl,
                            varlist = c("PCMloglik", "setParams", "loadParams", "q_new", "q", "S_p"),
                            envir = environment())

    # log-unnormalized posterior
    logp <- t(parallel::parApply(cl, q, 1, lupost,
                                 model$model, X, model$tree, priors_tr, tr)) # log-unnormalized posterior & likelihood
    logp <- as.vector(logp) # convert to vector

    # log proposal
    logq <- parallel::parSapply(cl, 1:nsample, function(i){mvtnorm::dmvnorm(q[i,], mean = q_new[i,], sigma = S_p, log = TRUE)})

    on.exit(parallel::stopCluster(cl))

  }else{
    # log-unnormalized posterior
    logp <- t(apply(q, 1, lupost, model$model, X, model$tree, priors_tr, tr)) # log-unnormalized posterior & likelihood
    logp <- as.vector(logp) # convert to vector
    logq <- sapply(1:nsample, function(i){mvtnorm::dmvnorm(q[i,], mean = q_new[i,], sigma = S_p, log = TRUE)}) # log density of normal
  }


  logW <- logp - logq # log of weights
  # W <- exp(logW) # weights

  logWst <- logW - logsumexp(logW) # log of normalized weights
  Wst <- exp(logWst) # normalized weights


  Q <- t(apply(q, 1, tr$g)) # proposed values on the original space

  return(list(Q = Q, W = Wst))
}


# function to estimate (marginal) quantiles given sample and weights
est_quantiles <- function(Q, W, probs = c(0.025, 0.5, 0.975)){

  # number of params
  d <- dim(Q)[2]

  # matrix of quantiles
  qM <- matrix(nrow = d, ncol = length(probs))

  for (i in 1:d){
    # sort the parameter values
    qSort <- sort(Q[,i], index.return = TRUE)
    # sorted parameter values
    q_sorted <- qSort$x
    # indices
    idx <- qSort$ix

    # compute quantiles
    qtls <- sapply(probs, function(p){min(which(cumsum(W[idx]) >= p))})

    # store to matrix
    qM[i,] <- q_sorted[qtls]
  }

  colnames(qM) <- probs

  return(qM)
}


# model evaluation using posterior predictive deviance
# takes the normalized weights (W) and samples from the proposal distribution (Q)
ppred_loss <- function(P, X, model){
  # P: particles with weights (W) and values from proposal distribution (Q)
  # X: observed values at the tips
  # model: bgphy_model class

  # simulate from the posterior predictive function
  ntips <- dim(X)[2]

  # sample from posterior
  S <- length(P$W) # number of samples
  post_id <- sample(1:S, S, replace = TRUE, prob = P$W)
  post_samples <- P$Q[post_id,]

  # and simulate using PCMBase
  PCMsimulate <- function(p){
    PCMBase::PCMSim(model$tree, setParams(p, model$model), p[1])[,1:ntips]
  }

  Xppred <- matrix(nrow = S, ncol = ncol(X))
  cl <- parallel::makeCluster(parallel::detectCores(),"PSOCK")
  parallel::clusterExport(cl, varlist = c("PCMloglik", "setParams", "loadParams",
                                          "model", "X", "ntips"), envir = environment())
  Xppred <- t(parallel::parApply(cl, post_samples, 1, PCMsimulate))
  on.exit(parallel::stopCluster(cl))

  # return a list of loss score and the posterior predictive samples
  loss <- sum((X - colMeans(Xppred))^2) + sum(diag(stats::cov(Xppred)))

  return(list(loss = loss, Xppred = Xppred))
}






# ------------------------- Helper functions -------------------------

# assign a vector of parameters to a model
setParams <- function(p, model){
  # p: all parameters, excluding Sigmae_x
  # p := c(X0, H1, Theta1, Sigma_x1, ..., Hr, Thetar, Sigma_xr, Sigmae_x)
  # Load parameters into the model

  # # create an empty vector of the size of total parameters count
  # v <- numeric(attr(model, "p"))
  #
  # # fill everything except the last with p
  # v[-length(v)] <- p

  # load p into the model
  model <- loadParams(model, p)

  return(model)
}


# log-likelihood given data, tree, and model with parameter values
PCMloglik <- function(X, tree, model, p){
  # model with parameters p
  prop_model <- setParams(p, model)

  # log likelihood calculation
  loglik <- PCMBase::PCMLik(X, tree, prop_model, metaI = PCMBaseCpp::PCMInfoCpp)

  return(loglik[1])
}

# log-unnormalized-posterior
lupost <- function(p, model, X, tree, priors_tr, tr){

  # p: vector of parameters, c(X0, alpha_1, theta_1, sigma_1, ..., alpha_r, theta_r, sigma_r)
  # model: PCM model
  # tree: phylo object
  # priors_tr: list of priors pdf on unbounded space
  # tr: transformation functions, obtained from trfunc()

  if (is.null(dim(X))){X <- matrix(X, nrow = 1)}

  # change parameters to p_b (bounded),
  # because PCMLik (the likelihood) needs parameter values from the original space
  p_b <- tr$g(p)

  # do not raise warning
  op <- options(PCMBase.Raise.Lik.Errors = FALSE)
  on.exit(options(op))

  # ===== log likelihood
  loglik <- PCMloglik(X, tree, model, p_b)

  # ===== log priors
  # calculated in the unbounded space
  log_p <- sum(mapply(function(f, x){f(x)}, priors_tr, p))

  # ===== sum of log-likelihood and log priors
  sum_log <- loglik + log_p

  # if infinite, return a very low value
  if (is.infinite(sum_log) || is.na(sum_log)){
    return(-1e20)
  }

  # return log unnormalized posterior

  return(sum_log)
}

