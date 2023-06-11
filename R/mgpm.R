
# @export
# mgpm <- function(model, X, tree, priors, iter = 10000, burn = 0.1*iter, method = "IS", initial = NULL, progress = TRUE, ...){
#   match.arg(method, c("AM"))
#
#   # burn-in period needs to be less than iter
#
#   # initial position, sample from priors if not provided
#   if (is.null(initial)){
#     initial <- prior_sampler(priors)(1)
#   }
#
#   if (progress){timestamp()}
#
#   P <- rwm(model = model, X = X, tree = tree, priors = priors,
#            initial = initial, iter = iter + burn, burn = burn, progress = progress, ...)
#
#   # acceptance rate
#   ar <- attr(P, "accept")
#
#
#   # cut the burn-in samples
#   P <- P[-c(1:burn),]
#
#
#   attr(P, "accept") <- ar
#   class(P) <- "mgpm_posterior"
#
#   if (progress){timestamp()}
#
#
#   return(P)
# }



#' @export
mgpm <- function(model, X, tree, priors, nsample = 10000, initial = NULL, scale = 1, parallel = TRUE){

  # initial position, sample from priors if not provided
  if (is.null(initial)){
    initial <- prior_sampler(priors)(1)
  }

  if (is.null(dim(X))){X <- matrix(X, nrow = 1)}

  res <- IS(model = model, X = X, tree = tree, priors = priors,
            initial = initial, nsample = nsample, scale = scale, parallel = parallel)

  class(res) <- "mgpm_posterior"

  return(res)
}
