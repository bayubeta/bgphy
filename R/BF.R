

#================== Naive Monte Carlo integral approximation

# # marginal likelihood
# logMargLik <- function(X, tree, model, priors, nsamples){
#
#   # create g(theta), priors sampler
#   g <- lapply(priors, priorsampler)
#
#   # parameter names
#   names <- names(g)
#
#   # vector of log-likelihood
#   loglikv <- numeric(nsamples)
#
#   op <- options(PCMBase.Raise.Lik.Errors = FALSE)
#   on.exit(options(op))
#
#   for (i in 1:nsamples){
#     repeat{
#       # sample prior
#       pr <- sapply(g, function(f){f(1)})
#       # likelihood from priors
#       loglik <- PCMloglik(X = X, tree = tree, model = model, pr)
#       if (!is.na(loglik)){
#         break()
#       }
#     }
#     loglikv[i] <- loglik
#   }
#
#
#   # apply logsumexp
#   lmax <- max(loglikv)
#   loglikv_shifted <- loglikv - lmax
#   return(-log(nsamples) + lmax + log(sum(exp(loglikv_shifted))))
# }
#
#
#
# # ================ Generalized Stepping Stone Method ================
#
# # reference distribution from posterior samples
# # takes prior distributions for the bounds info
# # takes posterior samples for the parameters
#
# refDist <- function(P, priors){
#
#   # calculate sample mean and variance
#   mu_hat <- colMeans(P)
#   sigma2_hat <- apply(P, 2, var)
#
#   # set reference dist. as a prior object similar to priors
#   rd <- priors
#
#   # check bounds of each parameter
#   for (name in names(rd)){
#     # if unbounded, use normal
#     if ( all(is.infinite(attr(priors[[name]], "bounds"))) ){
#       rd[[name]] <- prior_normal(mean = mu_hat[[name]], sd = sqrt(sigma2_hat[[name]]))
#     }else if (attr(priors[[name]], "bounds")[1] == 0 & is.infinite(attr(priors[[name]], "bounds")[2])){
#       # if > 0, use gamma
#       mh <- mu_hat[[name]]
#       s2 <- sigma2_hat[[name]]
#       rd[[name]] <- prior_gamma(shape = (mh^2)/s2, rate = mh/s2)
#     }
#   }
#
#   return(rd)
# }
#
#
#
#
# # # calculate log ratio, rk
# # logrk <- function(X, tree, model, priors, P, nsamples = 1000, b = NULL){
# #
# #   # create a reference distribution
# #   rd <- refDist(P, priors)
# #
# #   # run importance sampler sequentially, from b = 0 to b = 1
# #   if (is.null(b)){
# #     b <- seq(0, 1, length.out = 11)
# #   }
# #
# #   # sample parameters from reference dist on k = 0
# #   Th <- Prio
# #
# #   for (k in 2:length(b)){
# #     # calculate likelihood of the samples
# #
# #
# #     # prior pdf of samples
# #     prior_sampler()
# #   }
# # }
#
#
#
# # ================ Importance Sampling ================
# # use importance sampling with reference distribution as the proposal
#
# logMargLik_IS <- function(model, X, tree, priors, posterior, nsample = 10000, progress = TRUE){
#
#   if (progress){timestamp()}
#
#   # create a reference distribution
#   rd <- refDist(posterior, priors)
#
#   # sample values from reference dist. as proposals
#   Pr <- prior_sampler(rd)(nsample)
#
#   # log p_ref
#   log_rd <- colSums(apply(Pr, 1, function(p){mapply(function(f, x){f(x)}, rd, p)}))
#
#   # log priors
#   log_p <- colSums(apply(Pr, 1, function(p){mapply(function(f, x){f(x)}, priors, p)}))
#
#   # log posteriors
#   log_post <- apply(Pr, 1, function(p){PCMloglik(X, tree, model, p)})
#
#   # aggregate
#   q <- log_post + log_p - log_rd
#
#   N <- dim(Pr)[1]
#
#   if (progress){timestamp()}
#
#   return(-log(N) + logsumexp(q))
# }



