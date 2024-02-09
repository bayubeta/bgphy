test_that("setParams works", {
  # global BM
  expect_silent(BM1 <- setParams(c(1,1), PCMBase::PCM("BM")))
  expect_s3_class(BM1, c("BM", "GaussianPCM", "PCM"))
  expect_equal(BM1$X0[1], 1)
  expect_equal(BM1$Sigma_x[[1]], 1)

  # global OU
  expect_silent(newModel <- setParams(c(1,1,1,1), modelOU))
  expect_s3_class(newModel, c("OU", "GaussianPCM", "PCM"))
  expect_equal(newModel$X0[1], 1)
  expect_equal(newModel$H[[1]], 1)
  expect_equal(newModel$Sigma_x[[1]], 1)
  expect_equal(newModel$Theta[[1]], 1)

  # OUOU
  expect_silent(newMixed <- setParams(c(1,2,2,2,3,3,3), modelMixedOU))
  expect_s3_class(newModel, c("MixedGaussian", "GaussianPCM", "PCM"))
  expect_equal(newMixed$X0[1], 1)
  expect_equal(newMixed$"1"$H[[1]], 2)
  expect_equal(newMixed$"1"$Sigma_x[[1]], 2)
  expect_equal(newMixed$"1"$Theta[[1]], 2)
  expect_equal(newMixed$"2"$H[[1]], 3)
  expect_equal(newMixed$"2"$Sigma_x[[1]], 3)
  expect_equal(newMixed$"2"$Theta[[1]], 3)

  # BMOU
  BMOU <- setModel(tree = lizardTree,
                   regime_names = c("Ancestral", "New"),
                   modeltypes = c("BM", "OU"),
                   startNodes = list(Ancestral = c(101), New = c(135)))
  expect_silent(newBMOU <- setParams(c(1,2,3,3,3), BMOU$model))
  expect_s3_class(newBMOU, c("MixedGaussian", "GaussianPCM", "PCM"))
  expect_equal(newBMOU$X0[1], 1)
  expect_equal(newBMOU$Ancestral$Sigma_x[[1]], 2)
  expect_equal(newBMOU$New$H[[1]], 3)
  expect_equal(newBMOU$New$Sigma_x[[1]], 3)
  expect_equal(newBMOU$New$Theta[[1]], 3)
})


test_that("likelihood and posterior works", {
  # check that likelihoods are the same when using star tree
  mytr <- ape::read.tree(text = paste0("(",
                                       paste0(sapply(1:100,
                                                     function(i){paste0("t", i, ":2")}),
                                              collapse = ","), ");"))
  set.seed(1234, kind = "L'Ecuyer-CMRG"); X <- rnorm(100, 0, sqrt(8))
  # BM model
  BM2 <- setParams(c(0,2), PCMBase::PCM("BM"))
  expect_equal(PCMloglik(matrix(X, nrow = 1, dimnames = list(NULL, mytr$tip.label)), mytr, BM2, p = c(0,2)),
               sum(stats::dnorm(X, 0, sqrt(8), log = TRUE)))
  # OU model
  OU2 <- setParams(c(0,2,2,2), PCMBase::PCM("OU"))
  expect_equal(PCMloglik(matrix(X, nrow = 1, dimnames = list(NULL, mytr$tip.label)), mytr, OU2, p = c(0,2,2,2)),
               sum(stats::dnorm(X, (1-exp(-2*2))*2, sqrt(1-exp(-2*2*2)), log = TRUE)))

  # check that likelihoods are the same when using a simple bifurcated tree
  mytr2 <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
  # BM
  C <- matrix(c(2,1,0,0,
                1,2,0,0,
                0,0,2,1,
                0,0,1,2), nrow = 4, byrow = TRUE)
  set.seed(1234, kind = "L'Ecuyer-CMRG"); X2 <- rnorm(4, 0, sqrt(8))
  expect_equal(PCMloglik(matrix(X2, nrow = 1, dimnames = list(NULL, mytr2$tip.label)), mytr2, BM2, p = c(0,2)),
               mvtnorm::dmvnorm(X2, mean = rep(0,4), sigma = 2^2 * C, log = TRUE))
  #OU
  V <- diag(1-exp(-2*2*2), 4)
  V[1,2] <- V[2,1] <- V[3,4] <- V[4,3] <- exp(-4) - exp(-8)
  expect_equal(PCMloglik(matrix(X2, nrow = 1, dimnames = list(NULL, mytr2$tip.label)), mytr2, OU2, p = c(0,2,2,2)),
               mvtnorm::dmvnorm(X2, mean = rep((1-exp(-2*2))*2, 4), sigma = V, log = TRUE))


  # check log unnormalized posteriors
  OU1 <- setModel(tree = mytr2, regime_names = "Regime1", modeltypes = "OU")
  priors_tr_OU <- lapply(OU1$priors, prior_transform)
  tr_OU <- trfunc(priors_tr_OU)
  expect_equal(lupost(p = tr_OU$f(c(0,2,2,2)), model = OU1$model, X = matrix(X2, nrow = 1, dimnames = list(NULL, mytr2$tip.label)),
                      tree = OU1$tree, priors_tr = priors_tr_OU, tr = tr_OU),
               c(stats::dnorm(0, mean = 0, sd = 10, log = TRUE) + stats::dnorm(2, mean = 0, sd = 10, log = TRUE) +
                 extraDistr::dhnorm(2, sigma = 3.465736, log = TRUE) + log(2) +
                 extraDistr::dht(2, nu = 1, sigma = 6, log = TRUE) + log(2) +
                 mvtnorm::dmvnorm(X2, mean = rep((1-exp(-2*2))*2, 4), sigma = V, log = TRUE),
                 mvtnorm::dmvnorm(X2, mean = rep((1-exp(-2*2))*2, 4), sigma = V, log = TRUE)))
})






