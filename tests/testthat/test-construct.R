
test_that("Model is correctly constructed", {
  ### Global model
  expect_silent(OU1 <- setModel(tree = lizardTree, regime_names = "Regime1",
                                modeltypes = "OU"))
  # test the model
  expect_s3_class(OU1$model, "OU")
  expect_contains(names(OU1$model), c("X0", "H", "Theta", "Sigma_x"))
  # test the tree
  expect_s3_class(OU1$tree, "PCMTree")
  expect_named(OU1$tree$part.regime, "101")
  expect_equal(OU1$tree$part.regime[[1]], "Regime1")
  # test the priors
  expect_s3_class(OU1$priors, "bgphy_priors")
  expect_equal(names(OU1$priors), c("X0", "alpha", "theta", "sigma"))
  expect_s3_class(OU1$priors$X0, c("priorpdf", "normal"))
  expect_s3_class(OU1$priors$alpha, c("priorpdf", "halfnormal"))
  expect_s3_class(OU1$priors$theta, c("priorpdf", "normal"))
  expect_s3_class(OU1$priors$sigma, c("priorpdf", "halft"))

  ### Mixed model
  # BM to OU
  expect_silent(BMOU <- setModel(tree = lizardTree,
                                 regime_names = c("Ancestral", "New"),
                                 modeltypes = c("BM", "OU"),
                                 startNodes = list(Ancestral = c(101),
                                                   New = c(135))))
  # test the model
  expect_s3_class(BMOU$model, "MixedGaussian")
  expect_s3_class(BMOU$model$X0, "_Global")
  expect_s3_class(BMOU$model$Ancestral, "BM")
  expect_s3_class(BMOU$model$New, "OU")
  # test the tree
  expect_s3_class(BMOU$tree, "PCMTree")
  expect_named(BMOU$tree$part.regime, c("101", "135"))
  expect_equal(as.vector(BMOU$tree$part.regime), c("Ancestral", "New"))
  # test the priors
  expect_s3_class(BMOU$priors, "bgphy_priors")
  expect_equal(names(BMOU$priors), c("X0", "sigma_1", "alpha_2",
                                        "theta_2", "sigma_2"))
  expect_s3_class(BMOU$priors$X0, c("priorpdf", "normal"))
  expect_s3_class(BMOU$priors$sigma_1, c("priorpdf", "halft"))
  expect_s3_class(BMOU$priors$alpha_2, c("priorpdf", "halfnormal"))
  expect_s3_class(BMOU$priors$theta_2, c("priorpdf", "normal"))
  expect_s3_class(BMOU$priors$sigma_2, c("priorpdf", "halft"))


  # BM to BM
  expect_silent(BMBM <- setModel(tree = lizardTree,
                                 regime_names = c("R1", "R2"),
                                 modeltypes = c("BM", "BM"),
                                 startNodes = list(R1 = c(101),
                                                   R2 = c(185))))
  # test the model
  expect_s3_class(BMBM$model, "MixedGaussian")
  expect_s3_class(BMBM$model$X0, "_Global")
  expect_s3_class(BMBM$model$R1, "BM")
  expect_s3_class(BMBM$model$R2, "BM")
  # test the tree
  expect_s3_class(BMBM$tree, "PCMTree")
  expect_named(BMBM$tree$part.regime, c("101", "185"))
  expect_equal(as.vector(BMBM$tree$part.regime), c("R1", "R2"))
  # test the priors
  expect_s3_class(BMBM$priors, "bgphy_priors")
  expect_equal(names(BMBM$priors), c("X0", "sigma_1", "sigma_2"))
  expect_s3_class(BMBM$priors$X0, c("priorpdf", "normal"))
  expect_s3_class(BMBM$priors$sigma_1, c("priorpdf", "halft"))
  expect_s3_class(BMBM$priors$sigma_2, c("priorpdf", "halft"))


  # OU to BM
  expect_silent(OUBM <- setModel(tree = lizardTree,
                                 regime_names = c("A", "B"),
                                 modeltypes = c("OU", "BM"),
                                 startNodes = list(B = c(160),
                                                   A = c(101))))
  # test the model
  expect_s3_class(OUBM$model, "MixedGaussian")
  expect_s3_class(OUBM$model$X0, "_Global")
  expect_s3_class(OUBM$model$A, "OU")
  expect_s3_class(OUBM$model$B, "BM")
  # test the tree
  expect_s3_class(OUBM$tree, "PCMTree")
  expect_named(OUBM$tree$part.regime, c("101", "160"))
  expect_equal(as.vector(OUBM$tree$part.regime), c("A", "B"))
  # test the priors
  expect_s3_class(OUBM$priors, "bgphy_priors")
  expect_equal(names(OUBM$priors), c("X0", "alpha_1", "theta_1",
                                     "sigma_1", "sigma_2"))
  expect_s3_class(OUBM$priors$X0, c("priorpdf", "normal"))
  expect_s3_class(OUBM$priors$alpha_1, c("priorpdf", "halfnormal"))
  expect_s3_class(OUBM$priors$theta_1, c("priorpdf", "normal"))
  expect_s3_class(OUBM$priors$sigma_1, c("priorpdf", "halft"))
  expect_s3_class(OUBM$priors$sigma_2, c("priorpdf", "halft"))



  # OU to OU
  expect_silent(OUOU <- setModel(tree = lizardTree,
                                 regime_names = c("1", "2"),
                                 modeltypes = c("OU", "OU"),
                                 startNodes = list("1" = c(101, 162, 143, "ahli"),
                                                   "2" = c(134, 105))))
  # test the model
  expect_s3_class(OUOU$model, "MixedGaussian")
  expect_s3_class(OUOU$model$X0, "_Global")
  expect_s3_class(OUOU$model$"1", "OU")
  expect_s3_class(OUOU$model$"2", "OU")
  # test the tree
  expect_s3_class(OUOU$tree, "PCMTree")
  expect_contains(names(OUOU$tree$part.regime), c(101, 162, 143, "ahli", 134, 105))
  expect_equal(as.vector(OUOU$tree$part.regime[c(101, 162, 143, "ahli", 134, 105)]),
               c("1", "1", "1", "1", "2", "2"))
  # test the priors
  expect_s3_class(OUOU$priors, "bgphy_priors")
  expect_equal(names(OUOU$priors),
                  c("X0", "alpha_1", "theta_1", "sigma_1",
                    "alpha_2", "theta_2", "sigma_2"))
  expect_s3_class(OUOU$priors$X0, c("priorpdf", "normal"))
  expect_s3_class(OUOU$priors$alpha_1, c("priorpdf", "halfnormal"))
  expect_s3_class(OUOU$priors$theta_1, c("priorpdf", "normal"))
  expect_s3_class(OUOU$priors$sigma_1, c("priorpdf", "halft"))
  expect_s3_class(OUOU$priors$alpha_2, c("priorpdf", "halfnormal"))
  expect_s3_class(OUOU$priors$theta_2, c("priorpdf", "normal"))
  expect_s3_class(OUOU$priors$sigma_2, c("priorpdf", "halft"))
  })




test_that("Throw an error for wrong inputs", {
  # tree
  expect_error(setModel(tree = "Tree", regime_names = "Regime1", modeltypes = "OU"))

  # model types: reject anything other than "BM" and "OU"
  expect_error(setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "ou"))
  expect_error(setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "bm"))
  expect_error(setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "type"))
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "oU"),
                        startNodes = list(Ancestral = c(101), New = c(135))))

  # startNodes:
  # needs to be a list
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "OU"),
                        startNodes = c("101", "135")))
  # names match regime_names
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "OU"),
                        startNodes = list(R1 = c(101), R2 = c(135))))
  # number of regimes match
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "OU"),
                        startNodes = list(Ancestral = c(101))))
  expect_error(setModel(tree = lizardTree, regime_names = "Ancestral",
                        modeltypes = c("BM", "OU"),
                        startNodes = list(Ancestral = c(101), New = c(135))))
  # values in startNodes must be in the tree
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "OU"),
                        startNodes = list(Ancestral = c(101), New = c(135, "abc"))))
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "OU"),
                        startNodes = list(Ancestral = c(101, 99), New = c(135))))
  # ancestral nodes must be included in startNodes
  expect_error(setModel(tree = lizardTree, regime_names = c("Ancestral", "New"),
                        modeltypes = c("BM", "OU"),
                        startNodes = list(Ancestral = c(102), New = c(135))))
  # number of regimes must match the length of modeltypes
  expect_error(setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = c("OU", "OU")))
})


test_that("Outputs are correct",{
  OU1 <- setModel(tree = lizardTree, regime_names = "R1", modeltypes = "OU")
  expect_output(print(OU1), "R1        OU\\(alpha, theta, sigma\\)")
  expect_output(print(OU1), "X0 ~ normal\\(mean = 0, sd = 10\\)")
  expect_output(print(OU1), "alpha ~ halfnormal\\(sigma = 6.931472\\)")
  expect_output(print(OU1), "theta ~ normal\\(mean = 0, sd = 10\\)")
  expect_output(print(OU1), "sigma ~ halft\\(nu = 1, sigma = 3\\)")

  BM1 <- setModel(tree = lizardTree, regime_names = "Ancestral", modeltypes = "BM")
  expect_output(print(BM1), "Ancestral   BM\\(sigma\\)")
  expect_output(print(BM1), "X0 ~ normal\\(mean = 0, sd = 10\\)")
  expect_output(print(BM1), "sigma ~ halft\\(nu = 1, sigma = 3\\)")

  expect_silent(OUOU <- setModel(tree = lizardTree,
                                 regime_names = c("1", "2"),
                                 modeltypes = c("OU", "OU"),
                                 startNodes = list("1" = c(101, 162, 143, "ahli"),
                                                   "2" = c(134, 105))))
  expect_output(print(OUOU), "1         OU\\(alpha_1, theta_1, sigma_1\\)")
  expect_output(print(OUOU), "2         OU\\(alpha_2, theta_2, sigma_2\\)")
  expect_output(print(OUOU), "alpha_1 ~ halfnormal\\(sigma = 6.931472\\)")
  expect_output(print(OUOU), "theta_1 ~ normal\\(mean = 0, sd = 10\\)")
  expect_output(print(OUOU), "sigma_1 ~ halft\\(nu = 1, sigma = 3\\)")
  expect_output(print(OUOU), "alpha_2 ~ halfnormal\\(sigma = 6.931472\\)")
  expect_output(print(OUOU), "theta_2 ~ normal\\(mean = 0, sd = 10\\)")
  expect_output(print(OUOU), "sigma_2 ~ halft\\(nu = 1, sigma = 3\\)")
})
