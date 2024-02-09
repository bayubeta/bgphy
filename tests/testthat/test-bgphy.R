test_that("bgphy() throws an error for invalid inputs", {
  OU1 <- setModel(tree = lizardTree, regime_names = "Regime1", modeltypes = "OU")
  # model and data needs to be proper
  expect_error(bgphy(OU1$model, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU)))))
  expect_error(bgphy(OU1, XOU[1,]))

  # names on the tips match the column names
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, sort(colnames(XOU))))))

  # valid value for nsample
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     nsample = "a"))
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     nsample = -1))
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     nsample = 1.4))

  # valid value for scale
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     scale = 0))
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     scale = -1))
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     scale = "asdf"))

  # valid value for parallel
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU))),
                     parallel = "TRUE"))


  # number of tips on the tree match number of data points
  mytr2 <- ape::read.tree(text = "((t1:1,t2:1):1,(t3:1,t4:1):1);")
  OU1 <- setModel(tree = mytr2, regime_names = "Regime1", modeltypes = "OU")
  expect_error(bgphy(OU1, matrix(XOU[1,], nrow = 1, dimnames = list(NULL, colnames(XOU)))))

})
